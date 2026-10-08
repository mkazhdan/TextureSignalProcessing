/*
Copyright (c) 2019, Michael Kazhdan and Fabian Prada
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

////////////////
// SquaredEDT //
////////////////

template< typename Real , unsigned int Dim >
template< typename BinaryType >
RegularGrid< Dim , unsigned int > SquaredEDT< Real , Dim >::Saito( const RegularGrid< Dim , BinaryType > &binaryGrid , bool verbose )
{
	Miscellany::PerformanceMeter pMeter( '.' );

	// Allocate the squaredEDT
	RegularGrid< Dim , unsigned int > squaredEDT;
	unsigned int threads = ThreadPool::NumThreads();
	{
		unsigned int res[Dim];
		for( int d=0 ; d<Dim ; d++ ) res[d] = binaryGrid.res(d);
		squaredEDT.resize( res );
		unsigned int max = 0;
		for( int d=0 ; d<Dim ; d++ ) max += res[d] * res[d];
		ThreadPool::ParallelFor( 0 , squaredEDT.resolution() , [&]( unsigned int , size_t i ){ squaredEDT[i] = max; } );
	}

	// Scan along the first axis
	pMeter.reset();
	{
		unsigned int dim = 0;
		unsigned int lineSize = binaryGrid.res(dim);
		size_t preLineCount=1 , postLineCount=1;
		for( unsigned int d=0 ; d<dim ; d++ ) preLineCount *= binaryGrid.res(d);
		for( unsigned int d=dim+1 ; d<Dim ; d++ ) postLineCount *= binaryGrid.res(d);
		size_t lineCount = preLineCount * postLineCount;

		// In parallel process the different lines
		ThreadPool::ParallelFor
			(
				0 , lineCount ,
				[&]( unsigned int , size_t i )
				{
					size_t pre = i % preLineCount , post = i / preLineCount;
					size_t offset = pre + post * ( preLineCount*lineSize );
					Pointer( unsigned int ) squaredEDTPtr = squaredEDT() + pre + post * ( preLineCount*lineSize );

					// Forward scan
					{
						bool first = true;
						unsigned int dist = 0;
						for( int j=0 ; j<(int)lineSize ; j++ )
						{
							if( binaryGrid[ offset + j*preLineCount ] )
							{
								squaredEDTPtr[j*preLineCount] = 0;
								dist = 0;
								first = false;
							}
							else if( !first )
							{
								dist++;
								squaredEDTPtr[j*preLineCount] = dist*dist;
							}
						}
					}

					// backward scan
					{
						bool first = true;
						unsigned int dist = 0;
						for( int j=(int)lineSize-1 ; j>=0 ; j-- )
						{
							if( binaryGrid[ offset + j*preLineCount ] )
							{
								squaredEDTPtr[j*preLineCount] = 0;
								dist = 0;
								first = false;
							}
							else if( !first )
							{
								dist++;
								unsigned int square = dist*dist;
								if( square<squaredEDTPtr[j*preLineCount] ) squaredEDTPtr[j*preLineCount] = square;
							}
						}
					}
				}
			);
	}
	if( verbose ) std::cout << pMeter( "Initial propagation" ) << std::endl;

	for( unsigned int dim=1 ; dim<Dim ; dim++ )
	{
		unsigned int lineSize = binaryGrid.res(dim);
		size_t preLineCount=1 , postLineCount=1;
		for( unsigned int d=0 ; d<dim ; d++ ) preLineCount *= binaryGrid.res(d);
		for( unsigned int d=dim+1 ; d<Dim ; d++ ) postLineCount *= binaryGrid.res(d);
		size_t lineCount = preLineCount * postLineCount;

		std::vector< Pointer( unsigned int ) > oldBuffer( threads ) , newBuffer( threads );
		for( unsigned int i=0 ; i<threads ; i++ ) oldBuffer[i] = NewPointer< unsigned int >( lineSize ) , newBuffer[i] = NewPointer< unsigned int >( lineSize );

		// In parallel process the different lines
		ThreadPool::ParallelFor
			(
				0 , lineCount ,
				[&]( unsigned int t , size_t i )
				{
					size_t pre = i % preLineCount , post = i / preLineCount;
					Pointer( unsigned int ) squaredEDTPtr = squaredEDT() + pre + post * ( preLineCount*lineSize );

					Pointer( unsigned int ) _oldBuffer = oldBuffer[t];
					Pointer( unsigned int ) _newBuffer = newBuffer[t];

					// Forward scan
					{
						unsigned int s=0;
						for( int j=0 ; j<(int)lineSize ; j++ )
						{
							_oldBuffer[j] = squaredEDTPtr[j*preLineCount];
							unsigned int dist = _oldBuffer[j];
							bool foundCloser = false;
							if( dist )
							{
								for( int t=s ; t<=j ; t++ )
								{
									unsigned int new_dist = _oldBuffer[t] + (j - t) * (j - t);
									if( new_dist<=dist ) dist = new_dist , s=t , foundCloser = true;
								}
							}
							if( !foundCloser ) s = j;
							_newBuffer[j] = dist;
						}
					}
					// Backward scan
					{
						unsigned int s = lineSize-1;
						for( int j=(int)lineSize-1 ; j>=0 ; j-- )
						{
							unsigned int dist = _newBuffer[j];
							bool foundCloser = false;
							if( dist )
							{
								for( int t=s ; t>j ; t-- )
								{
									unsigned int new_dist = _oldBuffer[t] + (j - t) * (j - t);
									if( new_dist<=dist ) dist = new_dist , s=t , foundCloser = true;
								}
								squaredEDTPtr[j*preLineCount] = dist;
							}
							if( !foundCloser ) s = j;
						}
					}
				}
			);
		for( unsigned int i=0 ; i<threads ; i++ ){ DeletePointer( oldBuffer[i] ) ; DeletePointer( newBuffer[i] ); }
	}
	if( verbose ) std::cout << pMeter( "Remaining propagation" ) << std::endl;

	return squaredEDT;
}

template< typename Real , unsigned int Dim >
template< typename BinaryType >
RegularGrid< Dim , std::pair< unsigned int , size_t > > SquaredEDT< Real , Dim >::FullSaito( const RegularGrid< Dim , BinaryType > &binaryGrid , bool verbose )
{
	Miscellany::PerformanceMeter pMeter( '.' );

	// Allocate the squaredEDT
	RegularGrid< Dim , std::pair< unsigned int , size_t > > squaredEDT;
	unsigned int threads = ThreadPool::NumThreads();
	{
		unsigned int res[Dim];
		for( int d=0 ; d<Dim ; d++ ) res[d] = binaryGrid.res(d);
		squaredEDT.resize( res );
		unsigned int max = 0;
		for( int d=0 ; d<Dim ; d++ ) max += res[d] * res[d];
		ThreadPool::ParallelFor( 0 , squaredEDT.resolution() , [&]( unsigned int , size_t i ){ squaredEDT[i].first = max; } );
	}

	// Scan along the first axis
	pMeter.reset();
	{
		unsigned int dim = 0;
		unsigned int lineSize = binaryGrid.res(dim);
		size_t preLineCount=1 , postLineCount=1;
		for( unsigned int d=0 ; d<dim ; d++ ) preLineCount *= binaryGrid.res(d);
		for( unsigned int d=dim+1 ; d<Dim ; d++ ) postLineCount *= binaryGrid.res(d);
		size_t lineCount = preLineCount * postLineCount;

		// In parallel process the different lines
		ThreadPool::ParallelFor
			(
				0 , lineCount ,
				[&]( unsigned int , size_t i )
				{
					size_t pre = i % preLineCount , post = i / preLineCount;
					size_t offset = pre + post * ( preLineCount*lineSize );
					Pointer( std::pair< unsigned int , size_t > ) squaredEDTPtr = squaredEDT() + pre + post * ( preLineCount*lineSize );

					// Forward scan
					{
						bool first = true;
						unsigned int dist = 0;
						size_t idx;
						for( int j=0 ; j<(int)lineSize ; j++ )
						{
							if( binaryGrid[ offset + j*preLineCount ] )
							{
								idx = offset + j*preLineCount;
								squaredEDTPtr[j*preLineCount].first = 0;
								squaredEDTPtr[j*preLineCount].second = idx;
								dist = 0;
								first = false;
							}
							else if( !first )
							{
								dist++;
								squaredEDTPtr[j*preLineCount].first = dist*dist;
								squaredEDTPtr[j*preLineCount].second = idx;
							}
						}
					}

					// backward scan
					{
						bool first = true;
						unsigned int dist = 0;
						size_t idx;
						for( int j=(int)lineSize-1 ; j>=0 ; j-- )
						{
							if( binaryGrid[ offset + j*preLineCount ] )
							{
								idx = offset + j*preLineCount;
								squaredEDTPtr[j*preLineCount].first = 0;
								squaredEDTPtr[j*preLineCount].second = idx;
								dist = 0;
								first = false;
							}
							else if( !first )
							{
								dist++;
								unsigned int square = dist*dist;
								if( square<squaredEDTPtr[j*preLineCount].first ) squaredEDTPtr[j*preLineCount].first = square , squaredEDTPtr[j*preLineCount].second = idx;
							}
						}
					}
				}
			);
	}
	if( verbose ) std::cout << pMeter( "Initial propagation" ) << std::endl;

	for( unsigned int dim=1 ; dim<Dim ; dim++ )
	{
		unsigned int lineSize = binaryGrid.res(dim);
		size_t preLineCount=1 , postLineCount=1;
		for( unsigned int d=0 ; d<dim ; d++ ) preLineCount *= binaryGrid.res(d);
		for( unsigned int d=dim+1 ; d<Dim ; d++ ) postLineCount *= binaryGrid.res(d);
		size_t lineCount = preLineCount * postLineCount;

		std::vector< Pointer( std::pair< unsigned int , size_t > ) > oldBuffer( threads ) , newBuffer( threads );
		for( unsigned int i=0 ; i<threads ; i++ ) oldBuffer[i] = NewPointer< std::pair< unsigned int , size_t > >( lineSize ) , newBuffer[i] = NewPointer< std::pair< unsigned int , size_t > >( lineSize );

		// In parallel process the different lines
		ThreadPool::ParallelFor
			( 
				0 , lineCount ,
				[&]( unsigned int t , size_t i )
				{
					size_t pre = i % preLineCount , post = i / preLineCount;
					Pointer( std::pair< unsigned int , size_t > ) squaredEDTPtr = squaredEDT() + pre + post * ( preLineCount*lineSize );

					Pointer( std::pair< unsigned int , size_t > ) _oldBuffer = oldBuffer[t];
					Pointer( std::pair< unsigned int , size_t > ) _newBuffer = newBuffer[t];

					// Forward scan
					{
						unsigned int s=0;
						for( int j=0 ; j<(int)lineSize ; j++ )
						{
							_oldBuffer[j] = squaredEDTPtr[j*preLineCount];
							unsigned int dist = _oldBuffer[j].first;
							size_t idx = _oldBuffer[j].second;
							bool foundCloser = false;
							if( dist )
							{
								for( int t=s ; t<=j ; t++ )
								{
									unsigned int new_dist = _oldBuffer[t].first + (j - t) * (j - t);
									if( new_dist<=dist ) dist = new_dist , s=t , foundCloser = true , idx = _oldBuffer[t].second;
								}
							}
							if( !foundCloser ) s = j;
							_newBuffer[j].first = dist;
							_newBuffer[j].second = idx;
						}
					}
					// Backward scan
					{
						unsigned int s = lineSize-1;
						for( int j=(int)lineSize-1 ; j>=0 ; j-- )
						{
							unsigned int dist = _newBuffer[j].first;
							size_t idx = _newBuffer[j].second;
							bool foundCloser = false;
							if( dist )
							{
								for( int t=s ; t>j ; t-- )
								{
									unsigned int new_dist = _oldBuffer[t].first + (j - t) * (j - t);
									if( new_dist<=dist ) dist = new_dist , s=t , foundCloser = true , idx = _oldBuffer[t].second;
								}
								squaredEDTPtr[j*preLineCount].first = dist;
								squaredEDTPtr[j*preLineCount].second = idx;
							}
							if( !foundCloser ) s = j;
						}
					}
				}
			);

		for( unsigned int i=0 ; i<threads ; i++ ){ DeletePointer( oldBuffer[i] ) ; DeletePointer( newBuffer[i] ); }
	}
	if( verbose ) std::cout << pMeter( "Remaining propagation" ) << std::endl;

	return squaredEDT;
}

template< typename Real , unsigned int Dim >
template< typename IndexType , unsigned int K >
RegularGrid< Dim , unsigned int > SquaredEDT< Real , Dim >::Saito( const SimplicialComplex< Real , Dim , K > &simplicialComplex , unsigned int depth , unsigned int lockDepth , XForm< Real , Dim+1 > &gridToModel , Real bBoxScale , bool verbose )
{
	Miscellany::PerformanceMeter pMeter( '.' );

	// Set the transformation from the unit cube to the grid;
	XForm< Real , Dim+1 > unitCubeToGrid;
	{
		Real res = (Real)( 1<< depth );
		Point< Real , Dim > offset;
		for( unsigned int d=0 ; d<Dim ; d++ ) offset[d] = (Real)-0.5;
		unitCubeToGrid = XForm< Real , Dim+1 >( SquareMatrix< Real , Dim >::Identity() * (Real)(1<<depth) , offset );
	}

	// Rasterize the simplices into a regular grid
	pMeter.reset();
	XForm< Real , Dim+1 > unitCubeToModel;
	typename Rasterizer< Real , Dim >::template SimplexRasterizationGrid< IndexType , K > sRaster = Rasterizer< Real , Dim >::template Rasterize< IndexType >( simplicialComplex , depth , lockDepth , unitCubeToModel , bBoxScale );
	gridToModel = unitCubeToModel * unitCubeToGrid.inverse();
	if( verbose ) std::cout << pMeter( "Rasterized" ) << std::endl;

	// Mark all voxels containing geometry
	RegularGrid< Dim , bool > raster;
	{
		unsigned int res[Dim];
		for( int d=0 ; d<Dim ; d++ ) res[d] = sRaster.res(d);
		raster.resize( res );
		ThreadPool::ParallelFor( 0 , raster.resolution() , [&]( unsigned int , size_t i ){ raster[i] = sRaster[i].size()!=0; } );
	}
	if( verbose ) std::cout << pMeter( "Marked geometry voxels" ) << std::endl;

#if 1
	return Saito( raster , verbose );
#else
	// Compute the squaredEDT
	RegularGrid< Dim , unsigned int > squaredEDT;
	unsigned int threads = ThreadPool::NumThreads();
	{
		unsigned int res[Dim];
		for( int d=0 ; d<Dim ; d++ ) res[d] = raster.res(d);
		squaredEDT.resize( res );
		unsigned int max = 0;
		for( int d=0 ; d<Dim ; d++ ) max += res[d] * res[d];
		ThreadPool::ParallelFor( 0 , squaredEDT.resolution() , [&]( unsigned int , size_t i ){ squaredEDT[i] = max; } );
	}

	// Scan along the first axis
	pMeter.reset();
	{
		unsigned int dim = 0;
		unsigned int lineSize = raster.res(dim);
		size_t preLineCount=1 , postLineCount=1;
		for( unsigned int d=0 ; d<dim ; d++ ) preLineCount *= raster.res(d);
		for( unsigned int d=dim+1 ; d<Dim ; d++ ) postLineCount *= raster.res(d);
		size_t lineCount = preLineCount * postLineCount;

		// In parallel process the different lines
		ThreadPool::ParallelFor
			(
				0 , lineCount ,
				[&]( unsigned int , size_t i )
				{
					size_t pre = i % preLineCount , post = i / preLineCount;
					size_t offset = pre + post * ( preLineCount*lineSize );
					Pointer( unsigned int ) squaredEDTPtr = squaredEDT() + pre + post * ( preLineCount*lineSize );

					// Forward scan
					{
						bool first = true;
						unsigned int dist = 0;
						for( int j=0 ; j<(int)lineSize ; j++ )
						{
							if( raster[ offset + j*preLineCount ] )
							{
								squaredEDTPtr[j*preLineCount] = 0;
								dist = 0;
								first = false;
							}
							else if( !first )
							{
								dist++;
								squaredEDTPtr[j*preLineCount] = dist*dist;
							}
						}
					}

					// backward scan
					{
						bool first = true;
						unsigned int dist = 0;
						for( int j=(int)lineSize-1 ; j>=0 ; j-- )
						{
							if( raster[ offset + j*preLineCount ] )
							{
								squaredEDTPtr[j*preLineCount] = 0;
								dist = 0;
								first = false;
							}
							else if( !first )
							{
								dist++;
								unsigned int square = dist*dist;
								if( square<squaredEDTPtr[j*preLineCount] ) squaredEDTPtr[j*preLineCount] = square;
							}
						}
					}
				}
			);
	}
	if( verbose ) std::cout << pMeter( "Initial propagation" ) << std::endl;

	for( unsigned int dim=1 ; dim<Dim ; dim++ )
	{
		unsigned int lineSize = raster.res(dim);
		size_t preLineCount=1 , postLineCount=1;
		for( unsigned int d=0 ; d<dim ; d++ ) preLineCount *= raster.res(d);
		for( unsigned int d=dim+1 ; d<Dim ; d++ ) postLineCount *= raster.res(d);
		size_t lineCount = preLineCount * postLineCount;

		std::vector< Pointer( unsigned int ) > oldBuffer( threads ) , newBuffer( threads );
		for( unsigned int i=0 ; i<threads ; i++ ) oldBuffer[i] = NewPointer< unsigned int >( lineSize ) , newBuffer[i] = NewPointer< unsigned int >( lineSize );

		// In parallel process the different lines
		ThreadPool::ParallelFor
			(
				0 , lineCount ,
				[&]( unsigned int t , size_t i )
				{
					size_t pre = i % preLineCount , post = i / preLineCount;
					Pointer( unsigned int ) squaredEDTPtr = squaredEDT() + pre + post * ( preLineCount*lineSize );

					Pointer( unsigned int ) _oldBuffer = oldBuffer[t];
					Pointer( unsigned int ) _newBuffer = newBuffer[t];

					// Forward scan
					{
						unsigned int s=0;
						for( int j=0 ; j<(int)lineSize ; j++ )
						{
							_oldBuffer[j] = squaredEDTPtr[j*preLineCount];
							unsigned int dist = _oldBuffer[j];
							bool foundCloser = false;
							if( dist )
							{
								for( int t=s ; t<=j ; t++ )
								{
									unsigned int new_dist = _oldBuffer[t] + (j - t) * (j - t);
									if( new_dist<=dist ) dist = new_dist , s=t , foundCloser = true;
								}
							}
							if( !foundCloser ) s = j;
							_newBuffer[j] = dist;
						}
					}
					// Backward scan
					{
						unsigned int s = lineSize-1;
						for( int j=(int)lineSize-1 ; j>=0 ; j-- )
						{
							unsigned int dist = _newBuffer[j];
							bool foundCloser = false;
							if( dist )
							{
								for( int t=s ; t>j ; t-- )
								{
									unsigned int new_dist = _oldBuffer[t] + (j - t) * (j - t);
									if( new_dist<=dist ) dist = new_dist , s=t , foundCloser = true;
								}
								squaredEDTPtr[j*preLineCount] = dist;
							}
							if( !foundCloser ) s = j;
						}
					}
				}
			);
		for( unsigned int i=0 ; i<threads ; i++ ){ DeletePointer( oldBuffer[i] ) ; DeletePointer( newBuffer[i] ); }
	}
	if( verbose ) std::cout << pMeter( "Remaining propagation" ) << std::endl;

	return squaredEDT;
#endif
}

////////////////
// Danielsson //
////////////////
template< typename Real , unsigned int Dim >
template< typename IndexType , unsigned int K >
RegularGrid< Dim , Real > SquaredEDT< Real , Dim >::Danielsson( const SimplicialComplex< Real , Dim , K > &simplicialComplex , unsigned int depth , unsigned int lockDepth , unsigned int radius , XForm< Real , Dim+1 > &gridToModel , Real bBoxScale , bool verbose )
{
	Miscellany::PerformanceMeter pMeter( '.' );

	// Set the transformation from the unit cube to the grid;
	XForm< Real , Dim+1 > unitCubeToGrid;
	{
		Real res = (Real)( 1<< depth );
		Point< Real , Dim > offset;
		for( unsigned int d=0 ; d<Dim ; d++ ) offset[d] = (Real)-0.5;
		unitCubeToGrid = XForm< Real , Dim+1 >( SquareMatrix< Real , Dim >::Identity() * (Real)(1<<depth) , offset );
	}

	pMeter.reset();
	// Rasterize the simplices into a regular grid
	XForm< Real , Dim+1 > unitCubeToModel;
	typename Rasterizer< Real , Dim >::template SimplexRasterizationGrid< IndexType , K > sRaster = Rasterizer< Real , Dim >::template Rasterize< IndexType >( simplicialComplex , depth , lockDepth , unitCubeToModel , bBoxScale );
	gridToModel = unitCubeToModel * unitCubeToGrid.inverse();
	if( verbose ) std::cout << pMeter( "Rasterized" ) << std::endl;

	// Construct the nearest keys
	std::vector< typename Simplex< Real , Dim , K >::NearestKey > nearestKeys( simplicialComplex.size() );
	{
		TransformedSimplicialComplex< Real , Dim , K > tSimplicialComplex( simplicialComplex , unitCubeToModel.inverse() );
		ThreadPool::ParallelFor( 0 , nearestKeys.size() , [&]( unsigned int , size_t i ){ nearestKeys[i].init( tSimplicialComplex[i] ); } );
	}
	if( verbose ) std::cout << pMeter( "Set keys" ) << std::endl;

	// Compute neighbor data info
	struct NeighborData
	{
		int offset[Dim];
		Real minDistance2;

		NeighborData( void ) : minDistance2(0) { memset( offset , 0 , sizeof(offset) ); }
		void init( const int off[Dim] , unsigned int  depth )
		{
			Real width = (Real)( 1./(1<<depth) );
			memcpy( offset , off , sizeof(int)*Dim );
			minDistance2 = 0;
			for( unsigned int d=0 ; d<Dim ; d++ ) if( off[d] ) minDistance2 += (Real)( ( -0.5 + fabs( off[d] ) ) * ( -0.5 + fabs( off[d] ) ) );
			minDistance2 /= width*width;
		};
		static bool Compare( const NeighborData &n1 , const NeighborData &n2 ){ return n1.minDistance2 < n2.minDistance2; }
	};
	std::vector< NeighborData > neighborData;
	{
		unsigned int neighborCount = 1;
		for( unsigned int d=0 ; d<Dim ; d++ ) neighborCount *= (2*radius+1);
		neighborData.resize( neighborCount );
		for( unsigned int n=0 ; n<neighborData.size() ; n++ )
		{
			int idx[Dim];
			unsigned int _n = n;
			for( unsigned int d=0 ; d<Dim ; d++ ){ idx[d] = (int)( _n % (2*radius+1) ) - radius ; _n /= (2*radius+1); }
			neighborData[n].init( idx , depth );
		}
		std::sort( neighborData.begin() , neighborData.end() , NeighborData::Compare );
	}

	// For every voxel, compute the index of and distance to the nearest geometry within the specified radius
	unsigned int res[Dim];
	RegularGrid< Dim , std::pair< Real , IndexType > > nearest;
	{
		for( int d=0 ; d<Dim ; d++ ) res[d] = sRaster.res(d);
		nearest.resize( res );

		// For each voxel, compute the index of, and distance to, the nearest simplex with the prescribed radius
		ThreadPool::ParallelFor
			(
				0 , nearest.resolution() ,
				[&]( unsigned int , size_t i )
				{
					// Get the index and center of the voxel
					unsigned int idx[Dim];
					Point< Real , Dim > c;
					{
						long long ii = i;
						for( unsigned int d=0 ; d<Dim ; d++ )
						{
							idx[d] = ii % nearest.res(d);
							ii /= nearest.res(d);
							c[d] = (Real)( idx[d]+0.5 ) / nearest.res( d );
						}
					}

					// Iterate over all neighbors and compute the distance from the center to the neighbors' geometry
					IndexType index = -1;
					Real l2 = std::numeric_limits< Real >::infinity();
					for( int i=0 ; i<neighborData.size() ; i++ )
					{
						if( neighborData[i].minDistance2>l2 ) break;
						int _idx[Dim];
						bool inBounds = true;
						for( unsigned int d=0 ; d<Dim ; d++ )
						{
							_idx[d] = idx[d] + neighborData[i].offset[d];
							if( _idx[d]<0 || _idx[d]>=(int)nearest.res(d) ) inBounds = false;
						}
						if( inBounds )
						{
							const std::vector< std::pair< IndexType , Simplex< Real , Dim , K > > > &indexedSimplices = sRaster( _idx );
							for( size_t j=0 ; j<indexedSimplices.size() ; j++ )
							{
								Real _l2 = Point< Real , Dim >::SquareNorm( c - nearestKeys[ indexedSimplices[j].first ].nearest( c ) );
								if( _l2<l2 ) l2 = _l2 , index = indexedSimplices[j].first;
							}
						}
					}
					nearest[i] = std::pair< Real , IndexType >( l2 , index );
				}
			);
	}
	if( verbose ) std::cout << pMeter( "Computed initial distances" ) << std::endl;

	_Danielsson< IndexType , K , Dim , true >( nearestKeys , nearest() , res , Point< Real , Dim >() );
	if( verbose ) std::cout << pMeter( "Propagated distances" ) << std::endl;

	// Transform the nearest geometry information into distances
	RegularGrid< Dim , Real > edt;
	edt.resize( res );
	ThreadPool::ParallelFor( 0 , edt.resolution() , [&]( unsigned int , size_t i ){ edt[i] = nearest[i].first; } );
	return edt;
}

template< typename Real , unsigned int Dim >
template< typename IndexType , unsigned int K , unsigned int SliceDim , bool MultiThreaded >
typename std::enable_if< SliceDim!=0 >::type SquaredEDT< Real , Dim >::_Danielsson( const std::vector< typename Simplex< Real , Dim , K >::NearestKey > &nearestKeys , Pointer( std::pair< Real , IndexType > ) sliceNearest , const unsigned int res[SliceDim] , Point< Real , Dim > center )
{
	size_t sliceSize = 1;
	for( unsigned int d=0 ; d<SliceDim-1 ; d++ ) sliceSize *= res[d];

	auto UpdateCenter = []( long long i , const unsigned int res[] , Point< Real , Dim > &center ){ for( unsigned int d=0 ; d<SliceDim-1 ; d++ ){ center[d] = (Real)( (i%res[d]) + 0.5 ) / res[d] ; i /= res[d]; } };
	auto Update = [&]( long long i , int j , int offset )
	{
		if( sliceNearest[ sliceSize*(j+offset) + i ].second!=-1 )
		{
			Point< Real , Dim > c = center;
			UpdateCenter( i , res , c );
			Point< Real , Dim > n = nearestKeys[ sliceNearest[ sliceSize*(j+offset) + i ].second ].nearest( c );
			Real l2 = Point< Real , Dim >::SquareNorm( n - c );
			if( l2<sliceNearest[ sliceSize*j + i ].first ) sliceNearest[ sliceSize*j + i ] = std::pair< Real , IndexType >( l2 , sliceNearest[ sliceSize*(j+offset) + i ].second );
		}
	};

	for( int j=1 ; j<(int)res[SliceDim-1] ; j++ )
	{
		center[ SliceDim-1 ] = (Real)( j + 0.5 ) / res[ SliceDim-1 ];
		// Update from the previous slice
		if( MultiThreaded ) ThreadPool::ParallelFor( 0 , sliceSize , [&]( unsigned int , size_t i ){ Update( i , j , -1 ); } );
		else for( long long i=0 ; i<(long long)sliceSize ; i++ ) Update( i , j , -1 );
		// Update the current slice
		_Danielsson< IndexType , K , SliceDim-1 , false >( nearestKeys , sliceNearest + sliceSize*j , res , center );
	}
	for( int j=(int)res[SliceDim-1]-2 ; j>=0 ; j-- )
	{
		center[ SliceDim-1 ] = (Real)( j + 0.5 ) / res[ SliceDim-1 ];
		// Update from the next slice
		if( MultiThreaded ) ThreadPool::ParallelFor( 0 , sliceSize , [&]( unsigned int , size_t i ){ Update( i , j , 1 ); } );
		else
			for( long long i=0 ; i<(long long)sliceSize ; i++ )  Update( i , j , 1 );
		// Update the current slice
		_Danielsson< IndexType , K , SliceDim-1 , false >( nearestKeys , sliceNearest + sliceSize*j , res , center );
	}
}
