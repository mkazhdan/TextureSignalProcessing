/*
Copyright (c) 2025, Michael Kazhdan
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

///////////////////////
// Texels::TexelInfo //
///////////////////////

template< bool NodeAtCellCenter , typename Index , unsigned int Dim >
template< unsigned int K >
template< unsigned int EmbeddingDim , typename SimplexEmbeddingFunctor /* = std::function< Simplex< double , EmbeddingDim , K > ( size_t ) > */ >
Point< double , EmbeddingDim >
Texels< NodeAtCellCenter , Index , Dim >::TexelInfo< K >::position
(
	const SimplexEmbeddingFunctor & SEF
)
const
{
	static_assert( _IsValidSimplexFunctor< SimplexEmbeddingFunctor , EmbeddingDim , K >() , "[ERROR] SimplexEmbeddingFunctor poorly formed" );
	return Point< double , EmbeddingDim >( SEF( sIdx )( bc ) );
}

////////////
// Texels //
////////////


//////////////////////
// Helper functions //
//////////////////////

template< bool NodeAtCellCenter , typename Index , unsigned int Dim >
template< typename Real , unsigned int EmbeddingDim , unsigned int K , typename SimplexEmbeddingFunctor /* = std::function< Simplex< double , EmbeddingDim , K > > ( size_t ) > */ >
RegularGrid< Dim , Point< Real , EmbeddingDim > >
Texels< NodeAtCellCenter , Index , Dim >::GetTexelPositions
(
	size_t simplexNum ,
	SimplexEmbeddingFunctor && SEF ,
	const RegularGrid< Dim , TexelInfo< K > > &texelInfo
)
{
	static_assert( _IsValidSimplexFunctor< SimplexEmbeddingFunctor , EmbeddingDim , K >() , "[ERROR] SimplexEmbeddingFunctor poorly formed" );

	RegularGrid< Dim , Point< Real , EmbeddingDim > > texturePositions( texelInfo.res() );

	Point< double , EmbeddingDim > badPosition;
	for( unsigned int d=0 ; d<EmbeddingDim ; d++ ) badPosition[d] = std::numeric_limits< double >::infinity();

	ThreadPool::ParallelFor
	(
		0 , texelInfo.resolution() ,
		[&]( unsigned int , size_t i )
		{
			if( texelInfo[i].sIdx==-1 ) texturePositions[i] = badPosition;
			else texturePositions[i] = Point< Real , EmbeddingDim >( texelInfo[i].template position< EmbeddingDim >( SEF ) );
		}
	);
	return texturePositions;
}

template< bool NodeAtCellCenter , typename Index , unsigned int Dim >
Point< double , Dim >
Texels< NodeAtCellCenter , Index , Dim >::NodePosition
(
	typename RegularGrid< Dim >::Index I
)
{
	Point< double , Dim > p;
	if constexpr( NodeAtCellCenter ) for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = I[d] + 0.5;
	else                             for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = I[d] + 0.0;
	return p;
}

template< bool NodeAtCellCenter , typename Index , unsigned int Dim >
template< unsigned int K >
Simplex< double , Dim , K >
Texels< NodeAtCellCenter , Index , Dim >::TexelSimplex
(
	Simplex< double , Dim , K > simplex ,
	const unsigned int res[/*Dim*/]
)
{
	for( unsigned int k=0 ; k<=K ; k++ )
		if constexpr( NodeAtCellCenter ) for( unsigned int d=0 ; d<Dim ; d++ ) simplex[k][d] *= res[d];
		else                             for( unsigned int d=0 ; d<Dim ; d++ ) simplex[k][d] *= res[d]-1;
	return simplex;
}

//////////////////////////////
// Get simplices in support //
//////////////////////////////

template< bool NodeAtCellCenter , typename Index , unsigned int Dim >
template< bool Nearest , unsigned int K , typename SimplexFunctor /* = std::function< Simplex< double , Dim , K > ( size_t ) > */ >
RegularGrid< Dim , std::vector< Index > >
Texels< NodeAtCellCenter , Index , Dim >::GetSupportedSimplexIndices
(
	size_t simplexNum ,
	SimplexFunctor && SF ,
	const unsigned int res[/*Dim*/]
)
{
	static_assert( _IsValidSimplexFunctor< SimplexFunctor , Dim , K >() , "[ERROR] SimplexFunctor poorly formed" );

	typename RegularGrid< Dim >::Range range;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.second[d] = res[d];

	// Get the list of texel-index/simplex-index pairs
	std::vector< std::vector< std::pair< typename RegularGrid< Dim >::Index , Index > > > _tIndices( ThreadPool::NumThreads() );
	ThreadPool::ParallelFor
	(
		0 , simplexNum ,
		[&]( unsigned int t , size_t sIdx )
		{
			Rasterizer2D::RasterizeSupports< Nearest , NodeAtCellCenter >
				(
					TexelSimplex( SF(sIdx) , res ) ,
					[&]( typename RegularGrid< Dim >::Index I ){ _tIndices[t].push_back( std::make_pair( I , static_cast< Index >(sIdx) ) ); } ,
					range
				);
		}
	);

	// Add the simplex indices to the associated texels in the grid
	RegularGrid< Dim , std::vector< Index > > tIndices( res );
	for( unsigned int i=0 ; i<_tIndices.size() ; i++ ) for( unsigned int j=0 ; j<_tIndices[i].size() ; j++ ) tIndices( _tIndices[i][j].first ).push_back( _tIndices[i][j].second );

	// Remove redundant indices
	ThreadPool::ParallelFor
	(
		0 , tIndices.size() ,
		[&]( size_t i )
		{
			std::set< Index > _tIndices;
			for( unsigned int j=0 ; j<tIndices[i].size() ; j++ ) _tIndices.insert( tIndices[i][j] );
			if( _tIndices.size()!=tIndices[i].size() )
			{
				tIndices[i].resize(0);
				tIndices[i].reserve( _tIndices.size() );
				for( auto iter=_tIndices.begin() ; iter!=_tIndices.end() ; iter++ ) tIndices[i].push_back( *iter );
			}
		}
	);

	return tIndices;
}

template< bool NodeAtCellCenter , typename Index , unsigned int Dim >
template< bool Nearest , unsigned int K , typename SimplexFunctor >
RegularGrid< Dim , typename Texels< NodeAtCellCenter , Index , Dim >::template TexelInfo< K > >
Texels< NodeAtCellCenter , Index , Dim >::GetSupportedTexelInfo
(
	size_t simplexNum ,
	SimplexFunctor && SF ,
	const unsigned int res[/*Dim*/]
)
{
	static_assert( _IsValidSimplexFunctor< SimplexFunctor , Dim , K >() , "[ERROR] SimplexFunctor poorly formed" );

	typename RegularGrid< Dim >::Range range;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.second[d] = res[d];

	std::vector< std::vector< std::pair< typename RegularGrid< Dim >::Index , TexelInfo< K > > > > _texelInfo( ThreadPool::NumThreads() );
	ThreadPool::ParallelFor
	(
		0 , simplexNum ,
		[&]( unsigned int t , size_t sIdx )
		{
			Simplex< double , Dim , K > s = TexelSimplex( SF(sIdx) , res );
			auto RasterizationFunctor = [&]( typename RegularGrid< Dim >::Index I )
				{
					TexelInfo< K > ti;
					ti.sIdx = static_cast< Index >(sIdx);
					ti.bc = s.barycentricCoordinates( NodePosition( I ) );
					_texelInfo[t].push_back( std::make_pair( I , ti ) );
				};
			Rasterizer2D::RasterizeSupports< Nearest , NodeAtCellCenter >( s , RasterizationFunctor , range );
		}
	);

	RegularGrid< Dim , TexelInfo< K > > texelInfo( res );
	for( unsigned int i=0 ; i<_texelInfo.size() ; i++ ) for( unsigned int j=0 ; j<_texelInfo[i].size() ; j++ )
		texelInfo( _texelInfo[i][j].first ) = _texelInfo[i][j].second;

	return texelInfo;
}

////////////////////////////////////
// Get node-overlapping simplices //
////////////////////////////////////

template< bool NodeAtCellCenter , typename Index , unsigned int Dim >
template< typename SimplexFunctor /* = std::function< Simplex< double , Dim , Dim > ( size_t ) > */ >
RegularGrid< Dim , std::vector< Index > >
Texels< NodeAtCellCenter , Index , Dim >::GetNodeSimplexIndices
(
	size_t simplexNum ,
	SimplexFunctor && SF ,
	const unsigned int res[/*Dim*/]
)
{
	static_assert( _IsValidSimplexFunctor< SimplexFunctor , Dim , Dim >() , "[ERROR] SimplexFunctor is poorly formed" );

	typename RegularGrid< Dim >::Range range;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.second[d] = res[d];

	std::vector< std::vector< std::pair< typename RegularGrid< Dim >::Index , Index > > > _tIndices( ThreadPool::NumThreads() );
	ThreadPool::ParallelFor
	(
		0 , simplexNum ,
		[&]( unsigned int t , size_t sIdx )
		{
			Rasterizer2D::RasterizeNodes< NodeAtCellCenter >
				(
					TexelSimplex( SF(sIdx) , res ) ,
					[&]( typename RegularGrid< Dim >::Index I ){ _tIndices[t].push_back( std::make_pair( I , static_cast< Index >(sIdx) ) ); },
					range
				);
		}
	);

	RegularGrid< Dim , std::vector< Index > > tIndices( res );

	for( unsigned int i=0 ; i<_tIndices.size() ; i++ ) for( unsigned int j=0 ; j<_tIndices[i].size() ; j++ ) tIndices( _tIndices[i][j].first ).push_back( _tIndices[i][j].second );

	ThreadPool::ParallelFor
	(
		0 , tIndices.size() ,
		[&]( size_t i )
		{
			std::set< Index > _tIndices;
			for( unsigned int j=0 ; j<tIndices[i].size() ; j++ ) _tIndices.insert( tIndices[i][j] );
			if( _tIndices.size()!=tIndices[i].size() )
			{
				tIndices[i].resize(0);
				tIndices[i].reserve( _tIndices.size() );
				for( auto iter=_tIndices.begin() ; iter!=_tIndices.end() ; iter++ ) tIndices[i].push_back( *iter );
			}
		}
	);

	return tIndices;
}


template< bool NodeAtCellCenter , typename Index , unsigned int Dim >
template< typename SimplexFunctor /* = std::function< Simplex< double , Dim , Dim > ( size_t ) > */ >
RegularGrid< Dim , typename Texels< NodeAtCellCenter , Index , Dim >::template TexelInfo< Dim > >
Texels< NodeAtCellCenter , Index , Dim >::GetNodeTexelInfo
(
	size_t simplexNum ,
	SimplexFunctor && SF ,
	const unsigned int res[/*Dim*/]
)
{
	static_assert( _IsValidSimplexFunctor< SimplexFunctor , Dim , Dim >() , "[ERROR] SimplexFunctor is poorly formed" );

	typename RegularGrid< Dim >::Range range;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.second[d] = res[d];

	std::vector< std::vector< std::pair< typename RegularGrid< Dim >::Index , Index > > > _tIndices( ThreadPool::NumThreads() );
	ThreadPool::ParallelFor
	(
		0 , simplexNum ,
		[&]( unsigned int t , size_t sIdx )
		{
			Rasterizer2D::RasterizeNodes< NodeAtCellCenter >
				(
					TexelSimplex( SF(sIdx) , res ) ,
					[&]( typename RegularGrid< Dim >::Index I ){ _tIndices[t].push_back( std::make_pair( I , static_cast< Index >(sIdx) ) ); },
					range
				);
		}
	);

	RegularGrid< Dim , TexelInfo< Dim > > tInfo( res );

	ThreadPool::ParallelFor
		(
			0 , _tIndices.size() ,
			[&]( size_t i )
			{
				for( unsigned int j=0 ; j<_tIndices[i].size() ; j++ ) Atomic< Index >::Set( tInfo( _tIndices[i][j].first ).sIdx , _tIndices[i][j].second );
			}
		);

	ThreadPool::ParallelFor
		(
			0 , tInfo.size() ,
			[&]( size_t i )
			{
				if( tInfo[i].sIdx!=-1 )
				{
					typename RegularGrid< Dim >::Index I = RegularGrid< Dim >::FactorIndex( i , tInfo.res() );
					Point< double , Dim > p = NodePosition( I );
					Simplex< double , Dim , Dim > s = TexelSimplex( SF( tInfo[i].sIdx ) , tInfo.res() );
					tInfo[i].bc = s.barycentricCoordinates( p );	
				}
			}
		);

	return tInfo;
}

//////////////////////////////////////////////
// Dilate assignment of simplices to texels //
//////////////////////////////////////////////

template< bool NodeAtCellCenter , typename Index , unsigned int Dim >
size_t 
Texels< NodeAtCellCenter , Index , Dim >::DilateSimplexIndices
(
	RegularGrid< Dim , std::vector< Index > > &tIndices ,
	unsigned int dilationRadius ,
	bool verbose
)
{
	typename RegularGrid< Dim >::Range range;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.second[d] = tIndices.res(d);

	std::atomic< size_t > boundaryCount = 0;

	for( unsigned int r=0 ; r<dilationRadius ; r++ )
	{
		Miscellany::PerformanceMeter pMeter;

		RegularGrid< Dim , std::vector< Index > > _tIndices( tIndices.res() );

		// Get the new boundary simplex indices
		ThreadPool::ParallelFor
		(
			0 , tIndices.size() ,
			[&]( size_t i )
			{
				// If this is not an interior texel
				if( tIndices[i].size()==0 )
				{
					typename RegularGrid< Dim >::Index I = RegularGrid< Dim >::FactorIndex( i, tIndices.res() );

					std::set< Index > neighbors;

					bool isBoundary = false;
					RegularGrid< Dim >::Range::Intersect( typename RegularGrid< Dim >::Range(I).dilate(1) , range ).process( [&]( typename RegularGrid< Dim >::Index J ){ if( tIndices(J).size() ) isBoundary = true; } );
					if( isBoundary )
					{
						auto Kernel = [&]( typename RegularGrid< Dim >::Index J ){ for( unsigned int j=0 ; j<tIndices(J).size() ; j++ ) neighbors.insert( tIndices(J)[j] ); };
						RegularGrid< Dim >::Range::Intersect( typename RegularGrid< Dim >::Range(I).dilate(2) , range ).process( Kernel );
					}

					if( neighbors.size() )
					{
						boundaryCount++;
						_tIndices[i].reserve( neighbors.size() );
						for( auto iter=neighbors.begin() ; iter!=neighbors.end() ; iter++ ) _tIndices[i].push_back( *iter );
					}
				}
			}
		);

		// Copy the new boundary simplex indices in
		ThreadPool::ParallelFor( 0 , tIndices.size() , [&]( size_t i ){ if( _tIndices[i].size() ) std::swap( tIndices[i] , _tIndices[i] ); } );

		if( verbose ) std::cout << pMeter( "Radius: " + std::to_string( r+1 ) + " / " + std::to_string( dilationRadius ) ) << std::endl;
	}
	return boundaryCount;
}

template< bool NodeAtCellCenter , typename Index , unsigned int Dim >
template< unsigned int K , typename SimplexFunctor /* = std::function< Simplex< double , Dim , Dim > ( size_t ) > */ >
size_t
Texels< NodeAtCellCenter , Index , Dim >::DilateTexelInfo( SimplexFunctor && SF , RegularGrid< Dim , TexelInfo< K > > &tInfo , unsigned int dilationRadius , bool verbose )
{
	static_assert( std::is_convertible_v< SimplexFunctor , std::function< Simplex< double , Dim , K > ( size_t ) > > , "[ERROR] SimplexFunctor poorly formed" );

	typename RegularGrid< Dim >::Range range;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.second[d] = tInfo.res(d);

	std::atomic< size_t > boundaryCount = 0;
	for( unsigned int r=0 ; r<dilationRadius ; r++ )
	{
		Miscellany::PerformanceMeter pMeter;

		RegularGrid< Dim , std::vector< Index > > sIndices( tInfo.res() );

		// For the (exterior) boundary texels, set the set of active 2-ring neighbor simplex indices
		ThreadPool::ParallelFor
		(
			0 , tInfo.size() , 
			[&]( size_t i )
			{
				// If this is not an interior texel
				if( tInfo[i].sIdx==-1 )
				{
					typename RegularGrid< Dim >::Index I = RegularGrid< Dim >::FactorIndex( i , tInfo.res() );

					bool isBoundary = false;
					RegularGrid< Dim >::Range::Intersect( typename RegularGrid< Dim >::Range(I).dilate(1) , range ).process( [&]( typename RegularGrid< Dim >::Index J ){ if( tInfo(J).sIdx!=-1 ) isBoundary = true; } );
					if( isBoundary )
					{
						boundaryCount++;
						std::set< Index > nbrs;
						RegularGrid< Dim >::Range::Intersect( typename RegularGrid< Dim >::Range(I).dilate(2) , range ).process( [&]( typename RegularGrid< Dim >::Index J ){ if( tInfo(J).sIdx!=-1 ) nbrs.insert( tInfo(J).sIdx ); } );

						sIndices[i].reserve( nbrs.size() );
						for( auto iter=nbrs.begin() ; iter!=nbrs.end() ; iter++ ) sIndices[i].push_back( *iter );
					}
				}
			}
		);

		// For the (exterior) boundary texels, set the texel info from the active neighbor simplex indices set
		ThreadPool::ParallelFor
		(
			0 , tInfo.size() , 
			[&]( size_t i )
			{
				// If this is not an interior texel
				if( sIndices[i].size() )
				{
					if( tInfo[i].sIdx!=-1 ) MK_ERROR_OUT( "Exepcted exterior texel" );

					typename RegularGrid< Dim >::Index I = RegularGrid< Dim >::FactorIndex( i , tInfo.res() );
					Point< double , Dim > p = NodePosition( I );

					unsigned int sIdx =-1;
					double squreDist = std::numeric_limits< double >::infinity();
					for( auto iter=sIndices[i].begin() ; iter!=sIndices[i].end() ; iter++ )
					{
						Simplex< double , Dim , Dim > s = TexelSimplex( SF( *iter ) , tInfo.res() );
						double d2 = Point< double , Dim >::SquareDistance( p , s( s.nearestBC( p ) ) );
						if( d2<squreDist ) squreDist = d2 , sIdx = *iter;
					}
					if( sIdx==-1 ) MK_ERROR_OUT( "Could not find neighboring simplex" );

					Simplex< double , Dim , Dim > s = TexelSimplex( SF( sIdx ) , tInfo.res() );

					tInfo[i].sIdx = sIdx;
					tInfo[i].bc = s.barycentricCoordinates( p );	
				}
			}
		);

		if( verbose ) std::cout << pMeter( "Radius: " + std::to_string( r+1 ) + " / " + std::to_string( dilationRadius ) ) << std::endl;
	}

	return boundaryCount;
}

//////////////////////////////
// Compute nearest triangle //
//////////////////////////////

template< bool NodeAtCellCenter , typename Index , unsigned int Dim >
template< unsigned int K , typename SimplexFunctor /* = std::function< Simplex< double , Dim , Dim > ( size_t ) > */ >
RegularGrid< Dim , typename Texels< NodeAtCellCenter , Index , Dim >::template TexelInfo< K > >
Texels< NodeAtCellCenter , Index , Dim >::GetNearestTexelInfo
(
	SimplexFunctor && SF ,
	const RegularGrid< Dim , std::vector< Index > > & tIndices
)
{
	RegularGrid< Dim , TexelInfo< K > > texelInfo;
	texelInfo.resize( tIndices.res() );

	ThreadPool::ParallelFor
	(
		0 , texelInfo.size() ,
		[&]( size_t i )
		{
			if( tIndices[i].size()!=0 )
			{
				TexelInfo< K > tInfo;

				// Get the center of the texel
				typename RegularGrid< Dim >::Index I = RegularGrid< Dim >::FactorIndex( i , tIndices.res() );
				Point< double , Dim > p = NodePosition( I );

				// Finding the neighboring simplex closest to the texel center
				double dist2 = std::numeric_limits< double >::infinity();
				for( unsigned int j=0 ; j<tIndices[i].size() ; j++ )
				{
					Simplex< double , Dim , K > s = TexelSimplex( SF( tIndices[i][j] ) , tIndices.res() );
					double d2 = Point< double , Dim >::SquareDistance( p , s( s.nearestBC( p ) ) );
					if( d2<dist2 ) dist2 = d2 , tInfo.sIdx = tIndices[i][j];
				}

				Simplex< double , Dim , K > s = TexelSimplex( SF( tInfo.sIdx ) , tIndices.res() );
				tInfo.bc = s.barycentricCoordinates( p );	
				texelInfo[i] = tInfo;
			}
		}
	);

	return texelInfo;
}


template< bool NodeAtCellCenter , typename Index , unsigned int Dim >
template< bool Nearest , unsigned int K , bool UseVector , typename SimplexFunctor /* = std::function< Simplex< double , Dim , K > ( size_t ) > */ >
RegularGrid< Dim , typename Texels< NodeAtCellCenter , Index , Dim >::template TexelInfo< K > >
Texels< NodeAtCellCenter , Index , Dim >::GetSupportedTexelInfo
(
	size_t simplexNum ,
	SimplexFunctor && SF ,
	const unsigned int res[/*Dim*/] ,
	unsigned int dilationRadius ,
	bool verbose
)
{
	static_assert( _IsValidSimplexFunctor< SimplexFunctor , Dim , K >() , "[ERROR] SimplexEmbeddingFunctor poorly formed" );

	Miscellany::PerformanceMeter pMeter;

	RegularGrid< Dim , TexelInfo< K > > texelInfo;
	if constexpr( UseVector )
	{
		RegularGrid< Dim , std::vector< Index > > tIndices = GetSupportedSimplexIndices< Nearest , K >( simplexNum , SF , res );
		if( verbose ) std::cout << pMeter( "Supported simplex indices" ) << std::endl;
		if( dilationRadius )
		{
			DilateSimplexIndices( tIndices , dilationRadius , verbose );
			if( verbose ) std::cout << pMeter( "Dilated" ) << std::endl;
		}
		texelInfo = GetNearestTexelInfo< K >( SF , tIndices );
		if( verbose ) std::cout << pMeter( "Simplex indices to nearest" ) << std::endl;
	}
	else
	{
		texelInfo = GetSupportedTexelInfo< Nearest , K >( simplexNum , SF , res );
		if( verbose ) std::cout << pMeter( "Supported texel info" ) << std::endl;
		if( dilationRadius )
		{
			DilateTexelInfo( SF , texelInfo , dilationRadius , verbose );
			if( verbose ) std::cout << pMeter( "Dilated" ) << std::endl;
		}
	}

	return texelInfo;
}

template< bool NodeAtCellCenter , typename Index , unsigned int Dim >
template< bool UseVector , typename SimplexFunctor /* = std::function< Simplex< double , Dim , Dim > ( size_t ) > */ >
RegularGrid< Dim , typename Texels< NodeAtCellCenter , Index , Dim >::template TexelInfo< Dim > >
Texels< NodeAtCellCenter , Index , Dim >::GetNodeTexelInfo
(
	size_t simplexNum ,
	SimplexFunctor && SF ,
	const unsigned int res[/*Dim*/] ,
	unsigned int dilationRadius ,
	bool verbose
)
{
	static_assert( _IsValidSimplexFunctor< SimplexFunctor , Dim , Dim >() , "[ERROR] SimplexEmbeddingFunctor poorly formed" );

	Miscellany::PerformanceMeter pMeter;

	RegularGrid< Dim , TexelInfo< Dim > > texelInfo;
	if constexpr( UseVector )
	{
		RegularGrid< Dim , std::vector< Index > > tIndices = GetNodeSimplexIndices( simplexNum , SF , res );
		if( verbose ) std::cout << pMeter( "Node simplex indices" ) << std::endl;
		if( dilationRadius )
		{
			DilateSimplexIndices( tIndices , dilationRadius , verbose );
			if( verbose ) std::cout << pMeter( "Dilated" ) << std::endl;
		}
		texelInfo = GetNearestTexelInfo( SF , tIndices );
		if( verbose ) std::cout << pMeter( "Simplex indices to nearest" ) << std::endl;
	}
	else
	{
		texelInfo = GetNodeTexelInfo( simplexNum , SF , res );
		if( verbose ) std::cout << pMeter( "Node texel info" ) << std::endl;
		if( dilationRadius )
		{
			DilateTexelInfo( SF , texelInfo , dilationRadius , verbose );
			if( verbose ) std::cout << pMeter( "Dilated" ) << std::endl;
		}
	}

	return texelInfo;
}

template< bool NodeAtCellCenter , typename Index , unsigned int Dim >
template< unsigned int EmbeddingDim , bool Nearest , unsigned int K , bool UseVector , typename VertexEmbeddingFunctor /* = std::function< Point< double , EmbeddingDim > ( size_t ) > */ , typename SimplexIndexFunctor /* = std::function< SimplexIndex< Dim >( size_t ) > */ , typename SimplexFunctor /* = std::function< Simplex< double , Dim , K > ( size_t ) > */ >
RegularGrid< Dim , typename Texels< NodeAtCellCenter , Index , Dim >::template TexelInfo< K > >
Texels< NodeAtCellCenter , Index , Dim >::GetSupportedTexelInfo
(
	size_t simplexNum ,
	VertexEmbeddingFunctor && VF ,
	SimplexIndexFunctor && SIF ,
	SimplexFunctor && SF ,
	const unsigned int res[/*Dim*/] ,
	unsigned int dilationRadius ,
	bool verbose
)
{
	static_assert( _IsValidVertexFunctor< VertexEmbeddingFunctor , EmbeddingDim >() , "[ERROR] VertexEmbeddingFunctor poorly formed" );
	static_assert( _IsValidSimplexIndexFunctor< SimplexIndexFunctor , Dim >()       , "[ERROR] SimplexIndexFunctor poorly formed" );
	static_assert( _IsValidSimplexFunctor< SimplexFunctor , Dim , K >()             , "[ERROR] SimplexEmbeddingFunctor poorly formed" );

	RegularGrid< Dim , TexelInfo< K > > texelInfo = GetSupportedTexelInfo< Nearest , K >( simplexNum , std::forward< SimplexFunctor >( SF ) , res , dilationRadius , verbose );

	Miscellany::PerformanceMeter pMeter;
	FlowTexelInfoToInterior< EmbeddingDim >( simplexNum , std::forward< VertexEmbeddingFunctor >( VF ) , std::forward< SimplexIndexFunctor >( SIF ) , texelInfo );
	if( verbose ) std::cout << pMeter( "Flowed" ) << std::endl;

	return texelInfo;
}

template< bool NodeAtCellCenter , typename Index , unsigned int Dim >
template< unsigned int EmbeddingDim , bool UseVector , typename VertexEmbeddingFunctor /* = std::function< Point< double , EmbeddingDim > ( size_t ) > */ , typename SimplexIndexFunctor /* = std::function< SimplexIndex< Dim >( size_t ) > */ , typename SimplexFunctor /* = std::function< Simplex< double , Dim , Dim > ( size_t ) > */ >
RegularGrid< Dim , typename Texels< NodeAtCellCenter , Index , Dim >::template TexelInfo< Dim > >
Texels< NodeAtCellCenter , Index , Dim >::GetNodeTexelInfo
(
	size_t simplexNum ,
	VertexEmbeddingFunctor && VF ,
	SimplexIndexFunctor && SIF ,
	SimplexFunctor && SF ,
	const unsigned int res[/*Dim*/] ,
	unsigned int dilationRadius ,
	bool verbose
)
{
	static_assert( _IsValidVertexFunctor< VertexEmbeddingFunctor , EmbeddingDim >() , "[ERROR] VertexEmbeddingFunctor poorly formed" );
	static_assert( _IsValidSimplexIndexFunctor< SimplexIndexFunctor , Dim >()       , "[ERROR] SimplexIndexFunctor poorly formed" );
	static_assert( _IsValidSimplexFunctor< SimplexFunctor , Dim , Dim >()           , "[ERROR] SimplexEmbeddingFunctor poorly formed" );

	RegularGrid< Dim , TexelInfo< Dim > > texelInfo = GetNodeTexelInfo( simplexNum , std::forward< SimplexFunctor >( SF ) , res , dilationRadius , verbose );

	Miscellany::PerformanceMeter pMeter;
	FlowTexelInfoToInterior< EmbeddingDim >( simplexNum , std::forward< VertexEmbeddingFunctor >( VF ) , std::forward< SimplexIndexFunctor >( SIF ) , texelInfo );
	if( verbose ) std::cout << pMeter( "Flowed" ) << std::endl;

	return texelInfo;
}

template< bool NodeAtCellCenter , typename Index , unsigned int Dim >
template< unsigned int EmbeddingDim , typename VertexEmbeddingFunctor /* = std::function< Point< double , EmbeddingDim > ( size_t ) > */ , typename SimplexIndexFunctor /* = std::function< SimplexIndex< Dim >( size_t ) > */ >
void
Texels< NodeAtCellCenter , Index , Dim >::FlowTexelInfoToInterior
(
	size_t simplexNum ,
	VertexEmbeddingFunctor && VF ,
	SimplexIndexFunctor && SIF ,
	RegularGrid< Dim , TexelInfo< Dim > > &texelInfo ,
	double eps ,
	bool noWarning
)
{
	static_assert( _IsValidVertexFunctor< VertexEmbeddingFunctor , EmbeddingDim >() , "[ERROR] VertexEmbeddingFunctor poorly formed" );
	static_assert( _IsValidSimplexIndexFunctor< SimplexIndexFunctor , Dim >() , "[ERROR] SimplexIndexFunctor poorly formed" );

	std::vector< SimplexIndex< Dim > > simplices( simplexNum );
	for( size_t i=0 ; i<simplexNum ; i++ ) simplices[i] = SIF( i );
	FEM::RiemannianMesh< double , unsigned int > rMesh( GetPointer( simplices ) , simplexNum );
	rMesh.template setMetricFromEmbedding< EmbeddingDim >( [&]( unsigned int v ){ return VF(v); } );
	Pointer( FEM::CoordinateXForm< double > ) xForms = rMesh.getCoordinateXForms();

	ThreadPool::ParallelFor
	(
		0 , texelInfo.resolution() ,
		[&]( unsigned int , size_t i )
		{
			if( texelInfo[i].sIdx!=-1 )
			{
				bool isInterior = true;
				for( unsigned int k=0 ; k<=2 ; k++ ) isInterior &= texelInfo[i].bc[k]>=0;
				if( !isInterior )
				{
					// Get the representation of the texel center in the tangent frame of the closest simplex
					FEM::HermiteSamplePoint< double > h;
					h.tIdx = texelInfo[i].sIdx;
					for( unsigned int d=0 ; d<Dim ; d++ ) h.p[d] = 1./(Dim+1) , h.v[d] = texelInfo[i].bc[d+1] - 1./(Dim+1);

					// Walk along the tangent direction
					rMesh.exp( xForms , h , eps , noWarning );
					texelInfo[i].sIdx = h.tIdx;
					texelInfo[i].bc[0] = 1.;
					for( unsigned int d=0 ; d<Dim ; d++ ) texelInfo[i].bc[d+1] = h.p[d] , texelInfo[i].bc[0] -= h.p[d];
				}
			}
		}
	);
	DeletePointer( xForms );
}