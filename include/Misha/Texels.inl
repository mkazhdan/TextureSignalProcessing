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

///////////////
// TexelInfo //
///////////////


template< typename Real , unsigned int EmbeddingDim , typename SimplexEmbeddingFunctor /* = std::function< Simplex< double , EmbeddingDim , Dim > ( size_t ) > */ >
Point< Real , EmbeddingDim > TexelInfo::position( const SimplexEmbeddingFunctor & SEF ) const
{
	return Point< double , EmbeddingDim >( SEF( sIdx )( bc ) );
}


///////////////
///////////////

inline double DistanceToEdge( Point< double , Dim > p , Point< double , Dim > v0 , Point< double , Dim > v1 )
{
	// E(s) = || p - ( v0*(1-s) + v1*s ) ||^2
	// 0 = E'(s) = < p - ( v0*(1-s) + v1*s ) , v0 - v1 >
	// =>   < p - v0 , v0 - v1 > = - s * || v0 - v1 ||^2
	// =>   s = < p - v0  , v1-v0 > / || v1 - v0 ||^2
	return Point< double , Dim >::Dot( p-v0 , v1-v0 ) / Point< double , Dim >::SquareNorm( v1-v0 );
}

inline double DistanceToTriangle( Point< double , Dim > p , Simplex< double , Dim , Dim > s )
{
	Point< double , Dim+1 > bc = s.barycentricCoordinates( p );
	if( bc[0]>=0 && bc[1]>=0 && bc[2]>=0 ) return 0.;

	if     ( bc[0]<0 && bc[1]<0 ) return Point< double , Dim >::Length( p - s[0] );
	else if( bc[0]<0 && bc[2]<0 ) return Point< double , Dim >::Length( p - s[2] );
	else if( bc[1]<0 && bc[2]<0 ) return Point< double , Dim >::Length( p - s[1] );
	else if( bc[0]<0 ) return DistanceToEdge( p , s[0] , s[2] );
	else if( bc[1]<0 ) return DistanceToEdge( p , s[1] , s[0] );
	else if( bc[2]<0 ) return DistanceToEdge( p , s[2] , s[1] );
	else MK_ERROR_OUT( "At least one of the barycentric coordinates should be positive" );
	return 0.;
}

template< typename Real , unsigned int EmbeddingDim , typename SimplexEmbeddingFunctor /* = std::function< Simplex< double , EmbeddingDim , Dim > > ( size_t ) > */ >
RegularGrid< Dim , Point< Real , EmbeddingDim > > GetTexelPositions( size_t simplexNum , SimplexEmbeddingFunctor && SEF , const RegularGrid< Dim , TexelInfo > &texelInfo )
{
	static_assert( std::is_convertible_v< SimplexEmbeddingFunctor , std::function< Simplex< double , EmbeddingDim , Dim >( size_t ) > > , "[ERROR] SimplexEmbeddingFunctor poorly formed" );

	RegularGrid< Dim , Point< float , EmbeddingDim > > texturePositions( texelInfo.res() );

	Point< Real , EmbeddingDim > badPosition;
	for( unsigned int d=0 ; d<EmbeddingDim ; d++ ) badPosition[d] = std::numeric_limits< Real >::infinity();

	ThreadPool::ParallelFor
	(
		0 , texelInfo.resolution() ,
		[&]( unsigned int , size_t i )
		{
			if( texelInfo[i].sIdx==-1 ) texturePositions[i] = badPosition;
			else texturePositions[i] = texelInfo[i].template position< Real , EmbeddingDim >( SEF );
		}
	);
	return texturePositions;
}

template< bool NodeAtCellCenter >
Point< double , Dim > TexelNodePosition( typename RegularGrid< Dim >::Index I )
{
	Point< double , Dim > p;
	if constexpr( NodeAtCellCenter ) for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = I[d] + 0.5;
	else                             for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = I[d] + 0.0;
	return p;
}

template< bool NodeAtCellCenter , unsigned int K >
Simplex< double , Dim , K > GetTexelSpaceSimplex( Simplex< double , Dim , K > simplex , const unsigned int res[Dim] )
{
	for( unsigned int k=0 ; k<=K ; k++ )
		if constexpr( NodeAtCellCenter ) for( unsigned int d=0 ; d<Dim ; d++ ) simplex[k][d] *= res[d];
		else                             for( unsigned int d=0 ; d<Dim ; d++ ) simplex[k][d] *= res[d]-1;
	return simplex;
}

template< bool NodeAtCellCenter , typename SimplexFunctor /* = std::function< Simplex< double , Dim , Dim > ( size_t ) > */ >
RegularGrid< Dim , TexelInfo > GetNodeTexelInfo( size_t simplexNum , SimplexFunctor && SF , const unsigned int res[Dim] , bool forceThreadSafe )
{
	static_assert( std::is_convertible_v< SimplexFunctor , std::function< Simplex< double , Dim , Dim > ( size_t ) > > , "[ERROR] SimplexFunctor is poorly formed" );

	RegularGrid< Dim , TexelInfo > texelInfo( res );

	RegularGrid< Dim >::Range range;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.second[d] = res[d];

	if( forceThreadSafe )
	{
		std::vector< std::vector< std::pair< RegularGrid< Dim >::Index , TexelInfo > > > _texelInfo( ThreadPool::NumThreads() );
		ThreadPool::ParallelFor
		(
			0 , simplexNum ,
			[&]( unsigned int t , size_t sIdx )
			{
				Simplex< double , Dim , Dim > s = GetTexelSpaceSimplex< NodeAtCellCenter >( SF(sIdx) , res );
				auto RasterizationFunctor = [&]( RegularGrid< Dim >::Index I )
					{
						TexelInfo ti;
						ti.sIdx = (unsigned int)sIdx;
						ti.bc = s.barycentricCoordinates( TexelNodePosition< NodeAtCellCenter >( I ) );
						_texelInfo[t].push_back( std::make_pair( I , ti ) );
					};
				Rasterizer2D::RasterizeNodes< NodeAtCellCenter >( s , RasterizationFunctor , range );
			}
		);
		for( unsigned int i=0 ; i<_texelInfo.size() ; i++ ) for( unsigned int j=0 ; j<_texelInfo[i].size() ; j++ )
			texelInfo( _texelInfo[i][j].first ) = _texelInfo[i][j].second;
	}
	else
	{
		ThreadPool::ParallelFor
		(
			0 , simplexNum , 
			[&]( unsigned int , size_t sIdx )
			{
				Simplex< double , Dim , Dim > s = GetTexelSpaceSimplex< NodeAtCellCenter >( SF(sIdx) , res );
				auto RasterizationFunctor = [&]( typename RegularGrid< Dim >::Index I )
					{
						texelInfo( I ).sIdx = (unsigned int)sIdx;
						texelInfo( I ).bc = s.barycentricCoordinates( TexelNodePosition< NodeAtCellCenter >( I ) );
					};
				Rasterizer2D::RasterizeNodes< NodeAtCellCenter >( s , RasterizationFunctor , range );
			}
		);
	}

	return texelInfo;
}


template< bool Nearest , bool NodeAtCellCenter , typename SimplexFunctor /* = std::function< Simplex< double , Dim , Dim > ( size_t ) > */ >
RegularGrid< Dim , TexelInfo > GetSupportedTexelInfo( size_t simplexNum , SimplexFunctor && SF , const unsigned int res[Dim] )
{
	static_assert( IsValidSimplexFunctor< SimplexFunctor , Dim , Dim >() , "[ERROR] SimplexFunctor poorly formed" );
	RegularGrid< Dim , TexelInfo > texelInfo( res );
	RegularGrid< Dim >::Range range;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.second[d] = res[d];

	std::vector< std::vector< std::pair< RegularGrid< Dim >::Index , TexelInfo > > > _texelInfo( ThreadPool::NumThreads() );
	ThreadPool::ParallelFor
		(
			0 , simplexNum ,
			[&]( unsigned int t , size_t sIdx )
			{
				Simplex< double , Dim , Dim > s = GetTexelSpaceSimplex< NodeAtCellCenter >( SF(sIdx) , res );
				auto RasterizationFunctor = [&]( RegularGrid< Dim >::Index I )
					{
						TexelInfo ti;
						ti.sIdx = (unsigned int)sIdx;
						ti.bc = s.barycentricCoordinates( TexelNodePosition< NodeAtCellCenter >( I ) );
						_texelInfo[t].push_back( std::make_pair( I , ti ) );
					};
				Rasterizer2D::RasterizeSupports< Nearest , NodeAtCellCenter >( s , RasterizationFunctor , range );
			}
		);
	for( unsigned int i=0 ; i<_texelInfo.size() ; i++ ) for( unsigned int j=0 ; j<_texelInfo[i].size() ; j++ )
		texelInfo( _texelInfo[i][j].first ) = _texelInfo[i][j].second;

	return texelInfo;
}

#ifdef NEW_TEXEL_CODE
template< bool Nearest , bool NodeAtCellCenter , unsigned int K=Dim , typename SimplexFunctor = std::function< Simplex< double , Dim , K > ( size_t ) > >
RegularGrid< Dim , std::vector< size_t > > GetSupportedSimplexIndices( size_t simplexNum , SimplexFunctor && SF , const unsigned int res[Dim] )
#else // !NEW_TEXEL_CODE
template< bool Nearest , bool NodeAtCellCenter , typename SimplexFunctor /* = std::function< Simplex< double , Dim , Dim > ( size_t ) > */ >
RegularGrid< Dim , std::vector< size_t > > GetSupportedSimplexIndices( size_t simplexNum , SimplexFunctor && SF , const unsigned int res[Dim] )
#endif // NEW_TEXEL_CODE
{
#ifdef NEW_TEXEL_CODE
	static_assert( IsValidSimplexFunctor< SimplexFunctor , Dim , K >() , "[ERROR] SimplexFunctor poorly formed" );
#else // !NEW_TEXEL_CODE
	static_assert( IsValidSimplexFunctor< SimplexFunctor , Dim , Dim >() , "[ERROR] SimplexFunctor poorly formed" );
#endif // NEW_TEXEL_CODE

	RegularGrid< Dim , std::vector< size_t > > simplexIndices( res );
	RegularGrid< Dim >::Range range;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.second[d] = res[d];

	std::vector< std::vector< std::pair< RegularGrid< Dim >::Index , size_t > > > _simplexIndices( ThreadPool::NumThreads() );
	ThreadPool::ParallelFor
	(
		0 , simplexNum ,
		[&]( unsigned int t , size_t sIdx )
		{
			Simplex< double , Dim , K > s = GetTexelSpaceSimplex< NodeAtCellCenter >( SF( sIdx) , res );
			auto RasterizationFunctor = [&]( RegularGrid< Dim >::Index I ){ _simplexIndices[t].push_back( std::make_pair( I , sIdx ) ); };
			Rasterizer2D::RasterizeSupports< Nearest , NodeAtCellCenter >( s , RasterizationFunctor , range );
		}
	);
	for( unsigned int i=0 ; i<_simplexIndices.size() ; i++ ) for( unsigned int j=0 ; j<_simplexIndices[i].size() ; j++ )
		simplexIndices( _simplexIndices[i][j].first ).push_back( _simplexIndices[i][j].second );

	return simplexIndices;
}

template< bool Nearest , bool NodeAtCellCenter , typename SimplexFunctor /* = std::function< Simplex< double , Dim , Dim > ( size_t ) > */ >
size_t DilateTexelInfo( SimplexFunctor && SF , RegularGrid< Dim , TexelInfo > &texelInfo )
{
	static_assert( std::is_convertible_v< SimplexFunctor , std::function< Simplex< double , Dim , Dim > ( size_t ) > > , "[ERROR] SimplexFunctor poorly formed" );

	using Index = RegularGrid< Dim >::Index;
	using Range = RegularGrid< Dim >::Range;

	Range range;
	for( unsigned int d=0 ; d<Dim ; d++ ) range.second[d] = texelInfo.res(d);

	struct BoundaryInfo
	{
		Index I;
		std::vector< unsigned int > neighbors;
	};

	// Get the list of boundary texels and (interior) neighboring
	std::vector< std::vector< BoundaryInfo > > boundary( ThreadPool::NumThreads() );
	auto F = [&]( unsigned int t , Index I )
		{
			// If this is not an interior texel
			if( texelInfo(I).sIdx==-1 )
			{
				BoundaryInfo nInfo;
				Range::Intersect( Range(I).dilate(1) , range ).process( [&]( Index I ){ if( texelInfo(I).sIdx!=-1 ) nInfo.neighbors.push_back( texelInfo(I).sIdx ); } );
				if( nInfo.neighbors.size() )
				{
					nInfo.I = I;
					boundary[t].push_back( nInfo );
				}
			}
		};
	range.processParallel( F );

	for( unsigned int t=0 ; t<boundary.size() ; t++ )
	{
		std::vector< BoundaryInfo > &_boundary = boundary[t];

		// For each boundary texel
		ThreadPool::ParallelFor
		(
			0 , _boundary.size() ,
			[&]( unsigned int , size_t i )
			{
				Index I = _boundary[i].I;

				// Get the center of the texel
				Point< double , Dim > p = TexelNodePosition< NodeAtCellCenter >( I );

				// Finding the neighboring simplex closest to the texel center
				unsigned int sIdx = -1;
				double dist = std::numeric_limits< double >::infinity();
				for( unsigned int j=0 ; j<_boundary[i].neighbors.size() ; j++ )
				{
					Simplex< double , Dim , Dim > s = GetTexelSpaceSimplex< NodeAtCellCenter >( SF( _boundary[i].neighbors[j] ) , texelInfo.res() );
					double d = DistanceToTriangle( p , s );
					if( d<dist ) dist = d , sIdx = _boundary[i].neighbors[j];
				}
				if( sIdx==-1 ) MK_ERROR_OUT( "Could not find neighboring simplex" );

				Simplex< double , Dim , Dim > s = GetTexelSpaceSimplex< NodeAtCellCenter >( SF( sIdx ) , texelInfo.res() );

				texelInfo( I ).sIdx = sIdx;
				texelInfo( I ).bc = s.barycentricCoordinates( p );
			}
		);
	}

	return boundary.size();
}


template< unsigned int EmbeddingDim , typename VertexEmbeddingFunctor /* = std::function< Point< double , EmbeddingDim > ( size_t ) > */ , typename SimplexIndexFunctor /* = std::function< SimplexIndex< Dim >( size_t ) > */ >
void FinalizeTexelInfo( size_t simplexNum , VertexEmbeddingFunctor && VF , SimplexIndexFunctor && SIF , RegularGrid< Dim , TexelInfo > &texelInfo )
{
	static_assert( std::is_convertible_v< VertexEmbeddingFunctor , std::function< Point< double , EmbeddingDim > ( size_t ) > > , "[ERROR] VertexEmbeddingFunctor poorly formed" );
	static_assert( std::is_convertible_v< SimplexIndexFunctor , std::function< SimplexIndex< Dim > ( size_t ) > > , "[ERROR] SimplexIndexFunctor poorly formed" );

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
					// Get the representation of the texel center in the tangent frame of the closest triangle
					FEM::HermiteSamplePoint< double > h;
					h.tIdx = texelInfo[i].sIdx;
					for( unsigned int d=0 ; d<Dim ; d++ ) h.p[d] = 1./(Dim+1) , h.v[d] = texelInfo[i].bc[d+1] - 1./(Dim+1);

					// Walk along the tangent direction
					rMesh.exp( xForms , h );
					texelInfo[i].sIdx = h.tIdx;
					texelInfo[i].bc[0] = 1.;
					for( unsigned int d=0 ; d<Dim ; d++ ) texelInfo[i].bc[d+1] = h.p[d] , texelInfo[i].bc[0] -= h.p[d];
				}
			}
		}
	);
	DeletePointer( xForms );
}


template< unsigned int EmbeddingDim , bool Nearest , bool NodeAtCellCenter , typename VertexEmbeddingFunctor /* = std::function< Point< double , EmbeddingDim > ( size_t ) > */ , typename SimplexIndexFunctor /* = std::function< SimplexIndex< Dim >( size_t ) > */ , typename SimplexFunctor /* = std::function< Simplex< double , Dim , Dim > ( size_t ) > */ >
RegularGrid< Dim , TexelInfo > GetSupportedTexelInfo( size_t simplexNum , VertexEmbeddingFunctor && VF , SimplexIndexFunctor && SIF , SimplexFunctor && SF , const unsigned int res[Dim] , unsigned int dilationRadius , bool finalize , bool verbose )
{
	static_assert( IsValidVertexFunctor< VertexEmbeddingFunctor , EmbeddingDim >() , "[ERROR] VertexEmbeddingFunctor poorly formed" );
	static_assert( IsValidSimplexIndexFunctor< SimplexIndexFunctor , Dim >()       , "[ERROR] SimplexFunctor poorly formed" );
	static_assert( IsValidSimplexFunctor< SimplexFunctor , Dim , Dim >()           , "[ERROR] SimplexEmbeddingFunctor poorly formed" );

	Miscellany::PerformanceMeter pMeter;

	auto ActiveTexelCount = []( const RegularGrid< Dim , TexelInfo > &texelInfo )
		{
			size_t count = 0;
			for( size_t i=0 ; i<texelInfo.resolution() ; i++ ) if( texelInfo[i].sIdx!=-1 ) count++;
			return count;
		};

	RegularGrid< Dim , TexelInfo > texelInfo = GetSupportedTexelInfo< Nearest , NodeAtCellCenter >( simplexNum , SF , res );

	if( dilationRadius )
	{
		for( unsigned int i=0 ; i<dilationRadius ; i++ )
		{
			Miscellany::PerformanceMeter pMeter;
			DilateTexelInfo< Nearest , NodeAtCellCenter >( SF , texelInfo );
			if( verbose ) std::cout << pMeter( "Radius: " + std::to_string( i+1 ) + " / " + std::to_string( dilationRadius ) ) << std::endl;
		}
		if( verbose )
		{
			std::stringstream sStream;
			sStream << ActiveTexelCount( texelInfo );
			std::cout << pMeter( "Dilation" , sStream.str() ) << std::endl;
		}
	}

	if( finalize )
	{
		FinalizeTexelInfo< EmbeddingDim >( simplexNum , std::forward< VertexEmbeddingFunctor >( VF ) , std::forward< SimplexIndexFunctor >( SIF ) , texelInfo );
		if( verbose )
		{
			std::stringstream sStream;
			sStream << ActiveTexelCount( texelInfo );
			std::cout << pMeter( "Finalized" , sStream.str() ) << std::endl;
		}
	}

	return texelInfo;
}
