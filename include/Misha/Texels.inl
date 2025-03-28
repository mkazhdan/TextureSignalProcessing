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


template< typename Real , unsigned int Dim , typename VertexEmbeddingFunctor /* = std::function< Point< double , Dim > ( size_t ) > */ , typename SimplexFunctor /* = std::function< SimplexIndex< K >( size_t ) > */ >
Point< Real , Dim > TexelInfo::position( const VertexEmbeddingFunctor & VF , const SimplexFunctor & SF ) const
{
	SimplexIndex< K > si = SF( sIdx );
	Simplex< double , Dim , K > s;
	for( unsigned int k=0 ; k<=K ; k++ ) s[k] = VF( si[k] );
	return Point< Real , Dim >( s( bc ) );
}


///////////////
///////////////

inline double DistanceToEdge( Point< double , K > p , Point< double , K > v0 , Point< double , K > v1 )
{
	// E(s) = || p - ( v0*(1-s) + v1*s ) ||^2
	// 0 = E'(s) = < p - ( v0*(1-s) + v1*s ) , v0 - v1 >
	// =>   < p - v0 , v0 - v1 > = - s * || v0 - v1 ||^2
	// =>   s = < p - v0  , v1-v0 > / || v1 - v0 ||^2
	return Point< double , K >::Dot( p-v0 , v1-v0 ) / Point< double , K >::SquareNorm( v1-v0 );
}

inline double DistanceToTriangle( Point< double , K > p , Simplex< double , K , K > s )
{
	Point< double , K+1 > bc = s.barycentricCoordinates( p );
	if( bc[0]>=0 && bc[1]>=0 && bc[2]>=0 ) return 0.;

	if     ( bc[0]<0 && bc[1]<0 ) return Point< double , K >::Length( p - s[0] );
	else if( bc[0]<0 && bc[2]<0 ) return Point< double , K >::Length( p - s[2] );
	else if( bc[1]<0 && bc[2]<0 ) return Point< double , K >::Length( p - s[1] );
	else if( bc[0]<0 ) return DistanceToEdge( p , s[0] , s[2] );
	else if( bc[1]<0 ) return DistanceToEdge( p , s[1] , s[0] );
	else if( bc[2]<0 ) return DistanceToEdge( p , s[2] , s[1] );
	else MK_ERROR_OUT( "At least one of the barycentric coordinates should be positive" );
	return 0.;
}

template< typename Real , unsigned int Dim , typename VertexEmbeddingFunctor /* = std::function< Point< double , Dim > ( size_t ) > */ , typename SimplexFunctor /* = std::function< SimplexIndex< K >( size_t ) > */ >
RegularGrid< K , Point< Real , Dim > > GetTexelPositions( size_t simplexNum , VertexEmbeddingFunctor && VF , SimplexFunctor && SF , const RegularGrid< K , TexelInfo > &texelInfo )
{
	static_assert( std::is_convertible_v< VertexEmbeddingFunctor , std::function< Point< double , Dim >( size_t ) > > , "[ERROR] VertexEmbeddingFunctor poorly formed" );
	static_assert( std::is_convertible_v< SimplexFunctor , std::function< SimplexIndex< K >( size_t ) > > , "[ERROR] SimplexFunctor poorly formed" );

	RegularGrid< K , Point< float , Dim > > texturePositions( texelInfo.res() );

	Point< Real , Dim > badPosition;
	for( unsigned int d=0 ; d<Dim ; d++ ) badPosition[d] = std::numeric_limits< Real >::infinity();

	ThreadPool::ParallelFor
	(
		0 , texelInfo.resolution() ,
		[&]( unsigned int , size_t i )
		{
			if( texelInfo[i].sIdx==-1 ) texturePositions[i] = badPosition;
			else texturePositions[i] = texelInfo[i].template position< Real , Dim >( VF , SF );
		}
	);
	return texturePositions;
}

template< bool NodeAtCellCenter >
Point< double , K > TexelNodePosition( typename RegularGrid< K >::Index I )
{
	Point< double , K > p;
	if constexpr( NodeAtCellCenter ) for( unsigned int k=0 ; k<K ; k++ ) p[k] = I[k] + 0.5;
	else                             for( unsigned int k=0 ; k<K ; k++ ) p[k] = I[k] + 0.0;
	return p;
}

template< bool NodeAtCellCenter , typename UVCoordinateFunctor /* = std::function< Point< double , K > ( unsigned int , unsigned int ) > */ >
Simplex< double , K , K > GetSimplex( unsigned int sIdx , const UVCoordinateFunctor & UV , const unsigned int res[K] )
{
	Simplex< double , K , K > s;
	for( unsigned int k=0 ; k<=K ; k++ )
	{
		s[k] = UV( sIdx , k );
		if constexpr( NodeAtCellCenter ) for( unsigned int _k=0 ; _k<K ; _k++ ) s[k][_k] *= res[_k];
		else                             for( unsigned int _k=0 ; _k<K ; _k++ ) s[k][_k] *= res[_k]-1;
	}
	return s;
}

template< bool NodeAtCellCenter , typename UVCoordinateFunctor /* = std::function< Point< double , K > ( unsigned int , unsigned int ) > */ >
RegularGrid< K , TexelInfo > GetNodeTexelInfo( size_t simplexNum , UVCoordinateFunctor && UV , unsigned int width , unsigned int height )
{
	static_assert( std::is_convertible_v< UVCoordinateFunctor , std::function< Point< double , K >( unsigned int , unsigned int ) > > , "[ERROR] UVCoordinateFunctor is poorly formed" );

	const unsigned int res[] = { width , height };

	RegularGrid< K , TexelInfo > texelInfo;
	texelInfo.resize( res );

	RegularGrid< K >::Range range;
	for( unsigned int k=0 ; k<K ; k++ ) range.second[k] = res[k];

	ThreadPool::ParallelFor
	(
		0 , simplexNum , 
		[&]( unsigned int , size_t sIdx )
		{
			Simplex< double , K , K > s = GetSimplex< NodeAtCellCenter >( (unsigned int)sIdx , UV , res );
			auto RasterizationFunctor = [&]( typename RegularGrid< K >::Index I )
				{
					texelInfo( I ).sIdx = (unsigned int)sIdx;
					texelInfo( I ).bc = s.barycentricCoordinates( TexelNodePosition< NodeAtCellCenter >( I ) );
//					for( unsigned int k=0 ; k<=K ; k++ ) if( texelInfo( I ).bc[k]<0 || texelInfo(I).bc[k]>1 ) MK_ERROR_OUT( "Bad barycentric coordinate:" , texelInfo(I).bc );
				};
			Rasterizer2D::RasterizeNodes< NodeAtCellCenter >( s , RasterizationFunctor , range );
		}
	);

	return texelInfo;
}


template< bool Nearest , bool NodeAtCellCenter , typename UVCoordinateFunctor /* = std::function< Point< double , K > ( unsigned int , unsigned int ) > */ >
RegularGrid< K , TexelInfo > GetSupportedTexelInfo( size_t simplexNum , UVCoordinateFunctor && UV , unsigned int width , unsigned int height )
{
	const unsigned int res[] = { width , height };
	RegularGrid< K , Texels::TexelInfo > texelInfo( res );
	RegularGrid< K >::Range range;
	for( unsigned int k=0 ; k<K ; k++ ) range.second[k] = res[k];

	for( size_t sIdx=0 ; sIdx<simplexNum ; sIdx++ )
	{
		Simplex< double , K , K > s = GetSimplex< NodeAtCellCenter >( (unsigned int)sIdx , UV , res );
		auto RasterizationFunctor = [&]( RegularGrid< K >::Index I )
			{
				texelInfo( I ).sIdx = (unsigned int)sIdx;
				texelInfo( I ).bc = s.barycentricCoordinates( TexelNodePosition< NodeAtCellCenter >( I ) );
			};
		Rasterizer2D::RasterizeSupports< Nearest , NodeAtCellCenter >( s , RasterizationFunctor , range );
	}
	return texelInfo;
}

template< bool Nearest , bool NodeAtCellCenter , typename UVCoordinateFunctor /* = std::function< Point< double , K > ( unsigned int , unsigned int ) > */ >
size_t DilateTexelInfo( UVCoordinateFunctor && UV , RegularGrid< K , TexelInfo > &texelInfo )
{
	using Index = RegularGrid< K >::Index;
	using Range = RegularGrid< K >::Range;

	Range range;
	for( unsigned int k=0 ; k<K ; k++ ) range.second[k] = texelInfo.res(k);

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
				Point< double , K > p = TexelNodePosition< NodeAtCellCenter >( I );

				// Finding the neighboring simplex closest to the texel center
				unsigned int sIdx = -1;
				double dist = std::numeric_limits< double >::infinity();
				for( unsigned int j=0 ; j<_boundary[i].neighbors.size() ; j++ )
				{
					Simplex< double , K , K > s = GetSimplex< NodeAtCellCenter >( _boundary[i].neighbors[j] , UV , texelInfo.res() );
					double d = DistanceToTriangle( p , s );
					if( d<dist ) dist = d , sIdx = _boundary[i].neighbors[j];
				}
				if( sIdx==-1 ) MK_ERROR_OUT( "Could not find neighboring simplex" );

				Simplex< double , K , K > s = GetSimplex< NodeAtCellCenter >( sIdx , UV , texelInfo.res() );

				texelInfo( I ).sIdx = sIdx;
				texelInfo( I ).bc = s.barycentricCoordinates( p );
			}
		);
	}

	return boundary.size();
}


template< unsigned int Dim , typename VertexEmbeddingFunctor /* = std::function< Point< double , Dim > ( size_t ) > */ , typename SimplexFunctor /* = std::function< SimplexIndex< K >( size_t ) > */ >
void FinalizeTexelInfo( size_t simplexNum , VertexEmbeddingFunctor && VF , SimplexFunctor && SF , RegularGrid< K , TexelInfo > &texelInfo )
{
	static_assert( std::is_convertible_v< VertexEmbeddingFunctor , std::function< Point< double , Dim > ( size_t ) > > , "[ERROR] VertexEmbeddingFunctor poorly formed" );
	static_assert( std::is_convertible_v< SimplexFunctor , std::function< SimplexIndex< K > ( size_t ) > > , "[ERROR] SimplexFunctor poorly formed" );

	std::vector< SimplexIndex< K > > simplices( simplexNum );
	for( size_t i=0 ; i<simplexNum ; i++ ) simplices[i] = SF( i );
	FEM::RiemannianMesh< double , unsigned int > rMesh( GetPointer( simplices ) , simplexNum );
	rMesh.template setMetricFromEmbedding< Dim >( [&]( unsigned int v ){ return VF(v); } );
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
					for( unsigned int k=0 ; k<K ; k++ ) h.p[k] = 1./(K+1) , h.v[k] = texelInfo[i].bc[k+1] - 1./(K+1);

					// Walk along the tangent direction
					rMesh.exp( xForms , h );
					texelInfo[i].sIdx = h.tIdx;
					texelInfo[i].bc[0] = 1.;
					for( unsigned int k=0 ; k<K ; k++ ) texelInfo[i].bc[k+1] = h.p[k] , texelInfo[i].bc[0] -= h.p[k];
				}
			}
		}
	);
	DeletePointer( xForms );
}


template< unsigned int Dim , bool Nearest , bool NodeAtCellCenter , typename VertexEmbeddingFunctor /* = std::function< Point< double , Dim > ( size_t ) > */ , typename SimplexFunctor /* = std::function< SimplexIndex< K >( size_t ) > */ , typename UVCoordinateFunctor /* = std::function< Point< double , K > ( unsigned int , unsigned int ) > */ >
RegularGrid< K , TexelInfo > GetSupportedTexelInfo( size_t simplexNum , VertexEmbeddingFunctor && VF , SimplexFunctor && SF , UVCoordinateFunctor && UV , unsigned int width , unsigned int height , unsigned int dilationRadius , bool finalize , bool verbose )
{
	static_assert( std::is_convertible_v< VertexEmbeddingFunctor , std::function< Point< double , Dim > ( size_t ) > >           , "[ERROR] VertexEmbeddingFunctor poorly formed" );
	static_assert( std::is_convertible_v< SimplexFunctor , std::function< SimplexIndex< K > ( size_t ) > >                       , "[ERROR] SimplexFunctor poorly formed" );
	static_assert( std::is_convertible_v< UVCoordinateFunctor , std::function< Point< double , K > ( size_t , unsigned int ) > > , "[ERROR] UVCoordinateFunctor poorly formed" );

	Miscellany::PerformanceMeter pMeter;

	auto ActiveTexelCount = []( const RegularGrid< K , Texels::TexelInfo > &texelInfo )
		{
			size_t count = 0;
			for( size_t i=0 ; i<texelInfo.resolution() ; i++ ) if( texelInfo[i].sIdx!=-1 ) count++;
			return count;
		};

	RegularGrid< K , Texels::TexelInfo > texelInfo = GetSupportedTexelInfo< Nearest , NodeAtCellCenter >( simplexNum , UV , width , height );

	if( dilationRadius )
	{
		for( unsigned int i=0 ; i<dilationRadius ; i++ )
		{
			Miscellany::PerformanceMeter pMeter;
			DilateTexelInfo< Nearest , NodeAtCellCenter >( UV , texelInfo );
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
		FinalizeTexelInfo< Dim >( simplexNum , std::forward< VertexEmbeddingFunctor >( VF ) , std::forward< SimplexFunctor >( SF ) , texelInfo );
		if( verbose )
		{
			std::stringstream sStream;
			sStream << ActiveTexelCount( texelInfo );
			std::cout << pMeter( "Finalized" , sStream.str() ) << std::endl;
		}
	}

	return texelInfo;
}
