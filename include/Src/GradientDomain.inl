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


template< typename Real >
template< typename Functor >
constexpr bool GradientDomain< Real >::_IsTriangleCornerFunctor( void ){ return std::is_convertible_v< Functor , std::function< size_t ( size_t , unsigned int ) > >; }

template< typename Real >
template< typename Functor >
constexpr bool GradientDomain< Real >::_IsSurfaceMetricFunctor( void ){ return std::is_convertible_v< Functor , std::function< SquareMatrix< Real , 2 > ( size_t ) > >; }

template< typename Real >
template< typename Functor >
constexpr bool GradientDomain< Real >::_IsSurfaceVertexFunctor( void ){ return std::is_convertible_v< Functor , std::function< Point< Real , 3 > ( size_t ) > >; }

template< typename Real >
template< typename Functor >
constexpr bool GradientDomain< Real >::_IsSurfaceVertexOrMetricFunctor( void ){ return _IsSurfaceVertexFunctor< Functor >() || _IsSurfaceMetricFunctor< Functor >(); }

template< typename Real >
template< typename Functor >
constexpr bool GradientDomain< Real >::_IsTextureVertexFunctor( void ){ return std::is_convertible_v< Functor , std::function< Point< Real , 2 > ( size_t ) > >; }

template< typename Real >
template
<
	typename SurfaceCornerFunctor ,         /* = std::function< size_t ( size_t , unsigned int ) > */
	typename SurfaceVertexOrMetricFunctor , /* = std::function< Point< Real , 3 > ( size_t ) > || std::function< SquareMatrix< Real , 2 > ( size_t ) > */
	typename TextureCornerFunctor ,         /* = std::function< size_t ( size_t , unsigned int ) > */
	typename TextureVertexFunctor           /* = std::function< Point< Real , 2 > ( size_t ) > */
>
GradientDomain< Real >::GradientDomain
(
	unsigned quadraturePointsPerTriangle ,
	size_t numTriangles ,
	size_t numSurfaceVertices ,
	size_t numTextureVertices ,
	SurfaceCornerFunctor && surfaceCornerFunctor ,
	SurfaceVertexOrMetricFunctor && surfaceVertexOrMetricFunctor ,
	TextureCornerFunctor && textureCornerFunctor ,
	TextureVertexFunctor && textureVertexFunctor ,
	unsigned int width ,
	unsigned int height ,
	bool normalize
)
{
	static_assert( _IsTriangleCornerFunctor< SurfaceCornerFunctor >()                , "[ERROR] SurfaceCornerFunctor poorly formed" );
	static_assert( _IsSurfaceVertexOrMetricFunctor< SurfaceVertexOrMetricFunctor >() , "[ERROR] SurfaceVertexOrMetricFunctor poorly formed" );
	static_assert( _IsTriangleCornerFunctor< TextureCornerFunctor >()                , "[ERROR] TextureCornerFunctor poorly formed" );
	static_assert( _IsTextureVertexFunctor< TextureVertexFunctor >()                 , "[ERROR] TextureVertexFunctor poorly formed" );

	static const bool HasSurfaceMetric = _IsSurfaceMetricFunctor< SurfaceVertexOrMetricFunctor >();
	static const bool HasSurfaceVertex = _IsSurfaceVertexFunctor< SurfaceVertexOrMetricFunctor >();

	TexturedTriangleMesh< Real > mesh;
	{
		mesh.surface.triangles.resize( numTriangles );
		mesh.texture.triangles.resize( numTriangles );
		if constexpr( HasSurfaceVertex ) mesh.surface.vertices.resize( numSurfaceVertices );
		mesh.texture.vertices.resize( numTextureVertices );
		for( size_t t=0 ; t<numTriangles ; t++ ) for( unsigned int k=0 ; k<=2 ; k++ )
			mesh.surface.triangles[t][k] = surfaceCornerFunctor( t , k ) , mesh.texture.triangles[t][k] = textureCornerFunctor( t , k );
		if constexpr( HasSurfaceVertex ) for( size_t i=0 ; i<numSurfaceVertices ; i++ ) mesh.surface.vertices[i] = surfaceVertexOrMetricFunctor( i );
		for( size_t i=0 ; i<numTextureVertices ; i++ ) mesh.texture.vertices[i] = textureVertexFunctor( i );

		// Flip the vertical axis
		for( size_t i=0 ; i<mesh.texture.vertices.size() ; i++ ) mesh.texture.vertices[i][1] = (Real)1. - mesh.texture.vertices[i][1];
		for( size_t i=0 ; i<mesh.texture.triangles.size() ; i++ ) if( mesh.texture.triangle( static_cast< unsigned int >(i) ).measure()==0 ) MK_WARN( "Zero area texture triangle: " , i );
	}

	HierarchicalSystem< Real , Real > hierarchy;
	ExplicitIndexVector< ChartIndex , AtlasChart< Real > > atlasCharts;
	MultigridBlockInfo multigridBlockInfo;

	InitializeHierarchy( mesh , width , height , 1 , _textureNodes , hierarchy , atlasCharts , multigridBlockInfo );
	ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< Real , 2 > > > parameterMetric;
	if constexpr( HasSurfaceVertex ) InitializeMetric( mesh , EMBEDDING_METRIC , atlasCharts , parameterMetric , normalize );
	else if constexpr( HasSurfaceMetric )
	{
		ExplicitIndexVector< AtlasMeshTriangleIndex , SquareMatrix< Real , 2 > > surfaceMetric( numTriangles );
		for( size_t t=0 ; t<numTriangles ; t++ ) surfaceMetric[t] = surfaceVertexOrMetricFunctor(t);
		if( normalize )
		{
			Real totalArea = 0;
			ThreadPool::ParallelFor( 0 , surfaceMetric.size() , [&]( size_t t ){ Atomic< Real >::Add( totalArea , (Real)sqrt( fabs( surfaceMetric[t].determinant() ) / 2. ) ); } );
			ThreadPool::ParallelFor( 0 , surfaceMetric.size() , [&]( size_t t ){ surfaceMetric[t] /= totalArea; } );
		}
		InitializeParameterMetric( mesh , surfaceMetric , atlasCharts , parameterMetric );
	}
	OperatorInitializer::Initialize( quadraturePointsPerTriangle , _massAndStiffnessOperators , hierarchy.gridAtlases[0] , parameterMetric , atlasCharts , _divergenceOperator );
}

template< typename Real >
size_t GradientDomain< Real >::numNodes( void ) const
{
	return _massAndStiffnessOperators.indexConverter.numCombined();
}

template< typename Real >
size_t GradientDomain< Real >::numEdges( void ) const
{
	return _divergenceOperator.edges.size();
}

template< typename Real >
std::pair< unsigned int , unsigned int > GradientDomain< Real >::node( size_t n ) const
{
	return std::pair< unsigned int , unsigned int >( _textureNodes[n].ci , _textureNodes[n].cj );
}

template< typename Real >
bool GradientDomain< Real >::isCovered( size_t n ) const
{
	return _textureNodes[n].isInterior;
}

template< typename Real >
bool GradientDomain< Real >::isChartCrossing( size_t e ) const
{
	return _textureNodes[ _divergenceOperator.edges[e][0] ].chartID != _textureNodes[ _divergenceOperator.edges[e][1] ].chartID;
}

template< typename Real >
std::pair< size_t , size_t > GradientDomain< Real >::edge( size_t e ) const
{
	return std::pair< size_t , size_t >( static_cast< size_t >( _divergenceOperator.edges[e][0] ) , static_cast< size_t >( _divergenceOperator.edges[e][1] ) );
}

template< typename Real >
Eigen::SparseMatrix< Real > GradientDomain< Real >::mass( void ) const
{
	return _massAndStiffnessOperators.mass();
}

template< typename Real >
Eigen::SparseMatrix< Real > GradientDomain< Real >::stiffness( void ) const
{
	return _massAndStiffnessOperators.stiffness();
}

template< typename Real >
Eigen::SparseMatrix< Real > GradientDomain< Real >::divergence( void ) const
{
	return _divergenceOperator();
}

template< typename Real >
Eigen::SparseMatrix< Real > GradientDomain< Real >::finiteDifferences( void ) const
{
	std::vector< Eigen::Triplet< Real > > triplets( 2*numEdges() );

	ThreadPool::ParallelFor
		(
			0 , numEdges() ,
			[&]( size_t e )
			{
				std::pair< size_t , size_t > edge = this->edge( e );
				triplets[2*e+0] = Eigen::Triplet< Real >( static_cast< int >( e ) , static_cast< int >( edge.first  ) , static_cast< Real >(-1) );
				triplets[2*e+1] = Eigen::Triplet< Real >( static_cast< int >( e ) , static_cast< int >( edge.second ) , static_cast< Real >( 1) );
			}
		);

	Eigen::SparseMatrix< Real > D( numEdges() , numNodes() );
	D.setFromTriplets( triplets.begin() , triplets.end() );
	return D;
}


template< typename Real >
template< typename T >
void GradientDomain< Real >::mass( const T * in , T * out ) const
{
	return _massAndStiffnessOperators.mass( in , out );
}

template< typename Real >
template< typename T >
void GradientDomain< Real >::stiffness( const T * in , T * out ) const
{
	return _massAndStiffnessOperators.stiffness( in , out );
}

template< typename Real >
template< typename T >
void GradientDomain< Real >::divergence( const T * in , T * out ) const
{
	return _divergenceOperator( in , out );
}

template< typename Real >
template< typename T >
void GradientDomain< Real >::finiteDifferences( const T * in , T * out ) const
{
	ThreadPool::ParallelFor
	(
		0 , numEdges() ,
		[&]( size_t e )
		{
			std::pair< size_t , size_t > edge = this->edge( e );
			out[e] = in[ edge.second ] - in[ edge.first ];
		}
	);
}

template< typename Real >
void GradientDomain< Real >::unitTests( unsigned int numTests , double eps , bool verbose ) const
{
	Eigen::SparseMatrix< Real > M = mass() , S = stiffness() , Div = divergence() , FD = finiteDifferences();

	auto NormalizedDifference = []( const Eigen::VectorXd &x1 , const Eigen::VectorXd &x2 )
		{
			double n = ( x1.squaredNorm() + x2.squaredNorm() ) / 2.;
			return sqrt( ( x1 - x2 ).squaredNorm() / n );
		};


	// Check that the matrices make sense
	{
		Eigen::VectorXd one( numNodes() );
		for( size_t i=0 ; i<numNodes() ; i++ ) one[i] = 1;

		// Check that the area is equal to one (assuming the mass is computed from the area-normalized mesh)
		{
			double area = ( M * one ).dot( one );
			if( verbose || fabs( area-1. )>eps ) std::cout << "Area: " << area << std::endl;
		}

		// Check that constant functions vanish
		{
			double n1 = ( S * one ).squaredNorm() , n2 = ( FD * one ).squaredNorm();
			n1 = sqrt( n1 / S.squaredNorm() );
			n2 = sqrt( n2 / FD.squaredNorm() );
			if( n1>eps || verbose ) std::cout << "Constant (stiffness): " << n1 << std::endl;
			if( n2>eps || verbose ) std::cout << "Constant (finite-differences): " << n2 << std::endl;
		}

		// Check that the stiffness matrix is factored as the product of the finite-difference matrix and the divergence operator
		{
			Eigen::SparseMatrix< Real > DivFD = Div * FD;
			double n = ( S.squaredNorm() + ( DivFD ).squaredNorm() ) / 2. , d = ( S - DivFD ).squaredNorm();
			d = sqrt( d / n );
			if( verbose || d>eps ) std::cout << "Stiffness factorization: " << d << std::endl;
		}
	}

	// Check that the matrices are consistent with the evaluation
	for( unsigned int n=0 ; n<numTests ; n++ )
	{
		Eigen::VectorXd x( numNodes() ) , e( numEdges() );
		for( size_t i=0 ; i<numNodes() ; i++ ) x[i] = Random< Real >();
		for( size_t i=0 ; i<numEdges() ; i++ ) e[i] = Random< Real >();
		{
			Eigen::VectorXd b( numNodes() );
			mass( &x[0] , &b[0] );
			double diff = NormalizedDifference( b , M * x );
			if( verbose || diff>eps ) std::cout << "\t" << n << "] Mass: " << diff << std::endl;
		}
		{
			Eigen::VectorXd b( numNodes() );
			stiffness( &x[0] , &b[0] );
			double diff = NormalizedDifference( b , S * x );
			if( verbose || diff>eps ) std::cout << "\t" << n << "] Stiffness: " << diff << std::endl;
		}
		{
			Eigen::VectorXd b( numEdges() );
			finiteDifferences( &x[0] , &b[0] );
			double diff = NormalizedDifference( b , FD * x );
			if( verbose || diff>eps ) std::cout << "\t" << n << "] Finite-differences: " << diff << std::endl;
		}
		{
			Eigen::VectorXd b( numNodes() );
			divergence( &e[0] , &b[0] );
			double diff = NormalizedDifference( b , Div * e );
			if( verbose || diff>eps ) std::cout << "\t" << n << "] Divergence: " << diff << std::endl;
		}
	}
}