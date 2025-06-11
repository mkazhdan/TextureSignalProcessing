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
constexpr bool GradientDomain< Real >::_IsTriangleFunctor( void ){ return std::is_convertible_v< Functor , std::function< SimplexIndex< 2 > ( size_t ) > >; }

template< typename Real >
template< typename Functor >
constexpr bool GradientDomain< Real >::_IsSurfaceMetricFunctor( void ){ return std::is_convertible_v< Functor , std::function< SquareMatrix< Real , 2 > ( size_t ) > >; }

template< typename Real >
template< typename Functor >
constexpr bool GradientDomain< Real >::_IsSurfaceVertexFunctor( void ){ return std::is_convertible_v< Functor , std::function< Point< Real , 3 > ( size_t ) > >; }

template< typename Real >
template< typename Functor >
constexpr bool GradientDomain< Real >::_IsSurfaceMetricOrVertexFunctor( void ){ return _IsSurfaceMetricFunctor< Functor >() || _IsSurfaceVertexFunctor< Functor >(); }

template< typename Real >
template< typename Functor >
constexpr bool GradientDomain< Real >::_IsTextureVertexFunctor( void ){ return std::is_convertible_v< Functor , std::function< Point< Real , 2 > ( size_t ) > >; }


template< typename Real >
template
<
	typename SurfaceTriangleFunctor ,       /* = std::function< SimplexIndex< 2 > ( size_t ) > */
	typename SurfaceMetricOrVertexFunctor , /* = std::function< SquareMatrix< Real , 2 > ( size_t ) > || std::function< Point< Real ,32 > ( size_t ) > */
	typename TextureTriangleFunctor ,       /* = std::function< SimplexIndex< 2 > ( size_t ) > */
	typename TextureVertexFunctor           /* = std::function< Point< Real , 2 > ( size_t ) > */
>
GradientDomain< Real >::GradientDomain
(
	size_t numTriangles ,
	size_t numSurfaceVertices ,
	size_t numTextureVertices ,
	unsigned quadraturePointsPerTriangle ,
	SurfaceTriangleFunctor && surfaceTriangleFunctor ,
	SurfaceMetricOrVertexFunctor && surfaceMetricOrVertexFunctor ,
	TextureTriangleFunctor && textureTriangleFunctor ,
	TextureVertexFunctor   && textureVertexFunctor ,
	unsigned int width ,
	unsigned int height ,
	bool normalize
)
{
	static_assert( _IsTriangleFunctor< SurfaceTriangleFunctor >()                    , "[ERROR] SurfaceTriangleFunctor poorly formed" );
	static_assert( _IsSurfaceMetricOrVertexFunctor< SurfaceMetricOrVertexFunctor >() , "[ERROR] SurfaceMetricOrVertexFunctor poorly formed" );
	static_assert( _IsTriangleFunctor< TextureTriangleFunctor >()                    , "[ERROR] TextureTriangleFunctor poorly formed" );
	static_assert( _IsTextureVertexFunctor< TextureVertexFunctor >()                 , "[ERROR] TextureVertexFunctor poorly formed" );

	static const bool UseSurfaceMetric = _IsSurfaceMetricFunctor< SurfaceMetricOrVertexFunctor >();
	static const bool HasSurfaceVertex = _IsSurfaceVertexFunctor< SurfaceMetricOrVertexFunctor >();

	TexturedTriangleMesh< Real > mesh;
	{
		mesh.surface.triangles.resize( numTriangles );
		mesh.texture.triangles.resize( numTriangles );
		if constexpr( HasSurfaceVertex ) mesh.surface.vertices.resize( numSurfaceVertices );
		mesh.texture.vertices.resize( numTextureVertices );
		for( size_t i=0 ; i<numTriangles ; i++ ) mesh.surface.triangles[i] = surfaceTriangleFunctor( i ) , mesh.texture.triangles[i] = textureTriangleFunctor( i );
		if constexpr( HasSurfaceVertex ) for( size_t i=0 ; i<numSurfaceVertices ; i++ ) mesh.surface.vertices[i] = surfaceMetricOrVertexFunctor( i );
		for( size_t i=0 ; i<numTextureVertices ; i++ ) mesh.texture.vertices[i] = textureVertexFunctor( i );

		// Flip the vertical axis
		for( size_t i=0 ; i<mesh.texture.vertices.size() ; i++ ) mesh.texture.vertices[i][1] = (Real)1. - mesh.texture.vertices[i][1];
		for( size_t i=0 ; i<mesh.texture.triangles.size() ; i++ ) if( mesh.texture.triangle(i).measure()==0 ) MK_WARN( "Zero area texture triangle: " , i );
	}

	HierarchicalSystem< Real , Real > hierarchy;
	ExplicitIndexVector< ChartIndex , AtlasChart< Real > > atlasCharts;
	MultigridBlockInfo multigridBlockInfo;

	InitializeHierarchy( mesh , width , height , 1 , _textureNodes , hierarchy , atlasCharts , multigridBlockInfo );

	ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< Real , 2 > > > parameterMetric;
	if constexpr( HasSurfaceVertex ) InitializeMetric( mesh , EMBEDDING_METRIC , atlasCharts , parameterMetric );
	else
	{
		ExplicitIndexVector< AtlasMeshTriangleIndex , SquareMatrix< Real , 2 > > surfaceMetric( numTriangles );
		for( size_t t=0 ; t<numTriangles ; t++ ) surfaceMetric[t] = surfaceMetricOrVertexFunctor[t];
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
