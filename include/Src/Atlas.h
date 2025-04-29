/*
Copyright (c) 2018, Fabian Prada and Michael Kazhdan
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

#ifndef ATLAS_MESH_INLCUDED
#define ATLAS_MESH_INLCUDED

#define DEBUG_ATLAS

#include <set>
#include <Misha/Miscellany.h>
#include <Misha/Exceptions.h>
#include "SimpleTriangleMesh.h"
#ifdef NEW_INDEXING
#include "Indices.h"
#endif // NEW_INDEXING

namespace MishaK
{
	template< typename GeometryReal > struct AtlasChart;

	template< typename GeometryReal >
#ifdef NEW_INDEXING
	struct AtlasMesh : protected SimpleTriangleMesh< GeometryReal , 2 >
#else // !NEW_INDEXING
	struct AtlasMesh : public SimpleTriangleMesh< GeometryReal , 2 >
#endif // NEW_INDEXING
	{
#ifdef NEW_INDEXING
		Point< GeometryReal , 2 > vertex( AtlasMeshVertexIndex v ) const { return SimpleTriangleMesh< GeometryReal , 2 >::vertices[ static_cast< unsigned int >(v) ]; }
		SimplexIndex< 2 , AtlasMeshVertexIndex > triangle( AtlasMeshTriangleIndex t ) const
		{
			SimplexIndex< 2 > _tri = SimpleTriangleMesh< GeometryReal , 2 >::triangles[ static_cast< unsigned int >(t) ]
			SimplexIndex< 2 , AtlasMeshVertexIndex > tri;
			for( unsigned int k=0 ; k<=2 ; k++ ) tri[k] = AtlasMeshVertexIndex( _tri[k] )
			return tri;
		}
		size_t numVertices( void ) const { return SimpleTriangleMesh< GeometryReal , 2 >::vertices.size(); }
		size_t numTriangles( void ) const { return SimpleTriangleMesh< GeometryReal , 2 >::triangles.size(); }
#else // !NEW_INDEXING
		using SimpleTriangleMesh< GeometryReal , 2 >::vertices;		// The chart coordinates
		using SimpleTriangleMesh< GeometryReal , 2 >::triangles;	// The chart triangles
		using SimpleTriangleMesh< GeometryReal , 2 >::edgeIndex;
#endif // NEW_INDEXING

		// Returns the index of the chart the triangle has been assigned to
#ifdef NEW_INDEXING
		ChartIndex triangleToChart( AtlasMeshTriangleIndex t ) const { return _triangleToChart[ static_cast< unsigned int >(t) ]; }
#else // !NEW_INDEXING
		unsigned int triangleToChart( unsigned int t ) const { return _triangleToChart[t]; }
#endif // NEW_INDEXING

		// Returns the index of the edge associated with the half-edge
#ifdef NEW_INDEXING
		AtlasMeshEdgeIndex halfEdgeToEdge( AtlasMeshHalfEdgeIndex he ) const { return _halfEdgeToEdge[ static_cast< unsigned int >(he) ]; }
#else // !NEW_INDEXING
		unsigned int halfEdgeToEdge( unsigned int he ) const { return _halfEdgeToEdge[he]; }
#endif // NEW_INDEXING

		// Returns the index of the chart vertex as an atlas vertex
#ifdef NEW_INDEXING
		AtlasMeshVertexIndex chartToAtlasVertex( ChartMeshVertexIndex v ) const { return _chartToAtlasVertex[ static_cast< unsigned int >(v) ]; }
#else // !NEW_INDEXING
		unsigned int chartToAtlasVertex( unsigned int v ) const { return _chartToAtlasVertex[v]; }
#endif // NEW_INDEXING

		unsigned int numCharts( void ) const { return _numCharts; }

		void initialize( const TexturedTriangleMesh< GeometryReal > &inputMesh );

		// Displace vertex positions if they are too close to the axes
		void jitter( unsigned int width , unsigned int height , GeometryReal epsilon=(GeometryReal)1e-6 );

		std::vector< AtlasChart< GeometryReal > > getCharts( const std::vector< bool > &isBoundaryHalfEdge , unsigned int width , unsigned int height ) const;

	protected:
		unsigned int _numCharts;
#ifdef NEW_INDEXING
		std::vector< ChartIndex > _triangleToChart;
		std::vector< AtlasMeshEdgeIndex > _halfEdgeToEdge;
		std::vector< AtlasMeshVertexIndex > _chartToAtlasVertex;
#else // !NEW_INDEXING
		std::vector< unsigned int > _triangleToChart;
		std::vector< unsigned int > _halfEdgeToEdge;
		std::vector< unsigned int > _chartToAtlasVertex;
#endif // NEW_INDEXING
	};

	template< typename GeometryReal >
#ifdef NEW_INDEXING
	struct AtlasChart : protected SimpleTriangleMesh< GeometryReal , 2 >
#else // !NEW_INDEXING
	struct AtlasChart : public SimpleTriangleMesh< GeometryReal , 2 >
#endif // NEW_INDEXING
	{
#ifdef NEW_INDEXING
		Point< GeometryReal , 2 > vertex( ChartMeshVertexIndex v ) const { return SimpleTriangleMesh< GeometryReal , 2 >::vertices[ static_cast< unsigned int >(v) ]; }
		SimplexIndex< 2 , ChartMeshVertexIndex > triangleIndex( ChartMeshTriangleIndex t ) const
		{
			SimplexIndex< 2 > _tri = SimpleTriangleMesh< GeometryReal , 2 >::triangles[ static_cast< unsigned int >(t) ];
			SimplexIndex< 2 , ChartMeshVertexIndex > tri;
			for( unsigned int k=0 ; k<=2 ; k++ ) tri[k] = ChartMeshVertexIndex( _tri[k] );
			return tri;
		}
		SimplexIndex< 1 , ChartMeshVertexIndex > edgeIndex( ChartMeshHalfEdgeIndex he ) const
		{
			SimplexIndex< 1 > _edge = SimpleTriangleMesh< GeometryReal , 2 >::edgeIndex( static_cast< unsigned int >(he) );
			SimplexIndex< 1 , ChartMeshVertexIndex > edge;
			for( unsigned int k=0 ; k<=1 ; k++ ) edge[k] = ChartMeshVertexIndex( _edge[k] );
			return edge;
		}
		size_t numVertices( void ) const { return SimpleTriangleMesh< GeometryReal , 2 >::vertices.size(); }
		size_t numTriangles( void ) const { return SimpleTriangleMesh< GeometryReal , 2 >::triangles.size(); }
#else // !NEW_INDEXING
		using SimpleTriangleMesh< GeometryReal , 2 >::vertices;		// The chart coordinates
		using SimpleTriangleMesh< GeometryReal , 2 >::triangles;	// The chart triangles
#endif // NEW_INDEXING

		Point2D< GeometryReal > minCorner;
		Point2D< GeometryReal > maxCorner;
		Point2D< GeometryReal > gridOrigin;
		unsigned int originCoords[2];

		// The list of half edges on the boundary of the chart
#ifdef NEW_INDEXING
		std::vector< ChartMeshHalfEdgeIndex > boundaryHalfEdges;
#else // !NEW_INDEXING
		std::vector< unsigned int > boundaryHalfEdges;
#endif // NEW_INDEXING

		// Returns the index of the atlas edge associated with the chart half-edge
#ifdef NEW_INDEXING
		AtlasMeshEdgeIndex atlasEdge( ChartMeshHalfEdgeIndex he ) const { return _chartHalfEdgeToAtlasEdge[ static_cast< unsigned int >(he) ]; }
#else // !NEW_INDEXING
		unsigned int atlasEdge( unsigned int he ) const { return _chartHalfEdgeToAtlasEdge[he]; }
#endif // NEW_INDEXING

		// Returns the index of the atlas half-edge associated with the chart half-edge
#ifdef NEW_INDEXING
		AtlasMeshHalfEdgeIndex atlasHalfEdge( ChartMeshHalfEdgeIndex he ) const { return AtlasMeshHalfEdgeIndex( static_cast< unsigned int >( _chartToAtlasTriangle[static_cast< unsigned int >(he)/3] )*3 + (static_cast< unsigned int >(he)%3) ); }
#else // !NEW_INDEXING
		unsigned int atlasHalfEdge( unsigned int he ) const { return _chartToAtlasTriangle[he/3]*3 + (he%3); }
#endif // NEW_INDEXING

		// Returns the index of the chart vertex as an atlas vertex
#ifdef NEW_INDEXING
		AtlasMeshVertexIndex atlasVertex( ChartMeshVertexIndex v ) const { return _chartToAtlasVertex[ static_cast< unsigned int >(v) ]; }
#else // !NEW_INDEXING
		unsigned int atlasVertex( unsigned int v ) const { return _chartToAtlasVertex[v]; }
#endif // NEW_INDEXING

		// Returns the index of the triangle within the atlas
#ifdef NEW_INDEXING
		AtlasMeshTriangleIndex atlasTriangle( ChartMeshTriangleIndex t ) const { return _chartToAtlasTriangle[ static_cast< unsigned int >(t) ]; }
#else // !NEW_INDEXING
		unsigned int atlasTriangle( unsigned int t ) const { return _chartToAtlasTriangle[t]; }
#endif // NEW_INDEXING

		struct AtlasInfo
		{
			// The opposite half-edges (in the atlas mesh)
#ifdef NEW_INDEXING
			std::vector< AtlasMeshHalfEdgeIndex > oppositeHalfEdges;
#else // !NEW_INDEXING
			std::vector< unsigned int > oppositeHalfEdges;
#endif // NEW_INDEXING

			// A map assigning an index to atlas boundary vertices
#ifdef NEW_INDEXING
			std::map< AtlasMeshVertexIndex , AtlasBoundaryNodeIndex > atlasBoundaryVertexToIndex;
#else // !NEW_INDEXING
			std::map< unsigned int , unsigned int > atlasBoundaryVertexToIndex;
#endif // NEW_INDEXING

			// Is the atlas mesh water-tight
			bool isClosed;
		};

		static std::vector< AtlasChart< GeometryReal > > GetCharts
			(
				const TexturedTriangleMesh< GeometryReal > &mesh ,
				unsigned int width ,
				unsigned int height ,
				AtlasInfo &atlasInfo
			);

	protected:
		friend AtlasMesh< GeometryReal >;

#ifdef NEW_INDEXING
		std::vector< AtlasMeshEdgeIndex > _chartHalfEdgeToAtlasEdge;
		std::vector< AtlasMeshVertexIndex > _chartToAtlasVertex;
		std::vector< AtlasMeshTriangleIndex > _chartToAtlasTriangle;
#else // !NEW_INDEXING
		std::vector< unsigned int > _chartHalfEdgeToAtlasEdge;
		std::vector< unsigned int > _chartToAtlasVertex;
		std::vector< unsigned int > _chartToAtlasTriangle;
#endif // NEW_INDEXING
	};

	template< typename GeometryReal >
	struct IndexedVector2D
	{
#ifdef NEW_INDEXING
		IndexedVector2D( Point2D< GeometryReal > p , ChartMeshVertexIndex index , AtlasMeshVertexIndex vertex ) : p(p) , index(index) , vertex(vertex){}
#else // !NEW_INDEXING
		IndexedVector2D( Point2D< GeometryReal > p , unsigned int index , unsigned int vertex ) : p(p) , index(index) , vertex(vertex){}
#endif // NEW_INDEXING
		Point2D< GeometryReal > p;
#ifdef NEW_INDEXING
		ChartMeshVertexIndex index;
		AtlasMeshVertexIndex vertex;
#else // !NEW_INDEXING
		unsigned int index;
		unsigned int vertex;
#endif // NEW_INDEXING
		bool operator < ( const IndexedVector2D &p2 ) const
		{
			for( int i=0 ; i<2 ; i++ )
			{
				if      ( p[i]<p2.p[i] ) return true;
				else if ( p2.p[i]<p[i] ) return false;
				else
				{
					if     ( vertex<p2.vertex ) return true;
					else if( p2.vertex<vertex ) return false;
				}
			}
			return false;
		}
	};

#ifdef PRE_CLIP_TRIANGLES
#else // !PRE_CLIP_TRIANGLES
	template< typename GeometryReal >
	struct EdgeEquation
	{
		EdgeEquation( void ) : _offset(0) {}
		EdgeEquation( Point2D< GeometryReal > v1 , Point2D< GeometryReal > v2 , bool normalize=false )
		{
			_n = v2 - v1;
			_n = Point2D< GeometryReal >( -_n[1] , _n[0] );
			if( normalize ) _n /= Point2D< GeometryReal >::Length( _n );
			_offset = -Point2D< GeometryReal >::Dot( v1 , _n );
		}
		GeometryReal operator()( Point2D< GeometryReal > p ) const { return Point2D< GeometryReal >::Dot( _n , p ) + _offset; }
		bool makePositive( Point2D< GeometryReal > p )
		{
			if( operator()( p )<0 )
			{
				_n = -_n;
				_offset = -_offset;
				return true;
			}
			else return false;
		}
	protected:
		Point2D< GeometryReal > _n;
		GeometryReal _offset;
	};
#endif // PRE_CLIP_TRIANGLES


#include "AtlasMesh.inl"
#include "AtlasCharts.inl"
}
#endif// ATLAS_MESH_INLCUDED