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
#include "Indices.h"

namespace MishaK
{
	template< typename GeometryReal > struct AtlasChart;

	template< typename GeometryReal >
	struct AtlasMesh : protected SimpleTriangleMesh< GeometryReal , 2 >
	{
		Point< GeometryReal , 2 > vertex( AtlasMeshVertexIndex v ) const { return SimpleTriangleMesh< GeometryReal , 2 >::vertices[ static_cast< unsigned int >(v) ]; }
		SimplexIndex< 2 , AtlasMeshVertexIndex > triangle( AtlasMeshTriangleIndex t ) const
		{
			SimplexIndex< 2 > _tri = SimpleTriangleMesh< GeometryReal , 2 >::triangles[ static_cast< unsigned int >(t) ];
			SimplexIndex< 2 , AtlasMeshVertexIndex > tri;
			for( unsigned int k=0 ; k<=2 ; k++ ) tri[k] = AtlasMeshVertexIndex( _tri[k] );
			return tri;
		}
		size_t numVertices( void ) const { return SimpleTriangleMesh< GeometryReal , 2 >::vertices.size(); }
		size_t numTriangles( void ) const { return SimpleTriangleMesh< GeometryReal , 2 >::triangles.size(); }

		// Returns the index of the chart the triangle has been assigned to
		ChartIndex triangleToChart( AtlasMeshTriangleIndex t ) const { return _triangleToChart[t]; }

		// Returns the index of the edge associated with the half-edge
		AtlasMeshEdgeIndex halfEdgeToEdge( AtlasMeshHalfEdgeIndex he ) const { return _halfEdgeToEdge[he]; }

		// Returns the index of the chart vertex as an atlas vertex
		AtlasMeshVertexIndex chartToAtlasVertex( ChartMeshVertexIndex v ) const { return _chartToAtlasVertex[v]; }

		unsigned int numCharts( void ) const { return _numCharts; }

		void initialize( const TexturedTriangleMesh< GeometryReal > &inputMesh );

		// Displace vertex positions if they are too close to the axes
		void jitter( unsigned int width , unsigned int height , GeometryReal epsilon=(GeometryReal)1e-6 );

		IndexVector< ChartIndex , AtlasChart< GeometryReal > > getCharts( const std::vector< bool > &isBoundaryHalfEdge , unsigned int width , unsigned int height ) const;

	protected:
		unsigned int _numCharts;
		IndexVector< AtlasMeshTriangleIndex , ChartIndex > _triangleToChart;
		IndexVector< AtlasMeshHalfEdgeIndex , AtlasMeshEdgeIndex > _halfEdgeToEdge;
		IndexVector< ChartMeshVertexIndex , AtlasMeshVertexIndex > _chartToAtlasVertex;
	};

	template< typename GeometryReal >
	struct AtlasChart : protected SimpleTriangleMesh< GeometryReal , 2 >
	{
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

		Point2D< GeometryReal > minCorner;
		Point2D< GeometryReal > maxCorner;
		Point2D< GeometryReal > gridOrigin;
		unsigned int originCoords[2];

		// The list of half edges on the boundary of the chart
		std::vector< ChartMeshHalfEdgeIndex > boundaryHalfEdges;

		// Returns the index of the atlas edge associated with the chart half-edge
		AtlasMeshEdgeIndex atlasEdge( ChartMeshHalfEdgeIndex he ) const { return _chartHalfEdgeToAtlasEdge[he]; }

		// Returns the index of the atlas half-edge associated with the chart half-edge
		AtlasMeshHalfEdgeIndex atlasHalfEdge( ChartMeshHalfEdgeIndex he ) const { auto f = FactorChartMeshHalfEdgeIndex( he ) ; return GetAtlasMeshHalfEdgeIndex( _chartToAtlasTriangle[f.first] , f.second ); }

		// Returns the index of the chart vertex as an atlas vertex
		AtlasMeshVertexIndex atlasVertex( ChartMeshVertexIndex v ) const { return _chartToAtlasVertex[v]; }

		// Returns the index of the triangle within the atlas
		AtlasMeshTriangleIndex atlasTriangle( ChartMeshTriangleIndex t ) const { return _chartToAtlasTriangle[t]; }

		struct AtlasInfo
		{
			// The opposite half-edges (in the atlas mesh)
			IndexVector< AtlasMeshHalfEdgeIndex , AtlasMeshHalfEdgeIndex > oppositeHalfEdges;

			// A map assigning an index to atlas boundary vertices
			std::map< AtlasMeshVertexIndex , AtlasInteriorOrBoundaryNodeIndex > atlasBoundaryVertexToIndex;

			// Is the atlas mesh water-tight
			bool isClosed;
		};

		static IndexVector< ChartIndex , AtlasChart< GeometryReal > > GetCharts
			(
				const TexturedTriangleMesh< GeometryReal > &mesh ,
				unsigned int width ,
				unsigned int height ,
				AtlasInfo &atlasInfo
			);

	protected:
		friend AtlasMesh< GeometryReal >;

		IndexVector< ChartMeshHalfEdgeIndex , AtlasMeshEdgeIndex > _chartHalfEdgeToAtlasEdge;
		IndexVector< ChartMeshVertexIndex , AtlasMeshVertexIndex > _chartToAtlasVertex;
		IndexVector< ChartMeshTriangleIndex , AtlasMeshTriangleIndex > _chartToAtlasTriangle;
	};

	template< typename GeometryReal >
	struct IndexedVector2D
	{
		IndexedVector2D( Point2D< GeometryReal > p , ChartMeshVertexIndex index , AtlasMeshVertexIndex vertex ) : p(p) , index(index) , vertex(vertex){}
		Point2D< GeometryReal > p;
		ChartMeshVertexIndex index;
		AtlasMeshVertexIndex vertex;
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