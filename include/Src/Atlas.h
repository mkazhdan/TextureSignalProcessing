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

namespace MishaK
{
	template< typename GeometryReal > struct AtlasChart;

	template< typename GeometryReal >
	struct AtlasMesh : public SimpleTriangleMesh< GeometryReal , 2 >
	{
		using SimpleTriangleMesh< GeometryReal , 2 >::vertices;		// The texture coordinates
		using SimpleTriangleMesh< GeometryReal , 2 >::triangles;	// The triangles

		// Returns the index of the chart the triangle has been assigned to
		unsigned int triangleChart( unsigned int t ) const { return _triangleToChart[t]; }

		// Returns the index of the edge associated with the half-edge
		unsigned int halfEdgeToEdge( unsigned int he ) const { return _halfEdgeToEdge[he]; }

		// Returns the index of the texture vertex as a surface vertex
		unsigned int textureToSurfaceVertex( unsigned int v ) const { return _textureToSurfaceVertex[v]; }

		unsigned int numCharts;

		void initialize( const TexturedTriangleMesh< GeometryReal > &inputMesh );

		// Displace vertex positions if they are too close to the axes
		void jitter( unsigned int width , unsigned int height , GeometryReal epsilon=(GeometryReal)1e-6 );

		std::vector< AtlasChart< GeometryReal > > getCharts( const std::vector< bool > &isBoundaryHalfEdge , unsigned int width , unsigned int height ) const;

	protected:
		std::vector< unsigned int > _triangleToChart;
		std::vector< unsigned int > _halfEdgeToEdge;
		std::vector< unsigned int > _textureToSurfaceVertex;
	};

	template< typename GeometryReal >
	struct AtlasChart : public SimpleTriangleMesh< GeometryReal , 2 >
	{
		using SimpleTriangleMesh< GeometryReal , 2 >::vertices;		// The texture coordinates
		using SimpleTriangleMesh< GeometryReal , 2 >::triangles;	// The triangles

		Point2D< GeometryReal > minCorner;
		Point2D< GeometryReal > maxCorner;
		Point2D< GeometryReal > gridOrigin;
		int originCoords[2];

		// The list of half edges on the boundary of the chart
		std::vector< unsigned int > boundaryHalfEdges;

		// Returns the index of the atlas edge associated with the chart half-edge
		unsigned int atlasEdge( unsigned int he ) const { return _chartHalfEdgeToAtlasEdge[he]; }

		// Returns the index of the atlas half-edge associated with the chart half-edge
		unsigned int atlasHalfEdge( unsigned int he ) const { return _chartToAtlasTriangle[he/3]*3 + (he%3); }

		// Returns the index of the texture vertex as a surface vertex
		unsigned int surfaceVertex( unsigned int v ) const { return _chartToSurfaceVertex[v]; }

		// Returns the index of the triangle within the atlas
		unsigned int atlasTriangle( unsigned int t ) const { return _chartToAtlasTriangle[t]; }

		struct AtlasInfo
		{
			// The opposite half-edges (in the atlas mesh)
			std::vector< unsigned int > oppositeHalfEdges;

			// A map assigning an index to surface boundary verticess
			std::map< unsigned int , unsigned int > surfaceBoundaryVertexToIndex;

			// Is the surface mesh water-tight
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

		std::vector< unsigned int > _chartHalfEdgeToAtlasEdge;
		std::vector< unsigned int > _chartToSurfaceVertex;
		std::vector< unsigned int > _chartToAtlasTriangle;
	};

	template< typename GeometryReal >
	class IndexedVector2D
	{
	public:
		IndexedVector2D( Point2D< GeometryReal > p , int index , int vertex ) : p(p) , index(index) , vertex(vertex){}
		Point2D< GeometryReal > p;
		int index;
		int vertex;
	};

	template< typename GeometryReal >
	class IndexedVector2DComparison
	{
	public:
		bool operator()( const IndexedVector2D< GeometryReal > &p1 , const IndexedVector2D< GeometryReal > &p2 ) const
		{
			for( int i=0 ; i<2 ; i++ )
			{
				if      ( p1.p[i]<p2.p[i] ) return true;
				else if ( p2.p[i]<p1.p[i] ) return false;
				else
				{
					if     ( p1.vertex<p2.vertex ) return true;
					else if( p2.vertex<p1.vertex ) return false;
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