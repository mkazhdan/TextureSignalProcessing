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
#pragma once

#include <array>
#include <optional>
#include "Basis.h"
#include "Indices.h"

namespace MishaK
{
#pragma message( "[WARNING] How is this different from an AuxiliaryNode?" )	// It is used to index both the auxiliary nodes and the covered texels
	template< typename GeometryReal , typename IndexType >
	struct NodeInfo
	{
		NodeInfo( void ) : index(-1){}
		NodeInfo( IndexType index , Point2D< GeometryReal > position ) : index(index) , position(position) {}

		IndexType index;
		Point2D< GeometryReal > position;
	};


	// A characterization of the intersection of grid and mesh edges
	// If the mesh/grid index is undefined, it indexes a grid/mesh node.
	struct GridMeshIntersectionKey
	{
		GridMeshIntersectionKey( void ) : _gridIndex(-1) , _meshIndex(-1){}
		GridMeshIntersectionKey( AtlasGridEdgeIndex g , AtlasMeshEdgeIndex m ) : _gridIndex( static_cast< unsigned int >(g) ) , _meshIndex( static_cast< unsigned int >(m) ){}
		std::optional< ChartMeshVertexIndex > chartVertex( void ) const { if( _gridIndex==-1 && _meshIndex!=-1 ) return ChartMeshVertexIndex( _meshIndex ) ; else return std::nullopt; }
		std::optional<    AtlasGridVertexIndex > gridNode( void ) const { if( _meshIndex==-1 && _gridIndex!=-1 ) return    AtlasGridVertexIndex( _gridIndex ) ; else return std::nullopt; }
		std::optional< std::pair< AtlasGridEdgeIndex , AtlasMeshEdgeIndex > > intersection( void ) const { if( _meshIndex!=-1 && _gridIndex!=-1 ) return std::pair< AtlasGridEdgeIndex , AtlasMeshEdgeIndex >( AtlasGridEdgeIndex( _gridIndex ) , AtlasMeshEdgeIndex( _meshIndex ) ) ; else return std::nullopt; }

		bool operator <  ( const GridMeshIntersectionKey &key ) const { return _gridIndex< key._gridIndex || ( _gridIndex==key._gridIndex && _meshIndex<key._meshIndex ); }
		bool operator == ( const GridMeshIntersectionKey &key ) const { return _gridIndex==key._gridIndex && _meshIndex==key._meshIndex; }
		bool operator != ( const GridMeshIntersectionKey &key ) const { return _gridIndex!=key._gridIndex || _meshIndex!=key._meshIndex; }

		static GridMeshIntersectionKey    GridNodeKey(    AtlasGridVertexIndex g ){ return _Init( static_cast< unsigned int >(g) , -1 ); }
		static GridMeshIntersectionKey ChartVertexKey( ChartMeshVertexIndex m ){ return _Init( -1 , static_cast< unsigned int >(m) ); }

		friend std::ostream &operator << ( std::ostream &s ,  const GridMeshIntersectionKey &iKey ){ return s << "( " << iKey._gridIndex << " , " << iKey._meshIndex << " )"; }
	protected:
		static GridMeshIntersectionKey _Init( unsigned int g , unsigned int m )
		{
			GridMeshIntersectionKey key;
			key._gridIndex = g;
			key._meshIndex = m;
			return key;
		}
		unsigned int _gridIndex , _meshIndex;
	};

	template< typename GeometryReal >
	struct BoundaryIndexedTriangle
	{
		unsigned int id;
		Point2D< GeometryReal > vertices[3];
		AtlasMeshEdgeIndex atlasVertexParentEdge[3];
		ChartMeshVertexIndex vertexIndices[3];
		AtlasMeshEdgeIndex atlasEdgeIndices[3];
		QuadraticElementIndex indices;
		Point2D< GeometryReal >& operator [] ( size_t idx ){ return vertices[idx]; }
		const Point2D< GeometryReal >& operator [] ( size_t idx ) const { return vertices[idx]; }
	};

	template< typename GeometryReal , unsigned int N >
	struct _IndexedPolygon
	{
		template< typename T > using Array = std::conditional_t< N==-1 , std::vector< T > , std::array< T , N > >;
		Array< Point2D < GeometryReal > > vertices;			// The positions of the vertices within the chart
		Array< AtlasInteriorOrBoundaryNodeIndex > indices;	// The index of the boundary vertex
		Array< ChartMeshVertexIndex > vertexIndices;		// The index of the atlas/chart vertex (or -1 if it is not an atlas vertex)
		Array< AtlasMeshEdgeIndex > atlasVertexParentEdge;	// If this is a not an original mesh vertex, the index of the associated polygon edge 
		Array< AtlasMeshEdgeIndex > atlasEdgeIndices;		// If this is a boundary segment, the index of the associated polygon edge
		size_t size( void ) const { return vertices.size(); }
		Point2D< GeometryReal >& operator [] ( size_t idx ){ return vertices[idx]; }
		const Point2D< GeometryReal >& operator [] ( size_t idx ) const { return vertices[idx]; }
	};
	template< typename GeometryReal > using IndexedPolygon  = _IndexedPolygon< GeometryReal , static_cast< unsigned int >( -1 ) >;
	template< typename GeometryReal > using IndexedTriangle = _IndexedPolygon< GeometryReal , 3 >;


	template< typename GeometryReal , unsigned int N >
	struct _IndexedIntersectionPolygon
	{
		template< typename T > using Array = std::conditional_t< N==-1 , std::vector< T > , std::array< T , N > >;
		Array< Point2D< GeometryReal > > vertices;
		Array< GridMeshIntersectionKey > cornerKeys;
		Array< AtlasGridOrMeshEdgeIndex > outgoingEdgeIndices; // Could be an AtlasMeshEdgeIndex or a AtlasGridEdgeIndex
	};
	template< typename GeometryReal > using IndexedIntersectionPolygon  = _IndexedIntersectionPolygon< GeometryReal , static_cast< unsigned int >( -1 ) >;
	template< typename GeometryReal > using IndexedIntersectionTriangle = _IndexedIntersectionPolygon< GeometryReal , 3 >;


	template< typename GeometryReal >
	struct IntersectionInfo
	{
		GridMeshIntersectionKey intersectionKey;
		Point2D< GeometryReal > position;
		GeometryReal time;
		AtlasRefinedBoundaryVertexIndex index;

		IntersectionInfo( void ) : index(-1){}
		IntersectionInfo( GridMeshIntersectionKey intersectionKey , Point2D< GeometryReal > position , GeometryReal time , AtlasRefinedBoundaryVertexIndex index=AtlasRefinedBoundaryVertexIndex(-1) )
			: intersectionKey(intersectionKey) , position(position) , time(time) , index(index){}

		static bool CompareByTime( const IntersectionInfo< GeometryReal > &i0 , const IntersectionInfo< GeometryReal > &i1 ){ return i0.time < i1.time; };
	};

	template< typename GeometryReal >
	struct BoundarySegmentInfo
	{
		GeometryReal startTime;
		GeometryReal endTime;
		ChartMeshHalfEdgeIndex chartHalfEdge;

		BoundarySegmentInfo( GeometryReal startTime=0 , GeometryReal endTime=1 , ChartMeshHalfEdgeIndex chartHalfEdge=ChartMeshHalfEdgeIndex(-1) )
			: startTime(startTime) , endTime(endTime) , chartHalfEdge(chartHalfEdge){}
	};
}