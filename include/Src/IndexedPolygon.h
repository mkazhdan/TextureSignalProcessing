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
#include "Basis.h"

namespace MishaK
{
	template< typename GeometryReal >
	struct NodeInfo
	{
		NodeInfo( void ) : index(-1){}
		NodeInfo( unsigned int index , Point2D< GeometryReal > position ) : index(index) , position(position) {}

		unsigned int index;
		Point2D< GeometryReal > position;
	};

	struct GridMeshIntersectionKey
	{
		unsigned int gridIndex , meshIndex;
		GridMeshIntersectionKey( void ) : gridIndex(-1) , meshIndex(-1){}
		GridMeshIntersectionKey( unsigned int g , unsigned int m ) : gridIndex(g) , meshIndex(m){}
		bool isMeshVertex  ( void ) const { return gridIndex==-1 && meshIndex!=-1; }
		bool isGridVertex  ( void ) const { return meshIndex==-1 && gridIndex!=-1; }
		bool isIntersection( void ) const { return meshIndex!=-1 && gridIndex!=-1; }

		bool operator < ( const GridMeshIntersectionKey &key ) const { return gridIndex<key.gridIndex || ( gridIndex==key.gridIndex && meshIndex<key.meshIndex ); }
		bool operator == ( const GridMeshIntersectionKey &key ) const { return gridIndex==key.gridIndex && meshIndex==key.meshIndex; }
		bool operator != ( const GridMeshIntersectionKey &key ) const { return gridIndex!=key.gridIndex || meshIndex!=key.meshIndex; }

		static GridMeshIntersectionKey GridKey( unsigned int g ){ return GridMeshIntersectionKey( g , -1 ); }
		static GridMeshIntersectionKey MeshKey( unsigned int m ){ return GridMeshIntersectionKey( -1 , m ); }

		friend std::ostream &operator << ( std::ostream &s ,  const GridMeshIntersectionKey &iKey ){ return s << "( " << iKey.gridIndex << " , " << iKey.meshIndex << " )"; }
	};

	template< typename GeometryReal >
	struct BoundaryIndexedTriangle
	{
		int id;
		Point2D< GeometryReal > vertices[3];
		int atlasVertexParentEdge[3];
		int atlasVertexIndices[3];
		int atlasEdgeIndices[3];
		QuadraticElementIndex indices;
		Point2D< GeometryReal >& operator [] ( size_t idx ){ return vertices[idx]; }
		const Point2D< GeometryReal >& operator [] ( size_t idx ) const { return vertices[idx]; }
	};

	template< typename GeometryReal , unsigned int N >
	struct _AtlasIndexedPolygon
	{
		template< typename T > using Array = std::conditional_t< N==-1 , std::vector< T > , std::array< T , N > >;
		Array< Point2D < GeometryReal > > vertices;
		Array< int > indices;
		Array< int > atlasVertexIndices;
		Array< int > atlasVertexParentEdge;
		Array< int > atlasEdgeIndices;
		size_t size( void ) const { return vertices.size(); }
		Point2D< GeometryReal >& operator [] ( size_t idx ){ return vertices[idx]; }
		const Point2D< GeometryReal >& operator [] ( size_t idx ) const { return vertices[idx]; }
	};
	template< typename GeometryReal > using AtlasIndexedPolygon  = _AtlasIndexedPolygon< GeometryReal , static_cast< unsigned int >( -1 ) >;
	template< typename GeometryReal > using AtlasIndexedTriangle = _AtlasIndexedPolygon< GeometryReal , 3 >;

	template< typename GeometryReal , unsigned int N >
	struct _IndexedIntersectionPolygon
	{
		template< typename T > using Array = std::conditional_t< N==-1 , std::vector< T > , std::array< T , N > >;
		Array< Point2D< GeometryReal > > vertices;
		Array< GridMeshIntersectionKey > indices;
		Array< int > edgeIndices;	// [???]
	};
	template< typename GeometryReal > using IndexedIntersectionPolygon  = _IndexedIntersectionPolygon< GeometryReal , static_cast< unsigned int >( -1 ) >;
	template< typename GeometryReal > using IndexedIntersectionTriangle = _IndexedIntersectionPolygon< GeometryReal , 3 >;

	template< typename GeometryReal >
	struct IntersectionInfo
	{
		GridMeshIntersectionKey intersectionKey;
		int intersectionIndex;
		Point2D< GeometryReal > position;
		GeometryReal time;

		static bool CompareByTime( const IntersectionInfo< GeometryReal > &i0 , const IntersectionInfo< GeometryReal > &i1 ){ return i0.time < i1.time; };
	};

	template< typename GeometryReal >
	struct BoundarySegmentInfo
	{
		GeometryReal startTime;
		GeometryReal endTime;
		unsigned int chartHalfEdge;
	};
}