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

#include "Basis.h"

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

template< typename GeometryReal >
struct AtlasIndexedPolygon
{
	std::vector < Point2D < GeometryReal > > vertices;
	std::vector < int > indices;
	std::vector < int > atlasVertexIndices;
	std::vector < int > atlasVertexParentEdge;
	std::vector < int > atlasEdgeIndices;
	size_t size( void ) const { return vertices.size(); }
	Point2D< GeometryReal >& operator [] ( size_t idx ){ return vertices[idx]; }
	const Point2D< GeometryReal >& operator [] ( size_t idx ) const { return vertices[idx]; }
};

template< typename GeometryReal >
struct AtlasIndexedTriangle
{
	int id;
	Point2D< GeometryReal > vertices[3];
	int indices[3];
	int atlasVertexParentEdge[3];
	int atlasVertexIndices[3];
	int atlasEdgeIndices[3];
};

template< typename GeometryReal >
struct IndexedIntersectionPolygon
{
	std::vector< Point2D < GeometryReal > > vertices;
	std::vector< unsigned long long > indices;
	std::vector< int > edgeIndices;
};

template< typename GeometryReal >
struct IndexedIntersectionTriangle
{
	Point2D< GeometryReal > vertices[3];
	unsigned long long indices[3];
	int edgeIndices[3];
};

unsigned long long SetIntersectionKey(const unsigned long i0, const unsigned long i1) {
	return ( ( (static_cast<unsigned long long>(i0) << 32) & 0xFFFFFFFF00000000) | (static_cast<unsigned long long>(i1) & 0x00000000FFFFFFFF));
}

void GetIntersectionKey(unsigned long long key, unsigned long & i0, unsigned long & i1) {
	i1 = static_cast<unsigned long>(key & 0x00000000FFFFFFFF);
	i0 = static_cast<unsigned long>((key >> 32) & 0x00000000FFFFFFFF);
}

template< typename GeometryReal >
struct IntersectionInfo
{
	unsigned long long intersectionKey;
	int intersectionIndex;
	Point2D< GeometryReal > position;
	GeometryReal time;
};

template< typename GeometryReal >
bool IntersectionComparison( const IntersectionInfo< GeometryReal > &i0 , const IntersectionInfo< GeometryReal > &i1 ){ return i0.time < i1.time; };

template< typename GeometryReal >
struct BoundarySegmentInfo
{
	GeometryReal startTime;
	GeometryReal endTime;
	int halfEdge;
};
