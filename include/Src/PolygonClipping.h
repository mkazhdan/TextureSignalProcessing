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
#ifndef POLYGON_CLIPPING_INCLUDED
#define POLYGON_CLIPPING_INCLUDED

#include <Misha/Miscellany.h>
#include <Misha/Geometry.h>
#include "IndexedPolygon.h"

template< typename GeometryReal >
struct CellClippedTriangle
{
	unsigned int sz;
	Point2D< GeometryReal > vertices[7];

	CellClippedTriangle( void ) : sz(0){}
	CellClippedTriangle( const std::vector< Point2D< GeometryReal > > &triangle )
	{
		sz = (int)triangle.size();
		for( int i=0 ; i<triangle.size() ; i++ ) vertices[i] = triangle[i];
	}
	int size( void ) const { return (int)sz; }
	void push_back( Point2D< GeometryReal > p ){ vertices[sz++] = p; }
	Point2D< GeometryReal >& operator[] ( unsigned int idx ) { return vertices[idx]; }
	const Point2D< GeometryReal >& operator[] ( unsigned int idx ) const { return vertices[idx]; }
};

template< typename GeometryReal >
void ClipConvexPolygon( CellClippedTriangle< GeometryReal > &vertices , const Point2D< GeometryReal > &normal , const GeometryReal &bOff , const GeometryReal &fOff )
{
	enum
	{
		BACK=-1 ,
		INSIDE ,
		FRONT
	};
	int vCount = (int)vertices.size();


	CellClippedTriangle< GeometryReal > _vertices;

	Point2D< GeometryReal > pVertex = vertices[vCount - 1];
	GeometryReal pDot = Point2D< GeometryReal >::Dot( pVertex , normal );
	auto GetLabel = [&]( GeometryReal dot )
	{
		if     ( dot< bOff ) return BACK;
		else if( dot>=fOff ) return FRONT;
		else                 return INSIDE;
	};
	int pLabel = GetLabel( pDot );
	for( int i=0 ; i<vCount ; i++ )
	{
		Point2D< GeometryReal > cVertex = vertices[i];
		GeometryReal cDot = Point2D< GeometryReal >::Dot( cVertex , normal );
		int cLabel = GetLabel( cDot );
		GeometryReal bAlpha , fAlpha;
		switch( cLabel+pLabel )
		{
		case -1:
			bAlpha = (cDot-bOff)/(cDot-pDot);
			_vertices.push_back( pVertex*bAlpha + cVertex*( (GeometryReal)1.-bAlpha ) );
			break;
		case 1:
			fAlpha = (cDot-fOff)/(cDot-pDot);
			_vertices.push_back( pVertex*fAlpha + cVertex*( (GeometryReal)1.-fAlpha ) );
			break;
		case 0:
			if( pLabel==BACK )
			{
				bAlpha = (cDot-bOff)/(cDot-pDot);
				fAlpha = (cDot-fOff)/(cDot-pDot);
				_vertices.push_back( pVertex*bAlpha + cVertex*( (GeometryReal)1.-bAlpha ) );
				_vertices.push_back( pVertex*fAlpha + cVertex*( (GeometryReal)1.-fAlpha ) );
			}
			else if( pLabel==FRONT )
			{
				fAlpha = (cDot-fOff)/(cDot-pDot);
				bAlpha = (cDot-bOff)/(cDot-pDot);
				_vertices.push_back( pVertex*fAlpha + cVertex*( (GeometryReal)1.-fAlpha ) );
				_vertices.push_back( pVertex*bAlpha + cVertex*( (GeometryReal)1.-bAlpha ) );
			}
		}
		if( cLabel==INSIDE ) _vertices.push_back( cVertex );
		pVertex = cVertex , pDot = cDot , pLabel = cLabel;
	}

	vertices = _vertices;
}
template< typename GeometryReal >
void ClipConvexPolygon( std::vector< Point2D< GeometryReal > > &vertices , const Point2D< GeometryReal > &normal , const GeometryReal &offset )
{
	int verticesCount = (int)vertices.size();

	std::vector< Point2D< GeometryReal > > outputVertices;
	outputVertices.reserve( verticesCount+1 );

	Point2D< GeometryReal > previousVertex = vertices[ verticesCount-1 ];
	GeometryReal previousLevel = Point2D< GeometryReal >::Dot( previousVertex , normal ) - offset;
	bool isPreviousInterior =  previousLevel > 0;
	for( int i=0 ; i<vertices.size() ; i++ )
	{
		Point2D< GeometryReal > currentVertex = vertices[i];
		GeometryReal currentLevel = Point2D< GeometryReal >::Dot( currentVertex , normal ) - offset;
		bool isInterior = currentLevel > 0;
		if( isInterior!=isPreviousInterior )
		{
			GeometryReal alpha = currentLevel / (currentLevel - previousLevel);
			Point2D< GeometryReal > intersection = previousVertex*alpha + currentVertex*(GeometryReal)( 1.-alpha );
			outputVertices.push_back( intersection );
		}
		if( isInterior ) outputVertices.push_back( currentVertex );
		previousVertex = currentVertex;
		previousLevel = currentLevel;
		isPreviousInterior = isInterior;
	}

	vertices = outputVertices;
}

template< typename GeometryReal >
int ClipTriangleToPrimalCell( CellClippedTriangle< GeometryReal >&tri , int i , int j , GeometryReal cellSizeW , GeometryReal cellSizeH )
{
	ClipConvexPolygon( tri , Point2D< GeometryReal >(0,1) , cellSizeH*j , cellSizeH*(j+1) );
	ClipConvexPolygon( tri , Point2D< GeometryReal >(1,0) , cellSizeW*i , cellSizeW*(i+1) );
	return (int)tri.size();
}

//Vertex type 
// -2 non initialized
// -1 exterior
// 0 on edge
// 1 interior

template< typename GeometryReal >
void ClipPartiallyIndexedPolygonToIndexedEdge( AtlasIndexedPolygon< GeometryReal > &polygon , const Point2D< GeometryReal > &normal , GeometryReal &offset , int edgeIndex, int atlasVertexIndices[2] )
{
	std::vector< Point2D< GeometryReal > > outputVertices;
	std::vector< int > outputVertexIndices; 
	std::vector< int > outputEdgeIndices;
	std::vector< int > outputParentVertexEdgeIndices;
	outputVertices.reserve( polygon.vertices.size()+3 );
	outputVertexIndices.reserve( polygon.vertices.size()+3 );
	outputEdgeIndices.reserve( polygon.vertices.size()+3 );
	outputParentVertexEdgeIndices.reserve( polygon.vertices.size()+3 );

	Point2D< GeometryReal > previousVertex = polygon.vertices[polygon.vertices.size() - 1];
	int previousVertexIndex = polygon.atlasVertexIndices[polygon.vertices.size() - 1];
	int previousEdgeIndex = polygon.atlasEdgeIndices[polygon.vertices.size() - 2];
	int nextEdgeIndex = polygon.atlasEdgeIndices[polygon.vertices.size() - 1];
	int previousVertexEdgeSupport = polygon.atlasVertexParentEdge[polygon.vertices.size() - 1];

	int previousVertexType = -2;
	GeometryReal previousLevel;

	bool emptyPolygon = true;

	bool isCurrentVertexACorner = previousVertexIndex != -1 && (previousVertexIndex == atlasVertexIndices[0] || previousVertexIndex == atlasVertexIndices[1]);
	bool isOnTheEdge = previousVertexEdgeSupport != -1 && previousVertexEdgeSupport == edgeIndex;
	if (isOnTheEdge || isCurrentVertexACorner)
	{
		previousVertexType = 0;
		previousLevel = 0;
	}
	else
	{
		previousLevel = Point2D< GeometryReal >::Dot( previousVertex , normal ) - offset;
		previousVertexType = previousLevel > 0 ? 1 : -1;
	}

	for( int i=0 ; i<polygon.vertices.size() ; i++ )
	{
		Point2D< GeometryReal > currentVertex = polygon.vertices[i];
		int currentVertexIndex = polygon.atlasVertexIndices[i];
		int currentVertexType = -2;
		int currentVertexEdgeSupport = polygon.atlasVertexParentEdge[i];
		GeometryReal currentLevel;

		isCurrentVertexACorner = currentVertexIndex != -1 && (currentVertexIndex == atlasVertexIndices[0] || currentVertexIndex == atlasVertexIndices[1]);
		isOnTheEdge = currentVertexEdgeSupport != -1 && currentVertexEdgeSupport == edgeIndex;
		//isPreviousEdgeColinear = isNextEdgeColinear;
		previousEdgeIndex = nextEdgeIndex;

		nextEdgeIndex = polygon.atlasEdgeIndices[i];
		//isNextEdgeColinear = nextEdgeIndex != -1 && nextEdgeIndex == edgeIndex;

		if (isOnTheEdge || isCurrentVertexACorner)
		{
			currentVertexType = 0;
			currentLevel = 0;
		}
		else
		{
			currentLevel = Point2D< GeometryReal >::Dot( currentVertex , normal ) - offset;
			currentVertexType = currentLevel > 0 ? 1 : -1;
		}

		if (previousVertexType == -1){
			if (currentVertexType == -1){
				//Do nothing
			}
			else if (currentVertexType == -0){
				//Do nothing
			}
			else //Entrying edge
			{
				GeometryReal alpha = -currentLevel / (previousLevel - currentLevel);
				Point2D< GeometryReal > intersection = previousVertex*alpha + currentVertex*(GeometryReal)( 1.-alpha );
				outputVertices.push_back(intersection);
				outputVertexIndices.push_back(-1);
				outputEdgeIndices.push_back(previousEdgeIndex);
				outputParentVertexEdgeIndices.push_back(edgeIndex);
				//lastAddedVertexType = 0;
			}
		}
		else if (previousVertexType == 0){
			outputVertices.push_back(previousVertex);
			outputVertexIndices.push_back(previousVertexIndex);
			outputParentVertexEdgeIndices.push_back(previousVertexEdgeSupport);
			//lastAddedVertexType = 0;
			if (currentVertexType < 0){
				outputEdgeIndices.push_back(edgeIndex);
			}
			else if (currentVertexType == 0){
				outputEdgeIndices.push_back(edgeIndex);
			}
			else{
				outputEdgeIndices.push_back(previousEdgeIndex);
			}
		}
		else{
			outputVertices.push_back(previousVertex);
			outputVertexIndices.push_back(previousVertexIndex);
			outputParentVertexEdgeIndices.push_back(previousVertexEdgeSupport);
			emptyPolygon = false;
			//lastAddedVertexType = 1;
			if( currentVertexType<0 ) //Exiting edge
			{
				outputEdgeIndices.push_back(previousEdgeIndex);

				GeometryReal alpha = -currentLevel / (previousLevel - currentLevel);
				Point2D< GeometryReal > intersection = previousVertex*alpha + currentVertex*(GeometryReal)( 1.-alpha );
				outputVertices.push_back(intersection);
				outputVertexIndices.push_back(-1);
				outputEdgeIndices.push_back(edgeIndex);
				outputParentVertexEdgeIndices.push_back(edgeIndex);
				//lastAddedVertexType = 0;
			}
			else if (currentVertexType == 0){
				outputEdgeIndices.push_back(previousEdgeIndex);
			}
			else{
				outputEdgeIndices.push_back(previousEdgeIndex);
			}
		}

		previousVertex = currentVertex;
		previousLevel = currentLevel;
		previousVertexType = currentVertexType;
		previousVertexIndex = currentVertexIndex;
		previousVertexEdgeSupport = currentVertexEdgeSupport;
	}

	if (emptyPolygon){
		polygon.vertices.clear();
		polygon.atlasVertexIndices.clear();
		polygon.atlasEdgeIndices.clear();
		polygon.atlasVertexParentEdge.clear();
	}
	else{
		polygon.vertices = outputVertices;
		polygon.atlasVertexIndices = outputVertexIndices;
		polygon.atlasEdgeIndices = outputEdgeIndices;
		polygon.atlasVertexParentEdge = outputParentVertexEdgeIndices;

		if( polygon.vertices.size()!=polygon.atlasVertexIndices.size() || polygon.vertices.size()!=polygon.atlasEdgeIndices.size() || polygon.vertices.size()!=polygon.atlasVertexParentEdge.size() )
			THROW( "Polygon array size does not match" );

		//Check for non consecutive colinear edges
		for (int i = 0; i < polygon.atlasEdgeIndices.size(); i++){
			if( polygon.atlasEdgeIndices[i]!=-1 && polygon.atlasEdgeIndices[i]==polygon.atlasEdgeIndices[ (i+1)%polygon.atlasEdgeIndices.size() ] )
				THROW( "Unexpected consecutive colinear edges" );
		}
	}
}

//Only for convex polygons
template< typename GeometryReal >
int ClipPartiallyIndexedPolygonToIndexedTriangle(AtlasIndexedPolygon< GeometryReal > &polygon , const AtlasIndexedTriangle< GeometryReal > &triangle , bool verbose=false )
{
	Point2D< GeometryReal > triangleCenter = (triangle.vertices[0] + triangle.vertices[1] + triangle.vertices[2]) / 3;

	for (int k = 0; k < 3; k++){

		//Compute edge normal (pointing inside)
		Point2D< GeometryReal > dir =   triangle.vertices[ (k+1)%3 ] - triangle.vertices[k];
		Point2D< GeometryReal > mid = ( triangle.vertices[ (k+1)%3 ] + triangle.vertices[k] ) / 2;
		Point2D< GeometryReal > normal(-dir[1], dir[0]);
		normal /= (GeometryReal)Point2D< GeometryReal >::Length(normal);
		GeometryReal offset = Point2D< GeometryReal >::Dot( mid , normal );
		if( Point2D< GeometryReal >::Dot( triangleCenter , normal )<offset )
		{
			normal *= -(GeometryReal)1.;
			offset *= -(GeometryReal)1.;
		}

		int edgeIndex = triangle.atlasEdgeIndices[k];
		int atlasVertexIndices[2] = { triangle.atlasVertexIndices[k], triangle.atlasVertexIndices[(k + 1) % 3] };
		ClipPartiallyIndexedPolygonToIndexedEdge( polygon , normal , offset , edgeIndex , atlasVertexIndices );
			
		if( verbose )
		{
			printf("polygon after clipping against edge %d, corners %d %d, normal  %f %f, offset %f. \n", edgeIndex, atlasVertexIndices[0], atlasVertexIndices[1],normal[0],normal[1],offset);
			printf("%d \n", (int)polygon.vertices.size());
			for (int v = 0; v < polygon.vertices.size(); v++) printf("%f %f %f \n", polygon.vertices[v][0], polygon.vertices[v][1], 0.);
			for (int v = 0; v < polygon.vertices.size(); v++) printf("%d ", polygon.atlasVertexIndices[v]);
			printf("\n");
			for (int v = 0; v < polygon.vertices.size(); v++) printf("%d ", polygon.atlasEdgeIndices[v]);
			printf("\n");
		}
	}

	return (int)polygon.vertices.size();
}

template< typename GeometryReal >
void SetAtlasIndexedPolygonFromTriangle( const AtlasIndexedTriangle< GeometryReal > &triangle , AtlasIndexedPolygon< GeometryReal > &polygon )
{
	for( int k=0 ; k<3 ; k++ )
	{
		polygon.vertices.push_back( triangle.vertices[k] );
		polygon.indices.push_back( triangle.indices[k] );
		polygon.atlasVertexIndices.push_back( triangle.atlasVertexIndices[k] );
		polygon.atlasEdgeIndices.push_back( triangle.atlasEdgeIndices[k] );
		polygon.atlasVertexParentEdge.push_back( triangle.atlasVertexParentEdge[k]) ;
	}
}

//Points in general positions
template< typename GeometryReal >
void ClipIndexedIntersectionPolygonToIndexedIntersectionEdge( IndexedIntersectionPolygon< GeometryReal > &polygon , const Point2D< GeometryReal > &normal , GeometryReal &offset , int edgeIndex )
{
	std::vector< Point2D< GeometryReal > > outputVertices;
	std::vector< unsigned long long > outputIndices;
	std::vector< int > outputEdgeIndices;

	Point2D< GeometryReal > previousVertex = polygon.vertices[ polygon.vertices.size()-1 ];
	unsigned long long previousVertexIndex = polygon.indices[polygon.vertices.size() - 1];
	int previousEdgeIndex = polygon.edgeIndices[polygon.vertices.size() - 2];
	int nextEdgeIndex = polygon.edgeIndices[polygon.vertices.size() - 1];
	
	GeometryReal previousLevel = Point2D< GeometryReal >::Dot( previousVertex , normal ) - offset;
	bool isPreviousInterior = previousLevel > 0;

	for (int i = 0; i < polygon.vertices.size(); i++)
	{
		Point2D< GeometryReal > currentVertex = polygon.vertices[i];
		unsigned long long currentVertexIndex = polygon.indices[i];
		previousEdgeIndex = nextEdgeIndex;
		nextEdgeIndex = polygon.edgeIndices[i];
		
		GeometryReal currentLevel = Point2D< GeometryReal >::Dot( currentVertex , normal ) - offset;
		bool isCurrentInterior = currentLevel > 0;
		
		if (!isPreviousInterior){
			if (isCurrentInterior){//Entrying edge
				GeometryReal alpha = -currentLevel / (previousLevel - currentLevel);
				Point2D< GeometryReal > intersection = previousVertex*alpha + currentVertex*(GeometryReal)( 1.-alpha );
				outputVertices.push_back(intersection);

				unsigned long long intersectionVertexKey = SetIntersectionKey(edgeIndex, previousEdgeIndex);
				outputIndices.push_back(intersectionVertexKey);

				outputEdgeIndices.push_back(previousEdgeIndex);
			}
		}
		else{
			outputVertices.push_back(previousVertex);
			outputIndices.push_back(previousVertexIndex);
			outputEdgeIndices.push_back(previousEdgeIndex);

			if (!isCurrentInterior) //Exiting edge
			{
				GeometryReal alpha = -currentLevel / (previousLevel - currentLevel);
				Point2D< GeometryReal > intersection = previousVertex*alpha + currentVertex*(GeometryReal)( 1.-alpha );
				outputVertices.push_back(intersection);
				unsigned long long intersectionVertexKey = SetIntersectionKey(edgeIndex, previousEdgeIndex);
				outputIndices.push_back(intersectionVertexKey);
				outputEdgeIndices.push_back(edgeIndex);
			}
		}
		previousVertex = currentVertex;
		previousLevel = currentLevel;
		previousVertexIndex = currentVertexIndex;
		isPreviousInterior = isCurrentInterior;
	}

	polygon.vertices = outputVertices;
	polygon.indices = outputIndices;
	polygon.edgeIndices = outputEdgeIndices;
}

template< typename GeometryReal >
int ClipIndexedIntersectionPolygonToIndexedIntersectionTriangle( IndexedIntersectionPolygon< GeometryReal > &polygon , const IndexedIntersectionTriangle< GeometryReal > &triangle , bool verbose=false )
{
	Point2D< GeometryReal > triangleCenter = ( triangle.vertices[0] + triangle.vertices[1] + triangle.vertices[2] ) / 3;
	unsigned long long cornerKeys[6];
	for (int k = 0; k < 3; k++){
		cornerKeys[2 * k] = SetIntersectionKey(triangle.edgeIndices[k], triangle.edgeIndices[(k + 2) % 3]);
		cornerKeys[2 * k + 1] = SetIntersectionKey(triangle.edgeIndices[(k + 2) % 3], triangle.edgeIndices[k]);
	}
	unsigned long long cornerIndices[6];
	for (int k = 0; k < 3; k++){
		cornerIndices[2 * k] = cornerIndices[2 * k + 1] = triangle.indices[k];
	}
	
	if (verbose){
		printf("Triangle vertices\n");
		for (int v = 0; v < 3; v++)printf("%f %f %f \n", triangle.vertices[v][0], triangle.vertices[v][1], 0.);
		printf("Triangle corner keys \n");
		printf("%llu %llu %llu\n", cornerKeys[0], cornerKeys[1], cornerKeys[2]);
		printf("Triangle corner indices \n");
		printf("%llu %llu %llu\n", cornerIndices[0], cornerIndices[1], cornerIndices[2]);
	}

	bool reverseOrientation = false;

	for( int k=0 ; k<3 ; k++ )
	{
		//Compute edge normal (pointing inside)
		Point2D< GeometryReal > dir =  triangle.vertices[ (k+1)%3 ] - triangle.vertices[k];
		Point2D< GeometryReal > mid = (triangle.vertices[ (k+1)%3 ] + triangle.vertices[k]) / 2.0;
		Point2D< GeometryReal > normal(-dir[1], dir[0]);
		normal /= Point2D< GeometryReal >::Length( normal );
		GeometryReal offset = Point2D< GeometryReal  >::Dot( mid , normal );
		if( Point2D< GeometryReal >::Dot( triangleCenter , normal )<offset )
		{
			normal *= -(GeometryReal)1.;
			offset *= -(GeometryReal)1.;
			reverseOrientation = true;
		}

		int edgeIndex = triangle.edgeIndices[k];
		ClipIndexedIntersectionPolygonToIndexedIntersectionEdge( polygon , normal , offset , edgeIndex );

		if( verbose )
		{
			printf("polygon after clipping against edge %d, normal  %f %f, offset %f. \n",edgeIndex,normal[0], normal[1], offset);
			printf("%d \n", (int)polygon.vertices.size());
			for (int v = 0; v < polygon.vertices.size(); v++)printf("%f %f %f \n", polygon.vertices[v][0], polygon.vertices[v][1], 0.);
			printf(" Vertex Indices \n");
			for (int v = 0; v < polygon.vertices.size(); v++)printf("%llu ", polygon.indices[v]);
			printf("\n");
			printf(" Edge Indices \n");
			for (int v = 0; v < polygon.vertices.size(); v++) printf("%d ", polygon.edgeIndices[v]);
			printf("\n");
		}
	}

	//Identify triangle corners
	for (int i = 0; i < polygon.indices.size(); i++) for (int k = 0; k < 6; k++) if (polygon.indices[i] == cornerKeys[k]) polygon.indices[i] = cornerIndices[k];

	if (reverseOrientation){
		int n = (int)polygon.vertices.size();
		std::vector< Point2D < GeometryReal > > reversedVertices(n);
		std::vector < unsigned long long > reversedIndices(n);
		std::vector < int > reversedEdges(n);
		for (int k = 0; k < n; k++){
			reversedVertices[k] = polygon.vertices[ n-1-k ];
			reversedIndices[k] = polygon.indices[ n-1-k ];
			reversedEdges[k] = polygon.edgeIndices[ (2*n-k-2)%n ];
		}
		polygon.vertices = reversedVertices;
		polygon.indices = reversedIndices;
		polygon.edgeIndices = reversedEdges;
	}

	return (int)polygon.vertices.size();
}

template< typename GeometryReal >
void GetTriangleIntegerBBox( Point2D< GeometryReal > tPos[3] , const GeometryReal invCellSizeW , const GeometryReal invCellSizeH,  int minCorner[2], int maxCorner[2] )
{
	GeometryReal fminx = std::min< GeometryReal >( std::min< GeometryReal >( tPos[0][0] , tPos[1][0] ) , tPos[2][0]) ;
	fminx = std::max< GeometryReal >( fminx , (GeometryReal)0. );
	GeometryReal fminy = std::min< GeometryReal >( std::min< GeometryReal >( tPos[0][1] , tPos[1][1] ) , tPos[2][1] );
	fminy = std::max< GeometryReal >( fminy , (GeometryReal)0. );
	GeometryReal fmaxx = std::max< GeometryReal >( std::max< GeometryReal >( tPos[0][0] , tPos[1][0] ) , tPos[2][0] );
	fmaxx = std::min< GeometryReal >( fmaxx , (GeometryReal)1. );
	GeometryReal fmaxy = std::max< GeometryReal >( std::max< GeometryReal >( tPos[0][1] , tPos[1][1] ) , tPos[2][1] );
	fmaxy = std::min< GeometryReal >( fmaxy , (GeometryReal)1. );

	minCorner[0] = static_cast< int >( floor( fminx*invCellSizeW ) );
	minCorner[1] = static_cast< int >( floor( fminy*invCellSizeH ) );
	maxCorner[0] = static_cast< int >( ceil ( fmaxx*invCellSizeW ) );
	maxCorner[1] = static_cast< int >( ceil ( fmaxy*invCellSizeH ) );
}

template< typename GeometryReal >
void GetEdgeIntegerBBox( Point2D< GeometryReal > tPos[3] , const GeometryReal invCellSizeW , const GeometryReal invCellSizeH , int minCorner[2] , int maxCorner[2] )
{
	GeometryReal fminx = std::min< GeometryReal >( tPos[0][0] , tPos[1][0] );
	fminx = std::max< GeometryReal >( fminx , (GeometryReal)0. );
	GeometryReal fminy = std::min< GeometryReal >( tPos[0][1] , tPos[1][1] );
	fminy = std::max< GeometryReal >( fminy , (GeometryReal)0. );
	GeometryReal fmaxx = std::max< GeometryReal >( tPos[0][0] , tPos[1][0] );
	fmaxx = std::min< GeometryReal >( fmaxx , (GeometryReal)1. );
	GeometryReal fmaxy = std::max< GeometryReal >( tPos[0][1] , tPos[1][1] );
	fmaxy = std::min< GeometryReal >( fmaxy , (GeometryReal)1. );

	minCorner[0] = static_cast< int >( floor( fminx*invCellSizeW ) );
	minCorner[1] = static_cast< int >( floor( fminy*invCellSizeH ) );
	maxCorner[0] = static_cast< int >( ceil ( fmaxx*invCellSizeW ) );
	maxCorner[1] = static_cast< int >( ceil ( fmaxy*invCellSizeH ) );
}

template< typename GeometryReal >
SquareMatrix< GeometryReal, 2 > GetBarycentricMap( Point2D< GeometryReal > tPos[3] )
{
	SquareMatrix< GeometryReal , 2 > parametrizationMap;
	parametrizationMap.coords[0][0] = tPos[1][0] - tPos[0][0];
	parametrizationMap.coords[0][1] = tPos[1][1] - tPos[0][1];
	parametrizationMap.coords[1][0] = tPos[2][0] - tPos[0][0];
	parametrizationMap.coords[1][1] = tPos[2][1] - tPos[0][1];
	return parametrizationMap.inverse();
}

#endif // POLYGON_CLIPPING_INCLUDED

