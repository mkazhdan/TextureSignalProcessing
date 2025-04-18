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

namespace MishaK
{
#ifdef PRE_CLIP_TRIANGLES
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

	template< typename GeometryReal >
	struct CellClippedTriangle
	{
		unsigned int sz;
		Point2D< GeometryReal > vertices[7];

		CellClippedTriangle( void ) : sz(0){}
		CellClippedTriangle( const std::vector< Point2D< GeometryReal > > &triangle )
		{
			sz = (unsigned int)triangle.size();
			for( unsigned int i=0 ; i<triangle.size() ; i++ ) vertices[i] = triangle[i];
		}
		CellClippedTriangle( const Simplex< GeometryReal , 2 , 2 > &triangle ) : sz(3) { for( unsigned int i=0 ; i<3 ; i++ ) vertices[i] = triangle[i]; }

		unsigned int size( void ) const { return sz; }
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
	unsigned int ClipTriangleToPrimalCell( CellClippedTriangle< GeometryReal >&tri , int i , int j , GeometryReal cellSizeW , GeometryReal cellSizeH )
	{
		ClipConvexPolygon( tri , Point2D< GeometryReal >(0,1) , cellSizeH*j , cellSizeH*(j+1) );
		ClipConvexPolygon( tri , Point2D< GeometryReal >(1,0) , cellSizeW*i , cellSizeW*(i+1) );
		return (unsigned int)tri.size();
	}

	//Vertex type 
	// -2 non initialized
	// -1 exterior
	// 0 on edge
	// 1 interior

	template< typename GeometryReal >
	void ClipPartiallyIndexedPolygonToIndexedEdge
	(
		AtlasIndexedPolygon< GeometryReal > &polygon ,
		const EdgeEquation< GeometryReal > &edgeEquation , 
		unsigned int edgeIndex,
		unsigned int atlasVertexIndices[2]
	)
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
			previousLevel = edgeEquation( previousVertex );
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
				currentLevel = edgeEquation( currentVertex );
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
				MK_THROW( "Polygon array size does not match" );

			//Check for non consecutive colinear edges
			for (int i = 0; i < polygon.atlasEdgeIndices.size(); i++){
				if( polygon.atlasEdgeIndices[i]!=-1 && polygon.atlasEdgeIndices[i]==polygon.atlasEdgeIndices[ (i+1)%polygon.atlasEdgeIndices.size() ] )
					MK_THROW( "Unexpected consecutive colinear edges" );
			}
		}
	}

	//Only for convex polygons
	template< typename GeometryReal >
	unsigned int ClipPartiallyIndexedPolygonToIndexedTriangle( AtlasIndexedPolygon< GeometryReal > &polygon , const AtlasIndexedTriangle< GeometryReal > &triangle )
	{
		Point2D< GeometryReal > triangleCenter = (triangle.vertices[0] + triangle.vertices[1] + triangle.vertices[2]) / 3;

		for( int k=0 ; k<3 ; k++ )
		{
			EdgeIndex eIndex = CornerEdgeIndex( k );

			EdgeEquation< GeometryReal > edgeEquation( triangle.vertices[ eIndex[0] ] , triangle.vertices[ eIndex[1] ] );
			edgeEquation.makePositive( triangleCenter );

			unsigned int edgeIndex = triangle.atlasEdgeIndices[k];
			unsigned int atlasVertexIndices[2] = { triangle.atlasVertexIndices[ eIndex[0] ], triangle.atlasVertexIndices[ eIndex[1] ] };
			ClipPartiallyIndexedPolygonToIndexedEdge( polygon , edgeEquation , edgeIndex , atlasVertexIndices );
		}

		return (unsigned int)polygon.vertices.size();
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

	// Points in general positions
	template< typename GeometryReal >
	void ClipIndexedIntersectionPolygonToIndexedIntersectionEdge
	(
		IndexedIntersectionPolygon< GeometryReal > &polygon ,
		const EdgeEquation< GeometryReal > &edgeEquation , 
		int edgeIndex
	)
	{
		std::vector< Point2D< GeometryReal > > outputVertices;
		std::vector< GridMeshIntersectionKey > outputIndices;
		std::vector< int > outputEdgeIndices;

		Point2D< GeometryReal > previousVertex = polygon.vertices[ polygon.vertices.size()-1 ];
		GridMeshIntersectionKey previousVertexIndex = polygon.indices[polygon.vertices.size() - 1];
		int previousEdgeIndex = polygon.edgeIndices[ polygon.vertices.size()-2 ];
		int nextEdgeIndex = polygon.edgeIndices[ polygon.vertices.size()-1 ];

		GeometryReal previousLevel = edgeEquation( previousVertex );
		bool isPreviousInterior = previousLevel > 0;

		// Iterate over the vertices of the triangle:
		// -- If the vertex is interior:
		// ----- If the previous vertex was exterior, create and add the crossing edge
		// ----- Add the crossing edge
		// -- If the vertex is exterior
		// ----- If the previous vertex was interior, create and add the crossing edge

		for( int i=0 ; i<polygon.vertices.size() ; i++ )
		{
			Point2D< GeometryReal > currentVertex = polygon.vertices[i];
			GridMeshIntersectionKey currentVertexIndex = polygon.indices[i];
			previousEdgeIndex = nextEdgeIndex;
			nextEdgeIndex = polygon.edgeIndices[i];

			GeometryReal currentLevel = edgeEquation( currentVertex );
			bool isCurrentInterior = currentLevel > 0;

			if( isPreviousInterior )
			{
				outputVertices.push_back( previousVertex );
				outputIndices.push_back( previousVertexIndex );
				outputEdgeIndices.push_back( previousEdgeIndex );
			}

			if( isPreviousInterior!=isCurrentInterior )
			{
				GeometryReal alpha = -currentLevel / (previousLevel - currentLevel);
				outputVertices.push_back( previousVertex * alpha + currentVertex * ( 1-alpha ) );
				outputIndices.push_back( GridMeshIntersectionKey( previousEdgeIndex , edgeIndex ) );

				outputEdgeIndices.push_back( isPreviousInterior ? edgeIndex : previousEdgeIndex );
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

	// A function clipping the edges of a (convex) cell to the sides of a triangle
	template< typename GeometryReal >
	int ClipIndexedIntersectionPolygonToIndexedIntersectionTriangle
	(
		IndexedIntersectionPolygon< GeometryReal > &polygon ,
		const IndexedIntersectionTriangle< GeometryReal > &triangle
	)
	{
		Point2D< GeometryReal > triangleCenter = ( triangle.vertices[0] + triangle.vertices[1] + triangle.vertices[2] ) / 3;
		GridMeshIntersectionKey cornerKeys[6];
		for( int k=0 ; k<3 ; k++ )
		{
			// [NOTE] In sequentially clipping the polygon to the triangle's edges, 
			//        we can construct a vertex corresponding to a triangle's corner
			//        that is marked with the two edges of the triangle that meet there.
			//        In that case we need to clean up the indexing and mark it is a triangle corner
			//        so that the GridMeshIntersectionKey becomes consistent
			// [WARNING] This assumes that there isn't a grid edge with the same index as the triangle edge
			cornerKeys[ 2*k+0 ] = GridMeshIntersectionKey( triangle.edgeIndices[(k + 2) % 3] , triangle.edgeIndices[k] );
			cornerKeys[ 2*k+1 ] = GridMeshIntersectionKey( triangle.edgeIndices[k] , triangle.edgeIndices[(k + 2) % 3] );
		}
		GridMeshIntersectionKey cornerIndices[6];
		for( int k=0 ; k<3 ; k++ ) cornerIndices[2*k] = cornerIndices[2*k+1] = triangle.indices[k];

		bool reverseOrientation = false;

		// Clip the polygon to each edge of the triangle
		for( int k=0 ; k<3 ; k++ )
		{
			// Compute edge normal (pointing inside)
			EdgeIndex eIndex = CornerEdgeIndex( k );
			EdgeEquation< GeometryReal > edgeEquation( triangle.vertices[ eIndex[0] ] , triangle.vertices[ eIndex[1] ] , true );
			reverseOrientation = edgeEquation.makePositive( triangleCenter );

			ClipIndexedIntersectionPolygonToIndexedIntersectionEdge( polygon , edgeEquation , triangle.edgeIndices[k] );
		}

		// Identify triangle corners
		for( int i=0 ; i<polygon.indices.size() ; i++ ) for( int k=0 ; k<6 ; k++ ) if( polygon.indices[i]==cornerKeys[k] ) polygon.indices[i] = cornerIndices[k];

		if( reverseOrientation )
		{
			int n = (int)polygon.vertices.size();
			std::vector< Point2D < GeometryReal > > reversedVertices(n);
			std::vector< GridMeshIntersectionKey > reversedIndices(n);
			std::vector < int > reversedEdges(n);
			for( int k=0 ; k<n ; k++ )
			{
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

#ifdef NEW_RASTERIZER
#else // !NEW_RASTERIZER
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
#endif // NEW_RASTERIZER

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
}
#endif // POLYGON_CLIPPING_INCLUDED

