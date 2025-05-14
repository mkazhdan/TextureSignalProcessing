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
		IndexedPolygon< GeometryReal > &polygon ,
		const EdgeEquation< GeometryReal > &edgeEquation , 
		AtlasMeshEdgeIndex edgeIndex ,
		ChartMeshVertexIndex atlasVertexIndices[2]
	)
	{
		std::vector< Point2D< GeometryReal > > outputVertices;
		std::vector< ChartMeshVertexIndex > outputVertexIndices; 
		std::vector< AtlasMeshEdgeIndex > outputEdgeIndices;
		std::vector< AtlasMeshEdgeIndex > outputParentVertexEdgeIndices;
		outputVertices.reserve( polygon.vertices.size()+3 );
		outputVertexIndices.reserve( polygon.vertices.size()+3 );
		outputEdgeIndices.reserve( polygon.vertices.size()+3 );
		outputParentVertexEdgeIndices.reserve( polygon.vertices.size()+3 );

		Point2D< GeometryReal > previousVertex = polygon.vertices[polygon.vertices.size() - 1];
		ChartMeshVertexIndex previousVertexIndex = polygon.vertexIndices[ polygon.vertices.size()-1 ];
		AtlasMeshEdgeIndex previousEdgeIndex = polygon.atlasEdgeIndices[ polygon.vertices.size()-2 ];
		AtlasMeshEdgeIndex     nextEdgeIndex = polygon.atlasEdgeIndices[ polygon.vertices.size()-1 ];
		AtlasMeshEdgeIndex previousVertexEdgeSupport = polygon.atlasVertexParentEdge[polygon.vertices.size() - 1];

		int previousVertexType = -2;
		GeometryReal previousLevel;

		bool emptyPolygon = true;

		bool isCurrentVertexACorner = previousVertexIndex!=ChartMeshVertexIndex(-1) && ( previousVertexIndex==atlasVertexIndices[0] || previousVertexIndex==atlasVertexIndices[1] );
		bool isOnTheEdge = previousVertexEdgeSupport!=AtlasMeshEdgeIndex(-1) && previousVertexEdgeSupport==edgeIndex;
		if( isOnTheEdge || isCurrentVertexACorner )
		{
			previousVertexType = 0;
			previousLevel = 0;
		}
		else
		{
			previousLevel = edgeEquation( previousVertex );
			previousVertexType = previousLevel > 0 ? 1 : -1;
		}

		for( unsigned int i=0 ; i<polygon.vertices.size() ; i++ )
		{
			Point2D< GeometryReal > currentVertex = polygon.vertices[i];
			ChartMeshVertexIndex currentVertexIndex = polygon.vertexIndices[i];
			int currentVertexType = -2;
			AtlasMeshEdgeIndex currentVertexEdgeSupport = polygon.atlasVertexParentEdge[i];
			GeometryReal currentLevel;

			isCurrentVertexACorner = currentVertexIndex!=ChartMeshVertexIndex(-1) && ( currentVertexIndex==atlasVertexIndices[0] || currentVertexIndex==atlasVertexIndices[1] );
			isOnTheEdge = currentVertexEdgeSupport!=AtlasMeshEdgeIndex(-1) && currentVertexEdgeSupport==edgeIndex;
			//isPreviousEdgeColinear = isNextEdgeColinear;
			previousEdgeIndex = nextEdgeIndex;

			nextEdgeIndex = polygon.atlasEdgeIndices[i];
			//isNextEdgeColinear = nextEdgeIndex != -1 && nextEdgeIndex == edgeIndex;

			if( isOnTheEdge || isCurrentVertexACorner )
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
					outputVertexIndices.push_back( ChartMeshVertexIndex(-1) );
					outputEdgeIndices.push_back( previousEdgeIndex );
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
					outputVertexIndices.push_back( ChartMeshVertexIndex(-1) );
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
			polygon.vertexIndices.clear();
			polygon.atlasEdgeIndices.clear();
			polygon.atlasVertexParentEdge.clear();
		}
		else{
			polygon.vertices = outputVertices;
			polygon.vertexIndices = outputVertexIndices;
			polygon.atlasEdgeIndices = outputEdgeIndices;
			polygon.atlasVertexParentEdge = outputParentVertexEdgeIndices;

			if( polygon.vertices.size()!=polygon.vertexIndices.size() || polygon.vertices.size()!=polygon.atlasEdgeIndices.size() || polygon.vertices.size()!=polygon.atlasVertexParentEdge.size() )
				MK_THROW( "Polygon array size does not match" );

			//Check for non consecutive colinear edges
			for( unsigned int i=0 ; i<polygon.atlasEdgeIndices.size() ; i++ )
				if( polygon.atlasEdgeIndices[i]!=AtlasMeshEdgeIndex(-1) && polygon.atlasEdgeIndices[i]==polygon.atlasEdgeIndices[ (i+1)%polygon.atlasEdgeIndices.size() ] )
					MK_THROW( "Unexpected consecutive colinear edges" );
		}
	}

	// Only for convex polygons
	template< typename GeometryReal >
	unsigned int ClipPartiallyIndexedPolygonToIndexedTriangle( IndexedPolygon< GeometryReal > &polygon , const IndexedTriangle< GeometryReal > &triangle )
	{
		Point2D< GeometryReal > triangleCenter = ( triangle.vertices[0] + triangle.vertices[1] + triangle.vertices[2] ) / 3;

		for( unsigned int k=0 ; k<3 ; k++ )
		{
			SimplexIndex< 1 > eIndex = OutgoingEdgeIndex( k );

			EdgeEquation< GeometryReal > edgeEquation( triangle.vertices[ eIndex[0] ] , triangle.vertices[ eIndex[1] ] );
			edgeEquation.makePositive( triangleCenter );

			ChartMeshVertexIndex vertexIndices[2] = { triangle.vertexIndices[ eIndex[0] ] , triangle.vertexIndices[ eIndex[1] ] };
			ClipPartiallyIndexedPolygonToIndexedEdge( polygon , edgeEquation , triangle.atlasEdgeIndices[k] , vertexIndices );
		}

		return (unsigned int)polygon.vertices.size();
	}

	template< typename GeometryReal >
	void SetIndexedPolygonFromTriangle( const IndexedTriangle< GeometryReal > &triangle , IndexedPolygon< GeometryReal > &polygon )
	{
		for( int k=0 ; k<3 ; k++ )
		{
			polygon.vertices.push_back( triangle.vertices[k] );
			polygon.indices.push_back( triangle.indices[k] );
			polygon.vertexIndices.push_back( triangle.vertexIndices[k] );
			polygon.atlasEdgeIndices.push_back( triangle.atlasEdgeIndices[k] );
			polygon.atlasVertexParentEdge.push_back( triangle.atlasVertexParentEdge[k]) ;
		}
	}

	// Points in general positions
	// Clipping a convex polygon with an edge equation associated to the given edge index
	template< typename GeometryReal , typename MeshEdgesToKey /* = std::function< GridMeshIntersectionKey ( AtlasMeshEdgeIndex , AtlasMeshIndex ) > */ >
	void ClipIndexedIntersectionPolygonToIndexedIntersectionEdge
	(
		IndexedIntersectionPolygon< GeometryReal > &polygon ,
		const EdgeEquation< GeometryReal > &edgeEquation , 
		AtlasGridOrMeshEdgeIndex edgeIndex ,
		MeshEdgesToKey m2k
	)
	{
		static_assert( std::is_convertible_v< MeshEdgesToKey , std::function< GridMeshIntersectionKey ( AtlasMeshEdgeIndex , AtlasMeshEdgeIndex ) > > , "MeshEdgesToKey poorly formed" );

		std::vector< Point2D< GeometryReal > > outputVertices;
		std::vector< GridMeshIntersectionKey > outputCornerKeys;
		std::vector< AtlasGridOrMeshEdgeIndex > outputEdgeIndices;

		Point2D< GeometryReal > previousVertex = polygon.vertices[ polygon.vertices.size()-1 ];
		GridMeshIntersectionKey previousCornerKey = polygon.cornerKeys[ polygon.vertices.size()-1 ];
		AtlasGridOrMeshEdgeIndex previousEdgeIndex = polygon.outgoingEdgeIndices[ polygon.vertices.size()-1 ];

		GeometryReal previousLevel = edgeEquation( previousVertex );
		bool isPreviousInterior = previousLevel > 0;

		// Iterate over the vertices of the triangle:
		// -- If the vertex is interior:
		// ----- If the previous vertex was exterior, create and add the crossing edge
		// ----- Add the crossing edge
		// -- If the vertex is exterior
		// ----- If the previous vertex was interior, create and add the crossing edge

		for( unsigned int i=0 ; i<polygon.vertices.size() ; i++ )
		{
			Point2D< GeometryReal > currentVertex = polygon.vertices[i];
			GridMeshIntersectionKey currentCornerKey = polygon.cornerKeys[i];

			GeometryReal currentLevel = edgeEquation( currentVertex );
			bool isCurrentInterior = currentLevel > 0;

			if( isPreviousInterior )
			{
				outputVertices.push_back( previousVertex );
				outputCornerKeys.push_back( previousCornerKey );
				outputEdgeIndices.push_back( previousEdgeIndex );
			}

			if( isPreviousInterior!=isCurrentInterior )
			{
				GeometryReal alpha = -currentLevel / (previousLevel - currentLevel);
				outputVertices.push_back( previousVertex * alpha + currentVertex * ( 1-alpha ) );
				if( std::optional< AtlasMeshEdgeIndex > m = edgeIndex.mesh() )
				{
					// [WARNING] In the case that the previous edge index is of atlas type, we are creating an invalid GridMeshIntersectionKey
					if     ( std::optional< AtlasGridEdgeIndex >  g = previousEdgeIndex.grid() ) outputCornerKeys.push_back( GridMeshIntersectionKey( * g , *m ) );
					else if( std::optional< AtlasMeshEdgeIndex > _m = previousEdgeIndex.mesh() ) outputCornerKeys.push_back( m2k( *_m , *m ) );
					else MK_THROW( "Bad previous edge index" );
				}
				else MK_THROW( "Expected mesh edge type" );

				// If the previous is interior, the edge emenating from the new vertex follows the introduced edge
				outputEdgeIndices.push_back( isPreviousInterior ? edgeIndex : previousEdgeIndex );
			}

			previousVertex = currentVertex;
			previousLevel = currentLevel;
			previousCornerKey = currentCornerKey;
			isPreviousInterior = isCurrentInterior;

			previousEdgeIndex = polygon.outgoingEdgeIndices[i];
		}

		polygon.vertices = outputVertices;
		polygon.cornerKeys = outputCornerKeys;
		polygon.outgoingEdgeIndices = outputEdgeIndices;
	}

	// A function clipping the edges of a (convex) cell to the sides of a triangle
	template< typename GeometryReal >
	unsigned int ClipIndexedIntersectionPolygonToIndexedIntersectionTriangle
	(
		IndexedIntersectionPolygon< GeometryReal > &polygon ,
		const IndexedIntersectionTriangle< GeometryReal > &triangle
	)
	{
		Point2D< GeometryReal > triangleCenter = ( triangle.vertices[0] + triangle.vertices[1] + triangle.vertices[2] ) / 3;

		bool reverseOrientation = false;

		std::pair< AtlasMeshEdgeIndex , AtlasMeshEdgeIndex > cornerEdges[6];
		for( unsigned int k=0 ; k<=2 ; k++ )
		{
			std::optional< AtlasMeshEdgeIndex > m1 = triangle.outgoingEdgeIndices[k].mesh();
			std::optional< AtlasMeshEdgeIndex > m2 = triangle.outgoingEdgeIndices[(k+2)%3].mesh();
			if( !m1 || !m2 ) MK_THROW( "Expected incident edges to be mesh edges" );
			cornerEdges[2*k+0] = std::make_pair( *m1 , *m2 );
			cornerEdges[2*k+1] = std::make_pair( *m2 , *m1 );
		}

		auto IntersectingMeshEdgesToCornerKey = [&]( AtlasMeshEdgeIndex m1 , AtlasMeshEdgeIndex m2 )
			{
				std::pair< AtlasMeshEdgeIndex , AtlasMeshEdgeIndex > m = std::make_pair( m1 , m2 );
				for( unsigned int i=0 ; i<6 ; i++ ) if( m==cornerEdges[i] ) return triangle.cornerKeys[i/2];
				MK_THROW( "Could not match mesh edges: " , m1 , " : " , m2 );
				return GridMeshIntersectionKey();
			};

		// Clip the polygon to each edge of the triangle
		for( unsigned int k=0 ; k<3 ; k++ )
		{
			// Compute edge normal (pointing inside)
			SimplexIndex< 1 > eIndex = OutgoingEdgeIndex( k );
			EdgeEquation< GeometryReal > edgeEquation( triangle.vertices[ eIndex[0] ] , triangle.vertices[ eIndex[1] ] , true );
#ifdef SANITY_CHECK
			bool _reverseOrientation = edgeEquation.makePositive( triangleCenter );
			if( !k ) reverseOrientation = _reverseOrientation;
			else if( reverseOrientation!=_reverseOrientation ) MK_THROW( "Inconsistent orientation" );
#else // !SANITY_CHECK
			reverseOrientation = edgeEquation.makePositive( triangleCenter );
#endif // SANITY_CHECK

			ClipIndexedIntersectionPolygonToIndexedIntersectionEdge( polygon , edgeEquation , triangle.outgoingEdgeIndices[k] , IntersectingMeshEdgesToCornerKey );
		}

		if( reverseOrientation )
		{
			unsigned int n = (unsigned int)polygon.vertices.size();
			std::vector< Point2D < GeometryReal > > reversedVertices(n);
			std::vector< GridMeshIntersectionKey > reversedCornerKeys(n);
			std::vector < AtlasGridOrMeshEdgeIndex > reversedEdges(n);
			for( unsigned int k=0 ; k<n ; k++ )
			{
				reversedVertices[k] = polygon.vertices[ n-1-k ];
				reversedCornerKeys[k] = polygon.cornerKeys[ n-1-k ];
				// Need to grab the edge from the vertex before
				reversedEdges[k] = polygon.outgoingEdgeIndices[ ( n-1-k + (n-1) )%n ];
			}
			polygon.vertices = reversedVertices;
			polygon.cornerKeys = reversedCornerKeys;
			polygon.outgoingEdgeIndices = reversedEdges;
		}

		return (unsigned int)polygon.vertices.size();
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
}
#endif // POLYGON_CLIPPING_INCLUDED

