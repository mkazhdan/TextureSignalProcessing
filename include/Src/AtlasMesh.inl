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

#include <Misha/Miscellany.h>
#include <Misha/Exceptions.h>

template< typename GeometryReal >
void InitializeAtlasMesh( const OrientedTexturedTriangleMesh< GeometryReal > &inputMesh , AtlasMesh< GeometryReal > &outputMesh , const int width , const int height , bool verbose )
{
	// Compute the mapping from triangles to charts
	InitializeTriangleChartIndexing( inputMesh , outputMesh.triangleChartIndex , outputMesh.numCharts );
	if( verbose ) printf( "Num Charts %d \n" , outputMesh.numCharts );

	// Compute the 2D mesh(es)
	{
		int lastVertexIndex = 0;
		std::set< IndexedVector2D< GeometryReal > , IndexedVector2DComparison< GeometryReal > > IndexedPointSet;
		typename std::set< IndexedVector2D< GeometryReal > , IndexedVector2DComparison< GeometryReal > >::iterator it;

		for( int t=0 ; t<inputMesh.triangles.size() ; t++ )
		{
			int cornerIndices[3];
			for( int k=0 ; k<3 ; k++ )
			{
				int currentCorner = -1;
				IndexedVector2D< GeometryReal > idxP( inputMesh.textureCoordinates[ 3*t+k ] , lastVertexIndex , inputMesh.triangles[t][k] );
				it = IndexedPointSet.find( idxP );
				if( it==IndexedPointSet.end() )
				{
					IndexedPointSet.insert( idxP );
					outputMesh.vertexMap.push_back( inputMesh.triangles[t][k] );
					currentCorner = lastVertexIndex;
					outputMesh.vertices.push_back( inputMesh.textureCoordinates[ 3*t+k ] );
					lastVertexIndex++;
				}
				else
				{
					IndexedVector2D< GeometryReal > indexPoint = *it;
					currentCorner = indexPoint.index;
				}
				cornerIndices[k] = currentCorner;
			}
			outputMesh.triangles.push_back( TriangleIndex( cornerIndices[0] , cornerIndices[1] , cornerIndices[2] ) );
		}
	}

	if( true ) // Jitter vertices lying on the grid to avoid degeneracies
	{
		GeometryReal precision = (GeometryReal)1e-6;
		for( int i=0 ; i<outputMesh.vertices.size() ; i++ ) for( int c=0 ; c<2 ; c++ )
		{
			GeometryReal dimSize = c == 0 ? (GeometryReal)width : (GeometryReal)height;
			GeometryReal scaled = outputMesh.vertices[i][c] * dimSize - (GeometryReal)0.5;
			GeometryReal offset = scaled - (GeometryReal)round(scaled);
			if (fabs(offset) < precision)
			{
				if( offset>0 ) scaled = round(scaled) + precision*( (GeometryReal)1. + Random< GeometryReal >() );
				else           scaled = round(scaled) - precision*( (GeometryReal)1. + Random< GeometryReal >() );
				outputMesh.vertices[i][c] = (GeometryReal)(scaled + 0.5) / dimSize;
			}
		}
	}

	// Set the mapping from half-edges to edges
	{
		// A map from vertex to pairs to half-edge indices
		std::unordered_map< unsigned long long , int > halfEdgeIndex;
		for( int i=0 ; i<outputMesh.triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ )
		{
			unsigned long long edgeKey = outputMesh.edgeKey(i,k);
			if( halfEdgeIndex.find(edgeKey)==halfEdgeIndex.end() ) halfEdgeIndex[edgeKey] = 3*i+k;
			else THROW( "Non oriented manifold mesh" );	// If the same half-edge appears twice
		}

		int lastEdgeIndex = 0;
		std::vector< int > halfEdgeToEdgeIndex( 3*outputMesh.triangles.size() , -1 );
		for( int i=0 ; i<outputMesh.triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ )
		{
			int currentHalfEdgeIndex = 3*i+k;
			unsigned long long oppositeEdgeKey = outputMesh.edgeKey(i,k,true);

			// Set the edge associated to both halves of the edge (once)
			if( halfEdgeIndex.find(oppositeEdgeKey)!=halfEdgeIndex.end() ) // If the opposite edge exists
			{
				int oppositeHalfEdgeIndex = halfEdgeIndex[oppositeEdgeKey];

				if( currentHalfEdgeIndex<oppositeHalfEdgeIndex ) halfEdgeToEdgeIndex[currentHalfEdgeIndex] = halfEdgeToEdgeIndex[oppositeHalfEdgeIndex] = lastEdgeIndex++;
			}
			else halfEdgeToEdgeIndex[currentHalfEdgeIndex] = lastEdgeIndex++;
		}
		for( int i=0 ; i<outputMesh.triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ ) if( halfEdgeToEdgeIndex[3*i+k]==-1 ) THROW( "Non indexed half edge" );

		outputMesh.halfEdgeToEdgeIndex = halfEdgeToEdgeIndex;
	}
}

template< typename GeometryReal >
void InitializeBoundaryHalfEdges( const OrientedTexturedTriangleMesh< GeometryReal > &mesh , std::vector< int > &boundaryHalfEdges , std::vector< int > &oppositeHalfEdge , std::vector< bool > &isBoundaryHalfEdge , bool &isClosedMesh )
{
	oppositeHalfEdge = mesh.getOppositeHalfEdges( isClosedMesh );

	isBoundaryHalfEdge.resize( 3*mesh.triangles.size() , false );

	for( int i=0 ; i<mesh.triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ )
	{
		int heIndex = 3*i+k;
		if( oppositeHalfEdge[heIndex]!=-1 )
		{
			int _heIndex = oppositeHalfEdge[heIndex];

			int _i = _heIndex / 3 , _k = _heIndex % 3;
			if( mesh.textureCoordinates[ 3*i+(k+1)%3 ][0] == mesh.textureCoordinates[ 3*_i+_k ][0] &&
				mesh.textureCoordinates[ 3*i+(k+1)%3 ][1] == mesh.textureCoordinates[ 3*_i+_k ][1] &&
				mesh.textureCoordinates[ 3*i+k ][0] == mesh.textureCoordinates[ 3*_i+(_k+1)%3 ][0] &&
				mesh.textureCoordinates[ 3*i+k ][1] == mesh.textureCoordinates[ 3*_i+(_k+1)%3 ][1] ) ;
			else
			{
				boundaryHalfEdges.push_back( heIndex );
				isBoundaryHalfEdge[heIndex] = true;
			}
		}
		else
		{
			isClosedMesh = false;
			boundaryHalfEdges.push_back( heIndex );
			isBoundaryHalfEdge[heIndex] = true;
		}
	}
}

template< typename GeometryReal >
void InitiallizeBoundaryVertices( const OrientedTexturedTriangleMesh< GeometryReal > &mesh , const std::vector< int > &boundaryHalfEdges , std::unordered_map< int , int > &boundaryVerticesIndices , int &lastBoundaryIndex )
{
	lastBoundaryIndex = 0;

	for( int b=0 ; b< boundaryHalfEdges.size() ; b++ )
	{
		int he = boundaryHalfEdges[b];
		int i = he / 3 , k = he % 3;
		if( boundaryVerticesIndices.find( mesh.triangles[i][k] )==boundaryVerticesIndices.end() ) boundaryVerticesIndices[ mesh.triangles[i][k] ] = lastBoundaryIndex++;
		if( boundaryVerticesIndices.find( mesh.triangles[i][(k+1)%3] )==boundaryVerticesIndices.end() ) boundaryVerticesIndices[ mesh.triangles[i][(k+1)%3] ] = lastBoundaryIndex++;
	}
}
