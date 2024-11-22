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

template< typename GeometryReal >
void InitializeAtlasCharts( AtlasMesh< GeometryReal > &atlasMesh , const std::vector< bool > &isBoundaryHalfEdge , int width , int height , std::vector< AtlasChart< GeometryReal > > &atlasCharts )
{
	atlasCharts.resize( atlasMesh.numCharts );

	for( int i=0 ; i<atlasMesh.numCharts ; i++ )
	{
		atlasCharts[i].minCorner = Point2D< GeometryReal >(1,1);
		atlasCharts[i].maxCorner = Point2D< GeometryReal >(0,0);
	}

	std::vector< int > lastVertexID( atlasMesh.numCharts , 0 );
	std::vector< int > chartVertexID( atlasMesh.vertices.size() , -1 );

	atlasMesh.triangleIndexInChart.resize( atlasMesh.triangles.size() );
	std::vector< int > lastTriangleID( atlasMesh.numCharts , 0 );

	for( int t=0 ; t<atlasMesh.triangles.size() ; t++ )
	{
		int chartID = atlasMesh.triangleChartIndex[t];

		atlasMesh.triangleIndexInChart[t] = lastTriangleID[chartID];
		lastTriangleID[chartID]++;

		atlasCharts[chartID].meshTriangleIndices.push_back(t);
		int vertexID[3];
		for( int k=0 ; k<3 ; k++ )
		{
			atlasCharts[chartID].atlasEdgeIndices.push_back( atlasMesh.halfEdgeToEdgeIndex[ 3*t+k ] );
			int atlasVertexID = atlasMesh.triangles[t][k];
			if( chartVertexID[atlasVertexID]==-1 )
			{
				chartVertexID[atlasVertexID] = lastVertexID[chartID];
				lastVertexID[chartID]++;
				Point2D< GeometryReal > vertexPos = atlasMesh.vertices[atlasVertexID];
				for( int c=0 ; c<2 ; c++ )
				{
					atlasCharts[chartID].minCorner[c] = std::min< GeometryReal >( vertexPos[c] , atlasCharts[chartID].minCorner[c] );
					atlasCharts[chartID].maxCorner[c] = std::max< GeometryReal >( vertexPos[c] , atlasCharts[chartID].maxCorner[c] );
				}

				atlasCharts[chartID].vertices.push_back( vertexPos );
				atlasCharts[chartID].meshVertexIndices.push_back( atlasMesh.vertexMap[atlasVertexID] );
			}
			vertexID[k] = chartVertexID[atlasVertexID];
		}

#ifdef DEBUG_ATLAS
		atlasCharts[chartID].triangles.push_back( _TriangleIndex( TriangleIndex( vertexID[0] , vertexID[1] , vertexID[2] ) , t ) );
#else // !DEBUG_ATLAS
		atlasCharts[chartID].triangles.push_back( TriangleIndex( vertexID[0] , vertexID[1] , vertexID[2] ) );
#endif // DEBUG_ATLAS
	}

	for( int i=0 ; i<atlasCharts.size() ; i++ )
	{
		Point2D< GeometryReal > midBBox = ( atlasCharts[i].minCorner+atlasCharts[i].maxCorner ) / 2;
		midBBox[0] *= (GeometryReal)width;
		midBBox[1] *= (GeometryReal)height;
		midBBox -= Point2D< GeometryReal >( (GeometryReal)0.5 , (GeometryReal)0.5 );
		midBBox = Point2D< GeometryReal >( (GeometryReal)floor( midBBox[0] ) , (GeometryReal)floor( midBBox[1] ) );
		atlasCharts[i].originCoords[0] = (int)round( midBBox[0] );
		atlasCharts[i].originCoords[1] = (int)round( midBBox[1] );

		midBBox += Point2D< GeometryReal >( (GeometryReal)0.5 , (GeometryReal)0.5 );
		midBBox[0] /= (GeometryReal)width;
		midBBox[1] /= (GeometryReal)height;
		atlasCharts[i].gridOrigin = midBBox;
	}

	for( int i=0 ; i<atlasCharts.size() ; i++ )
	{
		std::vector< int > &boundaryHalfEdges = atlasCharts[i].boundaryHalfEdges;
		for( int t=0 ; t<atlasCharts[i].meshTriangleIndices.size() ; t++ )
		{
			int tIndex = atlasCharts[i].meshTriangleIndices[t];
			for( int k=0 ; k<3 ; k++ ) if( isBoundaryHalfEdge[ 3*tIndex+k ] ) boundaryHalfEdges.push_back( 3*t+k );
		}
	}
}


template< typename GeometryReal >
void InitializeAtlasMesh
(
	const OrientedTexturedTriangleMesh< GeometryReal > &mesh ,
	int width ,
	int height ,
	AtlasMesh< GeometryReal > &atlasMesh ,
	std::vector< AtlasChart< GeometryReal > > &atlasCharts ,
	std::vector< int > &oppositeHalfEdge ,
	std::unordered_map< int , int > &boundaryVerticesIndices ,
	int &numBoundaryVertices ,
	bool &isClosedMesh ,
	bool verbose
)
{
	// Compute the 2D mesh and set the half-edge to edge map
	InitializeAtlasMesh( mesh , atlasMesh , width , height , verbose );

	// 1. Compute the opposite half-edges
	// 2. Identify the half-edges that are chart boundaries
	std::vector< int > boundaryHalfEdges;
	std::vector< bool > isBoundaryHalfEdge;
	InitializeBoundaryHalfEdges( mesh , boundaryHalfEdges , oppositeHalfEdge , isBoundaryHalfEdge , isClosedMesh );

	// Set the map taking boundary vertices to boundary vertex indices
	InitiallizeBoundaryVertices( mesh , boundaryHalfEdges , boundaryVerticesIndices , numBoundaryVertices );

	InitializeAtlasCharts( atlasMesh , isBoundaryHalfEdge , width , height , atlasCharts );
}
