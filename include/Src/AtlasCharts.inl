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
ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > >
AtlasMesh< GeometryReal >::getCharts
(
	const std::vector< bool > &isBoundaryHalfEdge ,
	unsigned int width ,
	unsigned int height
)
const
{
	ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > atlasCharts( _numCharts );

	for( unsigned int i=0 ; i<_numCharts ; i++ )
	{
		atlasCharts[ ChartIndex(i) ].minCorner = Point2D< GeometryReal >(1,1);
		atlasCharts[ ChartIndex(i) ].maxCorner = Point2D< GeometryReal >(0,0);
	}

	// The map taking a vertex index in the atlas and giving the index of the corresponding vertex within the chart
	std::vector< unsigned int > chartVertexID( SimpleTriangleMesh< GeometryReal , 2 >::vertices.size() , static_cast< unsigned int >(-1) );

	for( unsigned int t=0 ; t<SimpleTriangleMesh< GeometryReal , 2 >::triangles.size() ; t++ )
	{
		AtlasChart< GeometryReal > &atlasChart = atlasCharts[ triangleToChart( AtlasMeshTriangleIndex( t ) ) ];

		SimplexIndex< 2 > tri;
		for( unsigned int k=0 ; k<3 ; k++ )
		{
			atlasChart._chartHalfEdgeToAtlasEdge.push_back( halfEdgeToEdge( AtlasMeshHalfEdgeIndex( 3*t+k ) ) );
			unsigned int v = SimpleTriangleMesh< GeometryReal , 2 >::triangles[t][k];
			if( chartVertexID[v]==-1 )
			{
				chartVertexID[v] = (unsigned int)atlasChart.vertices.size();
				Point2D< GeometryReal > vertexPos = SimpleTriangleMesh< GeometryReal , 2 >::vertices[v];
				for( unsigned int c=0 ; c<2 ; c++ )
				{
					atlasChart.minCorner[c] = std::min< GeometryReal >( vertexPos[c] , atlasChart.minCorner[c] );
					atlasChart.maxCorner[c] = std::max< GeometryReal >( vertexPos[c] , atlasChart.maxCorner[c] );
				}

				atlasChart._chartToAtlasVertex.push_back( chartToAtlasVertex( ChartMeshVertexIndex(v) ) );
				atlasChart.vertices.push_back( vertexPos );
			}
			tri[k] = chartVertexID[v];
		}

		atlasChart.triangles.push_back(tri);
		atlasChart._chartToAtlasTriangle.push_back( AtlasMeshTriangleIndex(t) );
	}

	for( unsigned int i=0 ; i<atlasCharts.size() ; i++ )
	{
		Point2D< GeometryReal > midBBox = ( atlasCharts[ ChartIndex(i) ].minCorner+atlasCharts[ ChartIndex(i) ].maxCorner ) / 2;
		midBBox[0] *= (GeometryReal)width;
		midBBox[1] *= (GeometryReal)height;

		midBBox -= Point2D< GeometryReal >( (GeometryReal)0.5 , (GeometryReal)0.5 );
		midBBox = Point2D< GeometryReal >( (GeometryReal)floor( midBBox[0] ) , (GeometryReal)floor( midBBox[1] ) );
		atlasCharts[ ChartIndex(i) ].originCoords[0] = (int)round( midBBox[0] );
		atlasCharts[ ChartIndex(i) ].originCoords[1] = (int)round( midBBox[1] );

		midBBox += Point2D< GeometryReal >( (GeometryReal)0.5 , (GeometryReal)0.5 );
		midBBox[0] /= (GeometryReal)width;
		midBBox[1] /= (GeometryReal)height;
		atlasCharts[ ChartIndex(i) ].gridOrigin = midBBox;
	}

	// Add the boundary edges to the chart
	for( unsigned int i=0 ; i<atlasCharts.size() ; i++ )
	{
		AtlasChart< GeometryReal > &atlasChart = atlasCharts[ ChartIndex(i) ];
		std::vector< ChartMeshHalfEdgeIndex > & boundaryHalfEdges = atlasChart.boundaryHalfEdges;
		for( unsigned int he=0 ; he<atlasChart.numTriangles()*3 ; he++ )
			if( isBoundaryHalfEdge[ static_cast< unsigned int >( atlasChart.atlasHalfEdge( ChartMeshHalfEdgeIndex(he) ) ) ] )
				boundaryHalfEdges.push_back( ChartMeshHalfEdgeIndex(he) );
	}
	return atlasCharts;
}


template< typename GeometryReal >
ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > >
AtlasChart< GeometryReal >::GetCharts
(
	const TexturedTriangleMesh< GeometryReal > &mesh ,
	unsigned int width ,
	unsigned int height ,
	AtlasInfo &atlasInfo
)
{
	// Compute the 2D atlas
	AtlasMesh< GeometryReal > atlasMesh;
	atlasMesh.initialize( mesh );
	atlasMesh.jitter( width , height );

	// Compute the opposite half-edges on the surface mesh and mark whether half-edges are boundaries on the texture mesh
	std::vector< unsigned int > textureBoundaryHalfEdges;
	{
		std::vector< unsigned int > oppositeHalfEdges;
		mesh.setBoundaryHalfEdgeInfo( textureBoundaryHalfEdges , oppositeHalfEdges );
		atlasInfo.oppositeHalfEdges.resize( oppositeHalfEdges.size() );
		for( unsigned int i=0 ; i<oppositeHalfEdges.size() ; i++ ) atlasInfo.oppositeHalfEdges[ AtlasMeshHalfEdgeIndex(i) ] = AtlasMeshHalfEdgeIndex( oppositeHalfEdges[i] );
	}
	std::vector< bool > isTextureBoundaryHalfEdge( mesh.numTriangles()*3 , false );
	for( unsigned int i=0 ; i<textureBoundaryHalfEdges.size() ; i++ ) isTextureBoundaryHalfEdge[ textureBoundaryHalfEdges[i] ] = true;

	// Determine if the mesh is water-tight
	atlasInfo.isClosed = true;
	for( unsigned int i=0 ; i<atlasInfo.oppositeHalfEdges.size() ; i++ ) if( atlasInfo.oppositeHalfEdges[ AtlasMeshHalfEdgeIndex(i) ]==AtlasMeshHalfEdgeIndex(-1) ) atlasInfo.isClosed = false;

	// Set the map taking the indices of surface vertices lying on the texture boundary to vertex indices
	{
		std::map< unsigned int , unsigned int > atlasBoundaryVertexToIndex;
		mesh.setBoundaryVertexInfo( textureBoundaryHalfEdges , atlasBoundaryVertexToIndex );
		for( auto iter=atlasBoundaryVertexToIndex.begin() ; iter!=atlasBoundaryVertexToIndex.end() ; iter++ )
			atlasInfo.atlasMeshVertexToBoundaryVertex[ AtlasMeshVertexIndex( iter->first ) ] = AtlasMeshBoundaryVertexIndex( iter->second );
	}

	return atlasMesh.getCharts( isTextureBoundaryHalfEdge , width , height );
}
