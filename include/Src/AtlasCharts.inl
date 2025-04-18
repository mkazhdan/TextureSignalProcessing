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
std::vector< AtlasChart< GeometryReal > >
AtlasMesh< GeometryReal >::getCharts
(
	const std::vector< bool > &isBoundaryHalfEdge ,
	unsigned int width ,
	unsigned int height
)
const
{
	std::vector< AtlasChart< GeometryReal > > atlasCharts( numCharts );

	for( unsigned int i=0 ; i<numCharts ; i++ )
	{
		atlasCharts[i].minCorner = Point2D< GeometryReal >(1,1);
		atlasCharts[i].maxCorner = Point2D< GeometryReal >(0,0);
	}

	// The map taking a vertex index in the atlas and giving the index of the corresponding vertex within the chart
	std::vector< unsigned int > chartVertexID( vertices.size() , static_cast< unsigned int >(-1) );

	for( unsigned int t=0 ; t<triangles.size() ; t++ )
	{
		AtlasChart< GeometryReal > &atlasChart = atlasCharts[ triangleChart(t) ];

		SimplexIndex< 2 > tri;
		for( unsigned int k=0 ; k<3 ; k++ )
		{
			atlasChart._chartHalfEdgeToAtlasEdge.push_back( halfEdgeToEdge( 3*t+k ) );
			unsigned int v = triangles[t][k];
			if( chartVertexID[v]==-1 )
			{
				chartVertexID[v] = atlasChart.vertices.size();
				Point2D< GeometryReal > vertexPos = vertices[v];
				for( unsigned int c=0 ; c<2 ; c++ )
				{
					atlasChart.minCorner[c] = std::min< GeometryReal >( vertexPos[c] , atlasChart.minCorner[c] );
					atlasChart.maxCorner[c] = std::max< GeometryReal >( vertexPos[c] , atlasChart.maxCorner[c] );
				}

				atlasChart._chartToSurfaceVertex.push_back( textureToSurfaceVertex(v) );
				atlasChart.vertices.push_back( vertexPos );
			}
			tri[k] = chartVertexID[v];
		}

		atlasChart.triangles.push_back(tri);
		atlasChart._chartToAtlasTriangle.push_back(t);
	}

	for( unsigned int i=0 ; i<atlasCharts.size() ; i++ )
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

	// Add the boundary edges to the chart
	for( int i=0 ; i<atlasCharts.size() ; i++ )
	{
		AtlasChart< GeometryReal > &atlasChart = atlasCharts[i];
		std::vector< unsigned int > &boundaryHalfEdges = atlasChart.boundaryHalfEdges;
		for( unsigned int t=0 ; t<atlasChart.triangles.size() ; t++ )
		{
			unsigned int tIndex = atlasChart.atlasTriangle(t);
			for( unsigned int k=0 ; k<3 ; k++ ) if( isBoundaryHalfEdge[ 3*tIndex+k ] ) boundaryHalfEdges.push_back( 3*t+k );
		}

#ifdef REORDER_BOUNDARY
		// Re-orer the boundary edges so that they are in sequence
		struct Edge
		{
			EdgeIndex edge;
			unsigned int index;
			bool processed;
		};
		std::vector< Edge > _boundaryHalfEdges( boundaryHalfEdges.size() );
		for( unsigned int i=0 ; i<boundaryHalfEdges.size() ; i++ )
		{
			unsigned he = boundaryHalfEdges[i];
			unsigned int t = he / 3 , k = he % 3;
			_boundaryHalfEdges[i].index = he;
			_boundaryHalfEdges[i].edge = atlasChart.edgeIndex( he );
			_boundaryHalfEdges[i].processed = false;
		}
		std::sort( _boundaryHalfEdges.begin() , _boundaryHalfEdges.end() , []( const Edge &e1 , const Edge &e2 ){ return e1.edge[0]<e2.edge[0]; } );

		std::function< unsigned int ( unsigned int , const Edge * , unsigned int , unsigned int ) > FindEdge = [&]( unsigned int v , const Edge *edges , unsigned int sz , unsigned int off )
			{
				if( sz==1 )
				{
					if( edges[0].edge[0]==v ) return off;
					else
					{
						MK_ERROR_OUT( "Could not find vertex: " , v );
						return static_cast< unsigned int >(-1);
					}
				}
				else if( v<edges[sz/2].edge[0] ) return FindEdge( v , edges , sz/2 , off );
				else return FindEdge( v , edges+(sz/2) , sz-(sz/2) , off+(sz/2) );
			};
		unsigned int e=0;

		for( unsigned int i=0 ; i<boundaryHalfEdges.size() ; i++ )
		{
			e = FindEdge( _boundaryHalfEdges[e].edge[1] , &_boundaryHalfEdges[0] , (unsigned int)_boundaryHalfEdges.size() , 0 );
			if( _boundaryHalfEdges[e].processed ) MK_ERROR_OUT( "Edge already processed" );
			_boundaryHalfEdges[e].processed = true;
			boundaryHalfEdges[e] = _boundaryHalfEdges[e].index;
		}
#endif // REORDER_BOUNDARY
	}
	return atlasCharts;
}


template< typename GeometryReal >
std::vector< AtlasChart< GeometryReal > >
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
	mesh.setBoundaryHalfEdgeInfo( textureBoundaryHalfEdges , atlasInfo.oppositeHalfEdges );
	std::vector< bool > isTextureBoundaryHalfEdge( mesh.numTriangles()*3 , false );
	for( unsigned int i=0 ; i<textureBoundaryHalfEdges.size() ; i++ ) isTextureBoundaryHalfEdge[ textureBoundaryHalfEdges[i] ] = true;

	// Determine if the mesh is water-tight
	atlasInfo.isClosed = true;
	for( unsigned int i=0 ; i<atlasInfo.oppositeHalfEdges.size() ; i++ ) if( atlasInfo.oppositeHalfEdges[i]==-1 ) atlasInfo.isClosed = false;

	// Set the map taking the indices of surface vertices lying on the texture boundary to vertex indices
#ifdef NEW_CODE
	mesh.setBoundaryVertexInfo( textureBoundaryHalfEdges , atlasInfo.surfaceBoundaryVertexToIndex );
#else // !NEW_CODE
	mesh.setBoundaryVertexInfo( textureBoundaryHalfEdges , atlasInfo.chartToSurfaceBoundaryVertex );
#endif // NEW_CODE

	return atlasMesh.getCharts( isTextureBoundaryHalfEdge , width , height );
}
