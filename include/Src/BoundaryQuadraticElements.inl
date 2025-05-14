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


template< typename GeometryReal >
void SetIndexedPolygonFromBoundaryTriangle
(
	const BoundaryIndexedTriangle< GeometryReal > &triangle ,
	IndexedPolygon< GeometryReal > &polygon
)
{
	for( unsigned int k=0 ; k<3 ; k++ )
	{
		polygon.vertices.push_back( triangle.vertices[k] );
		polygon.indices.push_back( AtlasInteriorOrBoundaryNodeIndex( triangle.indices[k] ) );
		polygon.vertexIndices.push_back( triangle.vertexIndices[k] );
		polygon.atlasEdgeIndices.push_back( triangle.atlasEdgeIndices[k] );
		polygon.atlasVertexParentEdge.push_back( triangle.atlasVertexParentEdge[k] );
	}
}

template< typename GeometryReal >
IndexedIntersectionPolygon< GeometryReal > IndexedPolygonFromCell( unsigned int i , unsigned int j , const GridChart< GeometryReal > &gridChart )
{
	IndexedIntersectionPolygon< GeometryReal > polygon;

	polygon.vertices.resize(4);
	polygon.cornerKeys.resize(4);
	polygon.outgoingEdgeIndices.resize(4);

	polygon.vertices[0] = gridChart.nodePosition( i+0 , j+0 );
	polygon.vertices[1] = gridChart.nodePosition( i+1 , j+0 );
	polygon.vertices[2] = gridChart.nodePosition( i+1 , j+1 );
	polygon.vertices[3] = gridChart.nodePosition( i+0 , j+1 );

	polygon.cornerKeys[0] = GridMeshIntersectionKey::GridNodeKey( gridChart.nodeIndex( i+0 , j+0 ) );
	polygon.cornerKeys[1] = GridMeshIntersectionKey::GridNodeKey( gridChart.nodeIndex( i+1 , j+0 ) );
	polygon.cornerKeys[2] = GridMeshIntersectionKey::GridNodeKey( gridChart.nodeIndex( i+1 , j+1 ) );
	polygon.cornerKeys[3] = GridMeshIntersectionKey::GridNodeKey( gridChart.nodeIndex( i+0 , j+1 ) );

	polygon.outgoingEdgeIndices[0] = AtlasGridOrMeshEdgeIndex::FromGrid( gridChart.edgeIndex( i+0 , j+0 , 0 ) );
	polygon.outgoingEdgeIndices[1] = AtlasGridOrMeshEdgeIndex::FromGrid( gridChart.edgeIndex( i+1 , j+0 , 1 ) );
	polygon.outgoingEdgeIndices[2] = AtlasGridOrMeshEdgeIndex::FromGrid( gridChart.edgeIndex( i+0 , j+1 , 0 ) );
	polygon.outgoingEdgeIndices[3] = AtlasGridOrMeshEdgeIndex::FromGrid( gridChart.edgeIndex( i+0 , j+0 , 1 ) );

	return polygon;
}

// Computes the intersections of the edge with the grid lines and adds to the vector of IntersectionInfo
template< typename GeometryReal >
void AddEdgeGridIntersection
(
	Simplex< GeometryReal , 2 , 1 > edge ,
	const GridChart< GeometryReal > &gridChart ,
	AtlasMeshEdgeIndex edgeIndex ,
	std::vector< IntersectionInfo< GeometryReal > > &intersections
)
{
	GeometryReal fmin[2] , fmax[2];
	for( unsigned int dir=0 ; dir<2 ; dir++ )
	{
		fmin[dir] = std::min< GeometryReal >( 1.f , std::max< GeometryReal >( 0.f , std::min< GeometryReal >( edge[0][dir] , edge[1][dir] ) ) );
		fmax[dir] = std::max< GeometryReal >( 0.f , std::min< GeometryReal >( 1.f , std::max< GeometryReal >( edge[0][dir] , edge[1][dir] ) ) );
	}

	unsigned int imin[2] , imax[2];
	imin[0] = static_cast< unsigned int >( floor( fmin[0] / gridChart.cellSizeW ) );
	imin[1] = static_cast< unsigned int >( floor( fmin[1] / gridChart.cellSizeH ) );
	imax[0] = static_cast< unsigned int >( floor( fmax[0] / gridChart.cellSizeW ) );
	imax[1] = static_cast< unsigned int >( floor( fmax[1] / gridChart.cellSizeH ) );

	for( unsigned int dir=0 ; dir<2 ; dir++ ) for( unsigned int s=imin[dir] ; s<=imax[dir] ; s++ )
	{
		GeometryReal    cellSize = dir == 0 ? gridChart.cellSizeW : gridChart.cellSizeH;
		GeometryReal oppCellSize = dir == 0 ? gridChart.cellSizeH : gridChart.cellSizeW;
		GeometryReal level = (GeometryReal)s*cellSize;
		GeometryReal alpha = ( level-edge[0][dir] ) / ( edge[1][dir]-edge[0][dir] );

		if( alpha>0 && alpha<1 )
		{
			Point2D< GeometryReal > position = edge[0] * (GeometryReal)( 1.-alpha ) + edge[1] * alpha;
			int oppDimIndex = (int)floor(position[1-dir] / oppCellSize);
			int cellIntersectionIndices[2];
			cellIntersectionIndices[dir] = s;
			cellIntersectionIndices[1-dir] = oppDimIndex;

			IntersectionInfo< GeometryReal > info;

			// [NOTE] When checking dir=0, we are intersecting with the vertical axis (integer x-values) so the direction is 1
			info.intersectionKey = GridMeshIntersectionKey( gridChart.edgeIndex( cellIntersectionIndices[0] , cellIntersectionIndices[1] , 1-dir ) , edgeIndex );
			info.time = alpha;
			info.position = position;

			intersections.push_back( info );
		}
	}
}

template< typename GeometryReal >
void
InitializeChartBoundaryTriangleGridIntersections
(
	const AtlasChart< GeometryReal > & atlasChart ,
	const std::map< AtlasMeshVertexIndex , AtlasMeshBoundaryVertexIndex > & atlasMeshVertexToBoundaryVertex ,
	GridChart< GeometryReal > & gridChart ,
	AtlasCoveredTexelIndex endCoveredTexelIndex ,
	AtlasRefinedBoundaryVertexIndex & endBoundaryVertex ,
	std::map< GridMeshIntersectionKey , NodeInfo< GeometryReal , AtlasRefinedBoundaryVertexIndex > > & gridMeshIntersectionKeyToNodeInfo
)
{
	// Assuming the boundary edges have processed. Now processing intersections of triangles with boundary cells.
	auto GetSimplex = [&]( ChartMeshTriangleIndex t )
		{
			Simplex< double , 2 , 2 > simplex;
			SimplexIndex< 2 , ChartMeshVertexIndex > tri = atlasChart.triangleIndex( t );
			for( unsigned int k=0 ; k<=2 ; k++ )
			{
				simplex[k] = atlasChart.vertex( tri[k] ) - gridChart.corner;
				simplex[k][0] /= gridChart.cellSizeW;
				simplex[k][1] /= gridChart.cellSizeH;
			}
			return simplex;
		};

	auto GetIndexedTriangle = [&]( ChartMeshTriangleIndex t )
		{
			IndexedIntersectionTriangle< GeometryReal > indexedTriangle;
			SimplexIndex< 2 , ChartMeshVertexIndex > tri = atlasChart.triangleIndex( t );

			for( unsigned int k=0 ; k<3 ; k++ )
			{
				indexedTriangle.vertices[k] = atlasChart.vertex( tri[k] ) - gridChart.corner;
				indexedTriangle.cornerKeys[k] = GridMeshIntersectionKey::ChartVertexKey( tri[k] );
				indexedTriangle.outgoingEdgeIndices[k] = AtlasGridOrMeshEdgeIndex::FromMesh( atlasChart.atlasEdge( GetChartMeshHalfEdgeIndex( t , k ) ) );
			}
			return indexedTriangle;
		};

	auto GetIndexedIntersectionPolygon = [&]( unsigned int i , unsigned int j ){ return IndexedPolygonFromCell( i , j , gridChart ); };

	using Range = RegularGrid< 2 >::Range;
	using Index = RegularGrid< 2 >::Index;
	Range nodeRange , cellRange;
	nodeRange.second[0] = gridChart.width  , cellRange.second[0] = gridChart.width-1;
	nodeRange.second[1] = gridChart.height , cellRange.second[1] = gridChart.height-1;

	for( unsigned int t=0 ; t<atlasChart.numTriangles() ; t++ )
	{
		IndexedIntersectionTriangle< GeometryReal > indexedTriangle = GetIndexedTriangle( ChartMeshTriangleIndex(t) );

		// Compute the associated triangle in (shifted) texel coordinates
		Simplex< double , 2 , 2 > simplex = GetSimplex( ChartMeshTriangleIndex(t) );

		auto Kernel = [&]( Index I )
			{
				if( gridChart.cellType(I)==CellType::Boundary )
				{
					ChartBoundaryCellIndex cellID = gridChart.cellIndices(I).boundary;
#ifdef SANITY_CHECK
					if( cellID==ChartBoundaryCellIndex(-1) ) MK_THROW( "Boundary cell invalid ID" );
#endif // SANITY_CHECK
					IndexedIntersectionPolygon< GeometryReal > poly = GetIndexedIntersectionPolygon( I[0] , I[1] );
					ClipIndexedIntersectionPolygonToIndexedIntersectionTriangle( poly , indexedTriangle );
#ifdef SANITY_CHECK
					if( !poly.vertices.size() ) MK_THROW( "Expected triangle to intersect cell" );
#endif // SANITY_CHECK

					// Assign indices to the boundary-edge/grid intersections
					for( unsigned int i=0 ; i<poly.cornerKeys.size() ; i++ )
					{
						// If we haven't seen this vertex before...
						if( gridMeshIntersectionKeyToNodeInfo.find( poly.cornerKeys[i] )==gridMeshIntersectionKeyToNodeInfo.end() )
						{
#ifdef SANITY_CHECK
							if( std::optional< ChartMeshVertexIndex > v = poly.cornerKeys[i].chartVertex() )
								if( atlasMeshVertexToBoundaryVertex.find( atlasChart.atlasVertex( *v ) )!=atlasMeshVertexToBoundaryVertex.end() )
									MK_THROW( "Boundary vertices should have already been processed" );
#endif // SANITY_CHECK
							AtlasRefinedBoundaryVertexIndex boundaryVertex = endBoundaryVertex++;

							gridMeshIntersectionKeyToNodeInfo[ poly.cornerKeys[i] ] = NodeInfo< GeometryReal , AtlasRefinedBoundaryVertexIndex >( boundaryVertex , poly.vertices[i] );

#ifdef SANITY_CHECK
#pragma message( "[WARNING] Converting AtlasRefinedBoundaryVertexIndex -> AtlasInteriorOrBoundaryNodeIndex" )
#endif // SANITY_CHECK
							gridChart.auxiliaryNodes.emplace_back( AtlasInteriorOrBoundaryNodeIndex( endCoveredTexelIndex ) + static_cast< unsigned int >(boundaryVertex) , poly.vertices[i] );
						}
					}
				}
			};
		Rasterizer2D::RasterizeSupports< true , true >( simplex , Kernel , cellRange );
	}
}

template< typename GeometryReal >
void InitializeChartBoundaryEdgeGridIntersections
(
	const AtlasChart< GeometryReal > & atlasChart ,
	const std::map< AtlasMeshVertexIndex , AtlasMeshBoundaryVertexIndex > & atlasMeshVertexToBoundaryVertex ,
	GridChart< GeometryReal > & gridChart ,
	AtlasCoveredTexelIndex endCoveredTexelIndex ,
	AtlasRefinedBoundaryVertexIndex & endBoundaryVertex ,
	std::map< AtlasMeshHalfEdgeIndex , std::vector< IntersectionInfo< GeometryReal > > > &atlasBoundaryHalfEdgeToIntersectionInfos ,
	std::map< SimplexIndex< 1 , AtlasRefinedBoundaryVertexIndex > , BoundarySegmentInfo< GeometryReal > > &segmentToBoundarySegmentInfo ,
	std::map< GridMeshIntersectionKey , NodeInfo< GeometryReal , AtlasRefinedBoundaryVertexIndex > > & gridMeshIntersectionKeyToNodeInfo
)
{
	for( unsigned int b=0 ; b<atlasChart.boundaryHalfEdges.size() ; b++ )
	{
		ChartMeshHalfEdgeIndex chartHalfEdge = atlasChart.boundaryHalfEdges[b];
		AtlasMeshHalfEdgeIndex atlasHalfEdge = atlasChart.atlasHalfEdge( chartHalfEdge );

		// The intersections of the boundary edge with the grid edges
		std::vector< IntersectionInfo< GeometryReal > > boundaryHalfEdgeIntersectionsInfo;

		SimplexIndex< 1 , ChartMeshVertexIndex > eIndex = atlasChart.edgeIndex( atlasChart.boundaryHalfEdges[b] );
		Simplex< GeometryReal , 2 , 1 > edge( atlasChart.vertex( ChartMeshVertexIndex( eIndex[0] ) ) - gridChart.corner , atlasChart.vertex( ChartMeshVertexIndex( eIndex[1] ) ) - gridChart.corner );

		// Add the end-points of the edge
		boundaryHalfEdgeIntersectionsInfo.emplace_back( GridMeshIntersectionKey::ChartVertexKey( eIndex[0] ) , edge[0] , 0 );
		boundaryHalfEdgeIntersectionsInfo.emplace_back( GridMeshIntersectionKey::ChartVertexKey( eIndex[1] ) , edge[1] , 1 );

		// Add the points where the edge intersects the grid
		AddEdgeGridIntersection( edge , gridChart , atlasChart.atlasEdge( chartHalfEdge ) , boundaryHalfEdgeIntersectionsInfo );

		// Sort intersections (required to define segments)
		std::sort( boundaryHalfEdgeIntersectionsInfo.begin() , boundaryHalfEdgeIntersectionsInfo.end() , IntersectionInfo< GeometryReal >::CompareByTime );

		// Assign indices to the boundary-edge/grid intersections
		for( unsigned int i=0 ; i<boundaryHalfEdgeIntersectionsInfo.size() ; i++ )
		{
			// If we haven't seen this vertex before...
			if( gridMeshIntersectionKeyToNodeInfo.find( boundaryHalfEdgeIntersectionsInfo[i].intersectionKey )==gridMeshIntersectionKeyToNodeInfo.end() )
			{
				AtlasRefinedBoundaryVertexIndex boundaryVertex;
				if( std::optional< ChartMeshVertexIndex > v = boundaryHalfEdgeIntersectionsInfo[i].intersectionKey.chartVertex() )	// Start/end vertex
				{
#ifdef SANITY_CHECK
					if( i!=0 && i!=boundaryHalfEdgeIntersectionsInfo.size()-1 ) MK_THROW( "Expected boundary vertex" );
#endif // SANITY_CHECK
					auto iter = atlasMeshVertexToBoundaryVertex.find( atlasChart.atlasVertex( *v ) );
#ifdef SANITY_CHECK
					if( iter==atlasMeshVertexToBoundaryVertex.end() ) MK_THROW( "Boundary vertex not found" );
#endif // SANITY_CHECK
#ifdef SANITY_CHECK
#pragma message( "[WARNING] Converting AtlasMeshBoundaryVertexIndex -> AtlasRefinedBoundaryVertexIndex" )
#endif // SANITY_CHECK
					boundaryVertex = AtlasRefinedBoundaryVertexIndex( iter->second );
				}
				else // Interior (intersection) vertex
				{
#ifdef SANITY_CHECK
					if( i==0 || i==boundaryHalfEdgeIntersectionsInfo.size()-1 ) MK_THROW( "Expected interior vertex" );
#endif // SANITY_CHECK
					boundaryVertex = endBoundaryVertex++;
				}

				gridMeshIntersectionKeyToNodeInfo[ boundaryHalfEdgeIntersectionsInfo[i].intersectionKey ]
					= NodeInfo< GeometryReal , AtlasRefinedBoundaryVertexIndex >( boundaryVertex , boundaryHalfEdgeIntersectionsInfo[i].position );

#ifdef SANITY_CHECK
#pragma message( "[WARNING] Converting AtlasRefinedBoundaryVertexIndex -> AtlasInteriorOrBoundaryNodeIndex" )
#endif // SANITY_CHECK
				gridChart.auxiliaryNodes.emplace_back( AtlasInteriorOrBoundaryNodeIndex( endCoveredTexelIndex ) + static_cast< unsigned int >(boundaryVertex) , boundaryHalfEdgeIntersectionsInfo[i].position );
			}
			boundaryHalfEdgeIntersectionsInfo[i].index = gridMeshIntersectionKeyToNodeInfo[ boundaryHalfEdgeIntersectionsInfo[i].intersectionKey ].index;
		}

		// Add the segment information
		for( unsigned int i=0 ; i<boundaryHalfEdgeIntersectionsInfo.size()-1 ; i++ )
		{
			SimplexIndex< 1 , AtlasRefinedBoundaryVertexIndex > eIndex( gridMeshIntersectionKeyToNodeInfo[ boundaryHalfEdgeIntersectionsInfo[i].intersectionKey ].index , gridMeshIntersectionKeyToNodeInfo[ boundaryHalfEdgeIntersectionsInfo[i+1].intersectionKey ].index );
			BoundarySegmentInfo< GeometryReal > segmentInfo( boundaryHalfEdgeIntersectionsInfo[i].time , boundaryHalfEdgeIntersectionsInfo[i+1].time , chartHalfEdge );

#ifdef SANITY_CHECK
			if( segmentToBoundarySegmentInfo.find( eIndex )!=segmentToBoundarySegmentInfo.end() ) MK_THROW( "Replicated segment key" );
#endif // SANITY_CHECK
			segmentToBoundarySegmentInfo[ eIndex ] = segmentInfo;
		}

		atlasBoundaryHalfEdgeToIntersectionInfos[ atlasHalfEdge ] = boundaryHalfEdgeIntersectionsInfo;
	}
}

template< typename GeometryReal >
ExplicitIndexVector< ChartBoundaryCellIndex , std::vector< std::pair< ChartMeshTriangleIndex , IndexedPolygon< GeometryReal > > > >
GetChartBoundaryPolygons
(
	const ExplicitIndexVector< AtlasMeshHalfEdgeIndex , AtlasMeshHalfEdgeIndex > & oppositeHalfEdge ,
	const AtlasChart< GeometryReal > &atlasChart ,
	GridChart< GeometryReal > &gridChart ,
	AtlasCoveredTexelIndex endCoveredTexelIndex ,
	unsigned int numBoundaryVertices ,
#ifdef SANITY_CHECK
	AtlasRefinedBoundaryVertexIndex endBoundaryVertex ,
#endif // SANITY_CHECK
	const std::map< AtlasMeshHalfEdgeIndex , std::vector< IntersectionInfo< GeometryReal > > > & atlasBoundaryHalfEdgeToIntersectionInfos ,
	const std::map< SimplexIndex< 1 , AtlasRefinedBoundaryVertexIndex > , BoundarySegmentInfo< GeometryReal > > & segmentToBoundarySegmentInfo ,
	const std::map< GridMeshIntersectionKey , NodeInfo< GeometryReal , AtlasRefinedBoundaryVertexIndex > > & gridMeshIntersectionKeyToNodeInfo
)
{
	ExplicitIndexVector< ChartBoundaryCellIndex , std::vector< std::pair< ChartMeshTriangleIndex , IndexedPolygon< GeometryReal > > > > boundaryPolygons;
	// A mapping giving the clipped polygons, keyed off of cell index
	std::map< ChartBoundaryCellIndex , std::vector< std::pair< ChartMeshTriangleIndex , std::vector< GridMeshIntersectionKey > > > > cellPolygons;

	auto GetIndexedTriangle = [&]( ChartMeshTriangleIndex t )
		{
			IndexedIntersectionTriangle< GeometryReal > indexedTriangle;
			SimplexIndex< 2 , ChartMeshVertexIndex > tri = atlasChart.triangleIndex( t );

			for( unsigned int k=0 ; k<3 ; k++ )
			{
				indexedTriangle.vertices[k] = atlasChart.vertex( tri[k] ) - gridChart.corner;
				indexedTriangle.cornerKeys[k] = GridMeshIntersectionKey::ChartVertexKey( tri[k] );
				indexedTriangle.outgoingEdgeIndices[k] = AtlasGridOrMeshEdgeIndex::FromMesh( atlasChart.atlasEdge( GetChartMeshHalfEdgeIndex( t , k ) ) );
			}
			return indexedTriangle;
		};

	auto GetIndexedIntersectionPolygon = [&]( unsigned int i , unsigned int j ){ return IndexedPolygonFromCell( i , j , gridChart ); };

	using Range = RegularGrid< 2 >::Range;
	using Index = RegularGrid< 2 >::Index;
	Range nodeRange , cellRange;
	nodeRange.second[0] = gridChart.width  , cellRange.second[0] = gridChart.width-1;
	nodeRange.second[1] = gridChart.height , cellRange.second[1] = gridChart.height-1;

	auto GetSimplex = [&]( ChartMeshTriangleIndex t )
		{
			Simplex< double , 2 , 2 > simplex;
			SimplexIndex< 2 , ChartMeshVertexIndex > tri = atlasChart.triangleIndex( t );
			for( unsigned int k=0 ; k<=2 ; k++ )
			{
				simplex[k] = atlasChart.vertex( tri[k] ) - gridChart.corner;
				simplex[k][0] /= gridChart.cellSizeW;
				simplex[k][1] /= gridChart.cellSizeH;
			}
			return simplex;
		};

	//(1) Rasterize all triangles and, if they fall into the support of a boundary cell, add them to the list of triangles associated with that cell
	for( unsigned int t=0 ; t<atlasChart.numTriangles() ; t++ )
	{
		IndexedIntersectionTriangle< GeometryReal > indexedTriangle = GetIndexedTriangle( ChartMeshTriangleIndex(t) );
		// Compute the associated triangle in (shifted) texel coordinates
		Simplex< double , 2 , 2 > simplex = GetSimplex( ChartMeshTriangleIndex(t) );

		auto Kernel = [&]( Index I )
			{
				if( gridChart.cellType(I)==CellType::Boundary )
				{
					ChartBoundaryCellIndex cellID = gridChart.cellIndices(I).boundary;
#ifdef SANITY_CHECK
					if( cellID==ChartBoundaryCellIndex(-1) ) MK_THROW( "Boundary cell invalid ID" );
#endif // SANITY_CHECK
					IndexedIntersectionPolygon< GeometryReal > poly = GetIndexedIntersectionPolygon( I[0] , I[1] );
					ClipIndexedIntersectionPolygonToIndexedIntersectionTriangle( poly , indexedTriangle );
#ifdef SANITY_CHECK
					if( !poly.vertices.size() ) MK_THROW( "Expected triangle to intersect cell" );
#endif // SANITY_CHECK

					std::vector< std::pair< ChartMeshTriangleIndex , std::vector< GridMeshIntersectionKey > > > & polygons = cellPolygons[cellID];
					polygons.emplace_back( ChartMeshTriangleIndex(t) , poly.cornerKeys );
				}
			};
		Rasterizer2D::RasterizeSupports< true , true >( simplex , Kernel , cellRange );
	}

	//(2) Process cells
	boundaryPolygons.resize( gridChart.numBoundaryCells() );
	for( auto iter=cellPolygons.begin() ; iter!=cellPolygons.end() ; iter++ )
	{
		ChartBoundaryCellIndex cellID = iter->first;
		std::vector< std::pair< ChartMeshTriangleIndex , std::vector< GridMeshIntersectionKey > > > & loopKeys = iter->second;

		// The position and index of boundary nodes/vertices
		std::vector< std::vector< NodeInfo< GeometryReal , AtlasInteriorOrBoundaryNodeIndex > > > loopNodes( loopKeys.size() );

		// For vertices on the original mesh, the index of the associated mesh vertex (-1 otherwise)
		std::vector< std::vector< ChartMeshVertexIndex > > loopAtlasVertexIndices( loopKeys.size() );

		// Remove vertices generated when non-boundary edges intersect the grid
		for( unsigned int i=0 ; i<loopKeys.size() ; i++ )
		{
			loopNodes[i].reserve( loopKeys[i].second.size() );
			loopAtlasVertexIndices[i].reserve( loopKeys[i].second.size() );

			// For each polygon vertex
			for( unsigned int j=0 ; j<loopKeys[i].second.size() ; j++ )
			{
				GridMeshIntersectionKey currentVertexKey = loopKeys[i].second[j];

				if( std::optional< AtlasGridVertexIndex > g = currentVertexKey.gridNode() )	// Texel node
				{
					unsigned int pi , pj;
					if( !gridChart.factorNodeIndex( *g , pi , pj ) ) MK_THROW( "Could not factor node index: " , *g );

					AtlasCoveredTexelIndex coveredTexelIndex = gridChart.texelIndices(pi,pj).covered;

#ifdef SANITY_CHECK
					// Confirm that the texel is inside the chart
					if( coveredTexelIndex==AtlasCoveredTexelIndex(-1) ) MK_THROW( "Invalid texel: " , Point2D< int >( pi , pj ) , " -> " , gridChart.texelIndices(pi,pj).covered , " : " , gridChart.texelIndices(pi, pj).interior );
#endif // SANITY_CHECK

#ifdef SANITY_CHECK
#pragma message( "[WARNING] Converting AtlasCoveredTexelIndex -> AtlasInteriorOrBoundaryNodeIndex" )
#endif // SANITY_CHECK
					loopNodes[i].emplace_back( AtlasInteriorOrBoundaryNodeIndex( coveredTexelIndex ) , gridChart.nodePosition(pi,pj) );
					loopAtlasVertexIndices[i].push_back( ChartMeshVertexIndex(-1) );
				}
				else
				{
					auto iter = gridMeshIntersectionKeyToNodeInfo.find( currentVertexKey );
					if( iter!=gridMeshIntersectionKeyToNodeInfo.end() )
					{
#ifdef SANITY_CHECK
#pragma message( "[WARNING] Converting AtlasRefinedBoundaryVertexIndex -> AtlasInteriorOrBoundaryNodeIndex" )
#endif // SANITY_CHECK
						loopNodes[i].emplace_back( AtlasInteriorOrBoundaryNodeIndex( iter->second.index ) + static_cast< unsigned int >(endCoveredTexelIndex) , iter->second.position );

						if( std::optional< ChartMeshVertexIndex > v = currentVertexKey.chartVertex() ) loopAtlasVertexIndices[i].push_back( *v );
						else                                                                           loopAtlasVertexIndices[i].push_back( ChartMeshVertexIndex(-1) );
					}
#ifdef SANITY_CHECK
					else MK_THROW( "Could not find vertex" );
#endif // SANITY_CHECK
				}
			}
		}

		// Add the intermediate vertices for boundary segments
		std::vector< std::vector< AtlasMeshEdgeIndex > > loopAtlasEdges( loopKeys.size() );
		std::vector< std::vector< AtlasMeshEdgeIndex > > loopAtlasVertexParentEdges( loopKeys.size() );
		for( unsigned int i=0 ; i<loopKeys.size() ; i++ )
		{
			std::vector< std::vector< AtlasRefinedBoundaryVertexIndex > > indicesToInsert;
			std::vector< std::vector< GeometryReal > > timesToInsert;
			std::vector< bool > isInsertionPos( loopKeys[i].second.size() , false );
			std::vector< AtlasMeshEdgeIndex > polygonAtlasEdgeIndex( loopKeys[i].second.size() , AtlasMeshEdgeIndex(-1) );
			std::vector< AtlasMeshEdgeIndex > polygonAtlasVertexParentEdges( loopKeys[i].second.size() , AtlasMeshEdgeIndex(-1) );	// The mesh edge the boundary vertex/node comes from

			for( unsigned int j=0 ; j<loopKeys[i].second.size() ; j++ )
			{
				GridMeshIntersectionKey currentVertexKey = loopKeys[i].second[j];
				GridMeshIntersectionKey    nextVertexKey = loopKeys[i].second[(j+1) % loopKeys[i].second.size()];

				if( std::optional< std::pair< AtlasGridEdgeIndex , AtlasMeshEdgeIndex > > x = currentVertexKey.intersection() ) polygonAtlasVertexParentEdges[j] = x->second;

				// If this is part of a boundary edge
				auto currentBoundaryNodeIter = gridMeshIntersectionKeyToNodeInfo.find( currentVertexKey );
				auto    nextBoundaryNodeIter = gridMeshIntersectionKeyToNodeInfo.find(    nextVertexKey );
				if( currentBoundaryNodeIter!=gridMeshIntersectionKeyToNodeInfo.end() && nextBoundaryNodeIter!=gridMeshIntersectionKeyToNodeInfo.end() )
				{
					SimplexIndex< 1 , AtlasRefinedBoundaryVertexIndex > eIndex( currentBoundaryNodeIter->second.index , nextBoundaryNodeIter->second.index );
					auto segmentToBoundarySegmentInfoIter = segmentToBoundarySegmentInfo.find( eIndex );
					if( segmentToBoundarySegmentInfoIter!=segmentToBoundarySegmentInfo.end() )
					{
						BoundarySegmentInfo< GeometryReal > segmentInfo = segmentToBoundarySegmentInfoIter->second;
						polygonAtlasEdgeIndex[j] = atlasChart.atlasEdge( segmentInfo.chartHalfEdge );

						AtlasMeshHalfEdgeIndex oppHalfEdge = oppositeHalfEdge[ atlasChart.atlasHalfEdge( segmentInfo.chartHalfEdge ) ];
						if( oppHalfEdge!=AtlasMeshHalfEdgeIndex(-1) )
						{
							GeometryReal startTime = segmentInfo.startTime;
							GeometryReal endTime = segmentInfo.endTime;

							// Find the information for the opposite half-edge
							auto atlasBoundaryHalfEdgeToIntersectionInfosIter = atlasBoundaryHalfEdgeToIntersectionInfos.find( oppHalfEdge );
							if( atlasBoundaryHalfEdgeToIntersectionInfosIter==atlasBoundaryHalfEdgeToIntersectionInfos.end() )
								MK_THROW( "Opposite edge intersections not found. Current  edge " , atlasChart.atlasHalfEdge( segmentInfo.chartHalfEdge ) , ". Opposite " , oppHalfEdge );

							std::vector< AtlasRefinedBoundaryVertexIndex > segmentIndicesToInsert;
							std::vector< GeometryReal > segmentTimesToInsert;

							// Iterate through the partitions of the opposite half-edge and, if they happen
							// within the time interval associated with this segment, insert them
							std::vector< IntersectionInfo< GeometryReal > > oppositeIntersection = atlasBoundaryHalfEdgeToIntersectionInfosIter->second;
							for( unsigned int oppIter=(unsigned int)oppositeIntersection.size()-2 ; oppIter!=0 ; oppIter-- )
							{
								GeometryReal reversedTime = (GeometryReal)1. - oppositeIntersection[oppIter].time;
								GeometryReal normalizedTime = ( reversedTime-startTime ) / ( endTime-startTime );
#ifdef INSERTION_EPSILON
								if( normalizedTime>INSERTION_EPSILON && normalizedTime<(1.-INSERTION_EPSILON) )
#else // !INSERTION_EPSILON
								if (normalizedTime > 0 && normalizedTime < 1)
#endif // INSERTION_EPSILON
								{
									segmentIndicesToInsert.push_back( oppositeIntersection[oppIter].index );
									segmentTimesToInsert.push_back( normalizedTime );
								}
							}


							if( segmentIndicesToInsert.size()>0 )
							{
								isInsertionPos[j] = true;
								indicesToInsert.push_back( segmentIndicesToInsert );
								timesToInsert.push_back( segmentTimesToInsert );
							}
						}
					}
					// [NOTE] Both end-point may be on the boundary (e.g. when two sides of a triangle pass through the same grid edge)
					//        but that may not make this a boundary segment
				}
			}

			// Do insertions
			unsigned int insertionCount = 0;
			std::vector< NodeInfo< GeometryReal , AtlasInteriorOrBoundaryNodeIndex > > expandedLoopNodes;
			std::vector< ChartMeshVertexIndex > expandedLoopAtlasVertexIndices;
			std::vector< AtlasMeshEdgeIndex > expandedLoopAtlasEdgeIndex;
			std::vector< AtlasMeshEdgeIndex > expandedVertexParentEdgeIndex;
			for( unsigned int j=0 ; j<loopKeys[i].second.size() ; j++ )
			{
				// Add information for original loop vertex

				expandedLoopNodes.push_back( loopNodes[i][j] );
				expandedLoopAtlasVertexIndices.push_back( loopAtlasVertexIndices[i][j] );

				AtlasMeshEdgeIndex currentSegmentAtlasEdgeIndex = polygonAtlasEdgeIndex[j];
				expandedLoopAtlasEdgeIndex.push_back( currentSegmentAtlasEdgeIndex );
				expandedVertexParentEdgeIndex.push_back( polygonAtlasVertexParentEdges[j] );

				// If there are vertices to insert from the opposite boundary, do that now
				if( isInsertionPos[j] )
				{
					const std::vector< AtlasRefinedBoundaryVertexIndex > & _indices = indicesToInsert[ insertionCount ];
					const std::vector< GeometryReal > & _times = timesToInsert[insertionCount];
					Point2D< GeometryReal > currentPos = loopNodes[i][j].position;
					Point2D< GeometryReal >    nextPos = loopNodes[i][ (j+1)%loopKeys[i].second.size() ].position;

					if( currentSegmentAtlasEdgeIndex==AtlasMeshEdgeIndex(-1) ) MK_THROW( "Invalid atlas edge index" );

					for( unsigned int k=0 ; k<_indices.size() ; k++ )
					{
						GeometryReal alpha = _times[k];
						Point2D< GeometryReal > interpolatedPos = (GeometryReal)( 1.-alpha )*currentPos + nextPos*alpha;
						if( true ) // Add orthogonal perturbation to avoid colinearity. Delanauy triangulation library failed for colinear inputs on some systems.
						{
							Point2D< GeometryReal > segmentDir = nextPos - currentPos;
							interpolatedPos += Point2D< GeometryReal >( -segmentDir[1] , segmentDir[0] ) * (GeometryReal)( ( Random< GeometryReal >()*2 - 1. ) * 1e-10 );
						}

						// Validate that the added node is not an original atlas vertex (and has been processed)
#ifdef SANITY_CHECK
						if( static_cast< unsigned int >(_indices[k])<numBoundaryVertices || static_cast< unsigned int >(_indices[k])>static_cast< unsigned int >(endBoundaryVertex) )
							MK_THROW( "Out of bounds index: " , _indices[k] , " not in " , numBoundaryVertices , " " , endBoundaryVertex );				
#endif // SANITY_CHECK

#ifdef SANITY_CHECK
#pragma message( "[WARNING] Converting AtlasRefinedBoundaryVertexIndex -> AtlasInteriorOrBoundaryNodeIndex" )
#endif // SANITY_CHECK
						gridChart.auxiliaryNodes.emplace_back( AtlasInteriorOrBoundaryNodeIndex( _indices[k] ) + static_cast< unsigned int >(endCoveredTexelIndex) , interpolatedPos );
						expandedLoopNodes.emplace_back( AtlasInteriorOrBoundaryNodeIndex( static_cast< unsigned int >(_indices[k]) ) + static_cast< unsigned int >(endCoveredTexelIndex) , interpolatedPos );

						expandedLoopAtlasVertexIndices.push_back( ChartMeshVertexIndex(-1) );
						expandedLoopAtlasEdgeIndex.push_back( currentSegmentAtlasEdgeIndex );
						expandedVertexParentEdgeIndex.push_back( currentSegmentAtlasEdgeIndex );
					}
					insertionCount++;
				}
			}
			if( insertionCount!=indicesToInsert.size() ) MK_THROW( "Intersection chains count does not match" );
			loopNodes[i] = expandedLoopNodes;
			loopAtlasVertexIndices[i] = expandedLoopAtlasVertexIndices;
			loopAtlasEdges[i] = expandedLoopAtlasEdgeIndex;
			loopAtlasVertexParentEdges[i] = expandedVertexParentEdgeIndex;
		}

		for( unsigned int i=0 ; i<loopNodes.size() ; i++ )
		{
			IndexedPolygon< GeometryReal > poly;
			poly.vertices.resize( loopNodes[i].size() ) , poly.indices.resize( loopNodes[i].size() );
			for( unsigned int j=0 ; j<loopNodes[i].size() ; j++ ) poly.indices[j] = loopNodes[i][j].index , poly.vertices[j] = loopNodes[i][j].position;
			poly.vertexIndices = loopAtlasVertexIndices[i];
			poly.atlasEdgeIndices = loopAtlasEdges[i];
			poly.atlasVertexParentEdge = loopAtlasVertexParentEdges[i];
			boundaryPolygons[ cellID ].emplace_back( loopKeys[i].first , poly );
		}
	}
	return boundaryPolygons;
}


template< typename GeometryReal , typename MatrixReal >
ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartBoundaryCellIndex , std::vector< std::pair< ChartMeshTriangleIndex , IndexedPolygon< GeometryReal > > > > >
GetBoundaryPolygons
(
	GridAtlas< GeometryReal , MatrixReal > &gridAtlas ,
	const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > &atlasCharts ,
	const typename AtlasChart< GeometryReal >::AtlasInfo &atlasInfo
)
{
	ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartBoundaryCellIndex , std::vector< std::pair< ChartMeshTriangleIndex , IndexedPolygon< GeometryReal > > > > > boundaryPolygons( gridAtlas.gridCharts.size() );

	ExplicitIndexVector< ChartIndex , GridChart< GeometryReal > > &gridCharts = gridAtlas.gridCharts;
	gridAtlas.endBoundaryVertexIndex = AtlasRefinedBoundaryVertexIndex( atlasInfo.atlasMeshVertexToBoundaryVertex.size() );

	// A map taking (atlas) boundary half-edge indices to their decomposition by grid edges
	std::map< AtlasMeshHalfEdgeIndex , std::vector< IntersectionInfo< GeometryReal > > > atlasBoundaryHalfEdgeToIntersectionInfos;

	// Maps taking a segment to the associated segment information
	std::vector< std::map< SimplexIndex< 1 , AtlasRefinedBoundaryVertexIndex > , BoundarySegmentInfo< GeometryReal > > > segmentToBoundarySegmentInfo( gridCharts.size() );

	// Maps taking an intersection key to the associated node information
	std::vector< std::map< GridMeshIntersectionKey , NodeInfo< GeometryReal , AtlasRefinedBoundaryVertexIndex > > > gridMeshIntersectionKeyToNodeInfo( gridCharts.size() );

	for( unsigned int i=0 ; i<gridCharts.size() ; i++ )
	{
		// -- Split the chart's boundary half-edges to the edges of the grid
		// -- Add associated auxiliary nodes
		// -- Add associated segments
		InitializeChartBoundaryEdgeGridIntersections
		(
			atlasCharts[ ChartIndex(i) ] ,
			atlasInfo.atlasMeshVertexToBoundaryVertex ,
			gridCharts[ ChartIndex(i) ] ,
			gridAtlas.endCoveredTexelIndex ,
			gridAtlas.endBoundaryVertexIndex ,
			atlasBoundaryHalfEdgeToIntersectionInfos ,
			segmentToBoundarySegmentInfo[i] ,
			gridMeshIntersectionKeyToNodeInfo[i]
		);

		// -- Split boundary-cell-inersecting triangles to the edges of the grid
		// -- Add associated auxiliary nodes
		InitializeChartBoundaryTriangleGridIntersections
		(
			atlasCharts[ ChartIndex(i) ] ,
			atlasInfo.atlasMeshVertexToBoundaryVertex ,
			gridCharts[ ChartIndex(i) ] ,
			gridAtlas.endCoveredTexelIndex ,
			gridAtlas.endBoundaryVertexIndex ,
			gridMeshIntersectionKeyToNodeInfo[i]
		);
	}

	for( unsigned int i=0 ; i<gridCharts.size() ; i++ )
	{
		try
		{
			boundaryPolygons[ ChartIndex(i) ] = GetChartBoundaryPolygons
			(
				atlasInfo.oppositeHalfEdges ,
				atlasCharts[ ChartIndex(i) ] ,
				gridCharts[ ChartIndex(i) ] ,
				gridAtlas.endCoveredTexelIndex ,
				(unsigned int)atlasInfo.atlasMeshVertexToBoundaryVertex.size() ,
#ifdef SANITY_CHECK
				gridAtlas.endBoundaryVertexIndex ,
#endif // SANITY_CHECK
				atlasBoundaryHalfEdgeToIntersectionInfos ,
				segmentToBoundarySegmentInfo[i] ,
				gridMeshIntersectionKeyToNodeInfo[i]
			);
		}
		catch( const Exception & e )
		{
			std::cerr << e.what() << std::endl;
			MK_THROW( "While processing chart: " , i );
		}
	}
	return boundaryPolygons;
}

template< typename GeometryReal >
void InitializeChartQuadraticElements
(
	GridChart< GeometryReal > &gridChart ,
	std::map< SimplexIndex< 1 > , BoundaryMidPointIndex > &midPointMap ,
	BoundaryMidPointIndex & endMidPointIndex ,
	const ExplicitIndexVector< ChartBoundaryCellIndex , std::vector< std::pair< ChartMeshTriangleIndex , IndexedPolygon< GeometryReal > > > > & boundaryPolygons ,
	unsigned int previouslyAddedNodes
)
{
	ExplicitIndexVector< ChartBoundaryCellIndex , std::vector< BoundaryIndexedTriangle< GeometryReal > > > & boundaryTriangles = gridChart.boundaryTriangles;
	boundaryTriangles.resize( gridChart.numBoundaryCells() );
	gridChart.endBoundaryTriangleIndex = ChartBoundaryTriangleIndex(0);
	for( unsigned int i=0 ; i<boundaryPolygons.size() ; i++ )
	{
		for( unsigned int j=0 ; j<boundaryPolygons[ ChartBoundaryCellIndex(i) ].size() ; j++ )
		{
			const IndexedPolygon< GeometryReal > &currentPolygon = boundaryPolygons[ ChartBoundaryCellIndex(i) ][j].second;
			std::vector< SimplexIndex< 2 > > polyTriangles;
			TriangulateConvexPolygon( currentPolygon.vertices , polyTriangles );
			if( polyTriangles.size()!=currentPolygon.vertices.size()-2 )
			{
				unsigned int localCellPos[] = { static_cast< unsigned int >(-1) , static_cast< unsigned int >(-1) };
				for( unsigned int li=0 ; li<gridChart.cellIndices.res(0) ; li++ ) for( unsigned int lj=0 ; lj<gridChart.cellIndices.res(1) ; lj++ ) if( gridChart.cellIndices(li,lj).boundary==static_cast< ChartBoundaryCellIndex >(i) ) localCellPos[0] = li , localCellPos[1] = lj;
				MK_WARN( "Unexpected number of triangles produced by convex polygon triangulation at combined cell ( " , gridChart.cornerCoords[0] + localCellPos[0] , " , " , gridChart.cornerCoords[1] + localCellPos[1] , " ). Polygon may self intersect! " , polyTriangles.size() , " != " , currentPolygon.vertices.size()-2 );
			}

			for( unsigned int k=0 ; k<polyTriangles.size() ; k++ )
			{
				BoundaryIndexedTriangle< GeometryReal > triangle;
				for( unsigned int v=0 ; v<3 ; v++ )
				{
					unsigned int currentPolygonVertex = polyTriangles[k][v];
					triangle.indices[v] = currentPolygon.indices[ currentPolygonVertex ];
					triangle.vertices[v] = currentPolygon.vertices[ currentPolygonVertex ];
					triangle.vertexIndices[v] = currentPolygon.vertexIndices[ currentPolygonVertex ];
					triangle.atlasVertexParentEdge[v] = currentPolygon.atlasVertexParentEdge[ currentPolygonVertex ];
					triangle.atlasEdgeIndices[v] = AtlasMeshEdgeIndex(-1);
					unsigned int nextPolygonVertex = (currentPolygonVertex + 1) % currentPolygon.vertices.size();
					if( polyTriangles[k][(v + 1) % 3] == nextPolygonVertex )
					{
						//preserve orientation
						triangle.atlasEdgeIndices[v] = currentPolygon.atlasEdgeIndices[currentPolygonVertex];
					}
					unsigned int previousPolygonVertex = (unsigned int)( (currentPolygonVertex + currentPolygon.vertices.size() - 1) % currentPolygon.vertices.size() );
					if( polyTriangles[k][(v + 1) % 3] == previousPolygonVertex )
					{
						//reverse orientation
						triangle.atlasEdgeIndices[v] = currentPolygon.atlasEdgeIndices[previousPolygonVertex];
					}

				}
				AtlasInteriorOrBoundaryNodeIndex midVertexIndex[3];

				for( unsigned int v=0 ; v<3 ; v++ )
				{
					Point2D< GeometryReal > edgeMidPoint = ( triangle.vertices[ (v+1)%3 ] + triangle.vertices[ (v+2)%3 ] ) / 2;
					AtlasInteriorOrBoundaryNodeIndex edgeCorners[2] = { triangle.indices[ (v+1)%3 ], triangle.indices[ (v+2)%3 ] };
					if( edgeCorners[1]<edgeCorners[0] ) std::swap( edgeCorners[0] , edgeCorners[1] );
					SimplexIndex< 1 > edgeKey( edgeCorners[0] , edgeCorners[1] );

					auto iter = midPointMap.find(edgeKey);
					if( iter==midPointMap.end() )
					{
						midPointMap[edgeKey] = endMidPointIndex;
						midVertexIndex[v] = AtlasInteriorOrBoundaryNodeIndex( previouslyAddedNodes + static_cast< unsigned int >(endMidPointIndex) );
						endMidPointIndex++;
					}
					else midVertexIndex[v] = AtlasInteriorOrBoundaryNodeIndex( previouslyAddedNodes + static_cast< unsigned int >(iter->second) );

					gridChart.auxiliaryNodes.emplace_back( midVertexIndex[v] , edgeMidPoint );
				}
				triangle.idx = gridChart.endBoundaryTriangleIndex++;
				triangle.sourceIdx = boundaryPolygons[ ChartBoundaryCellIndex(i) ][j].first;
				triangle.indices[3] = midVertexIndex[0];
				triangle.indices[4] = midVertexIndex[1];
				triangle.indices[5] = midVertexIndex[2];
				boundaryTriangles[ ChartBoundaryCellIndex(i) ].push_back( triangle );
			}
		}
	}
}

template< typename GeometryReal , typename MatrixReal >
void InitializeBoundaryTriangulation
(
	GridAtlas< GeometryReal , MatrixReal > &gridAtlas ,
	const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > &atlasCharts ,
	typename AtlasChart< GeometryReal >::AtlasInfo &atlasInfo
)
{
	// Add the vertices
	ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartBoundaryCellIndex , std::vector< std::pair< ChartMeshTriangleIndex , IndexedPolygon< GeometryReal > > > > > boundaryPolygons = GetBoundaryPolygons( gridAtlas , atlasCharts , atlasInfo );

	// Add the (mid-edge) nodes, fuse, and create triangles
	std::map< SimplexIndex< 1 > , BoundaryMidPointIndex > midPointMap;

	gridAtlas.endMidPointIndex = BoundaryMidPointIndex(0);
	for( unsigned int i=0 ; i<gridAtlas.gridCharts.size() ; i++ ) InitializeChartQuadraticElements( gridAtlas.gridCharts[ ChartIndex(i) ] , midPointMap , gridAtlas.endMidPointIndex , boundaryPolygons[ ChartIndex(i) ] , static_cast< unsigned int >(gridAtlas.endBoundaryVertexIndex) + static_cast< unsigned int >(gridAtlas.endCoveredTexelIndex) );
	gridAtlas.numFineNodes = static_cast< unsigned int >( gridAtlas.endCoveredTexelIndex ) + static_cast< unsigned int >(gridAtlas.endBoundaryVertexIndex) + static_cast< unsigned int >(gridAtlas.endMidPointIndex);
}

template< typename MatrixReal >
struct BoundaryProlongationData
{
	std::vector< AtlasInteriorOrBoundaryNodeIndex > fineBoundaryIndex;
	unsigned int numFineBoundaryNodes;
	SparseMatrix< MatrixReal , int > coarseBoundaryFineBoundaryProlongation;
	SparseMatrix< MatrixReal , int > fineBoundaryCoarseBoundaryRestriction;
};

template< typename GeometryReal , typename MatrixReal >
void InitializeCoarseBoundaryToFineBoundaryProlongation
(
	const GridAtlas< GeometryReal , MatrixReal > &gridAtlas ,
	SparseMatrix< MatrixReal , int > &coarseBoundaryFineBoundaryProlongation ,
	std::vector< AtlasInteriorOrBoundaryNodeIndex > &fineBoundaryIndex ,
	unsigned int &numFineBoundaryNodes ,
	bool verbose=false
)
{
	std::vector< unsigned int > boundaryFineToFullFine;
	std::vector< Eigen::Triplet< MatrixReal > > prolongationTriplets;
	const typename GridAtlas<>::IndexConverter & indexConverter = gridAtlas.indexConverter;

	fineBoundaryIndex.resize( gridAtlas.numFineNodes , AtlasInteriorOrBoundaryNodeIndex(-1) );
	numFineBoundaryNodes = 0;
	for( unsigned int i=0 ; i<gridAtlas.gridCharts.size() ; i++ )
	{
		const GridChart< GeometryReal > &gridChart = gridAtlas.gridCharts[ ChartIndex(i) ];
		for( unsigned int j=0 ; j<gridChart.texelIndices.size() ; j++ )
		{
			if( gridChart.texelIndices[j].covered!=AtlasCoveredTexelIndex(-1) && gridChart.texelIndices[j].interior==AtlasInteriorTexelIndex(-1) )
			{
				//////////////////////////
				// Boundary and covered //
				//////////////////////////

				AtlasTexelIndex coarseCombinedIndex = gridChart.texelIndices[j].combined;

				AtlasBoundaryTexelIndex boundaryIndex = indexConverter.combinedToBoundary( coarseCombinedIndex );
				if( boundaryIndex==AtlasBoundaryTexelIndex(-1) ) MK_THROW( "Coarse node is not boundary. Combined index " , coarseCombinedIndex , ". Boundary index " , boundaryIndex );
				prolongationTriplets.emplace_back( numFineBoundaryNodes , static_cast< unsigned int >(boundaryIndex) , (MatrixReal)1. );
				fineBoundaryIndex[ static_cast< unsigned int >(gridChart.texelIndices[j].covered) ] = static_cast< AtlasInteriorOrBoundaryNodeIndex >( numFineBoundaryNodes );
				numFineBoundaryNodes++;
				boundaryFineToFullFine.push_back( static_cast< unsigned int >(gridChart.texelIndices[j].covered) );
			}
		}
	}
	for( unsigned int i=static_cast< unsigned int >(gridAtlas.endCoveredTexelIndex) ; i<gridAtlas.numFineNodes ; i++ )
	{
		fineBoundaryIndex[i] = AtlasInteriorOrBoundaryNodeIndex( numFineBoundaryNodes++ );
		boundaryFineToFullFine.push_back(i);
	}

	if( verbose ) printf( "Fine boundary elements %d\n" , numFineBoundaryNodes );

	const ExplicitIndexVector< ChartIndex , GridChart< GeometryReal > > &gridCharts = gridAtlas.gridCharts;
	const ExplicitIndexVector< AtlasTexelIndex , TexelInfo > & texelInfo = gridAtlas.texelInfo;

	unsigned int numFineNodes = gridAtlas.numFineNodes;

	std::vector< unsigned int > auxiliaryNodesDegree( numFineNodes - static_cast< unsigned int >(gridAtlas.endCoveredTexelIndex) , 0);

	for( unsigned int i=0 ; i<gridCharts.size() ; i++ )
	{
		const GridChart< GeometryReal > &gridChart = gridCharts[ ChartIndex(i) ];
		for( unsigned int j=0 ; j<gridChart.auxiliaryNodes.size(); j++ )
		{
			unsigned int auxiliaryID = static_cast< unsigned int >( gridChart.auxiliaryNodes[j].index ) - static_cast< unsigned int >( gridAtlas.endCoveredTexelIndex );
			auxiliaryNodesDegree[auxiliaryID]++;
		}
	}

	GeometryReal precision_error = (GeometryReal)1e-10;

	std::vector< MatrixReal > auxiliaryNodesCumWeight( numFineNodes - static_cast< unsigned int >(gridAtlas.endCoveredTexelIndex) , 0 );

	for( unsigned int i=0 ; i<gridCharts.size() ; i++ )
	{
		const GridChart< GeometryReal > &gridChart = gridCharts[ ChartIndex(i) ];
		for( unsigned int j=0 ; j<gridChart.auxiliaryNodes.size() ; j++ )
		{
			unsigned int auxiliaryID = static_cast< unsigned int >(gridChart.auxiliaryNodes[j].index) - static_cast< unsigned int >( gridAtlas.endCoveredTexelIndex);
			AtlasInteriorOrBoundaryNodeIndex fineBoundaryID = fineBoundaryIndex[ static_cast< unsigned int >(gridChart.auxiliaryNodes[j].index) ];
			unsigned int nodeDegree = auxiliaryNodesDegree[auxiliaryID];
			Point2D< GeometryReal >nodePosition = gridChart.auxiliaryNodes[j].position;
			int corner[2] = { (int)floor(nodePosition[0] / gridChart.cellSizeW), (int)floor(nodePosition[1] / gridChart.cellSizeH) };
			ChartCellIndex cellID = gridChart.cellIndices( corner[0] , corner[1] ).combined;
			if( cellID==ChartCellIndex(-1) ) MK_THROW( "Invalid cell index. Node position " , nodePosition[0] / gridChart.cellSizeW , " " , nodePosition[1] / gridChart.cellSizeH );
			nodePosition[0] /= gridChart.cellSizeW;
			nodePosition[1] /= gridChart.cellSizeH;
			nodePosition[0] -= (GeometryReal)corner[0];
			nodePosition[1] -= (GeometryReal)corner[1];
			if( nodePosition[0] < 0-precision_error || nodePosition[0] > 1+precision_error || nodePosition[1] < 0-precision_error || nodePosition[1] > 1+precision_error )
				MK_THROW( "Sample out of unit box! (" , nodePosition[0] , " " , nodePosition[0] , ")" );
			for( unsigned int k=0 ; k<4 ; k++ )
			{
				MatrixReal texelWeight = (MatrixReal)BilinearElementValue(k, nodePosition) / nodeDegree;
				if( fabs(texelWeight)>1e-11 )
				{
					auxiliaryNodesCumWeight[auxiliaryID] += texelWeight;
					AtlasTexelIndex texelIndex = gridChart.combinedCellCombinedTexelBilinearElementIndices[cellID][k];
					if( texelInfo[texelIndex].texelType==TexelType::InteriorSupported )
						MK_THROW( "Interior-supported texel cannot be in the support of an auxiliary node. Weight " , texelWeight , " (A)" );

					if( static_cast< unsigned int >(gridChart.auxiliaryNodes[j].index)<static_cast< unsigned int >(gridAtlas.endCoveredTexelIndex) || static_cast< unsigned int >(gridChart.auxiliaryNodes[j].index)>numFineNodes || texelIndex==AtlasTexelIndex(-1) || static_cast< unsigned int >(texelIndex)>static_cast< unsigned int >(gridAtlas.endCombinedTexelIndex) ) MK_THROW( "Out of bounds index" );

					AtlasBoundaryTexelIndex boundaryIndex = indexConverter.combinedToBoundary( AtlasTexelIndex(texelIndex) );
					if( boundaryIndex==AtlasBoundaryTexelIndex(-1) ) MK_THROW( "Coarse node is not boundary" );
					prolongationTriplets.emplace_back( static_cast< unsigned int >(fineBoundaryID) , static_cast< unsigned int >(boundaryIndex) , texelWeight );
				}
			}
		}
	}

	coarseBoundaryFineBoundaryProlongation = SetSparseMatrix( prolongationTriplets , numFineBoundaryNodes , static_cast< unsigned int >(gridAtlas.endBoundaryTexelIndex) , false );
}

template< typename GeometryReal , typename MatrixReal >
void InitializeBoundaryProlongationData
(
	const GridAtlas< GeometryReal , MatrixReal > &gridAtlas ,
	BoundaryProlongationData< MatrixReal > &boundaryProlongation
)
{
	InitializeCoarseBoundaryToFineBoundaryProlongation( gridAtlas , boundaryProlongation.coarseBoundaryFineBoundaryProlongation , boundaryProlongation.fineBoundaryIndex , boundaryProlongation.numFineBoundaryNodes );
	boundaryProlongation.fineBoundaryCoarseBoundaryRestriction = boundaryProlongation.coarseBoundaryFineBoundaryProlongation.transpose();
}