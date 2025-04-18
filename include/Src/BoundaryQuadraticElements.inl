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
void SetAtlasIndexedPolygonFromBoundaryTriangle( const BoundaryIndexedTriangle< GeometryReal > &triangle , AtlasIndexedPolygon< GeometryReal > &polygon )
{
	for( int k=0 ; k<3 ; k++ )
	{
		polygon.vertices.push_back( triangle.vertices[k] );
		polygon.indices.push_back( triangle.indices[k] );
		polygon.atlasVertexIndices.push_back( triangle.atlasVertexIndices[k] );
		polygon.atlasEdgeIndices.push_back( triangle.atlasEdgeIndices[k] );
		polygon.atlasVertexParentEdge.push_back( triangle.atlasVertexParentEdge[k] );
	}
}

template< typename GeometryReal >
void IndexedPolygonFromCell( int i , int j , const GridChart< GeometryReal > &gridChart , IndexedIntersectionPolygon< GeometryReal > &polygon )
{
	polygon.vertices.resize(4);
	polygon.indices.resize(4);
	polygon.edgeIndices.resize(4);

	polygon.vertices[0] = gridChart.nodePosition( i+0 , j+0 );
	polygon.vertices[1] = gridChart.nodePosition( i+1 , j+0 );
	polygon.vertices[2] = gridChart.nodePosition( i+1 , j+1 );
	polygon.vertices[3] = gridChart.nodePosition( i+0 , j+1 );

	polygon.indices[0] = GridMeshIntersectionKey::GridKey( gridChart.nodeIndex( i+0 , j+0 ) );
	polygon.indices[1] = GridMeshIntersectionKey::GridKey( gridChart.nodeIndex( i+1 , j+0 ) );
	polygon.indices[2] = GridMeshIntersectionKey::GridKey( gridChart.nodeIndex( i+1 , j+1 ) );
	polygon.indices[3] = GridMeshIntersectionKey::GridKey( gridChart.nodeIndex( i+0 , j+1 ) );

	polygon.edgeIndices[0] = gridChart.edgeIndex( i+0 , j+0 , 1 );
	polygon.edgeIndices[1] = gridChart.edgeIndex( i+1 , j+0 , 0 );
	polygon.edgeIndices[2] = gridChart.edgeIndex( i+0 , j+1 , 1 );
	polygon.edgeIndices[3] = gridChart.edgeIndex( i+0 , j+0 , 0 );
}

// Computes the intersections of the edge with the grid lines and adds to the vector of IntersectionInfo
template< typename GeometryReal >
void AddEdgeGridIntersection
(
	Point2D< GeometryReal > edge[2] ,
	const GridChart< GeometryReal > &gridChart ,
	int edgeIndex ,
	std::vector< IntersectionInfo< GeometryReal > > &intersections
)
{
	GeometryReal fmin[2] , fmax[2];
	for( int dir=0 ; dir<2 ; dir++ )
	{
		fmin[dir] = std::max< GeometryReal >( 0.f , std::min< GeometryReal >( edge[0][dir] , edge[1][dir] ) );
		fmax[dir] = std::min< GeometryReal >( 1.f , std::max< GeometryReal >( edge[0][dir] , edge[1][dir] ) );
	}

	int imin[2] , imax[2];
	imin[0] = static_cast< int >( floor( fmin[0] / gridChart.cellSizeW ) );
	imin[1] = static_cast< int >( floor( fmin[1] / gridChart.cellSizeH ) );
	imax[0] = static_cast< int >( floor( fmax[0] / gridChart.cellSizeW ) );
	imax[1] = static_cast< int >( floor( fmax[1] / gridChart.cellSizeH ) );

	for( int dir=0 ; dir<2 ; dir++ ) for( int s=imin[dir] ; s<=imax[dir] ; s++ )
	{
		GeometryReal cellSize = dir == 0 ? gridChart.cellSizeW : gridChart.cellSizeH;
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

			info.intersectionKey = GridMeshIntersectionKey( gridChart.edgeIndex( cellIntersectionIndices[0] , cellIntersectionIndices[1] , dir ) , edgeIndex );
			info.time = alpha;
			info.position = position;

			intersections.push_back( info );
		}
	}
}

template< typename GeometryReal >
void InitializeChartBoundaryEdgeGridIntersections
(
	const AtlasChart< GeometryReal > & atlasChart ,
	const std::map< unsigned int , unsigned int > & surfaceBoundaryVertexToIndex ,
	GridChart< GeometryReal > & gridChart ,
	unsigned int & lastBoundaryIndex ,
	unsigned int numInteriorNodes ,
	std::unordered_map< unsigned int, std::vector< IntersectionInfo< GeometryReal > > > &boundaryEdgeIntersections ,
	std::map< EdgeIndex , BoundarySegmentInfo< GeometryReal > > &localBoundarySegmentsInfo ,
	std::map< GridMeshIntersectionKey , NodeInfo< GeometryReal > > & boundaryNodes 	// IntersectionKey -> boundary node index/position
)
{
	for( int b=0 ; b<atlasChart.boundaryHalfEdges.size() ; b++ )
	{
		unsigned int chartHalfEdge = atlasChart.boundaryHalfEdges[b];
		unsigned int tIndex = chartHalfEdge / 3;
		unsigned int kIndex = chartHalfEdge % 3;

		unsigned int atlasHalfEdge = atlasChart.atlasHalfEdge( chartHalfEdge );

		std::vector< IntersectionInfo< GeometryReal > > edgeIntersectionsInfo;
		Point2D< GeometryReal > edge[2];
		EdgeIndex eIndex = atlasChart.edgeIndex( atlasChart.boundaryHalfEdges[b] );
		edge[0] = atlasChart.vertices[ eIndex[0] ] - gridChart.corner;
		edge[1] = atlasChart.vertices[ eIndex[1] ] - gridChart.corner;

		// Add the end-points of the edge
		for( unsigned int v=0 ; v<2 ; v++ )
		{
			GridMeshIntersectionKey vertexKey = GridMeshIntersectionKey::MeshKey( eIndex[v] );

			IntersectionInfo< GeometryReal > cornerIntersection;
			cornerIntersection.intersectionKey = vertexKey;
			cornerIntersection.time = (GeometryReal)v;
			cornerIntersection.position = edge[v];
			edgeIntersectionsInfo.push_back( cornerIntersection );
		}

		// Add the points where the edge intersects the grid
		AddEdgeGridIntersection( edge , gridChart , atlasChart.atlasEdge( chartHalfEdge ) , edgeIntersectionsInfo );

		// Sort intersections
		std::sort( edgeIntersectionsInfo.begin() , edgeIntersectionsInfo.end() , IntersectionInfo< GeometryReal >::CompareByTime );

		for( int i=0 ; i<edgeIntersectionsInfo.size() ; i++ )
		{
			// If we haven't seen this vertex before...
			if( boundaryNodes.find( edgeIntersectionsInfo[i].intersectionKey )==boundaryNodes.end() )
			{
				int index;
				if( edgeIntersectionsInfo[i].intersectionKey.isMeshVertex() )	// Start/end vertex
				{
					auto iter = surfaceBoundaryVertexToIndex.find( atlasChart.surfaceVertex( edgeIntersectionsInfo[i].intersectionKey.meshIndex ) );
					if( iter==surfaceBoundaryVertexToIndex.end() ) MK_THROW( "Boundary vertex not found" );
					index = iter->second;
				}
				else index = lastBoundaryIndex++;

				boundaryNodes[ edgeIntersectionsInfo[i].intersectionKey ] = NodeInfo< GeometryReal >( index , edgeIntersectionsInfo[i].position );

				AuxiliaryNode< GeometryReal > auxNode;
				auxNode.index = numInteriorNodes + index;
				auxNode.position = edgeIntersectionsInfo[i].position;
				gridChart.auxiliaryNodes.push_back( auxNode );
			}
			edgeIntersectionsInfo[i].intersectionIndex = boundaryNodes[ edgeIntersectionsInfo[i].intersectionKey ].index;
		}

		for( int i=0 ; i<edgeIntersectionsInfo.size()-1 ; i++ )
		{
			unsigned int segmentCornerIndices[] = { boundaryNodes[ edgeIntersectionsInfo[i].intersectionKey ].index , boundaryNodes[ edgeIntersectionsInfo[i+1].intersectionKey ].index };
			EdgeIndex eIndex( segmentCornerIndices[0] , segmentCornerIndices[1] );
			BoundarySegmentInfo< GeometryReal > segmentInfo;
			segmentInfo.startTime = edgeIntersectionsInfo[i].time;
			segmentInfo.endTime = edgeIntersectionsInfo[i+1].time;
			segmentInfo.chartHalfEdge = chartHalfEdge;

			if( localBoundarySegmentsInfo.find( eIndex )!=localBoundarySegmentsInfo.end() ) MK_THROW( "Replicated segment key" );
			localBoundarySegmentsInfo[ eIndex ] = segmentInfo;
		}
		boundaryEdgeIntersections[ atlasHalfEdge ] = edgeIntersectionsInfo;
	}
}

template< typename GeometryReal >
void InitializeChartBoundaryPolygons
(
	const std::vector< unsigned int > &oppositeHalfEdge ,
	const AtlasChart< GeometryReal > &atlasChart ,
	GridChart< GeometryReal > &gridChart ,
	unsigned int numInteriorNodes ,
	unsigned int numBoundaryVertices ,
	unsigned int numBoundaryNodes ,
	std::unordered_map< unsigned int , std::vector< IntersectionInfo< GeometryReal > > > &boundaryEdgeIntersections ,
	std::map< EdgeIndex , BoundarySegmentInfo< GeometryReal > > &boundarySegmentsInfo ,
	std::map< GridMeshIntersectionKey , NodeInfo< GeometryReal > > & boundaryNodes ,
	std::vector< unsigned int > &coveredOppositeBoundaryNode ,
	unsigned int chartID
)
{
#ifdef SEPARATE_POLYGONS
	// A mapping giving the clipped polygons, keyed off of cell index
	std::map< int , std::vector< std::vector< unsigned long long > > > cellPolygons;
#else // !SEPARATE_POLYGONS
	// An association of segments to individual cells of a chart
	std::unordered_map< int , std::pair< Point< int , 2 > , std::vector< std::pair< GridMeshIntersectionKey , GridMeshIntersectionKey > > > > cellSegments;
#endif // SEPARATE_POLYGONS

	auto GetIndexedTriangle = [&]( unsigned int t )
		{
			IndexedIntersectionTriangle< GeometryReal > indexedTriangle;
			Point2D< GeometryReal > tPos[3];
			for( unsigned int i=0 ; i<3 ; i++ ) tPos[i] = atlasChart.vertices[ atlasChart.triangles[t][i] ] - gridChart.corner;

			for( unsigned int k=0 ; k<3 ; k++ )
			{
				indexedTriangle.vertices[k] = tPos[k];
				indexedTriangle.indices[k] = GridMeshIntersectionKey::MeshKey( atlasChart.triangles[t][k] );
				indexedTriangle.edgeIndices[k] = atlasChart.atlasEdge( 3*t+k );
			}
			return indexedTriangle;
		};

	auto GetIndexedIntersectionPolygon = [&]( unsigned int i , unsigned int j )
		{
			IndexedIntersectionPolygon< GeometryReal > cellPolygon;
			IndexedPolygonFromCell( i , j , gridChart , cellPolygon );
			return cellPolygon;
		};

#ifdef USE_RASTERIZER
	using Range = RegularGrid< 2 >::Range;
	using Index = RegularGrid< 2 >::Index;
	Range nodeRange , cellRange;
	nodeRange.second[0] = gridChart.width  , cellRange.second[0] = gridChart.width-1;
	nodeRange.second[1] = gridChart.height , cellRange.second[1] = gridChart.height-1;

	auto GetSimplex = [&]( unsigned int t )
		{
			Simplex< double , 2 , 2 > simplex;
			for( unsigned int i=0 ; i<=2 ; i++ )
			{
				simplex[i] = atlasChart.vertices[ atlasChart.triangles[t][i] ] - gridChart.corner;
				simplex[i][0] /= gridChart.cellSizeW , simplex[i][1] /= gridChart.cellSizeH;
			}
			return simplex;
		};
#endif // USE_RASTERIZER

	//(1) Rasterize all triangles and, if they fall into the support of a boundary cell, add them to the list of triangles associated with that cell
	for( int t=0 ; t<atlasChart.triangles.size() ; t++ )
	{
		IndexedIntersectionTriangle< GeometryReal > indexedTriangle = GetIndexedTriangle( t );
#ifdef USE_RASTERIZER
		// Compute the associated triangle in (shifted) texel coordinates
		Simplex< double , 2 , 2 > simplex = GetSimplex( t );

		auto Kernel = [&]( Index I )
			{
				if( gridChart.cellType(I)==CellType::Boundary )
				{
					unsigned int cellID = gridChart.cellIndices(I).boundary;
					if( cellID==-1 ) MK_THROW( "Boundary cell invalid ID" );
#ifdef SEPARATE_POLYGONS
					IndexedIntersectionPolygon< GeometryReal > poly = GetIndexedIntersectionPolygon( I[0] , I[1] );
					if( ClipIndexedIntersectionPolygonToIndexedIntersectionTriangle( poly , indexedTriangle ) ) cellPolygons[ cellID ].push_back( poly.indices );
					else MK_THROW( "Expected triangle to intersect cell" );
#else // !SEPARATE_POLYGONS
					cellSegments[cellID].first = Point< int , 2 >( I[0]+gridChart.cornerCoords[0] , I[1]+gridChart.cornerCoords[1] );

					IndexedIntersectionPolygon< GeometryReal > poly = GetIndexedIntersectionPolygon( I[0] , I[1] );

					if( ClipIndexedIntersectionPolygonToIndexedIntersectionTriangle( poly , indexedTriangle ) )
						for( int s=0 ; s<poly.indices.size() ; s++ )
							cellSegments[cellID].second.push_back( std::make_pair( poly.indices[s] , poly.indices[ (s+1) % poly.indices.size() ] ) );
					else MK_THROW( "Expected triangle to intersect cell" );
#endif // SEPARATE_POLYGONS
				}
			};
		Rasterizer2D::RasterizeSupports< true , true >( simplex , Kernel , cellRange );
#else // !USE_RASTERIZER
		int minCorner[2] , maxCorner[2];
//		GetTriangleIntegerBBox( indexedTriangle.vertices , (GeometryReal)1./gridChart.cellSizeW , (GeometryReal)1./gridChart.cellSizeH , minCorner , maxCorner );
		GetTriangleIntegerBBox( &indexedTriangle.vertices[0] , (GeometryReal)1./gridChart.cellSizeW , (GeometryReal)1./gridChart.cellSizeH , minCorner , maxCorner );

		for( int j=minCorner[1] ; j<maxCorner[1] ; j++ ) for( int i=minCorner[0] ; i<maxCorner[0] ; i++ )
		{
			if( gridChart.cellType(i,j)==CellType::Boundary )
			{
				int cellID = gridChart.cellIndices(i,j).boundary;
				if( cellID==-1 ) MK_THROW( "Boundary cell invalid ID" );
				{
					Point< int , 2 > idx;
					idx[0] = i+gridChart.cornerCoords[0];
					idx[1] = j+gridChart.cornerCoords[1];
					cellSegments[cellID].first = idx;
				}

				IndexedIntersectionPolygon< GeometryReal > cellPolygon = GetIndexedIntersectionPolygon( i , j );

				if( ClipIndexedIntersectionPolygonToIndexedIntersectionTriangle( cellPolygon , indexedTriangle ) )
				{
					for( int s=0 ; s<cellPolygon.indices.size() ; s++ )
					{
						{
							Point2D< GeometryReal > p[] = { cellPolygon.vertices[s] , cellPolygon.vertices[ (s+1)%cellPolygon.vertices.size() ] };
							if( Point2D< GeometryReal >::SquareNorm( p[0] - p[1] )<=1e-16 )
								MK_WARN( "Short clipped edge @ texel: " , cellSegments[cellID].first[0] , " , " , cellSegments[cellID].first[1] , " : " , sqrt( Point2D< GeometryReal >::SquareNorm( p[0] - p[1] ) ) );
						}
						std::pair< GridMeshIntersectionKey , GridMeshIntersectionKey > edge( cellPolygon.indices[s] , cellPolygon.indices[ (s+1) % cellPolygon.indices.size() ] );
						cellSegments[cellID].second.push_back( edge );
					}
				}
			}
		}
#endif // USE_RASTERIZER
	}

	//(2) Process cells
	gridChart.boundaryPolygons.resize( gridChart.numBoundaryCells() );
#ifdef SEPARATE_POLYGONS
	for( auto iter=cellPolygons.begin() ; iter!=cellPolygons.end() ; iter++ )
	{
		int cellID = iter->first;
		std::vector< std::vector< unsigned long long > > loopKeys = iter->second;
#else // !SEPARATE_POLYGONS
	for( auto cellIter=cellSegments.begin() ; cellIter!=cellSegments.end() ; cellIter++ )
	{
		int cellID = cellIter->first;
		std::vector< std::vector< GridMeshIntersectionKey > > loopKeys;

		// Get the polygons corresponding to the intersection of the cell with the the texture triangles
		{
			const std::vector< std::pair< GridMeshIntersectionKey , GridMeshIntersectionKey > > & segments = cellIter->second.second;

			// Extract the subset of boundary edges
			std::map< GridMeshIntersectionKey , GridMeshIntersectionKey > forwardMap;
			{
				// Get a set of all edges
				std::set< std::pair< GridMeshIntersectionKey , GridMeshIntersectionKey > > edgeSet;
				for( int k=0 ; k<segments.size() ; k++ ) edgeSet.insert( segments[k] );

				// Keep just the ones that don't have an opposite
				for( int k=0 ; k<segments.size() ; k++ )
				{
					std::pair< GridMeshIntersectionKey , GridMeshIntersectionKey > oppositeKey = std::make_pair( segments[k].second , segments[k].first );
					if( edgeSet.find( oppositeKey )==edgeSet.end() ) forwardMap[ segments[k].first ] = segments[k].second;
				}
			}

			// Transform the boundary edges to boundary loops
			try{ LoopVertices( forwardMap , loopKeys ); }
			catch( const Exception & e )
			{
				std::cerr << e.what() << std::endl;

				Point< int , 2 > idx = cellIter->second.first;
				{
					Point< int , 2 > _idx;
					_idx[0] = idx[0] - gridChart.cornerCoords[0];
					_idx[1] = idx[1] - gridChart.cornerCoords[1];
					std::stringstream sStream;
					for( unsigned int t=0 ; t<atlasChart.triangles.size() ; t++ )
					{
						IndexedIntersectionTriangle< GeometryReal > indexedTriangle = GetIndexedTriangle( t );
						IndexedIntersectionPolygon< GeometryReal > cellPolygon = GetIndexedIntersectionPolygon( _idx[0] , _idx[1] );
						if( ClipIndexedIntersectionPolygonToIndexedIntersectionTriangle( cellPolygon , indexedTriangle ) ) sStream << " " << atlasChart.atlasTriangle(t);
					}
					MK_WARN( "Overlapping triangles: " , sStream.str() );
				}
				MK_THROW( "While processing texel ( " , idx[0] , " , " , idx[1] , " ) -> " , cellID );
			}
		}
#endif // SEPARATE_POLYGONS

		// The position and index of boundary nodes/vertices
		std::vector< std::vector< NodeInfo< GeometryReal > > > loopNodes( loopKeys.size() );

		// For vertices on the original mesh, the index of the associated mesh vertex (-1 otherwise)
		std::vector< std::vector< int > > loopAtlasVertexIndices( loopKeys.size() );

		// Remove vertices generated when non-boundary edges intersect the grid
		for( int i=0 ; i<loopKeys.size() ; i++ )
		{
			std::vector< GridMeshIntersectionKey > reducedLoop;
			std::vector< NodeInfo< GeometryReal > > currentLoopNodes;
			std::vector< int > currentAtlasVertexIndices;

			// For each polygon vertex
			for( int j=0 ; j<loopKeys[i].size() ; j++ )
			{
				GridMeshIntersectionKey currentVertexKey = loopKeys[i][j];

				if( currentVertexKey.isGridVertex() )	// Texel node
				{
					unsigned int pi , pj;
					if( !gridChart.factorNodeIndex( currentVertexKey.gridIndex , pi , pj ) ) MK_ERROR_OUT( "Could not factor node index" );

					// Confirm that the texel is inside the chart
					int interiorTexelIndex = gridChart.texelIndices(pi,pj).interiorOrCovered;
					if( interiorTexelIndex==-1 ) MK_THROW( "Invalid texel: " , Point2D< int >( pi , pj ) , " -> " , gridChart.texelIndices(pi,pj).interiorOrCovered , " : " , gridChart.texelIndices(pi, pj).interior );

					reducedLoop.push_back( currentVertexKey );
					currentLoopNodes.emplace_back( interiorTexelIndex , gridChart.nodePosition(pi,pj) );
					currentAtlasVertexIndices.push_back(-1);
				}
				else if( boundaryNodes.find( currentVertexKey )!=boundaryNodes.end() )
				{
					reducedLoop.push_back( currentVertexKey );
					currentLoopNodes.push_back( boundaryNodes[currentVertexKey] );
					currentLoopNodes.back().index += numInteriorNodes;

					if( currentVertexKey.isMeshVertex() ) currentAtlasVertexIndices.push_back( currentVertexKey.meshIndex );
					else						          currentAtlasVertexIndices.push_back(-1);
				}
			}
			if( reducedLoop.size()==0 ) MK_THROW( "Reduced loop cannot be empty. Original loop size " , loopKeys[i].size() );
			loopKeys[i] = reducedLoop;
			loopNodes[i] = currentLoopNodes;
			loopAtlasVertexIndices[i] = currentAtlasVertexIndices;
		}

		// Add the intermediate vertices for boundary segments
		std::vector< std::vector< int > > loopAtlasEdges( loopKeys.size() );
		std::vector< std::vector< int > > loopAtlasVertexParentEdges( loopKeys.size() );
		for( int i=0 ; i<loopKeys.size() ; i++ )
		{
			std::vector< std::vector< int > > indicesToInsert;
			std::vector< std::vector< GeometryReal > > timesToInsert;
			std::vector< bool > isInsertionPos( loopKeys[i].size() , false );
			std::vector< int > polygonAtlasEdgeIndex( loopKeys[i].size() , -1 );
			std::vector< int > polygonAtlasVertexParentEdges( loopKeys[i].size() , -1 );	// The mesh edge the boundary vertex/node comes from

			for( int j=0 ; j<loopKeys[i].size() ; j++ )
			{
				GridMeshIntersectionKey currentVertexKey = loopKeys[i][j];
				GridMeshIntersectionKey nextVertexKey = loopKeys[i][(j + 1) % loopKeys[i].size()];


				if( currentVertexKey.isIntersection() ) polygonAtlasVertexParentEdges[j] = currentVertexKey.meshIndex;

				// If this is part of a boundary edge
				if( boundaryNodes.find( currentVertexKey )!=boundaryNodes.end() && boundaryNodes.find( nextVertexKey )!=boundaryNodes.end() )
				{
					EdgeIndex eIndex( boundaryNodes[currentVertexKey].index , boundaryNodes[nextVertexKey].index );
					if( boundarySegmentsInfo.find( eIndex )!=boundarySegmentsInfo.end() )
					{
						BoundarySegmentInfo< GeometryReal > segmentInfo = boundarySegmentsInfo[ eIndex ];
						polygonAtlasEdgeIndex[j] = atlasChart.atlasEdge( segmentInfo.chartHalfEdge );

						unsigned int oppHalfEdge = oppositeHalfEdge[ atlasChart.atlasHalfEdge( segmentInfo.chartHalfEdge ) ];
						if( oppHalfEdge!=-1 )
						{
							GeometryReal startTime = segmentInfo.startTime;
							GeometryReal endTime = segmentInfo.endTime;


							if( boundaryEdgeIntersections.find( oppHalfEdge )==boundaryEdgeIntersections.end() )
								MK_THROW( "Opposite edge intersections not found. Current  edge " , atlasChart.atlasHalfEdge( segmentInfo.chartHalfEdge ) , ". Opposite " , oppHalfEdge );

							std::vector< int > segmentIndicesToInsert;
							std::vector< GeometryReal > segmentTimesToInsert;

							std::vector< IntersectionInfo< GeometryReal > > oppositeIntersection = boundaryEdgeIntersections[ oppHalfEdge ];
							for( int oppIter = (int)oppositeIntersection.size()-2 ; oppIter>0 ; oppIter-- )
							{
								GeometryReal reversedTime = (GeometryReal)1. - oppositeIntersection[oppIter].time;
								GeometryReal normalizedTime = ( reversedTime-startTime ) / ( endTime-startTime );
#ifdef INSERTION_EPSILON
								if( normalizedTime>INSERTION_EPSILON && normalizedTime<(1.-INSERTION_EPSILON) )
#else // !INSERTION_EPSILON
								if (normalizedTime > 0 && normalizedTime < 1)
#endif // INSERTION_EPSILON
								{
									segmentIndicesToInsert.push_back( oppositeIntersection[oppIter].intersectionIndex );
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
				}
			}
			// Do insertions
			int insertionCount = 0;
			std::vector< NodeInfo< GeometryReal > > expandedLoopNodes;
			std::vector< int > expandedLoopAtlasVertexIndices;
			std::vector< int > expandedLoopAtlasEdgeIndex;
			std::vector< int > expandedVertexParentEdgeIndex;
			for( int j=0 ; j<loopKeys[i].size() ; j++ )
			{
				//expandedLoopKeys.push_back(loopKeys[i][j]);
				expandedLoopNodes.push_back( loopNodes[i][j] );
				expandedLoopAtlasVertexIndices.push_back(loopAtlasVertexIndices[i][j]);

				int currentSegmentAtlasEdgeIndex = polygonAtlasEdgeIndex[j];
				expandedLoopAtlasEdgeIndex.push_back(currentSegmentAtlasEdgeIndex);
				expandedVertexParentEdgeIndex.push_back(polygonAtlasVertexParentEdges[j]);

				if( isInsertionPos[j] )
				{
					std::vector< int > _indices = indicesToInsert[insertionCount];
					std::vector< GeometryReal > _times = timesToInsert[insertionCount];
					Point2D< GeometryReal > currentPos = loopNodes[i][j].position;
					Point2D< GeometryReal > nextPos = loopNodes[i][ (j+1)%loopKeys[i].size() ].position;

					if( currentSegmentAtlasEdgeIndex==-1 ) MK_THROW( "Invalid atlas edge index" );
					for( int k=0 ; k<_indices.size() ; k++ )
					{
						GeometryReal alpha = _times[k];
						Point2D< GeometryReal > interpolatedPos = (GeometryReal)( 1.-alpha )*currentPos + nextPos*alpha;
						if( true ) // Add orthogonal perturbation to avoid colinearity. Delanauy triangulation library failed for colinear inputs on some systems.
						{
							Point2D< GeometryReal > segmentDir = nextPos - currentPos;
							interpolatedPos += Point2D< GeometryReal >( -segmentDir[1] , segmentDir[0] ) * (GeometryReal)( ( Random< GeometryReal >()*2 - 1. ) * 1e-10 );
						}
						if( _indices[k]<numBoundaryVertices || _indices[k]>numBoundaryNodes ) MK_THROW( "Out of bounds index: " , _indices[k] , " not in " , numBoundaryVertices , " " , numBoundaryNodes );
						coveredOppositeBoundaryNode[_indices[k] - numBoundaryVertices] += 1;

						AuxiliaryNode< GeometryReal > auxNode;
						auxNode.index = _indices[k] + numInteriorNodes;
						auxNode.position = interpolatedPos;
						gridChart.auxiliaryNodes.push_back(auxNode);

						expandedLoopNodes.emplace_back( _indices[k] + numInteriorNodes , interpolatedPos );
						expandedLoopAtlasVertexIndices.push_back(-1);
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

		for( int i=0 ; i<loopNodes.size() ; i++ )
		{
			AtlasIndexedPolygon< GeometryReal > poly;
			poly.vertices.resize( loopNodes[i].size() ) , poly.indices.resize( loopNodes[i].size() );
			for( unsigned int j=0 ; j<loopNodes[i].size() ; j++ ) poly.indices[j] = loopNodes[i][j].index , poly.vertices[j] = loopNodes[i][j].position;
			poly.atlasVertexIndices = loopAtlasVertexIndices[i];
			poly.atlasEdgeIndices = loopAtlasEdges[i];
			poly.atlasVertexParentEdge = loopAtlasVertexParentEdges[i];
			gridChart.boundaryPolygons[ cellID ].push_back(poly);
		}
	}
}


template< typename GeometryReal >
void InitializeBoundaryPolygons
(
	const std::vector< AtlasChart< GeometryReal > > &atlasCharts ,
	const typename AtlasChart< GeometryReal >::AtlasInfo &atlasInfo ,
	std::vector< GridChart< GeometryReal > > &gridCharts ,
	unsigned int numInteriorNodes ,
	unsigned int &numBoundaryNodes
)
{ //Fine System

	unsigned int numBoundaryVertices = atlasInfo.surfaceBoundaryVertexToIndex.size();
	unsigned int lastBoundaryIndex = numBoundaryVertices;
	std::unordered_map< unsigned int , std::vector< IntersectionInfo< GeometryReal > > > boundaryEdgeIntersections;

	std::vector< std::map< EdgeIndex , BoundarySegmentInfo< GeometryReal > > > localBoundarySegmentsInfo( gridCharts.size() );
	std::vector< std::map< GridMeshIntersectionKey , NodeInfo< GeometryReal > > > boundaryNodes( gridCharts.size() );

	for( unsigned int i=0 ; i<gridCharts.size() ; i++ )
		InitializeChartBoundaryEdgeGridIntersections( atlasCharts[i] , atlasInfo.surfaceBoundaryVertexToIndex , gridCharts[i] , lastBoundaryIndex , numInteriorNodes , boundaryEdgeIntersections , localBoundarySegmentsInfo[i] , boundaryNodes[i] );	

	numBoundaryNodes = lastBoundaryIndex;
	std::vector< unsigned int > coveredOppositeBoundaryNode( numBoundaryNodes-numBoundaryVertices , 0 );
	for( unsigned int i=0 ; i<gridCharts.size() ; i++ )
	{
		try
		{
			InitializeChartBoundaryPolygons( atlasInfo.oppositeHalfEdges , atlasCharts[i] , gridCharts[i] , numInteriorNodes , numBoundaryVertices , numBoundaryNodes , boundaryEdgeIntersections , localBoundarySegmentsInfo[i] , boundaryNodes[i] , coveredOppositeBoundaryNode , i );
		}
		catch( const Exception & e )
		{
			std::cerr << e.what() << std::endl;
			MK_THROW( "While processing chart: " , i );
		}
	}

	if( atlasInfo.isClosed ) for ( int i=0 ; i<coveredOppositeBoundaryNode.size() ; i++ ) if( coveredOppositeBoundaryNode[i]!=1 ) MK_WARN( "Non-opposite boundary node at node " , i );
}

template< typename GeometryReal >
void InitializeChartQuadraticElements
(
	GridChart< GeometryReal > &gridChart ,
	std::map< EdgeIndex , unsigned int > &midPointIndex ,
	unsigned int &lastMidPointIndex ,
	unsigned int previouslyAddedNodes
)
{
	const std::vector< std::vector< AtlasIndexedPolygon< GeometryReal > > > & boundaryPolygons = gridChart.boundaryPolygons;
	std::vector< std::vector< BoundaryIndexedTriangle< GeometryReal > > > & boundaryTriangles = gridChart.boundaryTriangles;
	boundaryTriangles.resize( gridChart.numBoundaryCells() );
	int numBoundaryTriangles = 0;
	for( int i=0 ; i<boundaryPolygons.size() ; i++ )
	{
		for( int j=0 ; j<boundaryPolygons[i].size() ; j++ )
		{
			const AtlasIndexedPolygon< GeometryReal > &currentPolygon = boundaryPolygons[i][j];
			std::vector< SimplexIndex< 2 > > delanauyTriangles;
			TriangulatePolygon( currentPolygon.vertices , delanauyTriangles );
			if( delanauyTriangles.size()!=currentPolygon.vertices.size()-2 )
			{
				int localCellPos[2] = { -1,-1 };
				for( unsigned int li=0 ; li<gridChart.cellIndices.res(0) ; li++ ) for( unsigned int lj=0 ; lj<gridChart.cellIndices.res(1) ; lj++ ) if( gridChart.cellIndices(li,lj).boundary==i ) localCellPos[0] = li , localCellPos[1] = lj;
				MK_WARN( "Unexpected number of triangles produced by delaunay triangulation at global cell ( " , gridChart.cornerCoords[0] + localCellPos[0] , " , " , gridChart.cornerCoords[1] + localCellPos[1] , " ). Polygon may self intersect! " , delanauyTriangles.size() , " != " , currentPolygon.vertices.size()-2 );
			}

			for( int k=0 ; k<delanauyTriangles.size() ; k++ )
			{
				BoundaryIndexedTriangle< GeometryReal > triangle;
				for( int v=0 ; v<3 ; v++ )
				{
					int currentPolygonVertex = delanauyTriangles[k][v];
					triangle.indices[v] = currentPolygon.indices[currentPolygonVertex];
					triangle.vertices[v] = currentPolygon.vertices[currentPolygonVertex];
					triangle.atlasVertexIndices[v] = currentPolygon.atlasVertexIndices[currentPolygonVertex];
					triangle.atlasVertexParentEdge[v] = currentPolygon.atlasVertexParentEdge[currentPolygonVertex];
					triangle.atlasEdgeIndices[v] = -1;
					unsigned int nextPolygonVertex = (currentPolygonVertex + 1) % currentPolygon.vertices.size();
					if (delanauyTriangles[k][(v + 1) % 3] == nextPolygonVertex) { //preserve orientation
						triangle.atlasEdgeIndices[v] = currentPolygon.atlasEdgeIndices[currentPolygonVertex];
					}
					unsigned int previousPolygonVertex = (unsigned int)( (currentPolygonVertex + currentPolygon.vertices.size() - 1) % currentPolygon.vertices.size() );
					if( delanauyTriangles[k][(v + 1) % 3] == previousPolygonVertex) {//reverse orientation
						triangle.atlasEdgeIndices[v] = currentPolygon.atlasEdgeIndices[previousPolygonVertex];
					}

				}
				int midVexterIndex[3];

				for( int v=0 ; v<3 ; v++ )
				{
					Point2D< GeometryReal > edgeMidPoint = ( triangle.vertices[ (v+1)%3 ] + triangle.vertices[ (v+2)%3 ] ) / 2;
					int edgeCorners[2] = { (int)triangle.indices[(v + 1) % 3], (int)triangle.indices[(v + 2) % 3] };
					if (edgeCorners[0] > edgeCorners[1]) std::swap(edgeCorners[0], edgeCorners[1]);
					EdgeIndex edgeKey( edgeCorners[0] , edgeCorners[1] );

					if( midPointIndex.find(edgeKey)==midPointIndex.end() )
					{
						midPointIndex[edgeKey] = lastMidPointIndex;
						midVexterIndex[v] = lastMidPointIndex + previouslyAddedNodes;
						lastMidPointIndex++;
					}
					else midVexterIndex[v] = midPointIndex[edgeKey] + previouslyAddedNodes;

					AuxiliaryNode< GeometryReal > auxNode;
					auxNode.index = midVexterIndex[v];
					auxNode.position = edgeMidPoint;
					gridChart.auxiliaryNodes.push_back(auxNode);

				}
				triangle.id = numBoundaryTriangles++;
				triangle.indices[3] = midVexterIndex[0];
				triangle.indices[4] = midVexterIndex[1];
				triangle.indices[5] = midVexterIndex[2];
				boundaryTriangles[i].push_back( triangle );
			}
		}
	}
	gridChart.numBoundaryTriangles = numBoundaryTriangles;
}

template< typename GeometryReal >
void InitializeBoundaryQuadraticElements
(
	std::vector< GridChart< GeometryReal > > &gridCharts ,
	unsigned int previouslyAddedNodes ,
	unsigned int & numMidPoints )
{
	unsigned int lastMidPointIndex = 0;
	std::map< EdgeIndex , unsigned int > midPointIndex;
	for( unsigned int i=0 ; i<gridCharts.size() ; i++ ) InitializeChartQuadraticElements( gridCharts[i] , midPointIndex , lastMidPointIndex , previouslyAddedNodes );
	numMidPoints = lastMidPointIndex;
}

template< typename GeometryReal , typename MatrixReal >
void InitializeBoundaryTriangulation
(
	GridAtlas< GeometryReal , MatrixReal > &gridAtlas ,
	const std::vector< AtlasChart< GeometryReal > > &atlasCharts ,
	typename AtlasChart< GeometryReal >::AtlasInfo &atlasInfo
)
{
	InitializeBoundaryPolygons( atlasCharts , atlasInfo , gridAtlas.gridCharts , gridAtlas.numInteriorTexels , gridAtlas.numBoundaryNodes );
	InitializeBoundaryQuadraticElements( gridAtlas.gridCharts , gridAtlas.numBoundaryNodes + gridAtlas.numInteriorTexels , gridAtlas.numMidPoints );
	gridAtlas.numFineNodes = gridAtlas.numInteriorTexels + gridAtlas.numBoundaryNodes + gridAtlas.numMidPoints;
}

template< typename MatrixReal >
class BoundaryProlongationData
{
public:
	std::vector< unsigned int > fineBoundaryIndex;
	unsigned int numFineBoundaryNodes;
	SparseMatrix< MatrixReal , int > coarseBoundaryFineBoundaryProlongation;
	SparseMatrix< MatrixReal , int > fineBoundaryCoarseBoundaryRestriction;
};

template< typename GeometryReal , typename MatrixReal >
void InitializeCoarseBoundaryToFineBoundaryProlongation
(
	const GridAtlas< GeometryReal , MatrixReal > &gridAtlas ,
	SparseMatrix< MatrixReal , int > &coarseBoundaryFineBoundaryProlongation ,
	std::vector< unsigned int > &fineBoundaryIndex ,
	unsigned int &numFineBoundaryNodes ,
	bool verbose=false
)
{
	std::vector< unsigned int > boundaryFineToFullFine;
	std::vector< Eigen::Triplet< MatrixReal > > prolongationTriplets;
	const std::vector< unsigned int > & boundaryAndDeepIndex = gridAtlas.boundaryAndDeepIndex;

	fineBoundaryIndex.resize( gridAtlas.numFineNodes , static_cast< unsigned int >(-1) );
	unsigned int lastFineBoundaryIndex = 0;
	for (int i = 0; i < gridAtlas.gridCharts.size(); i++)
	{
		const GridChart< GeometryReal > &gridChart = gridAtlas.gridCharts[i];
		for( int j=0 ; j<gridChart.texelIndices.size() ; j++ )
		{
			if( gridChart.texelIndices[j].interiorOrCovered!=-1 && gridChart.texelIndices[j].interior==-1 )
			{ //Interior but not deep
				int coarseGlobalIndex = gridChart.texelIndices[j].combined;

				int coarseBondaryIndex = boundaryAndDeepIndex[coarseGlobalIndex] - 1;

				if( coarseBondaryIndex>=0 ) prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( lastFineBoundaryIndex , coarseBondaryIndex , (MatrixReal)1. ) );
				else MK_THROW( "Coarse node is not boundary. Global index " , coarseGlobalIndex , ". Boundary index " , coarseBondaryIndex );
				fineBoundaryIndex[ gridChart.texelIndices[j].interiorOrCovered ] = lastFineBoundaryIndex;
				lastFineBoundaryIndex++;
				boundaryFineToFullFine.push_back( gridChart.texelIndices[j].interiorOrCovered );
			}
		}
	}
	for( unsigned int i=gridAtlas.numInteriorTexels ; i<gridAtlas.numFineNodes ; i++ )
	{
		fineBoundaryIndex[i] = lastFineBoundaryIndex;
		lastFineBoundaryIndex++;
		boundaryFineToFullFine.push_back(i);
	}

	if( verbose ) printf( "Fine boundary elements %d\n" , lastFineBoundaryIndex );
	numFineBoundaryNodes = lastFineBoundaryIndex;

	const std::vector< GridChart< GeometryReal > > &gridCharts = gridAtlas.gridCharts;
	const std::vector<GridNodeInfo> & nodeInfo = gridAtlas.nodeInfo;

	int numInteriorTexels = gridAtlas.numInteriorTexels;
	int numFineNodes = gridAtlas.numFineNodes;
	int numCoarseNodes = gridAtlas.numTexels;

	int numAuxiliaryNodes = numFineNodes - numInteriorTexels;
	std::vector<int> auxiliaryNodesDegree(numAuxiliaryNodes, 0);

	for (int i = 0; i < gridCharts.size(); i++)
	{
		const GridChart< GeometryReal > &gridChart = gridCharts[i];
		for (int j = 0; j < gridChart.auxiliaryNodes.size(); j++) {
			int auxiliaryID = gridChart.auxiliaryNodes[j].index - numInteriorTexels;
			auxiliaryNodesDegree[auxiliaryID]++;
		}
	}

	GeometryReal precision_error = (GeometryReal)1e-10;

	std::vector< MatrixReal > auxiliaryNodesCumWeight( numAuxiliaryNodes , 0 );

	for( int i=0 ; i<gridCharts.size() ; i++ )
	{
		const GridChart< GeometryReal > &gridChart = gridCharts[i];
		for (int j = 0; j < gridChart.auxiliaryNodes.size(); j++) {
			int auxiliaryID = gridChart.auxiliaryNodes[j].index - numInteriorTexels;
			int fineBoundaryID = fineBoundaryIndex[gridChart.auxiliaryNodes[j].index];
			int nodeDegree = auxiliaryNodesDegree[auxiliaryID];
			Point2D< GeometryReal >nodePosition = gridChart.auxiliaryNodes[j].position;
			int corner[2] = { (int)floor(nodePosition[0] / gridChart.cellSizeW), (int)floor(nodePosition[1] / gridChart.cellSizeH) };
			unsigned int cellID = gridChart.cellIndices( corner[0] , corner[1] ).combined;
			if( cellID==-1 ) MK_THROW( "Invalid cell index. Node position " , nodePosition[0] / gridChart.cellSizeW , " " , nodePosition[1] / gridChart.cellSizeH );
			nodePosition[0] /= gridChart.cellSizeW;
			nodePosition[1] /= gridChart.cellSizeH;
			nodePosition[0] -= (GeometryReal)corner[0];
			nodePosition[1] -= (GeometryReal)corner[1];
			if( nodePosition[0] < 0-precision_error || nodePosition[0] > 1+precision_error || nodePosition[1] < 0-precision_error || nodePosition[1] > 1+precision_error )
				MK_THROW( "Sample out of unit box! (" , nodePosition[0] , " " , nodePosition[0] , ")" );
			for (int k = 0; k < 4; k++)
			{
				MatrixReal texelWeight = (MatrixReal)BilinearElementValue(k, nodePosition) / nodeDegree;
				if (fabs(texelWeight) > 1e-11) {
					auxiliaryNodesCumWeight[auxiliaryID] += texelWeight;
					int texelIndex = gridChart.bilinearElementIndices[cellID][k];
					if( nodeInfo[texelIndex].texelType==TexelType::InteriorSupported )
						MK_THROW( "Interior-supported texel cannot be in the support of an auxiliary node. Weight " , texelWeight , " (A)" );
					if( gridChart.auxiliaryNodes[j].index<numInteriorTexels || gridChart.auxiliaryNodes[j].index>numFineNodes || texelIndex<0 || texelIndex>numCoarseNodes ) MK_THROW( "Out of bounds index" );

					int coarseBoundaryID = boundaryAndDeepIndex[texelIndex] - 1;

					if( coarseBoundaryID<0 ) MK_THROW( "Coarse node is not boundary" );
					prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( fineBoundaryID , coarseBoundaryID , texelWeight ) );
				}
			}
		}
	}

	int numBoundaryTexels = gridAtlas.numBoundaryTexels;
	coarseBoundaryFineBoundaryProlongation = SetSparseMatrix( prolongationTriplets , numFineBoundaryNodes , numBoundaryTexels , false );
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