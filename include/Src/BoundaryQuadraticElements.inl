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
#ifdef NEW_INDEXING
		polygon.indices.push_back( AtlasBoundaryVertexIndex( triangle.indices[k] ) );
		polygon.vertexIndices.push_back( triangle.vertexIndices[k] );
#else // !NEW_INDEXING
		polygon.indices.push_back( triangle.indices[k] );
		polygon.vertexIndices.push_back( triangle.vertexIndices[k] );
#endif // NEW_INDEXING
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

#ifdef NEW_INDEXING
	polygon.outgoingEdgeIndices[0] = GridOrAtlasEdgeIndex::FromGridEdgeIndex( gridChart.edgeIndex( i+0 , j+0 , 0 ) );
	polygon.outgoingEdgeIndices[1] = GridOrAtlasEdgeIndex::FromGridEdgeIndex( gridChart.edgeIndex( i+1 , j+0 , 1 ) );
	polygon.outgoingEdgeIndices[2] = GridOrAtlasEdgeIndex::FromGridEdgeIndex( gridChart.edgeIndex( i+0 , j+1 , 0 ) );
	polygon.outgoingEdgeIndices[3] = GridOrAtlasEdgeIndex::FromGridEdgeIndex( gridChart.edgeIndex( i+0 , j+0 , 1 ) );
#else // !NEW_INDEXING
	polygon.outgoingEdgeIndices[0] = gridChart.edgeIndex( i+0 , j+0 , 0 );
	polygon.outgoingEdgeIndices[1] = gridChart.edgeIndex( i+1 , j+0 , 1 );
	polygon.outgoingEdgeIndices[2] = gridChart.edgeIndex( i+0 , j+1 , 0 );
	polygon.outgoingEdgeIndices[3] = gridChart.edgeIndex( i+0 , j+0 , 1 );
#endif // NEW_INDEXING

	return polygon;
}

// Computes the intersections of the edge with the grid lines and adds to the vector of IntersectionInfo
template< typename GeometryReal >
void AddEdgeGridIntersection
(
	Simplex< GeometryReal , 2 , 1 > edge ,
	const GridChart< GeometryReal > &gridChart ,
#ifdef NEW_INDEXING
	AtlasEdgeIndex edgeIndex ,
#else // !NEW_INDEXING
	int edgeIndex ,
#endif // NEW_INDEXING
	std::vector< IntersectionInfo< GeometryReal > > &intersections
)
{
	GeometryReal fmin[2] , fmax[2];
	for( unsigned int dir=0 ; dir<2 ; dir++ )
	{
		fmin[dir] = std::max< GeometryReal >( 0.f , std::min< GeometryReal >( edge[0][dir] , edge[1][dir] ) );
		fmax[dir] = std::min< GeometryReal >( 1.f , std::max< GeometryReal >( edge[0][dir] , edge[1][dir] ) );
	}

	int imin[2] , imax[2];
	imin[0] = static_cast< int >( floor( fmin[0] / gridChart.cellSizeW ) );
	imin[1] = static_cast< int >( floor( fmin[1] / gridChart.cellSizeH ) );
	imax[0] = static_cast< int >( floor( fmax[0] / gridChart.cellSizeW ) );
	imax[1] = static_cast< int >( floor( fmax[1] / gridChart.cellSizeH ) );

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
void InitializeChartBoundaryEdgeGridIntersections
(
	const AtlasChart< GeometryReal > & atlasChart ,
#ifdef NEW_INDEXING
	const std::map< AtlasVertexIndex , AtlasBoundaryVertexIndex > & atlasBoundaryVertexToIndex ,
#else // !NEW_INDEXING
	const std::map< unsigned int , unsigned int > & atlasBoundaryVertexToIndex ,
#endif // NEW_INDEXING
	GridChart< GeometryReal > & gridChart ,
	unsigned int & boundarySize ,
#ifdef NEW_INDEXING
	std::map< AtlasHalfEdgeIndex , std::vector< IntersectionInfo< GeometryReal > > > &atlasBoundaryHalfEdgeToIntersectionInfos ,
#else // !NEW_INDEXING
	std::map< unsigned int , std::vector< IntersectionInfo< GeometryReal > > > &atlasBoundaryHalfEdgeToIntersectionInfos ,
#endif // NEW_INDEXING
	std::map< SimplexIndex< 1 > , BoundarySegmentInfo< GeometryReal > > &segmentToBoundarySegmentInfo ,
	std::map< GridMeshIntersectionKey , NodeInfo< GeometryReal > > & gridMeshIntersectionKeyToNodeInfo
)
{
	for( unsigned int b=0 ; b<atlasChart.boundaryHalfEdges.size() ; b++ )
	{
#ifdef NEW_INDEXING
		ChartHalfEdgeIndex chartHalfEdge = atlasChart.boundaryHalfEdges[b];
		AtlasHalfEdgeIndex atlasHalfEdge = atlasChart.atlasHalfEdge( chartHalfEdge );
#else // !NEW_INDEXING
		unsigned int chartHalfEdge = atlasChart.boundaryHalfEdges[b];
		unsigned int atlasHalfEdge = atlasChart.atlasHalfEdge( chartHalfEdge );
#endif // NEW_INDEXING

		// The intersections of the boundary edge with the grid edges
		std::vector< IntersectionInfo< GeometryReal > > boundaryHalfEdgeIntersectionsInfo;

#ifdef NEW_INDEXING
		SimplexIndex< 1 , ChartVertexIndex > eIndex = atlasChart.edgeIndex( atlasChart.boundaryHalfEdges[b] );
		Simplex< GeometryReal , 2 , 1 > edge( atlasChart.vertex( ChartVertexIndex( eIndex[0] ) ) - gridChart.corner , atlasChart.vertex( ChartVertexIndex( eIndex[1] ) ) - gridChart.corner );
#else // !NEW_INDEXING
		SimplexIndex< 1 > eIndex = atlasChart.edgeIndex( atlasChart.boundaryHalfEdges[b] );
		Simplex< GeometryReal , 2 , 1 > edge( atlasChart.vertices[ eIndex[0] ] - gridChart.corner , atlasChart.vertices[ eIndex[1] ] - gridChart.corner );
#endif // NEW_INDEXING

		// Add the end-points of the edge
		boundaryHalfEdgeIntersectionsInfo.emplace_back( GridMeshIntersectionKey::ChartVertexKey( eIndex[0] ) , edge[0] , 0 );
		boundaryHalfEdgeIntersectionsInfo.emplace_back( GridMeshIntersectionKey::ChartVertexKey( eIndex[1] ) , edge[1] , 1 );

		// Add the points where the edge intersects the grid
		AddEdgeGridIntersection( edge , gridChart , atlasChart.atlasEdge( chartHalfEdge ) , boundaryHalfEdgeIntersectionsInfo );

		// Sort intersections
		std::sort( boundaryHalfEdgeIntersectionsInfo.begin() , boundaryHalfEdgeIntersectionsInfo.end() , IntersectionInfo< GeometryReal >::CompareByTime );

		// Assign indices to the boundary-edge/grid intersections
		for( unsigned int i=0 ; i<boundaryHalfEdgeIntersectionsInfo.size() ; i++ )
		{
			// If we haven't seen this vertex before...
			if( gridMeshIntersectionKeyToNodeInfo.find( boundaryHalfEdgeIntersectionsInfo[i].intersectionKey )==gridMeshIntersectionKeyToNodeInfo.end() )
			{
#ifdef NEW_INDEXING
				AtlasBoundaryVertexIndex index;
#else // !NEW_INDEXING
				unsigned int index;
#endif // NEW_INDEXING
#ifdef NEW_INDEXING
				if( std::optional< ChartVertexIndex > v = boundaryHalfEdgeIntersectionsInfo[i].intersectionKey.chartVertex() )	// Start/end vertex
#else // !NEW_INDEXING
				if( std::optional< unsigned int > v = boundaryHalfEdgeIntersectionsInfo[i].intersectionKey.meshVertex() )	// Start/end vertex
#endif // NEW_INDEXING
				{
					if( i!=0 && i!=boundaryHalfEdgeIntersectionsInfo.size()-1 ) MK_ERROR_OUT( "Expected boundary vertex" );
					auto iter = atlasBoundaryVertexToIndex.find( atlasChart.atlasVertex( *v ) );
					if( iter==atlasBoundaryVertexToIndex.end() ) MK_THROW( "Boundary vertex not found" );
					index = iter->second;
				}
				else
				{
					if( i==0 && i==boundaryHalfEdgeIntersectionsInfo.size()-1 ) MK_ERROR_OUT( "Expected interior vertex" );
#ifdef NEW_INDEXING
					index = AtlasBoundaryVertexIndex( boundarySize++ );
#else // !NEW_INDEXING
					index = boundarySize++;
#endif // NEW_INDEXING
				}

				gridMeshIntersectionKeyToNodeInfo[ boundaryHalfEdgeIntersectionsInfo[i].intersectionKey ]
					= NodeInfo< GeometryReal >( index , boundaryHalfEdgeIntersectionsInfo[i].position );

				gridChart.auxiliaryNodes.emplace_back( boundaryHalfEdgeIntersectionsInfo[i].position , index );
			}
			boundaryHalfEdgeIntersectionsInfo[i].index = gridMeshIntersectionKeyToNodeInfo[ boundaryHalfEdgeIntersectionsInfo[i].intersectionKey ].index;
		}

		// Add the segment information
		for( unsigned int i=0 ; i<boundaryHalfEdgeIntersectionsInfo.size()-1 ; i++ )
		{
			SimplexIndex< 1 > eIndex( gridMeshIntersectionKeyToNodeInfo[ boundaryHalfEdgeIntersectionsInfo[i].intersectionKey ].index , gridMeshIntersectionKeyToNodeInfo[ boundaryHalfEdgeIntersectionsInfo[i+1].intersectionKey ].index );
			BoundarySegmentInfo< GeometryReal > segmentInfo( boundaryHalfEdgeIntersectionsInfo[i].time , boundaryHalfEdgeIntersectionsInfo[i+1].time , chartHalfEdge );

			if( segmentToBoundarySegmentInfo.find( eIndex )!=segmentToBoundarySegmentInfo.end() ) MK_THROW( "Replicated segment key" );
			segmentToBoundarySegmentInfo[ eIndex ] = segmentInfo;
		}
		atlasBoundaryHalfEdgeToIntersectionInfos[ atlasHalfEdge ] = boundaryHalfEdgeIntersectionsInfo;
	}
}

template< typename GeometryReal >
void InitializeChartBoundaryPolygons
(
#ifdef NEW_INDEXING
	const std::vector< AtlasHalfEdgeIndex > & oppositeHalfEdge ,
#else // !NEW_INDEXING
	const std::vector< unsigned int > &oppositeHalfEdge ,
#endif // NEW_INDEXING
	const AtlasChart< GeometryReal > &atlasChart ,
	GridChart< GeometryReal > &gridChart ,
	unsigned int numInteriorNodes ,
	unsigned int numBoundaryVertices ,
	unsigned int numBoundaryNodes ,
#ifdef NEW_INDEXING
	const std::map< AtlasHalfEdgeIndex , std::vector< IntersectionInfo< GeometryReal > > > & atlasBoundaryHalfEdgeToIntersectionInfos ,
#else // !NEW_INDEXING
	const std::map< unsigned int , std::vector< IntersectionInfo< GeometryReal > > > & atlasBoundaryHalfEdgeToIntersectionInfos ,
#endif // NEW_INDEXING
	const std::map< SimplexIndex< 1 > , BoundarySegmentInfo< GeometryReal > > & segmentToBoundarySegmentInfo ,
	const std::map< GridMeshIntersectionKey , NodeInfo< GeometryReal > > & gridMeshIntersectionKeyToNodeInfo ,
	unsigned int chartID
)
{
#ifdef SEPARATE_POLYGONS
	// A mapping giving the clipped polygons, keyed off of cell index
	std::map< unsigned int , std::vector< std::vector< unsigned long long > > > cellPolygons;
#else // !SEPARATE_POLYGONS
	// An association of segments to individual cells of a chart
	std::map< unsigned int , std::pair< Point< unsigned int , 2 > , std::vector< std::pair< GridMeshIntersectionKey , GridMeshIntersectionKey > > > > cellSegments;
#endif // SEPARATE_POLYGONS

#ifdef NEW_INDEXING
	auto GetIndexedTriangle = [&]( ChartTriangleIndex t )
		{
			IndexedIntersectionTriangle< GeometryReal > indexedTriangle;
			SimplexIndex< 2 , ChartVertexIndex > tri = atlasChart.triangleIndex( t );

			for( unsigned int k=0 ; k<3 ; k++ )
			{
				indexedTriangle.vertices[k] = atlasChart.vertex( tri[k] ) - gridChart.corner;
				indexedTriangle.cornerKeys[k] = GridMeshIntersectionKey::ChartVertexKey( tri[k] );
				indexedTriangle.outgoingEdgeIndices[k] = GridOrAtlasEdgeIndex::FromAtlasEdgeIndex( atlasChart.atlasEdge( GetChartHalfEdgeIndex( t , k ) ) );
			}
			return indexedTriangle;
		};
#else // !NEW_INDEXING
	auto GetIndexedTriangle = [&]( unsigned int t )
		{
			IndexedIntersectionTriangle< GeometryReal > indexedTriangle;
			Point2D< GeometryReal > tPos[3];
			for( unsigned int i=0 ; i<3 ; i++ ) tPos[i] = atlasChart.vertices[ atlasChart.triangles[t][i] ] - gridChart.corner;

			for( unsigned int k=0 ; k<3 ; k++ )
			{
				indexedTriangle.vertices[k] = tPos[k];
				indexedTriangle.cornerKeys[k] = GridMeshIntersectionKey::ChartVertexKey( atlasChart.triangles[t][k] );
				indexedTriangle.outgoingEdgeIndices[k] = atlasChart.atlasEdge( 3*t+k );
			}
			return indexedTriangle;
		};
#endif // NEW_INDEXING

	auto GetIndexedIntersectionPolygon = [&]( unsigned int i , unsigned int j )
		{
			return IndexedPolygonFromCell( i , j , gridChart );
		};

#ifdef USE_RASTERIZER
	using Range = RegularGrid< 2 >::Range;
	using Index = RegularGrid< 2 >::Index;
	Range nodeRange , cellRange;
	nodeRange.second[0] = gridChart.width  , cellRange.second[0] = gridChart.width-1;
	nodeRange.second[1] = gridChart.height , cellRange.second[1] = gridChart.height-1;

#ifdef NEW_INDEXING
	auto GetSimplex = [&]( ChartTriangleIndex t )
		{
			Simplex< double , 2 , 2 > simplex;
			SimplexIndex< 2 , ChartVertexIndex > tri = atlasChart.triangleIndex( t );
			for( unsigned int k=0 ; k<=2 ; k++ )
			{
				simplex[k] = atlasChart.vertex( tri[k] ) - gridChart.corner;
				simplex[k][0] /= gridChart.cellSizeW;
				simplex[k][1] /= gridChart.cellSizeH;
			}
			return simplex;
		};
#else // !NEW_INDEXING
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
#endif // NEW_INDEXING
#endif // USE_RASTERIZER

	//(1) Rasterize all triangles and, if they fall into the support of a boundary cell, add them to the list of triangles associated with that cell
#ifdef NEW_INDEXING
	for( unsigned int t=0 ; t<atlasChart.numTriangles() ; t++ )
#else // !NEW_INDEXING
	for( unsigned int t=0 ; t<atlasChart.triangles.size() ; t++ )
#endif // NEW_INDEXING
	{
#ifdef NEW_INDEXING
		IndexedIntersectionTriangle< GeometryReal > indexedTriangle = GetIndexedTriangle( ChartTriangleIndex(t) );
#else // !NEW_INDEXING
		IndexedIntersectionTriangle< GeometryReal > indexedTriangle = GetIndexedTriangle( t );
#endif // NEW_INDEXING
#ifdef USE_RASTERIZER
		// Compute the associated triangle in (shifted) texel coordinates
#ifdef NEW_INDEXING
		Simplex< double , 2 , 2 > simplex = GetSimplex( ChartTriangleIndex(t) );
#else // !NEW_INDEXING
		Simplex< double , 2 , 2 > simplex = GetSimplex( t );
#endif // NEW_INDEXING

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
					cellSegments[cellID].first = Point< unsigned int , 2 >( I[0]+gridChart.cornerCoords[0] , I[1]+gridChart.cornerCoords[1] );

					IndexedIntersectionPolygon< GeometryReal > poly = GetIndexedIntersectionPolygon( I[0] , I[1] );

					if( ClipIndexedIntersectionPolygonToIndexedIntersectionTriangle( poly , indexedTriangle ) )
						for( int s=0 ; s<poly.cornerKeys.size() ; s++ )
							cellSegments[cellID].second.push_back( std::make_pair( poly.cornerKeys[s] , poly.cornerKeys[ (s+1) % poly.cornerKeys.size() ] ) );
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
		unsigned int cellID = cellIter->first;
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
#ifdef NEW_INDEXING
					for( unsigned int t=0 ; t<atlasChart.numTriangles() ; t++ )
#else // !NEW_INDEXING
					for( unsigned int t=0 ; t<atlasChart.triangles.size() ; t++ )
#endif // NEW_INDEXING
					{
#ifdef NEW_INDEXING
						IndexedIntersectionTriangle< GeometryReal > indexedTriangle = GetIndexedTriangle( ChartTriangleIndex(t) );
#else // !NEW_INDEXING
						IndexedIntersectionTriangle< GeometryReal > indexedTriangle = GetIndexedTriangle( t );
#endif // NEW_INDEXING
						IndexedIntersectionPolygon< GeometryReal > cellPolygon = GetIndexedIntersectionPolygon( _idx[0] , _idx[1] );
#ifdef NEW_INDEXING
						if( ClipIndexedIntersectionPolygonToIndexedIntersectionTriangle( cellPolygon , indexedTriangle ) ) sStream << " " << atlasChart.atlasTriangle( ChartTriangleIndex(t) );
#else // !NEW_INDEXING
						if( ClipIndexedIntersectionPolygonToIndexedIntersectionTriangle( cellPolygon , indexedTriangle ) ) sStream << " " << atlasChart.atlasTriangle(t);
#endif // NEW_INDEXING
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
#ifdef NEW_INDEXING
		std::vector< std::vector< ChartVertexIndex > > loopAtlasVertexIndices( loopKeys.size() );
#else // !NEW_INDEXING
		std::vector< std::vector< unsigned int > > loopAtlasVertexIndices( loopKeys.size() );
#endif // NEW_INDEXING

		// Remove vertices generated when non-boundary edges intersect the grid
		for( unsigned int i=0 ; i<loopKeys.size() ; i++ )
		{
			std::vector< GridMeshIntersectionKey > reducedLoop;		// The subset of vertices that are either grid nodes or boundary nodes
			std::vector< NodeInfo< GeometryReal > > currentLoopNodes;
#ifdef NEW_INDEXING
			std::vector< ChartVertexIndex > currentAtlasVertexIndices;
#else // !NEW_INDEXING
			std::vector< unsigned int > currentAtlasVertexIndices;
#endif // NEW_INDEXING

			// For each polygon vertex
			for( unsigned int j=0 ; j<loopKeys[i].size() ; j++ )
			{
				GridMeshIntersectionKey currentVertexKey = loopKeys[i][j];

#ifdef NEW_INDEXING
				if( std::optional< GridNodeIndex > g = currentVertexKey.gridNode() )	// Texel node
#else // !NEW_INDEXING
				if( std::optional< unsigned int > g = currentVertexKey.gridNode() )	// Texel node
#endif // NEW_INDEXING
				{
					unsigned int pi , pj;
					if( !gridChart.factorNodeIndex( *g , pi , pj ) ) MK_ERROR_OUT( "Could not factor node index" );

					// Confirm that the texel is inside the chart
					unsigned int interiorTexelIndex = gridChart.texelIndices(pi,pj).interiorOrCovered;
					if( interiorTexelIndex==-1 ) MK_THROW( "Invalid texel: " , Point2D< int >( pi , pj ) , " -> " , gridChart.texelIndices(pi,pj).interiorOrCovered , " : " , gridChart.texelIndices(pi, pj).interior );

					reducedLoop.push_back( currentVertexKey );
					currentLoopNodes.emplace_back( interiorTexelIndex , gridChart.nodePosition(pi,pj) );
#ifdef NEW_INDEXING
					currentAtlasVertexIndices.push_back( ChartVertexIndex(-1) );
#else // !NEW_INDEXING
					currentAtlasVertexIndices.push_back(-1);
#endif // NEW_INDEXING
				}
				else
				{
					auto iter = gridMeshIntersectionKeyToNodeInfo.find( currentVertexKey );
					if( iter!=gridMeshIntersectionKeyToNodeInfo.end() )
					{
						reducedLoop.push_back( currentVertexKey );
						currentLoopNodes.push_back( iter->second );
						currentLoopNodes.back().index += numInteriorNodes;

#ifdef NEW_INDEXING
						if( std::optional< ChartVertexIndex > v = currentVertexKey.chartVertex() ) currentAtlasVertexIndices.push_back( *v );
						else                                                                       currentAtlasVertexIndices.push_back( ChartVertexIndex(-1) );
#else // !NEW_INDEXING
#ifdef NEW_INDEXING
						if( std::optional< ChartVertexIndex > v = currentVertexKey.chartVertex() ) currentAtlasVertexIndices.push_back( *v );
#else // !NEW_INDEXING
						if( std::optional< unsigned int > v = currentVertexKey.meshVertex() ) currentAtlasVertexIndices.push_back( *v );
#endif // NEW_INDEXING
						else                                                                       currentAtlasVertexIndices.push_back(-1);
#endif // NEW_INDEXING
					}
				}
			}
			if( reducedLoop.size()==0 ) MK_THROW( "Reduced loop cannot be empty. Original loop size " , loopKeys[i].size() );
			loopKeys[i] = reducedLoop;
			loopNodes[i] = currentLoopNodes;
			loopAtlasVertexIndices[i] = currentAtlasVertexIndices;
		}

		// Add the intermediate vertices for boundary segments
#ifdef NEW_INDEXING
		std::vector< std::vector< AtlasEdgeIndex > > loopAtlasEdges( loopKeys.size() );
		std::vector< std::vector< AtlasEdgeIndex > > loopAtlasVertexParentEdges( loopKeys.size() );
#else // !NEW_INDEXING
		std::vector< std::vector< unsigned int > > loopAtlasEdges( loopKeys.size() );
		std::vector< std::vector< unsigned int > > loopAtlasVertexParentEdges( loopKeys.size() );
#endif // NEW_INDEXING
		for( unsigned int i=0 ; i<loopKeys.size() ; i++ )
		{
#ifdef NEW_INDEXING
			std::vector< std::vector< AtlasBoundaryVertexIndex > > indicesToInsert;
#else // !NEW_INDEXING
			std::vector< std::vector< unsigned int > > indicesToInsert;
#endif // NEW_INDEXING
			std::vector< std::vector< GeometryReal > > timesToInsert;
			std::vector< bool > isInsertionPos( loopKeys[i].size() , false );
#ifdef NEW_INDEXING
			std::vector< AtlasEdgeIndex > polygonAtlasEdgeIndex( loopKeys[i].size() , AtlasEdgeIndex(-1) );
			std::vector< AtlasEdgeIndex > polygonAtlasVertexParentEdges( loopKeys[i].size() , AtlasEdgeIndex(-1) );	// The mesh edge the boundary vertex/node comes from
#else // !NEW_INDEXING
			std::vector< unsigned int > polygonAtlasEdgeIndex( loopKeys[i].size() , -1 );
			std::vector< unsigned int > polygonAtlasVertexParentEdges( loopKeys[i].size() , -1 );	// The mesh edge the boundary vertex/node comes from
#endif // NEW_INDEXING

			for( unsigned int j=0 ; j<loopKeys[i].size() ; j++ )
			{
				GridMeshIntersectionKey currentVertexKey = loopKeys[i][j];
				GridMeshIntersectionKey    nextVertexKey = loopKeys[i][(j+1) % loopKeys[i].size()];

#ifdef NEW_INDEXING
				if( std::optional< std::pair< GridEdgeIndex , AtlasEdgeIndex > > x = currentVertexKey.intersection() ) polygonAtlasVertexParentEdges[j] = x->second;
#else // !NEW_INDEXING
				if( std::optional< std::pair< unsigned int , unsigned int > > x = currentVertexKey.intersection() ) polygonAtlasVertexParentEdges[j] = x->second;
#endif // NEW_INDEXING

				// If this is part of a boundary edge
				auto currentBoundaryNodeIter = gridMeshIntersectionKeyToNodeInfo.find( currentVertexKey );
				auto    nextBoundaryNodeIter = gridMeshIntersectionKeyToNodeInfo.find(    nextVertexKey );
				if( currentBoundaryNodeIter!=gridMeshIntersectionKeyToNodeInfo.end() && nextBoundaryNodeIter!=gridMeshIntersectionKeyToNodeInfo.end() )
				{
					SimplexIndex< 1 > eIndex( currentBoundaryNodeIter->second.index , nextBoundaryNodeIter->second.index );
					auto segmentToBoundarySegmentInfoIter = segmentToBoundarySegmentInfo.find( eIndex );
					if( segmentToBoundarySegmentInfoIter!=segmentToBoundarySegmentInfo.end() )
					{
						BoundarySegmentInfo< GeometryReal > segmentInfo = segmentToBoundarySegmentInfoIter->second;
						polygonAtlasEdgeIndex[j] = atlasChart.atlasEdge( segmentInfo.chartHalfEdge );

#ifdef NEW_INDEXING
						AtlasHalfEdgeIndex oppHalfEdge = oppositeHalfEdge[ static_cast< unsigned int >( atlasChart.atlasHalfEdge( segmentInfo.chartHalfEdge ) ) ];
						if( oppHalfEdge!=AtlasHalfEdgeIndex(-1) )
#else // !NEW_INDEXING
						unsigned int oppHalfEdge = oppositeHalfEdge[ atlasChart.atlasHalfEdge( segmentInfo.chartHalfEdge ) ];
						if( oppHalfEdge!=-1 )
#endif // NEW_INDEXING
						{
							GeometryReal startTime = segmentInfo.startTime;
							GeometryReal endTime = segmentInfo.endTime;

							// Find the information for the opposite half-edge
							auto atlasBoundaryHalfEdgeToIntersectionInfosIter = atlasBoundaryHalfEdgeToIntersectionInfos.find( oppHalfEdge );
							if( atlasBoundaryHalfEdgeToIntersectionInfosIter==atlasBoundaryHalfEdgeToIntersectionInfos.end() )
								MK_THROW( "Opposite edge intersections not found. Current  edge " , atlasChart.atlasHalfEdge( segmentInfo.chartHalfEdge ) , ". Opposite " , oppHalfEdge );

#ifdef NEW_INDEXING
							std::vector< AtlasBoundaryVertexIndex > segmentIndicesToInsert;
#else // !NEW_INDEXING
							std::vector< unsigned int > segmentIndicesToInsert;
#endif // NEW_INDEXING
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
			std::vector< NodeInfo< GeometryReal > > expandedLoopNodes;
#ifdef NEW_INDEXING
			std::vector< ChartVertexIndex > expandedLoopAtlasVertexIndices;
			std::vector< AtlasEdgeIndex > expandedLoopAtlasEdgeIndex;
			std::vector< AtlasEdgeIndex > expandedVertexParentEdgeIndex;
#else // !NEW_INDEXING
			std::vector< unsigned int > expandedLoopAtlasVertexIndices;
			std::vector< unsigned int > expandedLoopAtlasEdgeIndex;
			std::vector< unsigned int > expandedVertexParentEdgeIndex;
#endif // NEW_INDEXING
			for( unsigned int j=0 ; j<loopKeys[i].size() ; j++ )
			{
				// Add information for original loop vertex

				//expandedLoopKeys.push_back( loopKeys[i][j] );
				expandedLoopNodes.push_back( loopNodes[i][j] );
				expandedLoopAtlasVertexIndices.push_back( loopAtlasVertexIndices[i][j] );

#ifdef NEW_INDEXING
				AtlasEdgeIndex currentSegmentAtlasEdgeIndex = polygonAtlasEdgeIndex[j];
				expandedLoopAtlasEdgeIndex.push_back( currentSegmentAtlasEdgeIndex );
				expandedVertexParentEdgeIndex.push_back( polygonAtlasVertexParentEdges[j] );
#else // !NEW_INDEXING
				unsigned int currentSegmentAtlasEdgeIndex = polygonAtlasEdgeIndex[j];
				expandedLoopAtlasEdgeIndex.push_back( currentSegmentAtlasEdgeIndex );
				expandedVertexParentEdgeIndex.push_back( polygonAtlasVertexParentEdges[j] );
#endif // NEW_INDEXING

				// If there are vertices to insert from the opposite boundary, do that now
				if( isInsertionPos[j] )
				{
#ifdef NEW_INDEXING
					const std::vector< AtlasBoundaryVertexIndex > & _indices = indicesToInsert[insertionCount];
#else // !NEW_INDEXING
					const std::vector< unsigned int > & _indices = indicesToInsert[insertionCount];
#endif // NEW_INDEXING
					const std::vector< GeometryReal > & _times = timesToInsert[insertionCount];
					Point2D< GeometryReal > currentPos = loopNodes[i][j].position;
					Point2D< GeometryReal >    nextPos = loopNodes[i][ (j+1)%loopKeys[i].size() ].position;

#ifdef NEW_INDEXING
					if( currentSegmentAtlasEdgeIndex==AtlasEdgeIndex(-1) ) MK_THROW( "Invalid atlas edge index" );
#else // !NEW_INDEXING
					if( currentSegmentAtlasEdgeIndex==-1 ) MK_THROW( "Invalid atlas edge index" );
#endif // NEW_INDEXING

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
#ifdef NEW_INDEXING
						if( static_cast< unsigned int >(_indices[k])<numBoundaryVertices || static_cast< unsigned int >(_indices[k])>numBoundaryNodes ) MK_THROW( "Out of bounds index: " , _indices[k] , " not in " , numBoundaryVertices , " " , numBoundaryNodes );
#else // !NEW_INDEXING
						if( _indices[k]<numBoundaryVertices || _indices[k]>numBoundaryNodes ) MK_THROW( "Out of bounds index: " , _indices[k] , " not in " , numBoundaryVertices , " " , numBoundaryNodes );
#endif // NEW_INDEXING

						gridChart.auxiliaryNodes.emplace_back( interpolatedPos , _indices[k] + numInteriorNodes );

						expandedLoopNodes.emplace_back( _indices[k] + numInteriorNodes , interpolatedPos );
#ifdef NEW_INDEXING
						expandedLoopAtlasVertexIndices.push_back( ChartVertexIndex(-1) );
#else // !NEW_INDEXING
						expandedLoopAtlasVertexIndices.push_back(-1);
#endif // NEW_INDEXING
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
#ifdef NEW_INDEXING
			for( unsigned int j=0 ; j<loopNodes[i].size() ; j++ ) poly.indices[j] = loopNodes[i][j].index , poly.vertices[j] = loopNodes[i][j].position;
#else // !NEW_INDEXING
			for( unsigned int j=0 ; j<loopNodes[i].size() ; j++ ) poly.indices[j] = loopNodes[i][j].index , poly.vertices[j] = loopNodes[i][j].position;
#endif // NEW_INDEXING
			poly.vertexIndices = loopAtlasVertexIndices[i];
			poly.atlasEdgeIndices = loopAtlasEdges[i];
			poly.atlasVertexParentEdge = loopAtlasVertexParentEdges[i];
			gridChart.boundaryPolygons[ cellID ].push_back( poly );
		}
	}
}


template< typename GeometryReal >
unsigned int InitializeBoundaryPolygons
(
	const std::vector< AtlasChart< GeometryReal > > &atlasCharts ,
	const typename AtlasChart< GeometryReal >::AtlasInfo &atlasInfo ,
	std::vector< GridChart< GeometryReal > > &gridCharts ,
	unsigned int numInteriorNodes
)
{ //Fine System

	// Get the initial count of atlas vertices residing on the boundary
	unsigned int boundarySize = (unsigned int)atlasInfo.atlasBoundaryVertexToIndex.size();

	// A map taking (atlas) boundary half-edge indices to their decomposition by grid edges
#ifdef NEW_INDEXING
	std::map< AtlasHalfEdgeIndex , std::vector< IntersectionInfo< GeometryReal > > > atlasBoundaryHalfEdgeToIntersectionInfos;
	#else // !NEW_INDEXING
	std::map< unsigned int , std::vector< IntersectionInfo< GeometryReal > > > atlasBoundaryHalfEdgeToIntersectionInfos;
#endif // NEW_INDEXING

	// A map taking a segment to the associated segment information
	std::vector< std::map< SimplexIndex< 1 > , BoundarySegmentInfo< GeometryReal > > > segmentToBoundarySegmentInfo( gridCharts.size() );

	// A map taking an intersection key to the associated node information
	std::vector< std::map< GridMeshIntersectionKey , NodeInfo< GeometryReal > > > gridMeshIntersectionKeyToNodeInfo( gridCharts.size() );

	for( unsigned int i=0 ; i<gridCharts.size() ; i++ )
	{
		// -- Split the chart's boundary half-edges to the edges of the grid
		// -- Add associated auxiliary nodes
		// -- Add associated segments
		InitializeChartBoundaryEdgeGridIntersections( atlasCharts[i] , atlasInfo.atlasBoundaryVertexToIndex , gridCharts[i] , boundarySize , atlasBoundaryHalfEdgeToIntersectionInfos , segmentToBoundarySegmentInfo[i] , gridMeshIntersectionKeyToNodeInfo[i] );	
	}

	// Offset the auxiliary nodes' indices
	for( unsigned int i=0 ; i<gridCharts.size() ; i++ ) for( unsigned int j=0 ; j<gridCharts[i].auxiliaryNodes.size() ; j++ ) gridCharts[i].auxiliaryNodes[j].index += numInteriorNodes;

	for( unsigned int i=0 ; i<gridCharts.size() ; i++ )
	{
		try
		{
			InitializeChartBoundaryPolygons( atlasInfo.oppositeHalfEdges , atlasCharts[i] , gridCharts[i] , numInteriorNodes , (unsigned int)atlasInfo.atlasBoundaryVertexToIndex.size() , boundarySize , atlasBoundaryHalfEdgeToIntersectionInfos , segmentToBoundarySegmentInfo[i] , gridMeshIntersectionKeyToNodeInfo[i] , i );
		}
		catch( const Exception & e )
		{
			std::cerr << e.what() << std::endl;
			MK_THROW( "While processing chart: " , i );
		}
	}

	return boundarySize;
}

template< typename GeometryReal >
void InitializeChartQuadraticElements
(
	GridChart< GeometryReal > &gridChart ,
	std::map< SimplexIndex< 1 > , unsigned int > &midPointIndex ,
	unsigned int &lastMidPointIndex ,
	unsigned int previouslyAddedNodes
)
{
	const std::vector< std::vector< IndexedPolygon< GeometryReal > > > & boundaryPolygons = gridChart.boundaryPolygons;
	std::vector< std::vector< BoundaryIndexedTriangle< GeometryReal > > > & boundaryTriangles = gridChart.boundaryTriangles;
	boundaryTriangles.resize( gridChart.numBoundaryCells() );
	int numBoundaryTriangles = 0;
	for( int i=0 ; i<boundaryPolygons.size() ; i++ )
	{
		for( int j=0 ; j<boundaryPolygons[i].size() ; j++ )
		{
			const IndexedPolygon< GeometryReal > &currentPolygon = boundaryPolygons[i][j];
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
				for( unsigned int v=0 ; v<3 ; v++ )
				{
					int currentPolygonVertex = delanauyTriangles[k][v];
					triangle.indices[v] = currentPolygon.indices[ currentPolygonVertex ];
					triangle.vertices[v] = currentPolygon.vertices[ currentPolygonVertex ];
					triangle.vertexIndices[v] = currentPolygon.vertexIndices[ currentPolygonVertex ];
					triangle.atlasVertexParentEdge[v] = currentPolygon.atlasVertexParentEdge[ currentPolygonVertex ];
#ifdef NEW_INDEXING
					triangle.atlasEdgeIndices[v] = AtlasEdgeIndex(-1);
#else // !NEW_INDEXING
					triangle.atlasEdgeIndices[v] = -1;
#endif // NEW_INDEXING
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

				for( unsigned int v=0 ; v<3 ; v++ )
				{
					Point2D< GeometryReal > edgeMidPoint = ( triangle.vertices[ (v+1)%3 ] + triangle.vertices[ (v+2)%3 ] ) / 2;
					int edgeCorners[2] = { (int)triangle.indices[(v + 1) % 3], (int)triangle.indices[(v + 2) % 3] };
					if( edgeCorners[0] > edgeCorners[1]) std::swap( edgeCorners[0] , edgeCorners[1] );
					SimplexIndex< 1 > edgeKey( edgeCorners[0] , edgeCorners[1] );

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
	std::map< SimplexIndex< 1 > , unsigned int > midPointIndex;
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
	gridAtlas.numBoundaryNodes = InitializeBoundaryPolygons( atlasCharts , atlasInfo , gridAtlas.gridCharts , gridAtlas.numInteriorTexels );
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
	const typename GridAtlas<>::IndexConverter & indexConverter = gridAtlas.indexConverter;

	fineBoundaryIndex.resize( gridAtlas.numFineNodes , static_cast< unsigned int >(-1) );
	unsigned int lastFineBoundaryIndex = 0;
	for( unsigned int i=0 ; i<gridAtlas.gridCharts.size() ; i++ )
	{
		const GridChart< GeometryReal > &gridChart = gridAtlas.gridCharts[i];
		for( int j=0 ; j<gridChart.texelIndices.size() ; j++ )
		{
			if( gridChart.texelIndices[j].interiorOrCovered!=-1 && gridChart.texelIndices[j].interior==-1 )
			{ //Interior but not deep
				int coarseGlobalIndex = gridChart.texelIndices[j].combined;

				unsigned int boundaryIndex = indexConverter.supportedToBoundary( coarseGlobalIndex );
				if( boundaryIndex!=-1 ) prolongationTriplets.emplace_back( lastFineBoundaryIndex , boundaryIndex , (MatrixReal)1. );
				else MK_THROW( "Coarse node is not boundary. Global index " , coarseGlobalIndex , ". Boundary index " , boundaryIndex );
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

	for( unsigned int i=0 ; i<gridCharts.size() ; i++ )
	{
		const GridChart< GeometryReal > &gridChart = gridCharts[i];
		for( unsigned int j=0 ; j<gridChart.auxiliaryNodes.size() ; j++ )
		{
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
					int texelIndex = gridChart.combinedCellCombinedBilinearElementIndices[cellID][k];
					if( nodeInfo[texelIndex].texelType==TexelType::InteriorSupported )
						MK_THROW( "Interior-supported texel cannot be in the support of an auxiliary node. Weight " , texelWeight , " (A)" );
					if( gridChart.auxiliaryNodes[j].index<numInteriorTexels || gridChart.auxiliaryNodes[j].index>numFineNodes || texelIndex<0 || texelIndex>numCoarseNodes ) MK_THROW( "Out of bounds index" );

					unsigned int boundaryIndex = indexConverter.supportedToBoundary( texelIndex );
					if( boundaryIndex==-1 ) MK_THROW( "Coarse node is not boundary" );
					prolongationTriplets.emplace_back( fineBoundaryID , boundaryIndex , texelWeight );
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