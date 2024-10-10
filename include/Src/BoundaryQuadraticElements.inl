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

#include <Misha/Miscellany.h>
#include "PolygonClipping.h"
#include "MapLoop.h"


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
void SetIndexedPolygonFromCell( int i , int j , int width , int height , IndexedIntersectionPolygon< GeometryReal > &polygon )
{
	polygon.vertices.resize(4);
	polygon.indices.resize(4);
	polygon.edgeIndices.resize(4);

	polygon.vertices[0] = Point2D< GeometryReal >( (GeometryReal)( i+0.5 ) / width , (GeometryReal)( j+0.5 ) / height );
	polygon.vertices[1] = Point2D< GeometryReal >( (GeometryReal)( i+1.5 ) / width , (GeometryReal)( j+0.5 ) / height );
	polygon.vertices[2] = Point2D< GeometryReal >( (GeometryReal)( i+1.5 ) / width , (GeometryReal)( j+1.5 ) / height );
	polygon.vertices[3] = Point2D< GeometryReal >( (GeometryReal)( i+0.5 ) / width , (GeometryReal)( j+1.5 ) / height );

	polygon.indices[0] = SetIntersectionKey(4294967295, i + j*width);
	polygon.indices[1] = SetIntersectionKey(4294967295, i + 1 + j*width);
	polygon.indices[2] = SetIntersectionKey(4294967295, i + 1 + (j + 1)*width);
	polygon.indices[3] = SetIntersectionKey(4294967295, i + (j + 1)*width);

	polygon.edgeIndices[0] = width*height + i + j* width;
	polygon.edgeIndices[1] = i + 1 + j* width;
	polygon.edgeIndices[2] = width*height + i + (j + 1)* width;
	polygon.edgeIndices[3] = i + j* width;
}

template< typename GeometryReal >
void ComputeAtlasEdgeGridIntersection( Point2D< GeometryReal > p[2] , const GridChart< GeometryReal > &gridChart , int atlasEdgeIndex , std::vector< IntersectionInfo< GeometryReal > > &intersections )
{
	GeometryReal fmin[2] , fmax[2];
	for( int c=0 ; c<2 ; c++ )
	{
		fmin[c] = std::min< GeometryReal >( p[0][c] , p[1][c] );
		fmin[c] = std::max< GeometryReal >( fmin[c] , 0.f );
		fmax[c] = std::max< GeometryReal >( p[0][c] , p[1][c] );
		fmax[c] = std::min< GeometryReal >( fmax[c] , 1.f );
	}
	int imin[2];
	int imax[2];

	imin[0] = static_cast< int >( floor( fmin[0] / gridChart.cellSizeW ) );
	imin[1] = static_cast< int >( floor( fmin[1] / gridChart.cellSizeH ) );
	imax[0] = static_cast< int >( floor( fmax[0] / gridChart.cellSizeW ) );
	imax[1] = static_cast< int >( floor( fmax[1] / gridChart.cellSizeH ) );

	for (int c = 0; c < 2; c++)
	{
		for (int s = imin[c]; s <= imax[c]; s++)
		{
			GeometryReal cellSize = c == 0 ? gridChart.cellSizeW : gridChart.cellSizeH;
			GeometryReal oppCellSize = c == 0 ? gridChart.cellSizeH : gridChart.cellSizeW;
			GeometryReal level = (GeometryReal)s*cellSize;
			GeometryReal alpha = ( level-p[0][c] ) / ( p[1][c]-p[0][c] );
			if( alpha>0 && alpha<1 )
			{
				Point2D< GeometryReal > position = p[0] * (GeometryReal)( 1.-alpha ) + p[1] * alpha;
				int oppDimIndex = (int)floor(position[1 - c] / oppCellSize);
				int cellIntersectionIndices[2];
				cellIntersectionIndices[c] = s;
				cellIntersectionIndices[1 - c] = oppDimIndex;

				IntersectionInfo< GeometryReal > info;
				info.intersectionKey = SetIntersectionKey(atlasEdgeIndex, gridChart.gridIndexOffset + c*gridChart.width*gridChart.height + cellIntersectionIndices[0] + cellIntersectionIndices[1] * gridChart.width);

				info.time = alpha;
				info.position = position;
				intersections.push_back(info);
			}
		}
	}
}

template< typename GeometryReal >
void InitializeChartBoundaryEdgeGridIntersections
(
	const AtlasChart< GeometryReal > &atlasChart , GridChart< GeometryReal > &gridChart , std::unordered_map< int , int > &boundaryVerticesIndices , int &lastBoundaryIndex,  int numInteriorNodes,
	std::unordered_map< int, std::vector< IntersectionInfo< GeometryReal > > > &boundaryEdgeIntersections ,
	std::unordered_map< unsigned long long , BoundarySegmentInfo< GeometryReal > > &localBoundarySegmentsInfo ,
	std::unordered_map< unsigned long long , int > & localBoundaryNodeIndices ,
	std::unordered_map< unsigned long long , Point2D< GeometryReal > > &localBoundaryNodePosition
)
{

	for (int b = 0; b < atlasChart.boundaryHalfEdges.size(); b++) {
		int chartHalfEdgeIndex = atlasChart.boundaryHalfEdges[b];
		int tIndex = chartHalfEdgeIndex / 3;
		int kIndex = chartHalfEdgeIndex % 3;

		int halfEdgeIndex = atlasChart.meshTriangleIndices[tIndex] * 3 + kIndex;

		std::vector< IntersectionInfo< GeometryReal > > edgeIntersectionsInfo;
		Point2D< GeometryReal > edge[2];
		edge[0] = atlasChart.vertices[atlasChart.triangles[tIndex][kIndex]] - gridChart.corner;
		edge[1] = atlasChart.vertices[atlasChart.triangles[tIndex][(kIndex + 1) % 3]] - gridChart.corner;

		for (int v = 0; v < 2; v++) {
			//unsigned long long vertexKey = SetIntersectionKey(atlasMesh.triangles[tIndex][(kIndex + v) % 3], -1);
			unsigned long long vertexKey = SetIntersectionKey(atlasChart.triangles[tIndex][(kIndex + v) % 3], 4294967295);
			IntersectionInfo< GeometryReal > cornerIntersection;
			cornerIntersection.intersectionKey = vertexKey;
			cornerIntersection.time = (GeometryReal)v;
			cornerIntersection.position = edge[v];
			edgeIntersectionsInfo.push_back(cornerIntersection);
		}

		ComputeAtlasEdgeGridIntersection( edge , gridChart , atlasChart.atlasEdgeIndices[chartHalfEdgeIndex] , edgeIntersectionsInfo );

		//Sort intersections
		std::sort( edgeIntersectionsInfo.begin() , edgeIntersectionsInfo.end() , IntersectionComparison< GeometryReal > );

		for (int i = 0; i < edgeIntersectionsInfo.size(); i++) {
			if (localBoundaryNodeIndices.find(edgeIntersectionsInfo[i].intersectionKey) == localBoundaryNodeIndices.end()) {
				localBoundaryNodePosition[edgeIntersectionsInfo[i].intersectionKey] = edgeIntersectionsInfo[i].position;
				unsigned long tElement, vElement;
				GetIntersectionKey(edgeIntersectionsInfo[i].intersectionKey, tElement, vElement);
				if (vElement == 4294967295) {//boundary vertex
					int vertexIndex = atlasChart.meshVertexIndices[tElement];
					if( boundaryVerticesIndices.find(vertexIndex)==boundaryVerticesIndices.end() ) Miscellany::Throw( "Boundary vertex not found" );
					localBoundaryNodeIndices[edgeIntersectionsInfo[i].intersectionKey] = boundaryVerticesIndices[vertexIndex];

					AuxiliaryNode< GeometryReal > auxNode;
					auxNode.index = boundaryVerticesIndices[vertexIndex] + numInteriorNodes;
					auxNode.position = edgeIntersectionsInfo[i].position;
					gridChart.auxiliaryNodes.push_back(auxNode);

				}
				else {
					localBoundaryNodeIndices[edgeIntersectionsInfo[i].intersectionKey] = lastBoundaryIndex;
					AuxiliaryNode< GeometryReal > auxNode;
					auxNode.index = lastBoundaryIndex + numInteriorNodes;
					auxNode.position = edgeIntersectionsInfo[i].position;
					gridChart.auxiliaryNodes.push_back(auxNode);

					lastBoundaryIndex++;
				}
			}
			edgeIntersectionsInfo[i].intersectionIndex = localBoundaryNodeIndices[edgeIntersectionsInfo[i].intersectionKey];
		}

		for (int i = 0; i < edgeIntersectionsInfo.size() - 1; i++) {
			int segmentCornerIndices[2] = { localBoundaryNodeIndices[edgeIntersectionsInfo[i].intersectionKey], localBoundaryNodeIndices[edgeIntersectionsInfo[i + 1].intersectionKey] };
			unsigned long long segmentKey = SetMeshEdgeKey(segmentCornerIndices[0], segmentCornerIndices[1]);
			BoundarySegmentInfo< GeometryReal > segmentInfo;
			segmentInfo.startTime = edgeIntersectionsInfo[i].time;
			segmentInfo.endTime = edgeIntersectionsInfo[i + 1].time;
			segmentInfo.halfEdge = halfEdgeIndex;

			if( localBoundarySegmentsInfo.find(segmentKey)!=localBoundarySegmentsInfo.end() ) Miscellany::Throw( "Replicated segment key" );
			localBoundarySegmentsInfo[segmentKey] = segmentInfo;
		}
		boundaryEdgeIntersections[halfEdgeIndex] = edgeIntersectionsInfo;
	}
}

template< typename GeometryReal >
void IndexedPolygonFromCell( int i , int j , const GridChart< GeometryReal > &gridChart , IndexedIntersectionPolygon< GeometryReal > &polygon )
{
	int width = gridChart.width;
	int height = gridChart.height;
	int gridIndexOffset = gridChart.gridIndexOffset;

	polygon.vertices.resize(4);
	polygon.indices.resize(4);
	polygon.edgeIndices.resize(4);

	polygon.vertices[0] = Point2D< GeometryReal >( (GeometryReal)(i+0)*gridChart.cellSizeW , (GeometryReal)(j+0)*gridChart.cellSizeH );
	polygon.vertices[1] = Point2D< GeometryReal >( (GeometryReal)(i+1)*gridChart.cellSizeW , (GeometryReal)(j+0)*gridChart.cellSizeH );
	polygon.vertices[2] = Point2D< GeometryReal >( (GeometryReal)(i+1)*gridChart.cellSizeW , (GeometryReal)(j+1)*gridChart.cellSizeH );
	polygon.vertices[3] = Point2D< GeometryReal >( (GeometryReal)(i+0)*gridChart.cellSizeW , (GeometryReal)(j+1)*gridChart.cellSizeH );

	polygon.indices[0] = SetIntersectionKey( 4294967295 , gridIndexOffset + (i+0) + (j+0)*width );
	polygon.indices[1] = SetIntersectionKey( 4294967295 , gridIndexOffset + (i+1) + (j+0)*width );
	polygon.indices[2] = SetIntersectionKey( 4294967295 , gridIndexOffset + (i+1) + (j+1)*width );
	polygon.indices[3] = SetIntersectionKey( 4294967295 , gridIndexOffset + (i+0) + (j+1)*width );

	polygon.edgeIndices[0] = gridIndexOffset + width*height + i + j* width;
	polygon.edgeIndices[1] = gridIndexOffset + i + 1 + j* width;
	polygon.edgeIndices[2] = gridIndexOffset + width*height + i + (j + 1)* width;
	polygon.edgeIndices[3] = gridIndexOffset + i + j* width;
}


template< typename GeometryReal >
void InitializeChartBoundaryPolygons
(
	const std::vector< int > &atlasEdgeIndex , const std::vector< int > &oppositeHalfEdge ,
	const AtlasChart< GeometryReal > &atlasChart , GridChart< GeometryReal > &gridChart ,
	int numInteriorNodes , int numBoundaryVertices , int numBoundaryNodes ,
#if 1
//#pragma message( "[WARNING] making a fix here" )
	std::unordered_map< int , std::vector< IntersectionInfo< GeometryReal > > > &boundaryEdgeIntersections ,
#else
	std::unordered_map< int , std::vector< IntersectionInfo< GeometryReal > > > boundaryEdgeIntersections ,
#endif
	std::unordered_map< unsigned long long , BoundarySegmentInfo< GeometryReal > > &boundarySegmentsInfo ,
	std::unordered_map< unsigned long long , int > &boundaryNodeIndices ,
	std::unordered_map< unsigned long long , Point2D< GeometryReal > > & boundaryNodePosition ,
	std::vector< int > &coveredOppositeBoundaryNode
)
{

	std::unordered_map<int, std::vector<std::pair<unsigned long long, unsigned long long>>> cellSegments;

	//(1) Generate cell segments
	for (int t = 0; t < atlasChart.triangles.size(); t++)
	{
		Point2D< GeometryReal > tPos[3];
		for (int i = 0; i < 3; i++) tPos[i] = atlasChart.vertices[atlasChart.triangles[t][i]] - gridChart.corner;
		int minCorner[2];
		int maxCorner[2];
		GetTriangleIntegerBBox( tPos, (GeometryReal)1./gridChart.cellSizeW , (GeometryReal)1./gridChart.cellSizeH , minCorner , maxCorner );

		IndexedIntersectionTriangle< GeometryReal > indexedTriangle;
		for (int k = 0; k < 3; k++) {
			indexedTriangle.vertices[k] = tPos[k];
			indexedTriangle.indices[k] = SetIntersectionKey(atlasChart.triangles[t][k], 4294967295);
			indexedTriangle.edgeIndices[k] = atlasChart.atlasEdgeIndices[3 * t + k];
		}

		for (int j = minCorner[1]; j < maxCorner[1]; j++) for (int i = minCorner[0]; i < maxCorner[0]; i++) {
			if (gridChart.cellType(i, j) == 0) {
				int cellID = gridChart.localBoundaryCellIndex(i,j);
				if( cellID==-1 ) Miscellany::Throw( "Boundary cell invalid ID" );
				IndexedIntersectionPolygon< GeometryReal > cellPolygon;
				IndexedPolygonFromCell(i, j, gridChart, cellPolygon);
				if (ClipIndexedIntersectionPolygonToIndexedIntersectionTriangle(cellPolygon, indexedTriangle)) {
					for (int s = 0; s < cellPolygon.indices.size(); s++){
						std::pair<unsigned long long, unsigned long long> edge(cellPolygon.indices[s], cellPolygon.indices[(s + 1) % cellPolygon.indices.size()]);
						cellSegments[cellID].push_back(edge);
					}
				}
			}
		}
	}

	//(2) Process cells

	gridChart.boundaryPolygons.resize(gridChart.numBoundaryCells);
	for (auto cellIter = cellSegments.begin(); cellIter != cellSegments.end(); cellIter++) {
		std::vector<std::pair<unsigned long long, unsigned long long>> segments = (*cellIter).second;
		int cellID = (*cellIter).first;

		//Normalize vertex indices
		int lastCellVertex = 0;
		std::unordered_map<unsigned long long, int> cellVerticesIndex;
		std::vector<std::pair<int, int>> normalizedSegment;
		std::vector<unsigned long long> normalizedVertexKey;
		for (int k = 0; k < segments.size(); k++) {
			if (cellVerticesIndex.find(segments[k].first) == cellVerticesIndex.end()) {
				cellVerticesIndex[segments[k].first] = lastCellVertex;
				normalizedVertexKey.push_back(segments[k].first);
				lastCellVertex++;
			}
			if (cellVerticesIndex.find(segments[k].second) == cellVerticesIndex.end()) {
				cellVerticesIndex[segments[k].second] = lastCellVertex;
				normalizedVertexKey.push_back(segments[k].second);
				lastCellVertex++;
			}
			normalizedSegment.push_back(std::pair<int, int>(cellVerticesIndex[segments[k].first], cellVerticesIndex[segments[k].second]));
		}

		std::unordered_set<unsigned long long> segmentKeys;
		for (int k = 0; k < segments.size(); k++) segmentKeys.insert(SetMeshEdgeKey(normalizedSegment[k].first, normalizedSegment[k].second));

		std::unordered_map<unsigned long long, unsigned long long> forwardMap;
		for (int k = 0; k < segments.size(); k++) {
			unsigned long long oppositeKey = SetMeshEdgeKey(normalizedSegment[k].second, normalizedSegment[k].first);
			if (segmentKeys.find(oppositeKey) == segmentKeys.end()) {//Boundary segment
				forwardMap[segments[k].first] = segments[k].second;
			}
		}

		std::vector< std::vector< unsigned long long > > loopKeys;
		LoopVertices( forwardMap , loopKeys );

		//Keep only the grid nodes and the boundary vertices
		std::vector< std::vector< Point2D< GeometryReal > > > loopPositions(loopKeys.size());
		std::vector< std::vector< int> > loopIndices(loopKeys.size());
		std::vector< std::vector< int> > loopAtlasVertexIndices(loopKeys.size());
		for (int i = 0; i < loopKeys.size(); i++) {
			std::vector< unsigned long long > reduceLoop;
			std::vector< Point2D< GeometryReal > > currentLoopPositions;
			std::vector< int > currentLoopIndices;
			std::vector< int > currentAtlasVertexIndices;
			for (int j = 0; j < loopKeys[i].size(); j++) {
				unsigned long long currentVertexKey = loopKeys[i][j];
				unsigned long tElement, vElement;
				GetIntersectionKey(currentVertexKey, tElement, vElement);
				if (tElement == 4294967295) {
					vElement -= gridChart.gridIndexOffset;
					int pi = vElement % gridChart.width;
					int pj = vElement / gridChart.width;
					int interiorTexelIndex = gridChart.globalInteriorTexelIndex(pi, pj);
					if( interiorTexelIndex==-1 ) Miscellany::Throw( "Invalid texel" );
					reduceLoop.push_back( currentVertexKey );
					currentLoopPositions.push_back( Point2D< GeometryReal >( (GeometryReal)pi*gridChart.cellSizeW , (GeometryReal)pj*gridChart.cellSizeH ) );
					currentLoopIndices.push_back( interiorTexelIndex );
					currentAtlasVertexIndices.push_back( -1 );
				}
				else if (boundaryNodeIndices.find(currentVertexKey) != boundaryNodeIndices.end()) {
					reduceLoop.push_back(currentVertexKey);
					currentLoopPositions.push_back(boundaryNodePosition[currentVertexKey]);
					currentLoopIndices.push_back(boundaryNodeIndices[currentVertexKey] + numInteriorNodes);

					if (vElement == 4294967295) currentAtlasVertexIndices.push_back(tElement);
					else currentAtlasVertexIndices.push_back(-1);
				}
			}
			if( reduceLoop.size()==0 ) Miscellany::Throw( "Reduced loop cannot be empty. Original loop size %d" , (int)loopKeys[i].size() );
			loopKeys[i] = reduceLoop;
			loopPositions[i] = currentLoopPositions;
			loopIndices[i] = currentLoopIndices;
			loopAtlasVertexIndices[i] = currentAtlasVertexIndices;
		}

		//Add the intermediate vertices for boundary segments
		std::vector<std::vector<int>> loopAtlasEdges(loopKeys.size());
		std::vector<std::vector<int>> loopAtlasVertexParentEdges(loopKeys.size());
		for (int i = 0; i < loopKeys.size(); i++) {

			std::vector< std::vector< int > > indicesToInsert;
			std::vector< std::vector< GeometryReal > > timesToInsert;
			std::vector< bool > isInsertionPos( loopKeys[i].size() , false );
			std::vector< int > polygonAtlaEdgeIndex( loopKeys[i].size() , -1 );
			std::vector< int > polygonAtlaVertexParentEdges( loopKeys[i].size() , -1 );

			for (int j = 0; j < loopKeys[i].size(); j++) {
				unsigned long long currentVertexKey = loopKeys[i][j];
				unsigned long long nextVertexKey = loopKeys[i][(j + 1) % loopKeys[i].size()];

				unsigned long tElement, vElement;
				GetIntersectionKey(currentVertexKey, tElement, vElement);
				if (tElement != 4294967295 && vElement != 4294967295) {
					polygonAtlaVertexParentEdges[j] = tElement;
				}

				if (boundaryNodeIndices.find(currentVertexKey) != boundaryNodeIndices.end() && boundaryNodeIndices.find(nextVertexKey) != boundaryNodeIndices.end()) {
					unsigned long long segmentKey = SetMeshEdgeKey(boundaryNodeIndices[currentVertexKey], boundaryNodeIndices[nextVertexKey]);
					if (boundarySegmentsInfo.find(segmentKey) != boundarySegmentsInfo.end()) {
						BoundarySegmentInfo< GeometryReal > segmentInfo = boundarySegmentsInfo[segmentKey];
						polygonAtlaEdgeIndex[j] = atlasEdgeIndex[segmentInfo.halfEdge];

						int oppHalfEdge = oppositeHalfEdge[segmentInfo.halfEdge];
						if( oppHalfEdge!=-1 )
						{
							GeometryReal startTime = segmentInfo.startTime;
							GeometryReal endTime = segmentInfo.endTime;


							if( boundaryEdgeIntersections.find(oppHalfEdge)==boundaryEdgeIntersections.end() )
								Miscellany::Throw( "Opposite edge intersections not found. Current  edge %d. Opposite %d" , segmentInfo.halfEdge , oppHalfEdge );

							std::vector< int > segmentIndicesToInsert;
							std::vector< GeometryReal > segmentTimesToInsert;

							std::vector< IntersectionInfo< GeometryReal > > oppositeIntersection = boundaryEdgeIntersections[oppHalfEdge];
							for (int oppIter = (int)oppositeIntersection.size() - 2; oppIter > 0; oppIter--)
							{
								GeometryReal reversedTime = (GeometryReal)1. - oppositeIntersection[oppIter].time;
								GeometryReal normalizedTime = ( reversedTime-startTime ) / ( endTime-startTime );
								if (normalizedTime > 0 && normalizedTime < 1)
								{
									segmentIndicesToInsert.push_back(oppositeIntersection[oppIter].intersectionIndex);
									segmentTimesToInsert.push_back(normalizedTime);
								}
							}


							if (segmentIndicesToInsert.size() > 0)
							{
								isInsertionPos[j] = true;
								indicesToInsert.push_back(segmentIndicesToInsert);
								timesToInsert.push_back(segmentTimesToInsert);
							}
						}
					}
				}
			}
			//Do insertions
			int insertionCount = 0;
			std::vector< Point2D< GeometryReal > > expandedLoopPositions;
			std::vector< int > expandedLoopIndices;
			std::vector< int > expandedLoopAtlasVertexIndices;
			std::vector< int > expandedLoopAtlasEdgeIndex;
			std::vector< int > expandedVertexParentEdgeIndex;
			for (int j = 0; j < loopKeys[i].size(); j++) {

				//expandedLoopKeys.push_back(loopKeys[i][j]);
				expandedLoopPositions.push_back(loopPositions[i][j]);
				expandedLoopIndices.push_back(loopIndices[i][j]);
				expandedLoopAtlasVertexIndices.push_back(loopAtlasVertexIndices[i][j]);

				int currentSegmentAtlasEdgeIndex = polygonAtlaEdgeIndex[j];
				expandedLoopAtlasEdgeIndex.push_back(currentSegmentAtlasEdgeIndex);
				expandedVertexParentEdgeIndex.push_back(polygonAtlaVertexParentEdges[j]);

				if (isInsertionPos[j])
				{
					std::vector< int > _indices = indicesToInsert[insertionCount];
					std::vector< GeometryReal > _times = timesToInsert[insertionCount];
					Point2D< GeometryReal > currentPos = loopPositions[i][j];
					Point2D< GeometryReal > nextPos = loopPositions[i][ (j+1)%loopKeys[i].size() ];

					if( currentSegmentAtlasEdgeIndex==-1 ) Miscellany::Throw( "Invalid atlas edge index" );
					for (int k = 0; k < _indices.size(); k++)
					{
						GeometryReal alpha = _times[k];
						Point2D< GeometryReal > interpolatedPos = (GeometryReal)( 1.-alpha )*currentPos + nextPos*alpha;
						if( true ) // Add orthogonal perturbation to avoid colinearity. Delanauy triangulation library failed for colinear inputs on some systems.
						{
							Point2D< GeometryReal > segmentDir = nextPos - currentPos;
							interpolatedPos += Point2D< GeometryReal >( -segmentDir[1] , segmentDir[0] ) * (GeometryReal)( ( Random< GeometryReal >()*2 - 1. ) * 1e-10 );
						}
						if( _indices[k]<numBoundaryVertices || _indices[k]>numBoundaryNodes ) Miscellany::Throw( "Out of bounds index: %d not in %d %d" , _indices[k] , numBoundaryVertices , numBoundaryNodes );
						coveredOppositeBoundaryNode[_indices[k] - numBoundaryVertices] += 1;

						AuxiliaryNode< GeometryReal > auxNode;
						auxNode.index = _indices[k] + numInteriorNodes;
						auxNode.position = interpolatedPos;
						gridChart.auxiliaryNodes.push_back(auxNode);

						expandedLoopPositions.push_back(interpolatedPos);
						expandedLoopIndices.push_back(_indices[k] + numInteriorNodes);
						expandedLoopAtlasVertexIndices.push_back(-1);
						expandedLoopAtlasEdgeIndex.push_back(currentSegmentAtlasEdgeIndex);
						expandedVertexParentEdgeIndex.push_back(currentSegmentAtlasEdgeIndex);
					}
					insertionCount++;
				}
			}
			if( insertionCount!=indicesToInsert.size() ) Miscellany::Throw( "Intersection chains count does not match" );
			loopPositions[i] = expandedLoopPositions;
			loopIndices[i] = expandedLoopIndices;
			loopAtlasVertexIndices[i] = expandedLoopAtlasVertexIndices;
			loopAtlasEdges[i] = expandedLoopAtlasEdgeIndex;
			loopAtlasVertexParentEdges[i] = expandedVertexParentEdgeIndex;
		}

		for (int i = 0; i < loopPositions.size(); i++){
			AtlasIndexedPolygon< GeometryReal > poly;
			poly.vertices = loopPositions[i];
			poly.indices = loopIndices[i];
			poly.atlasVertexIndices = loopAtlasVertexIndices[i];
			poly.atlasEdgeIndices = loopAtlasEdges[i];
			poly.atlasVertexParentEdge = loopAtlasVertexParentEdges[i];
			gridChart.boundaryPolygons[cellID].push_back(poly);
		}
	}
}


template< typename GeometryReal >
void InitializeBoundaryPolygons
(
	const std::vector< int > &atlasEdgeIndex , const std::vector< int > &oppositeHalfEdge ,
	const std::vector< AtlasChart< GeometryReal > > &atlasCharts , std::vector< GridChart< GeometryReal > > &gridCharts ,
	std::unordered_map< int , int > &boundaryVerticesIndices ,
	int numInteriorNodes , int numBoundaryVertices , int &numBoundaryNodes , const bool isClosedMesh
)
{ //Fine System

	int lastBoundaryIndex = numBoundaryVertices;
	std::unordered_map< int , std::vector< IntersectionInfo< GeometryReal > > > boundaryEdgeIntersections;

	std::vector< std::unordered_map< unsigned long long , BoundarySegmentInfo< GeometryReal > > > localBoundarySegmentsInfo( gridCharts.size() );
	std::vector< std::unordered_map< unsigned long long , int > > localBoundaryNodeIndices( gridCharts.size() );
	std::vector< std::unordered_map< unsigned long long , Point2D< GeometryReal > > > localBoundaryNodePosition( gridCharts.size() );

	for( int i=0 ; i<gridCharts.size() ; i++ )
		InitializeChartBoundaryEdgeGridIntersections( atlasCharts[i] , gridCharts[i] , boundaryVerticesIndices , lastBoundaryIndex , numInteriorNodes , boundaryEdgeIntersections , localBoundarySegmentsInfo[i] , localBoundaryNodeIndices[i] , localBoundaryNodePosition[i] );

	numBoundaryNodes = lastBoundaryIndex;
	std::vector<int> coveredOppositeBoundaryNode(numBoundaryNodes - numBoundaryVertices, 0);
	for( int i=0 ; i<gridCharts.size() ; i++ )
		InitializeChartBoundaryPolygons( atlasEdgeIndex , oppositeHalfEdge , atlasCharts[i] , gridCharts[i] , numInteriorNodes , numBoundaryVertices , numBoundaryNodes , boundaryEdgeIntersections , localBoundarySegmentsInfo[i] , localBoundaryNodeIndices[i] , localBoundaryNodePosition[i] , coveredOppositeBoundaryNode );

	if( isClosedMesh ) for ( int i=0 ; i<coveredOppositeBoundaryNode.size() ; i++ ) if( coveredOppositeBoundaryNode[i]!=1 )
		Miscellany::Warn( "Non-opposite boundary node at node %d" , i );
}

#include "ConstrainedTriangulation.h"

template< typename GeometryReal >
void InitializeChartQuadraticElements( GridChart< GeometryReal > &gridChart , std::unordered_map< unsigned long long , int > &midPointIndex , int &lastMidPointIndex , int previouslyAddedNodes )
{
	const std::vector< std::vector< AtlasIndexedPolygon< GeometryReal > > > & boundaryPolygons = gridChart.boundaryPolygons;
	std::vector< std::vector< BoundaryIndexedTriangle< GeometryReal > > > & boundaryTriangles = gridChart.boundaryTriangles;
	boundaryTriangles.resize(gridChart.numBoundaryCells);
	int numBoundaryTriangles = 0;
	for (int i = 0; i < boundaryPolygons.size(); i++) {
		for (int j = 0; j < boundaryPolygons[i].size(); j++) {
			const AtlasIndexedPolygon< GeometryReal > &currentPolygon = boundaryPolygons[i][j];
			std::vector< TriangleIndex > delanauyTriangles;
			TriangulatePolygon( currentPolygon.vertices , delanauyTriangles );
			if( delanauyTriangles.size()!=currentPolygon.vertices.size()-2 )
			{
				int localCellPos[2] = { -1,-1 };
				for( int li=0 ; li<gridChart.localBoundaryCellIndex.width() ; li++ ) for( int lj=0 ; lj<gridChart.localBoundaryCellIndex.height() ; lj++ ) if( gridChart.localBoundaryCellIndex(li,lj) == i ) localCellPos[0] = li , localCellPos[1] = lj;
				Miscellany::Warn( "Unexpected number of triangles produced by delaunay triangulation at global cell (%d,%d). Polygon may self intersect!" , gridChart.cornerCoords[0] + localCellPos[0] , gridChart.cornerCoords[1] + localCellPos[1] );
			}

			for (int k = 0; k < delanauyTriangles.size(); k++) {
				BoundaryIndexedTriangle< GeometryReal > triangle;
				for (int v = 0; v < 3; v++) {
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

				for (int v = 0; v < 3; v++)
				{
					Point2D< GeometryReal > edgeMidPoint = ( triangle.vertices[ (v+1)%3 ] + triangle.vertices[ (v+2)%3 ] ) / 2;
					int edgeCorners[2] = { (int)triangle.indices[(v + 1) % 3], (int)triangle.indices[(v + 2) % 3] };
					if (edgeCorners[0] > edgeCorners[1]) std::swap(edgeCorners[0], edgeCorners[1]);
					unsigned long long edgeKey = SetMeshEdgeKey(edgeCorners[0], edgeCorners[1]);

					if (midPointIndex.find(edgeKey) == midPointIndex.end()) {
						midPointIndex[edgeKey] = lastMidPointIndex;
						midVexterIndex[v] = lastMidPointIndex + previouslyAddedNodes;
						lastMidPointIndex++;
					}
					else {
						midVexterIndex[v] = midPointIndex[edgeKey] + previouslyAddedNodes;
					}

					AuxiliaryNode< GeometryReal > auxNode;
					auxNode.index = midVexterIndex[v];
					auxNode.position = edgeMidPoint;
					gridChart.auxiliaryNodes.push_back(auxNode);

				}
				triangle.id = numBoundaryTriangles;
				triangle.indices[3] = midVexterIndex[0];
				triangle.indices[4] = midVexterIndex[1];
				triangle.indices[5] = midVexterIndex[2];
				boundaryTriangles[i].push_back(triangle);
				numBoundaryTriangles++;
			}
		}
	}
	gridChart.numBoundaryTriangles = numBoundaryTriangles;
}

template< typename GeometryReal >
void InitializeBoundaryQuadraticElements( std::vector< GridChart< GeometryReal > > &gridCharts , int previouslyAddedNodes , int & numMidPoints )
{
	int lastMidPointIndex = 0;
	std::unordered_map< unsigned long long , int > midPointIndex;
	for( int i=0 ; i<gridCharts.size() ; i++ ) InitializeChartQuadraticElements( gridCharts[i] , midPointIndex , lastMidPointIndex , previouslyAddedNodes );
	numMidPoints = lastMidPointIndex;
}

template< typename GeometryReal , typename MatrixReal >
void InitializeBoundaryTriangulation( GridAtlas< GeometryReal , MatrixReal > &gridAtlas , AtlasMesh< GeometryReal > &atlasMesh , std::vector< AtlasChart< GeometryReal > > &atlasCharts , std::vector< int > &oppositeHalfEdge , std::unordered_map< int , int > &boundaryVerticesIndices , int numBoundaryVertices , const bool &isClosedMesh , bool verbose=false )
{
	InitializeBoundaryPolygons( atlasMesh.halfEdgeToEdgeIndex , oppositeHalfEdge , atlasCharts , gridAtlas.gridCharts , boundaryVerticesIndices , gridAtlas.numInteriorTexels , numBoundaryVertices , gridAtlas.numBoundaryNodes , isClosedMesh );
	InitializeBoundaryQuadraticElements( gridAtlas.gridCharts , gridAtlas.numBoundaryNodes + gridAtlas.numInteriorTexels , gridAtlas.numMidPoints );

	if( verbose )
	{
		printf( "\t Boundary vertices %d \n" , numBoundaryVertices );
		printf( "\t Boundary nodes %d \n" , gridAtlas.numBoundaryNodes );
		printf( "\t Boundary cells %d \n" , gridAtlas.numBoundaryCells );
		printf( "\t Mid points %d \n" , gridAtlas.numMidPoints );
	}

	gridAtlas.numFineNodes = gridAtlas.numInteriorTexels + gridAtlas.numBoundaryNodes + gridAtlas.numMidPoints;

	if( verbose ) printf( "Num fine nodes %d =  Interior texels %d + Auxiliar boundary nodes %d\n" , gridAtlas.numFineNodes , gridAtlas.numInteriorTexels , gridAtlas.numBoundaryNodes + gridAtlas.numMidPoints );
}

template< typename MatrixReal >
class BoundaryProlongationData
{
public:
	std::vector< int > fineBoundaryIndex;
	int numFineBoundarNodes;
	SparseMatrix< MatrixReal , int > coarseBoundaryFineBoundaryProlongation;
	SparseMatrix< MatrixReal , int > fineBoundaryCoarseBoundaryRestriction;
};

template< typename GeometryReal , typename MatrixReal >
void InitializeCoarseBoundaryToFineBoundaryProlongation( const GridAtlas< GeometryReal , MatrixReal > &gridAtlas , SparseMatrix< MatrixReal , int > &coarseBoundaryFineBoundaryProlongation , std::vector< int > &fineBoundaryIndex , int &numFineBoundaryNodes , bool verbose=false )
{

	std::vector< int > boundaryFineToFullFine;
	std::vector< Eigen::Triplet< MatrixReal > > prolongationTriplets;
	const std::vector< int > & boundaryAndDeepIndex = gridAtlas.boundaryAndDeepIndex;


	fineBoundaryIndex.resize(gridAtlas.numFineNodes, -1);
	int lastFineBoundaryIndex = 0;
	for (int i = 0; i < gridAtlas.gridCharts.size(); i++)
	{
		const GridChart< GeometryReal > &gridChart = gridAtlas.gridCharts[i];
		for (int j = 0; j < gridChart.globalInteriorTexelIndex.size(); j++) {
			if (gridChart.globalInteriorTexelIndex[j] != -1 && gridChart.globalTexelDeepIndex[j] == -1) { //Interior but not deep
				int coarseGlobalIndex = gridChart.globalTexelIndex[j];

				int coarseBondaryIndex = boundaryAndDeepIndex[coarseGlobalIndex] - 1;

				if( coarseBondaryIndex>=0 ) prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( lastFineBoundaryIndex , coarseBondaryIndex , (MatrixReal)1. ) );
				else Miscellany::Throw( "Coarse node is not boundary. Global index %d. Boundary index %d" , coarseGlobalIndex , coarseBondaryIndex );
				fineBoundaryIndex[gridChart.globalInteriorTexelIndex[j]] = lastFineBoundaryIndex;
				lastFineBoundaryIndex++;
				boundaryFineToFullFine.push_back(gridChart.globalInteriorTexelIndex[j]);
			}
		}
	}
	for (int i = gridAtlas.numInteriorTexels; i < gridAtlas.numFineNodes;i++) {
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

	for (int i = 0; i < gridCharts.size(); i++){
		const GridChart< GeometryReal > &gridChart = gridCharts[i];
		for (int j = 0; j < gridChart.auxiliaryNodes.size(); j++) {
			int auxiliaryID = gridChart.auxiliaryNodes[j].index - numInteriorTexels;
			int fineBoundaryID = fineBoundaryIndex[gridChart.auxiliaryNodes[j].index];
			int nodeDegree = auxiliaryNodesDegree[auxiliaryID];
			Point2D< GeometryReal >nodePosition = gridChart.auxiliaryNodes[j].position;
			int corner[2] = { (int)floor(nodePosition[0] / gridChart.cellSizeW), (int)floor(nodePosition[1] / gridChart.cellSizeH) };
			int cellID = gridChart.localCellIndex(corner[0], corner[1]);
			if( cellID==-1 ) Miscellany::Throw( "Invalid cell index. Node position %f %f" , nodePosition[0] / gridChart.cellSizeW , nodePosition[1] / gridChart.cellSizeH );
			nodePosition[0] /= gridChart.cellSizeW;
			nodePosition[1] /= gridChart.cellSizeH;
			nodePosition[0] -= (GeometryReal)corner[0];
			nodePosition[1] -= (GeometryReal)corner[1];
			if( nodePosition[0] < 0-precision_error || nodePosition[0] > 1+precision_error || nodePosition[1] < 0-precision_error || nodePosition[1] > 1+precision_error )
				Miscellany::Throw( "Sample out of unit box! (%f %f)" , nodePosition[0] , nodePosition[1] );
			for (int k = 0; k < 4; k++)
			{
				MatrixReal texelWeight = (MatrixReal)BilinearElementValue(k, nodePosition) / nodeDegree;
				if (fabs(texelWeight) > 1e-11) {
					auxiliaryNodesCumWeight[auxiliaryID] += texelWeight;
					int texelIndex = gridChart.bilinearElementIndices[cellID][k];
					if( nodeInfo[texelIndex].nodeType==2 ) Miscellany::Throw( "Deep texel cannot be in the support of an auxiliary node. Weight %g (A)" , texelWeight );
					if( gridChart.auxiliaryNodes[j].index<numInteriorTexels || gridChart.auxiliaryNodes[j].index>numFineNodes || texelIndex<0 || texelIndex>numCoarseNodes ) Miscellany::Throw( "Out of bounds index" );

					int coarseBoundaryID = boundaryAndDeepIndex[texelIndex] - 1;

					if( coarseBoundaryID<0 ) Miscellany::Throw( "Coarse node is not boundary" );
					prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( fineBoundaryID , coarseBoundaryID , texelWeight ) );
				}
			}
		}
	}

	int numBoundaryTexels = gridAtlas.numBoundaryTexels;
	coarseBoundaryFineBoundaryProlongation = SetSparseMatrix( prolongationTriplets , numFineBoundaryNodes , numBoundaryTexels , false );
}

template< typename GeometryReal , typename MatrixReal >
void InitializeBoundaryProlongationData( const GridAtlas< GeometryReal , MatrixReal > &gridAtlas , BoundaryProlongationData< MatrixReal > &boundaryProlongation )
{
	InitializeCoarseBoundaryToFineBoundaryProlongation( gridAtlas , boundaryProlongation.coarseBoundaryFineBoundaryProlongation , boundaryProlongation.fineBoundaryIndex , boundaryProlongation.numFineBoundarNodes );
	boundaryProlongation.fineBoundaryCoarseBoundaryRestriction = boundaryProlongation.coarseBoundaryFineBoundaryProlongation.transpose();
}