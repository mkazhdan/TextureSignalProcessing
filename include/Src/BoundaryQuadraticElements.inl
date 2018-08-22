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
#include "PolygonClipping.h"
#include "MapLoop.h"


void SetAtlasIndexedPolygonFromBoundaryTriangle(const BoundaryIndexedTriangle & triangle, AtlasIndexedPolygon & polygon) {
	for (int k = 0; k < 3; k++) {
		polygon.vertices.push_back(triangle.vertices[k]);
		polygon.indices.push_back(triangle.indices[k]);
		polygon.atlasVertexIndices.push_back(triangle.atlasVertexIndices[k]);
		polygon.atlasEdgeIndices.push_back(triangle.atlasEdgeIndices[k]);
		polygon.atlasVertexParentEdge.push_back(triangle.atlasVertexParentEdge[k]);
	}
}
void SetIndexedPolygonFromCell(const int i, const int j, const int width, const int height, IndexedIntersectionPolygon & polygon) {
	//
	polygon.vertices.resize(4);
	polygon.indices.resize(4);
	polygon.edgeIndices.resize(4);

	polygon.vertices[0] = Point2D<double>((double(i) + 0.5) / double(width), (double(j) + 0.5) / double(height));
	polygon.vertices[1] = Point2D<double>((double(i) + 1.5) / double(width), (double(j) + 0.5) / double(height));
	polygon.vertices[2] = Point2D<double>((double(i) + 1.5) / double(width), (double(j) + 1.5) / double(height));
	polygon.vertices[3] = Point2D<double>((double(i) + 0.5) / double(width), (double(j) + 1.5) / double(height));

	polygon.indices[0] = SetIntersectionKey(4294967295, i + j*width);
	polygon.indices[1] = SetIntersectionKey(4294967295, i + 1 + j*width);
	polygon.indices[2] = SetIntersectionKey(4294967295, i + 1 + (j + 1)*width);
	polygon.indices[3] = SetIntersectionKey(4294967295, i + (j + 1)*width);

	polygon.edgeIndices[0] = width*height + i + j* width;
	polygon.edgeIndices[1] = i + 1 + j* width;
	polygon.edgeIndices[2] = width*height + i + (j + 1)* width;
	polygon.edgeIndices[3] = i + j* width;
}


int ComputeAtlasEdgeGridIntersection(Point2D<double> p[2], const GridChart & gridChart, int atlasEdgeIndex, std::vector<IntersectionInfo> & intersections) {
	double fmin[2];
	double fmax[2];
	for (int c = 0; c < 2; c++) {
		fmin[c] = std::min< double >(p[0][c], p[1][c]);
		fmin[c] = std::max< double >(fmin[c], 0.f);
		fmax[c] = std::max< double >(p[0][c], p[1][c]);
		fmax[c] = std::min< double >(fmax[c], 1.f);
	}
	int imin[2];
	int imax[2];

	imin[0] = static_cast<int>(floor(fmin[0] / gridChart.cellSizeW));
	imin[1] = static_cast<int>(floor(fmin[1] / gridChart.cellSizeH));
	imax[0] = static_cast<int>(floor(fmax[0] / gridChart.cellSizeW));
	imax[1] = static_cast<int>(floor(fmax[1] / gridChart.cellSizeH));

	for (int c = 0; c < 2; c++) {
		for (int s = imin[c]; s <= imax[c]; s++) {
			double cellSize = c == 0 ? gridChart.cellSizeW : gridChart.cellSizeH;
			double oppCellSize = c == 0 ? gridChart.cellSizeH : gridChart.cellSizeW;
			double level = double(s)*cellSize;
			double alpha = (level - p[0][c]) / (p[1][c] - p[0][c]);
			if (alpha > 0 && alpha < 1.0) {
				Point2D<double> position = p[0] * (1.0 - alpha) + p[1] * alpha;
				int oppDimIndex = (int)floor(position[1 - c] / oppCellSize);
				int cellIntersectionIndices[2];
				cellIntersectionIndices[c] = s;
				cellIntersectionIndices[1 - c] = oppDimIndex;

				IntersectionInfo info;
				info.intersectionKey = SetIntersectionKey(atlasEdgeIndex, gridChart.gridIndexOffset + c*gridChart.width*gridChart.height + cellIntersectionIndices[0] + cellIntersectionIndices[1] * gridChart.width);

				info.time = alpha;
				info.position = position;
				intersections.push_back(info);
			}
		}
	}

	return 1;
}



int InitializeChartBoundaryEdgeGridIntersections(const AtlasChart & atlasChart, GridChart & gridChart, std::unordered_map<int, int> & boundaryVerticesIndices, int & lastBoundaryIndex, const int numInteriorNodes,
	std::unordered_map<int, std::vector<IntersectionInfo>> & boundaryEdgeIntersections,
	std::unordered_map<unsigned long long, BoundarySegmentInfo> & localBoundarySegmentsInfo,
	std::unordered_map<unsigned long long, int> & localBoundaryNodeIndices,
	std::unordered_map<unsigned long long, Point2D<double>> & localBoundaryNodePosition) {

	for (int b = 0; b < atlasChart.boundaryHalfEdges.size(); b++) {
		int chartHalfEdgeIndex = atlasChart.boundaryHalfEdges[b];
		int tIndex = chartHalfEdgeIndex / 3;
		int kIndex = chartHalfEdgeIndex % 3;

		int halfEdgeIndex = atlasChart.meshTriangleIndices[tIndex] * 3 + kIndex;

		std::vector<IntersectionInfo> edgeIntersectionsInfo;
		Point2D<double> edge[2];
		edge[0] = atlasChart.vertices[atlasChart.triangles[tIndex][kIndex]] - gridChart.corner;
		edge[1] = atlasChart.vertices[atlasChart.triangles[tIndex][(kIndex + 1) % 3]] - gridChart.corner;

		for (int v = 0; v < 2; v++) {
			//unsigned long long vertexKey = SetIntersectionKey(atlasMesh.triangles[tIndex][(kIndex + v) % 3], -1);
			unsigned long long vertexKey = SetIntersectionKey(atlasChart.triangles[tIndex][(kIndex + v) % 3], 4294967295);
			IntersectionInfo cornerIntersection;
			cornerIntersection.intersectionKey = vertexKey;
			cornerIntersection.time = double(v);
			cornerIntersection.position = edge[v];
			edgeIntersectionsInfo.push_back(cornerIntersection);
		}

		if (!ComputeAtlasEdgeGridIntersection(edge, gridChart, atlasChart.atlasEdgeIndices[chartHalfEdgeIndex], edgeIntersectionsInfo)) {
			printf("Failed edge-grid intersection");
			return 0;
		}

		//Sort intersections
		std::sort(edgeIntersectionsInfo.begin(), edgeIntersectionsInfo.end(), IntersectionComparison);

		for (int i = 0; i < edgeIntersectionsInfo.size(); i++) {
			if (localBoundaryNodeIndices.find(edgeIntersectionsInfo[i].intersectionKey) == localBoundaryNodeIndices.end()) {
				localBoundaryNodePosition[edgeIntersectionsInfo[i].intersectionKey] = edgeIntersectionsInfo[i].position;
				unsigned long tElement, vElement;
				GetIntersectionKey(edgeIntersectionsInfo[i].intersectionKey, tElement, vElement);
				if (vElement == 4294967295) {//boundary vertex
					int vertexIndex = atlasChart.meshVertexIndices[tElement];
					if (boundaryVerticesIndices.find(vertexIndex) == boundaryVerticesIndices.end()) {
						printf("Boundary vertex not found! \n");
						return 0;
					}
					localBoundaryNodeIndices[edgeIntersectionsInfo[i].intersectionKey] = boundaryVerticesIndices[vertexIndex];

					AuxiliarNode auxNode;
					auxNode.index = boundaryVerticesIndices[vertexIndex] + numInteriorNodes;
					auxNode.position = edgeIntersectionsInfo[i].position;
					gridChart.auxiliarNodes.push_back(auxNode);

				}
				else {
					localBoundaryNodeIndices[edgeIntersectionsInfo[i].intersectionKey] = lastBoundaryIndex;
					AuxiliarNode auxNode;
					auxNode.index = lastBoundaryIndex + numInteriorNodes;
					auxNode.position = edgeIntersectionsInfo[i].position;
					gridChart.auxiliarNodes.push_back(auxNode);

					lastBoundaryIndex++;
				}
			}
			edgeIntersectionsInfo[i].intersectionIndex = localBoundaryNodeIndices[edgeIntersectionsInfo[i].intersectionKey];
		}

		for (int i = 0; i < edgeIntersectionsInfo.size() - 1; i++) {
			int segmentCornerIndices[2] = { localBoundaryNodeIndices[edgeIntersectionsInfo[i].intersectionKey], localBoundaryNodeIndices[edgeIntersectionsInfo[i + 1].intersectionKey] };
			unsigned long long segmentKey = SetMeshEdgeKey(segmentCornerIndices[0], segmentCornerIndices[1]);
			BoundarySegmentInfo segmentInfo;
			segmentInfo.startTime = edgeIntersectionsInfo[i].time;
			segmentInfo.endTime = edgeIntersectionsInfo[i + 1].time;
			segmentInfo.halfEdge = halfEdgeIndex;

			if (localBoundarySegmentsInfo.find(segmentKey) != localBoundarySegmentsInfo.end()) {
				printf("Error. Replicated segment key!\n");
				return 0;
			}
			localBoundarySegmentsInfo[segmentKey] = segmentInfo;
		}
		boundaryEdgeIntersections[halfEdgeIndex] = edgeIntersectionsInfo;
	}
	return 1;
}

void IndexedPolygonFromCell(const int i, const int j, const GridChart & gridChart, IndexedIntersectionPolygon & polygon) {

	int width = gridChart.width;
	int height = gridChart.height;
	int gridIndexOffset = gridChart.gridIndexOffset;

	polygon.vertices.resize(4);
	polygon.indices.resize(4);
	polygon.edgeIndices.resize(4);

	polygon.vertices[0] = Point2D<double>(double(i)*gridChart.cellSizeW, double(j)*gridChart.cellSizeH);
	polygon.vertices[1] = Point2D<double>(double(i + 1)*gridChart.cellSizeW, double(j)*gridChart.cellSizeH);
	polygon.vertices[2] = Point2D<double>(double(i + 1)*gridChart.cellSizeW, double(j + 1)*gridChart.cellSizeH);
	polygon.vertices[3] = Point2D<double>(double(i)*gridChart.cellSizeW, double(j + 1)*gridChart.cellSizeH);

	polygon.indices[0] = SetIntersectionKey(4294967295, gridIndexOffset + i + j*width);
	polygon.indices[1] = SetIntersectionKey(4294967295, gridIndexOffset + i + 1 + j*width);
	polygon.indices[2] = SetIntersectionKey(4294967295, gridIndexOffset + i + 1 + (j + 1)*width);
	polygon.indices[3] = SetIntersectionKey(4294967295, gridIndexOffset + i + (j + 1)*width);

	polygon.edgeIndices[0] = gridIndexOffset + width*height + i + j* width;
	polygon.edgeIndices[1] = gridIndexOffset + i + 1 + j* width;
	polygon.edgeIndices[2] = gridIndexOffset + width*height + i + (j + 1)* width;
	polygon.edgeIndices[3] = gridIndexOffset + i + j* width;
}


int InitializeChartBoundaryPolygons(const std::vector<int> & atlasEdgeIndex, const std::vector<int> & oppositeHalfEdge, const AtlasChart & atlasChart, GridChart & gridChart,
	const int numInteriorNodes, const int numBoundaryVertices, const int numBoundaryNodes,
	std::unordered_map<int, std::vector<IntersectionInfo>> boundaryEdgeIntersections,
	std::unordered_map<unsigned long long, BoundarySegmentInfo> & boundarySegmentsInfo,
	std::unordered_map<unsigned long long, int> & boundaryNodeIndices,
	std::unordered_map<unsigned long long, Point2D<double>> & boundaryNodePosition,
	std::vector<int> & coveredOppositeBoundaryNode) {

	std::unordered_map<int, std::vector<std::pair<unsigned long long, unsigned long long>>> cellSegments;

	//(1) Generate cell segments
	for (int t = 0; t < atlasChart.triangles.size(); t++) {
		Point2D< double > tPos[3];
		for (int i = 0; i < 3; i++) tPos[i] = atlasChart.vertices[atlasChart.triangles[t][i]] - gridChart.corner;
		int minCorner[2];
		int maxCorner[2];
		GetTriangleIntegerBBox(tPos, 1.0 / gridChart.cellSizeW, 1.0 / gridChart.cellSizeH , minCorner, maxCorner);

		IndexedIntersectionTriangle indexedTriangle;
		for (int k = 0; k < 3; k++) {
			indexedTriangle.vertices[k] = tPos[k];
			indexedTriangle.indices[k] = SetIntersectionKey(atlasChart.triangles[t][k], 4294967295);
			indexedTriangle.edgeIndices[k] = atlasChart.atlasEdgeIndices[3 * t + k];
		}

		for (int j = minCorner[1]; j < maxCorner[1]; j++) for (int i = minCorner[0]; i < maxCorner[0]; i++) {
			if (gridChart.cellType(i, j) == 0) {
				int cellId = gridChart.localBoundaryCellIndex(i, j);
				if (cellId == -1) {
					printf("Error: Boundary cell invalid id! \n");
					return 0;
				}
				IndexedIntersectionPolygon cellPolygon;
				IndexedPolygonFromCell(i, j, gridChart, cellPolygon);
				if (ClipIndexedIntersectionPolygonToIndexedIntersectionTriangle(cellPolygon, indexedTriangle)) {
					for (int s = 0; s < cellPolygon.indices.size(); s++){
						std::pair<unsigned long long, unsigned long long> edge(cellPolygon.indices[s], cellPolygon.indices[(s + 1) % cellPolygon.indices.size()]);
						cellSegments[cellId].push_back(edge);
					}
				}
			}
		}
	}

	//(2) Process cells

	gridChart.boundaryPolygons.resize(gridChart.numBoundaryCells);
	for (auto cellIter = cellSegments.begin(); cellIter != cellSegments.end(); cellIter++) {
		std::vector<std::pair<unsigned long long, unsigned long long>> segments = (*cellIter).second;
		int cellId = (*cellIter).first;

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

		std::vector<std::vector<unsigned long long>> loopKeys;
		if (!LoopVertices(forwardMap, loopKeys)) {
			printf("Cell segments does not form closed loops! \n");
			for (int i = 0; i < gridChart.width; i++)for (int j = 0; j < gridChart.height; j++) {
				if (gridChart.localBoundaryCellIndex(i, j) == cellId) {
					printf("Cell %d %d \n", gridChart.cornerCoords[0] + i, gridChart.cornerCoords[1] + j);
				}
			}
			return 0;
		}

		//Keep only the grid nodes and the boundary vertices
		std::vector<std::vector<Point2D<double>>> loopPositions(loopKeys.size());
		std::vector<std::vector<int>> loopIndices(loopKeys.size());
		std::vector<std::vector<int>> loopAtlasVertexIndices(loopKeys.size());
		for (int i = 0; i < loopKeys.size(); i++) {
			std::vector<unsigned long long> reduceLoop;
			std::vector<Point2D<double>> currentLoopPositions;
			std::vector<int> currentLoopIndices;
			std::vector<int> currentAtlasVertexIndices;
			for (int j = 0; j < loopKeys[i].size(); j++) {
				unsigned long long currentVertexKey = loopKeys[i][j];
				unsigned long tElement, vElement;
				GetIntersectionKey(currentVertexKey, tElement, vElement);
				if (tElement == 4294967295) {
					vElement -= gridChart.gridIndexOffset;
					int pi = vElement % gridChart.width;
					int pj = vElement / gridChart.width;
					int interiorTexelIndex = gridChart.globalInteriorTexelIndex(pi, pj);
					if (interiorTexelIndex == -1) {
						printf("Invalid texel! \n");
						return 0;
					}
					reduceLoop.push_back(currentVertexKey);
					currentLoopPositions.push_back(Point2D<double>(double(pi)*gridChart.cellSizeW, double(pj)*gridChart.cellSizeH));
					currentLoopIndices.push_back(interiorTexelIndex);
					currentAtlasVertexIndices.push_back(-1);
				}
				else if (boundaryNodeIndices.find(currentVertexKey) != boundaryNodeIndices.end()) {
					reduceLoop.push_back(currentVertexKey);
					currentLoopPositions.push_back(boundaryNodePosition[currentVertexKey]);
					currentLoopIndices.push_back(boundaryNodeIndices[currentVertexKey] + numInteriorNodes);

					if (vElement == 4294967295) currentAtlasVertexIndices.push_back(tElement);
					else currentAtlasVertexIndices.push_back(-1);
				}
			}
			if (reduceLoop.size() == 0) {
				printf("Reduced loop cannot be empty!. Original loop size %d \n", (int)loopKeys[i].size());
				for (int j = 0; j < loopKeys[i].size(); j++) {
					unsigned long long currentVertexKey = loopKeys[i][j];
					unsigned long tElement, vElement;
					GetIntersectionKey(currentVertexKey, tElement, vElement);
					printf("Mesh Index %lu. Grid Index %lu \n", tElement, vElement);
				}
				return 0;
			}
			loopKeys[i] = reduceLoop;
			loopPositions[i] = currentLoopPositions;
			loopIndices[i] = currentLoopIndices;
			loopAtlasVertexIndices[i] = currentAtlasVertexIndices;
		}

		//Add the intermediate vertices for boundary segments
		std::vector<std::vector<int>> loopAtlasEdges(loopKeys.size());
		std::vector<std::vector<int>> loopAtlasVertexParentEdges(loopKeys.size());
		for (int i = 0; i < loopKeys.size(); i++) {

			std::vector<std::vector<int>> indicesToInsert;
			std::vector<std::vector<double>> timesToInsert;
			std::vector<bool> isInsertionPos(loopKeys[i].size(), false);
			std::vector<int> polygonAtlaEdgeIndex(loopKeys[i].size(), -1);
			std::vector<int> polygonAtlaVertexParentEdges(loopKeys[i].size(), -1);

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
						BoundarySegmentInfo segmentInfo = boundarySegmentsInfo[segmentKey];
						polygonAtlaEdgeIndex[j] = atlasEdgeIndex[segmentInfo.halfEdge];

						int oppHalfEdge = oppositeHalfEdge[segmentInfo.halfEdge];
						if (oppHalfEdge != -1) {
							double startTime = segmentInfo.startTime;
							double endTime = segmentInfo.endTime;


							if (boundaryEdgeIntersections.find(oppHalfEdge) == boundaryEdgeIntersections.end()) {
								printf(" Opposite edge intersections not found! \n");
								printf(" Current  edge %d. Opposite %d \n", segmentInfo.halfEdge, oppHalfEdge);
								return 0;
							}

							std::vector<int> segmentIndicesToInsert;
							std::vector<double> segmentTimesToInsert;

							std::vector<IntersectionInfo> oppositeIntersection = boundaryEdgeIntersections[oppHalfEdge];
							for (int oppIter = (int)oppositeIntersection.size() - 2; oppIter > 0; oppIter--) {
								double reversedTime = 1.0 - oppositeIntersection[oppIter].time;
								double normalizedTime = (reversedTime - startTime) / (endTime - startTime);
								if (normalizedTime > 0 && normalizedTime < 1) {
									segmentIndicesToInsert.push_back(oppositeIntersection[oppIter].intersectionIndex);
									segmentTimesToInsert.push_back(normalizedTime);
								}
							}


							if (segmentIndicesToInsert.size() > 0) {
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
			std::vector<Point2D<double>> expandedLoopPositions;
			std::vector<int> expandedLoopIndices;
			std::vector<int> expandedLoopAtlasVertexIndices;
			std::vector<int> expandedLoopAtlasEdgeIndex;
			std::vector<int> expandedVertexParentEdgeIndex;
			for (int j = 0; j < loopKeys[i].size(); j++) {

				//expandedLoopKeys.push_back(loopKeys[i][j]);
				expandedLoopPositions.push_back(loopPositions[i][j]);
				expandedLoopIndices.push_back(loopIndices[i][j]);
				expandedLoopAtlasVertexIndices.push_back(loopAtlasVertexIndices[i][j]);

				int currentSegmentAtlasEdgeIndex = polygonAtlaEdgeIndex[j];
				expandedLoopAtlasEdgeIndex.push_back(currentSegmentAtlasEdgeIndex);
				expandedVertexParentEdgeIndex.push_back(polygonAtlaVertexParentEdges[j]);

				if (isInsertionPos[j]) {
					std::vector<int> _indices = indicesToInsert[insertionCount];
					std::vector<double> _times = timesToInsert[insertionCount];
					Point2D<double> currentPos = loopPositions[i][j];
					Point2D<double> nextPos = loopPositions[i][(j + 1) % loopKeys[i].size()];

					if (currentSegmentAtlasEdgeIndex == -1) {
						printf("Invalid atlas edge index! \n");
						return 0;
					}
					for (int k = 0; k < _indices.size(); k++) {
						double alpha = _times[k];
						Point2D<double> interpolatedPos = (1.0 - alpha)*currentPos + nextPos*alpha;
						if(1){ //Add orthogonal perturbation to avoid colinearity. Delanauy triangulation library failed for colinear inputs on some systems.
							Point2D<double> segmentDir = nextPos - currentPos;
							interpolatedPos += Point2D<double>(-segmentDir[1],segmentDir[0]) * ((2.0 * double(rand())/double(RAND_MAX)) - 1.0) * 1e-10;
						}
						if (_indices[k] < numBoundaryVertices || _indices[k] > numBoundaryNodes) {
							printf("Out of bounds index: %d not in %d %d \n", _indices[k], numBoundaryVertices, numBoundaryNodes);
							return 0;
						}
						coveredOppositeBoundaryNode[_indices[k] - numBoundaryVertices] += 1;

						AuxiliarNode auxNode;
						auxNode.index = _indices[k] + numInteriorNodes;
						auxNode.position = interpolatedPos;
						gridChart.auxiliarNodes.push_back(auxNode);

						expandedLoopPositions.push_back(interpolatedPos);
						expandedLoopIndices.push_back(_indices[k] + numInteriorNodes);
						expandedLoopAtlasVertexIndices.push_back(-1);
						expandedLoopAtlasEdgeIndex.push_back(currentSegmentAtlasEdgeIndex);
						expandedVertexParentEdgeIndex.push_back(currentSegmentAtlasEdgeIndex);
					}
					insertionCount++;
				}
			}
			if (insertionCount != indicesToInsert.size()) {
				printf("Intersection chains count does not match! \n");
				return 0;
			}
			loopPositions[i] = expandedLoopPositions;
			loopIndices[i] = expandedLoopIndices;
			loopAtlasVertexIndices[i] = expandedLoopAtlasVertexIndices;
			loopAtlasEdges[i] = expandedLoopAtlasEdgeIndex;
			loopAtlasVertexParentEdges[i] = expandedVertexParentEdgeIndex;
		}

		for (int i = 0; i < loopPositions.size(); i++) {
			AtlasIndexedPolygon poly;
			poly.vertices = loopPositions[i];
			poly.indices = loopIndices[i];
			poly.atlasVertexIndices = loopAtlasVertexIndices[i];
			poly.atlasEdgeIndices = loopAtlasEdges[i];
			poly.atlasVertexParentEdge = loopAtlasVertexParentEdges[i];
			gridChart.boundaryPolygons[cellId].push_back(poly);
		}
	}
	return 1;
}



int InitializeBoundaryPolygons(const std::vector<int> & atlasEdgeIndex, const std::vector<int> & oppositeHalfEdge, const std::vector<AtlasChart> & atlasCharts, std::vector<GridChart> & gridCharts, std::unordered_map<int, int> & boundaryVerticesIndices,
	const int numInteriorNodes, const int numBoundaryVertices, int & numBoundaryNodes, const bool isClosedMesh) { //Fine System

	int lastBoundaryIndex = numBoundaryVertices;
	std::unordered_map<int, std::vector<IntersectionInfo>> boundaryEdgeIntersections;

	std::vector<std::unordered_map<unsigned long long, BoundarySegmentInfo>> localBoundarySegmentsInfo(gridCharts.size());
	std::vector<std::unordered_map<unsigned long long, int>> localBoundaryNodeIndices(gridCharts.size());
	std::vector<std::unordered_map<unsigned long long, Point2D<double>>>  localBoundaryNodePosition(gridCharts.size());

	for (int i = 0; i < gridCharts.size(); i++) {
		if (!InitializeChartBoundaryEdgeGridIntersections(atlasCharts[i], gridCharts[i], boundaryVerticesIndices, lastBoundaryIndex, numInteriorNodes,
			boundaryEdgeIntersections, localBoundarySegmentsInfo[i], localBoundaryNodeIndices[i], localBoundaryNodePosition[i])) {
			printf("ERROR: Unable to initialize chart boundary edge grid intersections! \n");
		}
	}

	numBoundaryNodes = lastBoundaryIndex;
	std::vector<int> coveredOppositeBoundaryNode(numBoundaryNodes - numBoundaryVertices, 0);
	for (int i = 0; i < gridCharts.size(); i++) {
		if (!InitializeChartBoundaryPolygons(atlasEdgeIndex, oppositeHalfEdge, atlasCharts[i], gridCharts[i],
			numInteriorNodes, numBoundaryVertices, numBoundaryNodes,
			boundaryEdgeIntersections, localBoundarySegmentsInfo[i], localBoundaryNodeIndices[i], localBoundaryNodePosition[i], coveredOppositeBoundaryNode)) {
			printf("ERROR: Unable to initialize boundary polygons! \n");
			return 0;
		}
	}

	if (isClosedMesh) {
		for (int i = 0; i < coveredOppositeBoundaryNode.size(); i++) {
			if (coveredOppositeBoundaryNode[i] != 1) {
				printf("ERROR: non opposite boundary node at node %d !",i);
				return 0;
			}
		}
	}

	return 1;
}

#include "ConstrainedTriangulation.h"

int InitializeChartQuadraticElements(GridChart & gridChart, std::unordered_map<unsigned long long, int> & midPointIndex, int & lastMidPointIndex, const int & previouslyAddedNodes) {

	const std::vector<std::vector<AtlasIndexedPolygon>> & boundaryPolygons = gridChart.boundaryPolygons;
	std::vector<std::vector<BoundaryIndexedTriangle>> & boundaryTriangles = gridChart.boundaryTriangles;
	boundaryTriangles.resize(gridChart.numBoundaryCells);
	int numBoundaryTriangles = 0;
	for (int i = 0; i < boundaryPolygons.size(); i++) {
		for (int j = 0; j < boundaryPolygons[i].size(); j++) {
			const AtlasIndexedPolygon & currentPolygon = boundaryPolygons[i][j];
			std::vector<TriangleIndex> delanauyTriangles;
			TriangulatePolygon(currentPolygon.vertices, delanauyTriangles);
			if (delanauyTriangles.size() != currentPolygon.vertices.size() - 2) {
				int localCellPos[2] = { -1,-1 };
				for (int li = 0; li < gridChart.localBoundaryCellIndex.width(); li++)for (int lj = 0; lj < gridChart.localBoundaryCellIndex.height(); lj++) if (gridChart.localBoundaryCellIndex(li, lj) == i) localCellPos[0] = li, localCellPos[1] = lj;
				printf("WARNING: Unexpected number of triangles produced by delaunay triangulation at global cell (%d,%d). Polygon may self intersect!.\n", gridChart.cornerCoords[0] + localCellPos[0], gridChart.cornerCoords[1] + localCellPos[1]);
				printf("Num Vertices %d \n", (int)currentPolygon.vertices.size());
				printf("Num Face %d \n", (int)delanauyTriangles.size());
				for (int iv = 0; iv < currentPolygon.vertices.size(); iv++)printf("%f %f %f \n", currentPolygon.vertices[iv][0], currentPolygon.vertices[iv][1], 0.);
				for (int it = 0; it < delanauyTriangles.size(); it++)printf("3 %d %d %d \n", delanauyTriangles[it][0], delanauyTriangles[it][1], delanauyTriangles[it][2]);
				//return 0;
			}

			for (int k = 0; k < delanauyTriangles.size(); k++) {
				BoundaryIndexedTriangle triangle;
				for (int v = 0; v < 3; v++) {
					int currentPolygonVertex = delanauyTriangles[k][v];
					triangle.indices[v] = currentPolygon.indices[currentPolygonVertex];
					triangle.vertices[v] = currentPolygon.vertices[currentPolygonVertex];
					triangle.atlasVertexIndices[v] = currentPolygon.atlasVertexIndices[currentPolygonVertex];
					triangle.atlasVertexParentEdge[v] = currentPolygon.atlasVertexParentEdge[currentPolygonVertex];
					triangle.atlasEdgeIndices[v] = -1;
					bool isAtlasEdge = false;
					int nextPolygonVertex = (currentPolygonVertex + 1) % currentPolygon.vertices.size();
					if (delanauyTriangles[k][(v + 1) % 3] == nextPolygonVertex) { //preserve orientation
						triangle.atlasEdgeIndices[v] = currentPolygon.atlasEdgeIndices[currentPolygonVertex];
					}
					int previousPolygonVertex = (int)( (currentPolygonVertex + currentPolygon.vertices.size() - 1) % currentPolygon.vertices.size() );
					if (delanauyTriangles[k][(v + 1) % 3] == previousPolygonVertex) {//reverse orientation
						triangle.atlasEdgeIndices[v] = currentPolygon.atlasEdgeIndices[previousPolygonVertex];
					}

				}
				int midVexterIndex[3];

				for (int v = 0; v < 3; v++) {
					int mIndex = -1;

					Point2D<double> edgeMidPoint = (triangle.vertices[(v + 1) % 3] + triangle.vertices[(v + 2) % 3]) / 2.0;
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

					AuxiliarNode auxNode;
					auxNode.index = midVexterIndex[v];
					auxNode.position = edgeMidPoint;
					gridChart.auxiliarNodes.push_back(auxNode);

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
	return 1;
}

int InitializeBoundaryQuadraticElements(std::vector<GridChart> & gridCharts, const int previouslyAddedNodes, int & numMidPoints) {
	int lastMidPointIndex = 0;
	std::unordered_map<unsigned long long, int> midPointIndex;
	for (int i = 0; i < gridCharts.size(); i++) {
		if (!InitializeChartQuadraticElements(gridCharts[i], midPointIndex, lastMidPointIndex, previouslyAddedNodes)) {
			printf("Unable to initialize quadratic elements!\n");
			return 0;
		}
	}
	numMidPoints = lastMidPointIndex;
	return 1;
}


int InitializeBoundaryTriangulation(GridAtlas & gridAtlas, AtlasMesh & atlasMesh, std::vector<AtlasChart> & atlasCharts, std::vector<int> & oppositeHalfEdge, std::unordered_map<int, int> & boundaryVerticesIndices, const int & numBoundaryVertices, const bool & isClosedMesh, bool verbose = false) {
	if (!InitializeBoundaryPolygons(atlasMesh.halfEdgeToEdgeIndex, oppositeHalfEdge, atlasCharts, gridAtlas.gridCharts, boundaryVerticesIndices, gridAtlas.numInteriorTexels, numBoundaryVertices, gridAtlas.numBoundaryNodes, isClosedMesh)) {
		printf("ERROR: Unable to initialize boundary polygons! \n");
		return 0;
	}

	if (!InitializeBoundaryQuadraticElements(gridAtlas.gridCharts, gridAtlas.numBoundaryNodes + gridAtlas.numInteriorTexels, gridAtlas.numMidPoints)) {
		printf("ERROR: Unable to initialize boundary quadratic elements! \n");
		return 0;
	}

	if (verbose) {
		printf("\t Boundary vertices %d \n", numBoundaryVertices);
		printf("\t Boundary nodes %d \n", gridAtlas.numBoundaryNodes);
		printf("\t Boundary cells %d \n", gridAtlas.numBoundaryCells);
		printf("\t Mid points %d \n", gridAtlas.numMidPoints);
	}

	gridAtlas.numFineNodes = gridAtlas.numInteriorTexels + gridAtlas.numBoundaryNodes + gridAtlas.numMidPoints;

	if (verbose) printf("Num fine nodes %d =  Interior texels %d + Auxiliar boundary nodes %d \n", gridAtlas.numFineNodes, gridAtlas.numInteriorTexels, gridAtlas.numBoundaryNodes + gridAtlas.numMidPoints);

	return 1;
}

class BoundaryProlongationData {
public:
	std::vector<int> fineBoundaryIndex;
	int numFineBoundarNodes;
	SparseMatrix<double, int> coarseBoundaryFineBoundaryProlongation;
	SparseMatrix<double, int> fineBoundaryCoarseBoundaryRestriction;
};

int InitializeCoarseBoundaryToFineBoundaryProlongation(const GridAtlas & gridAtlas, SparseMatrix<double, int> & coarseBoundaryFineBoundaryProlongation, std::vector<int> & fineBoundaryIndex, int & numFineBoundaryNodes, bool verbose = false) {

	std::vector<int> boundaryFineToFullFine;
	std::vector<Eigen::Triplet<double>> prolongationTriplets;
	const std::vector<int> & boundaryAndDeepIndex = gridAtlas.boundaryAndDeepIndex;


	fineBoundaryIndex.resize(gridAtlas.numFineNodes, -1);
	int lastFineBoundaryIndex = 0;
	for (int i = 0; i < gridAtlas.gridCharts.size(); i++) {
		const GridChart & gridChart = gridAtlas.gridCharts[i];
		for (int j = 0; j < gridChart.globalInteriorTexelIndex.size(); j++) {
			if (gridChart.globalInteriorTexelIndex[j] != -1 && gridChart.globalTexelDeepIndex[j] == -1) { //Interior but not deep
				int coarseGlobalIndex = gridChart.globalTexelIndex[j];

				int coarseBondaryIndex = boundaryAndDeepIndex[coarseGlobalIndex] - 1;

				if (coarseBondaryIndex >= 0) {
					prolongationTriplets.push_back(Eigen::Triplet<double>(lastFineBoundaryIndex, coarseBondaryIndex, 1.0));
				}
				else {
					printf("ERROR: Coarse node is not boundary! \n");
					printf("Global index %d. Boundary index %d \n", coarseGlobalIndex, coarseBondaryIndex);
					return 0;
				}
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

	if (verbose)printf("Fine boundary elements %d \n", lastFineBoundaryIndex);
	numFineBoundaryNodes = lastFineBoundaryIndex;

	const std::vector<GridChart> & gridCharts = gridAtlas.gridCharts;
	const std::vector<GridNodeInfo> & nodeInfo = gridAtlas.nodeInfo;

	int numInteriorTexels = gridAtlas.numInteriorTexels;
	int numFineNodes = gridAtlas.numFineNodes;
	int numCoarseNodes = gridAtlas.numTexels;

	const int numAuxiliarNodes = numFineNodes - numInteriorTexels;
	std::vector<int> auxiliarNodesDegree(numAuxiliarNodes, 0);

	for (int i = 0; i < gridCharts.size(); i++) {
		const GridChart & gridChart = gridCharts[i];
		for (int j = 0; j < gridChart.auxiliarNodes.size(); j++) {
			int auxiliarId = gridChart.auxiliarNodes[j].index - numInteriorTexels;
			auxiliarNodesDegree[auxiliarId]++;
		}
	}

	double precision_error = 1e-10;

	std::vector<double> auxiliarNodesCumWeight(numAuxiliarNodes, 0);

	for (int i = 0; i < gridCharts.size(); i++) {
		const GridChart & gridChart = gridCharts[i];
		for (int j = 0; j < gridChart.auxiliarNodes.size(); j++) {
			int auxiliarId = gridChart.auxiliarNodes[j].index - numInteriorTexels;
			int fineBoundaryId = fineBoundaryIndex[gridChart.auxiliarNodes[j].index];
			int nodeDegree = auxiliarNodesDegree[auxiliarId];
			Point2D<double>nodePosition = gridChart.auxiliarNodes[j].position;
			int corner[2] = { (int)floor(nodePosition[0] / gridChart.cellSizeW), (int)floor(nodePosition[1] / gridChart.cellSizeH) };
			int cellId = gridChart.localCellIndex(corner[0], corner[1]);
			if (cellId == -1) {
				printf("ERROR : Invalid cell index!\n");
				printf("Node position %f %f \n", nodePosition[0] / gridChart.cellSizeW, nodePosition[1] / gridChart.cellSizeH);
				return 0;
			}
			nodePosition[0] /= gridChart.cellSizeW;
			nodePosition[1] /= gridChart.cellSizeH;
			nodePosition[0] -= double(corner[0]);
			nodePosition[1] -= double(corner[1]);
			if (nodePosition[0] < 0 - precision_error || nodePosition[0] > 1 + precision_error || nodePosition[1] < 0 - precision_error || nodePosition[1] > 1 + precision_error) {
				printf("Sample out of unit box! (%f %f)\n", nodePosition[0], nodePosition[1]);
				return 0;
			}
			for (int k = 0; k < 4; k++) {
				double texelWeight = BilinearElementValue(k, nodePosition) / double(nodeDegree);
				if (fabs(texelWeight) > 1e-11) {
					auxiliarNodesCumWeight[auxiliarId] += texelWeight;
					int texelIndex = gridChart.bilinearElementIndices[cellId][k];
					if (nodeInfo[texelIndex].nodeType == 2) {
						printf("ERROR: Deep texel cannot be in the support of an auxiliar node. Weight %g (A)\n", texelWeight);
						for (int _k = 0; _k < 4; _k++) printf("Neighbours weight %g \n", BilinearElementValue(_k, nodePosition) / double(nodeDegree));
						return 0;
					}
					if (gridChart.auxiliarNodes[j].index < numInteriorTexels || gridChart.auxiliarNodes[j].index > numFineNodes || texelIndex < 0 || texelIndex> numCoarseNodes) {
						printf("Out of bound index! \n");
						return 0;
					}

					int coarseBondaryId = boundaryAndDeepIndex[texelIndex] - 1;

					if (coarseBondaryId < 0) {
						printf("ERROR: Coarse node is not boundary! \n");
						return 0;
					}
					prolongationTriplets.push_back(Eigen::Triplet<double>(fineBoundaryId, coarseBondaryId, texelWeight));
				}
			}
		}
	}

	int numBoundaryTexels = gridAtlas.numBoundaryTexels;
	coarseBoundaryFineBoundaryProlongation = SetSparseMatrix( prolongationTriplets , numFineBoundaryNodes , numBoundaryTexels , false );
	return 1;
}


int InitializeBoundaryProlongationData(const GridAtlas & gridAtlas, BoundaryProlongationData & boundaryProlongation) {
	if (!InitializeCoarseBoundaryToFineBoundaryProlongation(gridAtlas, boundaryProlongation.coarseBoundaryFineBoundaryProlongation, boundaryProlongation.fineBoundaryIndex, boundaryProlongation.numFineBoundarNodes)) {
		printf("ERROR: Unable to initialize coarse to fine boundary prolongation\n");
		return 0;
	}
	boundaryProlongation.fineBoundaryCoarseBoundaryRestriction = boundaryProlongation.coarseBoundaryFineBoundaryProlongation.transpose();
	return 1;
}