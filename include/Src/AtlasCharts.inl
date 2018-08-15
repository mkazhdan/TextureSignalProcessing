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

int InitializeAtlasCharts(AtlasMesh & atlasMesh, const std::vector<bool> & isBoundaryHalfEdge, const int width, const int height, std::vector<AtlasChart> & atlasCharts) {
	atlasCharts.resize(atlasMesh.numCharts);

	for (int i = 0; i < atlasMesh.numCharts; i++) {
		atlasCharts[i].minCorner = Point2D<double>(1, 1);
		atlasCharts[i].maxCorner = Point2D<double>(0, 0);
	}

	std::vector<int> lastVertexId(atlasMesh.numCharts, 0);
	std::vector<int> chartVertexId(atlasMesh.vertices.size(), -1);

	atlasMesh.triangleIndexInChart.resize(atlasMesh.triangles.size());
	std::vector<int> lastTriangleId(atlasMesh.numCharts, 0);

	for (int t = 0; t < atlasMesh.triangles.size(); t++) {
		int chartId = atlasMesh.triangleChartIndex[t];
		//atlasCharts[chartId].atlasTriangleIndices.push_back(t);

		atlasMesh.triangleIndexInChart[t] = lastTriangleId[chartId];
		lastTriangleId[chartId]++;

		atlasCharts[chartId].meshTriangleIndices.push_back(t);
		//atlasCharts[chartId].metric.push_back(atlasMesh.metric[t]);
		int vertexId[3];
		for (int k = 0; k < 3; k++) {
			atlasCharts[chartId].atlasEdgeIndices.push_back(atlasMesh.halfEdgeToEdgeIndex[3 * t + k]);
			int atlasVertexId = atlasMesh.triangles[t][k];
			if (chartVertexId[atlasVertexId] == -1) {
				chartVertexId[atlasVertexId] = lastVertexId[chartId];
				lastVertexId[chartId]++;

				Point2D<double> vertexPos = atlasMesh.vertices[atlasVertexId];
				for (int c = 0; c < 2; c++) {
					atlasCharts[chartId].minCorner[c] = std::min<double>(vertexPos[c], atlasCharts[chartId].minCorner[c]);
					atlasCharts[chartId].maxCorner[c] = std::max<double>(vertexPos[c], atlasCharts[chartId].maxCorner[c]);
				}

				atlasCharts[chartId].vertices.push_back(vertexPos);
				//atlasCharts[chartId].atlasVertexIndices.push_back(atlasVertexId);
				atlasCharts[chartId].meshVertexIndices.push_back(atlasMesh.vertexMap[atlasVertexId]);
			}
			vertexId[k] = chartVertexId[atlasVertexId];
		}
		atlasCharts[chartId].triangles.push_back(TriangleIndex(vertexId[0], vertexId[1], vertexId[2]));
	}

	for (int i = 0; i < atlasCharts.size(); i++) {
		Point2D<double> midBBox = (atlasCharts[i].minCorner + atlasCharts[i].maxCorner) / 2.0;
		midBBox[0] *= double(width);
		midBBox[1] *= double(height);
		midBBox -= Point2D<double>(0.5, 0.5);
		midBBox = Point2D<double>(floor(midBBox[0]), floor(midBBox[1]));
		atlasCharts[i].originCoords[0] = int(round(midBBox[0]));
		atlasCharts[i].originCoords[1] = int(round(midBBox[1]));

		midBBox += Point2D<double>(0.5, 0.5);
		midBBox[0] /= double(width);
		midBBox[1] /= double(height);
		atlasCharts[i].gridOrigin = midBBox;
	}
	for (int i = 0; i < atlasCharts.size(); i++) {
		std::vector<int> & boundaryHalfEdges = atlasCharts[i].boundaryHalfEdges;
		for (int t = 0; t < atlasCharts[i].meshTriangleIndices.size(); t++) {
			int tIndex = atlasCharts[i].meshTriangleIndices[t];
			for (int k = 0; k < 3; k++) {
				if (isBoundaryHalfEdge[3 * tIndex + k]) {
					boundaryHalfEdges.push_back(3 * t + k);
				}
			}
		}
	}

	return 1;
}



int InitializeAtlasMesh(const TexturedMesh & mesh, const int width, const int height, AtlasMesh & atlasMesh, std::vector<AtlasChart> & atlasCharts, std::vector<int> & oppositeHalfEdge, std::unordered_map<int, int> & boundaryVerticesIndices, int & numBoundaryVertices, bool & isClosedMesh , bool verbose ) {
	if( !InitializeAtlasMesh( mesh , atlasMesh , width , height , verbose ) )
	{
		printf("Failed atlas mesh construction! \n");
		return 0;
	}

	std::vector<int> boundaryHalfEdges;
	std::vector<bool> isBoundaryHalfEdge;
	if (!InitializeBoundaryHalfEdges(mesh, boundaryHalfEdges, oppositeHalfEdge, isBoundaryHalfEdge, isClosedMesh)) {
		printf("Failed boundary halfedge construction! \n");
		return 0;
	}

	int lastBoundaryIndex;
	if (!InitiallizeBoundaryVertices(mesh, boundaryHalfEdges, boundaryVerticesIndices, lastBoundaryIndex)) {
		printf("Failed intialize boundary vertices! \n");
		return 0;
	}

	numBoundaryVertices = lastBoundaryIndex;

	if (!InitializeAtlasCharts(atlasMesh, isBoundaryHalfEdge, width, height, atlasCharts)) {
		printf("Failed intialize atlas charts! \n");
		return 0;
	}

	return 1;
}
