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

#ifndef ATLAS_MESH_INLCUDED
#define ATLAS_MESH_INLCUDED

#include "SimpleMesh.h"
#include "ChartDecomposition.h"
#include <set>

class AtlasMesh{
public:
	std::vector<Point2D<double>> vertices;
	std::vector<TriangleIndex> triangles;
	std::vector<int> triangleIndexInChart;
	std::vector<int> triangleChartIndex;
	std::vector<int> halfEdgeToEdgeIndex;
	std::vector<int> vertexMap;
	//std::vector<SquareMatrix<double, 2>> metric;
	int numCharts;
};

class AtlasChart {
public:
	Point2D<double> minCorner;
	Point2D<double> maxCorner;
	Point2D<double> gridOrigin;
	int originCoords[2];
	std::vector<TriangleIndex> triangles;
	std::vector<Point2D<double>> vertices;
	std::vector<int> boundaryHalfEdges;
	std::vector<int> atlasEdgeIndices;

	std::vector<int> meshVertexIndices;
	std::vector<int> meshTriangleIndices;
	//std::vector<SquareMatrix<double, 2>> metric;
};

class IndexedVector2D{
public:
	IndexedVector2D(Point2D<double> p_p, int p_index, int p_vertex){
		p = p_p;
		index = p_index;
		vertex = p_vertex;
	}
	Point2D<double> p;
	int index;
	int vertex;
};

class IndexedVector2DComparison
{
public:
	bool operator()(const IndexedVector2D & p1, const IndexedVector2D & p2) const{
		for (int i = 0; i < 2; i++){
			if (p1.p[i] < p2.p[i])
				return true;
			else if (p2.p[i] < p1.p[i])
				return false;
			else {
				if (p1.vertex < p2.vertex) return true;
				else if (p2.vertex < p1.vertex) return false;
			}
		}
		return false;
	}
};

#include "AtlasMesh.inl"
#include "AtlasCharts.inl"
#endif// ATLAS_MESH_INLCUDED