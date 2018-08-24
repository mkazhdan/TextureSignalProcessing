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
#ifndef SIMPLE_MESH_INCLUDED
#define SIMPLE_MESH_INCLUDED

#include <Eigen/Sparse>
#include <Misha/Ply.h>
#include <Misha/Image.h>
#include <Src/VectorIO.h>

class MeshSample{
public:
	MeshSample(){
		tId = -1;
		barycentricCoords = Point2D<double>(0, 0);
	}
	MeshSample(int _tId, Point2D<double> _barycentricCoords){
		tId = _tId;
		barycentricCoords = _barycentricCoords;
	}
	int tId;
	Point2D<double> barycentricCoords;
};

unsigned long long SetMeshEdgeKey(const unsigned long i0, const unsigned long i1){
	return ( ( (static_cast<unsigned long long>(i0) << 32) & 0xFFFFFFFF00000000) | (static_cast<unsigned long long>(i1) & 0x00000000FFFFFFFF) );
}

void GetMeshEdgeIndices(unsigned long long key, unsigned long & i0, unsigned long & i1){
	i1 = static_cast<unsigned long>(key & 0x00000000FFFFFFFF);
	i0 = static_cast<unsigned long>((key >> 32) & 0x00000000FFFFFFFF);
}

class SimpleMesh{
public:
	std::vector<Point3D<double>> vertices;
	std::vector<Point3D<double>> normals;
	std::vector<TriangleIndex> triangles;
};

class ColoredMesh : public SimpleMesh{
public:
	std::vector<Point3D<double>> colors;
};


class TexturedMesh : public SimpleMesh
{
public:
	std::vector< Point2D< double > > textureCoordinates;
	Image< Point3D< float > > texture;
};

void UpdateNormals(SimpleMesh & mesh){
	mesh.normals.clear();
	mesh.normals.resize(mesh.vertices.size(), Point3D<double>(0.0, 0.0, 0.0));
	for (int t = 0; t < mesh.triangles.size(); t++){
		Point3D<double> d01 = mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]];
		Point3D<double> d02 = mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]];
		Point3D<double> n = Point3D<double>::CrossProduct(d01,d02);
		for (int v = 0; v < 3; v++) mesh.normals[mesh.triangles[t][v]] += n;
	}

	for (int i = 0; i < mesh.normals.size(); i++){
		if (Point3D<double>::Length(mesh.normals[i])> 0) mesh.normals[i] /= Point3D<double>::Length(mesh.normals[i]);
	}
}

int ReadSimpleMesh(SimpleMesh & mesh, const char * fileName){
	mesh.vertices.clear();
	mesh.triangles.clear();
	int file_type;
	std::vector< PlyVertex< double > > ply_vertices;
	bool readFlags[PlyVertex< double >::ReadComponents];
	if (!PlyReadTriangles(fileName, ply_vertices, mesh.triangles, PlyVertex< double >::ReadProperties, readFlags, PlyVertex< double >::ReadComponents, file_type)) return 0;
	mesh.vertices.resize(ply_vertices.size());
	for (int i = 0; i < ply_vertices.size(); i++) mesh.vertices[i] = Point3D<double>(ply_vertices[i].point[0], ply_vertices[i].point[1], ply_vertices[i].point[2]);
	UpdateNormals(mesh);
	return 1;
}

void WriteSimpleMesh(SimpleMesh & mesh, const char * fileName){
	std::vector< PlyVertex< double > > ply_vertices(mesh.vertices.size());
	for (int i = 0; i<mesh.vertices.size(); i++) ply_vertices[i].point = Point3D<double>(mesh.vertices[i][0], mesh.vertices[i][1], mesh.vertices[i][2]);
	PlyWriteTriangles(fileName, ply_vertices, mesh.triangles, PlyVertex< double >::WriteProperties, PlyVertex< double >::WriteComponents, PLY_BINARY_NATIVE);
}

void WriteColoredMesh(ColoredMesh & mesh, const char * fileName){
	std::vector< PlyColorVertex< double > > ply_vertices(mesh.vertices.size());
	for (int i = 0; i<mesh.vertices.size(); i++) ply_vertices[i].point = Point3D<double>(mesh.vertices[i][0], mesh.vertices[i][1], mesh.vertices[i][2]), ply_vertices[i].color = Point3D<double>(mesh.colors[i][0], mesh.colors[i][1], mesh.colors[i][2]);
	PlyWriteTriangles(fileName, ply_vertices, mesh.triangles, PlyColorVertex< double >::WriteProperties, PlyColorVertex< double >::WriteComponents, PLY_BINARY_NATIVE);
}

void WriteTexturedMesh(TexturedMesh & mesh, const char * fileName, const char * atlasName = NULL)
{
	std::vector< PlyTexturedFace< double > > texturedTriangles;
	texturedTriangles.resize(mesh.triangles.size());
	for (int i = 0; i < mesh.triangles.size(); i++){
		texturedTriangles[i].resize(3);
		for (int j = 0; j < 3; j++){
			texturedTriangles[i][j] = mesh.triangles[i][j];
			texturedTriangles[i].texture(j) = Point2D<double>(mesh.textureCoordinates[3 * i + j][0], mesh.textureCoordinates[3 * i + j][1]);
		}
	}

	std::vector< PlyVertex< double > > vertices(mesh.vertices.size());
	for (int i = 0; i < mesh.vertices.size(); i++)vertices[i].point = Point3D<double>(mesh.vertices[i][0], mesh.vertices[i][1], mesh.vertices[i][2]);

	if (atlasName != NULL){
		char ** comments = new char *[1];
		char atlas_comment[256];
		sprintf(atlas_comment, "TextureFile %s", atlasName);
		comments[0] = atlas_comment;
		PlyWritePolygons(fileName, vertices, texturedTriangles, PlyVertex< double >::WriteProperties, PlyVertex< double >::WriteComponents, PlyTexturedFace< double >::WriteProperties, PlyTexturedFace< double >::WriteComponents, PLY_ASCII, comments, 1);
	}
	else{
		PlyWritePolygons(fileName, vertices, texturedTriangles, PlyVertex< double >::WriteProperties, PlyVertex< double >::WriteComponents, PlyTexturedFace< double >::WriteProperties, PlyTexturedFace< double >::WriteComponents, PLY_ASCII);
	}
}

int ReadTexturedMesh(TexturedMesh & mesh, const char * meshName, const char* atlasName , bool verbose )
{
	mesh.vertices.clear();
	mesh.triangles.clear();
	mesh.textureCoordinates.clear();
	int file_type;
	std::vector< PlyVertex< double > > ply_vertices;
	std::vector< PlyTexturedFace< double > > ply_faces;
	if (!PlyReadPolygons(meshName, ply_vertices, ply_faces, PlyVertex< double >::ReadProperties, NULL, PlyVertex< double >::ReadComponents, PlyTexturedFace< double >::ReadProperties, NULL, PlyTexturedFace< double >::ReadComponents, file_type)) return 0;
	
	mesh.vertices.resize(ply_vertices.size());
	for (int i = 0; i < ply_vertices.size(); i++) mesh.vertices[i] = ply_vertices[i].point;

	mesh.triangles.resize(ply_faces.size());
	mesh.textureCoordinates.resize(3*ply_faces.size());
	for (int i = 0; i < ply_faces.size(); i++){
		for (int j = 0; j < 3; j++){
			mesh.triangles[i][j] = ply_faces[i][j];
			mesh.textureCoordinates[3 * i + j] = ply_faces[i].texture(j);
		}
	}
	UpdateNormals(mesh);

	if (atlasName != NULL){

		char* ext = GetFileExtension( atlasName );
	
		if( !strcasecmp( ext , "normap" ) )
		{
			if( verbose ) printf( "Reading normal texture\n" );
			if( !ReadBinaryImage( mesh.texture , atlasName ) ){ printf( "Unable to read texture: %s\n" , atlasName ) ; delete[] ext ; return 0; }
		}
		else
		{
			if( verbose ) printf( "Reading color texture\n" );
			if( !mesh.texture.read( atlasName ) ){ printf( "Unable to read texture: %s\n" , atlasName ) ; delete[] ext ; return 0; }
		}
		delete[] ext;
	}

	return 1;
}

double GetMeshArea(const SimpleMesh & mesh){
	double meshArea = 0;
	for (int t = 0; t < mesh.triangles.size(); t++){
		Point3D<double> d01 = mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]];
		Point3D<double> d02 = mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]];
		Point3D<double> n = Point3D<double>::CrossProduct(d01, d02);
		meshArea += (Point3D<double>::Length(n) / 2.0);
	}
	return meshArea;
}

Point3D<double> GetMeshCentroid(const SimpleMesh & mesh) {
	double meshArea = 0;
	Point3D<double> centroid(0,0,0);
	for (int t = 0; t < mesh.triangles.size(); t++) {
		Point3D<double> d01 = mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]];
		Point3D<double> d02 = mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]];
		Point3D<double> n = Point3D<double>::CrossProduct(d01, d02);
		double tArea = (Point3D<double>::Length(n) / 2.0);
		Point3D<double> baricenter = (mesh.vertices[mesh.triangles[t][0]] + mesh.vertices[mesh.triangles[t][1]] + mesh.vertices[mesh.triangles[t][2]]) / 3.0;
		centroid += baricenter*tArea;
		meshArea += tArea;
	}
	return (centroid/meshArea);
}

double GetMeshRadius(const SimpleMesh & mesh, const Point3D<double> & centroid) {
	double radius = 0;
	for (int v = 0; v < mesh.vertices.size(); v++)radius = std::max<double>(radius, Point3D<double>::Length(mesh.vertices[v] - centroid));
	return radius;
}

void InitializeMeshMatrices(const SimpleMesh & mesh, Eigen::SparseMatrix<double> & meshMassMatrix, Eigen::SparseMatrix<double> & meshStiffnessMatrix){
	double meshMass = 0;
	for (int t = 0; t < mesh.triangles.size(); t++){
		Point3D<double> p[3] = { mesh.vertices[mesh.triangles[t][0]], mesh.vertices[mesh.triangles[t][1]], mesh.vertices[mesh.triangles[t][2]] };
		Point3D<double> d[2] = { mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]], mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]] };
		SquareMatrix<double, 2> g;
		for (int k = 0; k < 2; k++)for (int l = 0; l < 2; l++) g(k, l) = Point3D<double>::Dot(d[k], d[l]);
		double _mass = sqrt(g.determinant()) / 2.0;
		meshMass += _mass;
	}

	double cumMass = 0;
	std::vector<Eigen::Triplet<double>> massTriplets;
	std::vector<Eigen::Triplet<double>> stiffnessTriplets;
	for (int t = 0; t < mesh.triangles.size(); t++){
		Point3D<double> p[3] = { mesh.vertices[mesh.triangles[t][0]], mesh.vertices[mesh.triangles[t][1]], mesh.vertices[mesh.triangles[t][2]] };
		Point3D<double> d[2] = { mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]], mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]] };
		SquareMatrix<double, 2> g;
		for (int k = 0; k < 2; k++)for (int l = 0; l < 2; l++) g(k, l) = Point3D<double>::Dot(d[k], d[l]);
		g /= meshMass;

		double _mass = sqrt(g.determinant()) / 2.0;
		cumMass += _mass;
		for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) massTriplets.push_back(Eigen::Triplet<double>(mesh.triangles[t][k], mesh.triangles[t][l], k == l ? _mass / 6.0 : _mass / 12.0));
		Point2D<double> e[3] = { Point2D<double>(-1, -1), Point2D<double>(1, 0), Point2D<double>(0, 1) };
		SquareMatrix<double, 2> g_inverse = g.inverse();
		for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++){
			double stiffness_coeff = Point2D<double>::Dot(e[k], g_inverse*e[l]) * _mass;
			stiffnessTriplets.push_back(Eigen::Triplet<double>(mesh.triangles[t][k], mesh.triangles[t][l], stiffness_coeff));
		}
	}

	//printf("Mesh mass %f \n", cumMass);

	meshStiffnessMatrix.resize(mesh.vertices.size(), mesh.vertices.size());
	meshStiffnessMatrix.setFromTriplets(stiffnessTriplets.begin(), stiffnessTriplets.end());

	meshMassMatrix.resize(mesh.vertices.size(), mesh.vertices.size());
	meshMassMatrix.setFromTriplets(massTriplets.begin(), massTriplets.end());
}

int InitializeBoundaryEdges(const TexturedMesh & mesh, std::vector<int> & boundaryEdges) {

	bool isClosedMesh = true;

	std::unordered_map<unsigned long long, int> edgeIndex;
	for (int i = 0; i < mesh.triangles.size(); i++) {
		for (int k = 0; k < 3; k++) {
			unsigned long long  edgeKey = SetMeshEdgeKey(mesh.triangles[i][k], mesh.triangles[i][(k + 1) % 3]);
			if (edgeIndex.find(edgeKey) == edgeIndex.end()) {
				edgeIndex[edgeKey] = 3 * i + k;
			}
			else {
				printf("Non manifold mesh!! \n");
				return 0;
			}
		}
	}

	for (int i = 0; i < mesh.triangles.size(); i++) {
		for (int k = 0; k < 3; k++) {
			int currentEdgeIndex = 3 * i + k;
			unsigned long long edgeKey = SetMeshEdgeKey(mesh.triangles[i][(k + 1) % 3], mesh.triangles[i][k]);
			if (edgeIndex.find(edgeKey) != edgeIndex.end()) {
				int oppositeEdgeIndex = edgeIndex[edgeKey];
				int tIndex = oppositeEdgeIndex / 3;
				int kIndex = oppositeEdgeIndex % 3;
				if (mesh.textureCoordinates[3 * i + ((k + 1) % 3)][0] == mesh.textureCoordinates[3 * tIndex + kIndex][0] &&
					mesh.textureCoordinates[3 * i + ((k + 1) % 3)][1] == mesh.textureCoordinates[3 * tIndex + kIndex][1] &&
					mesh.textureCoordinates[3 * i + k][0] == mesh.textureCoordinates[3 * tIndex + ((kIndex + 1) % 3)][0] &&
					mesh.textureCoordinates[3 * i + k][1] == mesh.textureCoordinates[3 * tIndex + ((kIndex + 1) % 3)][1]
					) {
					//Do nothing
				}
				else {
					if (currentEdgeIndex < oppositeEdgeIndex) {
						boundaryEdges.push_back(currentEdgeIndex);
					}
				}
			}
			else {
				isClosedMesh = false;
				boundaryEdges.push_back(currentEdgeIndex);
			}
		}
	}
	return 1;
}

void SmoothSignal(const Eigen::SparseMatrix<double> & massMatrix, const Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > > & solver, std::vector<Point3D<double>> & signal, bool normalize = false){
	int vCount = (int)signal.size();

	for (int c = 0; c < 3; c++){
		Eigen::VectorXd x0;
		x0.resize(vCount);
		for (int i = 0; i < vCount; i++) x0[i] = signal[i][c];
		Eigen::VectorXd rhs = massMatrix*x0;
		Eigen::VectorXd sol = solver.solve(rhs);
		for (int i = 0; i < vCount; i++) signal[i][c] = sol[i];
	}

	if (normalize) for (int i = 0; i < signal.size(); i++) signal[i] /= Point3D<double>::Length(signal[i]);
}

void SmoothSignal(const Eigen::SparseMatrix<double> & massMatrix, const Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > > & solver,  std::vector<double> & signal){
	int vCount = (int)signal.size();
	Eigen::VectorXd x0;
	x0.resize(vCount);
	for (int i = 0; i < vCount; i++) x0[i] = signal[i];
	Eigen::VectorXd rhs = massMatrix*x0;
	Eigen::VectorXd sol = solver.solve(rhs);
	for (int i = 0; i < vCount; i++) signal[i] = sol[i];
}

template <class Vertex, class TCoord>
void SubdivideMeshAndTextureCoordinates(std::vector<Vertex> & vertices, std::vector<TriangleIndex> & triangles, std::vector<TCoord> & tCoordinates) {

#define EDGE_KEY( i1 , i2 ) ( (i1)>(i2) ? ( ( (long long) (i1) )<<32 ) | ( (long long) (i2) ) : ( ( (long long) (i2) )<<32 ) | ( (long long) (i1) ) )

	std::vector< TriangleIndex > _triangles(triangles.size() * 4);
	std::vector< Vertex > _vertices = vertices;
	std::vector< TCoord > _tCoordinates(triangles.size() * 12);

	std::unordered_map< long long, int > vMap;
	for (int i = 0; i<triangles.size(); i++)
	{
		long long keys[] = { EDGE_KEY(triangles[i][1], triangles[i][2]), EDGE_KEY(triangles[i][2], triangles[i][0]), EDGE_KEY(triangles[i][0], triangles[i][1]) };
		int eIndex[3];
		for (int j = 0; j<3; j++)
		{
			if (vMap.find(keys[j]) == vMap.end())
			{
				vMap[keys[j]] = eIndex[j] = (int)_vertices.size();
				_vertices.push_back((vertices[triangles[i][(j + 1) % 3]] + vertices[triangles[i][(j + 2) % 3]]) / 2);
			}
			else eIndex[j] = vMap[keys[j]];
		}

		TCoord cornerTCoordinates[3] = { tCoordinates[3 * i], tCoordinates[3 * i + 1], tCoordinates[3 * i + 2] };
		TCoord midTCoordinates[3] = { (cornerTCoordinates[1] + cornerTCoordinates[2]) / 2.0, (cornerTCoordinates[2] + cornerTCoordinates[0]) / 2.0, (cornerTCoordinates[0] + cornerTCoordinates[1]) / 2.0 };


		_triangles[4 * i + 0] = TriangleIndex(eIndex[0], eIndex[1], eIndex[2]);
		_tCoordinates[12 * i] = midTCoordinates[0];
		_tCoordinates[12 * i + 1] = midTCoordinates[1];
		_tCoordinates[12 * i + 2] = midTCoordinates[2];

		_triangles[4 * i + 1] = TriangleIndex(triangles[i][0], eIndex[2], eIndex[1]);
		_tCoordinates[12 * i + 3] = cornerTCoordinates[0];
		_tCoordinates[12 * i + 4] = midTCoordinates[2];
		_tCoordinates[12 * i + 5] = midTCoordinates[1];

		_triangles[4 * i + 2] = TriangleIndex(eIndex[2], triangles[i][1], eIndex[0]);
		_tCoordinates[12 * i + 6] = midTCoordinates[2];
		_tCoordinates[12 * i + 7] = cornerTCoordinates[1];
		_tCoordinates[12 * i + 8] = midTCoordinates[0];

		_triangles[4 * i + 3] = TriangleIndex(eIndex[1], eIndex[0], triangles[i][2]);
		_tCoordinates[12 * i + 9] = midTCoordinates[1];
		_tCoordinates[12 * i + 10] = midTCoordinates[0];
		_tCoordinates[12 * i + 11] = cornerTCoordinates[2];

	}
	triangles = _triangles;
	vertices = _vertices;
	tCoordinates = _tCoordinates;
#undef EDGE_KEY

}

void TexturedMeshSubdivide(TexturedMesh & mesh, int iters) {
	for (int i = 0; i < iters; i++) SubdivideMeshAndTextureCoordinates(mesh.vertices, mesh.triangles, mesh.textureCoordinates);
	UpdateNormals(mesh);
}

#endif//SIMPLE_MESH_INCLUDED
