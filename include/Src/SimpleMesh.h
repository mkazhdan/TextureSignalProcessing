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
#include <Misha/Miscellany.h>
#include <Src/VectorIO.h>

template< typename GeometryReal >
class MeshSample
{
public:
	MeshSample( void ) : tID(-1){}
	MeshSample( int _tID , Point2D< GeometryReal > _barycentricCoords ) : tID(_tID) , barycentricCoords(_barycentricCoords) { }
	int tID;
	Point2D< GeometryReal > barycentricCoords;
};

unsigned long long SetMeshEdgeKey( const unsigned long i0 , const unsigned long i1 )
{
	return ( ( ( static_cast< unsigned long long >(i0) << 32 ) & 0xFFFFFFFF00000000) | ( static_cast< unsigned long long >(i1) & 0x00000000FFFFFFFF ) );
}

void GetMeshEdgeIndices( unsigned long long key , unsigned long & i0 , unsigned long & i1 )
{
	i1 = static_cast< unsigned long >(  key      & 0x00000000FFFFFFFF );
	i0 = static_cast< unsigned long >( (key>>32) & 0x00000000FFFFFFFF );
}

template< typename GeometryReal >
class SimpleMesh
{
public:
	std::vector< Point3D< GeometryReal > > vertices;
	std::vector< Point3D< GeometryReal > > normals;
	std::vector< TriangleIndex > triangles;

	void updateNormals( void )
	{
		normals.clear();
		normals.resize( vertices.size() );
		for( int t=0 ; t<triangles.size() ; t++ )
		{
			Point3D< GeometryReal > d01 = vertices[ triangles[t][1] ] - vertices[ triangles[t][0] ];
			Point3D< GeometryReal > d02 = vertices[ triangles[t][2] ] - vertices[ triangles[t][0] ];
			Point3D< GeometryReal > n = Point3D< GeometryReal >::CrossProduct( d01 , d02 );
			for( int v=0 ; v<3 ; v++ ) normals[ triangles[t][v] ] += n;
		}

		for( int i=0 ; i<normals.size() ; i++ )
		{
			GeometryReal l = Point3D< GeometryReal >::Length( normals[i] );
			if( l>0 ) normals[i] /= l;
		}
	}
	void read( const char *fileName )
	{
		vertices.clear();
		triangles.clear();
		int file_type;
		std::vector< PlyVertex< GeometryReal > > ply_vertices;
		if ( !PlyReadTriangles( fileName , ply_vertices , triangles , PlyVertex< GeometryReal >::ReadProperties , NULL , PlyVertex< GeometryReal >::ReadComponents, file_type ) )
			Miscellany::Throw( "Failed to read ply file: %s" , fileName );
		vertices.resize( ply_vertices.size() );
		for( int i=0 ; i<ply_vertices.size() ; i++ ) vertices[i] = ply_vertices[i].point;
		updateNormals();
	}
	void write( const char *fileName ) const
	{
		std::vector< PlyVertex< GeometryReal > > ply_vertices( vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ ) ply_vertices[i].point = vertices[i];
		PlyWriteTriangles( fileName , ply_vertices , triangles , PlyVertex< GeometryReal >::WriteProperties , PlyVertex< GeometryReal >::WriteComponents , PLY_BINARY_NATIVE );
	}

	GeometryReal area( void ) const
	{
		GeometryReal meshArea = 0;
		for( int t=0 ; t<triangles.size() ; t++ )
		{
			Point3D< GeometryReal > d01 = vertices[ triangles[t][1] ] - vertices[ triangles[t][0] ];
			Point3D< GeometryReal > d02 = vertices[ triangles[t][2] ] - vertices[ triangles[t][0] ];
			Point3D< GeometryReal > n = Point3D< GeometryReal >::CrossProduct( d01 , d02 );
			meshArea += Point3D< GeometryReal >::Length(n) / 2;
		}
		return meshArea;
	}

	Point3D< GeometryReal > centroid( void ) const
	{
		GeometryReal meshArea = 0;
		Point3D< GeometryReal > centroid;
		for( int t=0 ; t<triangles.size() ; t++ )
		{
			Point3D< GeometryReal > d01 = vertices[ triangles[t][1] ] - vertices[ triangles[t][0] ];
			Point3D< GeometryReal > d02 = vertices[ triangles[t][2] ] - vertices[ triangles[t][0] ];
			Point3D< GeometryReal > n = Point3D< GeometryReal >::CrossProduct( d01 , d02 );
			GeometryReal tArea = Point3D< GeometryReal >::Length(n) / 2;
			Point3D< GeometryReal > baricenter = ( vertices[ triangles[t][0] ] + vertices[ triangles[t][1] ] + vertices[ triangles[t][2]] ) / 3;
			centroid += baricenter*tArea;
			meshArea += tArea;
		}
		return centroid/meshArea;
	}

	GeometryReal radius( const Point3D< GeometryReal > &centroid )
	{
		GeometryReal radius = 0;
		for( int v=0 ; v<vertices.size() ; v++ ) radius = std::max< GeometryReal >( radius , Point3D< GeometryReal >::Length( vertices[v] - centroid ) );
		return radius;
	}

	template< typename MatrixReal >
	void initializeMeshMatrices( Eigen::SparseMatrix< MatrixReal > &mass , Eigen::SparseMatrix< MatrixReal > &stiffness )
	{
		MatrixReal meshMass = 0;
		for( int t=0 ; t<triangles.size() ; t++ )
		{
			Point3D< MatrixReal > p[3] = { vertices[ triangles[t][0] ] , vertices[ triangles[t][1] ] , vertices[ triangles[t][2] ] };
			Point3D< MatrixReal > d[2] = { vertices[ triangles[t][1] ] - vertices[ triangles[t][0] ] , vertices[ triangles[t][2] ] - vertices[ triangles[t][0] ] };
			SquareMatrix< MatrixReal , 2 > g;
			for( int k=0 ; k<2 ; k++ ) for( int l=0 ; l<2 ; l++ ) g(k,l) = Point3D<MatrixReal>::Dot( d[k] , d[l] );
			MatrixReal _mass = sqrt(g.determinant()) / 2.0;
			meshMass += _mass;
		}

		MatrixReal cumMass = 0;
		std::vector< Eigen::Triplet< MatrixReal > > mTriplets;
		std::vector< Eigen::Triplet< MatrixReal > > sTriplets;
		for( int t=0 ; t<triangles.size() ; t++ )
		{
			Point3D< MatrixReal > p[3] = { vertices[ triangles[t][0] ] , vertices[ triangles[t][1] ] , vertices[ triangles[t][2] ] };
			Point3D< MatrixReal > d[2] = { vertices[ triangles[t][1] ] - vertices[ triangles[t][0] ] , vertices[ triangles[t][2] ] - vertices[ triangles[t][0] ] };
			SquareMatrix< MatrixReal , 2 > g;
			for( int k=0 ; k<2 ; k++ ) for( int l=0 ; l<2 ; l++ ) g(k,l) = Point3D< MatrixReal >::Dot( d[k] , d[l] );
			g /= meshMass;

			MatrixReal _mass = (MatrixReal)sqrt( g.determinant() ) / 2;
			cumMass += _mass;
			for( int k=0 ; k<3 ; k++ ) for( int l=0 ; l<3 ; l++ ) mTriplets.push_back( Eigen::Triplet< MatrixReal >( triangles[t][k] , triangles[t][l] , k == l ? _mass / 6.0 : _mass / 12.0 ) );
			Point2D< MatrixReal > e[3] = { Point2D< MatrixReal >(-1,-1) , Point2D< MatrixReal >(1,0) , Point2D< MatrixReal >(0,1) };
			SquareMatrix< MatrixReal , 2 > g_inverse = g.inverse();
			for( int k=0 ; k<3 ; k++ ) for( int l=0 ; l<3 ; l++ )
			{
				MatrixReal stiffness_coeff = Point2D< MatrixReal >::Dot( e[k] , g_inverse*e[l] ) * _mass;
				sTriplets.push_back( Eigen::Triplet< MatrixReal >( triangles[t][k] , triangles[t][l] , stiffness_coeff ) );
			}
		}

		stiffness.resize( vertices.size() , vertices.size() );
		stiffness.setFromTriplets( sTriplets.begin() , sTriplets.end() );

		mass.resize( vertices.size() , vertices.size() );
		mass.setFromTriplets( mTriplets.begin() , mTriplets.end() );
	}
};

template< typename GeometryReal >
class ColoredMesh : public SimpleMesh< GeometryReal >
{
public:
	using SimpleMesh< GeometryReal >::vertices;
	using SimpleMesh< GeometryReal >::triangles;
	std::vector< Point3D< GeometryReal > > colors;

	void write( const char *fileName ) const
	{
		std::vector< PlyColorVertex< GeometryReal > > ply_vertices( vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ ) ply_vertices[i].point = vertices[i] , ply_vertices[i].color = colors[i];
		PlyWriteTriangles( fileName , ply_vertices , triangles , PlyColorVertex< GeometryReal >::WriteProperties , PlyColorVertex< GeometryReal >::WriteComponents , PLY_BINARY_NATIVE );
	}

};

template< typename GeometryReal , typename ImageReal=float >
class TexturedMesh : public SimpleMesh< GeometryReal >
{
protected:
	static void _Subdivide( std::vector< Point3D< GeometryReal > > &vertices , std::vector< TriangleIndex > &triangles , std::vector< Point2D< GeometryReal > > & tCoordinates )
	{
#define EDGE_KEY( i1 , i2 ) ( (i1)>(i2) ? ( ( (long long) (i1) )<<32 ) | ( (long long) (i2) ) : ( ( (long long) (i2) )<<32 ) | ( (long long) (i1) ) )

		std::vector< TriangleIndex > _triangles(triangles.size() * 4);
		std::vector< Point3D< GeometryReal > > _vertices = vertices;
		std::vector< Point2D< GeometryReal > > _tCoordinates(triangles.size() * 12);

		std::unordered_map< long long , int > vMap;
		for( int i=0 ; i<triangles.size() ; i++)
		{
			long long keys[] = { EDGE_KEY( triangles[i][1] , triangles[i][2] ) , EDGE_KEY( triangles[i][2] , triangles[i][0] ) , EDGE_KEY( triangles[i][0] , triangles[i][1] ) };
			int eIndex[3];
			for( int j=0 ; j<3 ; j++ )
			{
				if( vMap.find( keys[j] )==vMap.end() )
				{
					vMap[ keys[j] ] = eIndex[j] = (int)_vertices.size();
					_vertices.push_back( ( vertices[ triangles[i][(j+1)%3] ] + vertices[ triangles[i][(j+2)%3] ] ) / 2 );
				}
				else eIndex[j] = vMap[ keys[j] ];
			}

			Point2D< GeometryReal > cornerTCoordinates[] = { tCoordinates[3*i] , tCoordinates[3*i+1] , tCoordinates[3*i+2] };
			Point2D< GeometryReal > midTCoordinates[] = { ( cornerTCoordinates[1] + cornerTCoordinates[2] ) / 2 , ( cornerTCoordinates[2] + cornerTCoordinates[0] ) / 2 , ( cornerTCoordinates[0] + cornerTCoordinates[1] ) / 2 };

			_triangles[4*i+0] = TriangleIndex( eIndex[0] , eIndex[1] , eIndex[2] );
			_tCoordinates[12*i+0] = midTCoordinates[0];
			_tCoordinates[12*i+1] = midTCoordinates[1];
			_tCoordinates[12*i+2] = midTCoordinates[2];

			_triangles[4*i+1] = TriangleIndex( triangles[i][0] , eIndex[2] , eIndex[1] );
			_tCoordinates[12*i+3] = cornerTCoordinates[0];
			_tCoordinates[12*i+4] = midTCoordinates[2];
			_tCoordinates[12*i+5] = midTCoordinates[1];

			_triangles[4*i+2] = TriangleIndex( eIndex[2] , triangles[i][1] , eIndex[0] );
			_tCoordinates[12*i+6] = midTCoordinates[2];
			_tCoordinates[12*i+7] = cornerTCoordinates[1];
			_tCoordinates[12*i+8] = midTCoordinates[0];

			_triangles[4*i+3] = TriangleIndex( eIndex[1] , eIndex[0] , triangles[i][2] );
			_tCoordinates[12*i+9] = midTCoordinates[1];
			_tCoordinates[12*i+10] = midTCoordinates[0];
			_tCoordinates[12*i+11] = cornerTCoordinates[2];

		}
		triangles = _triangles;
		vertices = _vertices;
		tCoordinates = _tCoordinates;
#undef EDGE_KEY
	}
public:
	using SimpleMesh< GeometryReal >::vertices;
	using SimpleMesh< GeometryReal >::triangles;
	std::vector< Point2D< GeometryReal > > textureCoordinates;
	Image< Point3D< ImageReal > > texture;

	void write( const char *fileName , const char *atlasName=NULL ) const
	{
		std::vector< PlyTexturedFace< GeometryReal > > texturedTriangles;
		texturedTriangles.resize( triangles.size() );
		for( int i=0 ; i<triangles.size() ; i++ )
		{
			texturedTriangles[i].resize(3);
			for( int j=0 ; j<3 ; j++ )
			{
				texturedTriangles[i][j] = triangles[i][j];
				texturedTriangles[i].texture(j) = textureCoordinates[3*i+j];
			}
		}

		std::vector< PlyVertex< GeometryReal > > vertices( vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ ) vertices[i].point = vertices[i];

		if( atlasName )
		{
			char **comments = new char*[1];
			char atlas_comment[256];
			sprintf( atlas_comment , "TextureFile %s" , atlasName );
			comments[0] = atlas_comment;
			PlyWritePolygons( fileName , vertices , texturedTriangles , PlyVertex< GeometryReal >::WriteProperties , PlyVertex< GeometryReal >::WriteComponents , PlyTexturedFace< GeometryReal >::WriteProperties , PlyTexturedFace< GeometryReal >::WriteComponents , PLY_BINARY_NATIVE , comments , 1 );
			delete[] comments;
		}
		else PlyWritePolygons( fileName , vertices , texturedTriangles , PlyVertex< GeometryReal >::WriteProperties , PlyVertex< GeometryReal >::WriteComponents , PlyTexturedFace< GeometryReal >::WriteProperties , PlyTexturedFace< GeometryReal >::WriteComponents , PLY_BINARY_NATIVE );
	}

	void read( const char *meshName , const char *atlasName , bool verbose )
	{
		vertices.clear();
		triangles.clear();
		textureCoordinates.clear();
#if 1 // NEW_CODE
		char * ext = GetFileExtension( meshName );
		if( !strcasecmp( ext , "ply" ) )
		{
			int file_type;
			std::vector< PlyVertex< GeometryReal > > ply_vertices;
			std::vector< PlyTexturedFace< GeometryReal > > ply_faces;
			if( !PlyReadPolygons( meshName , ply_vertices , ply_faces , PlyVertex< GeometryReal >::ReadProperties , NULL , PlyVertex< GeometryReal >::ReadComponents , PlyTexturedFace< GeometryReal >::ReadProperties , NULL , PlyTexturedFace< GeometryReal >::ReadComponents , file_type ) )
				Miscellany::Throw( "Failed to read ply file: %s" , meshName );

			vertices.resize( ply_vertices.size() );
			for( int i=0 ; i<ply_vertices.size() ; i++ ) vertices[i] = ply_vertices[i].point;

			triangles.resize( ply_faces.size() );
			textureCoordinates.resize( 3*ply_faces.size() );
			for( int i=0 ; i<ply_faces.size() ; i++ ) for( int j=0 ; j<3 ; j++ )
			{
				triangles[i][j] = ply_faces[i][j];
				textureCoordinates[3*i+j] = ply_faces[i].texture(j);
			}
		}
		else if( !strcasecmp( ext , "obj" ) )
		{
			struct ObjFaceIndex
			{
				int vIndex , tIndex;
			};
			std::vector< Point3D< GeometryReal > > obj_vertices;
			std::vector< Point2D< GeometryReal > > obj_textures;
			std::vector< std::vector< ObjFaceIndex > > obj_faces;
			std::ifstream in( meshName );
			if( !in.is_open() ) Miscellany::ErrorOut( "Could not open file for reading: " , std::string(meshName) );

			std::string( line );
			unsigned int count = 0;
			while( std::getline( in , line ) )
			{
				std::stringstream ss;
				if( line[0]=='v' && line[1]==' ' )
				{
					line = line.substr(2);
					std::stringstream ss( line );
					Point3D< GeometryReal > p;
					ss >> p[0] >> p[1] >> p[2];
					obj_vertices.push_back( p );
				}
				else if( line[0]=='v' && line[1]=='t' && line[2]==' ' )
				{
					line = line.substr(3);
					std::stringstream ss( line );
					Point2D< GeometryReal > p;
					ss >> p[0] >> p[1];
					obj_textures.push_back( p );
				}
				else if( line[0]=='f' && line[1]==' ' )
				{
					std::vector< ObjFaceIndex > face;
					line = line.substr(2);
					std::stringstream ss( line );
					std::string token;
					while( std::getline( ss , token , ' ' ) )
					{
						std::stringstream _ss(token);
						std::string _token;
						unsigned int count = 0;
						ObjFaceIndex idx;
						while( std::getline( _ss , _token , '/' ) )
						{
							if     ( count==0 ) idx.vIndex = std::stoi( _token );
							else if( count==1 ) idx.tIndex = std::stoi( _token );
							count++;
						}
						face.push_back( idx );
					}
					if( face.size()!=3 ) Miscellany::ErrorOut( "Expectred triangular face: " , face.size() );
					obj_faces.push_back( face );
				}
				count++;
			}

			vertices.resize( obj_vertices.size() );
			for( int i=0 ; i<obj_vertices.size() ; i++ ) vertices[i] = obj_vertices[i];

			triangles.resize( obj_faces.size() );
			textureCoordinates.resize( 3*obj_faces.size() );
#if 1
			for( int i=0 ; i<obj_faces.size() ; i++ ) for( int j=0 ; j<3 ; j++ )
			{
				if     ( obj_faces[i][j].vIndex>0 ) triangles[i][j] = obj_faces[i][j].vIndex-1;
				else if( obj_faces[i][j].vIndex<0 ) triangles[i][j] = (int)obj_vertices.size() + obj_faces[i][j].vIndex;
				else Miscellany::ErrorOut( "Zero vertex index unexpected in obj file" );

				if     ( obj_faces[i][j].tIndex>0 ) textureCoordinates[3*i+j] = obj_textures[ obj_faces[i][j].tIndex-1 ];
				else if( obj_faces[i][j].tIndex<0 ) textureCoordinates[3*i+j] = obj_textures[ (int)obj_textures.size() + obj_faces[i][j].tIndex ];
				else Miscellany::ErrorOut( "Zero texture index unexpected in obj file" );
			}
#else
			Miscellany::ErrorOut( "no obj support yet" );
#endif

		}
		else Miscellany::ErrorOut( "Unrecognized file extension: " , std::string( meshName ) );
		delete[] ext;

#else // !NEW_CODE
		int file_type;
		std::vector< PlyVertex< GeometryReal > > ply_vertices;
		std::vector< PlyTexturedFace< GeometryReal > > ply_faces;
		if( !PlyReadPolygons( meshName , ply_vertices , ply_faces , PlyVertex< GeometryReal >::ReadProperties , NULL , PlyVertex< GeometryReal >::ReadComponents , PlyTexturedFace< GeometryReal >::ReadProperties , NULL , PlyTexturedFace< GeometryReal >::ReadComponents , file_type ) )
			Miscellany::Throw( "Failed to read ply file: %s" , meshName );

		vertices.resize( ply_vertices.size() );
		for( int i=0 ; i<ply_vertices.size() ; i++ ) vertices[i] = ply_vertices[i].point;

		triangles.resize( ply_faces.size() );
		textureCoordinates.resize( 3*ply_faces.size() );
		for( int i=0 ; i<ply_faces.size() ; i++ ) for( int j=0 ; j<3 ; j++ )
		{
			triangles[i][j] = ply_faces[i][j];
			textureCoordinates[3*i+j] = ply_faces[i].texture(j);
		}
#endif // NEW_CODE
		SimpleMesh< GeometryReal >::updateNormals();

		if( atlasName )
		{

			char *ext = GetFileExtension( atlasName );

			if( !strcasecmp( ext , "normap" ) )
			{
				if( verbose ) printf( "Reading normal texture\n" );
				ReadBinaryImage( texture , atlasName );
			}
			else
			{
				if( verbose ) printf( "Reading color texture\n" );
				texture.read( atlasName );
			}
			delete[] ext;
		}
	}

	void initializeBoundaryEdges( std::vector< int > &boundaryEdges ) const
	{
		bool isClosedMesh = true;

		std::unordered_map< unsigned long long , int > edgeIndex;
		for( int i=0 ; i<triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ )
		{
			unsigned long long edgeKey = SetMeshEdgeKey( triangles[i][k] , triangles[i][(k+1)%3] );
			if( edgeIndex.find(edgeKey)==edgeIndex.end() ) edgeIndex[edgeKey] = 3*i+k;
			else Miscellany::Throw( "Non manifold mesh" );
		}

		for( int i=0 ; i<triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ )
		{
			int currentEdgeIndex = 3*i+k;
			unsigned long long edgeKey = SetMeshEdgeKey( triangles[i][(k+1)%3] , triangles[i][k] );
			if( edgeIndex.find(edgeKey)!=edgeIndex.end() )
			{
				int oppositeEdgeIndex = edgeIndex[edgeKey];
				int tIndex = oppositeEdgeIndex/3;
				int kIndex = oppositeEdgeIndex%3;
				if( textureCoordinates[ 3*i+(k+1)%3 ][0] == textureCoordinates[ 3*tIndex+kIndex ][0] &&
					textureCoordinates[ 3*i+(k+1)%3 ][1] == textureCoordinates[ 3*tIndex+kIndex ][1] &&
					textureCoordinates[ 3*i+k ][0] == textureCoordinates[ 3*tIndex+(kIndex+1)%3 ][0] &&
					textureCoordinates[ 3*i+k ][1] == textureCoordinates[ 3*tIndex+(kIndex+1)%3 ][1] )
				{
				}
				else
				{
					if( currentEdgeIndex< oppositeEdgeIndex ) boundaryEdges.push_back( currentEdgeIndex );
				}
			}
			else
			{
				isClosedMesh = false;
				boundaryEdges.push_back( currentEdgeIndex );
			}
		}
	}
	void subdivide( int iters )
	{
		for( int i=0 ; i<iters ; i++ ) _Subdivide( vertices , triangles , textureCoordinates );
		SimpleMesh< GeometryReal >::updateNormals();
	}
};
#endif//SIMPLE_MESH_INCLUDED
