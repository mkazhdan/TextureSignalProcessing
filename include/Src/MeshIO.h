/*
Copyright (c) 2025, Michael Kazhdan
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

#ifndef PLY_IO_INCLUDED
#define PLY_IO_INCLUDED

#include <string>
#include "Misha/Ply.h"
#include "Misha/PlyVertexData.h"
#include "Misha/Geometry.h"
#include "Misha/CmdLineParser.h"

namespace MishaK
{

	template< typename Index , typename Real , unsigned int Dim , unsigned int TDim >
	void ReadTexturedMesh( std::string fileName , std::vector< Point< Real , Dim > > &vertices , std::vector< Point< Real , TDim > > &textureCoordinates , std::vector< SimplexIndex< 2 , Index > > &simplices );

	template< typename Index , typename Real > struct PlyTexturedFace;

	/////////////////////////////////////////////////////////////////////////////////////

	template< typename Index , typename Real , unsigned int Dim , unsigned int TDim >
	void ReadTexturedMesh( std::string fileName , std::vector< Point< Real , Dim > > &vertices , std::vector< Point< Real , TDim > > &textureCoordinates , std::vector< SimplexIndex< 2 , Index > > &simplices )
	{
		static const unsigned int K=2;
		std::string ext = ToLower( GetFileExtension( fileName ) );
		if( ext==std::string( "ply" ) )
		{
			using VertexFactory = VertexFactory::PositionFactory< Real , Dim >;
			using Vertex = typename VertexFactory::VertexType;
			using Face = PlyTexturedFace< unsigned int , Real >;

			VertexFactory factory;

			std::vector< Vertex > inVertices;
			std::vector< Face > inFaces;
			bool *vFlags = new bool[ factory.plyReadNum() ];
			bool *fFlags = new bool[ Face::NumProperties ];
			int file_type;

			try{ file_type = PLY::ReadPolygons( fileName , factory , inVertices , inFaces , Face::Properties , Face::NumProperties , vFlags , fFlags ); }
			catch( const Exception & ){ MK_ERROR_OUT( "Failed to read ply file: " , fileName ); }

			if( !fFlags[0] ) MK_ERROR_OUT( "Failed to read face indices" );
			if( !fFlags[1] ) MK_ERROR_OUT( "Failed to read face textures" );
			delete[] vFlags;
			delete[] fFlags;

			if constexpr( K==2 )
			{
				size_t faceNum = inFaces.size();
				MinimalAreaTriangulation< Real , Dim > mat;
				for( unsigned int i=(unsigned int)inFaces.size() ; i!=0 ; i-- )
				{
					Face &face = inFaces[i-1];
					Face oldFace = face;

					if( face.size()>3 )
					{
						std::vector< Point< Real , Dim > > _vertices( face.size() );
						std::vector< SimplexIndex< K > > _triangles;
						for( unsigned int j=0 ; j<(unsigned int)face.size() ; j++ ) _vertices[j] = inVertices[ face[j] ];
						mat.GetTriangulation( _vertices , _triangles );

						auto TriangleToFace = [&]( SimplexIndex< K > si )
							{
								Face face;
								face.resize(K+1);
								for( unsigned int k=0 ; k<=K ; k++ ) face[k] = oldFace[ si[k] ] , face.texture(k) = oldFace.texture( si[k] );
								return face;
							};

						face = TriangleToFace( _triangles[0] );
						for( unsigned int j=1 ; j<_triangles.size() ; j++ ) inFaces.push_back( TriangleToFace( _triangles[j] ) );
					}
				}
				//		if( inFaces.size()!=faceNum ) WARN( "Triangulated: " , faceNum , " -> " , inFaces.size() );
			}

			for( unsigned int i=0 ; i<inFaces.size() ; i++ )
				if( inFaces[i].nr_vertices!=(K+1) ) MK_ERROR_OUT( "Face is not a simplex" );
				else if( inFaces[i].nr_uv_coordinates!=(K+1)*K ) MK_ERROR_OUT( "Unexpected number of texture coordinates: " , inFaces[i].nr_uv_coordinates , " != " , K+1 , " * " , K );

			simplices.resize( inFaces.size() );
			textureCoordinates.resize( inFaces.size()*(K+1) );
			for( unsigned int i=0 ; i<inFaces.size() ; i++ ) for( unsigned int k=0 ; k<=K ; k++ )
			{
				simplices[i][k] = inFaces[i][k];
				textureCoordinates[ i*(K+1)+k ] = inFaces[i].texture(k);
			}

			vertices.resize( inVertices.size() );
			for( unsigned int i=0 ; i<inVertices.size() ; i++ ) vertices[i] = inVertices[i];
		}
		else if( ext==std::string( "obj" ) )
		{
			struct ObjFaceIndex{ int vIndex , tIndex; };

			std::vector< Point3D< Real > > obj_vertices;
			std::vector< Point2D< Real > > obj_textures;
			std::vector< std::vector< ObjFaceIndex > > obj_faces;
			std::ifstream in( fileName );
			if( !in.is_open() ) MK_ERROR_OUT( "Could not open file for reading: " , fileName );

			std::string( line );
			unsigned int count = 0;
			while( std::getline( in , line ) )
			{
				std::stringstream ss;

				// Read vertex position
				if( line[0]=='v' && line[1]==' ' )
				{
					line = line.substr(2);
					std::stringstream ss( line );
					Point3D< Real > p;
					ss >> p[0] >> p[1] >> p[2];
					obj_vertices.push_back( p );
				}
				// Read texture coordinate
				else if( line[0]=='v' && line[1]=='t' && line[2]==' ' )
				{
					line = line.substr(3);
					std::stringstream ss( line );
					Point2D< Real > p;
					ss >> p[0] >> p[1];
					obj_textures.push_back( p );
				}
				// Read face
				else if( line[0]=='f' && line[1]==' ' )
				{
					std::vector< ObjFaceIndex > face;
					line = line.substr(1);
					while( line.size() && line[0]==' ' ) line = line.substr(1);
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
					obj_faces.push_back( face );
				}
				count++;
			}

			vertices.resize( obj_vertices.size() );
			for( int i=0 ; i<obj_vertices.size() ; i++ ) vertices[i] = obj_vertices[i];

			auto ObjIndexToArrayIndex = [&]( size_t sz , int index )
				{
					if( index>0 ) return index-1;
					else          return (int)sz + index;
				};

			// Triangulating polygonal faces
			{
				unsigned int tCount = 0;
				for( unsigned int i=0 ; i<obj_faces.size() ; i++ )
					if( obj_faces[i].size()<3 ) MK_ERROR_OUT( "Expected at least three vertices" );
					else tCount += (unsigned int)obj_faces[i].size()-2;
				if( tCount>obj_faces.size() )
				{
					MinimalAreaTriangulation< Real , 3 > mat;

					obj_faces.reserve( tCount );
					for( unsigned int i=(unsigned int)obj_faces.size() ; i!=0 ; i-- )
					{
						std::vector< ObjFaceIndex > &face = obj_faces[i-1];
						std::vector< ObjFaceIndex > oldFace = face;
						if( face.size()>3 )
						{
							std::vector< Point3D< Real > > _vertices( face.size() );
							std::vector< SimplexIndex< 2 > > triangles;
							for( unsigned int j=0 ; j<face.size() ; j++ ) _vertices[j] = vertices[ ObjIndexToArrayIndex( obj_vertices.size() , face[j].vIndex ) ];
							mat.GetTriangulation( _vertices , triangles );

							auto TriangleToOBJFace = [&]( SimplexIndex< 2 > tri )
								{
									std::vector< ObjFaceIndex > objFace(3);
									for( unsigned int i=0 ; i<3 ; i++ ) objFace[i] = oldFace[ tri[i] ];
									return objFace;
								};

							face = TriangleToOBJFace( triangles[0] );
							for( unsigned int j=1 ; j<triangles.size() ; j++ ) obj_faces.push_back( TriangleToOBJFace( triangles[j] ) );
						}
					}
				}
			}

			simplices.resize( obj_faces.size() );
			textureCoordinates.resize( 3*obj_faces.size() );
			for( int i=0 ; i<obj_faces.size() ; i++ ) for( int j=0 ; j<3 ; j++ )
			{
				if     ( obj_faces[i][j].vIndex>0 ) simplices[i][j] = obj_faces[i][j].vIndex-1;
				else if( obj_faces[i][j].vIndex<0 ) simplices[i][j] = (int)obj_vertices.size() + obj_faces[i][j].vIndex;
				else MK_ERROR_OUT( "Zero vertex index unexpected in .obj file: " , i );

				if     ( obj_faces[i][j].tIndex>0 ) textureCoordinates[3*i+j] = obj_textures[ obj_faces[i][j].tIndex-1 ];
				else if( obj_faces[i][j].tIndex<0 ) textureCoordinates[3*i+j] = obj_textures[ (int)obj_textures.size() + obj_faces[i][j].tIndex ];
				else MK_ERROR_OUT( "Zero texture index unexpected in .obj file: " , i );
			}
		}
		else MK_ERROR_OUT( "Unrecognized file type: " , fileName , " -> " , ext );
	}

	template< typename Index , typename Real >
	struct PlyTexturedFace
	{
		unsigned int nr_vertices , nr_uv_coordinates;
		Index *vertices;
		Real *uv_coordinates;
		PlyTexturedFace( void ){ vertices = nullptr , uv_coordinates = nullptr , nr_vertices = nr_uv_coordinates = 0; }
		~PlyTexturedFace( void ){ resize(0); }
		PlyTexturedFace( const PlyTexturedFace& face ) : vertices(nullptr) , uv_coordinates(nullptr) { *this = face; }
		PlyTexturedFace& operator = ( const PlyTexturedFace& face )
		{
			if( vertices ) free( vertices ) , vertices = nullptr;
			if( uv_coordinates ) free( uv_coordinates ) , uv_coordinates = nullptr;
			nr_vertices = face.nr_vertices , nr_uv_coordinates = face.nr_uv_coordinates;
			if( nr_vertices ) vertices = (Index*)malloc( sizeof(Index)*nr_vertices );
			else              vertices = nullptr;
			if( nr_uv_coordinates ) uv_coordinates = (Real*)malloc( sizeof(Real)*nr_uv_coordinates );
			else                    uv_coordinates = nullptr;
			memcpy( vertices , face.vertices , sizeof(Index)*nr_vertices );
			memcpy( uv_coordinates , face.uv_coordinates , sizeof(Real)*nr_uv_coordinates );
			return *this;
		}
		void resize( unsigned int count )
		{
			if( vertices ) free( vertices ) , vertices = nullptr;
			if( uv_coordinates ) free( uv_coordinates ) , uv_coordinates = nullptr;
			nr_vertices = nr_uv_coordinates = 0;
			if( count )
			{
				vertices = (Index*)malloc( sizeof(Index)*count ) , nr_vertices = count;
				uv_coordinates = (Real*)malloc( sizeof(Real)*2*count ) , nr_uv_coordinates = count*2;
			}
		}
		Index& operator[] ( unsigned int idx )       { return vertices[idx]; }
		Index  operator[] ( unsigned int idx ) const { return vertices[idx]; }
		Point2D< Real >  texture( unsigned int idx ) const { return Point2D< Real >( uv_coordinates[2*idx] , uv_coordinates[2*idx+1] ); }
		Point2D< Real >& texture( unsigned int idx )       { return *( (Point2D< Real >*)(uv_coordinates+2*idx) ); }
		int size( void ) const { return nr_vertices; }

		const static int NumProperties = 2;
		static GregTurk::PlyProperty Properties[];
	};

	template< typename Index , typename Real >
	GregTurk::PlyProperty PlyTexturedFace< Index , Real >::Properties[] =
	{
		{ "vertex_indices" , PLY::Type< Index >() , PLY::Type< Index >() , offsetof( PlyTexturedFace , vertices ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyTexturedFace , nr_vertices ) } ,
		{ "texcoord" , PLY::Type< Real >() , PLY::Type< Real >() , (int)offsetof( PlyTexturedFace , uv_coordinates ) , 1 , PLY_INT , PLY_INT , (int)offsetof( PlyTexturedFace , nr_uv_coordinates ) } ,
	};
}
#endif // TEXTURE_INCLUDED