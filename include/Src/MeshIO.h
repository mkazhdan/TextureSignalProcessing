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
#include "Misha/RegularGrid.h"

namespace MishaK
{
	namespace TSP
	{
		template< typename Index , typename Real , unsigned int Dim , unsigned int TDim >
		void ReadTexturedMesh( std::string fileName , std::vector< Point< Real , Dim > > &vertices , std::vector< Point< Real , TDim > > &textureCoordinates , std::vector< SimplexIndex< 2 , Index > > &simplices );

		template< typename Index , typename Real , unsigned int Dim , unsigned int TDim >
		void ReadTexturedMesh( std::string fileName , std::vector< Point< Real , Dim > > &vertices , std::vector< Point< Real , TDim > > &textureCoordinates , std::vector< SimplexIndex< 2 , Index > > &vSimplices , std::vector< SimplexIndex< 2 , Index > > &tSimplices , bool fuseTextureCoordinates=true );

		template< typename Index , typename Real , unsigned int Dim >
		void CollapseVertices( std::vector< Point< Real , Dim > > &vertices , std::vector< SimplexIndex< 2 , Index > > &simplices , double eps );

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
				catch( const Exception & ){ MK_THROW( "Failed to read ply file: " , fileName ); }

				if( !fFlags[0] ) MK_THROW( "Failed to read face indices" );
				if( !fFlags[1] ) MK_THROW( "Failed to read face textures" );
				delete[] vFlags;
				delete[] fFlags;

				if constexpr( K==2 )
				{
					size_t faceNum = inFaces.size();
					for( unsigned int i=(unsigned int)inFaces.size() ; i!=0 ; i-- )
					{
						Face &face = inFaces[i-1];
						Face oldFace = face;

						if( face.size()>3 )
						{
							std::vector< Point< Real , Dim > > _vertices( face.size() );
							std::vector< SimplexIndex< K > > _triangles;
							for( unsigned int j=0 ; j<(unsigned int)face.size() ; j++ ) _vertices[j] = inVertices[ face[j] ];
							MinimalAreaTriangulation::GetTriangulation( _vertices , _triangles );

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
					if( inFaces[i].nr_vertices!=(K+1) ) MK_THROW( "Face is not a simplex" );
					else if( inFaces[i].nr_uv_coordinates!=(K+1)*K ) MK_THROW( "Unexpected number of texture coordinates: " , inFaces[i].nr_uv_coordinates , " != " , K+1 , " * " , K );

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
				if( !in.is_open() ) MK_THROW( "Could not open file for reading: " , fileName );

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
						if( obj_faces[i].size()<3 ) MK_THROW( "Expected at least three vertices" );
						else tCount += (unsigned int)obj_faces[i].size()-2;
					if( tCount>obj_faces.size() )
					{
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
								MinimalAreaTriangulation::GetTriangulation( _vertices , triangles );

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
					else MK_THROW( "Zero vertex index unexpected in .obj file: " , i );

					if     ( obj_faces[i][j].tIndex>0 ) textureCoordinates[3*i+j] = obj_textures[ obj_faces[i][j].tIndex-1 ];
					else if( obj_faces[i][j].tIndex<0 ) textureCoordinates[3*i+j] = obj_textures[ (int)obj_textures.size() + obj_faces[i][j].tIndex ];
					else MK_THROW( "Zero texture index unexpected in .obj file: " , i );
				}
			}
			else MK_THROW( "Unrecognized file type: " , fileName , " -> " , ext );
		}

		template< typename Index , typename Real , unsigned int Dim , unsigned int TDim >
		void ReadTexturedMesh
		(
			std::string fileName ,
			std::vector< Point< Real ,  Dim > > &vertices ,
			std::vector< Point< Real , TDim > > &textures ,
			std::vector< SimplexIndex< 2 , Index > > &vSimplices ,
			std::vector< SimplexIndex< 2 , Index > > &tSimplices ,
			bool fuseTextureCoordinates 
		)
		{
			static const unsigned int K = 2;
			std::string ext = ToLower( GetFileExtension( fileName ) );
			if( ext==std::string( "ply" ) )
			{
				using VertexFactory = VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::TextureFactory< Real , TDim > >;
				using Vertex = typename VertexFactory::VertexType;
				using Face = PlyTexturedFace< unsigned int , Real >;

				VertexFactory factory;

				std::vector< Vertex > inVertices;
				std::vector< Face > inFaces;
				bool *vFlags = new bool[ factory.plyReadNum() ];
				bool *fFlags = new bool[ Face::NumProperties ];
				int file_type;

				try{ file_type = PLY::ReadPolygons( fileName , factory , inVertices , inFaces , Face::Properties , Face::NumProperties , vFlags , fFlags ); }
				catch( const Exception & ){ MK_THROW( "Failed to read ply file: " , fileName ); }

				bool perVertexTexture = factory.template plyValidReadProperties< 1 >( vFlags );
				if( !fFlags[0] ) MK_THROW( "Failed to read face indices" );
				if( !fFlags[1] && !perVertexTexture ) MK_THROW( "Failed to read face textures" );
				delete[] vFlags;
				delete[] fFlags;

				if constexpr( K==2 )
				{
					size_t faceNum = inFaces.size();
					for( unsigned int i=(unsigned int)inFaces.size() ; i!=0 ; i-- )
					{
						Face &face = inFaces[i-1];
						Face oldFace = face;

						if( face.size()>3 )
						{
							std::vector< Point< Real , Dim > > _vertices( face.size() );
							std::vector< SimplexIndex< K > > _triangles;
							for( unsigned int j=0 ; j<(unsigned int)face.size() ; j++ ) _vertices[j] = inVertices[ face[j] ].template get<0>();
							MinimalAreaTriangulation::GetTriangulation( _vertices , _triangles );

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

				vSimplices.resize( inFaces.size() );
				tSimplices.resize( inFaces.size() );
				vertices.resize( inVertices.size() );
				for( unsigned int i=0 ; i<inVertices.size() ; i++ ) vertices[i] = inVertices[i].template get<0>();

				if( perVertexTexture )
				{
					for( unsigned int i=0 ; i<inFaces.size() ; i++ ) for( unsigned int k=0 ; k<=K ; k++ ) vSimplices[i][k] = tSimplices[i][k] = inFaces[i][k];

					textures.resize( inVertices.size() );
					for( unsigned int i=0 ; i<inVertices.size() ; i++ ) textures[i] = inVertices[i].template get<1>();
				}
				else
				{
					for( unsigned int i=0 ; i<inFaces.size() ; i++ )
						if( inFaces[i].nr_vertices!=(K+1) ) MK_THROW( "Face is not a simplex" );
						else if( inFaces[i].nr_uv_coordinates!=(K+1)*K ) MK_THROW( "Unexpected number of texture coordinates: " , inFaces[i].nr_uv_coordinates , " != " , K+1 , " * " , K );

					for( unsigned int i=0 ; i<inFaces.size() ; i++ ) for( unsigned int k=0 ; k<=K ; k++ )
					{
						vSimplices[i][k] = inFaces[i][k];
						tSimplices[i][k] = i*(K+1)+k;
					}

					if( fuseTextureCoordinates )
					{
						struct IndexedTVertex
						{
							IndexedTVertex( void ) : vIndex(-1) , index(-1) {}
							IndexedTVertex( Point2D< Real > p , unsigned int vIndex , unsigned int index ) : p(p) , vIndex(vIndex) , index(index){}
							Point2D< Real > p;
							unsigned int vIndex , index;
							bool operator < ( const IndexedTVertex &v ) const
							{
								if( vIndex!=v.vIndex ) return vIndex<v.vIndex;
								if( p[0]!=v.p[0] ) return p[0]<v.p[0];
								else               return p[1]<v.p[1];
							}
							bool operator != ( const IndexedTVertex & v ) const { return vIndex!=v.vIndex || p[0]!=v.p[0] || p[1]!=v.p[1]; }
						};

						std::vector< IndexedTVertex > iVertices( inFaces.size()*3 );
						for( unsigned int i=0 ; i<inFaces.size() ; i++ ) for( unsigned int k=0 ; k<3 ; k++ )
							iVertices[3*i+k] = IndexedTVertex( inFaces[i].texture(k) , inFaces[i][k] , 3*i+k );
						std::sort( iVertices.begin() , iVertices.end() );
						std::vector< unsigned int > old2new( inFaces.size()*3 );
						unsigned int idx=0;
						for( unsigned int i=0 ; i<iVertices.size() ; i++ )
						{
							old2new[ iVertices[i].index ] = idx;
							if( i==iVertices.size()-1 || iVertices[i]!=iVertices[i+1] ) idx++;
						}
						textures.resize( idx );
						for( unsigned int i=0 , idx=0 ; i<inFaces.size() ; i++ ) for( unsigned int k=0 ; k<3 ; k++ , idx++ ) textures[ old2new[idx] ] = inFaces[i].texture(k);
						for( unsigned int i=0 , idx=0 ; i<tSimplices.size() ; i++ ) for( unsigned int k=0 ; k<3 ; k++ , idx++ ) tSimplices[i][k] = old2new[idx];
					}
					else
					{
						textures.resize( inFaces.size()*(K+1) );
						for( unsigned int i=0 ; i<inFaces.size() ; i++ ) for( unsigned int k=0 ; k<=K ; k++ ) textures[ i*(K+1)+k ] = inFaces[i].texture(k);
					}

				}
			}
			else if( ext==std::string( "obj" ) )
			{
				struct ObjFaceIndex{ int vIndex , tIndex; };

				std::vector< Point3D< Real > > obj_vertices;
				std::vector< Point2D< Real > > obj_textures;
				std::vector< std::vector< ObjFaceIndex > > obj_faces;
				std::ifstream in( fileName );
				if( !in.is_open() ) MK_THROW( "Could not open file for reading: " , fileName );

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
				textures.resize( obj_textures.size() );
				for( unsigned int i=0 ; i<obj_vertices.size() ; i++ ) vertices[i] = obj_vertices[i];
				for( unsigned int i=0 ; i<obj_textures.size() ; i++ ) textures[i] = obj_textures[i];

				auto ObjIndexToArrayIndex = [&]( size_t sz , int index )
					{
						if( index>0 ) return index-1;
						else          return (int)sz + index;
					};

				// Triangulating polygonal faces
				{
					unsigned int tCount = 0;
					for( unsigned int i=0 ; i<obj_faces.size() ; i++ )
						if( obj_faces[i].size()<3 ) MK_THROW( "Expected at least three vertices" );
						else tCount += (unsigned int)obj_faces[i].size()-2;
					if( tCount>obj_faces.size() )
					{
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
								MinimalAreaTriangulation::GetTriangulation( _vertices , triangles );

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

				vSimplices.resize( obj_faces.size() );
				tSimplices.resize( obj_faces.size() );
				for( int i=0 ; i<obj_faces.size() ; i++ ) for( int j=0 ; j<3 ; j++ )
				{
					if     ( obj_faces[i][j].vIndex>0 ) vSimplices[i][j] = obj_faces[i][j].vIndex-1;
					else if( obj_faces[i][j].vIndex<0 ) vSimplices[i][j] = (int)obj_vertices.size() + obj_faces[i][j].vIndex;
					else MK_THROW( "Zero vertex index unexpected in .obj file: " , i );

					if     ( obj_faces[i][j].tIndex>0 ) tSimplices[i][j] = obj_faces[i][j].tIndex-1;
					else if( obj_faces[i][j].tIndex<0 ) tSimplices[i][j] = (int)obj_textures.size() + obj_faces[i][j].tIndex;
					else MK_THROW( "Zero vertex index unexpected in .obj file: " , i );
				}
			}
			else MK_THROW( "Unrecognized file type: " , fileName , " -> " , ext );
		}

		template< typename Index , typename Real , unsigned int Dim >
		void CollapseVertices( std::vector< Point< Real , Dim > > &vertices , std::vector< SimplexIndex< 2 , Index > > &simplices , double eps )
		{
			auto SetMinMax = [&]( Point< Real , Dim > &min , Point< Real , Dim > &max )
				{
					min = max = Point< Real , Dim >( vertices[0] );
					for( unsigned int i=0 ; i<vertices.size() ; i++ ) for( unsigned int j=0 ; j<Dim ; j++ ) min[j] = std::min< Real >( min[j] , vertices[i][j] ) , max[j] = std::max< double >( max[j] , vertices[i][j] );
				};


			Point< Real , Dim > min , max;
			SetMinMax( min , max );

			std::vector< int > old2new( vertices.size() );
			for( unsigned int i=0 ; i<vertices.size() ; i++ ) old2new[i] = i;

			auto NewIndex = [&]( unsigned int i )
				{
					while( old2new[i]!=i ) i = old2new[i];
					return i;
				};

			// Need to choose a resolution so that:
			//	    1/res >= eps
			// <=>  res <= 1/eps
			// [NOTE] Distances are measured after rescaling to fit in the unit cube
			const long long res = static_cast< long long >( floor( 1./eps ) );
			auto VIndex = [&]( Point< Real , Dim > p )
				{
					p -= min;
					for( unsigned int j=0 ; j<Dim ; j++ ) p[j] /= max[j] - min[j];
					typename RegularGrid< Dim >::Index I;
					for( unsigned int d=0 ; d<Dim ; d++ ) I[d] = (int)( p[d] * res );
					return I;
				};

			std::map< typename RegularGrid< Dim >::Index , std::vector< unsigned int > > vMap;
			for( unsigned int i=0 ; i<vertices.size() ; i++ ) vMap[ VIndex( Point< Real , Dim >( vertices[i] ) ) ].push_back( i );

			for( unsigned int i=0 ; i<vertices.size() ; i++ ) if( old2new[i]==i )
			{
				Point< Real , Dim > p( vertices[i] );
				typename RegularGrid< Dim >::Index I = VIndex( p );
				typename RegularGrid< Dim >::Range range;
				std::vector< int > matches;
				for( unsigned int d=0 ; d<Dim ; d++ ) range.first[d] = I[d]-1 , range.second[d] = I[d]+1;


				auto Kernel = [&]( typename RegularGrid< Dim >::Index I )
					{
						auto iter = vMap.find( I );
						if( iter!=vMap.end() )
						{
							for( unsigned int j=0 ; j<iter->second.size() ; j++ )
							{
								Point< Real , Dim > q( vertices[ iter->second[j] ] );
								if( Point< Real , Dim >::Length( p - q )<eps ) matches.push_back( iter->second[j] );
							}
						}
					};
				range.process( Kernel );

				if( !matches.size() ) MK_THROW( "No matches found" );

				for( unsigned int j=0 ; j<matches.size() ; j++ ) old2new[ NewIndex( matches[j] ) ] = i;

			}
			std::map< int , int > _vMap;
			for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vMap[ NewIndex( i ) ] = 0;
			{
				int idx = 0;
				for( auto iter=_vMap.begin() ; iter!=_vMap.end() ; iter++ ) iter->second = idx++;
			}

			std::vector< Point< Real , Dim > > _vertices( _vMap.size() );
			std::vector< SimplexIndex< 2 , Index > > _simplices( simplices.size() );
			for( unsigned int i=0 ; i<vertices.size() ; i++ ) _vertices[ _vMap[ NewIndex(i) ] ] = vertices[i];
			for( unsigned int i=0 ; i<simplices.size() ; i++ ) for( unsigned int j=0 ; j<3 ; j++ ) _simplices[i][j] = _vMap[ NewIndex( simplices[i][j] ) ];

			vertices = _vertices;
			simplices = _simplices;
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
}
#endif // TEXTURE_INCLUDED