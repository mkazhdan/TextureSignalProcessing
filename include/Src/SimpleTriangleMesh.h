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
#ifndef SIMPLE_TRIANGLE_MESH_INCLUDED
#define SIMPLE_TRIANGLE_MESH_INCLUDED

#include <fstream>
#include <Eigen/Sparse>
#include <Misha/Ply.h>
#include <Misha/PlyVertexData.h>
#include <Misha/Image.h>
#include <Misha/Miscellany.h>
#include <Misha/Geometry.h>
#include <Src/VectorIO.h>
#include <Src/MeshIO.h>

namespace MishaK
{
	template< typename Real > using PlyVertexFactory = VertexFactory::PositionFactory< Real , 3 >;
	template< typename Real > using PlyVertex = typename PlyVertexFactory< Real >::VertexType;

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

	template< typename GeometryReal , unsigned int Dim >
	class SimpleTriangleMesh
	{
	public:
		template< typename Real=GeometryReal >
		using Position = std::conditional_t< Dim==2 , Point2D< Real > , std::conditional_t< Dim==3 , Point3D< Real > , Point< Real , Dim > > >;
		std::vector< Position<> > vertices;
		std::vector< SimplexIndex< 2 > > triangles;

		unsigned long long edgeKey( unsigned int t , unsigned int c , bool flip=false ) const
		{
			return flip ? SetMeshEdgeKey( triangles[t][(c+1)%3] , triangles[t][c] ) : SetMeshEdgeKey( triangles[t][c] , triangles[t][(c+1)%3] );
		}

		void read( std::string fileName )
		{
			vertices.clear();
			triangles.clear();
			int file_type;
			std::vector< PlyVertex< GeometryReal > > ply_vertices;
			if ( !PlyReadTriangles( fileName.c_str() , ply_vertices , triangles , PlyVertex< GeometryReal >::ReadProperties , NULL , PlyVertex< GeometryReal >::ReadComponents, file_type ) )
				MK_THROW( "Failed to read ply file: " , fileName );
			vertices.resize( ply_vertices.size() );
			for( int i=0 ; i<ply_vertices.size() ; i++ ) vertices[i] = ply_vertices[i].point;
		}

		void write( std::string fileName ) const
		{
			std::vector< PlyVertex< GeometryReal > > ply_vertices( vertices.size() );
			for( int i=0 ; i<vertices.size() ; i++ ) ply_vertices[i].point = vertices[i];
			PlyWriteTriangles( fileName.c_str() , ply_vertices , triangles , PlyVertex< GeometryReal >::WriteProperties , PlyVertex< GeometryReal >::WriteComponents , PLY_BINARY_NATIVE );
		}

		GeometryReal triangleArea( unsigned int t ) const
		{
			Position<> d[2] = { vertices[ triangles[t][1] ] - vertices[ triangles[t][0] ] , vertices[ triangles[t][2] ] - vertices[ triangles[t][0] ] };
			SquareMatrix< GeometryReal , 2 > M;
			for( unsigned int i=0 ; i<2 ; i++ ) for( unsigned int j=0 ; j<2 ; j++ ) M(i,j) = Position<>::Dot( d[i] , d[j] );
			return (GeometryReal)( sqrt( M.determinant() ) / 2. );
		}

		GeometryReal area( void ) const
		{
			GeometryReal meshArea = 0;
			for( int t=0 ; t<triangles.size() ; t++ ) meshArea += triangleArea( t );
			return meshArea;
		}

		Position<> centroid( void ) const
		{
			GeometryReal meshArea = 0;
			Position<> centroid;
			for( int t=0 ; t<triangles.size() ; t++ )
			{
				GeometryReal tArea = triangleArea(t);
				Position<> baricenter = ( vertices[ triangles[t][0] ] + vertices[ triangles[t][1] ] + vertices[ triangles[t][2]] ) / 3;
				centroid += baricenter*tArea;
				meshArea += tArea;
			}
			return centroid/meshArea;
		}

		GeometryReal radius( const Position<> &centroid )
		{
			GeometryReal radius = 0;
			for( int v=0 ; v<vertices.size() ; v++ ) radius = std::max< GeometryReal >( radius , Position<>::Length( vertices[v] - centroid ) );
			return radius;
		}

		template< typename MatrixReal >
		void initializeMeshMatrices( Eigen::SparseMatrix< MatrixReal > &mass , Eigen::SparseMatrix< MatrixReal > &stiffness )
		{
			MatrixReal meshMass = 0;
			for( int t=0 ; t<triangles.size() ; t++ )
			{
				Position< MatrixReal > p[3] = { vertices[ triangles[t][0] ] , vertices[ triangles[t][1] ] , vertices[ triangles[t][2] ] };
				Position< MatrixReal > d[2] = { vertices[ triangles[t][1] ] - vertices[ triangles[t][0] ] , vertices[ triangles[t][2] ] - vertices[ triangles[t][0] ] };
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
				Position< MatrixReal > p[3] = { vertices[ triangles[t][0] ] , vertices[ triangles[t][1] ] , vertices[ triangles[t][2] ] };
				Position< MatrixReal > d[2] = { vertices[ triangles[t][1] ] - vertices[ triangles[t][0] ] , vertices[ triangles[t][2] ] - vertices[ triangles[t][0] ] };
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

		std::vector< int > getOppositeHalfEdges( bool &isClosed ) const
		{
			std::vector< int > oppositeHalfEdges( 3*triangles.size() , -1 );

			isClosed = true;

			std::unordered_map< unsigned long long , int > edgeIndex;
			for( int i=0 ; i<triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ )
			{
				unsigned long long edgeKey = this->edgeKey(i,k);
				if( edgeIndex.find(edgeKey)==edgeIndex.end() ) edgeIndex[edgeKey] = 3*i+k;
				else MK_THROW( "Non manifold mesh" );
			}

			for( int i=0 ; i<triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ )
			{
				// Get an iterator to the opposite edge
				auto iter = edgeIndex.find( edgeKey( i , k , true ) );
				if( iter!=edgeIndex.end() ) oppositeHalfEdges[ 3*i+k ] = iter->second;
				else isClosed = false;
			}
			return oppositeHalfEdges;
		}
	};

	template< typename GeometryReal >
	class SimpleOrientedTriangleMesh : public SimpleTriangleMesh< GeometryReal , 3 >
	{
	public:
		using SimpleTriangleMesh< GeometryReal , 3 >::vertices;
		using SimpleTriangleMesh< GeometryReal , 3 >::triangles;

		std::vector< Point3D< GeometryReal > > normals;

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

		void read( std::string fileName )
		{
			SimpleTriangleMesh< GeometryReal , 3 >::read( fileName );
			updateNormals();
		}
	};

	template< typename GeometryReal , typename ImageReal=float >
	class OrientedTexturedTriangleMesh : public SimpleOrientedTriangleMesh< GeometryReal >
	{
	protected:
		static void _Subdivide( std::vector< Point3D< GeometryReal > > &vertices , std::vector< SimplexIndex< 2 > > &triangles , std::vector< Point2D< GeometryReal > > & tCoordinates )
		{
#define EDGE_KEY( i1 , i2 ) ( (i1)>(i2) ? ( ( (long long) (i1) )<<32 ) | ( (long long) (i2) ) : ( ( (long long) (i2) )<<32 ) | ( (long long) (i1) ) )

			std::vector< SimplexIndex< 2 > > _triangles(triangles.size() * 4);
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

				_triangles[4*i+0] = SimplexIndex< 2 >( eIndex[0] , eIndex[1] , eIndex[2] );
				_tCoordinates[12*i+0] = midTCoordinates[0];
				_tCoordinates[12*i+1] = midTCoordinates[1];
				_tCoordinates[12*i+2] = midTCoordinates[2];

				_triangles[4*i+1] = SimplexIndex< 2 >( triangles[i][0] , eIndex[2] , eIndex[1] );
				_tCoordinates[12*i+3] = cornerTCoordinates[0];
				_tCoordinates[12*i+4] = midTCoordinates[2];
				_tCoordinates[12*i+5] = midTCoordinates[1];

				_triangles[4*i+2] = SimplexIndex< 2 >( eIndex[2] , triangles[i][1] , eIndex[0] );
				_tCoordinates[12*i+6] = midTCoordinates[2];
				_tCoordinates[12*i+7] = cornerTCoordinates[1];
				_tCoordinates[12*i+8] = midTCoordinates[0];

				_triangles[4*i+3] = SimplexIndex< 2 >( eIndex[1] , eIndex[0] , triangles[i][2] );
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
		using SimpleOrientedTriangleMesh< GeometryReal >::vertices;
		using SimpleOrientedTriangleMesh< GeometryReal >::triangles;
#ifdef USE_TEXTURE_TRIANGLES
		SimplexTriangleMesh< GeometryReal , 2 > textureMesh;
#else // !USE_TEXTURE_TRIANGLES
		std::vector< Point2D< GeometryReal > > textureCoordinates;
#endif // USE_TEXTURE_TRIANGLES

		void write( std::string fileName , const char *atlasName=NULL ) const
		{
			std::vector< PlyTexturedFace< unsigned int , GeometryReal > > plyTexturedFaces( triangles.size() );
			for( int i=0 ; i<triangles.size() ; i++ )
			{
				plyTexturedFaces[i].resize(3);
				for( int j=0 ; j<3 ; j++ )
				{
					plyTexturedFaces[i][j] = triangles[i][j];
#ifdef USE_TEXTURE_TRIANGLES
					plyTexturedFaces[i].texture(j) = textureMesh.vertices[ textureTriangles[i][j] ];
#else // !USE_TEXTURE_TRIANGLES
					plyTexturedFaces[i].texture(j) = textureCoordinates[3*i+j];
#endif // USE_TEXTURE_TRIANGLES
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
				PlyWritePolygons( fileName.c_str() , vertices , plyTexturedFaces , PlyVertex< GeometryReal >::WriteProperties , PlyVertex< GeometryReal >::WriteComponents , PlyTexturedFace< unsigned int , GeometryReal >::WriteProperties , PlyTexturedFace< unsigned int , GeometryReal >::WriteComponents , PLY_BINARY_NATIVE , comments , 1 );
				delete[] comments;
			}
			else PlyWritePolygons( fileName.c_str() , vertices , plyTexturedFaces , PlyVertex< GeometryReal >::WriteProperties , PlyVertex< GeometryReal >::WriteComponents , PlyTexturedFace< unsigned int , GeometryReal >::WriteProperties , PlyTexturedFace< unsigned int , GeometryReal >::WriteComponents , PLY_BINARY_NATIVE );
		}

		void read( std::string meshName , bool verbose , double eps )
		{
			vertices.clear();
			triangles.clear();
			textureCoordinates.clear();
			ReadTexturedMesh( meshName , vertices , textureCoordinates , triangles );
			if( eps>0 ) CollapseVertices( vertices , triangles , eps );

			// Flip the vertical axis
#ifdef USE_TEXTURE_TRIANGLES
			for( int i=0 ; i<textureMesh.vertices.size() ; i++ ) textureMesh.vertices[i][1] = (GeometryReal)1. - textureMesh.vertices[i][1];
#else // !USE_TEXTURE_TRIANGLES
			for( int i=0 ; i<textureCoordinates.size() ; i++ ) textureCoordinates[i][1] = (GeometryReal)1. - textureCoordinates[i][1];
#endif // USE_TEXTURE_TRIANGLES

			SimpleOrientedTriangleMesh< GeometryReal >::updateNormals();

			for( unsigned int i=0 ; i<triangles.size() ; i++ )
			{
				Point3D< GeometryReal > v[3];
				for( unsigned int j=0 ; j<3 ; j++ ) for( unsigned int k=0 ; k<2 ; k++ ) v[j][k] = textureCoordinates[i*3+j][k];
				GeometryReal area2 = Point3D< GeometryReal >::SquareNorm( Point3D< GeometryReal >::CrossProduct( v[1]-v[0] , v[2]-v[0] ) );
				if( !area2 ) MK_WARN( "Zero area triangle: " , i );
			}
		}

		void initializeBoundaryEdges( std::vector< int > &boundaryEdges ) const
		{
			bool isClosedMesh = true;

			std::unordered_map< unsigned long long , int > edgeIndex;
			for( int i=0 ; i<triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ )
			{
				unsigned long long edgeKey = this->edgeKey(i,k);
				if( edgeIndex.find(edgeKey)==edgeIndex.end() ) edgeIndex[edgeKey] = 3*i+k;
				else MK_THROW( "Non manifold mesh" );
			}

			for( int i=0 ; i<triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ )
			{
				int currentEdgeIndex = 3*i+k;
				unsigned long long edgeKey = this->edgeKey(i,k,true);
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
			SimpleOrientedTriangleMesh< GeometryReal >::updateNormals();
		}
	};
}
#endif // SIMPLE_TRIANGLE_MESH_INCLUDED
