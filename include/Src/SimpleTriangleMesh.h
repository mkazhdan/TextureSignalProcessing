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
#include <queue>
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
	namespace TSP
	{
		template< typename Real >
		class MeshSample
		{
		public:
			MeshSample( void ) : tID(-1){}
			MeshSample( unsigned int tID , Point2D< Real > barycentricCoords ) : tID(tID) , barycentricCoords(barycentricCoords) { }
			unsigned int tID;
			Point2D< Real > barycentricCoords;
		};

		SimplexIndex< 1 > OutgoingEdgeIndex( unsigned int k , bool flip=false );

		template< typename Real , unsigned int Dim >
		struct SimpleTriangleMesh
		{
			std::vector< Point< Real , Dim > > vertices;
			std::vector< SimplexIndex< 2 > > triangles;

			// Returns the edge-index associated with the half-edge
			SimplexIndex< 1 > edgeIndex( unsigned int he , bool flip=false ) const;

			Simplex< Real , Dim , 2 > triangle( unsigned int t ) const;

			Real area( void ) const;

			Point< Real , Dim > operator()( MeshSample< Real > s ) const;

			std::vector< unsigned int > oppositeHalfEdges( void ) const;
			std::vector< unsigned int > boundaryHalfEdges( void ) const;
			std::vector< unsigned int > trianglesToComponents( unsigned int &numComponents ) const;

			Point< Real , Dim > centroid( void ) const;
			Real boundingRadius( Point< Real , Dim > center=Point< Real , Dim >() ) const;
		};

		template< typename Real >
		struct TexturedTriangleMesh
		{
			SimpleTriangleMesh< Real , 3 > surface;
			SimpleTriangleMesh< Real , 2 > texture;

			size_t numTriangles( void ) const;

			Simplex< Real , 3 , 2 > surfaceTriangle( unsigned int t ) const;
			Simplex< Real , 2 , 2 > textureTriangle( unsigned int t ) const;

			void read( std::string meshName , bool verbose , double eps , bool flip=true );

			// Sets the boundary half-edge information
			void setBoundaryHalfEdgeInfo( std::vector< unsigned int > &textureBoundaryHalfEdges , std::vector< unsigned int > &oppositeSurfaceHalfEdges ) const;

			// Sets the boundary vertex information
			void setBoundaryVertexInfo( const std::vector< unsigned int > &textureBoundaryHalfEdges , std::map< unsigned int , unsigned int > &surfaceBoundaryVertexToIndex ) const;		
		};
#include "SimpleTriangleMesh.inl"
	}
}
#endif // SIMPLE_TRIANGLE_MESH_INCLUDED
