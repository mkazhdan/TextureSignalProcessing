/*
Copyright (c) 2019, Michael Kazhdan
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


#ifndef PLY_INCLUDED
#define PLY_INCLUDED

#define NEW_PLY

#include <vector>
#include <string>
#include <functional>
#include <limits>
#include "PlyFile.h"
#include "Geometry.h"
#include "Exceptions.h"

namespace MishaK
{
	namespace PLY
	{
		// Converts from C-type to PLY type
		template< class Real > int Type( void );	

		// Converts from C-type to PLY name
		template< typename Integer > struct Traits{ static const std::string name; };

		// A structure representing a face
		template< typename Index >
		struct Face
		{
			unsigned int nr_vertices;
			Index *vertices;

			static GregTurk::PlyProperty Properties[];
		};

		int DefaultFileType( void );

		// PLY read functionality
		int ReadHeader( std::string fileName , const GregTurk::PlyProperty *properties , int propertyNum , bool *readFlags );

		std::vector< GregTurk::PlyProperty > ReadVertexHeader( std::string fileName );
		std::vector< GregTurk::PlyProperty > ReadVertexHeader( std::string fileName , int &file_type );

		template< typename VertexFactory , typename Index >
		int Read( std::string fileName , const VertexFactory &vFactory , std::vector< typename VertexFactory::VertexType > &vertices , std::vector< std::pair< Index , Index > > *edges , std::vector< std::vector< Index > > *polygons , bool *vertexPropertiesFlag=nullptr , std::vector< std::string > *comments=nullptr );

		template< typename VertexFactory >
		int ReadVertices( std::string fileName , const VertexFactory &vFactory , std::vector< typename VertexFactory::VertexType > &vertices , bool *vertexPropertiesFlag=nullptr , std::vector< std::string > *comments=nullptr );

		template< typename VertexFactory , typename Index >
		int ReadTriangles( std::string fileName , const VertexFactory &vFactory , std::vector< typename VertexFactory::VertexType > &vertices , std::vector< SimplexIndex< 2 , Index > > &triangles , bool *vertexPropertiesFlag=nullptr , std::vector< std::string > *comments=nullptr );

		template< typename VertexFactory , typename Real , unsigned int Dim , typename Index >
		int ReadTriangles( std::string fileName , const VertexFactory &vFactory , std::vector< typename VertexFactory::VertexType > &vertices , std::vector< SimplexIndex< 2 , Index > > &triangles , std::function< Point< Real , Dim > ( typename VertexFactory::VertexType ) > VertexToPointFunctor , bool *vertexPropertiesFlag=nullptr , std::vector< std::string > *comments=nullptr );

		template< typename VertexFactory , typename Index >
		int ReadPolygons( std::string fileName , const VertexFactory &vFactory , std::vector< typename VertexFactory::VertexType > &vertices , std::vector< std::vector< Index > > &polygons ,  bool *readFlags=nullptr , std::vector< std::string > *comments=nullptr );

		template< typename VertexFactory , typename Polygon >
		int ReadPolygons( std::string fileName , const VertexFactory &vFactory , std::vector< typename VertexFactory::VertexType >& vertices , std::vector< Polygon >& polygons , GregTurk::PlyProperty *polygonProperties , int polygonPropertyNum , bool *vertexPropertiesFlag=nullptr , bool *polygonPropertiesFlag=nullptr , std::vector< std::string > *comments=nullptr );

		template< typename VertexFactory , typename Index >
		int ReadTetrahedra( std::string fileName , const VertexFactory &vFactory , std::vector< typename VertexFactory::VertexType > &vertices , std::vector< SimplexIndex< 3 , Index > > &tetrahedra , bool *vertexPropertiesFlag=nullptr , std::vector< std::string > *comments=nullptr );

		template< typename VertexFactory , unsigned int Dim , typename Index >
		int ReadSimplices( std::string fileName , const VertexFactory &vFactory , std::vector< typename VertexFactory::VertexType > &vertices , std::vector< SimplexIndex< Dim , Index > > &simplices , bool *vertexPropertiesFlag=nullptr , std::vector< std::string > *comments=nullptr );

		// PLY write functionality
		template< typename VertexFactory >
		void WriteVertices( std::string fileName , const VertexFactory &vFactory , const std::vector< typename VertexFactory::VertexType > &vertices , int file_type , const std::vector< std::string > *comments=NULL );

		template< typename VertexFactory , typename Index >
		void Write( std::string fileName , const VertexFactory &vFactory , const std::vector< typename VertexFactory::VertexType > &vertices , const std::vector< std::pair< Index , Index > > *edges , const std::vector< std::vector< Index > > *polygons , int file_type , const std::vector< std::string > *comments=NULL );

		template< typename VertexFactory , typename Index >
		void WriteTriangles( std::string fileName , const VertexFactory &vFactory , const std::vector< typename VertexFactory::VertexType > &vertices , const std::vector< SimplexIndex< 2 , Index > > &triangles , int file_type , const std::vector< std::string > *comments=NULL );

		template< typename VertexFactory , typename Index >
		void WritePolygons( std::string fileName , const VertexFactory &vFactory , const std::vector< typename VertexFactory::VertexType > &vertices , const std::vector< std::vector< Index > > &polygons , int file_type , const std::vector< std::string > *comments=NULL );

		template< class VertexFactory , typename Polygon >
		void WritePolygons( std::string fileName , const VertexFactory &vFactory , const std::vector< typename VertexFactory::VertexType > &vertices , const std::vector< Polygon > &polygons , GregTurk::PlyProperty* polygonProperties , int polygonPropertyNum , int file_type , const std::vector< std::string > *comments=NULL );

		template< typename VertexFactory , typename Index >
		void WriteTetrahedra( std::string fileName , const VertexFactory &vFactory , const std::vector< typename VertexFactory::VertexType > &vertices , const std::vector< SimplexIndex< 3 , Index > > &tetrahedra , int file_type , const std::vector< std::string > *comments=NULL );

#include "Ply.inl"
	}
}
#endif // PLY_INCLUDED
