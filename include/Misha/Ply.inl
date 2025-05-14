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

template<> inline int Type<               char >( void ){ return PLY_CHAR      ; }
template<> inline int Type< unsigned      char >( void ){ return PLY_UCHAR     ; }
template<> inline int Type<                int >( void ){ return PLY_INT       ; }
template<> inline int Type< unsigned       int >( void ){ return PLY_UINT      ; }
template<> inline int Type<          long long >( void ){ return PLY_LONGLONG  ; }
template<> inline int Type< unsigned long long >( void ){ return PLY_ULONGLONG ; }
template<> inline int Type<             float  >( void ){ return PLY_FLOAT     ; }
template<> inline int Type<             double >( void ){ return PLY_DOUBLE    ; }
template< class Real > inline int Type( void )
{
	MK_ERROR_OUT( "Unrecognized type" );
	return -1;
}

template<> const std::string Traits<               char >::name="char";
template<> const std::string Traits< unsigned      char >::name="unsigned char";
template<> const std::string Traits<                int >::name="int";
template<> const std::string Traits< unsigned       int >::name="unsigned int";
template<> const std::string Traits<               long >::name="long";
template<> const std::string Traits< unsigned      long >::name="unsigned long";
template<> const std::string Traits<          long long >::name="long long";
template<> const std::string Traits< unsigned long long >::name="unsigned long long";
template<> const std::string Traits<              float >::name="float";
template<> const std::string Traits<             double >::name="double";

template< typename Index >
GregTurk::PlyProperty Face< Index >::Properties[] = { GregTurk::PlyProperty( "vertex_indices" , Type< Index >() , Type< Index >() , offsetof( Face , vertices ) , 1 , PLY_INT , PLY_INT , offsetof( Face , nr_vertices ) ) };

struct Edge{ int v1 , v2; };
const GregTurk::PlyProperty EdgeProps[] =
{ 
	{ "v1" , PLY_INT , PLY_INT , (int)offsetof( Edge , v1 ) , 0 , 0 , 0 , 0 },
	{ "v2" , PLY_INT , PLY_INT , (int)offsetof( Edge , v2 ) , 0 , 0 , 0 , 0 }
};

// Read
inline int ReadHeader( std::string fileName , const GregTurk::PlyProperty *properties , int propertyNum , bool *readFlags )
{
	int file_type;
	std::vector< std::string > elist;
	float version;

	GregTurk::PlyFile *ply = GregTurk::PlyFile::Read( fileName , elist , file_type , version );
	if( !ply ) MK_THROW( "could not read ply file: " , fileName );

	for( int i=0 ; i<(int)elist.size() ; i++ ) if( elist[i]=="vertex" ) for( int j=0 ; j<propertyNum ; j++ ) if( readFlags ) readFlags[j] = ply->get_property( elist[i] , &properties[j] )!=0;

	delete ply;
	return file_type;
}

inline std::vector< GregTurk::PlyProperty > ReadVertexHeader( std::string fileName , int &file_type )
{
	float version;
	std::vector< std::string > elist;
	std::vector< GregTurk::PlyProperty > properties;

	GregTurk::PlyFile *ply = GregTurk::PlyFile::Read( fileName , elist , file_type , version );
	if( !ply ) MK_THROW( "could not read ply file: " , fileName );

	for( int i=0 ; i<elist.size() ; i++ )
	{
		std::string &elem_name = elist[i];
		size_t num_elems;
		std::vector< GregTurk::PlyProperty * > plist = ply->get_element_description( elem_name , num_elems );
		if( !num_elems ) continue;
		else if( !plist.size() )
		{
			delete ply;
			MK_THROW( "could not get element description for: " , elem_name );
		}
		if( elem_name=="vertex" ) for( unsigned int i=0 ; i<plist.size() ; i++ ) properties.push_back( *plist[i] );
	}
	delete ply;
	return properties;
}

inline std::vector< GregTurk::PlyProperty > ReadVertexHeader( std::string fileName ){ int file_type; return ReadVertexHeader( fileName , file_type ); }


template< class VertexFactory , typename Index >
int Read
(
	std::string fileName ,
	const VertexFactory &vFactory ,
	std::vector< typename VertexFactory::VertexType > &vertices , 
	std::vector< std::pair< Index , Index > > *edges ,
	std::vector< std::vector< Index > > *polygons ,
	bool *vertexPropertiesFlag ,
	std::vector< std::string > *comments
)
{
	int file_type;
	float version;
	std::vector< std::string > elist;

	GregTurk::PlyFile *ply = GregTurk::PlyFile::Read( fileName , elist , file_type , version );
	if( !ply ) MK_THROW( "could not read ply file: " , fileName );

	if( comments )
	{
		comments->reserve( comments->size() + ply->comments.size() );
		for( int i=0 ; i<(int)ply->comments.size() ; i++ ) comments->push_back( ply->comments[i] );
	}

	for( int i=0 ; i<elist.size() ; i++ )
	{
		std::string &elem_name = elist[i];
		size_t num_elems;
		std::vector< GregTurk::PlyProperty * > plist = ply->get_element_description( elem_name , num_elems );
		if( !num_elems ) continue;
		else if( !plist.size() )
		{
			delete ply;
			MK_THROW( "could not get element description for: " , elem_name );
		}
#if 1
		if( elem_name=="vertex" )
		{
			for( unsigned int i=0 ; i<vFactory.plyReadNum() ; i++)
			{
#if 1
				GregTurk::PlyProperty prop;
				if constexpr( VertexFactory::IsStaticallyAllocated() ) prop = vFactory.plyStaticReadProperty(i);
				else                                                   prop = vFactory.plyReadProperty(i);
#else
				GregTurk::PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticReadProperty(i) : vFactory.plyReadProperty(i);
#endif
				int hasProperty = ply->get_property( elem_name , &prop );
				if( vertexPropertiesFlag ) vertexPropertiesFlag[i] = (hasProperty!=0);
			}
			vertices.resize( num_elems , vFactory() );

			char *buffer = new char[ vFactory.bufferSize() ];
			for( size_t j=0 ; j<num_elems ; j++ )
			{
#if 1
				if( VertexFactory::IsStaticallyAllocated() ) ply->get_element( (void *)&vertices[j] );
#else
				if( vFactory.isStaticallyAllocated() ) ply->get_element( (void *)&vertices[j] );
#endif
				else
				{
					ply->get_element( (void *)buffer );
					vFactory.fromBuffer( buffer , vertices[j] );
				}
			}
			delete[] buffer;
		}
#else
		if( elem_name=="vertex" )
		{
			for( unsigned int i=0 ; i<vFactory.plyReadNum() ; i++)
			{
#if 1
				GregTurk::PlyProperty prop;
				if constexpr( VertexFactory::IsStaticallyAllocated() ) prop = vFactory.plyStaticReadProperty(i);
				else                                                   prop = vFactory.plyReadProperty(i);
#else
				GregTurk::PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticReadProperty(i) : vFactory.plyReadProperty(i);
#endif
				int hasProperty = ply->get_property( elem_name , &prop );
				if( vertexPropertiesFlag ) vertexPropertiesFlag[i] = (hasProperty!=0);
			}
			vertices.resize( num_elems , vFactory() );

			char *buffer = new char[ vFactory.bufferSize() ];
			for( size_t j=0 ; j<num_elems ; j++ )
			{
#if 1
				if( VertexFactory::IsStaticallyAllocated() ) ply->get_element( (void *)&vertices[j] );
#else
				if( vFactory.isStaticallyAllocated() ) ply->get_element( (void *)&vertices[j] );
#endif
				else
				{
					ply->get_element( (void *)buffer );
					vFactory.fromBuffer( buffer , vertices[j] );
				}
			}
			delete[] buffer;
		}
#endif
		else if( elem_name=="face" && polygons )
		{
			ply->get_property( elem_name , &Face< Index >::Properties[0] );
			polygons->resize( num_elems );
			for( int j=0 ; j<num_elems ; j++ )
			{
				Face< Index > ply_face;
				ply->get_element( (void *)&ply_face );
				(*polygons)[j].resize( ply_face.nr_vertices );
				for( unsigned int k=0 ; k<ply_face.nr_vertices ; k++ ) (*polygons)[j][k] = ply_face.vertices[k];
				free( ply_face.vertices );
			}  // for, read faces
		}  // if face
		else if( elem_name=="edge" && edges )
		{
			ply->get_property( elem_name , &EdgeProps[0] );
			ply->get_property( elem_name , &EdgeProps[1] );
			edges->resize( num_elems );
			for( int j=0 ; j<num_elems ; j++ )
			{
				Edge ply_edge;
				ply->get_element( (void*)&ply_edge );
				(*edges)[j].first = ply_edge.v1 , (*edges)[j].second = ply_edge.v2;
			}
		}
		else ply->get_other_element( elem_name , num_elems );

		for( int j=0 ; j<plist.size() ; j++ ) delete plist[j];
	}  // for each type of element
	delete ply;
	return file_type;
}

template< class VertexFactory >
int ReadVertices
(
	std::string fileName ,
	const VertexFactory &vFactory ,
	std::vector< typename VertexFactory::VertexType > &vertices ,
	bool* vertexPropertiesFlag ,
	std::vector< std::string > *comments
)
{
	return Read< VertexFactory , unsigned int >( fileName , vFactory , vertices , nullptr , nullptr , vertexPropertiesFlag , comments );
}

template< typename VertexFactory , typename Real , unsigned int Dim , typename Index >
int ReadTriangles
(
	std::string fileName ,
	const VertexFactory &vFactory ,
	std::vector< typename VertexFactory::VertexType > &vertices ,
	std::vector< SimplexIndex< 2 , Index > > &triangles ,
	std::function< Point< Real , Dim > ( typename VertexFactory::VertexType ) > VertexToPointFunctor ,
	bool* vertexPropertiesFlag ,
	std::vector< std::string > *comments
)
{
	std::vector< std::vector< Index > > polygons;
	int file_type = ReadPolygons( fileName , vFactory , vertices , polygons , vertexPropertiesFlag , comments );
	std::vector< Point3D< Real > > poly;
	std::vector< SimplexIndex< 2 , Index > > tris;

	triangles.clear();
	for( unsigned int i=0 ; i<polygons.size() ; i++ )
	{
		poly.resize( polygons[i].size( ) );
		for( unsigned int j=0 ; j<polygons[i].size() ; j++ ) poly[j] = VertexToPointFunctor( vertices[ polygons[i][j] ] );
#ifdef NEW_MAT_CODE
		MinimalAreaTriangulation::GetTriangulation( poly , tris );
#else // !NEW_MAT_CODE
		MinimalAreaTriangulation< Real , Dim >::GetTriangulation( poly , tris );
#endif // NEW_MAT_CODE
		for( unsigned int j=0 ; j<tris.size() ; j++ )
		{
			SimplexIndex< 2 , Index > tri;
			tri[0] = polygons[i][ tris[j][0] ];
			tri[1] = polygons[i][ tris[j][1] ];
			tri[2] = polygons[i][ tris[j][2] ];
			triangles.push_back( tri );
		}
	}
	return file_type;
}


template< typename VertexFactory , typename Index >
int ReadTriangles
(
	std::string fileName ,
	const VertexFactory &vFactory ,
	std::vector< typename VertexFactory::VertexType > &vertices ,
	std::vector< SimplexIndex< 2 , Index > > &triangles ,
	bool* vertexPropertiesFlag ,
	std::vector< std::string > *comments
)
{
	std::vector< std::vector< Index > > polygons;
	int file_type = ReadPolygons( fileName , vFactory , vertices , polygons , vertexPropertiesFlag , comments );
	triangles.resize( polygons.size() );
	for( unsigned int i=0 ; i<polygons.size() ; i++ )
	{
		if( polygons[i].size()!=3 ) MK_ERROR_OUT( "Polygon is not a triangle: " , polygons[i].size() , " != " , 3 );
		for( int j=0 ; j<3 ; j++ ) triangles[i][j] = polygons[i][j];
	}
	return file_type;
}

template< typename VertexFactory , typename Index >
int ReadPolygons
(
	std::string fileName ,
	const VertexFactory &vFactory ,
	std::vector< typename VertexFactory::VertexType > &vertices ,
	std::vector< std::vector< Index > > &polygons ,
	bool *readFlags ,
	std::vector< std::string > *comments
)
{
	std::vector< std::string > elist;
	int file_type;
	float version;
	GregTurk::PlyFile *ply = GregTurk::PlyFile::Read( fileName , elist , file_type , version );
	if( !ply ) MK_THROW( "could not read ply file: " , fileName );

	if( comments )
	{
		comments->reserve( comments->size() + ply->comments.size() );
		for( int i=0 ; i<(int)ply->comments.size() ; i++ ) comments->push_back( ply->comments[i] );
	}

	for( int i=0 ; i<(int)elist.size() ; i++ )
	{
		std::string &elem_name = elist[i];
		size_t num_elems;
		std::vector< GregTurk::PlyProperty * > plist = ply->get_element_description( elem_name , num_elems );
		if( !num_elems ) continue;
		else if( !plist.size() )
		{
			delete ply;
			MK_THROW( "could not get element description for: " , elem_name );
		}
		if( elem_name=="vertex" )
		{
			for( unsigned int i=0 ; i<vFactory.plyReadNum() ; i++)
			{
#if 1
				GregTurk::PlyProperty prop;
				if constexpr( VertexFactory::IsStaticallyAllocated() ) prop = vFactory.plyStaticReadProperty(i);
				else                                                   prop = vFactory.plyReadProperty(i);
#else
				GregTurk::PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticReadProperty(i) : vFactory.plyReadProperty(i);
#endif
				int hasProperty = ply->get_property( elem_name , &prop );
				if( readFlags ) readFlags[i] = (hasProperty!=0);
			}
			vertices.resize( num_elems , vFactory() );

			char *buffer = new char[ vFactory.bufferSize() ];
			for( size_t j=0 ; j<num_elems ; j++ )
			{
#if 1
				if( VertexFactory::IsStaticallyAllocated() ) ply->get_element( (void *)&vertices[j] );
#else
				if( vFactory.isStaticallyAllocated() ) ply->get_element( (void *)&vertices[j] );
#endif
				else
				{
					ply->get_element( (void *)buffer );
					vFactory.fromBuffer( buffer , vertices[j] );
				}
			}
			delete[] buffer;
		}
		else if( elem_name=="face" )
		{
			ply->get_property( elem_name , &Face< Index >::Properties[0] );
			polygons.resize( num_elems );
			for( unsigned int j=0 ; j<num_elems ; j++ )
			{
				Face< Index > ply_face;
				ply->get_element( (void *)&ply_face );
				polygons[j].resize( ply_face.nr_vertices );
				for( unsigned int k=0 ; k<ply_face.nr_vertices ; k++ ) polygons[j][k] = ply_face.vertices[k];
				free( ply_face.vertices );
			}  // for, read faces
		}  // if face
		else ply->get_other_element( elem_name , num_elems );

		for( int j=0 ; j<(int)plist.size() ; j++ ) delete plist[j];
	}  // for each type of element

	delete ply;

	return file_type;
}

template< typename VertexFactory , typename Polygon >
int ReadPolygons
(
	std::string fileName ,
	const VertexFactory &vFactory ,
	std::vector< typename VertexFactory::VertexType >& vertices ,
	std::vector< Polygon >& polygons ,
	GregTurk::PlyProperty *polygonProperties ,
	int polygonPropertyNum ,
	bool *vertexPropertiesFlag ,
	bool *polygonPropertiesFlag ,
	std::vector< std::string > *comments
)
{
	std::vector< std::string > elist = { std::string( "vertex" ) , std::string( "face" ) };
	int file_type;
	float version;

	GregTurk::PlyFile *ply = GregTurk::PlyFile::Read( fileName , elist , file_type , version );
	if( !ply ) MK_THROW( "could not read ply file: " , fileName );

	if( comments )
	{
		comments->reserve( comments->size() + ply->comments.size() );
		for( int i=0 ; i<ply->comments.size() ; i++ ) comments->push_back( ply->comments[i] );
	}

	for( int i=0 ; i<elist.size() ; i++ )
	{
		std::string &elem_name = elist[i];
		size_t num_elems;
		std::vector< GregTurk::PlyProperty * > plist = ply->get_element_description( elem_name , num_elems );
		if( !plist.size() )
		{
			delete ply;
			MK_THROW( "Failed to read property list: " , elem_name );
		}		
		if( elem_name=="vertex" )
		{
			for( unsigned int i=0 ; i<vFactory.plyReadNum() ; i++ )
			{
#if 1
				GregTurk::PlyProperty prop;
				if constexpr( VertexFactory::IsStaticallyAllocated() ) prop = vFactory.plyStaticReadProperty(i);
				else                                                   prop = vFactory.plyReadProperty(i);
#else
				GregTurk::PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticReadProperty(i) : vFactory.plyReadProperty(i);
#endif
				int hasProperty = ply->get_property( elem_name , &prop );
				if( vertexPropertiesFlag ) vertexPropertiesFlag[i] = (hasProperty!=0);
			}
			vertices.resize( num_elems , vFactory() );

			char *buffer = new char[ vFactory.bufferSize() ];
			for( size_t j=0 ; j<num_elems ; j++ )
			{
#if 1
				if( VertexFactory::IsStaticallyAllocated() ) ply->get_element( (void *)&vertices[j] );
#else
				if( vFactory.isStaticallyAllocated() ) ply->get_element( (void *)&vertices[j] );
#endif
				else
				{
					ply->get_element( (void *)buffer );
					vFactory.fromBuffer( buffer , vertices[j] );
				}
			}
			delete[] buffer;
		}
		else if( elem_name=="face" )
		{
			for( int i=0 ; i<polygonPropertyNum ; i++ )
			{
				int hasProperty = ply->get_property( elem_name , &polygonProperties[i] );
				if( polygonPropertiesFlag ) polygonPropertiesFlag[i] = (hasProperty!=0);
			}
			polygons.resize( num_elems );
			for( size_t j=0 ; j<num_elems ; j++ ) ply->get_element( (void *)&polygons[j] );
		}
		else ply->get_other_element( elem_name , num_elems );

		for( int j=0 ; j<plist.size() ; j++ ) delete plist[j];
	}
	delete ply;
	return file_type;
}

template< class VertexFactory , typename Index >
int ReadTetrahedra
(
	std::string fileName ,
	const VertexFactory &vFactory ,
	std::vector< typename VertexFactory::VertexType > &vertices ,
	std::vector< SimplexIndex< 3 , Index > > &tetrahedra ,
	bool* vertexPropertiesFlag ,
	std::vector< std::string > *comments
)
{
	std::vector< std::vector< Index > > polygons;
	int file_type = ReadPolygons( fileName , vFactory , vertices , polygons , vertexPropertiesFlag , comments );

	for( int i=0 ; i<polygons.size() ; i++ ) if( polygons[i].size()!=4 ) MK_ERROR_OUT( "Expected polygon with four vertices" );
	tetrahedra.resize( polygons.size() );
	for( unsigned int i=0 ; i<polygons.size() ; i++ ) for( int j=0 ; j<4 ; j++ ) tetrahedra[i][j] = polygons[i][j];
	return file_type;
}

template< class VertexFactory , unsigned int K , typename Index >
int ReadSimplices
(
	std::string fileName ,
	const VertexFactory &vFactory ,
	std::vector< typename VertexFactory::VertexType > &vertices ,
	std::vector< SimplexIndex< K , Index > > &simplexIndices ,
	bool *vertexPropertiesFlag ,
	std::vector< std::string > *comments
)
{
	std::vector< std::vector< Index > > polygons;
	int file_type = ReadPolygons( fileName , vFactory , vertices , polygons , vertexPropertiesFlag , comments );

	for( int i=0 ; i<polygons.size() ; i++ ) if( polygons[i].size()!=K+1 ) MK_THROW( "Expected polygon with " , K+1 , " vertices" );
	simplexIndices.resize( polygons.size() );
	for( unsigned int i=0 ; i<polygons.size() ; i++ ) for( int j=0 ; j<=K ; j++ ) simplexIndices[i][j] = polygons[i][j];
	return file_type;
}

template< class VertexFactory , unsigned int K , typename Index >
void WriteSimplices
(
	std::string fileName ,
	const VertexFactory &vFactory ,
	const std::vector< typename VertexFactory::VertexType > &vertices ,
	const std::vector< SimplexIndex< K , Index > > &simplexIndices ,
	int file_type ,
	std::vector< std::string > *comments=NULL
)
{
	std::vector< std::vector< Index > > polygons( simplexIndices.size() );
	for( unsigned int i=0 ; i<simplexIndices.size() ; i++ )
	{
		polygons[i].resize( K+1 );
		for( unsigned int k=0 ; k<=K ; k++ ) polygons[i][k] = simplexIndices[i][k];
	}
	WritePolygons( fileName , vFactory , vertices , polygons , file_type , comments );
}

// Write
template< typename VertexFactory , typename Index >
void Write
(
	std::string fileName ,
	const VertexFactory &vFactory ,
	const std::vector< typename VertexFactory::VertexType > &vertices , 
	const std::vector< std::pair< Index , Index > > *edges , 
	const std::vector< std::vector< Index > > *polygons,
	int file_type ,
	const std::vector< std::string > *comments
)
{
	int nr_vertices =            (int) vertices.size()    ;
	int nr_edges    = edges    ? (int)   edges->size() : 0;
	int nr_faces    = polygons ? (int)polygons->size() : 0;
	float version;
	std::vector< std::string > elist = { std::string( "vertex" ) , std::string( "edge" ) , std::string( "face" ) };

	GregTurk::PlyFile *ply = GregTurk::PlyFile::Write( fileName , elist , file_type , version );
	if( !ply ) MK_THROW( "could not write ply file: " , fileName );

	//
	// describe vertex, edge, and face properties
	//
	{
		ply->element_count( "vertex", nr_vertices );
		for( unsigned int i=0 ; i<vFactory.plyWriteNum() ; i++ )
		{
#if 1
			GregTurk::PlyProperty prop;
			if constexpr( VertexFactory::IsStaticallyAllocated() ) prop = vFactory.plyStaticWriteProperty(i);
			else                                                   prop = vFactory.plyWriteProperty(i);
#else
			GregTurk::PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticWriteProperty(i) : vFactory.plyWriteProperty(i);
#endif
			ply->describe_property( "vertex" , &prop );
		}
	}
	{
		ply->element_count( "edge" , nr_edges );
		ply->describe_property( "edge" , &EdgeProps[0] );
		ply->describe_property( "edge" , &EdgeProps[1] );
	}
	{
		ply->element_count( "face" , nr_faces );
		ply->describe_property( "face" , &Face< Index >::Properties[0] );
	}

	// Write in the comments
	if( comments ) for( int i=0 ; i<comments->size() ; i++ ) ply->put_comment( (*comments)[i] );
	ply->header_complete();

	// write vertices
	ply->put_element_setup( elist[0] );
	{
		char *buffer = new char[ vFactory.bufferSize() ];
		for( size_t j=0 ; j<(int)vertices.size() ; j++ )
		{
#if 1
			if( VertexFactory::IsStaticallyAllocated() ) ply->put_element( (void *)&vertices[j] );
#else
			if( vFactory.isStaticallyAllocated() ) ply->put_element( (void *)&vertices[j] );
#endif
			else
			{
				vFactory.toBuffer( vertices[j] , buffer );
				ply->put_element( (void *)buffer );
			}
		}
		delete[] buffer;
	}

	// write edges
	if( nr_edges )
	{
		Edge ply_edge;
		ply->put_element_setup( "edge" );
		for( int i=0 ; i<nr_edges ; i++ )
		{
			ply_edge.v1 = (*edges)[i].first , ply_edge.v2 = (*edges)[i].second;
			ply->put_element( (void*)&ply_edge );
		}
	}

	// write faces
	if( nr_faces )
	{
		Face< Index > ply_face;
		int maxFaceVerts=3;
		ply_face.nr_vertices = 3;
		ply_face.vertices = new Index[3];

		ply->put_element_setup( "face" );
		for( int i=0 ; i<nr_faces ; i++ )
		{
			int face_size = (int)(*polygons)[i].size();
			if( face_size>maxFaceVerts )
			{
				delete[] ply_face.vertices;
				maxFaceVerts = face_size;
				ply_face.vertices = new Index[face_size];
			}
			ply_face.nr_vertices = face_size;
			for( unsigned int j=0 ; j<ply_face.nr_vertices ; j++ ) ply_face.vertices[j] = (*polygons)[i][j];
			ply->put_element( (void*)&ply_face );
		}
		delete[] ply_face.vertices;
	}
	delete ply;
}

template< class VertexFactory >
void WriteVertices
(
	std::string fileName ,
	const VertexFactory &vFactory ,
	const std::vector< typename VertexFactory::VertexType > &vertices ,
	int file_type ,
	const std::vector< std::string > *comments
)
{
	int nr_vertices = int(vertices.size());
	float version;
	std::vector< std::string > elem_names = { std::string( "vertex" ) };
	GregTurk::PlyFile *ply = GregTurk::PlyFile::Write( fileName , elem_names , file_type , version );
	if( !ply ) MK_THROW( "could not write ply file: " , fileName );

	//
	// describe vertex and face properties
	//
	ply->element_count( "vertex", nr_vertices );
	for( unsigned int i=0 ; i<vFactory.plyWriteNum() ; i++ )
	{
#if 1
		GregTurk::PlyProperty prop;
		if constexpr( VertexFactory::IsStaticallyAllocated() ) prop = vFactory.plyStaticWriteProperty(i);
		else                                                   prop = vFactory.plyWriteProperty(i);
#else
		GregTurk::PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticWriteProperty(i) : vFactory.plyWriteProperty(i);
#endif
		ply->describe_property( "vertex" , &prop );
	}

	// Write in the comments
	if( comments ) for( int i=0 ; i<comments->size() ; i++ ) ply->put_comment( (*comments)[i] );
	ply->header_complete();

	// write vertices
	ply->put_element_setup( elem_names[0] );
	for( int i=0 ; i<(int)vertices.size() ; i++ ) ply->put_element( (void *)&vertices[i] );

	delete ply;
}

template< class VertexFactory , typename Index >
void WriteTriangles
(
	std::string fileName ,
	const VertexFactory &vFactory ,
	const std::vector< typename VertexFactory::VertexType > &vertices ,
	const std::vector< SimplexIndex< 2 , Index > > &triangles ,
	int file_type ,
	const std::vector< std::string > *comments
)
{
	std::vector< std::vector< Index > > polygons( triangles.size() );
	for( int i=0 ; i<triangles.size() ; i++ )
	{
		polygons[i].resize( 3 );
		for( int j=0 ; j<3 ; j++ ) polygons[i][j] = triangles[i][j];
	}
	WritePolygons( fileName , vFactory , vertices , polygons , file_type , comments );
}

template< class VertexFactory , typename Index >
void WritePolygons
(
	std::string fileName ,
	const VertexFactory &vFactory ,
	const std::vector< typename VertexFactory::VertexType > &vertices ,
	const std::vector< std::vector< Index > > &polygons ,
	int file_type ,
	const std::vector< std::string > *comments
)
{
	int nr_vertices = int(vertices.size());
	int nr_faces = int(polygons.size());
	float version;
	std::vector< std::string > elem_names = { std::string( "vertex" ) , std::string( "face" ) };
	GregTurk::PlyFile *ply = GregTurk::PlyFile::Write( fileName , elem_names , file_type , version );
	if( !ply ) MK_THROW( "could not write ply file: " , fileName );

	//
	// describe vertex and face properties
	//
	ply->element_count( "vertex", nr_vertices );
	for( unsigned int i=0 ; i<vFactory.plyWriteNum() ; i++)
	{
#if 1
		GregTurk::PlyProperty prop;
		if constexpr( VertexFactory::IsStaticallyAllocated() ) prop = vFactory.plyStaticWriteProperty(i);
		else                                                   prop = vFactory.plyWriteProperty(i);
#else
		GregTurk::PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticWriteProperty(i) : vFactory.plyWriteProperty(i);
#endif
		ply->describe_property( "vertex" , &prop );
	}
	ply->element_count( "face" , nr_faces );
	ply->describe_property( "face" , &Face< Index >::Properties[0] );

	// Write in the comments
	if( comments ) for( size_t i=0 ; i<comments->size() ; i++ ) ply->put_comment( (*comments)[i] );
	ply->header_complete();

	// write vertices
	ply->put_element_setup( elem_names[0] );
	char *buffer = new char[ vFactory.bufferSize() ];
	for( size_t j=0 ; j<(int)vertices.size() ; j++ )
	{
#if 1
		if constexpr( VertexFactory::IsStaticallyAllocated() ) ply->put_element( (void *)&vertices[j] );
#else
		if( vFactory.isStaticallyAllocated() ) ply->put_element( (void *)&vertices[j] );
#endif
		else
		{
			vFactory.toBuffer( vertices[j] , buffer );
			ply->put_element( (void *)buffer );
		}
	}
	delete[] buffer;

	// write faces
	Face< Index > ply_face;
	int maxFaceVerts = 3;
	ply_face.nr_vertices = maxFaceVerts;
	ply_face.vertices = new Index[ maxFaceVerts ];

	ply->put_element_setup( elem_names[1] );
	for( int i=0 ; i<nr_faces ; i++ )
	{
		if( (int)polygons[i].size()>maxFaceVerts )
		{
			delete[] ply_face.vertices;
			maxFaceVerts = (int)polygons[i].size();
			ply_face.vertices=new Index[ maxFaceVerts ];
		}
		ply_face.nr_vertices = (int)polygons[i].size();
		for( unsigned int j=0 ; j<ply_face.nr_vertices ; j++ ) ply_face.vertices[j] = polygons[i][j];
		ply->put_element( (void *)&ply_face );
	}
	delete[] ply_face.vertices;
	delete ply;
}

template< class VertexFactory , typename Polygon >
void WritePolygons
(
	std::string fileName ,
	const VertexFactory &vFactory ,
	const std::vector< typename VertexFactory::VertexType > &vertices ,
	const std::vector< Polygon > &polygons ,
	GregTurk::PlyProperty* polygonProperties , int polygonPropertyNum ,
	int file_type ,
	const std::vector< std::string > *comments
)
{
	int nr_vertices = int(vertices.size());
	int nr_faces = int(polygons.size());
	float version;
	std::vector< std::string > elem_names = { std::string( "vertex" ) , std::string( "face" ) };
	GregTurk::PlyFile *ply = GregTurk::PlyFile::Write( fileName , elem_names , file_type , version );
	if( !ply ) MK_THROW( "could not write ply file: " , fileName );

	//
	// describe vertex and face properties
	//
	ply->element_count( "vertex", nr_vertices );
	for( unsigned int i=0 ; i<vFactory.plyWriteNum() ; i++)
	{
#if 1
		GregTurk::PlyProperty prop;
		if constexpr( VertexFactory::IsStaticallyAllocated() ) prop = vFactory.plyStaticWriteProperty(i);
		else                                                   prop = vFactory.plyWriteProperty(i);
#else
		GregTurk::PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticWriteProperty(i) : vFactory.plyWriteProperty(i);
#endif
		ply->describe_property( "vertex" , &prop );
	}
	ply->element_count( "face" , nr_faces );
	for( unsigned int i=0 ; i<(unsigned int)polygonPropertyNum ; i++ ) ply->describe_property( "face" , polygonProperties + i );

	// Write in the comments
	if( comments ) for( size_t i=0 ; i<comments->size() ; i++ ) ply->put_comment( (*comments)[i] );
	ply->header_complete();

	// write vertices
	ply->put_element_setup( elem_names[0] );
	char *buffer = new char[ vFactory.bufferSize() ];
	for( size_t j=0 ; j<(int)vertices.size() ; j++ )
	{
#if 1
		if( VertexFactory::IsStaticallyAllocated() ) ply->put_element( (void *)&vertices[j] );
#else
		if( vFactory.isStaticallyAllocated() ) ply->put_element( (void *)&vertices[j] );
#endif
		else
		{
			vFactory.toBuffer( vertices[j] , buffer );
			ply->put_element( (void *)buffer );
		}
	}
	delete[] buffer;

	// write faces
	ply->put_element_setup( elem_names[1] );
	for( int i=0 ; i<nr_faces ; i++ ) ply->put_element( (void *)&polygons[i] );
	delete ply;
}

template< class VertexFactory >
void WritePoints
(
	std::string fileName ,
	const VertexFactory &vFactory ,
	const std::vector< typename VertexFactory::VertexType > &vertices ,
	int file_type ,
	const std::vector< std::string > *comments=NULL
)
{
	int nr_vertices = int(vertices.size());
	float version;
	std::vector< std::string > elem_names = { std::string( "vertex" ) };
	GregTurk::PlyFile *ply = GregTurk::PlyFile::Write( fileName , elem_names , file_type , version );
	if( !ply ) MK_THROW( "could not write ply file: " , fileName );

	//
	// describe vertex and face properties
	//
	ply->element_count( "vertex", nr_vertices );
	for( unsigned int i=0 ; i<vFactory.plyWriteNum() ; i++)
	{
#if 1
		GregTurk::PlyProperty prop;
		if constexpr( VertexFactory::IsStaticallyAllocated() ) prop = vFactory.plyStaticWriteProperty(i);
		else                                                   prop = vFactory.plyWriteProperty(i);
#else
		GregTurk::PlyProperty prop = vFactory.isStaticallyAllocated() ? vFactory.plyStaticWriteProperty(i) : vFactory.plyWriteProperty(i);
#endif
		ply->describe_property( "vertex" , &prop );
	}

	// Write in the comments
	if( comments ) for( size_t i=0 ; i<comments->size() ; i++ ) ply->put_comment( (*comments)[i] );
	ply->header_complete();

	// write vertices
	ply->put_element_setup( elem_names[0] );
	char *buffer = new char[ vFactory.bufferSize() ];
	for( size_t j=0 ; j<(int)vertices.size() ; j++ )
	{
#if 1
		if constexpr( VertexFactory::IsStaticallyAllocated() ) ply->put_element( (void *)&vertices[j] );
#else
		if( vFactory.isStaticallyAllocated() ) ply->put_element( (void *)&vertices[j] );
#endif
		else
		{
			vFactory.toBuffer( vertices[j] , buffer );
			ply->put_element( (void *)buffer );
		}
	}
	delete[] buffer;
	delete ply;
}

template< class VertexFactory , typename Index >
void WriteTetrahedra
(
	std::string fileName ,
	const VertexFactory &vFactory ,
	const std::vector< typename VertexFactory::VertexType > &vertices ,
	const std::vector< SimplexIndex< 3 , Index > > &tetrahedra ,
	int file_type ,
	const std::vector< std::string > *comments
)
{
	std::vector< std::vector< Index > > polygons( tetrahedra.size() );
	for( int i=0 ; i<tetrahedra.size() ; i++ )
	{
		polygons[i].resize( 4 );
		for( int j=0 ; j<4 ; j++ ) polygons[i][j] = tetrahedra[i][j];
	}
	WritePolygons( fileName , vFactory , vertices , polygons , file_type , comments );
}

inline int DefaultFileType( void ){ return PLY_ASCII; }
