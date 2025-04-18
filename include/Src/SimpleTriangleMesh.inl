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

////////////////////////
// SimpleTriangleMesh //
////////////////////////

inline EdgeIndex CornerEdgeIndex( unsigned int k , bool flip )
{
#ifdef NEW_HALF_EDGE
	if( flip ) return EdgeIndex( (k+2)%3 , (k+1)%3 );
	else       return EdgeIndex( (k+1)%3 , (k+2)%3 );
#else // !NEW_HALF_EDGE
	if( flip ) return EdgeIndex( (k+1)%3 , k );
	else       return EdgeIndex( k , (k+1)%3 );
#endif // NEW_HALF_EDGE
}

template< typename Real , unsigned int Dim >
EdgeIndex SimpleTriangleMesh< Real , Dim >::edgeIndex( unsigned int he , bool flip ) const
{
	unsigned int t = he/3 , k = he%3;
	EdgeIndex eIndex = CornerEdgeIndex( k , flip );
	return EdgeIndex( triangles[t][ eIndex[0] ] , triangles[t][ eIndex[1] ] );
}

template< typename Real , unsigned int Dim >
void SimpleTriangleMesh< Real , Dim >::read( std::string fileName )
{
	vertices.clear();
	triangles.clear();
	int file_type;
	std::vector< PlyVertex > ply_vertices;
	if ( !PlyReadTriangles( fileName.c_str() , ply_vertices , triangles , PlyVertex::ReadProperties , NULL , PlyVertex::ReadComponents, file_type ) )
		MK_THROW( "Failed to read ply file: " , fileName );
	vertices.resize( ply_vertices.size() );
	for( int i=0 ; i<ply_vertices.size() ; i++ ) vertices[i] = ply_vertices[i].point;
}

template< typename Real , unsigned int Dim >
void SimpleTriangleMesh< Real , Dim >::write( std::string fileName ) const
{
	std::vector< PlyVertex > ply_vertices( vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) ply_vertices[i].point = vertices[i];
	PlyWriteTriangles( fileName.c_str() , ply_vertices , triangles , PlyVertex::WriteProperties , PlyVertex::WriteComponents , PLY_BINARY_NATIVE );
}

template< typename Real , unsigned int Dim >
Simplex< Real , Dim , 2 > SimpleTriangleMesh< Real , Dim >::triangle( unsigned int t ) const
{
	Simplex< Real , Dim , 2 > s;
	for( unsigned int k=0 ; k<3 ; k++ ) s[k] = vertices[ triangles[t][k] ];
	return s;
}
template< typename Real , unsigned int Dim >
Real SimpleTriangleMesh< Real , Dim >::area( void ) const
{
	Real meshArea = 0;
	for( int t=0 ; t<triangles.size() ; t++ ) meshArea += triangle(t).measure();
	return meshArea;
}

template< typename Real , unsigned int Dim >
Point< Real , Dim > SimpleTriangleMesh< Real , Dim >::operator()( MeshSample< Real > s ) const
{
	return triangle( s.tID )( s.barycentricCoords );
}

template< typename Real , unsigned int Dim >
std::vector< unsigned int > SimpleTriangleMesh< Real , Dim >::oppositeHalfEdges( void ) const
{
	std::vector< unsigned int > oppositeHalfEdges( 3*triangles.size() , static_cast< unsigned int >(-1) );

	std::map< EdgeIndex , unsigned int > edgeMap;
	for( unsigned int he=0 ; he<triangles.size()*3 ; he++ )
	{
		EdgeIndex e = edgeIndex( he ); 
		if( edgeMap.find(e)==edgeMap.end() ) edgeMap[e] = he;
		else MK_THROW( "Non manifold mesh" );
	}

	for( unsigned int he=0 ; he<triangles.size()*3 ; he++ )
	{
		// Get an iterator to the opposite edge
		auto iter = edgeMap.find( edgeIndex( he , true ) );
		if( iter!=edgeMap.end() ) oppositeHalfEdges[he] = iter->second;
	}

	return oppositeHalfEdges;
}

template< typename Real , unsigned int Dim >
std::vector< unsigned int > SimpleTriangleMesh< Real , Dim >::boundaryHalfEdges( void ) const
{
	std::vector< unsigned int > boundaryHalfEdges;

	std::vector< unsigned int > oppositeHalfEdges = this->oppositeHalfEdges();
	for( unsigned int he=0 ; he<oppositeHalfEdges.size() ; he++ ) if( oppositeHalfEdges[he]==-1 ) boundaryHalfEdges.push_back( he );
	return boundaryHalfEdges;
}

template< typename Real , unsigned int Dim >
Point< Real , Dim > SimpleTriangleMesh< Real , Dim >::centroid( void ) const
{
	Point< Real , Dim > center;
	Real area = 0;
	for( unsigned int i=0 ; i<triangles.size() ; i++ )
	{
		Simplex< Real , 3 , 2 > s = triangle(i);
		Real a = s.measure();
		center += s.center() * a;
		area += a;
	}
	return center / area;
}

template< typename Real , unsigned int Dim >
Real SimpleTriangleMesh< Real , Dim >::boundingRadius( Point< Real , Dim > center ) const
{
	Real l = 0;
	for( unsigned int i=0 ; i<vertices.size() ; i++ ) l = std::max< Real >( l , Point< Real , Dim >::Length( vertices[i]-center ) );
	return l;
}


template< typename Real , unsigned int Dim >
std::vector< unsigned int > SimpleTriangleMesh< Real , Dim >::trianglesToComponents( unsigned int &numComponents ) const
{
	std::vector< std::vector< unsigned int > > neighbours( triangles.size() );
	{
		std::vector< unsigned int > oppositeHalfEdges = this->oppositeHalfEdges();
		for( unsigned int he=0 ; he<triangles.size()*3 ; he++ )
			if( oppositeHalfEdges[he]!=-1 ) neighbours[he/3].push_back( oppositeHalfEdges[he]/3 );
	}

	std::vector< unsigned int > components( triangles.size() , static_cast< unsigned int >(-1) );

	auto AddComponent = [&]( unsigned int t , unsigned int c )
		{
			components[t] = c;
			std::queue< unsigned int > visitingQueue;
			visitingQueue.push(t);

			while( !visitingQueue.empty() )
			{
				unsigned int currentVertex = visitingQueue.front();
				visitingQueue.pop();
				const std::vector< unsigned int > & _neighbors = neighbours[ currentVertex ];
				for( unsigned int i=0 ; i<_neighbors.size() ; i++ )
				{
					if( components[ _neighbors[i] ]==-1 )
					{
						components[ _neighbors[i] ] = c;
						visitingQueue.push( _neighbors[i] );
					}
					else if( components[ _neighbors[i] ]==c ) ;
					else MK_THROW( "Unexpected Condition on a connected component. Expected " , c , ". Obtained " , components[ _neighbors[i] ] , "." );
				}
			}
		};

	numComponents = 0;
	for( unsigned int t=0 ; t<triangles.size() ; t++ ) if( components[t]==-1 ) AddComponent( t , numComponents++ );
	return components;
}

//////////////////////////
// TexturedTriangleMesh //
//////////////////////////
template< typename Real >
size_t TexturedTriangleMesh< Real >::numTriangles( void ) const { return surface.triangles.size(); }

template< typename Real >
void TexturedTriangleMesh< Real >::write( std::string fileName ) const
{
	std::vector< PlyTexturedFace< unsigned int , Real > > plyTexturedFaces( surface.triangles.size() );
	for( int i=0 ; i<surface.triangles.size() ; i++ )
	{
		plyTexturedFaces[i].resize(3);
		for( int j=0 ; j<3 ; j++ )
		{
			plyTexturedFaces[i][j] = surface.triangles[i][j];
			plyTexturedFaces[i].texture(j) = texture.vertices[ texture.triangles[i][j] ];
		}
	}

	std::vector< PlyVertex > vertices( surface.vertices.size() );
	for( int i=0 ; i<vertices.size() ; i++ ) vertices[i] = surface.vertices[i];

	PlyWritePolygons( fileName.c_str() , vertices , plyTexturedFaces , PlyVertex::WriteProperties , PlyVertex::WriteComponents , PlyTexturedFace< unsigned int , Real >::WriteProperties , PlyTexturedFace< unsigned int , Real >::WriteComponents , PLY_BINARY_NATIVE );
}

template< typename Real >
void TexturedTriangleMesh< Real >::read( std::string meshName , bool verbose , double eps )
{
	surface.vertices.clear();
	texture.vertices.clear();
	surface.triangles.clear();
	texture.triangles.clear();
	ReadTexturedMesh( meshName , surface.vertices , texture.vertices , surface.triangles , texture.triangles );
	if( surface.triangles.size()!=texture.triangles.size() ) MK_ERROR_OUT( "Triangle counts differ: " , surface.triangles.size() , " != " , texture.triangles.size() );

	if( eps>0 ) CollapseVertices( surface.vertices , surface.triangles , eps );

	// Flip the vertical axis
	for( int i=0 ; i<texture.vertices.size() ; i++ ) texture.vertices[i][1] = (Real)1. - texture.vertices[i][1];
	for( unsigned int i=0 ; i<texture.triangles.size() ; i++ ) if( texture.triangle(i).measure()==0 ) MK_WARN( "Zero area texture triangle: " , i );
}

template< typename Real >
Simplex< Real , 3 , 2 > TexturedTriangleMesh< Real >::surfaceTriangle( unsigned int t ) const { return surface.triangle(t); }

template< typename Real >
Simplex< Real , 2 , 2 > TexturedTriangleMesh< Real >::textureTriangle( unsigned int t ) const { return texture.triangle(t); }

template< typename Real >
void TexturedTriangleMesh< Real >::setBoundaryHalfEdgeInfo
(
	std::vector< unsigned int > &textureBoundaryHalfEdges ,
	std::vector< unsigned int > &oppositeSurfaceHalfEdges
)
const
{
	oppositeSurfaceHalfEdges = surface.oppositeHalfEdges();

	for( int c=0 ; c<numTriangles()*3 ; c++ )
		if( oppositeSurfaceHalfEdges[c]!=-1 )
		{
			unsigned int _c = oppositeSurfaceHalfEdges[c];
			EdgeIndex  eIndex = texture.edgeIndex(  c );
			EdgeIndex _eIndex = texture.edgeIndex( _c );
			if( eIndex[0]!=_eIndex[1] || eIndex[1]!=_eIndex[0] ) textureBoundaryHalfEdges.push_back( c );
		}
		else textureBoundaryHalfEdges.push_back( c );
}

template< typename Real >
unsigned int TexturedTriangleMesh< Real >::setBoundaryVertexInfo
(
	const std::vector< unsigned int > &textureBoundaryHalfEdges ,
	std::unordered_map< unsigned int , unsigned int > &textureBoundarySurfaceVertexIndices
)
const
{
	unsigned int idx = 0;

	// Iterate over the boundary half-edges, get the associated edges, and add their (surface) vertices to the map if not already there.
	for( int b=0 ; b<textureBoundaryHalfEdges.size() ; b++ )
	{
		EdgeIndex e = surface.edgeIndex( textureBoundaryHalfEdges[b] );
		if( textureBoundarySurfaceVertexIndices.find( e[0] )==textureBoundarySurfaceVertexIndices.end() ) textureBoundarySurfaceVertexIndices[ e[0] ] = idx++;
		if( textureBoundarySurfaceVertexIndices.find( e[1] )==textureBoundarySurfaceVertexIndices.end() ) textureBoundarySurfaceVertexIndices[ e[1] ] = idx++;
	}
	return idx;
}