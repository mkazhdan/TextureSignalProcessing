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
#pragma once

template< typename GeometryReal >
void AtlasMesh< GeometryReal >::initialize
(
	const TexturedTriangleMesh< GeometryReal > &inputMesh
)
{
	// Compute the mapping from triangles to charts
	{
		std::vector< unsigned int > triangleToChart = inputMesh.texture.trianglesToComponents( _numCharts );
		_triangleToChart.resize( triangleToChart.size() );
		for( unsigned int i=0 ; i<triangleToChart.size() ; i++ ) _triangleToChart[ AtlasMeshTriangleIndex(i) ] = ChartIndex( triangleToChart[i] );
	}

	// Compute the 2D mesh(es)
	{
		std::set< IndexedVector2D< GeometryReal > > IndexedPointSet;
		typename std::set< IndexedVector2D< GeometryReal > >::iterator it;

		for( unsigned int t=0 ; t<inputMesh.numTriangles() ; t++ )
		{
			SimplexIndex< 2 > tri;
			Simplex< GeometryReal , 2 , 2 > tTriangle = inputMesh.textureTriangle( t );
			for( unsigned int k=0 ; k<3 ; k++ )
			{
				ChartMeshVertexIndex currentCorner = ChartMeshVertexIndex(-1);
				IndexedVector2D< GeometryReal > idxP( tTriangle[k] , ChartMeshVertexIndex( (unsigned int)SimpleTriangleMesh< GeometryReal , 2 >::vertices.size() ) , AtlasMeshVertexIndex( inputMesh.surface.triangles[t][k] ) );
				it = IndexedPointSet.find( idxP );
				if( it==IndexedPointSet.end() )
				{
					IndexedPointSet.insert( idxP );
					_chartToAtlasVertex.push_back( AtlasMeshVertexIndex( inputMesh.surface.triangles[t][k] ) );
					currentCorner = ChartMeshVertexIndex( (unsigned int)SimpleTriangleMesh< GeometryReal , 2 >::vertices.size() );
					SimpleTriangleMesh< GeometryReal , 2 >::vertices.push_back( tTriangle[k] );
				}
				else
				{
					IndexedVector2D< GeometryReal > indexPoint = *it;
					currentCorner = indexPoint.index;
				}
				tri[k] = static_cast< unsigned int >( currentCorner );
			}
			SimpleTriangleMesh< GeometryReal , 2 >::triangles.push_back( tri );
		}
	}

	// Set the mapping from half-edges to edges
	{
		// A map from vertex to pairs to half-edge indices
		std::map< SimplexIndex< 1 > , unsigned int > edgeMap;
		for( unsigned int he=0 ; he<SimpleTriangleMesh< GeometryReal , 2 >::triangles.size()*3 ; he++ )
		{
			SimplexIndex< 1 > e = SimpleTriangleMesh< GeometryReal , 2 >::edgeIndex(he);
			if( edgeMap.find(e)==edgeMap.end() ) edgeMap[e] = he;
			else MK_THROW( "Non oriented manifold mesh" );	// If the same half-edge appears twice
		}

		unsigned int lastEdgeIndex = 0;
		ExplicitIndexVector< AtlasMeshHalfEdgeIndex , AtlasMeshEdgeIndex > halfEdgeToEdgeIndex( 3 * SimpleTriangleMesh< GeometryReal , 2 >::triangles.size() , AtlasMeshEdgeIndex(-1) );
		for( unsigned int he=0 ; he<SimpleTriangleMesh< GeometryReal , 2 >::triangles.size()*3 ; he++ )
		{
			SimplexIndex< 1 > _e = SimpleTriangleMesh< GeometryReal , 2 >::edgeIndex( he , true );

			// Set the edge associated to both halves of the edge (once)
			if( edgeMap.find(_e)!=edgeMap.end() ) // If the opposite edge exists
			{
				unsigned int _he = edgeMap[ _e ];
				if( he<_he ) halfEdgeToEdgeIndex[ AtlasMeshHalfEdgeIndex(he) ] = halfEdgeToEdgeIndex[ AtlasMeshHalfEdgeIndex(_he) ] = AtlasMeshEdgeIndex( lastEdgeIndex++ );
			}
			else halfEdgeToEdgeIndex[ AtlasMeshHalfEdgeIndex(he) ] = AtlasMeshEdgeIndex( lastEdgeIndex++ );
		}
		for( unsigned int he=0 ; he<SimpleTriangleMesh< GeometryReal , 2 >::triangles.size()*3 ; he++ ) if( halfEdgeToEdgeIndex[ AtlasMeshHalfEdgeIndex(he) ]==AtlasMeshEdgeIndex(-1) ) MK_THROW( "Non indexed half edge" );

		_halfEdgeToEdge = halfEdgeToEdgeIndex;
	}
}

template< typename GeometryReal >
void AtlasMesh< GeometryReal >::jitter( unsigned int width , unsigned int height , GeometryReal epsilon )
{
	for( int i=0 ; i<SimpleTriangleMesh< GeometryReal , 2 >::vertices.size() ; i++ ) for( int c=0 ; c<2 ; c++ )
	{
		GeometryReal dimSize = c == 0 ? (GeometryReal)width : (GeometryReal)height;
		GeometryReal scaled = SimpleTriangleMesh< GeometryReal , 2 >::vertices[i][c] * dimSize - (GeometryReal)0.5;
		GeometryReal offset = scaled - (GeometryReal)round(scaled);
		if( fabs(offset)<epsilon )
		{
			if( offset>0 ) scaled = round(scaled) + epsilon*( (GeometryReal)1. + Random< GeometryReal >() );
			else           scaled = round(scaled) - epsilon*( (GeometryReal)1. + Random< GeometryReal >() );
			SimpleTriangleMesh< GeometryReal , 2 >::vertices[i][c] = (GeometryReal)(scaled + 0.5) / dimSize;
		}
	}
}
