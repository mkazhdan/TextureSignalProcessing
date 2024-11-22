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
#ifndef CHART_DECOMPOSITION_INCLUDED
#define CHART_DECOMPOSITION_INCLUDED

#include <Misha/Miscellany.h>
#include <Src/SimpleTriangleMesh.h>
#include <queue>

void AddComponent( std::vector< int > &vertexComponent , int vIndex , int currentComponent , const std::vector< std::vector< int > > &neighbours )
{
	vertexComponent[vIndex] = currentComponent;
	std::queue< int > visitingQueue;
	visitingQueue.push( vIndex );

	while( !visitingQueue.empty() )
	{
		int currentVertex = visitingQueue.front();
		visitingQueue.pop();
		const std::vector< int > & vertexNeighbours = neighbours[currentVertex];
		for( int i=0 ; i<vertexNeighbours.size() ; i++ )
		{
			if( vertexComponent[ vertexNeighbours[i] ]==-1 )
			{
				vertexComponent[ vertexNeighbours[i] ] = currentComponent;
				visitingQueue.push( vertexNeighbours[i] );
			}
			else if( vertexComponent[ vertexNeighbours[i] ]==currentComponent ) ;
			else Miscellany::Throw( "Unexpected Condition on a connected component. Expected %d. Obtained %d.\n" , currentComponent , vertexComponent[ vertexNeighbours[i] ] );
		}
	}
}

template< typename GeometryReal >
void InitializeTriangleChartIndexing( const OrientedTexturedTriangleMesh< GeometryReal > &mesh , std::vector< int > &chartIndex , int &numCharts )
{
	std::unordered_map< unsigned long long , int > edgeIndex;
	for( int i=0 ; i<mesh.triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ )
	{
		unsigned long long  edgeKey = mesh.edgeKey(i,k);
		if( edgeIndex.find(edgeKey)==edgeIndex.end() ) edgeIndex[edgeKey] = 3*i+k;
		else Miscellany::Throw( "Non manifold mesh" );
	}

	std::vector< std::vector< int > > neighbours( mesh.triangles.size() );
	for( int i=0 ; i<mesh.triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ )
	{
		unsigned long long edgeKey = mesh.edgeKey(i,k,true);
		if( edgeIndex.find(edgeKey)!=edgeIndex.end() )
		{
			int tIndex = edgeIndex[edgeKey] / 3;
			int kIndex = edgeIndex[edgeKey] % 3;
			if( mesh.textureCoordinates[ 3*i+(k+1)%3 ][0] == mesh.textureCoordinates[ 3*tIndex+kIndex ][0] &&
				mesh.textureCoordinates[ 3*i+(k+1)%3 ][1] == mesh.textureCoordinates[ 3*tIndex+kIndex ][1] &&
				mesh.textureCoordinates[ 3*i+k ][0] == mesh.textureCoordinates[ 3*tIndex+(kIndex+1)%3 ][0] &&
				mesh.textureCoordinates[ 3*i+k ][1] == mesh.textureCoordinates[ 3*tIndex+(kIndex+1)%3 ][1] )
			{
				neighbours[i].push_back( tIndex );
			}
		}
	}
	chartIndex.clear();
	chartIndex.resize( mesh.triangles.size() , -1 );
	int currentComponent = -1;
	for( int v=0 ; v<mesh.triangles.size() ; v++ ) if( chartIndex[v]==-1 )
	{
		currentComponent++;
		AddComponent( chartIndex , v , currentComponent , neighbours );
	}
	numCharts = currentComponent + 1;
}


#endif //CHART_DECOMPOSITION_INCLUDED