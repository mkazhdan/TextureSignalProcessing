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
#ifndef MAPLOOP_INCLUDED
#define MAPLOOP_INCLUDED

#include <unordered_map>
#include <unordered_set>
#include <Misha/Miscellany.h>


void LoopVertices( std::unordered_map< unsigned long long , unsigned long long > &forwardMap , std::vector< std::vector< unsigned long long > > &loopVertices )
{
	loopVertices.clear();
	std::unordered_set<unsigned long long> alreadyVisitedVertex;
	unsigned int loopCounter = 0;
	unsigned int edgeCounter = 0;
	for (auto iter = forwardMap.begin(); iter != forwardMap.end(); iter++){
		unsigned long long sourceVertex = (*iter).first;
		if (alreadyVisitedVertex.find(sourceVertex) == alreadyVisitedVertex.end()){
			std::vector< unsigned long long> currentLoop;
			unsigned long long currentVertex = sourceVertex;
			bool terminate = false;
			unsigned int startEdgeCounter = edgeCounter;
			do{
				if (alreadyVisitedVertex.find(currentVertex) == alreadyVisitedVertex.end()){
					alreadyVisitedVertex.insert(currentVertex);
					currentLoop.push_back(currentVertex);
					auto mappedVertex = forwardMap.find(currentVertex);
					if (mappedVertex != forwardMap.end()){
						unsigned long long nextVertex = (*mappedVertex).second;
						edgeCounter++;
						currentVertex = nextVertex;
					}
					else Miscellany::Throw( "Vertex to dead end" );
				}
				else
				{
					if( currentVertex!=sourceVertex ) Miscellany::Throw( "Non-simple loop, node %llu" , currentVertex );
					terminate = true;
				}
			} while (!terminate);
			loopCounter++;
			loopVertices.push_back(currentLoop);
		}
	}
}

void ListVerticesSimpleLoop( std::unordered_map< int , int > &forwardMap , std::vector< int > &vertexList )
{
	vertexList.clear();
	std::unordered_set<int> alreadyVisitedVertex;
	int loopCounter = 0;

	int sourceVertex = (*forwardMap.begin()).first;
	int currentVertex = sourceVertex;
	bool terminate = false;
	do
	{
		if (alreadyVisitedVertex.find(currentVertex) == alreadyVisitedVertex.end()) {
			alreadyVisitedVertex.insert(currentVertex);
			vertexList.push_back(currentVertex);
			auto mappedVertex = forwardMap.find(currentVertex);
			if (mappedVertex != forwardMap.end()) {
				unsigned long long nextVertex = (*mappedVertex).second;
				currentVertex = (int)nextVertex;
			}
			else Miscellany::Throw( "Vertex to dead end" );
		}
		else
		{
			if( currentVertex!=sourceVertex ) Miscellany::Throw( "Non-simple loop" );
			terminate = true;
		}
	}
	while (!terminate);
}

#endif //MAPLOOP_INCLUDED