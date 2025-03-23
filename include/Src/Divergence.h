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

#ifndef DIVERGENCE_INCLUDED
#define DIVERGENCE_INCLUDED

#include <Misha/Miscellany.h>
#include <Misha/MultiThreading.h>
#include "EdgeIndexing.h"

namespace MishaK
{
	class DivegenceRasterLine
	{
	public:
		int prevEdgeRowStart;
		int currEdgeRowStart;
		int nextEdgeRowStart;
		int deepCoefficientsStart;
		int texelStart;
		int texelEnd;
	};

	void InitializeDivergenceRasteLines( std::unordered_map< unsigned long long , int > &coarseEdgeIndex , const std::vector< RasterLine > &rasterLines , std::vector< DivegenceRasterLine > &divergenceRasterLines )
	{
		divergenceRasterLines.resize(rasterLines.size());
		for (int i = 0; i < rasterLines.size(); i++) {
			const RasterLine & line = rasterLines[i];
			DivegenceRasterLine & divLine = divergenceRasterLines[i];
			divLine.texelStart = line.lineStartIndex;
			divLine.texelEnd = line.lineEndIndex;
			divLine.deepCoefficientsStart = line.coeffStartIndex;
			unsigned long long prevEdgeKey = SetMeshEdgeKey(line.prevLineIndex - 1, line.prevLineIndex);
			if( coarseEdgeIndex.find(prevEdgeKey)==coarseEdgeIndex.end() ) MK_THROW( "Edge not found" );
			divLine.prevEdgeRowStart = coarseEdgeIndex[prevEdgeKey];

			unsigned long long currEdgeKey = SetMeshEdgeKey(line.lineStartIndex - 1, line.lineStartIndex);
			if( coarseEdgeIndex.find(currEdgeKey)==coarseEdgeIndex.end() ) MK_THROW( "Edge not found" );
			divLine.currEdgeRowStart = coarseEdgeIndex[currEdgeKey];

			unsigned long long nextEdgeKey = SetMeshEdgeKey(line.nextLineIndex - 1, line.nextLineIndex);
			if( coarseEdgeIndex.find(nextEdgeKey)==coarseEdgeIndex.end() ) MK_THROW( "Edge not found" );
			divLine.nextEdgeRowStart = coarseEdgeIndex[nextEdgeKey];
		}
	}

	template< class Real , class Data >
	void ComputeDivergence( const std::vector< Data > &edgeValues , std::vector< Data > &texelDivergence , const std::vector< Real > &deepDivergenceCoefficients , const SparseMatrix< Real , int > &boundaryDivergenceMatrix , const std::vector< DivegenceRasterLine > &divergenceRasterLines )
	{
		//Update Boundary Texels 
		boundaryDivergenceMatrix.Multiply(&edgeValues[0], &texelDivergence[0]);

		//Update Deep Texels

		auto UpdateRow = [&](int r)
			{
				Data* out = texelDivergence.data() + divergenceRasterLines[r].texelStart;
				const Data* previousRowEdges = edgeValues.data() + divergenceRasterLines[r].prevEdgeRowStart;
				const Data* currentRowEdges = edgeValues.data() + divergenceRasterLines[r].currEdgeRowStart;
				const Data* nextRowEdges = edgeValues.data() + divergenceRasterLines[r].nextEdgeRowStart;

				const Real * coeff = deepDivergenceCoefficients.data() + divergenceRasterLines[r].deepCoefficientsStart * 12;
				int lineLenght = divergenceRasterLines[r].texelEnd - divergenceRasterLines[r].texelStart + 1;
				for (int i = 0; i < lineLenght; coeff += 12, previousRowEdges += 2, currentRowEdges += 2, nextRowEdges += 2, i++) {
					out[i] = previousRowEdges[0] * coeff[0] + previousRowEdges[1] * coeff[1] + previousRowEdges[2] * coeff[2] + previousRowEdges[3] * coeff[3] +
						previousRowEdges[5] * coeff[4] +

						currentRowEdges[0] * coeff[5] + currentRowEdges[1] * coeff[6] + currentRowEdges[2] * coeff[7] + currentRowEdges[3] * coeff[8] +
						currentRowEdges[5] * coeff[9] +

						nextRowEdges[0] * coeff[10] +
						nextRowEdges[2] * coeff[11];
				}

				// Edge indexing
				//		 ---0---
				//		|  
				//		1  
				//		|  


				// Reduced interior texel neighbour edge indexing
				//		--0-----2-- 
				//		|    |    |
				//		1    3    4
				//		|    |    |
				//		--5----7--
				//		|    |    |
				//		6    8    9
				//		|    |    |
				//      --10---11--


			};


		ThreadPool::ParallelFor( 0 , divergenceRasterLines.size() , [&]( unsigned int , size_t r ){ UpdateRow((unsigned int)r); } );
	}
}
#endif // DIVERGENCE_INCLUDED

