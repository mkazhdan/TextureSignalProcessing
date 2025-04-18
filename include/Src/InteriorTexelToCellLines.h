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

namespace MishaK
{
	template< typename GeometryReal , typename MatrixReal >
	void InitializeInteriorTexelToCellLines( std::vector< InteriorTexelToCellLine > &interiorTexeltoCellLine , const GridAtlas< GeometryReal , MatrixReal > &gridAtlas )
	{
		const std::vector<RasterLine> & rasterLines = gridAtlas.rasterLines;
		const std::vector<GridNodeInfo> & nodeInfo = gridAtlas.nodeInfo;
		const std::vector< GridChart< GeometryReal > > &gridCharts = gridAtlas.gridCharts;
		interiorTexeltoCellLine.resize(rasterLines.size());
		for (int i = 0; i < rasterLines.size(); i++) {
			int interiorTexelStart = rasterLines[i].lineStartIndex;
			int ci = nodeInfo[interiorTexelStart].ci;
			int cj = nodeInfo[interiorTexelStart].cj;
			int chartID = nodeInfo[interiorTexelStart].chartID;

			interiorTexeltoCellLine[i].texelStartIndex = rasterLines[i].lineStartIndex;
			interiorTexeltoCellLine[i].texelEndIndex = rasterLines[i].lineEndIndex;
			interiorTexeltoCellLine[i].coeffOffset = rasterLines[i].coeffStartIndex;

			if( gridCharts[chartID].cellType(ci-1,cj-1)!=CellType::Interior ) MK_THROW( "Non interior cell" );
			interiorTexeltoCellLine[i].previousCellStartIndex = gridCharts[chartID].cellIndices(ci - 1, cj - 1).combined + gridCharts[chartID].combinedCellOffset;

			if( gridCharts[chartID].cellType(ci-1,cj)!=CellType::Interior ) MK_THROW( "Non interior cell" );
			interiorTexeltoCellLine[i].nextCellStartIndex = gridCharts[chartID].cellIndices(ci - 1, cj).combined + gridCharts[chartID].combinedCellOffset;
		}
	}
}