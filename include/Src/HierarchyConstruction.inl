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

#include <Misha/Miscellany.h>
#include <Misha/Exceptions.h>

template< typename GeometryReal >
void InitializeBoundaryAndDeepTexelIndexing( const std::vector< GridChart< GeometryReal > > &gridCharts , const int numTexels , std::vector< int > &boundaryAndDeepIndex , std::vector< int > &boundaryGlobalIndex , std::vector< int > &deepGlobalIndex ) 
{
	boundaryAndDeepIndex.resize(numTexels, 0);
	int lastGlobalIndex = 0;
	int lastBoundaryIndex = 1;
	int lastDeepIndex = -1;

	for (int c = 0; c < gridCharts.size(); c++)
	{
		const GridChart< GeometryReal > &gridChart = gridCharts[c];
		for (int j = 0; j < gridChart.nodeType.height(); j++)for (int i = 0; i < gridChart.nodeType.width(); i++) {
			if (gridChart.nodeType(i, j) == 0 || gridChart.nodeType(i, j) == 1) {
				boundaryGlobalIndex.push_back(lastGlobalIndex);
				boundaryAndDeepIndex[lastGlobalIndex] = lastBoundaryIndex;
				lastBoundaryIndex++;
				if( gridChart.globalTexelIndex(i,j)!=lastGlobalIndex ) THROW( "Unexpected global index: actual " , gridChart.globalTexelIndex(i,j) , " , expected " , lastGlobalIndex );
				lastGlobalIndex++;
			}
			else if (gridChart.nodeType(i, j) == 2) {
				deepGlobalIndex.push_back(lastGlobalIndex);
				boundaryAndDeepIndex[lastGlobalIndex] = lastDeepIndex;
				lastDeepIndex--;
				if( gridChart.globalTexelIndex(i,j)!=lastGlobalIndex ) THROW( "Unexpected global index: actual " , gridChart.globalTexelIndex(i,j) , " , expected " , lastGlobalIndex );
				lastGlobalIndex++;
			}
		}
	}
}

//Node type : inactive(-1) , exterior (0), interior boundary (1), interior deep (2) hybryd (both deep and boundary for the solver)(3).
//Cell type : inactive(-1) , boundary (0), interior (1).
template< typename GeometryReal >
void InitializeGridChartsActiveNodes( const int chartID, const AtlasChart< GeometryReal > &atlasChart , GridChart< GeometryReal > &gridChart , std::vector< GridNodeInfo > &nodeInfo , std::vector< RasterLine > &rasterLines , std::vector< SegmentedRasterLine > &segmentedLines , std::vector< ThreadTask > &threadTasks , int &lastGlobalTexelIndex , int &lastGlobalTexelInteriorIndex , int &lastGlobalTexelDeepIndex , int &lastGlobalTexelBoundaryIndex , int &lastGlobalCellIndex , int &lastGlobalBoundaryCellIndex , int &lastGlobalInteriorCellIndex , const MultigridBlockInfo &multigridBlockInfo )
{
	int width = gridChart.width;
	int height = gridChart.height;
	GeometryReal cellSizeW = gridChart.cellSizeW;
	GeometryReal cellSizeH = gridChart.cellSizeH;
	Image< int > &nodeType = gridChart.nodeType;
	nodeType.resize( width , height );
	for( int i=0 ; i<nodeType.size() ; i++ ) nodeType[i] = -1;
#ifdef DEBUG_ATLAS
	Image< int > nodeOwner;
	nodeOwner.resize( width , height );
#endif // DEBUG_ATLAS

	Image< int > &cellType = gridChart.cellType;
	cellType.resize( width-1 , height-1 );
	for( int i=0 ; i<cellType.size() ; i++ ) cellType[i] = -1;

	Image< int > &triangleID = gridChart.triangleID;
	triangleID.resize( width , height );
	for ( int i=0 ; i<triangleID.size() ; i++ ) triangleID[i] = -1;

	Image< Point2D< GeometryReal > > &barycentricCoords = gridChart.barycentricCoords;
	barycentricCoords.resize(width, height);

	//(1) Add interior texels
	for( int t=0 ; t<atlasChart.triangles.size() ; t++ )
	{
		Point2D< GeometryReal > tPos[3];
		for( int i=0 ; i<3 ; i++ ) tPos[i] = atlasChart.vertices[ atlasChart.triangles[t][i] ] - gridChart.corner;
		int minCorner[2] , maxCorner[2];
		GetTriangleIntegerBBox( tPos , (GeometryReal)1./cellSizeW , (GeometryReal)1./cellSizeH , minCorner , maxCorner );

		SquareMatrix< GeometryReal , 2 > barycentricMap = GetBarycentricMap(tPos);

		for( int j=minCorner[1] ; j<=maxCorner[1] ; j++ ) for( int i=minCorner[0] ; i<=maxCorner[0] ; i++ )
		{
			Point2D< GeometryReal > texel_pos = Point2D< GeometryReal >( (GeometryReal)i*cellSizeW , (GeometryReal)j*cellSizeH ) - tPos[0];
			Point2D< GeometryReal > barycentricCoord = barycentricMap*texel_pos;
			if( barycentricCoord[0]>=0 && barycentricCoord[1]>=0 && ( barycentricCoord[0]+barycentricCoord[1] )<=1 )
			{
#ifdef DEBUG_ATLAS
				if( nodeType(i,j)!=-1 ) WARN( "Node ( " , i , " , " , j , " ) covered by two triangles: " , atlasChart.triangles[t]() , " " , nodeOwner(i,j) );
#else // !DEBUG_ATLAS
				if( nodeType(i,j)!=-1 ) WARN( "Node ( " , i , " , " , j , " ) in chart " , chartID , " already covered [" , t , "]" );
#endif // DEBUG_ATLAS
				nodeType(i,j) = 1;
#ifdef DEBUG_ATLAS
				nodeOwner(i,j) = atlasChart.triangles[t]();
#endif // DEBUG_ATLAS
				triangleID(i,j) = atlasChart.meshTriangleIndices[t];
				barycentricCoords(i,j) = barycentricCoord;
			}
		}
	}

	//(2) Add texels adjacent to boundary cells
	int interiorCellTriangles = 0;

	for( int e=0 ; e<atlasChart.boundaryHalfEdges.size() ; e++ )
	{
		int tIndex = atlasChart.boundaryHalfEdges[e] / 3;
		int kIndex = atlasChart.boundaryHalfEdges[e] % 3;

		Point2D< GeometryReal > ePos[2];
		ePos[0] = atlasChart.vertices[atlasChart.triangles[tIndex][kIndex]] - gridChart.corner;
		ePos[1] = atlasChart.vertices[atlasChart.triangles[tIndex][(kIndex + 1) % 3]] - gridChart.corner;

		int minCorner[2] , maxCorner[2];
		GetEdgeIntegerBBox( ePos , (GeometryReal)1./cellSizeW , (GeometryReal)1./cellSizeH , minCorner , maxCorner);
		Point2D< GeometryReal > tPos[3];
		for (int k = 0; k < 3; k++) tPos[k] = atlasChart.vertices[atlasChart.triangles[tIndex][k]] - gridChart.corner;

		SquareMatrix< GeometryReal , 2 > barycentricMap = GetBarycentricMap(tPos);

		Point2D< GeometryReal > edgeNormal;
		GeometryReal edgeLevel;
		Point2D< GeometryReal > edgeDirection = ePos[1] - ePos[0];
		edgeNormal = Point2D< GeometryReal >( edgeDirection[1] , -edgeDirection[0] );
		edgeNormal /= Point2D< GeometryReal >::Length( edgeNormal );
		edgeLevel = ( Point2D< GeometryReal >::Dot( edgeNormal , ePos[0] ) + Point2D< GeometryReal >::Dot( edgeNormal , ePos[1] ) ) / 2;

		//(2.1) Add texels adjacent to cell intersecting boundary edges

		for( int c=0 ; c<2 ; c++ ) for( int j=minCorner[1] ; j<=maxCorner[1] ; j++ ) for( int i=minCorner[0] ; i<=maxCorner[0] ; i++ )
		{
			Point2D< GeometryReal > cellNode[2] = { Point2D< GeometryReal >( (GeometryReal)i*cellSizeW , (GeometryReal)j*cellSizeH ), Point2D< GeometryReal >( (GeometryReal)i*cellSizeW , (GeometryReal)j*cellSizeH ) };
			if( c==0 ) cellNode[1][c] += cellSizeW;
			else cellNode[1][c] += cellSizeH;
			Point2D< GeometryReal > cellSide = cellNode[1] - cellNode[0];
			Point2D< GeometryReal > cellSideNormal = Point2D< GeometryReal >(cellSide[1], -cellSide[0]);
			cellSideNormal /= Point2D< GeometryReal >::Length(cellSideNormal);
			GeometryReal cellLevel = ( Point2D< GeometryReal >::Dot( cellSideNormal , cellNode[0] ) + Point2D< GeometryReal >::Dot( cellSideNormal , cellNode[1] ) ) / 2;

			bool oppositeEdgeSide = ( Point2D< GeometryReal >::Dot( edgeNormal , cellNode[0] ) - edgeLevel ) * ( Point2D< GeometryReal >::Dot( edgeNormal , cellNode[1] ) - edgeLevel )<0;
			bool oppositeCellSide = ( Point2D< GeometryReal >::Dot( cellSideNormal , ePos[1] ) - cellLevel ) * ( Point2D< GeometryReal >::Dot( cellSideNormal , ePos[0] ) - cellLevel )<0;

			if( oppositeEdgeSide && oppositeCellSide )
			{
				if( c==0 ) cellType(i,j-1) = cellType(i,j) = 0;
				if( c==1 ) cellType(i-1,j) = cellType(i,j) = 0;

				for( int dn=-1 ; dn<=1 ; dn++ ) for( int dc=0 ; dc<2 ; dc++ ) //Update nodes on adjacent cells
				{
					int nIndices[2] = { i, j };
					nIndices[(1 - c)] += dn;
					nIndices[c] += dc;
					nIndices[0] = std::min< int >( std::max< int >( 0 , nIndices[0] ) , width  - 1 );
					nIndices[1] = std::min< int >( std::max< int >( 0 , nIndices[1] ) , height - 1 );
					if( nodeType( nIndices[0] , nIndices[1] )!=1 )
					{
						nodeType( nIndices[0] , nIndices[1] ) = 0;

						Point2D< GeometryReal > texel_pos = Point2D< GeometryReal >( (GeometryReal)nIndices[0]*cellSizeW , (GeometryReal)nIndices[1]*cellSizeH ) - tPos[0];
						Point2D< GeometryReal > barycentricCoord = barycentricMap*texel_pos;

						if( triangleID( nIndices[0] , nIndices[1] )==-1 )
						{
							triangleID(nIndices[0], nIndices[1]) = atlasChart.meshTriangleIndices[tIndex];
							barycentricCoords(nIndices[0], nIndices[1]) = barycentricCoord;
						}
						else //Update the position to the closest triangle
						{
							Point2D< GeometryReal > oldBarycentricCoord = barycentricCoords( nIndices[0] , nIndices[1] );
							Point3D< GeometryReal > oldBarycentricCoord3( (GeometryReal)1. - oldBarycentricCoord[0] - oldBarycentricCoord[1] , oldBarycentricCoord[0] , oldBarycentricCoord[1] );
							Point3D< GeometryReal > newBarycentricCoord3( (GeometryReal)1. -    barycentricCoord[0] -    barycentricCoord[1] ,    barycentricCoord[0] ,    barycentricCoord[1] );
							GeometryReal minOld = std::min< GeometryReal >( std::min< GeometryReal >( oldBarycentricCoord3[0] , oldBarycentricCoord3[1] ) , oldBarycentricCoord3[2] );
							GeometryReal minNew = std::min< GeometryReal >( std::min< GeometryReal >( newBarycentricCoord3[0] , newBarycentricCoord3[1] ) , newBarycentricCoord3[2] );
							if( minNew>minOld )
							{
								triangleID( nIndices[0] , nIndices[1] ) = atlasChart.meshTriangleIndices[tIndex];
								barycentricCoords( nIndices[0] , nIndices[1] ) = barycentricCoord;
							}
						}
					}
				}
			}
		}

		//(2.2) Add texels adjacent to cells that contain triangles
		if ((minCorner[0] + 1 == maxCorner[0]) && (minCorner[1] + 1 == maxCorner[1])) {
			cellType(minCorner[0], minCorner[1]) = 0;
			interiorCellTriangles++;
			for (int di = 0; di < 2; di++) for (int dj = 0; dj < 2; dj++) {
				int nIndices[2] = { minCorner[0] + di, minCorner[1] + dj };
				nIndices[0] = std::min<int>(std::max<int>(0, nIndices[0]), width - 1);
				nIndices[1] = std::min<int>(std::max<int>(0, nIndices[1]), height - 1);
				if (nodeType(nIndices[0], nIndices[1]) != 1) {
					nodeType(nIndices[0], nIndices[1]) = 0;

					Point2D< GeometryReal > texel_pos = Point2D< GeometryReal >( (GeometryReal)nIndices[0]*cellSizeW , (GeometryReal)nIndices[1]*cellSizeH ) - tPos[0];
					Point2D< GeometryReal > barycentricCoord = barycentricMap*texel_pos;

					if (triangleID(nIndices[0], nIndices[1]) == -1) {
						triangleID(nIndices[0], nIndices[1]) = atlasChart.meshTriangleIndices[tIndex];
						barycentricCoords(nIndices[0], nIndices[1]) = barycentricCoord;
					}
					else //Update the position to the closest triangle
					{
						Point2D< GeometryReal > oldBarycentricCoord = barycentricCoords( nIndices[0] , nIndices[1] );
						Point3D< GeometryReal > oldBarycentricCoord3( (GeometryReal)1. - oldBarycentricCoord[0] - oldBarycentricCoord[1] , oldBarycentricCoord[0] , oldBarycentricCoord[1] );
						Point3D< GeometryReal > newBarycentricCoord3( (GeometryReal)1. -    barycentricCoord[0] -    barycentricCoord[1] ,    barycentricCoord[0] ,    barycentricCoord[1] );
						GeometryReal minOld = std::min< GeometryReal >( std::min< GeometryReal >( oldBarycentricCoord3[0] , oldBarycentricCoord3[1] ) , oldBarycentricCoord3[2] );
						GeometryReal minNew = std::min< GeometryReal >( std::min< GeometryReal >( newBarycentricCoord3[0] , newBarycentricCoord3[1] ) , newBarycentricCoord3[2] );
						if( minNew>minOld )
						{
							triangleID( nIndices[0] , nIndices[1] ) = atlasChart.meshTriangleIndices[tIndex];
							barycentricCoords( nIndices[0] , nIndices[1] ) = barycentricCoord;
						}
					}

				}
			}
		}
	}

	//(3) Add interior cells
	for (int j = 0; j < height - 1; j++)for (int i = 0; i < width - 1; i++) {
		if (nodeType(i, j) == 1 && nodeType(i + 1, j) == 1 && nodeType(i, j + 1) == 1 && nodeType(i + 1, j + 1) == 1 && cellType(i, j) == -1) {
			cellType(i, j) = 1;
		}
	}

	int deepRadius = multigridBlockInfo.boundaryDilation + 1;
	for (int j = deepRadius; j < height - deepRadius; j++)for (int i = deepRadius; i < width - deepRadius; i++) {
		bool validDeepNode = true;
		for (int di = -deepRadius; di < deepRadius; di++) for (int dj = -deepRadius; dj < deepRadius; dj++) {
			if (cellType(i + di, j + dj) != 1) validDeepNode = false;
		}
		if (validDeepNode) nodeType(i, j) = 2;
	}

	gridChart.nodeType = nodeType;
	gridChart.cellType = cellType;

	// (5) Enumerate variables in raster order
	gridChart.globalIndexTexelOffset = lastGlobalTexelIndex;
	gridChart.globalIndexInteriorTexelOffset = lastGlobalTexelInteriorIndex;
	gridChart.globalIndexCellOffset = lastGlobalCellIndex;
	gridChart.globalIndexInteriorCellOffset = lastGlobalInteriorCellIndex;
	gridChart.globalTexelDeepOffset = lastGlobalTexelDeepIndex;
	gridChart.globalTexelBoundaryOffset = lastGlobalTexelBoundaryIndex;

	Image<int> & globalTexelIndex = gridChart.globalTexelIndex;
	globalTexelIndex.resize(width, height);
	for (int i = 0; i < globalTexelIndex.size(); i++)globalTexelIndex[i] = -1;

	Image<int> & globalTexelInteriorIndex = gridChart.globalInteriorTexelIndex;
	globalTexelInteriorIndex.resize(width, height);
	for (int i = 0; i < globalTexelInteriorIndex.size(); i++)globalTexelInteriorIndex[i] = -1;

	Image<int> & globalTexelDeepIndex = gridChart.globalTexelDeepIndex;
	globalTexelDeepIndex.resize(width, height);
	for (int i = 0; i < globalTexelDeepIndex.size(); i++)globalTexelDeepIndex[i] = -1;

	Image<int> & globalTexelBoundaryIndex = gridChart.globalTexelBoundaryIndex;
	globalTexelBoundaryIndex.resize(width, height);
	for (int i = 0; i < globalTexelBoundaryIndex.size(); i++)globalTexelBoundaryIndex[i] = -1;

	int lastLocalTexelIndex = 0;
	int lastLocalTexelInteiorIndex = 0;
	int lastLocalTexelDeepIndex = 0;
	int lastLocalTexelBoundaryIndex = 0;
	for (int j = 0; j < height; j++)for (int i = 0; i < width; i++) {
		if (gridChart.nodeType(i, j) != -1) {
			globalTexelIndex(i, j) = lastGlobalTexelIndex + lastLocalTexelIndex;
			lastLocalTexelIndex++;

			GridNodeInfo currentNodeInfo;
			currentNodeInfo.ci = i;
			currentNodeInfo.cj = j;
			currentNodeInfo.chartID = chartID;
			currentNodeInfo.nodeType = gridChart.nodeType(i, j);
			nodeInfo.push_back(currentNodeInfo);

		}
		if (gridChart.nodeType(i, j) >= 1) {
			globalTexelInteriorIndex(i, j) = lastGlobalTexelInteriorIndex + lastLocalTexelInteiorIndex;
			lastLocalTexelInteiorIndex++;
		}

		if (gridChart.nodeType(i, j) == 0 || gridChart.nodeType(i, j) == 1) {
			globalTexelBoundaryIndex(i, j) = lastGlobalTexelBoundaryIndex + lastLocalTexelBoundaryIndex;
			lastLocalTexelBoundaryIndex++;
		}

		if (gridChart.nodeType(i, j) == 2) {
			globalTexelDeepIndex(i, j) = lastGlobalTexelDeepIndex + lastLocalTexelDeepIndex;
			lastLocalTexelDeepIndex++;
		}
	}
	lastGlobalTexelIndex += lastLocalTexelIndex;
	lastGlobalTexelInteriorIndex += lastLocalTexelInteiorIndex;
	lastGlobalTexelDeepIndex += lastLocalTexelDeepIndex;
	lastGlobalTexelBoundaryIndex += lastLocalTexelBoundaryIndex;

	Image<int> & localCellIndex = gridChart.localCellIndex;
	std::vector< BilinearElementIndex > & bilinearElementIndices = gridChart.bilinearElementIndices;
	localCellIndex.resize(width - 1, height - 1);
	for (int i = 0; i < localCellIndex.size(); i++)localCellIndex[i] = -1;

	Image<int> & localInteriorCellIndex = gridChart.localInteriorCellIndex;
	std::vector< BilinearElementIndex > & interiorCellCorners = gridChart.interiorCellCorners;
	localInteriorCellIndex.resize(width - 1, height - 1);
	for (int i = 0; i < localInteriorCellIndex.size(); i++)localInteriorCellIndex[i] = -1;

	std::vector< BilinearElementIndex > & interiorCellGlobalCorners = gridChart.interiorCellGlobalCorners;

	Image<int> & localBoundaryCellIndex = gridChart.localBoundaryCellIndex;
	localBoundaryCellIndex.resize(width - 1, height - 1);
	for (int i = 0; i < localBoundaryCellIndex.size(); i++)localBoundaryCellIndex[i] = -1;

	std::vector<int> & interiorCellIndexToLocalCellIndex = gridChart.interiorCellIndexToLocalCellIndex;
	std::vector<int> & boundaryCellIndexToLocalCellIndex = gridChart.boundaryCellIndexToLocalCellIndex;

	int lastLocalCellIndex = 0;
	int lastLocalBoundaryCellIndex = 0;
	int lastLocalInteriorCellIndex = 0;

	for (int j = 0; j < height - 1; j++)for (int i = 0; i < width - 1; i++) {
		if (gridChart.cellType(i, j) != -1) {
			localCellIndex(i, j) = lastLocalCellIndex;
			lastLocalCellIndex++;
			int globalTexelIndices[4] = { globalTexelIndex(i, j), globalTexelIndex(i + 1, j), globalTexelIndex(i + 1, j + 1), globalTexelIndex(i, j + 1) };
			if (globalTexelIndices[0] != -1 && globalTexelIndices[1] != -1 && globalTexelIndices[2] != -1 && globalTexelIndices[3] != -1) {
				bilinearElementIndices.push_back( BilinearElementIndex( globalTexelIndices[0] , globalTexelIndices[1] , globalTexelIndices[2] , globalTexelIndices[3] ) );
			}
			else THROW( "Active cell adjacent to unactive node" );

			if (gridChart.cellType(i, j) == 0) {
				localBoundaryCellIndex(i, j) = lastLocalBoundaryCellIndex;
				boundaryCellIndexToLocalCellIndex.push_back(lastLocalCellIndex - 1);
				lastLocalBoundaryCellIndex++;
			}
			else {
				localInteriorCellIndex(i, j) = lastLocalInteriorCellIndex;
				interiorCellIndexToLocalCellIndex.push_back(lastLocalCellIndex - 1);

				int globalTexelInteriorIndices[4] = { globalTexelInteriorIndex(i, j), globalTexelInteriorIndex(i + 1, j), globalTexelInteriorIndex(i + 1, j + 1), globalTexelInteriorIndex(i, j + 1) };
				if (globalTexelInteriorIndices[0] != -1 && globalTexelInteriorIndices[1] != -1 && globalTexelInteriorIndices[2] != -1 && globalTexelInteriorIndices[3] != -1) {
					interiorCellCorners.push_back( BilinearElementIndex( globalTexelInteriorIndices[0] , globalTexelInteriorIndices[1] , globalTexelInteriorIndices[2] , globalTexelInteriorIndices[3] ) );
					interiorCellGlobalCorners.push_back( BilinearElementIndex( globalTexelIndices[0] , globalTexelIndices[1] , globalTexelIndices[2] , globalTexelIndices[3] ) );
				}
				else THROW( "Interior cell adjacent to non interior node" );
				lastLocalInteriorCellIndex++;
			}
		}
	}

	lastGlobalCellIndex += lastLocalCellIndex;
	lastGlobalBoundaryCellIndex += lastLocalBoundaryCellIndex;
	lastGlobalInteriorCellIndex += lastLocalInteriorCellIndex;

	gridChart.numInteriorCells = lastLocalInteriorCellIndex;
	gridChart.numBoundaryCells = lastLocalBoundaryCellIndex;

	// (6) Construct raster lines
	for (int j = 0; j < height; j++) {
		int offset = 0;
		bool firstSegment = true;
		bool previousDeep = false;
		int rasterStart = -1;
		while (offset < width) {
			bool currentIsDeep = nodeType(offset, j) == 2;
			if (currentIsDeep && !previousDeep) {
				rasterStart = offset; //Start raster line
				if (firstSegment) {
					firstSegment = false;
					SegmentedRasterLine newSegmentLine;
					segmentedLines.push_back(newSegmentLine);
				}
			}
			if (!currentIsDeep && previousDeep) { //Terminate raster line
				RasterLine newLine;
				newLine.lineStartIndex = globalTexelIndex(rasterStart, j);
				newLine.lineEndIndex = globalTexelIndex(offset - 1, j);
				newLine.prevLineIndex = globalTexelIndex(rasterStart, j - 1);
				newLine.nextLineIndex = globalTexelIndex(rasterStart, j + 1);
				newLine.coeffStartIndex = globalTexelDeepIndex(rasterStart, j);

				if( newLine.lineStartIndex==-1 || newLine.lineEndIndex==-1 || newLine.prevLineIndex==-1 || newLine.nextLineIndex==-1 ) THROW( "Inavlid Indexing" );
				rasterLines.push_back(newLine);

				SegmentedRasterLine & newSegmentLine = segmentedLines.back();
				newSegmentLine.segments.push_back(newLine);
			}
			previousDeep = currentIsDeep;
			offset++;
		}
	}

	//Initialize thread tasks
	int blockHorizontalOffset = multigridBlockInfo.blockWidth - multigridBlockInfo.paddingWidth;
	int blockVerticalOffset = multigridBlockInfo.blockHeight - multigridBlockInfo.paddingHeight;
	int numHorizontalBlocks = ((width - multigridBlockInfo.paddingWidth - 1) / blockHorizontalOffset) + 1;
	int numVerticalBlocks = ((height - multigridBlockInfo.paddingHeight - 1) / blockVerticalOffset) + 1;

	for (int bj = 0; bj < numVerticalBlocks; bj++) {
		ThreadTask threadTask;
		int taskDeepTexels = 0;

		int blockVerticalStart = bj*blockVerticalOffset;
		int blockVerticalEnd = std::min<int>((bj + 1)*blockVerticalOffset + multigridBlockInfo.paddingHeight + 2 - 1, height - 1);

		for (int bi = 0; bi < numHorizontalBlocks; bi++) {
			BlockTask blockTask;

			int blockHorizontalStart = bi*blockHorizontalOffset;
			int blockHorizontalEnd = std::min<int>((bi + 1)*blockHorizontalOffset + multigridBlockInfo.paddingWidth + 2 - 1, width - 1);

			//Deep texel within rows[blockVerticalStart + 1,blockVerticalEnd - 1]column [blockHorizontalStart + 1,blockHorizontalEnd - 1] 
			std::vector<BlockDeepSegmentedLine>  &  blockDeepSegmentedLines = blockTask.blockDeepSegmentedLines;
			for (int j = blockVerticalStart + 1; j <= blockVerticalEnd - 1; j++) {
				BlockDeepSegmentedLine segmentedLine;
				int offset = blockHorizontalStart + 1;
				int segmentStart = -1;
				bool previousDeep = false;
				while (offset <= blockHorizontalEnd - 1) {
					bool currentDeep = (globalTexelDeepIndex(offset, j) != -1);
					if (currentDeep && !previousDeep) {//Start segment
						segmentStart = offset;
					}

					if ((previousDeep && !currentDeep) || (currentDeep && offset == blockHorizontalEnd - 1)) {//Terminate segment
						BlockDeepSegment deepSegment;
						deepSegment.currentStart = globalTexelIndex(segmentStart, j);
						if (currentDeep)deepSegment.currentEnd = globalTexelIndex(offset, j);
						else deepSegment.currentEnd = globalTexelIndex(offset - 1, j);
						deepSegment.previousStart = globalTexelIndex(segmentStart, j - 1);
						deepSegment.nextStart = globalTexelIndex(segmentStart, j + 1);

						deepSegment.deepStart = globalTexelDeepIndex(segmentStart, j);

						segmentedLine.blockDeepSegments.push_back(deepSegment);
					}
					previousDeep = currentDeep;
					offset++;
				}
				if (segmentedLine.blockDeepSegments.size()) {
					blockDeepSegmentedLines.push_back(segmentedLine);
				}
			}

			if (blockDeepSegmentedLines.size()) {

				//Deep texel within rows[blockVerticalStart,blockVerticalEnd]column [blockHorizontalStart,blockHorizontalEnd] 
				std::vector<BlockDeepSegmentedLine>  &  blockPaddedSegmentedLines = blockTask.blockPaddedSegmentedLines;
				for (int j = blockVerticalStart; j <= blockVerticalEnd; j++) {
					BlockDeepSegmentedLine segmentedLine;
					int offset = blockHorizontalStart;
					int segmentStart = -1;
					bool previousDeep = false;
					while (offset <= blockHorizontalEnd) {
						bool currentDeep = (globalTexelDeepIndex(offset, j) != -1);
						if (currentDeep && !previousDeep) {//Start segment
							segmentStart = offset;
						}

						if ((previousDeep && !currentDeep) || (currentDeep && offset == blockHorizontalEnd)) {//Terminate segment
							BlockDeepSegment deepSegment;
							deepSegment.currentStart = globalTexelIndex(segmentStart, j);
							if (currentDeep)deepSegment.currentEnd = globalTexelIndex(offset, j);
							else deepSegment.currentEnd = globalTexelIndex(offset - 1, j);
							deepSegment.previousStart = globalTexelIndex(segmentStart, j - 1);
							deepSegment.nextStart = globalTexelIndex(segmentStart, j + 1);

							deepSegment.deepStart = globalTexelDeepIndex(segmentStart, j);

							segmentedLine.blockDeepSegments.push_back(deepSegment);
						}
						previousDeep = currentDeep;
						offset++;
					}
					if (segmentedLine.blockDeepSegments.size()) {
						blockPaddedSegmentedLines.push_back(segmentedLine);
					}
				}
				threadTask.blockTasks.push_back(blockTask);
			}
		}
		if (threadTask.blockTasks.size()) {
			threadTask.taskDeepTexels = taskDeepTexels;
			threadTasks.push_back(threadTask);
		}
	}
}

template< typename GeometryReal >
void InitializeGridCharts( const std::vector< AtlasChart< GeometryReal > > &atlasCharts , const GeometryReal &cellSizeW , const GeometryReal &cellSizeH , std::vector< GridNodeInfo > &nodeInfo , std::vector< GridChart< GeometryReal > > &gridCharts , std::vector< RasterLine > &rasterLines , std::vector< SegmentedRasterLine > &segmentedLines , std::vector< ThreadTask > &threadTasks , int &numTexels , int &numInteriorTexels , int &numDeepTexels , int &numBoundaryTexels , int &numCells , int &numBoundaryCells , int &numInteriorCells , const MultigridBlockInfo &multigridBlockInfo )
{
	gridCharts.resize( atlasCharts.size() );

	numTexels = 0;
	numInteriorTexels = 0;
	numDeepTexels = 0;
	numBoundaryTexels = 0;
	numCells = 0;
	numBoundaryCells = 0;
	numInteriorCells = 0;

	for( int i=0 ; i<atlasCharts.size() ; i++ )
	{
		int halfSize[2][2];
		for( int c=0 ; c<2 ; c++ )
		{
			if( c==0 )
			{
				halfSize[c][0] = (int)ceil( ( atlasCharts[i].gridOrigin[c] - atlasCharts[i].minCorner[c] ) / cellSizeW );
				halfSize[c][1] = (int)ceil( ( atlasCharts[i].maxCorner[c] - atlasCharts[i].gridOrigin[c] ) / cellSizeW );
				gridCharts[i].corner[c] = atlasCharts[i].gridOrigin[c] - cellSizeW * (GeometryReal)halfSize[c][0];
			}
			else
			{
				halfSize[c][0] = (int)ceil( ( atlasCharts[i].gridOrigin[c] - atlasCharts[i].minCorner[c] ) / cellSizeH );
				halfSize[c][1] = (int)ceil( ( atlasCharts[i].maxCorner[c] - atlasCharts[i].gridOrigin[c] ) / cellSizeH );
				gridCharts[i].corner[c] = atlasCharts[i].gridOrigin[c] - cellSizeH * (GeometryReal)halfSize[c][0];
			}
			gridCharts[i].centerOffset[c] = halfSize[c][0];
			gridCharts[i].cornerCoords[c] = atlasCharts[i].originCoords[c] - halfSize[c][0];
		}

		gridCharts[i].width = halfSize[0][0] + halfSize[0][1] + 1;
		gridCharts[i].height = halfSize[1][0] + halfSize[1][1] + 1;
		gridCharts[i].cellSizeW = cellSizeW;
		gridCharts[i].cellSizeH = cellSizeH;
		int approxGridResolution = std::max<int>((int)ceil(1.0 / cellSizeW), (int)ceil(1.0 / cellSizeH));
		gridCharts[i].gridIndexOffset = 2 * approxGridResolution*approxGridResolution*i;
		InitializeGridChartsActiveNodes(i, atlasCharts[i], gridCharts[i], nodeInfo, rasterLines, segmentedLines, threadTasks, numTexels, numInteriorTexels, numDeepTexels, numBoundaryTexels, numCells, numBoundaryCells, numInteriorCells, multigridBlockInfo);
	}
}

template< typename GeometryReal >
void InitializeTextureNodes( const std::vector< GridChart< GeometryReal > > &gridCharts , std::vector< TextureNodeInfo< GeometryReal > > &textureNodes )
{
	for (int c = 0; c < gridCharts.size(); c++)
	{
		const GridChart< GeometryReal > &gridChart = gridCharts[c];
		const Image<int> globalTexelIndex = gridChart.globalTexelIndex;
		for (int j = 0; j < gridChart.height; j++)for (int i = 0; i < gridChart.width; i++) {
			int texelIndex = globalTexelIndex(i, j);
			if (texelIndex != -1)
			{
				TextureNodeInfo< GeometryReal > textureNode;
				textureNode.ci = gridChart.cornerCoords[0] + i;
				textureNode.cj = gridChart.cornerCoords[1] + j;
				textureNode.chartID = c;
				textureNode.tID = gridChart.triangleID(i, j);
				textureNode.barycentricCoords = gridChart.barycentricCoords(i, j);
				textureNode.isInterior = gridChart.nodeType(i, j) > 0;
				textureNodes.push_back(textureNode);
			}
		}
	}
}

template< typename GeometryReal >
void InitializeCellNodes( const std::vector< GridChart< GeometryReal > > &gridCharts , std::vector< BilinearElementIndex > &cellNodes )
{
	for( int c=0 ; c<gridCharts.size() ; c++ )
	{
		const GridChart< GeometryReal >& gridChart = gridCharts[c];
		const std::vector< BilinearElementIndex >& localBilinearElementIndices = gridChart.bilinearElementIndices;
		for( int i=0 ; i<localBilinearElementIndices.size() ; i++ ) cellNodes.push_back( localBilinearElementIndices[i] );
	}
}

template< typename GeometryReal , typename MatrixReal >
void InitializeAtlasHierachicalBoundaryCoefficients( const GridAtlas< GeometryReal , MatrixReal > &fineAtlas , GridAtlas< GeometryReal , MatrixReal > &coarseAtlas , SparseMatrix< MatrixReal , int > &boundaryCoarseFineProlongation , std::vector< BoundaryDeepIndex > &boundaryDeepIndices , std::vector< BoundaryBoundaryIndex< MatrixReal > > &boundaryBoundaryIndices )
{
	std::vector< Eigen::Triplet< MatrixReal > > prolongationTriplets;

	const std::vector<GridNodeInfo> & coarseNodeInfo = coarseAtlas.nodeInfo;
	for (int i = 0; i < coarseAtlas.boundaryGlobalIndex.size(); i++) {
		const GridNodeInfo & currentNode = coarseNodeInfo[coarseAtlas.boundaryGlobalIndex[i]];
		const GridChart< GeometryReal > &coarseChart = coarseAtlas.gridCharts[currentNode.chartID];
		const GridChart< GeometryReal > &fineChart = fineAtlas.gridCharts[currentNode.chartID];

		int coarseChartWidth = coarseChart.nodeType.width();
		int coarseChartHeight = coarseChart.nodeType.height();

		int fineChartWidth = fineChart.nodeType.width();
		int fineChartHeight = fineChart.nodeType.height();

		int ci = currentNode.ci;
		int cj = currentNode.cj;

		int numBoundaryNeighbours = 0;// -1: not defined, 0: boundary, 1: interior 
		int boundary_id[9];
		int boundary_offset_i[9];
		int boundary_offset_j[9];
		for (int li = -1; li < 2; li++)for (int lj = -1; lj < 2; lj++) {
			int pi = ci + li;
			int pj = cj + lj;
			if (pi >= 0 && pi < coarseChartWidth && pj >= 0 && pj < coarseChartHeight) {
				int coarseGlobalIndex = coarseChart.globalTexelIndex(pi, pj);
				if (coarseGlobalIndex != -1) {
					int coarseBoundaryIndex = coarseAtlas.boundaryAndDeepIndex[coarseGlobalIndex];
					if (coarseBoundaryIndex < 0) {//Deep
						int coarseDeepIndex = -coarseBoundaryIndex - 1;
						BoundaryDeepIndex bdIndex;
						bdIndex.boundaryIndex = i;
						bdIndex.deepGlobalIndex = coarseGlobalIndex;
						bdIndex.deepIndex = coarseDeepIndex;
						bdIndex.offset = (1 - lj) * 3 + (1 - li);
						boundaryDeepIndices.push_back(bdIndex);
					}
					else { //Boundary
						coarseBoundaryIndex = coarseBoundaryIndex - 1;
						boundary_id[numBoundaryNeighbours] = coarseBoundaryIndex;
						boundary_offset_i[numBoundaryNeighbours] = 2 * li;
						boundary_offset_j[numBoundaryNeighbours] = 2 * lj;
						numBoundaryNeighbours++;
					}

				}
			}
		}

		int ri = ci - coarseChart.centerOffset[0];
		int rj = cj - coarseChart.centerOffset[1];

		int fi = fineChart.centerOffset[0] + 2 * ri;
		int fj = fineChart.centerOffset[1] + 2 * rj;


		for (int di = -1; di < 2; di++)for (int dj = -1; dj < 2; dj++) {
			int pi = fi + di;
			int pj = fj + dj;
			if (pi >= 0 && pi < fineChartWidth && pj >= 0 && pj < fineChartHeight) {
				int fineGlobalIndex = fineChart.globalTexelIndex(pi, pj);
				if (fineGlobalIndex != -1) {
					int fineBoundaryIndex = fineAtlas.boundaryAndDeepIndex[fineGlobalIndex];
					MatrixReal weight = (MatrixReal)( 1.0 / (1.0 + std::abs(di) ) / ( 1.0 + std::abs(dj) ) );
					if (fineBoundaryIndex > 0) {//Boundary
						fineBoundaryIndex -= 1;
						prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( fineBoundaryIndex , i , weight ) );

						for (int ki = -1; ki < 2; ki++)for (int kj = -1; kj < 2; kj++) {
							int qi = pi + ki;
							int qj = pj + kj;
							if (qi >= 0 && qi < fineChartWidth && qj >= 0 && qj < fineChartHeight) {
								int neighboutFineGlobalIndex = fineChart.globalTexelIndex(qi, qj);
								if (neighboutFineGlobalIndex != -1) {
									int neighbourFineBoundaryIndex = fineAtlas.boundaryAndDeepIndex[neighboutFineGlobalIndex];
									if (neighbourFineBoundaryIndex < 0) {//Deep
										int neighbourFineDeepIndex = -neighbourFineBoundaryIndex - 1;
										int oi = di + ki;
										int oj = dj + kj;
										for (int n = 0; n < numBoundaryNeighbours; n++) {
											MatrixReal diff_i = (MatrixReal)std::abs( oi - boundary_offset_i[n] );
											MatrixReal diff_j = (MatrixReal)std::abs( oj - boundary_offset_j[n] );
											if (diff_i < 1.5 && diff_j < 1.5) {
												MatrixReal weight2 = (MatrixReal)( 1.0 / (1.0 + diff_i ) / ( 1.0 + diff_j ) );
												BoundaryBoundaryIndex< MatrixReal > bbIndex;
												bbIndex.coarsePrincipalBoundaryIndex = i;
												bbIndex.coarseSecondaryBoundaryIndex = boundary_id[n];
												bbIndex.fineDeepIndex = neighbourFineDeepIndex;
												bbIndex.offset = (1 - kj) * 3 + (1 - ki);
												bbIndex.weight = weight * weight2;
												boundaryBoundaryIndices.push_back(bbIndex);
											}
										}
									}
								}
							}
						}
					}
					else {//deep
						int fineDeepIndex = -fineBoundaryIndex - 1;
						for (int ki = -1; ki < 2; ki++)for (int kj = -1; kj < 2; kj++) {
							int oi = di + ki;
							int oj = dj + kj;

							for (int n = 0; n < numBoundaryNeighbours; n++)
							{
								MatrixReal diff_i = (MatrixReal)std::abs( oi - boundary_offset_i[n] );
								MatrixReal diff_j = (MatrixReal)std::abs( oj - boundary_offset_j[n] );
								if (diff_i < 1.5 && diff_j < 1.5)
								{
									MatrixReal weight2 = (MatrixReal)( 1.0 / ( 1.0 + diff_i ) / ( 1.0 + diff_j ) );
									BoundaryBoundaryIndex< MatrixReal > bbIndex;
									bbIndex.coarsePrincipalBoundaryIndex = i;
									bbIndex.coarseSecondaryBoundaryIndex = boundary_id[n];
									bbIndex.fineDeepIndex = fineDeepIndex;
									bbIndex.offset = (1 + kj) * 3 + (1 + ki);
									bbIndex.weight = weight * weight2;
									boundaryBoundaryIndices.push_back(bbIndex);
								}
							}
						}
					}

				}
			}
		}
	}

	int numCoarseBoundaryTexels = coarseAtlas.numBoundaryTexels;
	int numFineBoundaryTexels = fineAtlas.numBoundaryTexels;
	boundaryCoarseFineProlongation = SetSparseMatrix( prolongationTriplets , numFineBoundaryTexels , numCoarseBoundaryTexels , false );
}


template< typename GeometryReal , typename MatrixReal >
void InitializeHierarchy( const int width , const int height , HierarchicalSystem< GeometryReal , MatrixReal > &hierarchy , std::vector< AtlasChart< GeometryReal > > &atlasCharts , const int levels , const MultigridBlockInfo &multigridBlockInfo , bool verbose=false )
{
	Miscellany::Timer timer;

	std::vector< GridAtlas< GeometryReal , MatrixReal > > &gridAtlases = hierarchy.gridAtlases;
	gridAtlases.resize(levels);

	for( int i=0 ; i<levels ; i++ )
	{
		int reductionFactor = 1 << (i);
		GeometryReal cellSizeW = (GeometryReal)reductionFactor / width;
		GeometryReal cellSizeH = (GeometryReal)reductionFactor / height;
		//gridAtlases[i].resolution = resolution;
		InitializeGridCharts( atlasCharts , cellSizeW , cellSizeH , gridAtlases[i].nodeInfo , gridAtlases[i].gridCharts , gridAtlases[i].rasterLines , gridAtlases[i].segmentedLines , gridAtlases[i].threadTasks , gridAtlases[i].numTexels , gridAtlases[i].numInteriorTexels , gridAtlases[i].numDeepTexels , gridAtlases[i].numBoundaryTexels , gridAtlases[i].numCells , gridAtlases[i].numBoundaryCells , gridAtlases[i].numInteriorCells , multigridBlockInfo );

		if( gridAtlases[i].numTexels!=gridAtlases[i].numBoundaryTexels+gridAtlases[i].numDeepTexels )
			THROW( "Boundary and deep texels does not form a partition: " , gridAtlases[i].numTexels , " != " , gridAtlases[i].numBoundaryTexels , " + " , gridAtlases[i].numDeepTexels );

		InitializeBoundaryAndDeepTexelIndexing( gridAtlases[i].gridCharts , gridAtlases[i].numTexels , gridAtlases[i].boundaryAndDeepIndex , gridAtlases[i].boundaryGlobalIndex , gridAtlases[i].deepGlobalIndex );

		if( verbose )
		{
			printf( "Level %d:\n" , i );
			printf( "\t Active texels %d\n", gridAtlases[i].numTexels );
			printf( "\t Deep texels %d\n" , gridAtlases[i].numDeepTexels );
			printf( "\t Boundary texels %d\n" , gridAtlases[i].numBoundaryTexels );
			printf( "\t Active cells %d\n" , gridAtlases[i].numCells );
			printf( "\t Interior texture nodes %d\n" , gridAtlases[i].numInteriorTexels );
			printf( "\t Deep lines %d\n" , (int)gridAtlases[i].segmentedLines.size() );
			printf( "\t Deep segments %d\n" , (int)gridAtlases[i].rasterLines.size() );
		}
	}
	if( verbose ) printf( "Hierarchy structure =  %.4f \n" , timer.elapsed() );

	if( verbose ) timer.reset();
	hierarchy.boundaryRestriction.resize( levels-1 );
	for( int i=0 ; i<levels-1 ; i++ ) InitializeAtlasHierachicalRestriction( gridAtlases[i] , gridAtlases[i+1] , hierarchy.boundaryRestriction[i] );
	for( int i=0 ; i<levels-1 ; i++ ) InitializeAtlasHierachicalProlongation( gridAtlases[i] , gridAtlases[i+1] );

	if( verbose ) printf( "Hierarchy prolongation and restriction =  %.4f \n" , timer.elapsed() );

	hierarchy.boundaryCoarseFineProlongation.resize(levels);
	hierarchy.boundaryFineCoarseRestriction.resize(levels);
	hierarchy.boundaryDeepIndices.resize(levels);
	hierarchy.boundaryBoundaryIndices.resize(levels);

	if( verbose ) timer.reset();

	for( int i=1 ; i<levels ; i++ )
	{
		InitializeAtlasHierachicalBoundaryCoefficients( hierarchy.gridAtlases[i-1] , hierarchy.gridAtlases[i] , hierarchy.boundaryCoarseFineProlongation[i] , hierarchy.boundaryDeepIndices[i] , hierarchy.boundaryBoundaryIndices[i] );
		hierarchy.boundaryFineCoarseRestriction[i-1] = hierarchy.boundaryCoarseFineProlongation[i].transpose();
	}
	if( verbose ) printf( "Boundary coefficient initialization =  %.4f \n" , timer.elapsed() );
}

template< typename GeometryReal , typename MatrixReal >
void InitializeHierarchy
(
	OrientedTexturedTriangleMesh< GeometryReal > & mesh , const int width , const int height , const int levels , std::vector< TextureNodeInfo< GeometryReal > > &textureNodes , std::vector< BilinearElementIndex > &cellNodes ,
	HierarchicalSystem< GeometryReal , MatrixReal > &hierarchy , std::vector< AtlasChart< GeometryReal > > &atlasCharts , const MultigridBlockInfo &multigridBlockInfo , bool verbose=false , bool detailVerbose=false , bool computeProlongation=false
)
{

	Miscellany::Timer timer;

	//(1) Initialize Atlas Mesh
	AtlasMesh< GeometryReal > atlasMesh;
	std::vector< int > oppositeHalfEdge;
	std::unordered_map< int , int > boundaryVerticesIndices;
	int numBoundaryVertices;
	bool isClosedMesh;

	InitializeAtlasMesh( mesh , width , height , atlasMesh , atlasCharts , oppositeHalfEdge , boundaryVerticesIndices , numBoundaryVertices , isClosedMesh , detailVerbose );

	if( verbose ) timer.reset();

	//(2) Initialize Hierarchy
	InitializeHierarchy( width , height , hierarchy , atlasCharts , levels , multigridBlockInfo , detailVerbose );
	if( detailVerbose ) printf( "Hierarchy construction =  %.4f \n" , timer.elapsed() );

	//(3) Initialize fine level texture nodes and cells
	InitializeTextureNodes( hierarchy.gridAtlases[0].gridCharts , textureNodes );
	InitializeCellNodes( hierarchy.gridAtlases[0].gridCharts , cellNodes );

	//(4) Initialize boundary triangulation
	InitializeBoundaryTriangulation( hierarchy.gridAtlases[0] , atlasMesh , atlasCharts , oppositeHalfEdge , boundaryVerticesIndices , numBoundaryVertices , isClosedMesh );

	//(5) Unnnecesary
	if( computeProlongation )
	{
		Eigen::SparseMatrix< MatrixReal > __prolongation;
		InitializeProlongation( hierarchy.gridAtlases[0].numInteriorTexels , hierarchy.gridAtlases[0].numFineNodes , hierarchy.gridAtlases[0].numTexels , hierarchy.gridAtlases[0].gridCharts , hierarchy.gridAtlases[0].nodeInfo , __prolongation );
		hierarchy.gridAtlases[0].coarseToFineNodeProlongation = __prolongation;
	}
}