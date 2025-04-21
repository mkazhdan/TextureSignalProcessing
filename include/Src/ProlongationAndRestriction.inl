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
#include <Misha/MultiThreading.h>

template< typename GeometryReal , typename MatrixReal >
void InitializeProlongation( int numInteriorTexels , int numFineNodes , int numCoarseNodes , const std::vector< GridChart< GeometryReal > > &gridCharts , const std::vector< GridNodeInfo > &nodeInfo , Eigen::SparseMatrix< MatrixReal > &prolongation )
{
	std::vector< Eigen::Triplet< MatrixReal > > prolongationTriplets;
	std::unordered_set<int> coveredNodes;

	std::vector< int > interiorTexelIndices(numCoarseNodes, -1);
	for (int i = 0; i < gridCharts.size(); i++)
	{
		const GridChart< GeometryReal > &gridChart = gridCharts[i];
		for (int j = 0; j < gridChart.texelIndices.size(); j++) {
#if 1
			MK_ERROR_OUT( "Method disabled" );
#else
			if (gridChart.texelIndices[j].cobmiend != -1 && gridChart.texelIndices[j].interiorOrCovered != -1) {
				interiorTexelIndices[gridChart.texelIndices[j].combined] = gridChart.texelIndices[j].interiorOrCovered;
			}
#endif
		}
	}

	for( int i=0 ; i<interiorTexelIndices.size() ; i++ )
	{
		if( interiorTexelIndices[i]!=-1 )
		{
			if( interiorTexelIndices[i]<0 || interiorTexelIndices[i]>numFineNodes ) MK_THROW( "Out of bound index" );
			prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( interiorTexelIndices[i] , i , (MatrixReal)1. ) );
			coveredNodes.insert(i);
		}
	}

	int numAuxiliaryNodes = numFineNodes - numInteriorTexels;
	std::vector<int> auxiliaryNodesDegree(numAuxiliaryNodes, 0);

	for (int i = 0; i < gridCharts.size(); i++)
	{
		const GridChart< GeometryReal > &gridChart = gridCharts[i];
		for (int j = 0; j < gridChart.auxiliaryNodes.size(); j++) {
			int auxiliaryID = gridChart.auxiliaryNodes[j].index - numInteriorTexels;
			auxiliaryNodesDegree[auxiliaryID]++;
		}
	}

	GeometryReal precision_error = (GeometryReal)1e-10;

	std::vector<GeometryReal> auxiliaryNodesCumWeight(numAuxiliaryNodes, 0);

	for (int i = 0; i < gridCharts.size(); i++) {
		const GridChart< GeometryReal > &gridChart = gridCharts[i];
		for (int j = 0; j < gridChart.auxiliaryNodes.size(); j++) {
			int auxiliaryID = gridChart.auxiliaryNodes[j].index - numInteriorTexels;
			int nodeDegree = auxiliaryNodesDegree[auxiliaryID];
			Point2D< GeometryReal > nodePosition = gridChart.auxiliaryNodes[j].position;
			int corner[2] = { (int)floor(nodePosition[0] / gridChart.cellSizeW), (int)floor(nodePosition[1] / gridChart.cellSizeH) };
			unsigned int cellId = gridChart.cellIndices( corner[0] , corner[1] ).combined;

			nodePosition[0] /= gridChart.cellSizeW;
			nodePosition[1] /= gridChart.cellSizeH;
			nodePosition[0] -= (GeometryReal)corner[0];
			nodePosition[1] -= (GeometryReal)corner[1];
			if( nodePosition[0] < 0-precision_error || nodePosition[0] > 1+precision_error || nodePosition[1] < 0-precision_error || nodePosition[1] > 1+precision_error )
				MK_THROW( "Sample out of unit box: (" , nodePosition[0] , " " , nodePosition[1] , ")" );
			for (int k = 0; k < 4; k++)
			{
				GeometryReal texelWeight = BilinearElementValue( k , nodePosition ) / nodeDegree;
				if( fabs(texelWeight)>1e-11 )
				{
					auxiliaryNodesCumWeight[auxiliaryID] += texelWeight;
					int texelIndex = gridChart.combinedCellCombinedBilinearElementIndices[cellId][k];
					if( nodeInfo[texelIndex].texelType==TexelType::InteriorSupported )
						MK_THROW( "Interior-supported texel cannot be in the support of an auxiliary node. Weight " , texelWeight , " (B)" );
					coveredNodes.insert(texelIndex);
					if( gridChart.auxiliaryNodes[j].index<numInteriorTexels || gridChart.auxiliaryNodes[j].index>numFineNodes || texelIndex<0 || texelIndex>numCoarseNodes )
						MK_THROW( "Out of bounds index" );
					prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( gridChart.auxiliaryNodes[j].index , texelIndex , (MatrixReal)texelWeight ) );
				}
			}
		}
	}

	for( int i=0 ; i<numAuxiliaryNodes ; i++ ) if( fabs( auxiliaryNodesCumWeight[i]-1.0 )>precision_error ) MK_THROW( "Cum weight out of precision " , auxiliaryNodesCumWeight[i] );

	if( coveredNodes.size()!=numCoarseNodes ) MK_THROW( "Total active texels does not match total texels" );

	printf("Prolongation operator dimensions %d x %d \n", numFineNodes, numCoarseNodes);
	prolongation.resize(numFineNodes, numCoarseNodes);
	prolongation.setFromTriplets(prolongationTriplets.begin(), prolongationTriplets.end());
}

template< typename GeometryReal , typename MatrixReal >
void InitializeAtlasHierachicalProlongation( GridAtlas< GeometryReal , MatrixReal > &fineAtlas , const GridAtlas< GeometryReal , MatrixReal > &coarseAtlas )
{
	std::vector<ProlongationLine> & prolongationLines = fineAtlas.prolongationLines;

	for (int k = 0; k < fineAtlas.gridCharts.size(); k++)
	{
		const GridChart< GeometryReal > &fineChart = fineAtlas.gridCharts[k];
		const GridChart< GeometryReal > &coarseChart = coarseAtlas.gridCharts[k];
		unsigned int width = fineChart.texelIndices.res(0);
		for( unsigned int j=0 ; j<fineChart.texelIndices.res(1) ; j++ )
		{
			unsigned int offset = 0;
			bool previousIsValid = false;
			unsigned int rasterStart = -1;
			while (offset < width) {
				bool currentIsValid = (fineChart.texelIndices( offset , j).combined != -1);
				if (currentIsValid && !previousIsValid) rasterStart = offset; //Start raster line
				if ((!currentIsValid && previousIsValid) || (currentIsValid && offset == (width - 1))) { //Terminate raster line
					if (currentIsValid && offset == (width - 1)) offset++;
					ProlongationLine newLine;
					newLine.startIndex = fineChart.texelIndices( rasterStart , j).combined;
					newLine.length = offset - rasterStart;

					int fi = rasterStart;
					int fj = j;

					int _fi = fi - fineChart.centerOffset[0];
					int _di = ((_fi % 2) == 0) ? 0 : 1;
					int _fj = fj - fineChart.centerOffset[1];
					int _dj = ((_fj % 2) == 0) ? 0 : 1;
					int ci = (_fi - _di) / 2 + coarseChart.centerOffset[0];
					int cj = (_fj - _dj) / 2 + coarseChart.centerOffset[1];
					//int sign_fj = _fj < 0 ? -1 : (_fj > 0 ? 1 : 0);

					newLine.alignedStart = ((_fi % 2) == 0);

					if (_fj % 2 == 0) {
						newLine.prevLineIndex = -1;
						newLine.nextLineIndex = -1;
						int globalIndex = coarseChart.texelIndices(ci, cj).combined;
						if( globalIndex==-1 ) MK_THROW( "Coarse texel is inactive!(A)" );
						else newLine.centerLineIndex = globalIndex;
					}
					else
					{
						newLine.centerLineIndex = -1;

						int globalIndex = coarseChart.texelIndices(ci, cj).combined;
						if( globalIndex==-1 ) MK_THROW( "Coarse texel is inactive!(B)" );
						else newLine.prevLineIndex = globalIndex;

						globalIndex = coarseChart.texelIndices(ci, cj + 1).combined;

						if( globalIndex==-1 ) MK_THROW( "Coarse texel is inactive!(C)" );
						else newLine.nextLineIndex = globalIndex;
					}
					prolongationLines.push_back(newLine);
				}
				previousIsValid = currentIsValid;
				offset++;
			}
		}
	}

#ifdef SANITY_CHECK
	{
		std::vector< Eigen::Triplet< MatrixReal > > prolongationTriplets;

		for (int r = 0; r < prolongationLines.size(); r++) {
			int startIndex = prolongationLines[r].startIndex;
			int lineLength = prolongationLines[r].length;
			int centerLineStart = prolongationLines[r].centerLineIndex;
			int previousLineStart = prolongationLines[r].prevLineIndex;
			int nextLineStart = prolongationLines[r].nextLineIndex;
			int offset = prolongationLines[r].alignedStart ? 0 : 1;
			if (centerLineStart != -1) {
				for (int i = 0; i < lineLength; i++)
				{
					prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( startIndex+i , centerLineStart+ (i+offset+0) / 2 , (MatrixReal)0.5 ) );
					prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( startIndex+i , centerLineStart+ (i+offset+1) / 2 , (MatrixReal)0.5 ) );
				}
			}
			else {
				for (int i = 0; i < lineLength; i++)
				{
					prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( startIndex+i , previousLineStart + (i+offset+0)/2 , (MatrixReal)0.25 ) );
					prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( startIndex+i , previousLineStart + (i+offset+1)/2 , (MatrixReal)0.25 ) );
					prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( startIndex+i ,     nextLineStart + (i+offset+0)/2 , (MatrixReal)0.25 ) );
					prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( startIndex+i ,     nextLineStart + (i+offset+1)/2 , (MatrixReal)0.25 ) );
				}
			}
		}

		Eigen::SparseMatrix< MatrixReal > prolongation;
		prolongation.resize(fineAtlas.numTexels, coarseAtlas.numTexels);
		prolongation.setFromTriplets(prolongationTriplets.begin(), prolongationTriplets.end());
		typedef Eigen::Matrix< MatrixReal , Eigen::Dynamic , 1 > EVector;
		EVector ones = EVector::Ones( coarseAtlas.numTexels );
		EVector prolongedOnes = prolongation*ones;
		for( unsigned int i=0 ; i<fineAtlas.numTexels ; i++ ) if( fabs( prolongedOnes[i]-1.0 )>1e-10 ) MK_THROW( "Prolongation does not add up to one! " , i , " -> " , prolongedOnes[i] );
	}
#endif // SANITY_CHECK
}

template< typename GeometryReal , typename MatrixReal >
void InitializeProlongationMatrix( const GridAtlas< GeometryReal , MatrixReal > &fineAtlas , GridAtlas< GeometryReal , MatrixReal > &coarseAtlas , SparseMatrix< MatrixReal , int > &__prolongation )
{
	const std::vector<ProlongationLine> & prolongationLines = fineAtlas.prolongationLines;

	std::vector< Eigen::Triplet< MatrixReal > > prolongationTriplets;

	for (int r = 0; r < prolongationLines.size(); r++) {
		int startIndex = prolongationLines[r].startIndex;
		int lineLenght = prolongationLines[r].length;
		int centerLineStart = prolongationLines[r].centerLineIndex;
		int previousLineStart = prolongationLines[r].prevLineIndex;
		int nextLineStart = prolongationLines[r].nextLineIndex;
		int offset = prolongationLines[r].alignedStart ? 0 : 1;
		if (centerLineStart != -1) {
			for (int i = 0; i < lineLenght; i++)
			{
				prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( startIndex+i , centerLineStart + (i+offset+0) / 2 , (MatrixReal)0.5 ) );
				prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( startIndex+i , centerLineStart + (i+offset+1) / 2 , (MatrixReal)0.5 ) );
			}
		}
		else {
			for (int i = 0; i < lineLenght; i++)
			{
				prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( startIndex+i , previousLineStart + (i+offset+0) / 2 , (MatrixReal)0.25 ) );
				prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( startIndex+i , previousLineStart + (i+offset+1) / 2 , (MatrixReal)0.25 ) );
				prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( startIndex+i ,     nextLineStart + (i+offset+0) / 2 , (MatrixReal)0.25 ) );
				prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( startIndex+i ,     nextLineStart + (i+offset+1) / 2 , (MatrixReal)0.25 ) );
			}
		}
	}

	__prolongation = SetSparseMatrix( prolongationTriplets , fineAtlas.numTexels , coarseAtlas.numTexels , false );
	ThreadPool::ParallelFor
		(
			0 , __prolongation.rows ,
			[&]( unsigned int , size_t i )
			{
				MatrixReal sum = 0;
				for( int j=0 ; j<__prolongation.rowSizes[i] ; j++ ) sum += __prolongation[i][j].Value;
				if( fabs(sum-1.0)>1e-10 ) MK_THROW( "Prolongation does not add up to one! " , i , " -> " , sum );
			}
		);
}

//Coarse restriction
template< typename GeometryReal , typename MatrixReal  >
void InitializeAtlasHierachicalRestriction
(
	const GridAtlas< GeometryReal , MatrixReal > &fineAtlas ,
	GridAtlas< GeometryReal , MatrixReal > &coarseAtlas ,
	SparseMatrix< MatrixReal , int > &boundaryRestriction
)
{
	std::vector< Eigen::Triplet< MatrixReal > > boundaryRestrictionTriplets;

	const std::vector< GridNodeInfo > & coarseNodeInfo = coarseAtlas.nodeInfo;
	const typename GridAtlas<>::IndexConverter & indexConverter = coarseAtlas.indexConverter;

	const std::vector< RasterLine > &coarseRasterLines = coarseAtlas.rasterLines;
	std::vector<RasterLine> & restrictionLines = coarseAtlas.restrictionLines;
	std::vector<DeepLine> & deepLines = coarseAtlas.deepLines;
	//Initialize restrictionLines

	restrictionLines.resize(coarseRasterLines.size());
	deepLines.resize(coarseRasterLines.size());

	ThreadPool::ParallelFor
		(
			0 , coarseRasterLines.size() ,
			[&]( unsigned int , size_t i )
			{
				int lineCoarseStartIndex = coarseRasterLines[i].lineStartIndex;
				int lineCoarseLength = coarseRasterLines[i].lineEndIndex - lineCoarseStartIndex + 1;

				restrictionLines[i].coeffStartIndex = lineCoarseStartIndex; //global (NOT DEEP) variable index in the current level

				deepLines[i].coarseLineStartIndex = coarseRasterLines[i].coeffStartIndex;
				deepLines[i].coarseLineEndIndex = coarseRasterLines[i].coeffStartIndex + lineCoarseLength - 1;

				GridNodeInfo startNodeInfo = coarseNodeInfo[lineCoarseStartIndex];
				if( startNodeInfo.texelType!=TexelType::InteriorSupported ) MK_THROW( "Not an interior-supported teel" );

				int ci = startNodeInfo.ci;
				int cj = startNodeInfo.cj;
				int chartID = startNodeInfo.chartID;

				const GridChart< GeometryReal > &fineChart = fineAtlas.gridCharts[chartID];
				const GridChart< GeometryReal > &coarseChart = coarseAtlas.gridCharts[chartID];

				int fi = (int)( fineChart.centerOffset[0] + 2.0 *(ci - coarseChart.centerOffset[0]) );
				int fj = (int)( fineChart.centerOffset[1] + 2.0 *(cj - coarseChart.centerOffset[1]) );
				if( fi-1 < 0 || fi+1>(int)fineChart.width-1 || fj-1 < 0 || fj+1>(int)fineChart.height-1 ) MK_THROW( "Out of bounds node position" );

				int fineCurrentLineStart = fineChart.texelIndices(fi, fj).combined;
				if (fineCurrentLineStart != -1) {
					restrictionLines[i].lineStartIndex = fineCurrentLineStart;
					restrictionLines[i].lineEndIndex = fineCurrentLineStart + 2 * (coarseRasterLines[i].lineEndIndex - lineCoarseStartIndex);

					int fineCurrentDeep = fineChart.texelIndices(fi, fj).interior;
					if (fineCurrentDeep != -1) {
						deepLines[i].fineCurrentLineIndex = fineCurrentDeep;
					}
					else MK_THROW( "Invalid fine line start index" );
				}
				else MK_THROW( "Invalid fine line start index" );

				int finePreviousLineStart = fineChart.texelIndices(fi, fj - 1).combined;
				if (finePreviousLineStart != -1) {
					restrictionLines[i].prevLineIndex = finePreviousLineStart;

					int finePreviousDeep = fineChart.texelIndices(fi, fj - 1).interior;
					if (finePreviousDeep != -1) {
						deepLines[i].finePrevLineIndex = finePreviousDeep;
					}
					else MK_THROW( "Invalid fine line start index" );
				}
				else MK_THROW( "Invalid fine previous line start index" );


				int fineNextLineStart = fineChart.texelIndices(fi, fj + 1).combined;
				if (fineNextLineStart != -1) {
					restrictionLines[i].nextLineIndex = fineNextLineStart;

					int fineNextDeep = fineChart.texelIndices(fi, fj + 1).interior;
					if (fineNextDeep != -1) {
						deepLines[i].fineNextLineIndex = fineNextDeep;
					}
					else MK_THROW( "Invalid fine line start index" );
				}
				else MK_THROW( "Invalid fine next line start index" );
			}
		);

	//Initialize boundary nodes prolongation

	for( unsigned int k=0 ; k<fineAtlas.gridCharts.size() ; k++ )
	{
		const GridChart< GeometryReal > &fineChart = fineAtlas.gridCharts[k];
		const GridChart< GeometryReal > &coarseChart = coarseAtlas.gridCharts[k];

		for( unsigned int fj=0 ; fj<fineChart.height ; fj++ ) for( unsigned int fi=0 ; fi<fineChart.width ; fi++ ) if( fineChart.texelIndices(fi,fj).combined!=-1 )
		{

			int _fi = fi - fineChart.centerOffset[0];
			int sign_fi = _fi < 0 ? -1 : (_fi > 0 ? 1 : 0);
			int _fj = fj - fineChart.centerOffset[1];
			int sign_fj = _fj < 0 ? -1 : (_fj > 0 ? 1 : 0);

			int ci[2];
			MatrixReal ci_weights[2];
			if (_fi % 2 == 0) {
				ci[0] = _fi / 2 + coarseChart.centerOffset[0];
				ci[1] = -1;
				ci_weights[0] = 1.0;
				ci_weights[1] = 0.0;
			}
			else {
				ci[0] = _fi / 2 + coarseChart.centerOffset[0];
				ci[1] = (_fi + sign_fi) / 2 + coarseChart.centerOffset[0];
				ci_weights[0] = 0.5;
				ci_weights[1] = 0.5;
			}

			int cj[2];
			MatrixReal cj_weights[2];
			if (_fj % 2 == 0)
			{
				cj[0] = _fj / 2 + coarseChart.centerOffset[1];
				cj[1] = -1;
				cj_weights[0] = 1.0;
				cj_weights[1] = 0.0;
			}
			else 
			{
				cj[0] = _fj / 2 + coarseChart.centerOffset[1];
				cj[1] = (_fj + sign_fj) / 2 + coarseChart.centerOffset[1];
				cj_weights[0] = 0.5;
				cj_weights[1] = 0.5;
			}
			for( int di=0 ; di<2 ; di++ ) for( int dj=0 ; dj<2 ; dj++ ) if( ci[di]!=-1 && cj[dj]!=-1)
			{
				unsigned int coarseNodeGlobalIndex = coarseChart.texelIndices(ci[di], cj[dj]).combined;
				if( coarseNodeGlobalIndex==-1 ) MK_THROW( "Coarse texel is unactive! (D)" );
				else
				{
					unsigned int coarseBoundaryIndex = indexConverter.supportedToBoundary( coarseNodeGlobalIndex );
					if( coarseBoundaryIndex!=-1 ) boundaryRestrictionTriplets.emplace_back( coarseBoundaryIndex , fineChart.texelIndices(fi,fj).combined , ci_weights[di] * cj_weights[dj] );
				}
			}
		}
	}

	boundaryRestriction = SetSparseMatrix( boundaryRestrictionTriplets , coarseAtlas.numBoundaryTexels , fineAtlas.numTexels , false );
}
