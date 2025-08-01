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
void InitializeAtlasHierachicalProlongation( GridAtlas< GeometryReal , MatrixReal > &fineAtlas , const GridAtlas< GeometryReal , MatrixReal > &coarseAtlas )
{
	std::vector<ProlongationLine> & prolongationLines = fineAtlas.prolongationLines;

	for( unsigned int k=0 ; k<fineAtlas.gridCharts.size() ; k++ )
	{
		const GridChart< GeometryReal > &fineChart = fineAtlas.gridCharts[ ChartIndex(k) ];
		const GridChart< GeometryReal > &coarseChart = coarseAtlas.gridCharts[ ChartIndex(k) ];
		unsigned int width = fineChart.texelIndices.res(0);
		for( unsigned int j=0 ; j<fineChart.texelIndices.res(1) ; j++ )
		{
			unsigned int offset = 0;
			bool previousIsValid = false;
			unsigned int rasterStart = -1;
			while( offset<width )
			{
				bool currentIsValid = fineChart.texelIndices( offset , j ).combined!=AtlasTexelIndex(-1);
				if( currentIsValid && !previousIsValid ) rasterStart = offset; //Start raster line
				if( ( !currentIsValid && previousIsValid)  || ( currentIsValid && offset==(width-1) ) )
				{ //Terminate raster line
					if( currentIsValid && offset==(width-1) ) offset++;
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

					if( _fj%2==0 )
					{
						newLine.prevLineIndex = AtlasTexelIndex(-1);
						newLine.nextLineIndex = AtlasTexelIndex(-1);
						AtlasTexelIndex combinedIndex = coarseChart.texelIndices(ci,cj).combined;
#ifdef SANITY_CHECK
						if( combinedIndex==AtlasTexelIndex(-1) ) MK_THROW( "Coarse texel is inactive!(A)" );
						else newLine.centerLineIndex = combinedIndex;
#else // !SANITY_CHECK
					newLine.centerLineIndex = combinedIndex;
#endif // SANITY_CHECK
					}
					else
					{
						newLine.centerLineIndex = AtlasTexelIndex(-1);

						AtlasTexelIndex combinedIndex = coarseChart.texelIndices(ci, cj).combined;
#ifdef SANITY_CHECK
						if( combinedIndex==AtlasTexelIndex(-1) ) MK_THROW( "Coarse texel is inactive!(B)" );
						else newLine.prevLineIndex = combinedIndex;
#else // !SANITY_CHECK
						newLine.prevLineIndex = combinedIndex;
#endif // SANITY_CHECK

						combinedIndex = coarseChart.texelIndices(ci, cj + 1).combined;
#ifdef SANITY_CHECK
						if( combinedIndex==AtlasTexelIndex(-1) ) MK_THROW( "Coarse texel is inactive!(C)" );
						else newLine.nextLineIndex = combinedIndex;
#else // !SANITY_CHECK
						newLine.nextLineIndex = combinedIndex;
#endif // SANITY_CHECK
					}
					prolongationLines.push_back( newLine );
				}
				previousIsValid = currentIsValid;
				offset++;
			}
		}
	}

#ifdef SANITY_CHECK
	{
		std::vector< Eigen::Triplet< MatrixReal > > prolongationTriplets;

		for( int r=0 ; r<prolongationLines.size() ; r++ )
		{
			int lineLength = prolongationLines[r].length;
			AtlasTexelIndex startIndex = prolongationLines[r].startIndex;
			AtlasTexelIndex centerLineStart = prolongationLines[r].centerLineIndex;
			AtlasTexelIndex previousLineStart = prolongationLines[r].prevLineIndex;
			AtlasTexelIndex nextLineStart = prolongationLines[r].nextLineIndex;
			int offset = prolongationLines[r].alignedStart ? 0 : 1;
			if( centerLineStart!=AtlasTexelIndex(-1) )
			{
				for( int i=0 ; i<lineLength ; i++ )
				{
					prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( static_cast< unsigned int >(startIndex) + i , static_cast< unsigned int >(centerLineStart) + (i+offset+0) / 2 , (MatrixReal)0.5 ) );
					prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( static_cast< unsigned int >(startIndex) + i , static_cast< unsigned int >(centerLineStart) + (i+offset+1) / 2 , (MatrixReal)0.5 ) );
				}
			}
			else
			{
				for( int i=0 ; i<lineLength ; i++ )
				{
					prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( static_cast< unsigned int >(startIndex) + i , static_cast< unsigned int >(previousLineStart) + (i+offset+0)/2 , (MatrixReal)0.25 ) );
					prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( static_cast< unsigned int >(startIndex) + i , static_cast< unsigned int >(previousLineStart) + (i+offset+1)/2 , (MatrixReal)0.25 ) );
					prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( static_cast< unsigned int >(startIndex) + i , static_cast< unsigned int >    (nextLineStart) + (i+offset+0)/2 , (MatrixReal)0.25 ) );
					prolongationTriplets.push_back( Eigen::Triplet< MatrixReal >( static_cast< unsigned int >(startIndex) + i , static_cast< unsigned int >    (nextLineStart) + (i+offset+1)/2 , (MatrixReal)0.25 ) );
				}
			}
		}

		Eigen::SparseMatrix< MatrixReal > prolongation;
		prolongation.resize( static_cast< unsigned int >(fineAtlas.endCombinedTexelIndex) , static_cast< unsigned int >(coarseAtlas.endCombinedTexelIndex) );
		prolongation.setFromTriplets( prolongationTriplets.begin() , prolongationTriplets.end() );
		typedef Eigen::Matrix< MatrixReal , Eigen::Dynamic , 1 > EVector;
		EVector ones = EVector::Ones( static_cast< unsigned int >(coarseAtlas.endCombinedTexelIndex) );
		EVector prolongedOnes = prolongation*ones;
		for( unsigned int i=0 ; i<static_cast< unsigned int >(fineAtlas.endCombinedTexelIndex) ; i++ ) if( fabs( prolongedOnes[i]-1.0 )>PROLONGATION_EPSILON ) MK_THROW( "Prolongation does not add up to one! " , i , " -> " , prolongedOnes[i] );
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
#ifdef SANITY_CHECK
	ThreadPool::ParallelFor
		(
			0 , __prolongation.rows ,
			[&]( unsigned int , size_t i )
			{
				MatrixReal sum = 0;
				for( int j=0 ; j<__prolongation.rowSizes[i] ; j++ ) sum += __prolongation[i][j].Value;
				if( fabs(sum-1.0)>PROLONGATION_EPSILON ) MK_THROW( "Prolongation does not add up to one! " , i , " -> " , sum );
			}
		);
#endif // SANITY_CHECK
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

	const ExplicitIndexVector< AtlasTexelIndex , TexelInfo > & coarseTexelInfo = coarseAtlas.texelInfo;
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
				AtlasTexelIndex lineCoarseStartIndex = coarseRasterLines[i].lineStartIndex;
				int lineCoarseLength = coarseRasterLines[i].lineEndIndex - lineCoarseStartIndex + 1;

#ifdef SANITY_CHECK
#pragma message( "[WARNING] Converting AtlasTexelIndex -> AtlasInteriorTexelIndex" )
#endif // SANITY_CHECK
				restrictionLines[i].coeffStartIndex = AtlasInteriorTexelIndex( lineCoarseStartIndex ); //combined (NOT DEEP) variable index in the current level

				deepLines[i].coarseLineStartIndex = coarseRasterLines[i].coeffStartIndex;
				deepLines[i].coarseLineEndIndex = coarseRasterLines[i].coeffStartIndex + lineCoarseLength - 1;

				TexelInfo startTexelInfo = coarseTexelInfo[lineCoarseStartIndex];
#ifdef SANITY_CHECK
				if( startTexelInfo.texelType!=TexelType::InteriorSupported ) MK_THROW( "Not an interior-supported texel" );
#endif // SANITY_CHECK
				int ci = startTexelInfo.ci;
				int cj = startTexelInfo.cj;

				ChartIndex chartID = startTexelInfo.chartID;
				const GridChart< GeometryReal > &fineChart = fineAtlas.gridCharts[ chartID ];
				const GridChart< GeometryReal > &coarseChart = coarseAtlas.gridCharts[ chartID ];

				int fi = (int)( fineChart.centerOffset[0] + 2.0 *(ci - coarseChart.centerOffset[0]) );
				int fj = (int)( fineChart.centerOffset[1] + 2.0 *(cj - coarseChart.centerOffset[1]) );
#ifdef SANITY_CHECK
				if( fi-1 < 0 || fi+1>(int)fineChart.width-1 || fj-1 < 0 || fj+1>(int)fineChart.height-1 ) MK_THROW( "Out of bounds node position" );
#endif // SANITY_CHECK

				AtlasTexelIndex fineCurrentLineStart = fineChart.texelIndices(fi,fj).combined;
				if( fineCurrentLineStart!=AtlasTexelIndex(-1) )
				{
					restrictionLines[i].lineStartIndex = fineCurrentLineStart;
					restrictionLines[i].lineEndIndex = fineCurrentLineStart + 2 * (coarseRasterLines[i].lineEndIndex - lineCoarseStartIndex);

					AtlasInteriorTexelIndex fineCurrentDeep = fineChart.texelIndices(fi,fj).interior;
#ifdef SANITY_CHECK
					if( fineCurrentDeep!=AtlasInteriorTexelIndex(-1) ) deepLines[i].fineCurrentLineIndex = fineCurrentDeep;
					else MK_THROW( "Invalid fine line start index" );
#else // !SANITY_CHECK
					deepLines[i].fineCurrentLineIndex = fineCurrentDeep;
#endif // SANITY_CHECK
				}
#ifdef SANITY_CHECK
				else MK_THROW( "Invalid fine line start index" );
#endif // SANITY_CHECK

				AtlasTexelIndex finePreviousLineStart = fineChart.texelIndices(fi,fj-1).combined;
				if( finePreviousLineStart!=AtlasTexelIndex(-1) )
				{
					restrictionLines[i].prevLineIndex = finePreviousLineStart;

					AtlasInteriorTexelIndex finePreviousDeep = fineChart.texelIndices(fi,fj-1).interior;
#ifdef SANITY_CHECK
					if( finePreviousDeep!=AtlasInteriorTexelIndex(-1) ) deepLines[i].finePrevLineIndex = finePreviousDeep;
					else MK_THROW( "Invalid fine line start index" );
#else // !SANITY_CHECK
					deepLines[i].finePrevLineIndex = finePreviousDeep;
#endif // SANITY_CHECK
				}
#ifdef SANITY_CHECK
				else MK_THROW( "Invalid fine previous line start index" );
#endif // SANITY_CHECK


				AtlasTexelIndex fineNextLineStart = fineChart.texelIndices(fi,fj+1).combined;
				if( fineNextLineStart!=AtlasTexelIndex(-1) )
				{
					restrictionLines[i].nextLineIndex = fineNextLineStart;

					AtlasInteriorTexelIndex fineNextDeep = fineChart.texelIndices(fi,fj+1).interior;
#ifdef SANITY_CHECK
					if( fineNextDeep!=AtlasInteriorTexelIndex(-1) ) deepLines[i].fineNextLineIndex = fineNextDeep;
					else MK_THROW( "Invalid fine line start index" );
#else // !SANITY_CHECK
					deepLines[i].fineNextLineIndex = fineNextDeep;
#endif // SANITY_CHECK
				}
#ifdef SANITY_CHECK
				else MK_THROW( "Invalid fine next line start index" );
#endif // SANITY_CHECK
			}
		);

	//Initialize boundary nodes prolongation

	for( unsigned int k=0 ; k<fineAtlas.gridCharts.size() ; k++ )
	{
		const GridChart< GeometryReal > &fineChart = fineAtlas.gridCharts[ ChartIndex(k) ];
		const GridChart< GeometryReal > &coarseChart = coarseAtlas.gridCharts[ ChartIndex(k) ];

		for( unsigned int fj=0 ; fj<fineChart.height ; fj++ ) for( unsigned int fi=0 ; fi<fineChart.width ; fi++ ) if( fineChart.texelIndices(fi,fj).combined!=AtlasTexelIndex(-1) )
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
				AtlasTexelIndex coarseCombinedTexel = coarseChart.texelIndices(ci[di], cj[dj]).combined;
#ifdef SANITY_CHECK
				if( coarseCombinedTexel==AtlasTexelIndex(-1) ) MK_THROW( "Coarse texel is unactive! (D)" );
				else
				{
					AtlasBoundaryTexelIndex coarseBoundaryIndex = indexConverter.combinedToBoundary( coarseCombinedTexel );
					if( coarseBoundaryIndex!=AtlasBoundaryTexelIndex(-1) ) boundaryRestrictionTriplets.emplace_back( static_cast< unsigned int >(coarseBoundaryIndex) , static_cast< unsigned int >(fineChart.texelIndices(fi,fj).combined) , ci_weights[di] * cj_weights[dj] );
				}
#else // !SANITY_CHECK
				AtlasBoundaryTexelIndex coarseBoundaryIndex = indexConverter.combinedToBoundary( coarseCombinedTexel );
				if( coarseBoundaryIndex!=AtlasBoundaryTexelIndex(-1) ) boundaryRestrictionTriplets.emplace_back( static_cast< unsigned int >(coarseBoundaryIndex) , static_cast< unsigned int >(fineChart.texelIndices(fi,fj).combined) , ci_weights[di] * cj_weights[dj] );
#endif // SANITY_CHECK
			}
		}
	}

	boundaryRestriction = SetSparseMatrix( boundaryRestrictionTriplets , static_cast< unsigned int >(coarseAtlas.endBoundaryTexelIndex) , static_cast< unsigned int >(fineAtlas.endCombinedTexelIndex) , false );
}
