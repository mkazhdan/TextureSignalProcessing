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

#ifndef EDGE_INDEXING_INCLUDED
#define EDGE_INDEXING_INCLUDED

#include <Misha/Miscellany.h>
#include <Misha/MultiThreading.h>

namespace MishaK
{
	template< typename GeometryReal >
	void InitializeIntraChartEdgeIndexing
	(
		std::map< EdgeIndex , unsigned int > &boundaryCoarseEdgeIndex ,
		const GridChart< GeometryReal > &gridChart ,
		unsigned int &lastAddedEdgeIndex
	)
	{	
		int edgesPerCell = 2;
		int pairsToAdd[4] = { 0,1,0,3 };
		// Node indexing
		//		0 ---- 1
		//		|	   |
		//		|	   |
		//		3------2

		// Edge indexing
		//		 ---0---
		//		|  
		//		1  
		//		|  

		for (int i = 0; i < gridChart.combinedCellCombinedBilinearElementIndices.size(); i++){
			const BilinearElementIndex & indices = gridChart.combinedCellCombinedBilinearElementIndices[i];
			for (int k = 0; k < edgesPerCell; k++) {
				int vIndices[2] = { (int)indices[ pairsToAdd[2*k] ] , (int)indices[ pairsToAdd[2*k+1] ] };
				EdgeIndex edgeKey( vIndices[0] , vIndices[1] );
				if (boundaryCoarseEdgeIndex.find(edgeKey) == boundaryCoarseEdgeIndex.end()) {
					boundaryCoarseEdgeIndex[edgeKey] = lastAddedEdgeIndex;
					lastAddedEdgeIndex++;
				}
				else MK_THROW( "Edge already added" );
			}
		}
	}

	template< typename GeometryReal >
	void InitializeIntraChartEdgeIndexing
	(
		const std::vector< GridChart< GeometryReal > > &gridCharts ,
		std::map< EdgeIndex , unsigned int > &boundaryCoarseEdgeIndex
	)
	{
		//Add edges within charts
		unsigned int lastAddedEdgeIndex = 0;
		for( int i=0 ; i<gridCharts.size() ; i++ ) InitializeIntraChartEdgeIndexing( boundaryCoarseEdgeIndex , gridCharts[i] , lastAddedEdgeIndex );
	}

	template< typename MatrixReal >
	void InitializeBoundaryEdgeIndexing
	(
		const SparseMatrix< MatrixReal , int > &boundaryAdjancencyMatrix ,
		const typename GridAtlas<>::IndexConverter & indexConverter ,
		std::map< EdgeIndex , unsigned int > &coarseEdgeIndex ,
		std::vector< unsigned int > &boundaryEdgeToGlobalEdge ,
		std::map< EdgeIndex , unsigned int > &boundaryEdgeIndex
	)
	{
		int lastAddedCoarseEdgeIndex = (int)coarseEdgeIndex.size();
		int lastAddedBoundaryEdgeIndex = 0;
		for (int r = 0; r < boundaryAdjancencyMatrix.Rows(); r++){
			for(int c = 0; c< boundaryAdjancencyMatrix.RowSize(r); c++){
				int bIndices[2] = {r,boundaryAdjancencyMatrix[r][c].N};
				if (bIndices[0] != bIndices[1])
				{
					EdgeIndex bminKey( bIndices[0] , bIndices[1] );
					EdgeIndex bmaxKey( bIndices[1] , bIndices[0] );
					if( boundaryEdgeIndex.find(bminKey)==boundaryEdgeIndex.end() && boundaryEdgeIndex.find(bmaxKey)==boundaryEdgeIndex.end() )
					{
						unsigned int gIndices[2] = { indexConverter.boundaryToSupported( bIndices[0] ) , indexConverter.boundaryToSupported( bIndices[1] ) };
						EdgeIndex minKey( gIndices[0] , gIndices[1] );
						EdgeIndex maxKey( gIndices[1] , gIndices[0] );
						int globalEdgeIndex = -1;
						if( coarseEdgeIndex.find(minKey)!=coarseEdgeIndex.end() )
						{
							globalEdgeIndex = coarseEdgeIndex[minKey];
							boundaryEdgeIndex[bminKey] = lastAddedBoundaryEdgeIndex;
							lastAddedBoundaryEdgeIndex++;
						}
						else if (coarseEdgeIndex.find(maxKey) != coarseEdgeIndex.end()) {
							globalEdgeIndex = coarseEdgeIndex[maxKey];
							boundaryEdgeIndex[bmaxKey] = lastAddedBoundaryEdgeIndex;
							lastAddedBoundaryEdgeIndex++;
						}
						else {
							coarseEdgeIndex[minKey] = lastAddedCoarseEdgeIndex;
							globalEdgeIndex = lastAddedCoarseEdgeIndex;
							lastAddedCoarseEdgeIndex++;

							boundaryEdgeIndex[bminKey] = lastAddedBoundaryEdgeIndex;
							lastAddedBoundaryEdgeIndex++;
						}
						boundaryEdgeToGlobalEdge.push_back(globalEdgeIndex);
					}
				}
			}
		}
	}

	template< typename GeometryReal >
	void InitializeFineBoundaryEdgeChartIndexing
	(
		const std::vector< unsigned int > &fineBoundaryNodeIndex ,
		std::map< EdgeIndex , unsigned int > &fineBoundaryEdgeIndex ,
		const GridChart< GeometryReal > &gridChart ,
		unsigned int &lastAddedEdgeIndex
	)
	{

		for (int c = 0; c < gridChart.boundaryTriangles.size(); c++) {
			for (int b = 0; b < gridChart.boundaryTriangles[c].size(); b++) {
				const QuadraticElementIndex & indices = gridChart.boundaryTriangles[c][b].indices;
				int fineHexIndices[6];
				for (int k = 0; k < 6; k++) fineHexIndices[k] = fineBoundaryNodeIndex[indices[k]];
				for (int k = 1; k < 6; k++) for (int l = 0; l < k; l++) {
					int vIndices[2] = { fineHexIndices[k],fineHexIndices[l] };
					if (vIndices[0] > vIndices[1]) std::swap(vIndices[0], vIndices[1]);
					EdgeIndex edgeKey( vIndices[0] , vIndices[1] );
					if (fineBoundaryEdgeIndex.find(edgeKey) == fineBoundaryEdgeIndex.end()) {
						fineBoundaryEdgeIndex[edgeKey] = lastAddedEdgeIndex;
						lastAddedEdgeIndex++;
					}
				}
			}
		}
	}

	template< typename GeometryReal >
	void InitializeFineBoundaryEdgeIndexing
	(
		const std::vector< unsigned int > &fineBoundaryNodeIndex ,
		std::map< EdgeIndex , unsigned int > &fineBoundaryEdgeIndex ,
		const std::vector< GridChart< GeometryReal > > &gridCharts
	)
	{
		unsigned int lastAddedEdgeIndex = 0;
		for( int i=0 ; i<gridCharts.size() ; i++ ) InitializeFineBoundaryEdgeChartIndexing( fineBoundaryNodeIndex , fineBoundaryEdgeIndex , gridCharts[i] , lastAddedEdgeIndex );
	}

	template< typename MatrixReal >
	void InitializeBoundaryCoarseToFineBoundaryOneFormProlongation
	(
		const SparseMatrix< MatrixReal , int > &boundaryCoarseToFineNodeProlongation ,
		std::map< EdgeIndex , unsigned int > &boundaryCoarseEdgeIndex ,
		std::map< EdgeIndex , unsigned int > &boundaryFineEdgeIndex ,
		SparseMatrix< MatrixReal , int > &boundaryFineToBoundaryCoarseOneFormProlongation
	)
	{
		std::vector< Eigen::Triplet< MatrixReal > > coarseToFineOneFormProlongation;
		std::vector< std::vector< Eigen::Triplet< MatrixReal > > > _coarseToFineOneFormProlongation( ThreadPool::NumThreads() );
		std::vector< std::pair< EdgeIndex , int > > _boundaryFineEdgeIndex;

		// Transform the unordered_map into a vector of pairs for parallelization
		_boundaryFineEdgeIndex.reserve( boundaryFineEdgeIndex.size() );
		for( auto iter=boundaryFineEdgeIndex.begin() ; iter!=boundaryFineEdgeIndex.end() ; iter++ ) _boundaryFineEdgeIndex.push_back( std::pair< EdgeIndex , unsigned int >( iter->first , iter->second ) );

		ThreadPool::ParallelFor
		(
			0 , _boundaryFineEdgeIndex.size() ,
			[&]( unsigned int thread , size_t i )
			{
				unsigned int fineEdgeId = _boundaryFineEdgeIndex[i].second;
				EdgeIndex fineEdgeCorners = _boundaryFineEdgeIndex[i].first;

				for( int k=0 ; k<boundaryCoarseToFineNodeProlongation.RowSize( fineEdgeCorners[0] ) ; k++ )
				{
					int coarseIndex1 = boundaryCoarseToFineNodeProlongation[ fineEdgeCorners[0] ][k].N;
					MatrixReal coarseValue1 = boundaryCoarseToFineNodeProlongation[ fineEdgeCorners[0] ][k].Value;

					for( int l=0 ; l<boundaryCoarseToFineNodeProlongation.RowSize( fineEdgeCorners[1] ) ; l++ )
					{
						int coarseIndex2 = boundaryCoarseToFineNodeProlongation[ fineEdgeCorners[1] ][l].N;
						MatrixReal coarseValue2 = boundaryCoarseToFineNodeProlongation[ fineEdgeCorners[1] ][l].Value;

						if( coarseIndex1!=coarseIndex2 )
						{
							bool foundEdge = false;
							EdgeIndex coarseEdgeKey( coarseIndex1 , coarseIndex2 );
							auto coarseEdgePtr = boundaryCoarseEdgeIndex.find( coarseEdgeKey );
							if( coarseEdgePtr!=boundaryCoarseEdgeIndex.end() )
							{
								foundEdge = true;
								int coarseEdgeId = coarseEdgePtr->second;
								_coarseToFineOneFormProlongation[thread].push_back( Eigen::Triplet< MatrixReal >( fineEdgeId , coarseEdgeId , coarseValue1 * coarseValue2 ) );
							}
							else
							{
								coarseEdgeKey = EdgeIndex( coarseIndex2 , coarseIndex1 );
								coarseEdgePtr = boundaryCoarseEdgeIndex.find(coarseEdgeKey);
								if( coarseEdgePtr!=boundaryCoarseEdgeIndex.end() )
								{
									foundEdge = true;
									int coarseEdgeId = coarseEdgePtr->second;
									_coarseToFineOneFormProlongation[thread].push_back( Eigen::Triplet< MatrixReal >( fineEdgeId , coarseEdgeId , -coarseValue1 *coarseValue2 ) );
								}
							}
							if( !foundEdge ) MK_THROW( "Edge (" , coarseIndex1 , "," , coarseIndex2 , ") not found" );
						}
					}
				}
			}
		);

		// Merge the prolongation entries
		{
			size_t count = 0;
			for( int i=0 ; i<_coarseToFineOneFormProlongation.size() ; i++ ) count += _coarseToFineOneFormProlongation[i].size();
			coarseToFineOneFormProlongation.reserve( count );
			for( int i=0 ; i<_coarseToFineOneFormProlongation.size() ; i++ ) for( int j=0 ; j<_coarseToFineOneFormProlongation[i].size() ; j++ ) coarseToFineOneFormProlongation.push_back( _coarseToFineOneFormProlongation[i][j] );
		}

		int numCoarseOneForms = (int)boundaryCoarseEdgeIndex.size();
		int numFineOneForms = (int)boundaryFineEdgeIndex.size();

		boundaryFineToBoundaryCoarseOneFormProlongation = SetSparseMatrix( coarseToFineOneFormProlongation , numFineOneForms , numCoarseOneForms , false );
	}
}
#endif // EDGE_INDEXING_INCLUDED
