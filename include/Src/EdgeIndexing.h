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
		std::map< SimplexIndex< 1 , AtlasCombinedTexelIndex > , unsigned int > &boundaryCoarseEdgeIndex ,
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

		for( unsigned int i=0 ; i<gridChart.combinedCellCombinedTexelBilinearElementIndices.size() ; i++ )
		{
			const BilinearElementIndex< AtlasCombinedTexelIndex > & indices = gridChart.combinedCellCombinedTexelBilinearElementIndices[ ChartCombinedCellIndex(i) ];
			for( int k=0 ; k<edgesPerCell ; k++ )
			{
				AtlasCombinedTexelIndex vIndices[2] = { indices[ pairsToAdd[2*k] ] , indices[ pairsToAdd[2*k+1] ] };
				SimplexIndex< 1 , AtlasCombinedTexelIndex > edgeKey( vIndices[0] , vIndices[1] );
				if( boundaryCoarseEdgeIndex.find(edgeKey)==boundaryCoarseEdgeIndex.end() )
				{
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
		const ExplicitIndexVector< ChartIndex , GridChart< GeometryReal > > &gridCharts ,
		std::map< SimplexIndex< 1 , AtlasCombinedTexelIndex > , unsigned int > &boundaryCoarseEdgeIndex
	)
	{
		//Add edges within charts
		unsigned int lastAddedEdgeIndex = 0;
		for( unsigned int i=0 ; i<gridCharts.size() ; i++ ) InitializeIntraChartEdgeIndexing( boundaryCoarseEdgeIndex , gridCharts[ ChartIndex(i) ] , lastAddedEdgeIndex );
	}

	template< typename MatrixReal >
	void InitializeBoundaryEdgeIndexing
	(
		const SparseMatrix< MatrixReal , int > &boundaryAdjancencyMatrix ,
		const typename GridAtlas<>::IndexConverter & indexConverter ,
		std::map< SimplexIndex< 1 , AtlasCombinedTexelIndex > , unsigned int > &coarseEdgeIndex ,
		std::vector< unsigned int > &boundaryEdgeToGlobalEdge ,
#ifdef NEW_CODE
		std::map< SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > , AtlasRefinedBoundaryEdgeIndex > &boundaryEdgeIndex
#else // !NEW_CODE
		std::map< SimplexIndex< 1 > , AtlasInteriorOrBoundaryNodeIndex > &boundaryEdgeIndex
#endif // NEW_CODE
	)
	{
		int lastAddedCoarseEdgeIndex = (int)coarseEdgeIndex.size();
#ifdef NEW_CODE
		AtlasRefinedBoundaryEdgeIndex endRefinedBoundaryEdgeIndex(0);
#else // !NEW_CODE
		int lastAddedBoundaryEdgeIndex = 0;
#endif // NEW_CODE
		for( int r=0 ; r<boundaryAdjancencyMatrix.Rows() ; r++ )
		{
			for( int c=0 ; c<boundaryAdjancencyMatrix.RowSize(r) ; c++ )
			{
#ifdef NEW_CODE
				AtlasInteriorOrBoundaryNodeIndex bIndices[] = { AtlasInteriorOrBoundaryNodeIndex(r) , AtlasInteriorOrBoundaryNodeIndex(boundaryAdjancencyMatrix[r][c].N) };
#else // !NEW_CODE
				int bIndices[2] = { r , boundaryAdjancencyMatrix[r][c].N };
#endif // NEW_CODE
				if( bIndices[0]!=bIndices[1] )
				{
#ifdef NEW_CODE
					SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > bminKey( bIndices[0] , bIndices[1] );
					SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > bmaxKey( bIndices[1] , bIndices[0] );
#else // !NEW_CODE
					SimplexIndex< 1 > bminKey( bIndices[0] , bIndices[1] );
					SimplexIndex< 1 > bmaxKey( bIndices[1] , bIndices[0] );
#endif // NEW_CODE
					if( boundaryEdgeIndex.find(bminKey)==boundaryEdgeIndex.end() && boundaryEdgeIndex.find(bmaxKey)==boundaryEdgeIndex.end() )
					{
						AtlasCombinedTexelIndex gIndices[2] = { indexConverter.boundaryToCombined( AtlasBoundaryTexelIndex( bIndices[0] ) ) , indexConverter.boundaryToCombined( AtlasBoundaryTexelIndex( bIndices[1] ) ) };
						SimplexIndex< 1 , AtlasCombinedTexelIndex > minKey( gIndices[0] , gIndices[1] );
						SimplexIndex< 1 , AtlasCombinedTexelIndex > maxKey( gIndices[1] , gIndices[0] );
						int globalEdgeIndex = -1;
						if( coarseEdgeIndex.find(minKey)!=coarseEdgeIndex.end() )
						{
							globalEdgeIndex = coarseEdgeIndex[minKey];
#ifdef NEW_CODE
							boundaryEdgeIndex[bminKey] = endRefinedBoundaryEdgeIndex++;
#else // !NEW_CODE
							boundaryEdgeIndex[bminKey] = static_cast< AtlasInteriorOrBoundaryNodeIndex >( lastAddedBoundaryEdgeIndex );
							lastAddedBoundaryEdgeIndex++;
#endif // NEW_CODE
						}
						else if( coarseEdgeIndex.find(maxKey)!=coarseEdgeIndex.end() )
						{
							globalEdgeIndex = coarseEdgeIndex[maxKey];
#ifdef NEW_CODE
							boundaryEdgeIndex[bmaxKey] = endRefinedBoundaryEdgeIndex++;
#else // !NEW_CODE
							boundaryEdgeIndex[bmaxKey] = static_cast< AtlasInteriorOrBoundaryNodeIndex >( lastAddedBoundaryEdgeIndex );
							lastAddedBoundaryEdgeIndex++;
#endif // NEW_CODE
						}
						else
						{
							coarseEdgeIndex[minKey] = lastAddedCoarseEdgeIndex;
							globalEdgeIndex = lastAddedCoarseEdgeIndex;
							lastAddedCoarseEdgeIndex++;

#ifdef NEW_CODE
							boundaryEdgeIndex[bminKey] = endRefinedBoundaryEdgeIndex++;
#else // !NEW_CODE
							boundaryEdgeIndex[bminKey] = static_cast< AtlasInteriorOrBoundaryNodeIndex >( lastAddedBoundaryEdgeIndex );
							lastAddedBoundaryEdgeIndex++;
#endif // NEW_CODE
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
		const std::vector< AtlasInteriorOrBoundaryNodeIndex > &fineBoundaryNodeIndex ,
#ifdef NEW_CODE
		std::map< SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > , AtlasRefinedBoundaryEdgeIndex > &fineBoundaryEdgeIndex ,
#else // !NEW_CODE
		std::map< SimplexIndex< 1 > , AtlasInteriorOrBoundaryNodeIndex > &fineBoundaryEdgeIndex ,
#endif // NEW_CODE
		const GridChart< GeometryReal > &gridChart ,
#ifdef NEW_CODE
		AtlasRefinedBoundaryEdgeIndex & endRefinedBoundaryEdgeIndex
#else // !NEW_CODE
		unsigned int &lastAddedEdgeIndex
#endif // NEW_CODE
	)
	{
		for( unsigned int c=0 ; c<gridChart.boundaryTriangles.size() ; c++ )
		{
			const std::vector< BoundaryIndexedTriangle< GeometryReal > > & boundaryTriangles = gridChart.boundaryTriangles[ ChartBoundaryCellIndex(c) ];
			for( unsigned int b=0 ; b<boundaryTriangles.size() ; b++ )
			{
				const QuadraticElement::Index & indices = boundaryTriangles[b].indices;
				AtlasInteriorOrBoundaryNodeIndex fineHexIndices[6];
				for( unsigned int k=0 ; k<6 ; k++ ) fineHexIndices[k] = fineBoundaryNodeIndex[ static_cast< unsigned int >(indices[k]) ];
				for( unsigned int k=1 ; k<6 ; k++ ) for( unsigned int l=0 ; l<k ; l++ )
				{
					AtlasInteriorOrBoundaryNodeIndex vIndices[2] = { fineHexIndices[k] , fineHexIndices[l] };
					if( vIndices[1]<vIndices[0] ) std::swap( vIndices[0] , vIndices[1] );
#ifdef NEW_CODE
					SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > edgeKey( vIndices[0] , vIndices[1] );
					if( fineBoundaryEdgeIndex.find(edgeKey)==fineBoundaryEdgeIndex.end() ) fineBoundaryEdgeIndex[edgeKey] = endRefinedBoundaryEdgeIndex++;
#else // !NEW_CODE
					SimplexIndex< 1 > edgeKey( vIndices[0] , vIndices[1] );
					if( fineBoundaryEdgeIndex.find(edgeKey)==fineBoundaryEdgeIndex.end() ) fineBoundaryEdgeIndex[edgeKey] = static_cast< AtlasInteriorOrBoundaryNodeIndex >( lastAddedEdgeIndex++ );
#endif // NEW_CODE
				}
			}
		}
	}

	// Create the index for all edges between element nodes
	// -- Required for finite-differencing
	template< typename GeometryReal >
	void InitializeFineBoundaryEdgeIndexing
	(
		const std::vector< AtlasInteriorOrBoundaryNodeIndex > &fineBoundaryNodeIndex ,
#ifdef NEW_CODE
		std::map< SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > , AtlasRefinedBoundaryEdgeIndex > &fineBoundaryEdgeIndex ,
#else // !NEW_CODE
		std::map< SimplexIndex< 1 > , AtlasInteriorOrBoundaryNodeIndex > &fineBoundaryEdgeIndex ,
#endif // NEW_CODE
		const ExplicitIndexVector< ChartIndex , GridChart< GeometryReal > > &gridCharts
	)
	{
#ifdef NEW_CODE
		AtlasRefinedBoundaryEdgeIndex endRefinedBoundaryEdgeIndex(0);
		for( unsigned int i=0 ; i<gridCharts.size() ; i++ ) InitializeFineBoundaryEdgeChartIndexing( fineBoundaryNodeIndex , fineBoundaryEdgeIndex , gridCharts[ ChartIndex(i) ] , endRefinedBoundaryEdgeIndex );
#else // !NEW_CODE
		unsigned int lastAddedEdgeIndex = 0;
		for( unsigned int i=0 ; i<gridCharts.size() ; i++ ) InitializeFineBoundaryEdgeChartIndexing( fineBoundaryNodeIndex , fineBoundaryEdgeIndex , gridCharts[ ChartIndex(i) ] , lastAddedEdgeIndex );
#endif // NEW_CODE
	}

	template< typename MatrixReal >
	void InitializeBoundaryCoarseToFineBoundaryOneFormProlongation
	(
		const SparseMatrix< MatrixReal , int > &boundaryCoarseToFineNodeProlongation ,
#ifdef NEW_CODE
		std::map< SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > , AtlasRefinedBoundaryEdgeIndex > &boundaryCoarseEdgeIndex ,
		std::map< SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > , AtlasRefinedBoundaryEdgeIndex > &boundaryFineEdgeIndex ,
#else // !NEW_CODE
		std::map< SimplexIndex< 1 > , AtlasInteriorOrBoundaryNodeIndex > &boundaryCoarseEdgeIndex ,
		std::map< SimplexIndex< 1 > , AtlasInteriorOrBoundaryNodeIndex > &boundaryFineEdgeIndex ,
#endif // NEW_CODE
		SparseMatrix< MatrixReal , int > &boundaryFineToBoundaryCoarseOneFormProlongation
	)
	{
		std::vector< Eigen::Triplet< MatrixReal > > coarseToFineOneFormProlongation;
		std::vector< std::vector< Eigen::Triplet< MatrixReal > > > _coarseToFineOneFormProlongation( ThreadPool::NumThreads() );
#ifdef NEW_CODE
		std::vector< std::pair< SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > , AtlasRefinedBoundaryEdgeIndex > > _boundaryFineEdgeIndex;
#else // !NEW_CODE
		std::vector< std::pair< SimplexIndex< 1 > , AtlasInteriorOrBoundaryNodeIndex > > _boundaryFineEdgeIndex;
#endif // NEW_CODE

		// Transform the unordered_map into a vector of pairs for parallelization
		_boundaryFineEdgeIndex.reserve( boundaryFineEdgeIndex.size() );
#ifdef NEW_CODE
		for( auto iter=boundaryFineEdgeIndex.begin() ; iter!=boundaryFineEdgeIndex.end() ; iter++ ) _boundaryFineEdgeIndex.push_back( std::make_pair( iter->first , iter->second ) );
#else // !NEW_CODE
		for( auto iter=boundaryFineEdgeIndex.begin() ; iter!=boundaryFineEdgeIndex.end() ; iter++ ) _boundaryFineEdgeIndex.push_back( std::pair< SimplexIndex< 1 > , AtlasInteriorOrBoundaryNodeIndex >( iter->first , iter->second ) );
#endif // NEW_CODE

		ThreadPool::ParallelFor
		(
			0 , _boundaryFineEdgeIndex.size() ,
			[&]( unsigned int thread , size_t i )
			{
#ifdef NEW_CODE
				AtlasRefinedBoundaryEdgeIndex fineEdgeId = _boundaryFineEdgeIndex[i].second;
				SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > fineEdgeCorners = _boundaryFineEdgeIndex[i].first;

				for( unsigned int k=0 ; k<boundaryCoarseToFineNodeProlongation.RowSize( static_cast< unsigned int >(fineEdgeCorners[0]) ) ; k++ )
				{
					AtlasInteriorOrBoundaryNodeIndex coarseIndex1( boundaryCoarseToFineNodeProlongation[ static_cast< unsigned int >(fineEdgeCorners[0]) ][k].N );
					MatrixReal coarseValue1 = boundaryCoarseToFineNodeProlongation[ static_cast< unsigned int >(fineEdgeCorners[0]) ][k].Value;

					for( unsigned int l=0 ; l<boundaryCoarseToFineNodeProlongation.RowSize( static_cast< unsigned int >(fineEdgeCorners[1]) ) ; l++ )
					{
						AtlasInteriorOrBoundaryNodeIndex coarseIndex2( boundaryCoarseToFineNodeProlongation[ static_cast< unsigned int >(fineEdgeCorners[1]) ][l].N );
						MatrixReal coarseValue2 = boundaryCoarseToFineNodeProlongation[ static_cast< unsigned int >(fineEdgeCorners[1]) ][l].Value;
#else // !NEW_CODE
				AtlasInteriorOrBoundaryNodeIndex fineEdgeId = _boundaryFineEdgeIndex[i].second;
				SimplexIndex< 1 > fineEdgeCorners = _boundaryFineEdgeIndex[i].first;

				for( unsigned int k=0 ; k<boundaryCoarseToFineNodeProlongation.RowSize( fineEdgeCorners[0] ) ; k++ )
				{
					int coarseIndex1 = boundaryCoarseToFineNodeProlongation[ fineEdgeCorners[0] ][k].N;
					MatrixReal coarseValue1 = boundaryCoarseToFineNodeProlongation[ fineEdgeCorners[0] ][k].Value;

					for( unsigned int l=0 ; l<boundaryCoarseToFineNodeProlongation.RowSize( fineEdgeCorners[1] ) ; l++ )
					{
						int coarseIndex2 = boundaryCoarseToFineNodeProlongation[ fineEdgeCorners[1] ][l].N;
						MatrixReal coarseValue2 = boundaryCoarseToFineNodeProlongation[ fineEdgeCorners[1] ][l].Value;
#endif // NEW_CODE

						if( coarseIndex1!=coarseIndex2 )
						{
							bool foundEdge = false;
#ifdef NEW_CODE
							SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > coarseEdgeKey( coarseIndex1 , coarseIndex2 );
#else // !NEW_CODE
							SimplexIndex< 1 > coarseEdgeKey( coarseIndex1 , coarseIndex2 );
#endif // NEW_CODE
							auto coarseEdgePtr = boundaryCoarseEdgeIndex.find( coarseEdgeKey );
							if( coarseEdgePtr!=boundaryCoarseEdgeIndex.end() )
							{
								foundEdge = true;
#ifdef NEW_CODE
								AtlasRefinedBoundaryEdgeIndex coarseEdgeId = coarseEdgePtr->second;
#else // !NEW_CODE
								AtlasInteriorOrBoundaryNodeIndex coarseEdgeId = coarseEdgePtr->second;
#endif // NEW_CODE
								_coarseToFineOneFormProlongation[thread].push_back( Eigen::Triplet< MatrixReal >( static_cast< unsigned int >(fineEdgeId) , static_cast< unsigned int >(coarseEdgeId) , coarseValue1 * coarseValue2 ) );
							}
							else
							{
#ifdef NEW_CODE
								coarseEdgeKey = SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex >( coarseIndex2 , coarseIndex1 );
#else // !NEW_CODE
								coarseEdgeKey = SimplexIndex< 1 >( coarseIndex2 , coarseIndex1 );
#endif // NEW_CODE
								coarseEdgePtr = boundaryCoarseEdgeIndex.find(coarseEdgeKey);
								if( coarseEdgePtr!=boundaryCoarseEdgeIndex.end() )
								{
									foundEdge = true;
#ifdef NEW_CODE
									AtlasRefinedBoundaryEdgeIndex coarseEdgeId = coarseEdgePtr->second;
#else // !NEW_CODE
									AtlasInteriorOrBoundaryNodeIndex coarseEdgeId = coarseEdgePtr->second;
#endif // NEW_CODE
									_coarseToFineOneFormProlongation[thread].push_back( Eigen::Triplet< MatrixReal >( static_cast< unsigned int >(fineEdgeId) , static_cast< unsigned int >(coarseEdgeId) , -coarseValue1 *coarseValue2 ) );
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
