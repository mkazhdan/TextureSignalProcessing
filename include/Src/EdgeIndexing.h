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
	namespace TSP
	{
		template< typename GeometryReal >
		void InitializeIntraChartEdgeIndexing
		(
			std::map< SimplexIndex< 1 , AtlasTexelIndex > , unsigned int > &boundaryCoarseEdgeIndex ,
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
				const BilinearElementIndex< AtlasTexelIndex > & indices = gridChart.combinedCellCombinedTexelBilinearElementIndices[ ChartCellIndex(i) ];
				for( int k=0 ; k<edgesPerCell ; k++ )
				{
					AtlasTexelIndex vIndices[2] = { indices[ pairsToAdd[2*k] ] , indices[ pairsToAdd[2*k+1] ] };
					SimplexIndex< 1 , AtlasTexelIndex > edgeKey( vIndices[0] , vIndices[1] );
#ifdef SANITY_CHECK
					if( boundaryCoarseEdgeIndex.find(edgeKey)==boundaryCoarseEdgeIndex.end() )
#endif // SANITY_CHECK
					{
						boundaryCoarseEdgeIndex[edgeKey] = lastAddedEdgeIndex;
						lastAddedEdgeIndex++;
					}
#ifdef SANITY_CHECK
					else MK_THROW( "Edge already added" );
#endif // SANITY_CHECK
				}
			}
		}

		template< typename GeometryReal >
		void InitializeIntraChartEdgeIndexing
		(
			const ExplicitIndexVector< ChartIndex , GridChart< GeometryReal > > &gridCharts ,
			std::map< SimplexIndex< 1 , AtlasTexelIndex > , unsigned int > &boundaryCoarseEdgeIndex
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
			std::map< SimplexIndex< 1 , AtlasTexelIndex > , unsigned int > &coarseEdgeIndex ,
			std::vector< unsigned int > &boundaryEdgeToGlobalEdge ,
			std::map< SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > , AtlasRefinedBoundaryEdgeIndex > &boundaryEdgeIndex
		)
		{
			int lastAddedCoarseEdgeIndex = (int)coarseEdgeIndex.size();
			AtlasRefinedBoundaryEdgeIndex endRefinedBoundaryEdgeIndex(0);
			for( int r=0 ; r<boundaryAdjancencyMatrix.Rows() ; r++ )
			{
				for( int c=0 ; c<boundaryAdjancencyMatrix.RowSize(r) ; c++ )
				{
					AtlasInteriorOrBoundaryNodeIndex bIndices[] = { AtlasInteriorOrBoundaryNodeIndex(r) , AtlasInteriorOrBoundaryNodeIndex(boundaryAdjancencyMatrix[r][c].N) };
					if( bIndices[0]!=bIndices[1] )
					{
						SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > bminKey( bIndices[0] , bIndices[1] );
						SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > bmaxKey( bIndices[1] , bIndices[0] );
						if( boundaryEdgeIndex.find(bminKey)==boundaryEdgeIndex.end() && boundaryEdgeIndex.find(bmaxKey)==boundaryEdgeIndex.end() )
						{
							AtlasTexelIndex gIndices[2] = { indexConverter.boundaryToCombined( AtlasBoundaryTexelIndex( bIndices[0] ) ) , indexConverter.boundaryToCombined( AtlasBoundaryTexelIndex( bIndices[1] ) ) };
							SimplexIndex< 1 , AtlasTexelIndex > minKey( gIndices[0] , gIndices[1] );
							SimplexIndex< 1 , AtlasTexelIndex > maxKey( gIndices[1] , gIndices[0] );
							int globalEdgeIndex = -1;
							if( coarseEdgeIndex.find(minKey)!=coarseEdgeIndex.end() )
							{
								globalEdgeIndex = coarseEdgeIndex[minKey];
								boundaryEdgeIndex[bminKey] = endRefinedBoundaryEdgeIndex++;
							}
							else if( coarseEdgeIndex.find(maxKey)!=coarseEdgeIndex.end() )
							{
								globalEdgeIndex = coarseEdgeIndex[maxKey];
								boundaryEdgeIndex[bmaxKey] = endRefinedBoundaryEdgeIndex++;
							}
							else
							{
								coarseEdgeIndex[minKey] = lastAddedCoarseEdgeIndex;
								globalEdgeIndex = lastAddedCoarseEdgeIndex;
								lastAddedCoarseEdgeIndex++;

								boundaryEdgeIndex[bminKey] = endRefinedBoundaryEdgeIndex++;
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
			std::map< SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > , AtlasRefinedBoundaryEdgeIndex > &fineBoundaryEdgeIndex ,
			const GridChart< GeometryReal > &gridChart ,
			AtlasRefinedBoundaryEdgeIndex & endRefinedBoundaryEdgeIndex
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
						SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > edgeKey( vIndices[0] , vIndices[1] );
						if( fineBoundaryEdgeIndex.find(edgeKey)==fineBoundaryEdgeIndex.end() ) fineBoundaryEdgeIndex[edgeKey] = endRefinedBoundaryEdgeIndex++;
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
			std::map< SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > , AtlasRefinedBoundaryEdgeIndex > &fineBoundaryEdgeIndex ,
			const ExplicitIndexVector< ChartIndex , GridChart< GeometryReal > > &gridCharts
		)
		{
			AtlasRefinedBoundaryEdgeIndex endRefinedBoundaryEdgeIndex(0);
			for( unsigned int i=0 ; i<gridCharts.size() ; i++ ) InitializeFineBoundaryEdgeChartIndexing( fineBoundaryNodeIndex , fineBoundaryEdgeIndex , gridCharts[ ChartIndex(i) ] , endRefinedBoundaryEdgeIndex );
		}

		template< typename MatrixReal >
		void InitializeBoundaryCoarseToFineBoundaryOneFormProlongation
		(
			const SparseMatrix< MatrixReal , int > &boundaryCoarseToFineNodeProlongation ,
			std::map< SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > , AtlasRefinedBoundaryEdgeIndex > &boundaryCoarseEdgeIndex ,
			std::map< SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > , AtlasRefinedBoundaryEdgeIndex > &boundaryFineEdgeIndex ,
			SparseMatrix< MatrixReal , int > &boundaryFineToBoundaryCoarseOneFormProlongation
		)
		{
			std::vector< Eigen::Triplet< MatrixReal > > coarseToFineOneFormProlongation;
			std::vector< std::vector< Eigen::Triplet< MatrixReal > > > _coarseToFineOneFormProlongation( ThreadPool::NumThreads() );
			std::vector< std::pair< SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > , AtlasRefinedBoundaryEdgeIndex > > _boundaryFineEdgeIndex;

			// Transform the map into a vector of pairs for parallelization
			_boundaryFineEdgeIndex.reserve( boundaryFineEdgeIndex.size() );
			for( auto iter=boundaryFineEdgeIndex.begin() ; iter!=boundaryFineEdgeIndex.end() ; iter++ ) _boundaryFineEdgeIndex.push_back( std::make_pair( iter->first , iter->second ) );

			ThreadPool::ParallelFor
			(
				0 , _boundaryFineEdgeIndex.size() ,
				[&]( unsigned int thread , size_t i )
				{
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

							if( coarseIndex1!=coarseIndex2 )
							{
								bool foundEdge = false;
								SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > coarseEdgeKey( coarseIndex1 , coarseIndex2 );
								auto coarseEdgePtr = boundaryCoarseEdgeIndex.find( coarseEdgeKey );
								if( coarseEdgePtr!=boundaryCoarseEdgeIndex.end() )
								{
									foundEdge = true;
									AtlasRefinedBoundaryEdgeIndex coarseEdgeId = coarseEdgePtr->second;
									_coarseToFineOneFormProlongation[thread].push_back( Eigen::Triplet< MatrixReal >( static_cast< unsigned int >(fineEdgeId) , static_cast< unsigned int >(coarseEdgeId) , coarseValue1 * coarseValue2 ) );
								}
								else
								{
									coarseEdgeKey = SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex >( coarseIndex2 , coarseIndex1 );
									coarseEdgePtr = boundaryCoarseEdgeIndex.find(coarseEdgeKey);
									if( coarseEdgePtr!=boundaryCoarseEdgeIndex.end() )
									{
										foundEdge = true;
										AtlasRefinedBoundaryEdgeIndex coarseEdgeId = coarseEdgePtr->second;
										_coarseToFineOneFormProlongation[thread].push_back( Eigen::Triplet< MatrixReal >( static_cast< unsigned int >(fineEdgeId) , static_cast< unsigned int >(coarseEdgeId) , -coarseValue1 *coarseValue2 ) );
									}
								}
#ifdef SANITY_CHECK
								if( !foundEdge ) MK_THROW( "Edge (" , coarseIndex1 , "," , coarseIndex2 , ") not found" );
#endif // SANITY_CHECK
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
}
#endif // EDGE_INDEXING_INCLUDED
