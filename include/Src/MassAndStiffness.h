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
#include "Hierarchy.h"
#include "Divergence.h"
#include "PolygonClipping.h"

namespace MishaK
{
	template< unsigned int Samples , typename GeometryReal , typename MatrixReal >
	void InitializeChartMassAndStiffness
	(
		const std::vector< SquareMatrix< GeometryReal , 2 > > &texture_metrics ,
		const AtlasChart< GeometryReal > &atlasChart ,
		const GridChart< GeometryReal > &gridChart ,
		const std::vector< int >& boundaryAndDeepIndex ,
		const std::vector< int >& fineBoundaryIndex ,
		std::vector< MatrixReal >& deepMassCoefficients ,
		std::vector< MatrixReal >& deepStiffnessCoefficients ,
		std::vector< Eigen::Triplet< MatrixReal > >& boundaryBoundaryMassTriplets ,
		std::vector< Eigen::Triplet< MatrixReal > >& boundaryBoundaryStiffnessTriplets ,
		std::vector< Eigen::Triplet< MatrixReal > >& boundaryDeepMassTriplets ,
		std::vector< Eigen::Triplet< MatrixReal > >& boundaryDeepStiffnessTriplets ,
		bool computeCellBasedStiffness ,
		const std::vector< Point3D< MatrixReal > >& inputSignal ,
		const std::vector< Point3D< MatrixReal > >& boundarySignal ,
		std::vector< MatrixReal >& texelToCellCoeffs ,
		std::vector< Eigen::Triplet< Point3D< MatrixReal > > > &boundaryCellStiffnessTriplets ,
		bool computeDivergence,
		std::unordered_map< unsigned long long , int > & fineBoundaryEdgeIndex,
		std::unordered_map< unsigned long long , int > & coarseEdgeIndex,
		std::vector< Eigen::Triplet< MatrixReal > > & boundaryDeepDivergenceTriplets,
		std::vector< Eigen::Triplet< MatrixReal > > & boundaryBoundaryDivergenceTriplets,
		std::vector< MatrixReal > & deepDivergenceCoefficients
	)
	{
		std::vector< SquareMatrix< GeometryReal , 4 > > cellStiffness;
		std::vector< SquareMatrix< GeometryReal , 4 > > cellMass;
		cellMass.resize( gridChart.numInteriorCells );
		cellStiffness.resize( gridChart.numInteriorCells );

		std::vector< SquareMatrix< GeometryReal , 6 > > triangleElementStiffness;
		std::vector< SquareMatrix< GeometryReal , 6 > > triangleElementMass;
		triangleElementMass.resize( gridChart.numBoundaryTriangles );
		triangleElementStiffness.resize( gridChart.numBoundaryTriangles );

		std::vector< SquareMatrix< GeometryReal , 4 > > cellDivergence;
		if( computeDivergence ) cellDivergence.resize( gridChart.numInteriorCells );
		std::vector< Matrix< GeometryReal , 6 , 15 > > triangleElementDivergence;
		if( computeDivergence ) triangleElementDivergence.resize(gridChart.numBoundaryTriangles);

		//// Rasterize
		int zeroAreaElementCount = 0;
		GeometryReal PRECISION_ERROR = (GeometryReal)1e-3;

		// Node indexing
		//		0 ---- 1
		//		|	   |
		//		|	   |
		//		3------2

		int cellCornerPairs[12] = { 0,1,0,3,0,2,1,3,3,2,1,2 };

		// Cell edge indexing
		//		+---0--+
		//		|\    /|
		//		| 2  / |
		//		|  \/  5
		//		1  /\  |
		//		| 3  \ |
		//		|/    \|
		//      +---4--+

		int reducedCellCornerPairs[8] = { 0,1,0,3,3,2,1,2 };

		// Reduced cell edge indexing
		//		+--0---+
		//		|      |
		//		1      3
		//		|      |
		//		+--2---+


		auto InUnitSquare =   [&]( Point2D< GeometryReal > p ){ return !( p[0]<0-PRECISION_ERROR || p[1]<0-PRECISION_ERROR || p[0]>1+PRECISION_ERROR || p[1]>1+PRECISION_ERROR ); };
		auto InUnitTriangle = [&]( Point2D< GeometryReal > p ){ return !( p[0]<0-PRECISION_ERROR || p[1]<0-PRECISION_ERROR || ( p[0]+p[1] )>1+PRECISION_ERROR ); };
		auto CellInTriangle = [&]( int i , int j , const std::vector< Point2D< GeometryReal > >& vertices )
			{
				GeometryReal x1 = (GeometryReal)i*gridChart.cellSizeW , x2 = (GeometryReal)(i+1)*gridChart.cellSizeW;
				GeometryReal y1 = (GeometryReal)j*gridChart.cellSizeH , y2 = (GeometryReal)(j+1)*gridChart.cellSizeH;
				Point2D< GeometryReal > points[] = { Point2D< GeometryReal >(x1,y1) , Point2D< GeometryReal >(x2,y1) , Point2D< GeometryReal >(x2,y2) , Point2D< GeometryReal >(x1,y2) };
				Point2D< GeometryReal > normals[] = { vertices[1]-vertices[0] , vertices[2]-vertices[1] , vertices[0]-vertices[2] };
				for( int k=0 ; k<3 ; k++ ) normals[k] = Point2D< GeometryReal >( normals[k][1] , -normals[k][0] );
				GeometryReal offsets[] = { Point2D< GeometryReal >::Dot( normals[0] , vertices[0] ) , Point2D< GeometryReal >::Dot( normals[1] , vertices[1] ) , Point2D< GeometryReal >::Dot( normals[2] , vertices[2] ) };
				for( int k=0 ; k<4 ; k++ )
				{
					if( ( Point2D< GeometryReal >::Dot( points[k] , normals[0] ) - offsets[0] )<0 ) return false;
					if( ( Point2D< GeometryReal >::Dot( points[k] , normals[1] ) - offsets[1] )<0 ) return false;
					if( ( Point2D< GeometryReal >::Dot( points[k] , normals[2] ) - offsets[2] )<0 ) return false;
				}
				return true;
			};

		// Compute the square-root of the weights to make taking the weighted dot-product faster
		GeometryReal _integrator_sampleWeight[Samples];
		for( int s=0 ; s<Samples ; s++ ) _integrator_sampleWeight[s] = (GeometryReal)sqrt( TriangleIntegrator< Samples >::Weights[s] );
		SquareMatrix< GeometryReal , 2 > cell_to_texture_differential;
		cell_to_texture_differential(0,0) = (GeometryReal)gridChart.cellSizeW;
		cell_to_texture_differential(1,1) = (GeometryReal)gridChart.cellSizeH;
		cell_to_texture_differential(0,1) = cell_to_texture_differential(1,0) = 0;

		SquareMatrix< GeometryReal , 4 > interior_cell_mass , interior_cell_stiffnesses[2][2];
		SquareMatrix< GeometryReal , 4 > grad_edge_products[2][2];
		{
			std::vector< Point2D< GeometryReal > > polygon = { Point2D< GeometryReal >(0,0) , Point2D< GeometryReal >(1,0) , Point2D< GeometryReal >(1,1) , Point2D< GeometryReal >(0,1) };
			for( int p=2 ; p<polygon.size() ; p++ )
			{
				Point2D< GeometryReal > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };
				Point2D< GeometryReal > fragment_samples[Samples];
				for( int s=0 ; s<Samples ; s++ ) fragment_samples[s] = polygon[0] + dm[0] * (GeometryReal)TriangleIntegrator<Samples>::Positions[s][0] + dm[1] * (GeometryReal)TriangleIntegrator<Samples>::Positions[s][1];

				// Integrate scalar product and gradient field
				GeometryReal sampleValues[Samples][4];
				Point2D< GeometryReal > sampleGradients[Samples][4];
				Point2D< GeometryReal > __sampleGradients[Samples][4];

				for( int s=0 ; s<Samples ; s++ )
				{
					BilinearElementValuesAndGradients( fragment_samples[s] , sampleValues[s] , sampleGradients[s] );
					for( int k=0 ; k<4 ; k++ )
					{
						sampleValues   [s][k] *= _integrator_sampleWeight[s];
						__sampleGradients[s][k] = sampleGradients[s][k];
						sampleGradients[s][k] *= _integrator_sampleWeight[s];
					}
				}

				for( int k=0 ; k<4 ; k++ ) for( int l=0 ; l<4 ; l++ ) for( int s=0 ; s<Samples ; s++ )
				{
					interior_cell_mass(l,k) += sampleValues[s][k] * sampleValues[s][l] / 2;
					for( int m=0 ; m<2 ; m++ ) for( int n=0 ; n<2 ; n++ ) interior_cell_stiffnesses[m][n](l,k) += sampleGradients[s][l][m] * sampleGradients[s][k][n] / 2;
				}

				if( computeDivergence )
				{
					Point2D< GeometryReal > sampleVectorFields[Samples][4];
					for( int s=0 ; s<Samples ; s++ )
					{
						ReducedVectorFieldBasis( fragment_samples[s] , sampleVectorFields[s] );
						for( int k=0 ; k<4 ; k++ ) sampleVectorFields[s][k] *= _integrator_sampleWeight[s];
					}
					for (int k = 0; k < 4; k++)for (int l = 0; l < 4; l++) for (int s = 0; s<Samples; s++) {
						for (int m = 0; m < 2; m++) for (int n = 0; n < 2; n++) grad_edge_products[m][n](l, k) += sampleVectorFields[s][k][m] * sampleGradients[s][l][n] / 2;
					}
				}

			}
		}
		//#pragma omp parallel for 
		for( int t=0 ; t<atlasChart.triangles.size() ; t++ )
		{
			Point2D< GeometryReal > tPos[3];
			for( int i=0 ; i<3 ; i++ ) tPos[i] = atlasChart.vertices[ atlasChart.triangles[t][i] ] - gridChart.corner;

			SquareMatrix< GeometryReal , 2 > texture_metric = texture_metrics[t];
			SquareMatrix< GeometryReal , 2 > cell_metric = cell_to_texture_differential.transpose() * texture_metric * cell_to_texture_differential;
			SquareMatrix< GeometryReal , 2 > cell_metric_inverse = cell_metric.inverse();
			GeometryReal cell_area_scale_factor = (GeometryReal)sqrt( cell_metric.determinant() );

			//BBox
			int minCorner[2] , maxCorner[2];
			GetTriangleIntegerBBox( tPos , (GeometryReal)1./gridChart.cellSizeW , (GeometryReal)1./gridChart.cellSizeH , minCorner , maxCorner );

			std::vector< Point2D< GeometryReal > > parametricVertices(3);
			parametricVertices[0] = tPos[0] , parametricVertices[1] = tPos[1] , parametricVertices[2] = tPos[2];

			AtlasIndexedTriangle< GeometryReal > atlasTriangle;
			for( int k=0 ; k<3 ; k++ )
			{
				atlasTriangle.vertices[k] = tPos[k];
				atlasTriangle.atlasEdgeIndices[k] = atlasChart.atlasEdgeIndices[3 * t + k];
				atlasTriangle.atlasVertexIndices[k] = atlasChart.triangles[t][k];
				atlasTriangle.atlasVertexParentEdge[k] = -1;
			}

			for( int j=minCorner[1] ; j<maxCorner[1] ; j++ ) for( int i=minCorner[0] ; i<maxCorner[0] ; i++ )
			{
				auto TextureToCell = [&]( Point2D< GeometryReal > p ){ return Point2D< GeometryReal >( (GeometryReal)( p[0] / gridChart.cellSizeW ) - i , (GeometryReal)( p[1] / gridChart.cellSizeH ) - j ); };

				int localInteriorIndex = gridChart.localInteriorCellIndex(i,j) , localBoundaryIndex = gridChart.localBoundaryCellIndex(i,j);
				if( localInteriorIndex!=-1 && localBoundaryIndex!=-1 ) MK_THROW( "Cell simultaneosly interior and boundary" );

				// If the cell is entirely within the triangle...
				if( CellInTriangle( i , j , parametricVertices ) && localInteriorIndex!=-1 )
				{
					SquareMatrix< GeometryReal , 4 > polygonMass , polygonStiffness;

					polygonMass = interior_cell_mass * cell_area_scale_factor;
					for( int k=0 ; k<4 ; k++ ) for( int l=0 ; l<4 ; l++ )
					{
						GeometryReal v=0;
						for( int m=0 ; m<2 ; m++ ) for( int n=0 ; n<2 ; n++ ) v += cell_metric_inverse(m,n) * interior_cell_stiffnesses[m][n](k,l);
						polygonStiffness(k,l) = v * cell_area_scale_factor;
					}
					cellMass     [ localInteriorIndex ] += polygonMass;
					cellStiffness[ localInteriorIndex ] += polygonStiffness;

					if( computeDivergence )
					{
						SquareMatrix< GeometryReal , 4 > polygonDivergence;
						for( int k=0 ; k<4 ; k++ ) for( int l=0 ; l<4 ; l++ )
						{
							GeometryReal v=0;
							for (int m = 0; m<2; m++) for (int n = 0; n<2; n++) v += cell_metric_inverse(m, n) * grad_edge_products[m][n](l, k);
							polygonDivergence(l, k) = v * cell_area_scale_factor;
						}
						cellDivergence[localInteriorIndex] += polygonDivergence;
					}
				}
				else if( localInteriorIndex!=-1 )
				{
					// For interior cells, the cell and the element are the same thing
					auto TextureToElement = TextureToCell;

					SquareMatrix< GeometryReal , 2 > element_metric = cell_metric , element_metric_inverse = cell_metric_inverse;
					GeometryReal element_area_scale_factor = cell_area_scale_factor;
					CellClippedTriangle< GeometryReal > polygon = parametricVertices;

					// Clip the triangle to the cell
					if( ClipTriangleToPrimalCell( polygon , i , j , gridChart.cellSizeW , gridChart.cellSizeH ) )
					{
						// Transform the polygon vertices into the coordinate frame of the cell
						for( int ii=0 ; ii<polygon.size() ; ii++ ) polygon[ii] = TextureToElement( polygon[ii] );
						SquareMatrix< GeometryReal , 4 > polygonStiffness , polygonMass;
						SquareMatrix< GeometryReal , 4 > polygonDivergence;

						for( int p=2 ; p<polygon.size() ; p++ )
						{
							Point2D< GeometryReal > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };

							SquareMatrix< GeometryReal , 2 > fragment_to_element_differential;
							for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = dm[x][y];
							GeometryReal fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );

							Point2D< GeometryReal > fragment_samples[Samples];
							for( int s=0 ; s<Samples ; s++ )
							{
								fragment_samples[s] = polygon[0] + dm[0] * (GeometryReal)TriangleIntegrator<Samples>::Positions[s][0] + dm[1] * (GeometryReal)TriangleIntegrator<Samples>::Positions[s][1];
								if( !InUnitSquare( fragment_samples[s] ) ) MK_THROW( "Interior sample out of unit box! (" , fragment_samples[s][0] , " " , fragment_samples[s][1] , ")" );
							}

							// Integrate scalar product and gradient field
							// Make the code more efficient by:
							// -- pre-multiplying the values and gradients by the square-root of the quadrature weight
							// -- pre-multiplying the gradients by the inverse of the metric
							// so that the computation within the inner loop is faster.
							GeometryReal fragment_area = element_area_scale_factor * fragment_to_element_area_scale_factor / 2;
							GeometryReal sampleValues[Samples][4];
							Point2D< GeometryReal > sampleGradients[Samples][4] , _sampleGradients[Samples][4], __sampleGradients[Samples][4];
							for( int s=0 ; s<Samples ; s++ )
							{
								BilinearElementValuesAndGradients( fragment_samples[s] , sampleValues[s] , sampleGradients[s] );
								for( int k=0 ; k<4 ; k++ )
								{
									sampleValues[s][k] *= _integrator_sampleWeight[s];
									__sampleGradients[s][k] = sampleGradients[s][k];
									sampleGradients[s][k] *= _integrator_sampleWeight[s];
									_sampleGradients[s][k] = element_metric_inverse * sampleGradients[s][k];
								}
							}
							for( int k=0 ; k<4 ; k++ ) for( int l=0 ; l<4 ; l++ )
							{
								GeometryReal vIntegral=0 , gIntegral=0;
								for( int s=0 ; s<Samples ; s++ )
								{
									vIntegral += sampleValues[s][k] * sampleValues[s][l];
									gIntegral += Point2D< GeometryReal >::Dot( sampleGradients[s][l] , _sampleGradients[s][k] );
								}
								polygonMass(l,k) += vIntegral * fragment_area;
								polygonStiffness(l,k) += gIntegral * fragment_area;
							}
							if( computeDivergence )
							{
								Point2D< GeometryReal > sampleVectorFields[Samples][4];
								for (int s = 0; s < Samples; s++)
								{
									ReducedVectorFieldBasis(fragment_samples[s], sampleVectorFields[s]);
									for( int k=0 ; k<4 ; k++ ) sampleVectorFields[s][k] *= _integrator_sampleWeight[s];
								}
								for( int k=0 ; k<4 ; k++ ) for( int l=0 ; l<4 ; l++ )
								{
									GeometryReal dIntegral = 0;
									for( int s=0 ; s<Samples ; s++ ) dIntegral += Point2D< GeometryReal >::Dot( sampleVectorFields[s][k] , _sampleGradients[s][l] );
									polygonDivergence(l, k) += dIntegral * fragment_area;
								}
							}
						}
						cellMass[localInteriorIndex] += polygonMass;
						cellStiffness[localInteriorIndex] += polygonStiffness;
						if (computeDivergence) cellDivergence[localInteriorIndex] += polygonDivergence;
					}
				}
				else if( localBoundaryIndex!=-1 )
				{
					std::vector< BoundaryIndexedTriangle< GeometryReal > > cellBoundaryTriangles = gridChart.boundaryTriangles[localBoundaryIndex];

					// Iterate over all elements in the cell
					for( int bt=0 ; bt<cellBoundaryTriangles.size() ; bt++ )
					{
						BoundaryIndexedTriangle< GeometryReal > element = cellBoundaryTriangles[bt];
						std::vector< Point2D< GeometryReal > > element_vertices(3);
						for( int ii=0 ; ii<3 ; ii++ ) element_vertices[ii] = TextureToCell( element[ii] );
						int boundaryTriangleId = element.id;

						AtlasIndexedPolygon< GeometryReal > polygon;
						SetAtlasIndexedPolygonFromBoundaryTriangle( element , polygon );

						// Intersect the element with the atlas triangle
						int clippingResult = ClipPartiallyIndexedPolygonToIndexedTriangle( polygon , atlasTriangle );
						if( clippingResult>0 )
						{
							// Convert the polygon vertices from the texture frame to the cell frame
							for( int ii=0 ; ii<polygon.size() ; ii++ ) polygon[ii] = TextureToCell( polygon[ii] );

							SquareMatrix< GeometryReal , 2 > element_to_cell_differential , cell_to_element_differential;
							Point2D< GeometryReal > dm[] = { element_vertices[1]-element_vertices[0] , element_vertices[2]-element_vertices[0] };
							for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) element_to_cell_differential(x,y) = dm[x][y];
							cell_to_element_differential = element_to_cell_differential.inverse();
							auto CellToElement = [&]( Point2D< GeometryReal > p ){ return cell_to_element_differential * ( p - element_vertices[0] ); };
							// Convert the polygon vertices from the cell frame to the element frame
							for( int ii=0 ; ii<polygon.size() ; ii++ ) polygon[ii] = CellToElement( polygon[ii] );

							SquareMatrix< GeometryReal , 2 > element_metric = element_to_cell_differential.transpose() * cell_metric * element_to_cell_differential;
							SquareMatrix< GeometryReal , 2 > element_metric_inverse = element_metric.inverse();
							GeometryReal element_area_scale_factor = sqrt( element_metric.determinant() );
							SquareMatrix< GeometryReal , 6 > polygonStiffness , polygonMass;
							Matrix< GeometryReal , 6 , 15 > polygonDivergence;
							GeometryReal polygonArea = 0;

							for( int p=2 ; p<polygon.vertices.size() ; p++ )
							{
								Point2D< GeometryReal > d[] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };

								SquareMatrix< GeometryReal , 2 > fragment_to_element_differential;
								for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = d[x][y];
								GeometryReal fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );
								GeometryReal fragment_area = fragment_to_element_area_scale_factor * element_area_scale_factor / 2;
								if( fragment_area>0 )
								{
									polygonArea += fragment_area;

									Point2D< GeometryReal > fragment_samples[Samples];
									for( int s=0 ; s<Samples ; s++ )
									{
										fragment_samples[s] = polygon.vertices[0] + d[0] * (GeometryReal)TriangleIntegrator<Samples>::Positions[s][0] + d[1] * (GeometryReal)TriangleIntegrator<Samples>::Positions[s][1];
										if( !InUnitTriangle( fragment_samples[s] ) ) MK_THROW( "Boundary sample out of unit right triangle! (" , fragment_samples[s][0] , " " , fragment_samples[s][1] , ")" );
										else
										{
											fragment_samples[s][0] = std::max< GeometryReal >( fragment_samples[s][0] , 0 );
											fragment_samples[s][1] = std::max< GeometryReal >( fragment_samples[s][1] , 0 );
											GeometryReal excess = ( fragment_samples[s][0] + fragment_samples[s][1] ) - 1;
											if( excess>0 ) fragment_samples[s][0] -= excess/2 , fragment_samples[s][1] -= excess/2;
										}
									}
									GeometryReal sampleValues[Samples][6];
									Point2D< GeometryReal > sampleGradients[Samples][6] , _sampleGradients[Samples][6], __sampleGradients[Samples][6];
									for( int s=0 ; s<Samples ; s++ )
									{
										QuadraticElementValuesAndGradients( fragment_samples[s] , sampleValues[s] , sampleGradients[s] );
										for( int k=0 ; k<6 ; k++ )
										{
											sampleValues[s][k] *= _integrator_sampleWeight[s];
											__sampleGradients[s][k] = sampleGradients[s][k];
											sampleGradients[s][k] *= _integrator_sampleWeight[s];
											_sampleGradients[s][k] = element_metric_inverse * sampleGradients[s][k];										
										}
									}

									for( int k=0 ; k<6 ; k++ ) for( int l=0 ; l<6 ; l++ )
									{
										GeometryReal vIntegral=0 , gIntegral=0;
										for( int s=0 ; s<Samples ; s++ )
										{
											vIntegral += sampleValues[s][k] * sampleValues[s][l];
											gIntegral += Point2D< GeometryReal >::Dot( sampleGradients[s][k] , _sampleGradients[s][l] );
										}
										polygonMass( l , k ) += vIntegral * fragment_area;
										polygonStiffness( l , k ) += gIntegral * fragment_area;
									}

									if( computeDivergence )
									{
										Point2D< GeometryReal > sampleVectorFields[Samples][15];

										for( int k=1 ; k<6 ; k++ ) for( int l=0 ; l<k ; l++ )
										{
											int edgeId = (k*(k - 1)) / 2 + l;
											for (int s = 0; s < Samples; s++) {
												sampleVectorFields[s][edgeId] = sampleValues[s][k] * __sampleGradients[s][l] - sampleValues[s][l] * __sampleGradients[s][k];
											}
										}

										for( int k=0 ; k<15 ; k++ ) for( int l=0 ; l<6 ; l++ )
										{
											GeometryReal dIntegral = 0;
											for( int s=0 ; s<Samples ; s++ ) dIntegral += Point2D< GeometryReal >::Dot( sampleVectorFields[s][k] , _sampleGradients[s][l] );
											polygonDivergence(l,k) += dIntegral * fragment_area;
										}
									}
								}
								else
								{
									zeroAreaElementCount++;
									MK_WARN( "Zero area polygon at cell " , gridChart.cornerCoords[0] + i , " " , gridChart.cornerCoords[1] + j);
								}
							}

							GeometryReal integratedPolygonMass = 0;
							for( int k=0 ; k<6 ; k++ ) for( int l=0 ; l<6 ; l++ ) integratedPolygonMass += polygonMass(k,l);
							if( fabs( integratedPolygonMass - polygonArea )>PRECISION_ERROR ) MK_WARN( "Out of precision" );
							{
								for (int dk = 0; dk < 6; dk++)for (int dl = 0; dl < 6; dl++) {
									triangleElementStiffness[boundaryTriangleId](dk, dl) += (polygonStiffness(dk, dl) + polygonStiffness(dl, dk)) / 2;
									triangleElementMass[boundaryTriangleId](dk, dl) += (polygonMass(dk, dl) + polygonMass(dl, dk)) / 2;
								}
								if( computeDivergence ) triangleElementDivergence[boundaryTriangleId] += polygonDivergence;
							}
						}
					}
				}
			}
		}

		if( zeroAreaElementCount ) MK_WARN( "Element with zero area = " , zeroAreaElementCount );

		int offset_i[4] = { 0, 1 , 1 ,0 };
		int offset_j[4] = { 0, 0 , 1 ,1 };
		int cellOffset[4] = { 3, 2, 0, 1 };



		// Node indexing
		//		0 ---- 1
		//		|	   |
		//		|	   |
		//		3------2

		int cellEdgeCornerOffset[24] = { 13,  9,  0,  4,
			14, 10,  1,  5,
			15, 11,  2 , 6,
			16, 12,  3,  7,
			19, 18,  9, 13,
			17, 14,  5,  8 };

		// Interior texel neighbour edge indexing
		//		+--0---+---4--+ 
		//		|\    /|\    /|  
		//		| 2  / | 6  / |
		//		|  \/  |  \/  |
		//		1  /\  5  /\  8
		//		| 3  \ | 7  \ |
		//		|/    \|/    \|
		//		+--9---+--13--+
		//		|\    /|\    /|
		//		| 11 / | 15 / |
		//		|  \/  |  \/  |
		//		10 /\ 14  /\ 17
		//		| 12 \ | 16 \ |
		//		|/    \|/	 \|
		//		+--18--+--19--+


		int reducedCellEdgeCornerOffset[24] = {7, 5, 0, 2,
			8,  6,  1,  3,
			11, 10,  5,  7,
			9,  8,  3,  4
		};


		// Reduced interior texel neighbour edge indexing
		//		+-0--+--2-+ 
		//		|    |    |
		//		1    3    4
		//		|    |    |
		//		+-5--+-7--+
		//		|    |    |
		//		6    8    9
		//		|    |    |
		//      +-10-+-11-+

		auto NeighbourOffset = [&]( int k , int l ){ return ( offset_j[l] - offset_j[k] + 1 ) * 3 + ( offset_i[l] - offset_i[k] + 1 ); };
		for( int i=0 ; i<gridChart.interiorCellCorners.size() ; i++ )
		{
			const BilinearElementIndex & indicesGlobal = gridChart.interiorCellGlobalCorners[i];
			const BilinearElementIndex & indicesInterior = gridChart.interiorCellCorners[i];

			int localCellIndex = gridChart.interiorCellIndexToLocalCellIndex[i];
			int globalCellIndex = localCellIndex + gridChart.globalIndexCellOffset;

			int cellCoarseEdgeIndex[4];
			bool coarseEdgeIndexInitialized = false;

			Point< GeometryReal , 4 > prod[3];
			if( computeCellBasedStiffness )
			{
				Point3D< GeometryReal > values[4] = { inputSignal[ indicesGlobal[0] ] , inputSignal[ indicesGlobal[1] ] , inputSignal[ indicesGlobal[2] ] , inputSignal[ indicesGlobal[3] ] };
				for( int c=0 ; c<3 ; c++ )
				{
					Point< GeometryReal , 4 > v;
					v[0] = values[0][c];
					v[1] = values[1][c];
					v[2] = values[2][c];
					v[3] = values[3][c];

					prod[c] = cellStiffness[i] * v;
				}
			}

			for( int k=0 ; k<4 ; k++ )
			{
				int currentNode = indicesGlobal[k];
				int _currentBoundaryAndDeepIndex = boundaryAndDeepIndex[currentNode];
				if( _currentBoundaryAndDeepIndex<0 ) //Deep
				{
					int deepIndex = -_currentBoundaryAndDeepIndex - 1;
					for( int l=0 ; l<4 ; l++ )
					{
						deepMassCoefficients[ 10*deepIndex + NeighbourOffset(k,l) ] = (MatrixReal)( cellMass[i](k,l) + deepMassCoefficients[ 10*deepIndex + NeighbourOffset(k,l) ] );
						deepStiffnessCoefficients[ 10*deepIndex + NeighbourOffset(k,l) ] = (MatrixReal)( cellStiffness[i](k,l) + deepStiffnessCoefficients[ 10*deepIndex + NeighbourOffset(k,l) ] );
					}
					if( computeCellBasedStiffness )
					{
						//Add cell data
						texelToCellCoeffs[ 3*( 4*deepIndex + cellOffset[k] ) + 0 ] = (MatrixReal)prod[0][k];
						texelToCellCoeffs[ 3*( 4*deepIndex + cellOffset[k] ) + 1 ] = (MatrixReal)prod[1][k];
						texelToCellCoeffs[ 3*( 4*deepIndex + cellOffset[k] ) + 2 ] = (MatrixReal)prod[2][k];
					}
					if( computeDivergence ) for( int l=0 ; l<4 ; l++ ) deepDivergenceCoefficients[ 12*deepIndex + reducedCellEdgeCornerOffset[4*l+k] ] = (MatrixReal)( cellDivergence[i](k,l) + deepDivergenceCoefficients[ 12*deepIndex + reducedCellEdgeCornerOffset[4*l+k] ] );
				}
				else // Boundary
				{
					int boundaryIndex = _currentBoundaryAndDeepIndex - 1;
					for( int l=0 ; l<4 ; l++ )
					{
						int neighbourNode = indicesGlobal[l];
						int _neighbourBoundaryAndDeepIndex = boundaryAndDeepIndex[neighbourNode];
						if (_neighbourBoundaryAndDeepIndex < 0) {//Deep
							boundaryDeepMassTriplets.push_back( Eigen::Triplet< MatrixReal >( boundaryIndex , neighbourNode , (MatrixReal)cellMass[i](k,l) ) );
							boundaryDeepStiffnessTriplets.push_back(Eigen::Triplet< MatrixReal >( boundaryIndex , neighbourNode , (MatrixReal)cellStiffness[i](k,l) ) );
						}
						else // Boundary
						{
							boundaryBoundaryMassTriplets.push_back( Eigen::Triplet< MatrixReal >( fineBoundaryIndex[ indicesInterior[k] ] , fineBoundaryIndex[ indicesInterior[l] ] , (MatrixReal)cellMass[i](k,l) ) );
							boundaryBoundaryStiffnessTriplets.push_back( Eigen::Triplet< MatrixReal >( fineBoundaryIndex[ indicesInterior[k] ] , fineBoundaryIndex[ indicesInterior[l] ] , (MatrixReal)cellStiffness[i](k,l) ) );
						}
					}
					if( computeCellBasedStiffness )
					{
						// Add cell data
						int _fineBoundaryIndex = fineBoundaryIndex[indicesInterior[k]];
						Point3D< MatrixReal > p( (MatrixReal)prod[0][k] , (MatrixReal)prod[1][k] , (MatrixReal)prod[2][k] );
						boundaryCellStiffnessTriplets.push_back( Eigen::Triplet< Point3D< MatrixReal > >( _fineBoundaryIndex , globalCellIndex , p ) );
					}
					if( computeDivergence )
					{
						if( !coarseEdgeIndexInitialized )
						{
							for( int l=0 ; l<4 ; l++ )
							{
								int edgeSourceCoarseIndex = indicesGlobal[reducedCellCornerPairs[2 * l]];
								int edgeTargetCoarseIndex = indicesGlobal[reducedCellCornerPairs[2 * l + 1]];
								unsigned long long coarseEdgeKey = SetMeshEdgeKey(edgeSourceCoarseIndex, edgeTargetCoarseIndex);
								if( coarseEdgeIndex.find(coarseEdgeKey)==coarseEdgeIndex.end() ) MK_THROW( "Fine edge not found" );
								cellCoarseEdgeIndex[l] = coarseEdgeIndex[coarseEdgeKey];
							}
							coarseEdgeIndexInitialized = true;
						}
						for( int l=0 ; l<4 ; l++ ) boundaryDeepDivergenceTriplets.push_back( Eigen::Triplet< MatrixReal >( currentNode , cellCoarseEdgeIndex[l] , (MatrixReal)cellDivergence[i](k,l) ) );
					}
				}

			}
		}

		for( int c=0 ; c<gridChart.boundaryTriangles.size() ; c++ )
		{
			int localCellIndex = gridChart.boundaryCellIndexToLocalCellIndex[c];
			int globalCellIndex = localCellIndex + gridChart.globalIndexCellOffset;

			for (int b = 0; b < gridChart.boundaryTriangles[c].size(); b++) {
				int i = gridChart.boundaryTriangles[c][b].id;
				const QuadraticElementIndex & indices = gridChart.boundaryTriangles[c][b].indices;
				QuadraticElementIndex fineTriangleElementIndices;
				for (int k = 0; k < 6; k++) fineTriangleElementIndices[k] = fineBoundaryIndex[indices[k]];

				// Add cell data

				if( computeCellBasedStiffness )
				{
					Point< GeometryReal , 6 > prod[3];
					Point3D < GeometryReal > values[6] = { boundarySignal[fineTriangleElementIndices[0]],boundarySignal[fineTriangleElementIndices[1]],boundarySignal[fineTriangleElementIndices[2]],boundarySignal[fineTriangleElementIndices[3]],boundarySignal[fineTriangleElementIndices[4]],boundarySignal[fineTriangleElementIndices[5]] };
					for( int cc=0 ; cc<3 ; cc++ )
					{
						Point< GeometryReal , 6 > v;
						v[0] = values[0][cc];
						v[1] = values[1][cc];
						v[2] = values[2][cc];
						v[3] = values[3][cc];
						v[4] = values[4][cc];
						v[5] = values[5][cc];
						prod[cc] = triangleElementStiffness[i] * v;
					}

					for (int k = 0; k < 6; k++)
					{
						Point3D< MatrixReal > p( (MatrixReal)prod[0][k] , (MatrixReal)prod[1][k] , (MatrixReal)prod[2][k] );
						boundaryCellStiffnessTriplets.push_back( Eigen::Triplet< Point3D< MatrixReal > >( fineTriangleElementIndices[k] , globalCellIndex , p ) );
					}
				}
				for( int k=0 ; k<6 ; k++ ) for( int l=0 ; l<6 ; l++ )
				{
					boundaryBoundaryMassTriplets.push_back( Eigen::Triplet< MatrixReal >( fineTriangleElementIndices[k] , fineTriangleElementIndices[l] , (MatrixReal)triangleElementMass[i](l,k) ) );
					boundaryBoundaryStiffnessTriplets.push_back(Eigen::Triplet< MatrixReal >( fineTriangleElementIndices[k] , fineTriangleElementIndices[l] , (MatrixReal)triangleElementStiffness[i](l,k) ) );
				}

				if( computeDivergence )
				{
					for (int k = 1; k < 6; k++)for (int l = 0; l < k; l++) {
						int edgeId = (k*(k - 1)) / 2 + l;
						int edgeSourceFineIndex = fineTriangleElementIndices[k];
						int edgeTargetFineIndex = fineTriangleElementIndices[l];
						MatrixReal edgeSign = (MatrixReal)1.;
						if( edgeSourceFineIndex>edgeTargetFineIndex )
						{
							std::swap(edgeSourceFineIndex, edgeTargetFineIndex);
							edgeSign = (MatrixReal)-1.;
						}
						unsigned long long fineEdgeKey = SetMeshEdgeKey(edgeSourceFineIndex, edgeTargetFineIndex);
						int _fineEdgeIndex = -1;
						if( fineBoundaryEdgeIndex.find( fineEdgeKey )!=fineBoundaryEdgeIndex.end() ) _fineEdgeIndex = fineBoundaryEdgeIndex[fineEdgeKey];
						else MK_THROW( "Fine edge not found" );
						for( int n=0 ; n<6 ; n++ )
						{
							int fineNodeIndex = fineTriangleElementIndices[n];
							boundaryBoundaryDivergenceTriplets.push_back( Eigen::Triplet< MatrixReal >( fineNodeIndex , _fineEdgeIndex , (MatrixReal)triangleElementDivergence[i](n,edgeId) * edgeSign ) );
						}
					}
				}
			}
		}
	}

	template< unsigned int Samples , typename GeometryReal , typename MatrixReal >
	void InitializeMassAndStiffness
	(
		const std::vector< std::vector< SquareMatrix< GeometryReal , 2 > > > &parameterMetric ,
		const std::vector< AtlasChart< GeometryReal > > &atlasCharts ,
		const GridAtlas< GeometryReal , MatrixReal > &gridAtlas ,
		const std::vector< int > &fineBoundaryIndex ,
		int numFineBoundaryNodes ,
		std::vector< MatrixReal >& deepMassCoefficients ,
		std::vector< MatrixReal >& deepStiffnessCoefficients ,
		SparseMatrix< MatrixReal , int >& boundaryBoundaryMassMatrix ,
		SparseMatrix< MatrixReal , int >& boundaryBoundaryStiffnessMatrix ,
		SparseMatrix< MatrixReal , int >& boundaryDeepMassMatrix ,
		SparseMatrix< MatrixReal , int >& boundaryDeepStiffnessMatrix ,
		bool computeCellBasedStiffness ,
		const std::vector< Point3D< MatrixReal > >& inputSignal ,
		const std::vector< Point3D< MatrixReal > >& boundarySignal ,
		std::vector< MatrixReal >& texelToCellCoeffs ,
		SparseMatrix< MatrixReal , int > boundaryCellBasedStiffnessRHSMatrix[3] ,
		bool computeDivergence ,
		std::unordered_map< unsigned long long , int > & fineBoundaryEdgeIndex,
		std::unordered_map< unsigned long long , int > & coarseEdgeIndex,
		std::vector< Eigen::Triplet< MatrixReal > > & boundaryDeepDivergenceTriplets,
		std::vector< Eigen::Triplet< MatrixReal > > & boundaryBoundaryDivergenceTriplets,
		std::vector< MatrixReal >& deepDivergenceCoefficients
	)
	{
		auto MergeTriplets = [] ( const std::vector< std::vector< Eigen::Triplet< MatrixReal > > > &inTriplets , std::vector< Eigen::Triplet< MatrixReal > > &outTriplets )
			{
				size_t count = 0;
				for( int i=0 ; i<inTriplets.size() ; i++ ) count += inTriplets[i].size();
				outTriplets.reserve( count );
				for( int i=0 ; i<inTriplets.size() ; i++ ) for( int j=0 ; j<inTriplets[i].size() ; j++ ) outTriplets.push_back( inTriplets[i][j] );
			};

		const std::vector< GridChart< GeometryReal > > &gridCharts = gridAtlas.gridCharts;
		const std::vector< int > &boundaryAndDeepIndex = gridAtlas.boundaryAndDeepIndex;

		if( computeCellBasedStiffness ) texelToCellCoeffs.resize( 3*4*gridAtlas.numDeepTexels );

		std::vector< Eigen::Triplet< MatrixReal > > boundaryBoundaryMassTriplets;
		std::vector< Eigen::Triplet< MatrixReal > > boundaryBoundaryStiffnessTriplets;
		std::vector< Eigen::Triplet< MatrixReal > > boundaryDeepMassTriplets;
		std::vector< Eigen::Triplet< MatrixReal > > boundaryDeepStiffnessTriplets;

		std::vector< std::vector< Eigen::Triplet< MatrixReal > > > _boundaryBoundaryMassTriplets              ( ThreadPool::NumThreads() );
		std::vector< std::vector< Eigen::Triplet< MatrixReal > > > _boundaryBoundaryStiffnessTriplets         ( ThreadPool::NumThreads() );
		std::vector< std::vector< Eigen::Triplet< MatrixReal > > > _boundaryDeepMassTriplets                  ( ThreadPool::NumThreads() );
		std::vector< std::vector< Eigen::Triplet< MatrixReal > > > _boundaryDeepStiffnessTriplets             ( ThreadPool::NumThreads() );
		std::vector< std::vector< Eigen::Triplet< MatrixReal > > > _boundaryDeepDivergenceTriplets            ( ThreadPool::NumThreads() );
		std::vector< std::vector< Eigen::Triplet< MatrixReal > > > _boundaryBoundaryDivergenceTriplets        ( ThreadPool::NumThreads() );
		std::vector< std::vector< Eigen::Triplet< Point3D< MatrixReal > > > > _boundaryCellStiffnessTriplets  ( ThreadPool::NumThreads() );

		ThreadPool::ParallelFor
		(
			0 , gridCharts.size() ,
			[&]( unsigned int thread , size_t i )
			{
				InitializeChartMassAndStiffness< Samples >
					(
						parameterMetric[i] , atlasCharts[i] , gridCharts[i] , boundaryAndDeepIndex , fineBoundaryIndex , deepMassCoefficients , deepStiffnessCoefficients , 
						_boundaryBoundaryMassTriplets[thread] , _boundaryBoundaryStiffnessTriplets[thread] ,
						_boundaryDeepMassTriplets[thread] , _boundaryDeepStiffnessTriplets[thread] ,
						computeCellBasedStiffness , inputSignal , boundarySignal , texelToCellCoeffs ,
						_boundaryCellStiffnessTriplets[thread] ,
						computeDivergence , fineBoundaryEdgeIndex , coarseEdgeIndex ,
						_boundaryDeepDivergenceTriplets[thread] , _boundaryBoundaryDivergenceTriplets[thread] ,
						deepDivergenceCoefficients
					);			}
		);

		MergeTriplets( _boundaryBoundaryMassTriplets , boundaryBoundaryMassTriplets );
		MergeTriplets( _boundaryBoundaryStiffnessTriplets , boundaryBoundaryStiffnessTriplets );
		MergeTriplets( _boundaryDeepMassTriplets , boundaryDeepMassTriplets );
		MergeTriplets( _boundaryDeepStiffnessTriplets , boundaryDeepStiffnessTriplets );
		MergeTriplets( _boundaryDeepDivergenceTriplets , boundaryDeepDivergenceTriplets );
		MergeTriplets( _boundaryBoundaryDivergenceTriplets , boundaryBoundaryDivergenceTriplets );

		int numTexels = gridAtlas.numTexels;
		int numBoundaryTexels = gridAtlas.numBoundaryTexels;
		boundaryBoundaryMassMatrix      = SetSparseMatrix( boundaryBoundaryMassTriplets , numFineBoundaryNodes , numFineBoundaryNodes , true );
		boundaryBoundaryStiffnessMatrix = SetSparseMatrix( boundaryBoundaryStiffnessTriplets , numFineBoundaryNodes , numFineBoundaryNodes , true );
		boundaryDeepMassMatrix          = SetSparseMatrix( boundaryDeepMassTriplets , numBoundaryTexels , numTexels , false );
		boundaryDeepStiffnessMatrix     = SetSparseMatrix( boundaryDeepStiffnessTriplets , numBoundaryTexels , numTexels , false );

		if( computeCellBasedStiffness )
		{
			size_t count = 0;
			for( int i=0 ; i<_boundaryCellStiffnessTriplets.size() ; i++ ) count += _boundaryCellStiffnessTriplets[i].size();
			std::vector< Eigen::Triplet< MatrixReal > > boundaryCellStiffnessTriplets( count );

			for( int c=0 ; c<3 ; c++ )
			{
				for( int i=0 , idx=0 ; i<_boundaryCellStiffnessTriplets.size() ; i++ ) for( int j=0 ; j<_boundaryCellStiffnessTriplets[i].size() ; j++ , idx++ )
					boundaryCellStiffnessTriplets[idx] = Eigen::Triplet< MatrixReal >( _boundaryCellStiffnessTriplets[i][j].row() , _boundaryCellStiffnessTriplets[i][j].col() , _boundaryCellStiffnessTriplets[i][j].value()[c] );
				boundaryCellBasedStiffnessRHSMatrix[c] = SetSparseMatrix( boundaryCellStiffnessTriplets , numFineBoundaryNodes , gridAtlas.numCells , false );
			}
		}
	}

	template< unsigned int Samples , typename GeometryReal , typename MatrixReal >
	void InitializeMassAndStiffness
	(
		SystemCoefficients< MatrixReal > &mass ,
		SystemCoefficients< MatrixReal > &stiffness ,
		const HierarchicalSystem< GeometryReal , MatrixReal > &hierarchy ,
		const std::vector< std::vector< SquareMatrix< GeometryReal , 2 > > > &parameterMetric ,
		const std::vector< AtlasChart< GeometryReal > > &atlasCharts ,
		const BoundaryProlongationData< MatrixReal > &boundaryProlongation ,
		bool computeCellBasedStiffness ,
		const std::vector< Point3D< MatrixReal > > &inputSignal ,
		std::vector< MatrixReal > &texelToCellCoeffs ,
		SparseMatrix< MatrixReal , int > boundaryCellBasedStiffnessRHSMatrix[3] ,
		bool computeDivergence ,
		std::unordered_map< unsigned long long , int > & edgeIndex ,
		SparseMatrix< MatrixReal , int > & boundaryDivergenceMatrix ,
		std::vector< MatrixReal > & deepDivergenceCoefficients
	)
	{
		//(2) Initialize mass and stiffness
		mass.deepCoefficients.resize( 10*hierarchy.gridAtlases[0].numDeepTexels , 0 );
		stiffness.deepCoefficients.resize( 10*hierarchy.gridAtlases[0].numDeepTexels , 0 );
		if( computeDivergence ) deepDivergenceCoefficients.resize( 20*hierarchy.gridAtlases[0].numDeepTexels , 0 );

		SparseMatrix< MatrixReal , int > fineBoundaryBoundaryMassMatrix;
		SparseMatrix< MatrixReal , int > fineBoundaryBoundaryStiffnessMatrix;

		std::vector< Eigen::Triplet< MatrixReal > > boundaryDivergenceTriplets;
		std::vector< Eigen::Triplet< MatrixReal > > boundaryBoundaryDivergenceTriplets;
		std::unordered_map< unsigned long long , int > fineBoundaryEdgeIndex;

		if( computeDivergence )
		{
			InitializeFineBoundaryEdgeIndexing( boundaryProlongation.fineBoundaryIndex , fineBoundaryEdgeIndex , hierarchy.gridAtlases[0].gridCharts );
		}

		SparseMatrix< MatrixReal , int > fineBoundaryCellStiffnessRHSMatrix[3];
		std::vector< Point3D< MatrixReal > > fineBoundarySignal;

		if( computeCellBasedStiffness )
		{
			const std::vector< int > & boundaryGlobalIndex = hierarchy.gridAtlases[0].boundaryGlobalIndex;
			int numBoundaryTexels = (int)boundaryGlobalIndex.size();
			int numFineBoundaryNodes = boundaryProlongation.numFineBoundarNodes;
			std::vector< Point3D< MatrixReal > > coarseBoundarySignal;
			coarseBoundarySignal.resize( numBoundaryTexels );
			for( int i=0 ; i<numBoundaryTexels ; i++ ) coarseBoundarySignal[i] = inputSignal[ boundaryGlobalIndex[i] ];
			fineBoundarySignal.resize( numFineBoundaryNodes );
			boundaryProlongation.coarseBoundaryFineBoundaryProlongation.Multiply( &coarseBoundarySignal[0] , &fineBoundarySignal[0] );
		}

		InitializeMassAndStiffness< Samples >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0] , boundaryProlongation.fineBoundaryIndex , boundaryProlongation.numFineBoundarNodes , mass.deepCoefficients , stiffness.deepCoefficients , fineBoundaryBoundaryMassMatrix , fineBoundaryBoundaryStiffnessMatrix , mass.boundaryDeepMatrix , stiffness.boundaryDeepMatrix , computeCellBasedStiffness , inputSignal , fineBoundarySignal , texelToCellCoeffs , fineBoundaryCellStiffnessRHSMatrix , computeDivergence , fineBoundaryEdgeIndex , edgeIndex , boundaryDivergenceTriplets , boundaryBoundaryDivergenceTriplets , deepDivergenceCoefficients );

		{
			SparseMatrix< MatrixReal , int > temp = fineBoundaryBoundaryMassMatrix * boundaryProlongation.coarseBoundaryFineBoundaryProlongation;
			mass.boundaryBoundaryMatrix = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction * temp;
		}
		{
			SparseMatrix< MatrixReal , int > temp = fineBoundaryBoundaryStiffnessMatrix * boundaryProlongation.coarseBoundaryFineBoundaryProlongation;
			stiffness.boundaryBoundaryMatrix = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction * temp;
		}

		{
			std::vector< MatrixReal > in ( mass.boundaryBoundaryMatrix.Rows() , (MatrixReal)1. );
			std::vector< MatrixReal > out( mass.boundaryBoundaryMatrix.Rows() , (MatrixReal)0. );
			mass.boundaryBoundaryMatrix.Multiply( GetPointer(in) , GetPointer(out) );
			for( int i=0 ; i<out.size() ; i++ ) if( out[i]==0 ) MK_WARN( "Zero row at boundary index " , i , ". Try running with jittering." );
		}

		if( computeCellBasedStiffness ) for( int c=0 ; c<3 ; c++ ) boundaryCellBasedStiffnessRHSMatrix[c] = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction * fineBoundaryCellStiffnessRHSMatrix[c];

		if( computeDivergence )
		{
			SparseMatrix< MatrixReal , int > fineBoundaryBoundaryDivergenceMatrix = SetSparseMatrix( boundaryBoundaryDivergenceTriplets , boundaryProlongation.numFineBoundarNodes , (int)fineBoundaryEdgeIndex.size() , false );

			std::unordered_map< unsigned long long , int > boundaryCoarseEdgeIndex;
			std::vector< int > boundaryCoarseEdgeToGlobalEdge;

			InitializeBoundaryEdgeIndexing( mass.boundaryBoundaryMatrix , hierarchy.gridAtlases[0].boundaryGlobalIndex , edgeIndex , boundaryCoarseEdgeToGlobalEdge , boundaryCoarseEdgeIndex );

			SparseMatrix< MatrixReal , int > boundaryCoarseToFineBoundaryOneFormProlongation;
			InitializeBoundaryCoarseToFineBoundaryOneFormProlongation< MatrixReal >( boundaryProlongation.coarseBoundaryFineBoundaryProlongation , boundaryCoarseEdgeIndex , fineBoundaryEdgeIndex , boundaryCoarseToFineBoundaryOneFormProlongation );

			SparseMatrix< MatrixReal , int > temp = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction * fineBoundaryBoundaryDivergenceMatrix;
			SparseMatrix< MatrixReal , int > boundaryBoundaryDivergenceMatrix = temp *  boundaryCoarseToFineBoundaryOneFormProlongation;
			const std::vector< int > &boundaryGlobalIndex = hierarchy.gridAtlases[0].boundaryGlobalIndex;
			for( int i=0 ; i<boundaryBoundaryDivergenceMatrix.Rows() ; i++ )
			{
				int nodeIndex = boundaryGlobalIndex[i];
				for( int j=0 ; j<boundaryBoundaryDivergenceMatrix.RowSize(i) ; j++ )
					boundaryDivergenceTriplets.push_back( Eigen::Triplet< MatrixReal >( nodeIndex , boundaryCoarseEdgeToGlobalEdge[ boundaryBoundaryDivergenceMatrix[i][j].N ] , boundaryBoundaryDivergenceMatrix[i][j].Value ) );
			}
			boundaryDivergenceMatrix = SetSparseMatrix( boundaryDivergenceTriplets , hierarchy.gridAtlases[0].numTexels , (int)edgeIndex.size() , false );
		}
	}
	template< unsigned int Samples , typename GeometryReal , typename MatrixReal >
	void InitializeMassAndStiffness
	(
		SystemCoefficients< MatrixReal > &mass ,
		SystemCoefficients< MatrixReal > &stiffness ,
		const HierarchicalSystem< GeometryReal , MatrixReal >& hierarchy ,
		const std::vector< std::vector< SquareMatrix< GeometryReal , 2 > > > &parameterMetric ,
		const std::vector< AtlasChart< GeometryReal > > &atlasCharts ,
		const BoundaryProlongationData< MatrixReal > &boundaryProlongation ,
		bool computeCellBasedStiffness ,
		const std::vector< Point3D< MatrixReal > > &inputSignal ,
		std::vector< MatrixReal >& texelToCellCoeffs ,
		SparseMatrix< MatrixReal , int > boundaryCellBasedStiffnessRHSMatrix[3]
	)
	{
		std::unordered_map< unsigned long long , int > edgeIndex;
		SparseMatrix< MatrixReal , int > boundaryDivergenceMatrix;
		std::vector< MatrixReal > deepDivergenceCoefficients;
		InitializeMassAndStiffness< Samples >( mass , stiffness , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , computeCellBasedStiffness , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix , false , edgeIndex , boundaryDivergenceMatrix , deepDivergenceCoefficients );
	}
}
