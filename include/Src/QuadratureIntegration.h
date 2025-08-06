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

#include <Misha/Atomic.h>
#include <Misha/Atomic.Geometry.h>
#include <Misha/Array.h>
#include "PolygonClipping.h"

namespace MishaK
{
	namespace TSP
	{
		class InteriorCellLine
		{
		public:
			AtlasTexelIndex prevLineIndex , nextLineIndex;
			unsigned int length;
		};

		template< typename GeometryReal >
		void InitializeGridChartInteriorCellLines
		(
			const GridChart< GeometryReal > &gridChart ,
			std::vector< InteriorCellLine > &interiorCellLines ,
			ExplicitIndexVector< AtlasInteriorCellIndex , std::pair< unsigned int , unsigned int > > &interiorCellLineIndex
		)
		{
			const RegularGrid< 2 , CellType >& cellType = gridChart.cellType;
			unsigned int width = cellType.res(0) , height = cellType.res(1);

			const RegularGrid< 2 , TexelIndex > & texelIndices = gridChart.texelIndices;

			ChartInteriorCellIndex chartInteriorCellIndex(0);

			for( unsigned int j=0 ; j<height ; j++ )
			{
				unsigned int offset = 0;
				bool previousIsInterior = false;
				unsigned int rasterStart = -1;
				while( offset<width )
				{
					bool currentIsInterior = cellType(offset,j)==CellType::Interior;
#ifdef SANITY_CHECK
					if( ( offset==0 || offset==width-1 ) && currentIsInterior ) MK_THROW( "Unexpected interior cell" );
#endif // SANITY_CHECK
					if(  currentIsInterior && !previousIsInterior ) rasterStart = offset; //Start raster line
					if( !currentIsInterior &&  previousIsInterior ) //Terminate raster line
					{ 
						InteriorCellLine newLine;
						newLine.prevLineIndex = texelIndices( rasterStart , j+0 ).combined;
						newLine.nextLineIndex = texelIndices( rasterStart , j+1 ).combined;
						newLine.length = offset - rasterStart;

#ifdef SANITY_CHECK
						if( newLine.prevLineIndex==AtlasTexelIndex(-1) || newLine.nextLineIndex==AtlasTexelIndex(-1) ) MK_THROW( "Invalid indexing" );
#endif // SANITY_CHECK
						int currentLine = (int)interiorCellLines.size();

						for( unsigned int k=0 ; k<offset-rasterStart ; k++ )
						{
#ifdef SANITY_CHECK
							if( gridChart.interiorCellCoveredTexelBilinearElementIndices[chartInteriorCellIndex][0]!=texelIndices( rasterStart+k , j ).covered ) MK_THROW( "Unexpected corner ID" );
#endif // SANITY_CHECK

							interiorCellLineIndex.push_back( std::pair< int , int >( currentLine , k ) );
							chartInteriorCellIndex++;
						}

						interiorCellLines.push_back( newLine );
					}
					previousIsInterior = currentIsInterior;
					offset++;
				}
			}
		}

		template< typename GeometryReal >
		void InitializeGridAtlasInteriorCellLines
		(
			const ExplicitIndexVector< ChartIndex , GridChart< GeometryReal > >& gridCharts ,
			std::vector< InteriorCellLine >& interiorCellLines ,
			ExplicitIndexVector< AtlasInteriorCellIndex , std::pair< unsigned int , unsigned int > >& interiorCellLineIndex
		)
		{
			for( unsigned int i=0 ; i<gridCharts.size() ; i++ ) InitializeGridChartInteriorCellLines( gridCharts[ ChartIndex(i) ] , interiorCellLines , interiorCellLineIndex );
		}

		//////////////////////////////////////
		// Pre-compute the integration info //
		//////////////////////////////////////

		// In what follows we have the following mappings:
		// fragment -> element -> texture:
		//		fragment: a triangle in the decomposition of the clipping of the atlas triangle to the elements
		//		element:  either a square or a triangle (depending on whether the cell is interior or not)
		//		texture: the texture space grid
		template< unsigned int Samples , typename Real >
		void SetInteriorCellDuals
		(
			typename BilinearElementScalarSample< Real >::SampleData& sampleData ,
			Point2D< Real > pos ,
			Real fragment_to_element_quadrature_weight
		)
		{
			// Compute the integrated Scalar of the bilinear basis function in the frame of the cell, weighted by the area of the unit triangle
			for( int k=0 ; k<4 ; k++ ) sampleData.dualValues[k] += (Real)BilinearElementValue( k , pos ) * fragment_to_element_quadrature_weight;
		}

		template< unsigned int Samples , typename Real >
		void SetInteriorCellDuals
		(
			typename BilinearElementGradientSample< Real >::SampleData& sampleData ,
			Point2D< Real > pos ,
			Real fragment_to_element_quadrature_weight
		)
		{
			// Compute the integrated gradient of the bilinear basis function in the frame of the cell, weighted by the area of the unit triangle
			for( int k=0 ; k<4 ; k++ ) sampleData.dualGradients[k] += Point2D< Real >( BilinearElementGradient( k , pos ) * fragment_to_element_quadrature_weight );
		}

		template< unsigned int Samples , typename Real >
		void SetCellInTriangleDuals
		(
			typename BilinearElementScalarSample< Real >::SampleData& sampleData ,
			const BilinearElementScalarSample< Real > & ,
			const typename BilinearElementScalarSample< Real >::SampleData& interiorSampleData ,
			Real element_area
		)
		{
			for( unsigned int k=0 ; k<4 ; k++ ) sampleData.dualValues[k] = interiorSampleData.dualValues[k] * element_area;
		}

		template< unsigned int Samples , typename Real >
		void SetCellInTriangleDuals
		(
			typename BilinearElementGradientSample< Real >::SampleData& sampleData ,
			const BilinearElementGradientSample< Real >& sample ,
			const typename BilinearElementGradientSample< Real >::SampleData& interiorSampleData ,
			Real element_area
		)
		{
			for( unsigned int k=0 ; k<4 ; k++ ) sampleData.dualGradients[k] = sample.tensor * interiorSampleData.dualGradients[k] * element_area;
		}

		template< unsigned int Samples , typename Real >
		void SetInteriorDuals
		(
			typename BilinearElementScalarSample< Real >::SampleData& sampleData ,
			const BilinearElementScalarSample< Real > & ,
			Point2D< Real > pos , 
			Real fragment_quadrature_weight
		)
		{
			// Compute the integrated gradients of each bilinear basis functions in the frame of the cell, weighted by the area of the unit triangle
			// Dualize the samples so that integration is simple a dot-product
			for( unsigned int k=0 ; k<4 ; k++ ) sampleData.dualValues[k] += (Real)BilinearElementValue( k , pos ) * fragment_quadrature_weight;
		}

		template< unsigned int Samples , typename Real >
		void SetInteriorDuals
		(
			typename BilinearElementGradientSample< Real >::SampleData& sampleData ,
			const BilinearElementGradientSample< Real >& sample ,
			Point2D< Real > pos , 
			Real fragment_quadrature_weight
		)
		{
			// Compute the integrated gradients of each bilinear basis functions in the frame of the cell, weighted by the area of the unit triangle
			// Dualize the samples so that integration is simple a dot-product
			for( unsigned int k=0 ; k<4 ; k++ ) sampleData.dualGradients[k] += sample.tensor * Point2D< Real >( BilinearElementGradient( k , pos ) * fragment_quadrature_weight );
		}

		template< unsigned int Samples , typename Real >
		void SetBoundaryDuals
		(
			typename QuadraticElementScalarSample< Real >::SampleData& sampleData ,
			const QuadraticElementScalarSample< Real > & ,
			Point2D< Real > pos , 
			Real fragment_quadrature_weight
		)
		{
			// Compute the integrated gradients of each quadratic basis functions in the frame of the cell, weighted by the area of the unit triangle
			// Dualize the samples so that integration is simple a dot-product
			for( unsigned int k=0 ; k<6 ; k++ ) sampleData.dualValues[k] += (Real)QuadraticElement::Value( k , pos ) * fragment_quadrature_weight;
		}

		template< unsigned int Samples , typename Real >
		void SetBoundaryDuals
		(
			typename QuadraticElementGradientSample< Real >::SampleData& sampleData ,
			const QuadraticElementGradientSample< Real >& sample ,
			Point2D< Real > pos , 
			Real fragment_quadrature_weight
		)
		{
			// Compute the integrated gradients of each quadratic basis functions in the frame of the cell, weighted by the area of the unit triangle
			// Dualize the samples so that integration is simple a dot-product
			for( unsigned int k=0 ; k<6 ; k++ ) sampleData.dualGradients[k] += sample.tensor * Point2D< Real >( QuadraticElement::Differential( k , pos ) * fragment_quadrature_weight );
		}

		template< unsigned int Samples , typename GeometryReal , typename ElementSamples >
		void InitializeIntegration
		(
			const ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > &texture_metrics ,
			const AtlasChart< GeometryReal > &atlasChart ,
			const GridChart< GeometryReal > &gridChart ,
			const ExplicitIndexVector< AtlasInteriorCellIndex , std::pair< unsigned int , unsigned int > >& interiorCellLineIndex ,
			const std::vector< AtlasInteriorOrBoundaryNodeIndex > &fineBoundaryIndex ,
			ElementSamples &elementSamples ,
			bool fastIntegration
		)
		{
			// Rasterize
			const GeometryReal precision_error = PRECISION_EPSILON;

			auto InUnitSquare =   [&]( Point2D< GeometryReal > p ){ return !( p[0]<0-precision_error || p[1]<0-precision_error || p[0]>1+precision_error || p[1]>1+precision_error ); };
			auto InUnitTriangle = [&]( Point2D< GeometryReal > p ){ return !( p[0]<0-precision_error || p[1]<0-precision_error || ( p[0]+p[1] )>1+precision_error ); };
			auto CellInTriangle = [&]( int i , int j , const std::vector< Point2D< GeometryReal > >& vertices )
				{
					Point2D< GeometryReal > points[] = { gridChart.nodePosition(i,j) , gridChart.nodePosition(i+1,j) , gridChart.nodePosition(i+1,j+1) , gridChart.nodePosition(i,j+1) };
					EdgeEquation< GeometryReal > eq[3];
					Point2D< GeometryReal > c = ( vertices[0] + vertices[1] + vertices[2] ) / 3;
					for( unsigned int i=0 ; i<3 ; i++ )
					{
						SimplexIndex< 1 > eIndex = OutgoingEdgeIndex( i );
						eq[i] = EdgeEquation< GeometryReal >( vertices[ eIndex[0] ] , vertices[ eIndex[1] ] );
						eq[i].makePositive(c);
					}
					for( int k=0 ; k<4 ; k++ ) for( unsigned int i=0 ; i<3 ; i++ ) if( eq[i]( points[k] )<0 ) return false;
					return true;
				};

			SquareMatrix< GeometryReal , 2 > cell_to_texture_differential;

			cell_to_texture_differential(0,0) = gridChart.cellSizeW;
			cell_to_texture_differential(1,1) = gridChart.cellSizeH;
			cell_to_texture_differential(0,1) = cell_to_texture_differential(1,0) = 0.;

			typename ElementSamples::Bilinear::SampleData interior_cell_samples[2*Samples] , interior_cell_sample;
			{
				std::vector< Point2D< GeometryReal > > polygon = { Point2D< GeometryReal >(0,0) , Point2D< GeometryReal >(1,0) , Point2D< GeometryReal >(1,1) , Point2D< GeometryReal >(0,1) };
				for( int p=2 ; p<polygon.size() ; p++ )
				{
					Point2D< GeometryReal > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };
					SquareMatrix< GeometryReal , 2 > fragment_to_element_differential;
					for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = dm[x][y];
					GeometryReal fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );

					for( unsigned int s=0 ; s<Samples ; s++ )
					{
						Point2D< GeometryReal > pos = polygon[0] + dm[0] * (GeometryReal)TriangleIntegrator< Samples >::Positions[s][0] + dm[1] * (GeometryReal)TriangleIntegrator< Samples >::Positions[s][1];

						typename ElementSamples::Bilinear::SampleData& _interior_cell_sample = interior_cell_samples[ (p-2)*Samples + s ];
						SetInteriorCellDuals< Samples >( _interior_cell_sample , Point2D< typename ElementSamples::Real >( pos ) , (typename ElementSamples::Real)( TriangleIntegrator< Samples >::Weights[s] * fragment_to_element_area_scale_factor / 2 ) );
						SetInteriorCellDuals< Samples >(  interior_cell_sample , Point2D< typename ElementSamples::Real >( pos ) , (typename ElementSamples::Real)( TriangleIntegrator< Samples >::Weights[s] * fragment_to_element_area_scale_factor / 2 ) );
						_interior_cell_sample.init( pos );
					}
				}
				interior_cell_sample.init( Point2D< GeometryReal >( (GeometryReal)0.5 , (GeometryReal)0.5 ) );
			}

			auto PolygonCenter = []( const Point2D< GeometryReal >* polygon , unsigned int polygonSize )
				{
					GeometryReal area = 0;
					Point2D< GeometryReal > center;
					for( int p=2 ; p<(int)polygonSize ; p++ )
					{
						Point2D< GeometryReal > v[] = { polygon[0] , polygon[p-1] , polygon[p] };
						Point2D< GeometryReal > d[] = { v[1]-v[0] , v[2]-v[0] };
						GeometryReal a = (GeometryReal)fabs( d[0][0] * d[1][1] - d[0][1] * d[1][0] );
						area += a;
						Point2D< GeometryReal > c = ( v[0] + v[1] + v[2] ) / 3;
						center += c * a; 
					}
					return center / area;
				};

			using Range = RegularGrid< 2 >::Range;
			using Index = RegularGrid< 2 >::Index;
			Range nodeRange , cellRange;
			nodeRange.second[0] = gridChart.width , cellRange.second[0] = gridChart.width-1;
			nodeRange.second[1] = gridChart.height , cellRange.second[1] = gridChart.height-1;


			auto GetSimplex = [&]( unsigned int t )
				{
					Simplex< double , 2 , 2 > simplex;
					for( unsigned int i=0 ; i<=2 ; i++ )
					{
						simplex[i] = atlasChart.vertex( ChartMeshVertexIndex( atlasChart.triangleIndex( ChartMeshTriangleIndex(t) )[i] ) ) - gridChart.corner;
						simplex[i][0] /= gridChart.cellSizeW , simplex[i][1] /= gridChart.cellSizeH;
					}
					return simplex;
				};

			std::vector< std::vector< std::vector< typename ElementSamples::Bilinear > > > _bilinear( ThreadPool::NumThreads() );
			std::vector< std::vector< typename ElementSamples::Quadratic > > _quadratic( ThreadPool::NumThreads() );
			for( unsigned int i=0 ; i<ThreadPool::NumThreads() ; i++ ) _bilinear[i].resize( elementSamples.bilinear.size() );

			ThreadPool::ParallelFor
			(
				0 , atlasChart.numTriangles() ,
				[&]( unsigned int thread , size_t t )
				{
					std::vector< std::vector< typename ElementSamples::Bilinear > > & bilinear = _bilinear[thread];
					std::vector< typename ElementSamples::Quadratic > & quadratic = _quadratic[thread];

					SimplexIndex< 2 , ChartMeshVertexIndex > tri = atlasChart.triangleIndex( ChartMeshTriangleIndex(t) );
					Point2D< GeometryReal > tPos[3];
					for( unsigned int k=0 ; k<=2 ; k++ ) tPos[k] = atlasChart.vertex( tri[k] ) - gridChart.corner;

					SquareMatrix< GeometryReal , 2 > texture_metric = texture_metrics[ ChartMeshTriangleIndex(t) ];
					SquareMatrix< GeometryReal , 2 > cell_metric = cell_to_texture_differential.transpose() * texture_metric * cell_to_texture_differential;
					SquareMatrix< GeometryReal , 2 > cell_metric_inverse = cell_metric.inverse();
					GeometryReal cell_area = (GeometryReal)sqrt( cell_metric.determinant() );

					std::vector< Point2D< GeometryReal > > parametricVertices(3);
					parametricVertices[0] = tPos[0], parametricVertices[1] = tPos[1], parametricVertices[2] = tPos[2];

					IndexedTriangle< GeometryReal > atlasTriangle;
					for( unsigned int k=0 ; k<3 ; k++ )
					{
						atlasTriangle.vertices[k] = tPos[k];
						atlasTriangle.atlasEdgeIndices[k] = atlasChart.atlasEdge( GetChartMeshHalfEdgeIndex( ChartMeshTriangleIndex(t) , k ) );
						atlasTriangle.vertexIndices[k] = tri[k];
						atlasTriangle.atlasVertexParentEdge[k] = AtlasMeshEdgeIndex(-1);
					}

					// Iterate over the cells that overlap the triangle
					auto Kernel = [&]( Index I )
						{
							int i = I[0] , j = I[1];
							auto TextureToCell = [&]( Point2D< GeometryReal > p ){ return Point2D< GeometryReal >( ( p[0] / gridChart.cellSizeW ) - i , ( p[1] / gridChart.cellSizeH ) - j ); };

							ChartInteriorCellIndex chartInteriorIndex = gridChart.cellIndices(i,j).interior;
							ChartBoundaryCellIndex chartBoundaryIndex = gridChart.cellIndices(i,j).boundary;
#ifdef SANITY_CHECK
							if( chartInteriorIndex!=ChartInteriorCellIndex(-1) && chartBoundaryIndex!=ChartBoundaryCellIndex(-1) ) MK_THROW( "Cell simultaneosly interior and boundary" );
#endif // SANITY_CHECK

							// Interior cells
							// If the cell is entirely interior to the triangle
							if( CellInTriangle( i , j , parametricVertices ) && chartInteriorIndex!=ChartInteriorCellIndex(-1) )
							{
								// For interior cells, the cell and the element are the same thing
								auto TextureToElement = TextureToCell;
								SquareMatrix< GeometryReal , 2 > element_metric = cell_metric , element_metric_inverse = cell_metric_inverse;
								GeometryReal element_area = cell_area;

								AtlasInteriorCellIndex atlasInteriorIndex = gridChart.chartToAtlasInteriorCellIndex( chartInteriorIndex );
								unsigned int cellLineId = interiorCellLineIndex[atlasInteriorIndex].first;
								unsigned int cellLineOffset = interiorCellLineIndex[atlasInteriorIndex].second;

								typename ElementSamples::Bilinear bilinearElementSample( fastIntegration ? 1 : 2*Samples );
								bilinearElementSample.cellOffset = cellLineOffset;
								for( unsigned int x=0 ; x<2 ; x++ ) for( unsigned int y=0 ; y<2 ; y++ ) bilinearElementSample.tensor(x,y) = (typename ElementSamples::Real)element_metric_inverse(x,y);

								if( fastIntegration )
								{
									bilinearElementSample[0] = interior_cell_sample;
									SetCellInTriangleDuals< Samples >( bilinearElementSample[0] , bilinearElementSample , interior_cell_sample , (typename ElementSamples::Real)element_area/2 );
								}
								else
								{
									for( unsigned int s=0 ; s<2*Samples ; s++ )
									{
										bilinearElementSample[s] = interior_cell_samples[s];
										SetCellInTriangleDuals< Samples >( bilinearElementSample[s] , bilinearElementSample , interior_cell_samples[s] , (typename ElementSamples::Real)element_area/2 );
									}
								}
								bilinear[ cellLineId ].push_back( bilinearElementSample );
							}
							else if( chartInteriorIndex!=ChartInteriorCellIndex(-1) )
							{
								// For interior cells, the cell and the element are the same thing
								auto TextureToElement = TextureToCell;

								SquareMatrix< GeometryReal , 2 > element_metric = cell_metric , element_metric_inverse = cell_metric_inverse;
								GeometryReal element_area = cell_area;

								CellClippedTriangle< GeometryReal > polygon = parametricVertices;

								// Clip the triangle to the cell
								if( ClipTriangleToPrimalCell( polygon , i , j , gridChart.cellSizeW , gridChart.cellSizeH ) )
								{
									// Transform the polygon vertices into the coordinate frame of the cell
									for( unsigned int ii=0 ; ii<polygon.size() ; ii++ ) polygon[ii] = TextureToElement( polygon[ii] );

									AtlasInteriorCellIndex atlasInteriorIndex = gridChart.chartToAtlasInteriorCellIndex( chartInteriorIndex );
									unsigned int cellLineId = interiorCellLineIndex[atlasInteriorIndex].first;
									unsigned int cellLineOffset = interiorCellLineIndex[atlasInteriorIndex].second;

									// There is a single sample for the whole polygon
									typename ElementSamples::Bilinear bilinearElementSample( fastIntegration ? 1 : (polygon.size()-2)*Samples );
									bilinearElementSample.cellOffset = cellLineOffset;
									for( unsigned int x=0 ; x<2 ; x++ ) for( unsigned int y=0 ; y<2 ; y++ ) bilinearElementSample.tensor(x,y) = (typename ElementSamples::Real)element_metric_inverse(x,y);

									for( unsigned int p=2 ; p<polygon.size() ; p++ )
									{
										Point2D< GeometryReal > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };
										SquareMatrix< GeometryReal , 2 > fragment_to_element_differential;
										for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = dm[x][y];
										GeometryReal fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );
										GeometryReal fragment_area = element_area * fragment_to_element_area_scale_factor / 2;

										for( int s=0 ; s<Samples ; s++ )
										{
											Point2D< typename ElementSamples::Real > pos = polygon[0] + dm[0] * (typename ElementSamples::Real)TriangleIntegrator< Samples >::Positions[s][0] + dm[1] * (typename ElementSamples::Real)TriangleIntegrator< Samples >::Positions[s][1];
#ifdef SANITY_CHECK
											if( !InUnitSquare( pos ) ) MK_THROW( "Sample position out of unit square! (" , pos[0] , " " , pos[1] , ")" );
#endif // SANITY_CHECK

											typename ElementSamples::Bilinear::SampleData& sampleData = bilinearElementSample[ fastIntegration ? 0 : (p-2)*Samples+s ];
											SetInteriorDuals< Samples >( sampleData , bilinearElementSample , pos , (typename ElementSamples::Real)( TriangleIntegrator< Samples >::Weights[s] * fragment_area ) );
											sampleData.init( pos );
										}
									}
									if( fastIntegration ) bilinearElementSample[0].init( Point2D< GeometryReal >( PolygonCenter( &polygon[0] , (unsigned int)polygon.size() ) ) );
									bilinear[ cellLineId ].push_back( bilinearElementSample );
								}
							}
							// Boundary cell
							else if( chartBoundaryIndex!=ChartBoundaryCellIndex(-1) )
							{
								const std::vector< BoundaryIndexedTriangle< GeometryReal > > & cellBoundaryTriangles = gridChart.boundaryTriangles[chartBoundaryIndex];

								// Iterate over all elements in the cell
								for( unsigned int bt=0 ; bt<cellBoundaryTriangles.size() ; bt++ )
								{
									BoundaryIndexedTriangle< GeometryReal > element = cellBoundaryTriangles[bt];
									std::vector< Point2D< GeometryReal > > element_vertices(3);
									for( unsigned int ii=0 ; ii<3 ; ii++ ) element_vertices[ii] = TextureToCell( element[ii] );


									IndexedPolygon< GeometryReal > polygon;
									SetIndexedPolygonFromBoundaryTriangle( element , polygon );

									// Intersect the element with the atlas triangle
									if( ClipPartiallyIndexedPolygonToIndexedTriangle( polygon , atlasTriangle ) )
									{
										// Convert the polygon vertices from the texture frame to the cell frame
										for( unsigned int ii=0 ; ii<polygon.size() ; ii++ ) polygon[ii] = TextureToCell( polygon[ii] );

										SquareMatrix< GeometryReal , 2 > element_to_cell_differential , cell_to_element_differential;
										Point2D< GeometryReal > dm[] = { element_vertices[1]-element_vertices[0] , element_vertices[2]-element_vertices[0] };
										for( unsigned int x=0 ; x<2 ; x++ ) for( unsigned int y=0 ; y<2 ; y++ ) element_to_cell_differential(x,y) = dm[x][y];
										cell_to_element_differential = element_to_cell_differential.inverse();
										auto CellToElement = [&]( Point2D< GeometryReal > p ){ return cell_to_element_differential * ( p - element_vertices[0] ); };
										// Convert the polygon vertices from the cell frame to the element frame
										for( unsigned int ii=0 ; ii<polygon.size() ; ii++ ) polygon[ii] = CellToElement( polygon[ii] );

										SquareMatrix< GeometryReal , 2 > element_metric = element_to_cell_differential.transpose() * cell_metric * element_to_cell_differential;
										SquareMatrix< GeometryReal , 2 > element_metric_inverse = element_metric.inverse();
										GeometryReal element_area = sqrt( element_metric.determinant() );
										typename ElementSamples::Quadratic quadraticElementSample( fastIntegration ? 1 : (unsigned int)(polygon.size()-2)*Samples );
										for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) quadraticElementSample.tensor(x,y) = (typename ElementSamples::Real)element_metric_inverse(x,y);
										const QuadraticElement::Index& triangleElementIndices = gridChart.boundaryTriangles[chartBoundaryIndex][bt].indices;
										for( unsigned int k=0 ; k<6 ; k++ )
										{
											AtlasInteriorOrBoundaryNodeIndex _fineBoundaryIndex = fineBoundaryIndex[ static_cast< unsigned int >( triangleElementIndices[k] ) ];
#ifdef SANITY_CHECK
											if( _fineBoundaryIndex!=AtlasInteriorOrBoundaryNodeIndex(-1) ) quadraticElementSample.fineNodes[k] = _fineBoundaryIndex;
											else MK_THROW( "Invalid fine boundary index" );
#else // !SANITY_CHECK
											quadraticElementSample.fineNodes[k] = _fineBoundaryIndex;
#endif // SANITY_CHECK
										}

										for( unsigned int p=2 ; p<polygon.size() ; p++ )
										{
											Point2D< GeometryReal > d[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };
											SquareMatrix< GeometryReal , 2 > fragment_to_element_differential;
											for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = d[x][y];
											GeometryReal fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );
											GeometryReal fragment_area = element_area * fragment_to_element_area_scale_factor / 2;

											if( fragment_area>0 )
											{
												for( int s=0 ; s<Samples ; s++ )
												{
													Point2D< GeometryReal > pos = polygon[0] + d[0] * (GeometryReal)TriangleIntegrator< Samples >::Positions[s][0] + d[1] * (GeometryReal)TriangleIntegrator< Samples >::Positions[s][1];
#ifdef SANITY_CHECK
													if( !InUnitTriangle( pos ) ) MK_THROW( "Sample out of unit right triangle! (" , pos[0] , " " , pos[1] , ")" );
													else
													{
														pos[0] = std::max< GeometryReal >( pos[0] , 0 );
														pos[1] = std::max< GeometryReal >( pos[1] , 0 );
														GeometryReal excess = ( pos[0] + pos[1] ) - 1;
														if( excess>0 ) pos[0] -= excess/2 , pos[1] -= excess/2;
													}
#else // !SANITY_CHECK
													pos[0] = std::max< GeometryReal >( pos[0] , 0 );
													pos[1] = std::max< GeometryReal >( pos[1] , 0 );
													GeometryReal excess = ( pos[0] + pos[1] ) - 1;
													if( excess>0 ) pos[0] -= excess/2 , pos[1] -= excess/2;
#endif // SANITY_CHECK

													typename ElementSamples::Quadratic::SampleData& sampleData = quadraticElementSample[ fastIntegration ? 0 : (p-2)*Samples+s ];
													SetBoundaryDuals< Samples >( sampleData , quadraticElementSample , Point2D< typename ElementSamples::Real >( pos ) , (typename ElementSamples::Real)( TriangleIntegrator< Samples >::Weights[s] * fragment_area ) );
													sampleData.init( pos );
												}
											}
											else MK_WARN( "Element discarded due to zero mass. Triangle " , t , ". Boundary cell " , chartBoundaryIndex , ". Element " , bt , ". Sub triangle " , p-2 );
										}
										if( fastIntegration ) quadraticElementSample[0].init( Point2D< GeometryReal >( PolygonCenter( &polygon[0] , (unsigned int)polygon.size() ) ) );
										quadratic.push_back( quadraticElementSample );
									}
								}
							}
						};
					Rasterizer2D::RasterizeSupports< true , true >( GetSimplex( static_cast< unsigned int >(t) ) , Kernel , cellRange );
				}
			);
			for( unsigned int i=0 ; i<elementSamples.bilinear.size() ; i++ )
			{
				size_t sz = elementSamples.bilinear[i].size();
				for( unsigned int t=0 ; t<_bilinear.size() ; t++ ) sz += _bilinear[t][i].size();
				elementSamples.bilinear[i].reserve( sz );
				for( unsigned int t=0 ; t<_bilinear.size() ; t++ ) elementSamples.bilinear[i].insert( elementSamples.bilinear[i].end() , _bilinear[t][i].begin() , _bilinear[t][i].end() );
			}
			{
				size_t sz = elementSamples.quadratic.size();
				for( unsigned int t=0 ; t<_quadratic.size() ; t++ ) sz += _quadratic[t].size();
				elementSamples.quadratic.reserve( sz );
				for( unsigned int t=0 ; t<_quadratic.size() ; t++ ) elementSamples.quadratic.insert( elementSamples.quadratic.end() , _quadratic[t].begin() , _quadratic[t].end() );
			}
		}

		template< unsigned int Samples , typename GeometryReal , typename ElementSamples >
		void InitializeIntegration
		(
			const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > >& parameterMetric ,
			const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > &atlasCharts ,
			const ExplicitIndexVector< ChartIndex , GridChart< GeometryReal > > &gridCharts ,
			const ExplicitIndexVector< AtlasInteriorCellIndex , std::pair< unsigned int , unsigned int > > &interiorCellLineIndex ,
			const std::vector< AtlasInteriorOrBoundaryNodeIndex > &fineBoundaryIndex ,
			ElementSamples &elementSamples ,
			bool fastIntegration
		)
		{
			for( unsigned int i=0 ; i<gridCharts.size() ; i++ )
				InitializeIntegration< Samples >( parameterMetric[ ChartIndex(i) ] , atlasCharts[ ChartIndex(i) ] , gridCharts[ ChartIndex(i) ] , interiorCellLineIndex , fineBoundaryIndex , elementSamples , fastIntegration );
		}

		///////////////////////////
		// Compute the integrals //
		///////////////////////////
		template< class Real , typename OutT , typename InT , typename ValueFunctionType >
		void IntegrateBilinear
		(
			const BilinearElementScalarSample< Real > &sample ,
			const ValueFunctionType &ValueFunction ,
			const InT cornerValues[] ,
			OutT rhsValues[]
		)
		{
			for( unsigned int s=0 ; s<sample.size() ; s++ )
			{
				OutT scalar = ValueFunction( sample[s]( cornerValues ) , sample.tensor );
				for( unsigned int k=0 ; k<4 ; k++ ) rhsValues[k] += scalar * sample[s].dualValues[k];
			}
		}

		template< class Real , typename OutT , typename InT , typename GradientFunctionType >
		void IntegrateBilinear
		(
			const BilinearElementGradientSample< Real > &sample ,
			const GradientFunctionType &GradientFunction ,
			const InT cornerValues[] ,
			OutT rhsValues[]
		)
		{
			for( unsigned int s=0 ; s<sample.size() ; s++ )
			{
				Point2D< OutT > gradientVector = GradientFunction( sample[s]( cornerValues ) , sample.tensor );
				for( unsigned int k=0 ; k<4 ; k++ ) for( int d=0 ; d<2 ; d++ ) rhsValues[k] += gradientVector[d] * sample[s].dualGradients[k][d];
			}
		}

		template< class Real , typename OutT , typename InT , typename ValueFunctionType >
		void IntegrateQuadratic
		(
			const QuadraticElementScalarSample< Real >& sample ,
			const ValueFunctionType& ValueFunction ,
			const InT cornerValues[] ,
			OutT rhsValues[]
		)
		{
			for( unsigned int s=0 ; s<sample.size() ; s++ )
			{
				OutT scalar = ValueFunction( sample[s]( cornerValues ) , sample.tensor );
				for( unsigned int k=0 ; k<6 ; k++ ) rhsValues[k] += scalar * sample[s].dualValues[k];
			}
		}

		template< class Real , typename OutT , typename InT , typename GradientFunctionType >
		void IntegrateQuadratic
		(
			const QuadraticElementGradientSample< Real >& sample ,
			const GradientFunctionType& VectorFunction ,
			const InT cornerValues[] ,
			OutT rhsValues[]
		)
		{
			for( unsigned int s=0 ; s<sample.size() ; s++ )
			{
				Point2D< OutT > gradientVector = VectorFunction( sample[s]( cornerValues ) , sample.tensor );
				for( unsigned int k=0 ; k<6 ; k++ ) for( unsigned int d=0 ; d<2 ; d++ ) rhsValues[k] += gradientVector[d] * sample[s].dualGradients[k][d];
			}
		}

		template< typename Real , typename OutT , typename InT , typename ElementSamples , typename SampleFunctionType >
		void Integrate
		(
			const std::vector< InteriorCellLine >& interiorCellLines ,
			const ElementSamples &elementSamples ,
			const std::vector< InT >& potential ,
			const std::vector< InT >& boundary_potential ,
			const SampleFunctionType& SampleFunction ,
			std::vector< OutT >& rhs ,
			std::vector< OutT >& boundary_rhs ,
			bool verbose=false
		)
		{
			return Integrate< Real >( interiorCellLines , elementSamples , GetPointer( potential ) , GetPointer( boundary_potential ) , SampleFunction , GetPointer( rhs ) , GetPointer( boundary_rhs ) , verbose );
		}

		template< typename Real , typename OutT , typename InT , typename ElementSamples , typename SampleFunctionType >
		void Integrate
		(
			const std::vector< InteriorCellLine >& interiorCellLines ,
			const ElementSamples &elementSamples ,
			ConstPointer( InT ) potential ,
			ConstPointer( InT ) boundary_potential ,
			const SampleFunctionType& SampleFunction ,
			Pointer( OutT ) rhs ,
			Pointer( OutT ) boundary_rhs ,
			bool verbose=false
		)
		{
			Miscellany::Timer timer;
			auto UpdateRow = [&]( unsigned int r )
				{
					const InT * _inPrevious = &potential[0] + interiorCellLines[r].prevLineIndex;
					const InT *     _inNext = &potential[0] + interiorCellLines[r].nextLineIndex;

					OutT * _outPrevious = &rhs[0] + interiorCellLines[r].prevLineIndex;
					OutT *     _outNext = &rhs[0] + interiorCellLines[r].nextLineIndex;

					InT cornerValues[4];
					cornerValues[0] = *_inPrevious;
					_inPrevious++;
					cornerValues[1] = *_inPrevious;
					cornerValues[3] = *_inNext;
					_inNext++;
					cornerValues[2] = *_inNext;

					OutT rhsValues[] = { OutT() , OutT() , OutT() , OutT() };

					unsigned int numSamples = (unsigned int)elementSamples.bilinear[r].size();
					const typename ElementSamples::Bilinear *samplePtr = &elementSamples.bilinear[r][0];
					unsigned int currentOffset = 0;

					for( unsigned int j=0 ; j<numSamples ; j++ )
					{
						const typename ElementSamples::Bilinear &sample = *samplePtr;
						if( sample.cellOffset>currentOffset )
						{
							cornerValues[0] = cornerValues[1];
							_inPrevious++;
							cornerValues[1] = *_inPrevious;
							cornerValues[3] = cornerValues[2];
							_inNext++;
							cornerValues[2] = *_inNext;

							Atomic< OutT >::Add( *_outPrevious , rhsValues[0] );
							_outPrevious++;
							Atomic< OutT >::Add( *_outNext , rhsValues[3] );
							_outNext++;
							rhsValues[0] = rhsValues[1];
							rhsValues[1] = OutT();
							rhsValues[3] = rhsValues[2];
							rhsValues[2] = OutT();

							currentOffset++;
						}

						IntegrateBilinear< Real , OutT , InT >( sample , SampleFunction , cornerValues , rhsValues );
						samplePtr++;
					}

					Atomic< OutT >::Add( *_outPrevious , rhsValues[0] );
					_outPrevious++;
					Atomic< OutT >::Add( *_outPrevious , rhsValues[1] );
					Atomic< OutT >::Add( *_outNext     , rhsValues[3] );
					_outNext++;
					Atomic< OutT >::Add( *_outNext     , rhsValues[2] );
				};
			ThreadPool::ParallelFor( 0 , interiorCellLines.size() , [&]( unsigned int , size_t r ){ UpdateRow( static_cast< unsigned int >(r) ); } );

			if( verbose ) printf( "Integrated bilinear: %.2f(s)\n" , timer.elapsed() );

			timer.reset();
			for( unsigned int i=0 ; i<elementSamples.quadratic.size() ; i++ )
			{
				const typename ElementSamples::Quadratic &sample = elementSamples.quadratic[i];
				// The values of the potential at the vertices and edge mid-points
				InT cornerValues[] =
				{
					boundary_potential[ static_cast< unsigned int >(sample.fineNodes[0]) ] ,
					boundary_potential[ static_cast< unsigned int >(sample.fineNodes[1]) ] ,
					boundary_potential[ static_cast< unsigned int >(sample.fineNodes[2]) ] ,
					boundary_potential[ static_cast< unsigned int >(sample.fineNodes[3]) ] ,
					boundary_potential[ static_cast< unsigned int >(sample.fineNodes[4]) ] ,
					boundary_potential[ static_cast< unsigned int >(sample.fineNodes[5]) ]
				};
				OutT rhsValues[] = { OutT() , OutT() , OutT() , OutT() , OutT() , OutT() };
				IntegrateQuadratic< Real , OutT , InT >( sample , SampleFunction , cornerValues , rhsValues );
				for( unsigned int k=0 ; k<6 ; k++ ) boundary_rhs[ static_cast< unsigned int >( sample.fineNodes[k] ) ] += rhsValues[k];
			}
			if( verbose ) printf( "Integrated quadratic: %.2f(s)\n" , timer.elapsed() );
		}
	}
}