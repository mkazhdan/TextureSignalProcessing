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


#include <functional>
#include "Hierarchy.h"

class InteriorCellLine
{
public:
	int prevLineIndex;
	int nextLineIndex;
	int length;
};

int InitializeGridChartInteriorCellLines
(
	const AtlasChart& atlasChart ,
	const GridChart& gridChart ,
	std::vector< InteriorCellLine >& interiorCellLines ,
	std::vector< std::pair< int , int > >& interiorCellLineIndex
)
{
	const Image< int >& cellType = gridChart.cellType;
	int width = cellType.width();
	int height = cellType.height();

	const Image<int> & globalTexelIndex = gridChart.globalTexelIndex;
	const Image<int> & globalTexelInteriorIndex = gridChart.globalInteriorTexelIndex;

	int localInteriorCellIndex = 0;

	for( int j=0 ; j<height ; j++ )
	{
		int offset = 0;
		bool previousIsInterior = false;
		int rasterStart = -1;
		while (offset < width) {
			bool currentIsInterior = cellType(offset, j) == 1;
			if ((offset == 0 || offset == width - 1) && currentIsInterior) {
				printf("Unexpected interior cell! \n");
				return 0;
			}
			if (currentIsInterior && !previousIsInterior) rasterStart = offset; //Start raster line
			if (!currentIsInterior && previousIsInterior) { //Terminate raster line
				InteriorCellLine newLine;
				newLine.prevLineIndex = globalTexelIndex(rasterStart, j);
				newLine.nextLineIndex = globalTexelIndex(rasterStart, j + 1);
				newLine.length = offset - rasterStart;

				if (newLine.prevLineIndex == -1 || newLine.nextLineIndex == -1) {
					printf("Inavlid Indexing! \n");
					return 0;
				}
				int currentLine = (int)interiorCellLines.size();

				for (int k = 0; k < offset - rasterStart; k++) {
					if (gridChart.interiorCellCorners[localInteriorCellIndex][0] != globalTexelInteriorIndex(rasterStart + k, j)) {
						printf("Unexpected corner id! \n");
						return 0;
					}

					interiorCellLineIndex.push_back(std::pair<int, int>(currentLine, k));
					localInteriorCellIndex++;
				}

				interiorCellLines.push_back(newLine);
			}
			previousIsInterior = currentIsInterior;
			offset++;
		}
	}

	return 1;
}

int InitializeGridAtlasInteriorCellLines
(
	const std::vector< AtlasChart >& atlasCharts ,
	const std::vector< GridChart >& gridCharts ,
	std::vector< InteriorCellLine >& interiorCellLines ,
	std::vector< std::pair< int , int > >& interiorCellLineIndex
)
{
	for( int i=0 ; i<gridCharts.size() ; i++ ) InitializeGridChartInteriorCellLines( atlasCharts[i] , gridCharts[i] , interiorCellLines , interiorCellLineIndex );
	return 1;
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
	const BilinearElementScalarSample< Real >& sample ,
	const typename BilinearElementScalarSample< Real >::SampleData& interiorSampleData ,
	Real element_area
)
{
	for( int k=0 ; k<4 ; k++ ) sampleData.dualValues[k] = interiorSampleData.dualValues[k] * element_area;
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
	for( int k=0 ; k<4 ; k++ ) sampleData.dualGradients[k] = sample.tensor * interiorSampleData.dualGradients[k] * element_area;
}

template< unsigned int Samples , typename Real >
void SetInteriorDuals
(
	typename BilinearElementScalarSample< Real >::SampleData& sampleData ,
	const BilinearElementScalarSample< Real >& sample ,
	Point2D< Real > pos , 
	Real fragment_quadrature_weight
)
{
	// Compute the integrated gradients of each bilinear basis functions in the frame of the cell, weighted by the area of the unit triangle
	// Dualize the samples so that integration is simple a dot-product
	for( int k=0 ; k<4 ; k++ ) sampleData.dualValues[k] += (Real)BilinearElementValue( k , pos ) * fragment_quadrature_weight;
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
	for( int k=0 ; k<4 ; k++ ) sampleData.dualGradients[k] += sample.tensor * Point2D< Real >( BilinearElementGradient( k , pos ) * fragment_quadrature_weight );
}

template< unsigned int Samples , typename Real >
void SetBoundaryDuals
(
	typename QuadraticElementScalarSample< Real >::SampleData& sampleData ,
	const QuadraticElementScalarSample< Real >& sample ,
	Point2D< Real > pos , 
	Real fragment_quadrature_weight
)
{
	// Compute the integrated gradients of each quadratic basis functions in the frame of the cell, weighted by the area of the unit triangle
	// Dualize the samples so that integration is simple a dot-product
	for( int k=0 ; k<6 ; k++ ) sampleData.dualValues[k] += (Real)QuadraticElementValue( k , pos ) * fragment_quadrature_weight;
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
	for( int k=0 ; k<6 ; k++ ) sampleData.dualGradients[k] += sample.tensor * Point2D< Real >( QuadraticElementGradient( k , pos ) * fragment_quadrature_weight );
}


template< unsigned int Samples , typename Real , typename BilinearElementSample , typename QuadraticElementSample >
int InitializeIntegration
(
	const std::vector< SquareMatrix< double , 2 > >& texture_metrics ,
	const AtlasChart& atlasChart ,
	const GridChart& gridChart ,
	const std::vector< std::pair< int , int > >& interiorCellLineIndex ,
	const std::vector< int >& fineBoundaryIndex ,
	std::vector< std::vector< BilinearElementSample > >& bilinearElementSamples ,
	std::vector< QuadraticElementSample >& quadraticElementSamples ,
	bool fastIntegration
)
{
	////Rasterize
	const double PRECISION_ERROR = 1e-3;

	auto InUnitSquare =   [&]( Point2D< double > p ){ return !( p[0]<0-PRECISION_ERROR || p[1]<0-PRECISION_ERROR || p[0]>1+PRECISION_ERROR || p[1]>1+PRECISION_ERROR ); };
	auto InUnitTriangle = [&]( Point2D< double > p ){ return !( p[0]<0-PRECISION_ERROR || p[1]<0-PRECISION_ERROR || ( p[0]+p[1] )>1+PRECISION_ERROR ); };
	auto CellInTriangle = [&]( int i , int j , const std::vector< Point2D< double > >& vertices )
	{
		double x1 = (double)i*gridChart.cellSizeW , x2 = (double)(i+1)*gridChart.cellSizeW;
		double y1 = (double)j*gridChart.cellSizeH , y2 = (double)(j+1)*gridChart.cellSizeH;
		Point2D< double > points[] = { Point2D< double >(x1,y1) , Point2D< double >(x2,y1) , Point2D< double >(x2,y2) , Point2D< double >(x1,y2) };
		Point2D< double > normals[] = { vertices[1]-vertices[0] , vertices[2]-vertices[1] , vertices[0]-vertices[2] };
		for( int k=0 ; k<3 ; k++ ) normals[k] = Point2D< double >( normals[k][1] , -normals[k][0] );
		double offsets[] = { Point2D< double >::Dot( normals[0] , vertices[0] ) , Point2D< double >::Dot( normals[1] , vertices[1] ) , Point2D< double >::Dot( normals[2] , vertices[2] ) };
		for( int k=0 ; k<4 ; k++ )
		{
			if( ( Point2D< double >::Dot( points[k] , normals[0] ) - offsets[0] )<0 ) return false;
			if( ( Point2D< double >::Dot( points[k] , normals[1] ) - offsets[1] )<0 ) return false;
			if( ( Point2D< double >::Dot( points[k] , normals[2] ) - offsets[2] )<0 ) return false;
		}
		return true;
	};

	SquareMatrix< double , 2 > cell_to_texture_differential;
	cell_to_texture_differential(0,0) = gridChart.cellSizeW;
	cell_to_texture_differential(1,1) = gridChart.cellSizeH;
	cell_to_texture_differential(0,1) = cell_to_texture_differential(1,0) = 0.;

	typename BilinearElementSample::SampleData interior_cell_samples[2*Samples] , interior_cell_sample;
	{
		std::vector< Point2D< double > > polygon = { Point2D< double >(0,0) , Point2D< double >(1,0) , Point2D< double >(1,1) , Point2D< double >(0,1) };
		for( int p=2 ; p<polygon.size() ; p++ )
		{
			Point2D< double > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };
			SquareMatrix< double , 2 > fragment_to_element_differential;
			for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = dm[x][y];
			double fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );

			for( int s=0 ; s<Samples ; s++ )
			{
				Point2D< Real > pos = polygon[0] + dm[0] * TriangleIntegrator< Samples >::Positions[s][0] + dm[1] * TriangleIntegrator< Samples >::Positions[s][1];

				typename BilinearElementSample::SampleData& _interior_cell_sample = interior_cell_samples[ (p-2)*Samples + s ];
				SetInteriorCellDuals< Samples , Real >( _interior_cell_sample , pos , (Real)( TriangleIntegrator< Samples >::Weights[s] * fragment_to_element_area_scale_factor / 2.0 ) );
				SetInteriorCellDuals< Samples , Real >(  interior_cell_sample , pos , (Real)( TriangleIntegrator< Samples >::Weights[s] * fragment_to_element_area_scale_factor / 2.0 ) );
				_interior_cell_sample.pos = pos;
			}
		}
		interior_cell_sample.pos = Point2D< Real >( (Real)0.5 , (Real)0.5 );
	}

	auto PolygonCenter = []( const Point2D< double >* polygon , unsigned int polygonSize )
	{
		Real area = 0;
		Point2D< Real > center;
		for( int p=2 ; p<(int)polygonSize ; p++ )
		{
			Point2D< Real > v[] = { polygon[0] , polygon[p-1] , polygon[p] };
			Point2D< Real > d[] = { v[1]-v[0] , v[2]-v[0] };
			Real a = (Real)fabs( d[0][0] * d[1][1] - d[0][1] * d[1][0] );
			area += a;
			Point2D< Real > c = ( v[0] + v[1] + v[2] ) / (Real)3.;
			center += c * a; 
		}
		return center / area;
	};	

	for( int t=0 ; t<atlasChart.triangles.size() ; t++ )
	{
		Point2D< double > tPos[3];
		for( int i=0 ; i<3 ; i++ ) tPos[i] = atlasChart.vertices[ atlasChart.triangles[t][i] ] - gridChart.corner;

		SquareMatrix< double , 2 > texture_metric = texture_metrics[t];
		SquareMatrix< double , 2 > cell_metric = cell_to_texture_differential.transpose() * texture_metric * cell_to_texture_differential;
		SquareMatrix< double , 2 > cell_metric_inverse = cell_metric.inverse();
		double cell_area = sqrt( cell_metric.determinant() );

		//BBox
		int minCorner[2] , maxCorner[2];
		GetTriangleIntegerBBox( tPos , 1.0/gridChart.cellSizeW , 1.0/gridChart.cellSizeH , minCorner , maxCorner );

		std::vector< Point2D< double > > parametricVertices(3);
		parametricVertices[0] = tPos[0], parametricVertices[1] = tPos[1], parametricVertices[2] = tPos[2];

		AtlasIndexedTriangle atlasTriangle;
		for( int k=0 ; k<3 ; k++ )
		{
			atlasTriangle.vertices[k] = tPos[k];
			atlasTriangle.atlasEdgeIndices[k] = atlasChart.atlasEdgeIndices[ 3*t+k ];
			atlasTriangle.atlasVertexIndices[k] = atlasChart.triangles[t][k];
			atlasTriangle.atlasVertexParentEdge[k] = -1;
		}

		// Iterate over the cells that can overlap the triangle
		for( int j=minCorner[1] ; j<maxCorner[1] ; j++ ) for( int i=minCorner[0] ; i<maxCorner[0] ; i++ )
		{
			auto TextureToCell = [&]( Point2D< double > p ){ return Point2D< double >( ( p[0] / gridChart.cellSizeW ) - i , ( p[1] / gridChart.cellSizeH ) - j ); };

			int localInteriorIndex = gridChart.localInteriorCellIndex(i,j) , localBoundaryIndex = gridChart.localBoundaryCellIndex(i,j);
			if( localInteriorIndex!=-1 && localBoundaryIndex!=-1 ){ printf( "[ERROR] Cell simultaneosly interior and boundary!\n" ) ; return 0;	}

			// Interior cells
			// If the cell is entirely interior to the triangle
			if( CellInTriangle( i , j , parametricVertices ) && localInteriorIndex!=-1 )
			{
				// For interior cells, the cell and the element are the same thing
				auto TextureToElement = TextureToCell;
				SquareMatrix< double , 2 > element_metric = cell_metric , element_metric_inverse = cell_metric_inverse;
				double element_area = cell_area;

				int globalInteriorIndex = localInteriorIndex + gridChart.globalIndexInteriorCellOffset;
				int cellLineId = interiorCellLineIndex[globalInteriorIndex].first;
				int cellLineOffset = interiorCellLineIndex[globalInteriorIndex].second;

				BilinearElementSample bilinearElementSample( fastIntegration ? 1 : 2*Samples );
				bilinearElementSample.cellOffset = cellLineOffset;
				for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) bilinearElementSample.tensor(x,y) = (Real)element_metric_inverse(x,y);

				if( fastIntegration )
				{
					{
						bilinearElementSample[0].pos = interior_cell_sample.pos;
						SetCellInTriangleDuals< Samples , Real >( bilinearElementSample[0] , bilinearElementSample , interior_cell_sample , (Real)( element_area/2. ) );
					}
				}
				else
				{
					for( int s=0 ; s<2*Samples ; s++ )
					{
						bilinearElementSample[s].pos = interior_cell_samples[s].pos;
						SetCellInTriangleDuals< Samples , Real >( bilinearElementSample[s] , bilinearElementSample , interior_cell_samples[s] , (Real)( element_area/2. ) );
					}
				}
#pragma omp critical
				bilinearElementSamples[ cellLineId ].push_back( bilinearElementSample );
			}
			else if( localInteriorIndex!=-1 )
			{
				// For interior cells, the cell and the element are the same thing
				auto TextureToElement = TextureToCell;
				SquareMatrix< double , 2 > element_metric = cell_metric , element_metric_inverse = cell_metric_inverse;
				double element_area = cell_area;

				CellClippedTriangle polygon = parametricVertices;

				// Clip the triangle to the cell
				if( ClipTriangleToPrimalCell( polygon , i , j , gridChart.cellSizeW , gridChart.cellSizeH ) )
				{
					// Transform the polygon vertices into the coordinate frame of the cell
					for( int i=0 ; i<polygon.size() ; i++ ) polygon[i] = TextureToElement( polygon[i] );

					int globalInteriorIndex = localInteriorIndex + gridChart.globalIndexInteriorCellOffset;
					int cellLineId = interiorCellLineIndex[globalInteriorIndex].first;
					int cellLineOffset = interiorCellLineIndex[globalInteriorIndex].second;

					// There is a single sample for the whole polygon
					BilinearElementSample bilinearElementSample( fastIntegration ? 1 : (polygon.size()-2)*Samples );
					bilinearElementSample.cellOffset = cellLineOffset;
					for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) bilinearElementSample.tensor(x,y) = (Real)element_metric_inverse(x,y);

					for( int p=2 ; p<polygon.size() ; p++ )
					{
						Point2D< double > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };
						SquareMatrix< double , 2 > fragment_to_element_differential;
						for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = dm[x][y];
						double fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );
						double fragment_area = element_area * fragment_to_element_area_scale_factor / 2.;

						for( int s=0 ; s<Samples ; s++ )
						{
							Point2D< Real > pos = polygon[0] + dm[0] * TriangleIntegrator< Samples >::Positions[s][0] + dm[1] * TriangleIntegrator< Samples >::Positions[s][1];
							if( !InUnitSquare( pos ) ){ printf( "[ERROR] Sample position out of unit square! (%f %f)\n" , pos[0] , pos[1] ) ; return 0; }

							typename BilinearElementSample::SampleData& sampleData = bilinearElementSample[ fastIntegration ? 0 : (p-2)*Samples+s ];
							SetInteriorDuals< Samples , Real >( sampleData , bilinearElementSample , pos , (Real)( TriangleIntegrator< Samples >::Weights[s] * fragment_area ) );
							sampleData.pos = pos;
						}
					}
					if( fastIntegration ) bilinearElementSample[0].pos = Point2D< Real >( PolygonCenter( &polygon[0] , (unsigned int)polygon.size() ) );
#pragma omp critical
					bilinearElementSamples[ cellLineId ].push_back( bilinearElementSample );
				}
			}
			// Boundary cell
			else if( localBoundaryIndex!=-1 )
			{
				std::vector< BoundaryIndexedTriangle > cellBoundaryTriangles = gridChart.boundaryTriangles[localBoundaryIndex];

				// Iterate over all elements in the cell
				for( int bt=0 ; bt<cellBoundaryTriangles.size() ; bt++ )
				{
					BoundaryIndexedTriangle element = cellBoundaryTriangles[bt];
					std::vector< Point2D< double > > element_vertices(3);
					for( int i=0 ; i<3 ; i++ ) element_vertices[i] = TextureToCell( element[i] );
					int boundaryTriangleId = element.id;

					AtlasIndexedPolygon polygon;
					SetAtlasIndexedPolygonFromBoundaryTriangle( element , polygon );

					// Intersect the element with the atlas triangle
					int clippingResult = ClipPartiallyIndexedPolygonToIndexedTriangle( polygon , atlasTriangle );
					if( clippingResult>0 )
					{
						// Convert the polygon vertices from the texture frame to the cell frame
						for( int i=0 ; i<polygon.size() ; i++ ) polygon[i] = TextureToCell( polygon[i] );

						SquareMatrix< double , 2 > element_to_cell_differential , cell_to_element_differential;
						Point2D< double > dm[] = { element_vertices[1]-element_vertices[0] , element_vertices[2]-element_vertices[0] };
						for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) element_to_cell_differential(x,y) = dm[x][y];
						cell_to_element_differential = element_to_cell_differential.inverse();
						auto CellToElement = [&]( Point2D< double > p ){ return cell_to_element_differential * ( p - element_vertices[0] ); };
						// Convert the polygon vertices from the cell frame to the element frame
						for( int i=0 ; i<polygon.size() ; i++ ) polygon[i] = CellToElement( polygon[i] );

						SquareMatrix< double , 2 > element_metric = element_to_cell_differential.transpose() * cell_metric * element_to_cell_differential;
						SquareMatrix< double , 2 > element_metric_inverse = element_metric.inverse();
						double element_area = sqrt( element_metric.determinant() );

						QuadraticElementSample quadraticElementSample( fastIntegration ? 1 : (unsigned int)(polygon.size()-2)*Samples );
						for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) quadraticElementSample.tensor(x,y) = (Real)element_metric_inverse(x,y);
						const QuadraticElementIndex& triangleElementIndices = gridChart.boundaryTriangles[localBoundaryIndex][bt].indices;
						for( int k=0 ; k<6 ; k++ )
						{
							int _fineBoundaryIndex = fineBoundaryIndex[ triangleElementIndices[k] ];
							if( _fineBoundaryIndex!=-1 ) quadraticElementSample.fineNodes[k] = _fineBoundaryIndex;
							else{ printf( "[ERROR] Invalid fine boundary index!\n" ) ; return 0; }
						}

						for( int p=2 ; p<polygon.size() ; p++ )
						{
							Point2D< double > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };
							SquareMatrix< double , 2 > fragment_to_element_differential;
							for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = dm[x][y];
							double fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );
							double fragment_area = element_area * fragment_to_element_area_scale_factor / 2.;

							if( fragment_area>0 )
							{
								for( int s=0 ; s<Samples ; s++ )
								{
									Point2D< Real > pos = polygon[0] + dm[0] * TriangleIntegrator< Samples >::Positions[s][0] + dm[1] * TriangleIntegrator< Samples >::Positions[s][1];
									if( !InUnitTriangle( pos ) ){ printf( "[ERROR] Sample out of unit right triangle! (%f %f)\n", pos[0] , pos[1] ) ; return 0; }
									else
									{
										pos[0] = (Real)std::max< double >( pos[0] , 0 );
										pos[1] = (Real)std::max< double >( pos[1] , 0 );
										double excess = ( pos[0] + pos[1] ) - 1;
										if( excess>0 ) pos[0] -= (Real)excess/2 , pos[1] -= (Real)excess/2;
									}

									typename QuadraticElementSample::SampleData& sampleData = quadraticElementSample[ fastIntegration ? 0 : (p-2)*Samples+s ];
									SetBoundaryDuals< Samples , Real >( sampleData , quadraticElementSample , pos , (Real)( TriangleIntegrator< Samples >::Weights[s] * fragment_area ) );
									sampleData.pos = pos;
								}
							}
							else
							{
								printf( "[WARNING] Element discarded due to zero mass. Triangle %d. Boundary cell %d . Element %d. Sub triangle %d\n" , t , localBoundaryIndex , bt , p-2 );

								printf( "Atlas triangle\n" );
								printf( "%d \n" , 3 );
								for( int v=0 ; v<3 ; v++ ) printf( "%f %f %f\n" , atlasTriangle.vertices[v][0] , atlasTriangle.vertices[v][1] , 0. );
								printf( "%d %d %d\n" , atlasTriangle.atlasVertexIndices[0] , atlasTriangle.atlasVertexIndices[1] , atlasTriangle.atlasVertexIndices[2] );
								printf( "%d %d %d\n" , atlasTriangle.atlasEdgeIndices[0] , atlasTriangle.atlasEdgeIndices[1] , atlasTriangle.atlasEdgeIndices[2] );
								for( int _bt=0 ; _bt<cellBoundaryTriangles.size() ; _bt++ )
								{
									printf( "ELEMENT %d\n" , _bt );

									BoundaryIndexedTriangle _element = cellBoundaryTriangles[_bt];
									AtlasIndexedPolygon _polygon;
									SetAtlasIndexedPolygonFromBoundaryTriangle( _element , _polygon );

									printf( "Boundary triangle\n" );
									printf( "%d \n" , 3 );
									for( int v=0 ; v<3 ; v++ )printf( "%f %f %f\n" , _element.vertices[v][0] , _element.vertices[v][1] , 0. );
									printf( "%d %d %d\n" , _element.atlasVertexIndices[0] , _element.atlasVertexIndices[1] , _element.atlasVertexIndices[2] );
									printf( "%d %d %d\n" , _element.atlasEdgeIndices[0] , _element.atlasEdgeIndices[1] , _element.atlasEdgeIndices[2] );

									ClipPartiallyIndexedPolygonToIndexedTriangle( _polygon , atlasTriangle , _bt==bt );
								}
								return 0;
							}
						}
						if( fastIntegration ) quadraticElementSample[0].pos = Point2D< Real >( PolygonCenter( &polygon[0] , (unsigned int)polygon.size() ) );
#pragma omp critical
						quadraticElementSamples.push_back( quadraticElementSample );
					}
					else if( clippingResult<0 ) return 0;
				}
			}
		}
	}
	return 1;
}

template< unsigned int Samples , typename Real , typename BilinearElementSample , typename QuadraticElementSample >
int InitializeIntegration
(
	const std::vector<std::vector< SquareMatrix< double , 2 > > >& parameterMetric ,
	const std::vector< AtlasChart >& atlasCharts ,
	const std::vector< GridChart >& gridCharts ,
	const std::vector< std::pair< int , int > >& interiorCellLineIndex ,
	const std::vector< int >& fineBoundaryIndex ,
	std::vector<std::vector< BilinearElementSample > >& bilinearElementSamples ,
	std::vector< QuadraticElementSample >& quadraticElementSamples ,
	bool fastIntegration
)
{
#pragma omp parallel for
	for( int i=0 ; i<gridCharts.size() ; i++ )
		if( !InitializeIntegration< Samples , Real >( parameterMetric[i] , atlasCharts[i] , gridCharts[i] , interiorCellLineIndex , fineBoundaryIndex , bilinearElementSamples , quadraticElementSamples , fastIntegration ) )
			fprintf( stderr , "[ERROR] Failed to intialize integration: %d\n" , i ) , exit( 0 );
	return 1;
}

///////////////////////////
// Compute the integrals //
///////////////////////////
template< class Real , typename T , typename ValueFunctionType >
void IntegrateBilinear
(
	const BilinearElementScalarSample< Real >& sample ,
	const ValueFunctionType& ValueFunction ,
	const T cornerValues[] ,
	T rhsValues[]
)
{
	for( int s=0 ; s<(int)sample.size() ; s++ )
	{
		T scalar = ValueFunction( BilinearValue( cornerValues , sample[s].pos ) , sample.tensor );
		for( int k=0 ; k<4 ; k++ ) rhsValues[k] += scalar * sample[s].dualValues[k];
	}
}
template< class Real , typename T , typename GradientFunctionType >
void IntegrateBilinear
(
	const BilinearElementGradientSample< Real >& sample ,
	const GradientFunctionType& GradientFunction ,
	const T cornerValues[] ,
	T rhsValues[]
)
{
	for( int s=0 ; s<(int)sample.size() ; s++ )
	{
		Point2D< T > gradientVector = GradientFunction( BilinearGradient( cornerValues , sample[s].pos ) , sample.tensor );
		for( int k=0 ; k<4 ; k++ ) for( int d=0 ; d<2 ; d++ ) rhsValues[k] += gradientVector[d] * sample[s].dualGradients[k][d];
	}
}

template< class Real , typename T , typename ValueFunctionType >
void IntegrateQuadratic
(
	const QuadraticElementScalarSample< Real >& sample ,
	const ValueFunctionType& ValueFunction ,
	const T cornerValues[] ,
	T rhsValues[]
)
{
	for( int s=0 ; s<(int)sample.size() ; s++ )
	{
		T scalar = ValueFunction( QuadraticValue( cornerValues , sample[s].pos ) , sample.tensor );
		for( int k=0 ; k<6 ; k++ ) rhsValues[k] += scalar * sample[s].dualValues[k];
	}
}
template< class Real , typename T , typename GradientFunctionType >
void IntegrateQuadratic
(
	const QuadraticElementGradientSample< Real >& sample ,
	const GradientFunctionType& VectorFunction ,
	const T cornerValues[] ,
	T rhsValues[]
)
{
	for( int s=0 ; s<(int)sample.size() ; s++ )
	{
		Point2D< T > gradientVector = VectorFunction( QuadraticGradient( cornerValues , sample[s].pos ) , sample.tensor );
		for( int k=0 ; k<6 ; k++ ) for( int d=0 ; d<2 ; d++ ) rhsValues[k] += gradientVector[d] * sample[s].dualGradients[k][d];
	}
}

template< class Real , typename T , typename BilinearElementSample , typename QuadraticElementSample , typename SampleFunctionType >
int Integrate
(
	const std::vector< InteriorCellLine >& interiorCellLines ,
	const std::vector< std::vector< BilinearElementSample > >& bilinearElementSamples ,
	const std::vector< QuadraticElementSample >& quadraticElementSamples ,
	const std::vector< T >& potential ,
	const std::vector< T >& boundary_potential ,
	const SampleFunctionType& SampleFunction ,
	std::vector< T >& rhs ,
	std::vector< T >& boundary_rhs ,
	bool verbose=false
)
{
	clock_t begin = clock();
	auto UpdateRow = [&]( int r )
	{
		const T* _inPrevious = &potential[ interiorCellLines[r].prevLineIndex ];
		const T*     _inNext = &potential[ interiorCellLines[r].nextLineIndex ];

		T* _outPrevious = &rhs[ interiorCellLines[r].prevLineIndex ];
		T*     _outNext = &rhs[ interiorCellLines[r].nextLineIndex ];

		T cornerValues[4];
		cornerValues[0] = *_inPrevious;
		_inPrevious++;
		cornerValues[1] = *_inPrevious;
		cornerValues[3] = *_inNext;
		_inNext++;
		cornerValues[2] = *_inNext;

		T rhsValues[] = { T() , T() , T() , T() };

		int numSamples = (int)bilinearElementSamples[r].size();
		const BilinearElementSample* samplePtr = &bilinearElementSamples[r][0];
		int currentOffset = 0;

		for( int j=0 ; j<numSamples ; j++ )
		{
			const BilinearElementSample& sample = *samplePtr;
			if( sample.cellOffset>currentOffset )
			{
				cornerValues[0] = cornerValues[1];
				_inPrevious++;
				cornerValues[1] = *_inPrevious;
				cornerValues[3] = cornerValues[2];
				_inNext++;
				cornerValues[2] = *_inNext;

				*_outPrevious += rhsValues[0];
				_outPrevious++;
				*_outNext += rhsValues[3];
				_outNext++;
				rhsValues[0] = rhsValues[1];
				rhsValues[1] = T();
				rhsValues[3] = rhsValues[2];
				rhsValues[2] = T();

				currentOffset++;
			}

			IntegrateBilinear< Real , T >( sample , SampleFunction , cornerValues , rhsValues );
			samplePtr++;
		}

		*_outPrevious += rhsValues[0];
		_outPrevious++;
		*_outPrevious += rhsValues[1];
		*_outNext += rhsValues[3];
		_outNext++;
		*_outNext += rhsValues[2];
	};

	int threads = omp_get_max_threads();
	std::vector< int > lineRange( threads+1 );
	int blockSize = (int)interiorCellLines.size() / threads;
	for( int t=0 ; t<threads ; t++ ) lineRange[t] = t*blockSize;
	lineRange[threads] = (int)interiorCellLines.size();
#pragma omp parallel for
	for( int t=0 ; t<threads ; t++ )
	{
		const int tId = omp_get_thread_num();
		const int firstLine = lineRange[tId];
		const int lastLine = lineRange[tId+1];
		for( int r=firstLine ; r<lastLine ; r++ ) UpdateRow(r);
	}

	if( verbose ) printf( "Integrated bilinear: %.2f(s)\n" , double( clock()-begin ) / CLOCKS_PER_SEC );

	begin = clock();
	for( int i=0 ; i<quadraticElementSamples.size() ; i++ )
	{
		const QuadraticElementSample& sample = quadraticElementSamples[i];
		// The values of the potential at the vertices and edge mid-points
		T cornerValues[] = { boundary_potential[ sample.fineNodes[0] ] , boundary_potential[ sample.fineNodes[1] ] , boundary_potential[ sample.fineNodes[2] ] , boundary_potential[ sample.fineNodes[3] ] , boundary_potential[ sample.fineNodes[4] ] , boundary_potential[ sample.fineNodes[5] ] };
		T rhsValues[] = { T() , T() , T() , T() , T() , T() };
		IntegrateQuadratic< Real , T >( sample , SampleFunction , cornerValues , rhsValues );
		for( int k=0 ; k<6 ; k++ ) boundary_rhs[ sample.fineNodes[k] ] += rhsValues[k];
	}
	if( verbose ) printf( "Integrated quadratic: %.2f(s)\n" , double( clock()-begin ) / CLOCKS_PER_SEC );
	return 1;
}
