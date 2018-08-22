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

#ifdef MISHA_CODE
// Note that integration is implemented using a hybrid-single-point quadrature:
// -- In a pre-processing step we compute the integrated gradients of the four bilinear basis functions
// -- At run time we take a single dot-product with the gradient of the function sampled at the barycenter
// In what follows we have the following mappings:
// fragment -> element -> texture:
//		fragment: a triangle in the decomposition of the clipping of the atlas triangle to the elements
//		element:  either a square or a triangle (depending on whether the cell is interior or not)
//		texture: the texture space grid
template< unsigned int Samples , typename Real >
int InitializeVectorFieldIntegration
(
	const std::vector< SquareMatrix< double , 2 > >& texture_metrics ,
	const AtlasChart& atlasChart ,
	const GridChart& gridChart ,
	const std::vector< std::pair< int , int > >& interiorCellLineIndex ,
	const std::vector< int >& fineBoundaryIndex ,
	std::vector< std::vector< BilinearElementGradientSample< Real > > >& bilinearElementGradientSamples ,
	std::vector< QuadraticElementGradientSample< Real > >& quadraticElementGradientSamples ,
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


	typename BilinearElementGradientSample< Real >::SampleData interior_cell_samples[2*Samples] , interior_cell_sample;
	{
		std::vector< Point2D< double > > polygon = { Point2D< double >(0,0) , Point2D< double >(1,0) , Point2D< double >(1,1) , Point2D< double >(0,1) };

		interior_cell_sample.pos = Point2D< Real >( (Real)0.5 , (Real)0.5 );
		for( int p=2 ; p<polygon.size() ; p++ )
		{
			Point2D< double > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };
			SquareMatrix< double , 2 > fragment_to_element_differential;
			for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = dm[x][y];
			double fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );

			for( int s=0 ; s<Samples ; s++ )
			{
				typename BilinearElementGradientSample< Real >::SampleData& sample = interior_cell_samples[(p-2)*Samples+s];
				sample.pos = polygon[0] + dm[0] * TriangleIntegrator< Samples >::Positions[s][0] + dm[1] * TriangleIntegrator< Samples >::Positions[s][1];

				// Compute the integrated gradient of the bilinear basis function in the frame of the cell, weighted by the area of the unit triangle
				for( int k=0 ; k<4 ; k++ ) sample.dualGradients[k] = BilinearElementGradient( k , sample.pos ) * TriangleIntegrator< Samples >::Weights[s] * fragment_to_element_area_scale_factor / 2.0;
				// Compute the integrated gradient of the bilinear basis function in the frame of the cell, weighted by the area of the unit triangle
				for( int k=0 ; k<4 ; k++ ) interior_cell_sample.dualGradients[k] += sample.dualGradients[k];
			}
		}
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
		double cell_area_scale_factor = sqrt( cell_metric.determinant() );

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
				double element_area_scale_factor = cell_area_scale_factor;

				int globalInteriorIndex = localInteriorIndex + gridChart.globalIndexInteriorCellOffset;
				int cellLineId = interiorCellLineIndex[globalInteriorIndex].first;
				int cellLineOffset = interiorCellLineIndex[globalInteriorIndex].second;

				BilinearElementGradientSample< Real > bilinearElementGradientSample( fastIntegration ? 1 : 2*Samples );
				bilinearElementGradientSample.cellOffset = cellLineOffset;
				for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) bilinearElementGradientSample.tensor(x,y) = (Real)element_metric_inverse(x,y);

				if( fastIntegration )
				{
					{
						bilinearElementGradientSample[0].pos = interior_cell_sample.pos;
						for( int k=0 ; k<4 ; k++ ) bilinearElementGradientSample[0].dualGradients[k] = bilinearElementGradientSample.tensor * interior_cell_sample.dualGradients[k] * (Real)element_area_scale_factor;
					}
				}
				else
				{
					for( int s=0 ; s<2*Samples ; s++ )
					{
						bilinearElementGradientSample[s].pos = interior_cell_samples[ fastIntegration ? 0 : s ].pos;
						for( int k=0 ; k<4 ; k++ ) bilinearElementGradientSample[s].dualGradients[k] = bilinearElementGradientSample.tensor * interior_cell_samples[s].dualGradients[k] * (Real)element_area_scale_factor;
					}
				}
#pragma omp critical
				bilinearElementGradientSamples[ cellLineId ].push_back( bilinearElementGradientSample );
			}
			else if( localInteriorIndex!=-1 )
			{
				// For interior cells, the cell and the element are the same thing
				auto TextureToElement = TextureToCell;
				SquareMatrix< double , 2 > element_metric = cell_metric , element_metric_inverse = cell_metric_inverse;
				double element_area_scale_factor = cell_area_scale_factor;

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
					BilinearElementGradientSample< Real > bilinearElementGradientSample( fastIntegration ? 1 : (polygon.size()-2)*Samples );
					bilinearElementGradientSample.cellOffset = cellLineOffset;
					for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) bilinearElementGradientSample.tensor(x,y) = (Real)element_metric_inverse(x,y);

					for( int p=2 ; p<polygon.size() ; p++ )
					{
						Point2D< double > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };
						SquareMatrix< double , 2 > fragment_to_element_differential;
						for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = dm[x][y];
						double fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );

						for( int s=0 ; s<Samples ; s++ )
						{
							typename BilinearElementGradientSample< Real >::SampleData& sample = bilinearElementGradientSample[ fastIntegration ? 0 : (p-2)*Samples+s ];
							sample.pos = polygon[0] + dm[0] * TriangleIntegrator< Samples >::Positions[s][0] + dm[1] * TriangleIntegrator< Samples >::Positions[s][1];
							if( !InUnitSquare( sample.pos ) ){ printf( "[ERROR] Sample position out of unit square! (%f %f)\n" , sample.pos[0] , sample.pos[1] ) ; return 0; }

							// Compute the integrated gradients of each bilinear basis functions in the frame of the cell, weighted by the area of the unit triangle
							// Dualize the samples so that integration is simple a dot-product
							for( int k=0 ; k<4 ; k++ ) sample.dualGradients[k] += bilinearElementGradientSample.tensor * Point2D< Real >( BilinearElementGradient( k , sample.pos ) * (Real)( TriangleIntegrator< Samples >::Weights[s] * element_area_scale_factor * fragment_to_element_area_scale_factor / 2.0 ) );
						}
					}
					if( fastIntegration ) bilinearElementGradientSample[0].pos = Point2D< Real >( PolygonCenter( &polygon[0] , (unsigned int)polygon.size() ) );
#pragma omp critical
					bilinearElementGradientSamples[cellLineId].push_back( bilinearElementGradientSample );
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
						double element_area_scale_factor = sqrt( element_metric.determinant() );

						QuadraticElementGradientSample< Real > quadraticElementGradientSample( fastIntegration ? 1 : (unsigned int)(polygon.size()-2)*Samples );
						for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) quadraticElementGradientSample.tensor(x,y) = (Real)element_metric_inverse(x,y);
						const QuadraticElementIndex& triangleElementIndices = gridChart.boundaryTriangles[localBoundaryIndex][bt].indices;
						for( int k=0 ; k<6 ; k++ )
						{
							int _fineBoundaryIndex = fineBoundaryIndex[ triangleElementIndices[k] ];
							if( _fineBoundaryIndex!=-1 ) quadraticElementGradientSample.fineNodes[k] = _fineBoundaryIndex;
							else{ printf( "[ERROR] Invalid fine boundary index!\n" ) ; return 0; }
						}

						for( int p=2 ; p<polygon.size() ; p++ )
						{
							Point2D< double > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };
							SquareMatrix< double , 2 > fragment_to_element_differential;
							for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = dm[x][y];
							double fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );

							double fragment_area = element_area_scale_factor * fragment_to_element_area_scale_factor / 2.;
							if( fragment_area>0 )
							{
								for( int s=0 ; s<Samples ; s++ )
								{
									typename QuadraticElementGradientSample< Real >::SampleData& sample = quadraticElementGradientSample[ fastIntegration ? 0 : (p-2)*Samples+s ];
									sample.pos = polygon[0] + dm[0] * TriangleIntegrator< Samples >::Positions[s][0] + dm[1] * TriangleIntegrator< Samples >::Positions[s][1];
									if( !InUnitTriangle( sample.pos ) ){ printf( "[ERROR] Sample out of unit right triangle! (%f %f)\n", sample.pos[0] , sample.pos[1] ) ; return 0; }
									else
									{
										sample.pos[0] = (Real)std::max< double >( sample.pos[0] , 0 );
										sample.pos[1] = (Real)std::max< double >( sample.pos[1] , 0 );
										double excess = ( sample.pos[0] + sample.pos[1] ) - 1;
										if( excess>0 ) sample.pos[0] -= (Real)excess/2 , sample.pos[1] -= (Real)excess/2;
									}

									// Dualize the samples so that integration is simple a dot-product
									for( int k=0 ; k<6 ; k++ ) sample.dualGradients[k] += quadraticElementGradientSample.tensor * Point2D< Real >( QuadraticElementGradient( k , sample.pos ) * (Real)( TriangleIntegrator< Samples >::Weights[s] * fragment_area ) );
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
						if( fastIntegration ) quadraticElementGradientSample[0].pos = Point2D< Real >( PolygonCenter( &polygon[0] , (unsigned int)polygon.size() ) );
#pragma omp critical
						quadraticElementGradientSamples.push_back( quadraticElementGradientSample );
					}
					else if( clippingResult<0 ) return 0;
				}
			}
		}
	}
	return 1;
}

template< unsigned int Samples , typename Real >
int InitializeVectorFieldIntegration
(
	const std::vector<std::vector< SquareMatrix< double , 2 > > >& parameterMetric ,
	const std::vector< AtlasChart >& atlasCharts ,
	const std::vector< GridChart >& gridCharts ,
	const std::vector< std::pair< int , int > >& interiorCellLineIndex ,
	const std::vector< int >& fineBoundaryIndex ,
	std::vector<std::vector< BilinearElementGradientSample< Real > > >& bilinearElementGradientSamples ,
	std::vector< QuadraticElementGradientSample< Real > >& quadraticElementGradientSamples ,
	bool fastIntegration
)
{
#pragma omp parallel for
	for( int i=0 ; i<gridCharts.size() ; i++ )
		if( !InitializeVectorFieldIntegration< Samples >( parameterMetric[i] , atlasCharts[i] , gridCharts[i] , interiorCellLineIndex , fineBoundaryIndex , bilinearElementGradientSamples , quadraticElementGradientSamples , fastIntegration ) )
			fprintf( stderr , "[ERROR] Failed to intialize vector field integration: %d\n" , i ) , exit( 0 );
	return 1;
}

template< class Real >
int IntegrateVectorField
(
	const std::vector< InteriorCellLine >& interiorCellLines ,
	const std::vector<std::vector< BilinearElementGradientSample< Real > > >& bilinearElementGradientSamples ,
	const std::vector< QuadraticElementGradientSample< Real > >& quadraticElementGradientSamples ,
	const std::vector< Real >& potential ,
	const std::function< Point2D< Real > ( Point2D< Real > , SquareMatrix< Real , 2 >  ) >& VectorFunction ,
	std::vector< Real > & rhs ,
	const std::vector< Real >& boundary_potential ,
	std::vector< Real >& boundary_rhs ,
	bool verbose=false
)
{
	memset( &rhs[0] , 0 , rhs.size() * sizeof(Real) );
	memset( &boundary_rhs[0] , 0 , boundary_rhs.size() * sizeof(Real) );

	clock_t begin = clock();
	auto UpdateRow = [&]( int r )
	{
		const Real* _inPrevious = &potential[ interiorCellLines[r].prevLineIndex ];
		const Real*     _inNext = &potential[ interiorCellLines[r].nextLineIndex ];

		Real* _outPrevious = &rhs[ interiorCellLines[r].prevLineIndex ];
		Real*     _outNext = &rhs[ interiorCellLines[r].nextLineIndex ];

		Real cornerValues[4];
		cornerValues[0] = *_inPrevious;
		_inPrevious++;
		cornerValues[1] = *_inPrevious;
		cornerValues[3] = *_inNext;
		_inNext++;
		cornerValues[2] = *_inNext;

		Real rhsValues[4] = { 0 , 0 , 0 , 0 };

		int numSamples = (int)bilinearElementGradientSamples[r].size();
		const BilinearElementGradientSample< Real >* samplePtr = &bilinearElementGradientSamples[r][0];
		int currentOffset = 0;

		for( int j=0 ; j<numSamples ; j++ )
		{
			const BilinearElementGradientSample< Real >& sample = *samplePtr;
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
				rhsValues[1] = 0;
				rhsValues[3] = rhsValues[2];
				rhsValues[2] = 0;

				currentOffset++;
			}

			for( int s=0 ; s<(int)sample.size() ; s++ )
			{
				Point2D< Real > gradientVector = VectorFunction( BilinearGradient( cornerValues , sample[s].pos ) , sample.tensor );
				for( int k=0 ; k<4 ; k++ ) rhsValues[k] += Point2D< Real >::Dot( gradientVector , sample[s].dualGradients[k] );
			}
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

	if( verbose ) printf( "Integrated bilinear %.4f\n" , double( clock()-begin ) / CLOCKS_PER_SEC );

	begin = clock();
	for( int i=0 ; i<quadraticElementGradientSamples.size() ; i++ )
	{
		const QuadraticElementGradientSample< Real >& sample = quadraticElementGradientSamples[i];
		// The values of the potential at the vertices and edge mid-points
		Real cornerValues[] = { boundary_potential[ sample.fineNodes[0] ] , boundary_potential[ sample.fineNodes[1] ] , boundary_potential[ sample.fineNodes[2] ] , boundary_potential[ sample.fineNodes[3] ] , boundary_potential[ sample.fineNodes[4] ] , boundary_potential[ sample.fineNodes[5] ] };
		for( int s=0 ; s<(int)sample.size() ; s++ )
		{
			Point2D< Real > gradientVector = VectorFunction( QuadraticGradient( cornerValues , sample[s].pos ) , sample.tensor );
			for( int k=0 ; k<6 ; k++ ) boundary_rhs[ sample.fineNodes[k] ] += Point2D< Real >::Dot( gradientVector , sample[s].dualGradients[k] );
		}
	}
	if( verbose ) printf( "Integrated quadratic %.4f\n" , double( clock()-begin ) / CLOCKS_PER_SEC);
	return 1;
}
#else // !MISHA_CODE
// Note that integration is implemented using a hybrid-single-point quadrature:
// -- In a pre-processing step we compute the integrated gradients of the four bilinear basis functions
// -- At run time we take a single dot-product with the gradient of the function sampled at the barycenter
// In what follows we have the following mappings:
// fragment -> element -> texture:
//		fragment: a triangle in the decomposition of the clipping of the atlas triangle to the elements
//		element:  either a square or a triangle (depending on whether the cell is interior or not)
//		texture: the texture space grid
template< unsigned int Samples , typename Real >
int InitializeVectorFieldIntegration
(
	const std::vector< SquareMatrix< double , 2 > >& texture_metrics ,
	const AtlasChart& atlasChart ,
	const GridChart& gridChart ,
	const std::vector< std::pair< int , int > >& interiorCellLineIndex ,
	const std::vector< int >& fineBoundaryIndex ,
	std::vector< std::vector< BilinearElementGradientSample< Real > > >& bilinearElementGradientSamples ,
	std::vector< QuadraticElementGradientSample< Real > >& quadraticElementGradientSamples
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

	Point2D< double > interior_cell_barycenter(0.5,0.5) , interior_cell_v[4];
	{
		std::vector< Point2D< double > > polygon = { Point2D< double >(0,0) , Point2D< double >(1,0) , Point2D< double >(1,1) , Point2D< double >(0,1) };

		for( int p=2 ; p<polygon.size() ; p++ )
		{
			Point2D< double > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };
			SquareMatrix< double , 2 > fragment_to_element_differential;
			for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = dm[x][y];
			double fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );

			Point2D< double > fragment_samples[Samples];
			for( int s=0 ; s<Samples ; s++ ) fragment_samples[s] = polygon[0] + dm[0] * TriangleIntegrator<Samples>::Positions[s][0] + dm[1] * TriangleIntegrator<Samples>::Positions[s][1];

			// Compute the integrated gradients of each bilinear basis functions
			for( int k=0 ; k<4 ; k++ )
			{
				Point2D< double > integrated_gradient;
				// Compute the integrated gradient of the bilinear basis function in the frame of the cell, weighted by the area of the unit triangle
				for( int s=0 ; s<Samples ; s++ ) integrated_gradient += BilinearElementGradient( k , fragment_samples[s] ) * TriangleIntegrator<Samples>::Weights[s];
				interior_cell_v[k] += integrated_gradient * fragment_to_element_area_scale_factor / 2.0;
			}
		}
	}

	for( int t=0 ; t<atlasChart.triangles.size() ; t++ )
	{
		Point2D< double > tPos[3];
		for( int i=0 ; i<3 ; i++ ) tPos[i] = atlasChart.vertices[ atlasChart.triangles[t][i] ] - gridChart.corner;

		SquareMatrix< double , 2 > texture_metric = texture_metrics[t];
		SquareMatrix< double , 2 > cell_metric = cell_to_texture_differential.transpose() * texture_metric * cell_to_texture_differential;
		SquareMatrix< double , 2 > cell_metric_inverse = cell_metric.inverse();
		double cell_area_scale_factor = sqrt( cell_metric.determinant() );

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
				double element_area_scale_factor = cell_area_scale_factor;

				int globalInteriorIndex = localInteriorIndex + gridChart.globalIndexInteriorCellOffset;
				int cellLineId = interiorCellLineIndex[globalInteriorIndex].first;
				int cellLineOffset = interiorCellLineIndex[globalInteriorIndex].second;

				BilinearElementGradientSample< Real > bilinearElementGradientSample;
				bilinearElementGradientSample.cellOffset = cellLineOffset;

				bilinearElementGradientSample.integrationData.pos = Point2D< Real >( interior_cell_barycenter );
				for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) bilinearElementGradientSample.tensor(x,y) = (Real)element_metric_inverse(x,y);
				for( int k=0 ; k<4 ; k++ ) bilinearElementGradientSample.integrationData.dualGradients[k] = bilinearElementGradientSample.tensor * Point2D< Real >( interior_cell_v[k] * element_area_scale_factor );
#pragma omp critical
				bilinearElementGradientSamples[ cellLineId ].push_back( bilinearElementGradientSample );
			}
			else if( localInteriorIndex!=-1 )
			{
				// For interior cells, the cell and the element are the same thing
				auto TextureToElement = TextureToCell;
				SquareMatrix< double , 2 > element_metric = cell_metric , element_metric_inverse = cell_metric_inverse;
				double element_area_scale_factor = cell_area_scale_factor;

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
					BilinearElementGradientSample< Real > bilinearElementGradientSample;
					bilinearElementGradientSample.cellOffset = cellLineOffset;

					// Compute the center of the polygon
					Point2D< double > polygonBarycenter;
					{
						double polygonArea = 0;
						for( int p=2 ; p<polygon.size() ; p++ )
						{
							Point2D< double > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };
							double fragment_area = fabs( dm[0][0] * dm[1][1] - dm[0][1] * dm[1][0] ) / 2.0;
							polygonBarycenter += fragment_area * (polygon[0] + polygon[p-1] + polygon[p]) / 3.0;
							polygonArea += fragment_area;
						}
						polygonBarycenter /= polygonArea;
					}

					if( !InUnitSquare( polygonBarycenter ) ){ printf( "[ERROR] Center out of unit square! (%f %f)\n" , polygonBarycenter[0] , polygonBarycenter[1] ) ; return 0; }

					bilinearElementGradientSample.integrationData.pos = Point2D< Real >( (Real)polygonBarycenter[0] , (Real)polygonBarycenter[1] );
					for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) bilinearElementGradientSample.tensor(x,y) = (Real)element_metric_inverse(x,y);

					for( int p=2 ; p<polygon.size() ; p++ )
					{
						Point2D< double > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };
						SquareMatrix< double , 2 > fragment_to_element_differential;
						for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = dm[x][y];
						double fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );

						Point2D< double > fragment_samples[Samples];
						for( int s=0 ; s<Samples ; s++ )
						{
							fragment_samples[s] = polygon[0] + dm[0] * TriangleIntegrator<Samples>::Positions[s][0] + dm[1] * TriangleIntegrator<Samples>::Positions[s][1];
							if( !InUnitSquare( fragment_samples[s] ) ){ printf( "[ERROR] Sample out of unit square! (%f %f)\n" , fragment_samples[s][0] , fragment_samples[s][1]) ; return 0; }
						}

						// Compute the integrated gradients of each bilinear basis functions
						for( int k=0 ; k<4 ; k++ )
						{
							Point2D< double > integrated_gradient;
							// Compute the integrated gradient of the bilinear basis function in the frame of the cell, weighted by the area of the unit triangle
							for( int s=0 ; s<Samples ; s++ ) integrated_gradient += BilinearElementGradient( k , fragment_samples[s] ) * TriangleIntegrator<Samples>::Weights[s];
							// [MK] Why the division by 2?
							bilinearElementGradientSample.integrationData.dualGradients[k] += Point2D< Real >( integrated_gradient * element_area_scale_factor * fragment_to_element_area_scale_factor / 2.0 );
						}
					}
					// Dualize the samples so that integration is simple a dot-product
					for( int k=0 ; k<4 ; k++ ) bilinearElementGradientSample.integrationData.dualGradients[k] = bilinearElementGradientSample.tensor * bilinearElementGradientSample.integrationData.dualGradients[k];
#pragma omp critical
					bilinearElementGradientSamples[cellLineId].push_back( bilinearElementGradientSample );
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

						QuadraticElementGradientSample< Real > quadraticElementGradientSample;

						const QuadraticElementIndex& triangleElementIndices = gridChart.boundaryTriangles[localBoundaryIndex][bt].indices;
						for( int k=0 ; k<6 ; k++ )
						{
							int _fineBoundaryIndex = fineBoundaryIndex[ triangleElementIndices[k] ];
							if( _fineBoundaryIndex!=-1 ) quadraticElementGradientSample.fineNodes[k] = _fineBoundaryIndex;
							else{ printf( "[ERROR] Invalid fine boundary index!\n" ) ; return 0; }
						}

						SquareMatrix< double , 2 > element_to_cell_differential , cell_to_element_differential;
						Point2D< double > dm[] = { element_vertices[1]-element_vertices[0] , element_vertices[2]-element_vertices[0] };
						for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) element_to_cell_differential(x,y) = dm[x][y];
						cell_to_element_differential = element_to_cell_differential.inverse();
						auto CellToElement = [&]( Point2D< double > p ){ return cell_to_element_differential * ( p - element_vertices[0] ); };
						// Convert the polygon vertices from the cell frame to the element frame
						for( int i=0 ; i<polygon.size() ; i++ ) polygon[i] = CellToElement( polygon[i] );

						SquareMatrix< double , 2 > element_metric = element_to_cell_differential.transpose() * cell_metric * element_to_cell_differential;
						SquareMatrix< double , 2 > element_metric_inverse = element_metric.inverse();
						double element_area_scale_factor = sqrt( element_metric.determinant() );

						Point2D< double > polygonBarycenter;
						{
							double polygonArea = 0;
							for( int p=2 ; p<polygon.size() ; p++ )
							{
								Point2D< double > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };
								double fragment_area = fabs(dm[0][0] * dm[1][1] - dm[0][1] * dm[1][0]) / 2.0;
								polygonBarycenter += fragment_area * ( polygon[0] + polygon[p-1] + polygon[p] ) / 3.0;
								polygonArea += fragment_area;
							}
							polygonBarycenter /= polygonArea;
						}

						if( !InUnitTriangle( polygonBarycenter ) ){ printf( "[ERROR] Center out of unit right triangle! (%f %f)\n", polygonBarycenter[0] , polygonBarycenter[1] ) ; return 0; }

						quadraticElementGradientSample.integrationData.pos = Point2D< Real >( (Real)polygonBarycenter[0] , (Real)polygonBarycenter[1] );
						for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) quadraticElementGradientSample.tensor(x,y) = (Real)element_metric_inverse(x,y);
						// Iterate over all fragments
						for( int p=2 ; p<polygon.size() ; p++ )
						{
							Point2D< double > dm[] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };

							SquareMatrix< double , 2 > fragment_to_element_differential;
							for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = dm[x][y];
							double fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );

							double fragment_area = element_area_scale_factor * fragment_to_element_area_scale_factor / 2.;

							if( fragment_area>0 )
							{
								Point2D< double > fragment_samples[Samples];
								for( int s=0 ; s<Samples ; s++ )
								{
									fragment_samples[s] = polygon[0] + dm[0] * TriangleIntegrator<Samples>::Positions[s][0] + dm[1] * TriangleIntegrator<Samples>::Positions[s][1];
									if( !InUnitTriangle( fragment_samples[s] ) ){ printf( "[ERROR] Sample out of unit right triangle! (%f %f)\n" , fragment_samples[s][0] , fragment_samples[s][1] ) ; return 0; }
									else
									{
										fragment_samples[s][0] = std::max< double >( fragment_samples[s][0] , 0 );
										fragment_samples[s][1] = std::max< double >( fragment_samples[s][1] , 0 );
										double excess = ( fragment_samples[s][0] + fragment_samples[s][1] ) - 1;
										if( excess>0 ) fragment_samples[s][0] -= excess/2 , fragment_samples[s][1] -= excess/2;
									}
								}


								for( int k=0 ; k<6 ; k++ ) // which function
								{
									Point2D< double > integrated_gradient;
									for( int s=0 ; s<Samples ; s++ ) integrated_gradient += QuadraticElementGradient( k , fragment_samples[s] ) * TriangleIntegrator<Samples>::Weights[s];
									// [MK] Why the division by 2?
									quadraticElementGradientSample.integrationData.dualGradients[k] += Point2D< Real >( integrated_gradient * fragment_area );
								}
							}
							else
							{
								printf("WARNING: Element discarded due to zero mass. Triangle %d. Boundary cell %d . Element %d. Sub triangle %d\n", t, localBoundaryIndex, bt, p - 2);

								printf("Atlas triangle \n");
								printf("%d \n", 3);
								for (int v = 0; v < 3; v++)printf("%f %f %f \n", atlasTriangle.vertices[v][0], atlasTriangle.vertices[v][1], 0.);
								printf("%d %d %d \n", atlasTriangle.atlasVertexIndices[0], atlasTriangle.atlasVertexIndices[1], atlasTriangle.atlasVertexIndices[2]);
								printf("%d %d %d \n", atlasTriangle.atlasEdgeIndices[0], atlasTriangle.atlasEdgeIndices[1], atlasTriangle.atlasEdgeIndices[2]);
								for (int _bt = 0; _bt < cellBoundaryTriangles.size(); _bt++) {
									printf("ELEMENT %d \n", _bt);

									BoundaryIndexedTriangle _element = cellBoundaryTriangles[_bt];
									AtlasIndexedPolygon _polygon;
									SetAtlasIndexedPolygonFromBoundaryTriangle( _element , _polygon );

									printf("Boundary triangle \n");
									printf("%d \n", 3);
									for (int v = 0; v < 3; v++)printf("%f %f %f \n", _element.vertices[v][0], _element.vertices[v][1], 0.);
									printf("%d %d %d \n", _element.atlasVertexIndices[0], _element.atlasVertexIndices[1], _element.atlasVertexIndices[2]);
									printf("%d %d %d \n", _element.atlasEdgeIndices[0], _element.atlasEdgeIndices[1], _element.atlasEdgeIndices[2]);

									ClipPartiallyIndexedPolygonToIndexedTriangle( _polygon , atlasTriangle , _bt == bt );
								}
								return 0;
							}
						}
						// Dualize the samples so that integration is simple a dot-product
						for( int k=0 ; k<6 ; k++ ) quadraticElementGradientSample.integrationData.dualGradients[k] = quadraticElementGradientSample.tensor * quadraticElementGradientSample.integrationData.dualGradients[k];
#pragma omp critical
						quadraticElementGradientSamples.push_back( quadraticElementGradientSample );
					}
					else if( clippingResult<0 ) return 0;
				}
			}
		}
	}
	return 1;
}

template< unsigned int Samples , typename Real >
int InitializeVectorFieldIntegration
(
	const std::vector<std::vector< SquareMatrix< double , 2 > > >& parameterMetric ,
	const std::vector< AtlasChart >& atlasCharts ,
	const std::vector< GridChart >& gridCharts ,
	const std::vector< std::pair< int , int > >& interiorCellLineIndex ,
	const std::vector< int >& fineBoundaryIndex ,
	std::vector<std::vector< BilinearElementGradientSample< Real > > >& bilinearElementGradientSamples ,
	std::vector< QuadraticElementGradientSample< Real > >& quadraticElementGradientSamples
)
{
#pragma omp parallel for
	for( int i=0 ; i<gridCharts.size() ; i++ )
		if( !InitializeVectorFieldIntegration< Samples >( parameterMetric[i] , atlasCharts[i] , gridCharts[i] , interiorCellLineIndex , fineBoundaryIndex , bilinearElementGradientSamples , quadraticElementGradientSamples ) )
			fprintf( stderr , "[ERROR] Failed to intialize vector field integration: %d\n" , i ) , exit( 0 );
	return 1;
}

template< class Real >
int IntegrateVectorField
(
	const std::vector< InteriorCellLine >& interiorCellLines ,
	const std::vector<std::vector< BilinearElementGradientSample< Real > > >& bilinearElementGradientSamples ,
	const std::vector< QuadraticElementGradientSample< Real > >& quadraticElementGradientSamples ,
	const std::vector< Real >& potential ,
	const std::function< Point2D< Real > ( Point2D< Real > , SquareMatrix< Real , 2 >  ) >& VectorFunction ,
	std::vector< Real > & rhs ,
	const std::vector< Real >& boundary_potential ,
	std::vector< Real >& boundary_rhs ,
	bool verbose=false
)
{
	memset( &rhs[0] , 0 , rhs.size() * sizeof(Real) );
	memset( &boundary_rhs[0] , 0 , boundary_rhs.size() * sizeof(Real) );

	clock_t begin = clock();
	auto UpdateRow = [&]( int r )
	{
		const Real* _inPrevious = &potential[ interiorCellLines[r].prevLineIndex ];
		const Real*     _inNext = &potential[ interiorCellLines[r].nextLineIndex ];

		Real* _outPrevious = &rhs[ interiorCellLines[r].prevLineIndex ];
		Real*     _outNext = &rhs[ interiorCellLines[r].nextLineIndex ];

		Real cornerValues[4];
		cornerValues[0] = *_inPrevious;
		_inPrevious++;
		cornerValues[1] = *_inPrevious;
		cornerValues[3] = *_inNext;
		_inNext++;
		cornerValues[2] = *_inNext;

		Real rhsValues[4] = { 0 , 0 , 0 , 0 };

		int numSamples = (int)bilinearElementGradientSamples[r].size();
		const BilinearElementGradientSample< Real >* sample = &bilinearElementGradientSamples[r][0];
		int currentOffset = 0;

		for( int j=0 ; j<numSamples ; j++ )
		{
			if( sample->cellOffset>currentOffset )
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
				rhsValues[1] = 0;
				rhsValues[3] = rhsValues[2];
				rhsValues[2] = 0;

				currentOffset++;
			}

			Point2D< Real > gradientVector = VectorFunction( BilinearGradient( cornerValues , sample->integrationData.pos ) , sample->tensor );
			for( int k=0 ; k<4 ; k++ ) rhsValues[k] += Point2D<Real>::Dot( gradientVector , sample->integrationData.dualGradients[k] );
			sample++;
		}

		*_outPrevious += rhsValues[0];
		_outPrevious++;
		*_outPrevious += rhsValues[1];
		*_outNext += rhsValues[3];
		_outNext++;
		*_outNext += rhsValues[2];
	};

	//for (int r = 0; r < interiorCellLines.size(); r++) UpdateRow(r);

	int threads = omp_get_max_threads();
	std::vector<int> lineRange(threads + 1);
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

	if( verbose ) printf( "Integrating squares %.4f\n" , double( clock()-begin ) / CLOCKS_PER_SEC );

	begin = clock();
	for( int i=0 ; i<quadraticElementGradientSamples.size() ; i++ )
	{
		const QuadraticElementGradientSample< Real >& sample = quadraticElementGradientSamples[i];
		// The values of the potential at the vertices and edge mid-points
		Real cornerValues[] = { boundary_potential[ sample.fineNodes[0] ] , boundary_potential[ sample.fineNodes[1] ] , boundary_potential[ sample.fineNodes[2] ] , boundary_potential[ sample.fineNodes[3] ] , boundary_potential[ sample.fineNodes[4] ] , boundary_potential[ sample.fineNodes[5] ] };
		Point2D< Real > gradientVector = VectorFunction( QuadraticGradient( cornerValues , sample.integrationData.pos ) , sample.tensor );
		for( int k=0 ; k<6 ; k++ ) boundary_rhs[ sample.fineNodes[k] ] += Point2D< Real >::Dot( gradientVector , sample.integrationData.dualGradients[k] );
	}
	if( verbose ) printf( "Integrating triangles %.4f\n" , double( clock()-begin ) / CLOCKS_PER_SEC);
	return 1;
}
#endif // MISHA_CODE