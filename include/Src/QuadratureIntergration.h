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

namespace MishaK
{
	class InteriorCellLine
	{
	public:
		int prevLineIndex , nextLineIndex;
		int length;
	};

	template< typename GeometryReal >
	void InitializeGridChartInteriorCellLines
	(
		const GridChart< GeometryReal > &gridChart ,
		std::vector< InteriorCellLine > &interiorCellLines ,
		std::vector< std::pair< int , int > > &interiorCellLineIndex
	)
	{
		const Image< int >& cellType = gridChart.cellType;
		unsigned int width = cellType.res(0) , height = cellType.res(1);

		const Image< int > & globalTexelIndex = gridChart.globalTexelIndex;
		const Image< int > & globalTexelInteriorIndex = gridChart.globalInteriorTexelIndex;

		int localInteriorCellIndex = 0;

		for( unsigned int j=0 ; j<height ; j++ )
		{
			unsigned int offset = 0;
			bool previousIsInterior = false;
			unsigned int rasterStart = -1;
			while( offset<width )
			{
				bool currentIsInterior = cellType(offset, j) == 1;
				if( ( offset==0 || offset==width-1 ) && currentIsInterior ) MK_THROW( "Unexpected interior cell" );
				if(  currentIsInterior && !previousIsInterior ) rasterStart = offset; //Start raster line
				if( !currentIsInterior &&  previousIsInterior ) //Terminate raster line
				{ 
					InteriorCellLine newLine;
					newLine.prevLineIndex = globalTexelIndex( rasterStart , j+0 );
					newLine.nextLineIndex = globalTexelIndex( rasterStart , j+1 );
					newLine.length = offset - rasterStart;

					if( newLine.prevLineIndex==-1 || newLine.nextLineIndex==-1 ) MK_THROW( "Invalid indexing" );
					int currentLine = (int)interiorCellLines.size();

					for( unsigned int k=0 ; k<offset-rasterStart ; k++ )
					{
						if( (int)gridChart.interiorCellCorners[localInteriorCellIndex][0]!=globalTexelInteriorIndex( rasterStart+k , j ) ) MK_THROW( "Unexpected corner ID" );

						interiorCellLineIndex.push_back( std::pair< int , int >( currentLine , k ) );
						localInteriorCellIndex++;
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
		const std::vector< GridChart< GeometryReal > >& gridCharts ,
		std::vector< InteriorCellLine >& interiorCellLines ,
		std::vector< std::pair< int , int > >& interiorCellLineIndex
	)
	{
		for( int i=0 ; i<gridCharts.size() ; i++ ) InitializeGridChartInteriorCellLines( gridCharts[i] , interiorCellLines , interiorCellLineIndex );
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
		const BilinearElementScalarSample< Real > & ,
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
		const QuadraticElementScalarSample< Real > & ,
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

	template< unsigned int Samples , typename GeometryReal , typename ElementSamples >
	void InitializeIntegration
	(
		const std::vector< SquareMatrix< GeometryReal , 2 > > &texture_metrics ,
		const AtlasChart< GeometryReal > &atlasChart ,
		const GridChart< GeometryReal > &gridChart ,
		const std::vector< std::pair< int , int > >& interiorCellLineIndex ,
		const std::vector< int >& fineBoundaryIndex ,
		ElementSamples &elementSamples ,
		std::mutex &element_samples_bilinear_mutex ,
		std::mutex &element_samples_quadratic_mutex ,
		bool fastIntegration
	)
	{
		////Rasterize
		const GeometryReal PRECISION_ERROR = 1e-3;

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

				for( int s=0 ; s<Samples ; s++ )
				{
					Point2D< GeometryReal > pos = polygon[0] + dm[0] * (GeometryReal)TriangleIntegrator< Samples >::Positions[s][0] + dm[1] * (GeometryReal)TriangleIntegrator< Samples >::Positions[s][1];

					typename ElementSamples::Bilinear::SampleData& _interior_cell_sample = interior_cell_samples[ (p-2)*Samples + s ];
					SetInteriorCellDuals< Samples >( _interior_cell_sample , Point2D< typename ElementSamples::Real >( pos ) , (typename ElementSamples::Real)( TriangleIntegrator< Samples >::Weights[s] * fragment_to_element_area_scale_factor / 2 ) );
					SetInteriorCellDuals< Samples >(  interior_cell_sample , Point2D< typename ElementSamples::Real >( pos ) , (typename ElementSamples::Real)( TriangleIntegrator< Samples >::Weights[s] * fragment_to_element_area_scale_factor / 2 ) );
					_interior_cell_sample.pos = pos;
				}
			}
			interior_cell_sample.pos = Point2D< GeometryReal >( (GeometryReal)0.5 , (GeometryReal)0.5 );
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

		for( int t=0 ; t<atlasChart.triangles.size() ; t++ )
		{
			Point2D< GeometryReal > tPos[3];
			for( int i=0 ; i<3 ; i++ ) tPos[i] = atlasChart.vertices[ atlasChart.triangles[t][i] ] - gridChart.corner;

			SquareMatrix< GeometryReal , 2 > texture_metric = texture_metrics[t];
			SquareMatrix< GeometryReal , 2 > cell_metric = cell_to_texture_differential.transpose() * texture_metric * cell_to_texture_differential;
			SquareMatrix< GeometryReal , 2 > cell_metric_inverse = cell_metric.inverse();
			GeometryReal cell_area = (GeometryReal)sqrt( cell_metric.determinant() );

			//BBox
			int minCorner[2] , maxCorner[2];
			GetTriangleIntegerBBox( tPos , (GeometryReal)1./gridChart.cellSizeW , (GeometryReal)1./gridChart.cellSizeH , minCorner , maxCorner );

			std::vector< Point2D< GeometryReal > > parametricVertices(3);
			parametricVertices[0] = tPos[0], parametricVertices[1] = tPos[1], parametricVertices[2] = tPos[2];

			AtlasIndexedTriangle< GeometryReal > atlasTriangle;
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
				auto TextureToCell = [&]( Point2D< GeometryReal > p ){ return Point2D< GeometryReal >( ( p[0] / gridChart.cellSizeW ) - i , ( p[1] / gridChart.cellSizeH ) - j ); };

				int localInteriorIndex = gridChart.localInteriorCellIndex(i,j) , localBoundaryIndex = gridChart.localBoundaryCellIndex(i,j);
				if( localInteriorIndex!=-1 && localBoundaryIndex!=-1 ) MK_THROW( "Cell simultaneosly interior and boundary" );

				// Interior cells
				// If the cell is entirely interior to the triangle
				if( CellInTriangle( i , j , parametricVertices ) && localInteriorIndex!=-1 )
				{
					// For interior cells, the cell and the element are the same thing
					auto TextureToElement = TextureToCell;
					SquareMatrix< GeometryReal , 2 > element_metric = cell_metric , element_metric_inverse = cell_metric_inverse;
					GeometryReal element_area = cell_area;

					int globalInteriorIndex = localInteriorIndex + gridChart.globalIndexInteriorCellOffset;
					int cellLineId = interiorCellLineIndex[globalInteriorIndex].first;
					int cellLineOffset = interiorCellLineIndex[globalInteriorIndex].second;

					typename ElementSamples::Bilinear bilinearElementSample( fastIntegration ? 1 : 2*Samples );
					bilinearElementSample.cellOffset = cellLineOffset;
					for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) bilinearElementSample.tensor(x,y) = (typename ElementSamples::Real)element_metric_inverse(x,y);

					if( fastIntegration )
					{
						{
							bilinearElementSample[0].pos = interior_cell_sample.pos;
							SetCellInTriangleDuals< Samples >( bilinearElementSample[0] , bilinearElementSample , interior_cell_sample , (typename ElementSamples::Real)element_area/2 );
						}
					}
					else
					{
						for( int s=0 ; s<2*Samples ; s++ )
						{
							bilinearElementSample[s].pos = interior_cell_samples[s].pos;
							SetCellInTriangleDuals< Samples >( bilinearElementSample[s] , bilinearElementSample , interior_cell_samples[s] , (typename ElementSamples::Real)element_area/2 );
						}
					}
					{
						std::lock_guard< std::mutex > lock( element_samples_bilinear_mutex );
						elementSamples.bilinear[ cellLineId ].push_back( bilinearElementSample );
					}
				}
				else if( localInteriorIndex!=-1 )
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
						for( int ii=0 ; ii<polygon.size() ; ii++ ) polygon[ii] = TextureToElement( polygon[ii] );

						int globalInteriorIndex = localInteriorIndex + gridChart.globalIndexInteriorCellOffset;
						int cellLineId = interiorCellLineIndex[globalInteriorIndex].first;
						int cellLineOffset = interiorCellLineIndex[globalInteriorIndex].second;

						// There is a single sample for the whole polygon
						typename ElementSamples::Bilinear bilinearElementSample( fastIntegration ? 1 : (polygon.size()-2)*Samples );
						bilinearElementSample.cellOffset = cellLineOffset;
						for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) bilinearElementSample.tensor(x,y) = (typename ElementSamples::Real)element_metric_inverse(x,y);

						for( int p=2 ; p<polygon.size() ; p++ )
						{
							Point2D< GeometryReal > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };
							SquareMatrix< GeometryReal , 2 > fragment_to_element_differential;
							for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = dm[x][y];
							GeometryReal fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );
							GeometryReal fragment_area = element_area * fragment_to_element_area_scale_factor / 2;

							for( int s=0 ; s<Samples ; s++ )
							{
								Point2D< typename ElementSamples::Real > pos = polygon[0] + dm[0] * (typename ElementSamples::Real)TriangleIntegrator< Samples >::Positions[s][0] + dm[1] * (typename ElementSamples::Real)TriangleIntegrator< Samples >::Positions[s][1];
								if( !InUnitSquare( pos ) ) MK_THROW( "Sample position out of unit square! (" , pos[0] , " " , pos[1] , ")" );

								typename ElementSamples::Bilinear::SampleData& sampleData = bilinearElementSample[ fastIntegration ? 0 : (p-2)*Samples+s ];
								SetInteriorDuals< Samples >( sampleData , bilinearElementSample , pos , (typename ElementSamples::Real)( TriangleIntegrator< Samples >::Weights[s] * fragment_area ) );
								sampleData.pos = pos;
							}
						}
						if( fastIntegration ) bilinearElementSample[0].pos = Point2D< GeometryReal >( PolygonCenter( &polygon[0] , (unsigned int)polygon.size() ) );
						{
							std::lock_guard< std::mutex > lock( element_samples_bilinear_mutex );
							elementSamples.bilinear[ cellLineId ].push_back( bilinearElementSample );
						}
					}
				}
				// Boundary cell
				else if( localBoundaryIndex!=-1 )
				{
					std::vector< BoundaryIndexedTriangle< GeometryReal > > cellBoundaryTriangles = gridChart.boundaryTriangles[localBoundaryIndex];

					// Iterate over all elements in the cell
					for( int bt=0 ; bt<cellBoundaryTriangles.size() ; bt++ )
					{
						BoundaryIndexedTriangle< GeometryReal > element = cellBoundaryTriangles[bt];
						std::vector< Point2D< GeometryReal > > element_vertices(3);
						for( int ii=0 ; ii<3 ; ii++ ) element_vertices[ii] = TextureToCell( element[ii] );

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
							GeometryReal element_area = sqrt( element_metric.determinant() );

							typename ElementSamples::Quadratic quadraticElementSample( fastIntegration ? 1 : (unsigned int)(polygon.size()-2)*Samples );
							for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) quadraticElementSample.tensor(x,y) = (typename ElementSamples::Real)element_metric_inverse(x,y);
							const QuadraticElementIndex& triangleElementIndices = gridChart.boundaryTriangles[localBoundaryIndex][bt].indices;
							for( int k=0 ; k<6 ; k++ )
							{
								int _fineBoundaryIndex = fineBoundaryIndex[ triangleElementIndices[k] ];
								if( _fineBoundaryIndex!=-1 ) quadraticElementSample.fineNodes[k] = _fineBoundaryIndex;
								else MK_THROW( "Invalid fine boundary index" );
							}

							for( int p=2 ; p<polygon.size() ; p++ )
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
										if( !InUnitTriangle( pos ) ) MK_THROW( "Sample out of unit right triangle! (" , pos[0] , " " , pos[1] , ")" );
										else
										{
											pos[0] = std::max< GeometryReal >( pos[0] , 0 );
											pos[1] = std::max< GeometryReal >( pos[1] , 0 );
											GeometryReal excess = ( pos[0] + pos[1] ) - 1;
											if( excess>0 ) pos[0] -= excess/2 , pos[1] -= excess/2;
										}

										typename ElementSamples::Quadratic::SampleData& sampleData = quadraticElementSample[ fastIntegration ? 0 : (p-2)*Samples+s ];
										SetBoundaryDuals< Samples >( sampleData , quadraticElementSample , Point2D< typename ElementSamples::Real >( pos ) , (typename ElementSamples::Real)( TriangleIntegrator< Samples >::Weights[s] * fragment_area ) );
										sampleData.pos = pos;
									}
								}
								else
								{
									MK_WARN( "Element discarded due to zero mass. Triangle " , t , ". Boundary cell " , localBoundaryIndex , ". Element " , bt , ". Sub triangle " , p-2 );

									printf( "Atlas triangle\n" );
									printf( "%d \n" , 3 );
									for( int v=0 ; v<3 ; v++ ) printf( "%f %f %f\n" , atlasTriangle.vertices[v][0] , atlasTriangle.vertices[v][1] , 0. );
									printf( "%d %d %d\n" , atlasTriangle.atlasVertexIndices[0] , atlasTriangle.atlasVertexIndices[1] , atlasTriangle.atlasVertexIndices[2] );
									printf( "%d %d %d\n" , atlasTriangle.atlasEdgeIndices[0] , atlasTriangle.atlasEdgeIndices[1] , atlasTriangle.atlasEdgeIndices[2] );
									for( int _bt=0 ; _bt<cellBoundaryTriangles.size() ; _bt++ )
									{
										printf( "ELEMENT %d\n" , _bt );

										BoundaryIndexedTriangle< GeometryReal > _element = cellBoundaryTriangles[_bt];
										AtlasIndexedPolygon< GeometryReal > _polygon;
										SetAtlasIndexedPolygonFromBoundaryTriangle( _element , _polygon );

										printf( "Boundary triangle\n" );
										printf( "%d \n" , 3 );
										for( int v=0 ; v<3 ; v++ )printf( "%f %f %f\n" , _element.vertices[v][0] , _element.vertices[v][1] , 0. );
										printf( "%d %d %d\n" , _element.atlasVertexIndices[0] , _element.atlasVertexIndices[1] , _element.atlasVertexIndices[2] );
										printf( "%d %d %d\n" , _element.atlasEdgeIndices[0] , _element.atlasEdgeIndices[1] , _element.atlasEdgeIndices[2] );

										ClipPartiallyIndexedPolygonToIndexedTriangle( _polygon , atlasTriangle , _bt==bt );
									}
								}
							}
							if( fastIntegration ) quadraticElementSample[0].pos = Point2D< GeometryReal >( PolygonCenter( &polygon[0] , (unsigned int)polygon.size() ) );
							{
								std::lock_guard< std::mutex > lock( element_samples_quadratic_mutex );
								elementSamples.quadratic.push_back( quadraticElementSample );
							}
						}
						else if( clippingResult<0 ) MK_THROW( "Bad clipping result: " , clippingResult );
					}
				}
			}
		}
	}

	template< unsigned int Samples , typename GeometryReal , typename ElementSamples >
	void InitializeIntegration
	(
		const std::vector< std::vector< SquareMatrix< GeometryReal , 2 > > >& parameterMetric ,
		const std::vector< AtlasChart< GeometryReal > > &atlasCharts ,
		const std::vector< GridChart< GeometryReal > > &gridCharts ,
		const std::vector< std::pair< int , int > > &interiorCellLineIndex ,
		const std::vector< int > &fineBoundaryIndex ,
		ElementSamples &elementSamples ,
		bool fastIntegration
	)
	{
		std::mutex element_samples_bilinear_mutex , element_samples_quadratic_mutex;
		ThreadPool::ParallelFor
		(
			0 , gridCharts.size() ,
			[&]( unsigned int , size_t i )
			{
				InitializeIntegration< Samples >( parameterMetric[i] , atlasCharts[i] , gridCharts[i] , interiorCellLineIndex , fineBoundaryIndex , elementSamples , element_samples_bilinear_mutex , element_samples_quadratic_mutex , fastIntegration );
			}
		);
	}

	///////////////////////////
	// Compute the integrals //
	///////////////////////////
	template< class Real , typename T , typename ValueFunctionType >
	void IntegrateBilinear
	(
		const BilinearElementScalarSample< Real > &sample ,
		const ValueFunctionType &ValueFunction ,
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
		const BilinearElementGradientSample< Real > &sample ,
		const GradientFunctionType &GradientFunction ,
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

	template< typename Real , typename T , typename ElementSamples , typename SampleFunctionType >
	void Integrate
	(
		const std::vector< InteriorCellLine >& interiorCellLines ,
		const ElementSamples &elementSamples ,
		const std::vector< T >& potential ,
		const std::vector< T >& boundary_potential ,
		const SampleFunctionType& SampleFunction ,
		std::vector< T >& rhs ,
		std::vector< T >& boundary_rhs ,
		bool verbose=false
	)
	{
		Miscellany::Timer timer;
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

				int numSamples = (int)elementSamples.bilinear[r].size();
				const typename ElementSamples::Bilinear *samplePtr = &elementSamples.bilinear[r][0];
				int currentOffset = 0;

				for( int j=0 ; j<numSamples ; j++ )
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

		unsigned int threads = ThreadPool::NumThreads();
		std::vector< int > lineRange( threads+1 );
		int blockSize = (int)interiorCellLines.size() / threads;
		for( unsigned int t=0 ; t<threads ; t++ ) lineRange[t] = t*blockSize;
		lineRange[threads] = (int)interiorCellLines.size();
		ThreadPool::ParallelFor
		(
			0 , threads ,
			[&]( unsigned int , size_t t )
			{
				int firstLine = lineRange[t];
				int lastLine = lineRange[t+1];
				for( int r=firstLine ; r<lastLine ; r++ ) UpdateRow(r);
			}
		);

		if( verbose ) printf( "Integrated bilinear: %.2f(s)\n" , timer.elapsed() );

		timer.reset();
		for( int i=0 ; i<elementSamples.quadratic.size() ; i++ )
		{
			const typename ElementSamples::Quadratic &sample = elementSamples.quadratic[i];
			// The values of the potential at the vertices and edge mid-points
			T cornerValues[] = { boundary_potential[ sample.fineNodes[0] ] , boundary_potential[ sample.fineNodes[1] ] , boundary_potential[ sample.fineNodes[2] ] , boundary_potential[ sample.fineNodes[3] ] , boundary_potential[ sample.fineNodes[4] ] , boundary_potential[ sample.fineNodes[5] ] };
			T rhsValues[] = { T() , T() , T() , T() , T() , T() };
			IntegrateQuadratic< Real , T >( sample , SampleFunction , cornerValues , rhsValues );
			for( int k=0 ; k<6 ; k++ ) boundary_rhs[ sample.fineNodes[k] ] += rhsValues[k];
		}
		if( verbose ) printf( "Integrated quadratic: %.2f(s)\n" , timer.elapsed() );
	}
}