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

#include "Hierarchy.h"
#include "PolygonClipping.h"

template< typename Real >
int InitializeChartMassAndStiffness
(
	const std::vector< SquareMatrix< Real , 2 > >& texture_metrics ,
	const AtlasChart& atlasChart ,
	const GridChart& gridChart ,
	const std::vector< int > & boundaryAndDeepIndex , const std::vector< int >& fineBoundaryIndex ,
	std::vector< Real > & deepMassCoefficients ,
	std::vector< Real > & deepStiffnessCoefficients ,
	std::vector< Eigen::Triplet< Real > >& boundaryBoundaryMassTriplets ,
	std::vector< Eigen::Triplet< Real > >& boundaryBoundaryStiffnessTriplets ,
	std::vector< Eigen::Triplet< Real > >& boundaryDeepMassTriplets ,
	std::vector< Eigen::Triplet< Real > >& boundaryDeepStiffnessTriplets ,
	bool computeCellBasedStiffness ,
	const std::vector< Point3D< Real > >& inputSignal ,
	const std::vector< Point3D< Real > >& boundarySignal ,
	std::vector< Real >& texelToCellCoeffs ,
	std::vector< Eigen::Triplet< Real > > boundaryCellStiffnessTriplets[3]
)
{
	std::vector< SquareMatrix< Real , 4 > > cellStiffness;
	std::vector< SquareMatrix< Real , 4 > > cellMass;
	cellMass.resize(gridChart.numInteriorCells);
	cellStiffness.resize(gridChart.numInteriorCells);

	std::vector< SquareMatrix< Real , 6 > > triangleElementStiffness;
	std::vector< SquareMatrix< Real , 6 > > triangleElementMass;
	triangleElementMass.resize(gridChart.numBoundaryTriangles);
	triangleElementStiffness.resize(gridChart.numBoundaryTriangles);

	////Rasterize
	int zeroAreaElementCount = 0;
	Real precision_error = (Real)1e-3;

	auto InUnitSquare =   [&]( Point2D< Real > p ){ return !( p[0]<0-precision_error || p[1]<0-precision_error || p[0]>1+precision_error || p[1]>1+precision_error ); };
	auto InUnitTriangle = [&]( Point2D< Real > p ){ return !( p[0]<0-precision_error || p[1]<0-precision_error || ( p[0]+p[1] )>1+precision_error ); };
	auto CellInTriangle = [&]( int i , int j , const std::vector< Point2D< Real > >& vertices )
	{
		Real x1 = (Real)i*gridChart.cellSizeW , x2 = (Real)(i+1)*gridChart.cellSizeW;
		Real y1 = (Real)j*gridChart.cellSizeH , y2 = (Real)(j+1)*gridChart.cellSizeH;
		Point2D< Real > points[] = { Point2D< Real >(x1,y1) , Point2D< Real >(x2,y1) , Point2D< Real >(x2,y2) , Point2D< Real >(x1,y2) };
		Point2D< Real > normals[] = { vertices[1]-vertices[0] , vertices[2]-vertices[1] , vertices[0]-vertices[2] };
		for( int k=0 ; k<3 ; k++ ) normals[k] = Point2D< Real >( normals[k][1] , -normals[k][0] );
		Real offsets[] = { Point2D< Real >::Dot( normals[0] , vertices[0] ) , Point2D< Real >::Dot( normals[1] , vertices[1] ) , Point2D< Real >::Dot( normals[2] , vertices[2] ) };
		for( int k=0 ; k<4 ; k++ )
		{
			if( ( Point2D< Real >::Dot( points[k] , normals[0] ) - offsets[0] )<0 ) return false;
			if( ( Point2D< Real >::Dot( points[k] , normals[1] ) - offsets[1] )<0 ) return false;
			if( ( Point2D< Real >::Dot( points[k] , normals[2] ) - offsets[2] )<0 ) return false;
		}
		return true;
	};

	// Compute the square-root of the weights to make taking the weighted dot-product faster
	const Real _integrator4_sampleWeight[] =
	{
		(Real)sqrt( TriangleIntegrator< 6 >::Weights[0] ) ,
		(Real)sqrt( TriangleIntegrator< 6 >::Weights[1] ) ,
		(Real)sqrt( TriangleIntegrator< 6 >::Weights[2] ) ,
		(Real)sqrt( TriangleIntegrator< 6 >::Weights[3] ) ,
		(Real)sqrt( TriangleIntegrator< 6 >::Weights[4] ) ,
		(Real)sqrt( TriangleIntegrator< 6 >::Weights[5] )
	};
	SquareMatrix< Real , 2 > cell_to_texture_differential;
	cell_to_texture_differential(0,0) = gridChart.cellSizeW;
	cell_to_texture_differential(1,1) = gridChart.cellSizeH;
	cell_to_texture_differential(0,1) = cell_to_texture_differential(1,0) = 0.;

	SquareMatrix< Real , 4 > interior_cell_mass , interior_cell_stiffnesses[2][2];
	{
		std::vector< Point2D< Real > > polygon = { Point2D< Real >(0,0) , Point2D< Real >(1,0) , Point2D< Real >(1,1) , Point2D< Real >(0,1) };
		for( int p=2 ; p<polygon.size() ; p++ )
		{
			Point2D< Real > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };
			Point2D< Real > fragment_samples[6];
			for( int s=0 ; s<6 ; s++ )
			{
				fragment_samples[s] = polygon[0] + dm[0] * TriangleIntegrator<6>::Positions[s][0] + dm[1] * TriangleIntegrator<6>::Positions[s][1];
			}

			//Integrate scalar product and gradient field
			Real sampleValues[6][4];
			Point2D< Real > sampleGradients[6][4];
			for( int s=0 ; s<6 ; s++ )
			{
				BilinearElementValuesAndGradients( fragment_samples[s] , sampleValues[s] , sampleGradients[s] );
				for( int k=0 ; k<4 ; k++ )
				{
					sampleValues   [s][k] *= _integrator4_sampleWeight[s];
					sampleGradients[s][k] *= _integrator4_sampleWeight[s];
				}
			}

			for( int k=0 ; k<4 ; k++ ) for( int l=0 ; l<4 ; l++ ) for( int s=0 ; s<6 ; s++ )
			{
				interior_cell_mass(l,k) += sampleValues[s][k] * sampleValues[s][l] / 2.;
				for( int m=0 ; m<2 ; m++ ) for( int n=0 ; n<2 ; n++ ) interior_cell_stiffnesses[m][n](l,k) += sampleGradients[s][l][m] * sampleGradients[s][k][n] / 2.;
			}
		}
	}
	//#pragma omp parallel for 
	for( int t=0 ; t<atlasChart.triangles.size() ; t++ )
	{
		Point2D< Real > tPos[3];
		for( int i=0 ; i<3 ; i++ ) tPos[i] = atlasChart.vertices[ atlasChart.triangles[t][i] ] - gridChart.corner;

		SquareMatrix< Real , 2 > texture_metric = texture_metrics[t];
		SquareMatrix< Real , 2 > cell_metric = cell_to_texture_differential.transpose() * texture_metric * cell_to_texture_differential;
		SquareMatrix< Real , 2 > cell_metric_inverse = cell_metric.inverse();
		Real cell_area_scale_factor = sqrt( cell_metric.determinant() );

		//BBox
		int minCorner[2] , maxCorner[2];
		GetTriangleIntegerBBox( tPos , 1./gridChart.cellSizeW , 1./gridChart.cellSizeH , minCorner , maxCorner );

		std::vector< Point2D< Real > > parametricVertices(3);
		parametricVertices[0] = tPos[0], parametricVertices[1] = tPos[1], parametricVertices[2] = tPos[2];

		AtlasIndexedTriangle atlasTriangle;
		for( int k=0 ; k<3 ; k++ )
		{
			atlasTriangle.vertices[k] = tPos[k];
			atlasTriangle.atlasEdgeIndices[k] = atlasChart.atlasEdgeIndices[3 * t + k];
			atlasTriangle.atlasVertexIndices[k] = atlasChart.triangles[t][k];
			atlasTriangle.atlasVertexParentEdge[k] = -1;
		}

		for( int j=minCorner[1] ; j<maxCorner[1] ; j++ ) for( int i=minCorner[0] ; i<maxCorner[0] ; i++ )
		{
			auto TextureToCell = [&]( Point2D< Real > p ){ return Point2D< Real >( ( p[0] / gridChart.cellSizeW ) - i , ( p[1] / gridChart.cellSizeH ) - j ); };

			int localInteriorIndex = gridChart.localInteriorCellIndex(i,j) , localBoundaryIndex = gridChart.localBoundaryCellIndex(i,j);
			if( localInteriorIndex!=-1 && localBoundaryIndex!=-1 ){ printf( "[ERROR] Cell simultaneosly interior and boundary!\n" ) ; return 0; }

			// If the cell is entirely within the triangle...
			if( CellInTriangle( i , j , parametricVertices ) && localInteriorIndex!=-1 )
			{
				SquareMatrix< Real , 4 > polygonMass , polygonStiffness;

				polygonMass = interior_cell_mass * cell_area_scale_factor;
				for( int k=0 ; k<4 ; k++ ) for( int l=0 ; l<4 ; l++ )
				{
					Real v = 0;
					for( int m=0 ; m<2 ; m++ ) for( int n=0 ; n<2 ; n++ ) v += cell_metric_inverse(m,n) * interior_cell_stiffnesses[m][n](k,l);
					polygonStiffness(k,l) = v * cell_area_scale_factor;
				}
				cellMass     [ localInteriorIndex ] += polygonMass;
				cellStiffness[ localInteriorIndex ] += polygonStiffness;
			}
			else if( localInteriorIndex!=-1 )
			{
				// For interior cells, the cell and the element are the same thing
				auto TextureToElement = TextureToCell;
				SquareMatrix< Real , 2 > element_metric = cell_metric , element_metric_inverse = cell_metric_inverse;
				Real element_area_scale_factor = cell_area_scale_factor;

				CellClippedTriangle polygon = parametricVertices;

				// Clip the triangle to the cell
				if( ClipTriangleToPrimalCell( polygon , i , j , gridChart.cellSizeW , gridChart.cellSizeH ) )
				{
					// Transform the polygon vertices into the coordinate frame of the cell
					for( int i=0 ; i<polygon.size() ; i++ ) polygon[i] = TextureToElement( polygon[i] );
					SquareMatrix< Real , 4 > polygonStiffness , polygonMass;
					for( int p=2 ; p<polygon.size() ; p++ )
					{
						Point2D< Real > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };

						SquareMatrix< Real , 2 > fragment_to_element_differential;
						for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = dm[x][y];
						Real fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );

						Point2D< Real > fragment_samples[6];
						for( int s=0 ; s<6 ; s++ )
						{
							fragment_samples[s] = polygon[0] + dm[0] * TriangleIntegrator<6>::Positions[s][0] + dm[1] * TriangleIntegrator<6>::Positions[s][1];
							if( !InUnitSquare( fragment_samples[s] ) ){ printf( "[ERROR] Interior sample out of unit box! (%f %f)\n" , fragment_samples[s][0] , fragment_samples[s][1]) ; return 0; }
						}

						// Integrate scalar product and gradient field
						// Make the code more efficient by:
						// -- pre-multiplying the values and gradients by the square-root of the quadrature weight
						// -- pre-multiplying the gradients by the inverse of the metric
						// so that the computation within the inner loop is faster.
						Real fragment_area = element_area_scale_factor * fragment_to_element_area_scale_factor / 2.0;
						Real sampleValues[6][4];
						Point2D< Real > sampleGradients[6][4] , _sampleGradients[6][4];
						for( int s=0 ; s<6 ; s++ )
						{
							BilinearElementValuesAndGradients( fragment_samples[s] , sampleValues[s] , sampleGradients[s] );
							for( int k=0 ; k<4 ; k++ )
							{
								sampleValues[s][k] *= _integrator4_sampleWeight[s];
								sampleGradients[s][k] *= _integrator4_sampleWeight[s];
								_sampleGradients[s][k] = element_metric_inverse * sampleGradients[s][k];
							}
						}
						for( int k=0 ; k<4 ; k++ ) for( int l=0 ; l<4 ; l++ )
						{
							Real vIntegral=0 , gIntegral=0;
							for( int s=0 ; s<6 ; s++ )
							{
								vIntegral += sampleValues[s][k] * sampleValues[s][l];
								gIntegral += Point2D< Real >::Dot( sampleGradients[s][l] , _sampleGradients[s][k] );
							}
							polygonMass(l,k) += vIntegral * fragment_area;
							polygonStiffness(l,k) += gIntegral * fragment_area;
						}
					}
					cellMass[localInteriorIndex] += polygonMass;
					cellStiffness[localInteriorIndex] += polygonStiffness;
				}
			}
			else if( localBoundaryIndex!=-1 )
			{
				Real fine_precision_error = 1e-4;

				std::vector< BoundaryIndexedTriangle > cellBoundaryTriangles = gridChart.boundaryTriangles[localBoundaryIndex];

				// Iterate over all elements in the cell
				for( int bt=0 ; bt<cellBoundaryTriangles.size() ; bt++ )
				{
					BoundaryIndexedTriangle element = cellBoundaryTriangles[bt];
					std::vector< Point2D< Real > > element_vertices(3);
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

						SquareMatrix< Real , 2 > element_to_cell_differential , cell_to_element_differential;
						Point2D< Real > dm[] = { element_vertices[1]-element_vertices[0] , element_vertices[2]-element_vertices[0] };
						for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) element_to_cell_differential(x,y) = dm[x][y];
						cell_to_element_differential = element_to_cell_differential.inverse();
						auto CellToElement = [&]( Point2D< Real > p ){ return cell_to_element_differential * ( p - element_vertices[0] ); };
						// Convert the polygon vertices from the cell frame to the element frame
						for( int i=0 ; i<polygon.size() ; i++ ) polygon[i] = CellToElement( polygon[i] );

						SquareMatrix< Real , 2 > element_metric = element_to_cell_differential.transpose() * cell_metric * element_to_cell_differential;
						SquareMatrix< Real , 2 > element_metric_inverse = element_metric.inverse();
						Real element_area_scale_factor = sqrt( element_metric.determinant() );
						SquareMatrix< Real , 6 > polygonStiffness , polygonMass;

						Real polygonArea = 0;

						for( int p=2 ; p<polygon.vertices.size() ; p++ )
						{
							Point2D< Real > dm[] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };

							SquareMatrix< Real , 2 > fragment_to_element_differential;
							for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = dm[x][y];
							Real fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );
							Real fragment_area = fragment_to_element_area_scale_factor * element_area_scale_factor / 2.;
							if( fragment_area>0 )
							{
								polygonArea += fragment_area;

								Point2D< Real > fragment_samples[6]; //Belong to unit right triangles

								for( int s=0 ; s<6 ; s++ )
								{
									fragment_samples[s] = polygon.vertices[0] + dm[0] * TriangleIntegrator<6>::Positions[s][0] + dm[1] * TriangleIntegrator<6>::Positions[s][1];
									if( !InUnitTriangle( fragment_samples[s] ) ){ printf( "[ERROR] Boundary sample out of unit right triangle! (%f %f)\n" , fragment_samples[s][0] , fragment_samples[s][1] ) ; return 0; }
									else
									{
										fragment_samples[s][0] = std::max< Real >( fragment_samples[s][0] , 0 );
										fragment_samples[s][1] = std::max< Real >( fragment_samples[s][1] , 0 );
										Real excess = ( fragment_samples[s][0] + fragment_samples[s][1] ) - 1;
										if( excess>0 ) fragment_samples[s][0] -= excess/2 , fragment_samples[s][1] -= excess/2;
									}
								}

								Real sampleValues[6][6];
								Point2D< Real > sampleGradients[6][6] , _sampleGradients[6][6];

								for( int s=0 ; s<6 ; s++ )
								{
									QuadraticElementValuesAndGradients( fragment_samples[s] , sampleValues[s] , sampleGradients[s] );
									for( int k=0 ; k<6 ; k++ )
									{
										sampleValues[s][k] *= _integrator4_sampleWeight[s];
										sampleGradients[s][k] *= _integrator4_sampleWeight[s];
										_sampleGradients[s][k] = element_metric_inverse * sampleGradients[s][k];										
									}
								}

								for( int k=0 ; k<6 ; k++ ) for( int l=0 ; l<6 ; l++ )
								{
									Real vIntegral=0 , gIntegral=0;
									for( int s=0 ; s<6 ; s++ )
									{
										vIntegral += sampleValues[s][k] * sampleValues[s][l];
										gIntegral += Point2D< Real >::Dot( sampleGradients[s][k] , _sampleGradients[s][l] );
									}
									polygonMass( l , k ) += vIntegral * fragment_area;
									polygonStiffness( l , k ) += gIntegral * fragment_area;
								}
							}
							else
							{
								zeroAreaElementCount++;
								printf( "[WARNING] Zero area polygon at cell %d %d \n", gridChart.cornerCoords[0] + i, gridChart.cornerCoords[1] + j);
							}
						}

						Real integratedPolygonMass = 0;
						for( int k=0 ; k<6 ; k++ ) for( int l=0 ; l<6 ; l++ ) integratedPolygonMass += polygonMass(k,l);
						if( fabs( integratedPolygonMass - polygonArea )>precision_error )
						{
							printf( "[ERROR] out of precision! \n");
							//return 0;
						}
						//#pragma omp critical
						{
							for (int dk = 0; dk < 6; dk++)for (int dl = 0; dl < 6; dl++) {
								triangleElementStiffness[boundaryTriangleId](dk, dl) += (polygonStiffness(dk, dl) + polygonStiffness(dl, dk)) / 2.0;
								triangleElementMass[boundaryTriangleId](dk, dl) += (polygonMass(dk, dl) + polygonMass(dl, dk)) / 2.0;
							}
						}
					}
				}
			}
		}
	}
	if( zeroAreaElementCount ) printf( "[WARNING] Element with zero area = %d\n" , zeroAreaElementCount );

	int offset_i[4] = { 0, 1 , 1 ,0 };
	int offset_j[4] = { 0, 0 , 1 ,1 };
	int cellOffset[4] = { 3, 2, 0, 1 };

	auto NeighbourOffset = [&](int k, int l) {
		return  (offset_j[l] - offset_j[k] + 1) * 3 + (offset_i[l] - offset_i[k] + 1);
	};

	for (int i = 0; i < gridChart.interiorCellCorners.size(); i++) {
		const CellIndex & indicesGlobal = gridChart.interiorCellGlobalCorners[i];
		const CellIndex & indicesInterior = gridChart.interiorCellCorners[i];

		const int localCellIndex = gridChart.interiorCellIndexToLocalCellIndex[i];
		const int globalCellIndex = localCellIndex + gridChart.globalIndexCellOffset;

		Point< Real , 4 > prod[3];
		if (computeCellBasedStiffness){
			Point3D < Real > values[4] = { inputSignal[indicesGlobal[0]],inputSignal[indicesGlobal[1]], inputSignal[indicesGlobal[2]], inputSignal[indicesGlobal[3]] };
			for (int c = 0; c < 3; c++) {
				Point<Real, 4> v;
				v[0] = values[0][c];
				v[1] = values[1][c];
				v[2] = values[2][c];
				v[3] = values[3][c];

				prod[c] = cellStiffness[i] * v;
			}
		}

		for (int k = 0; k < 4; k++) {
			int currentNode = indicesGlobal[k];
			int _currentBoundaryAndDeepIndex = boundaryAndDeepIndex[currentNode];
			if (_currentBoundaryAndDeepIndex < 0) {//Deep
				int deepIndex = -_currentBoundaryAndDeepIndex - 1;
				for (int l = 0; l < 4; l++) {
					//offset_i(k,l);
					//offset_j(k,l);
					//int neighbourIndex = (neighbourNode.cj - currentNode.cj + 1) * 3 + (neighbourNode.ci - currentNode.ci + 1);
					deepMassCoefficients[10 * deepIndex + NeighbourOffset(k, l)] += cellMass[i](k, l);
					deepStiffnessCoefficients[10 * deepIndex + NeighbourOffset(k, l)] += cellStiffness[i](k, l);
				}
				if (computeCellBasedStiffness) {
					//Add cell data
					texelToCellCoeffs[3 * (4 * deepIndex + cellOffset[k]) + 0] = prod[0][k];
					texelToCellCoeffs[3 * (4 * deepIndex + cellOffset[k]) + 1] = prod[1][k];
					texelToCellCoeffs[3 * (4 * deepIndex + cellOffset[k]) + 2] = prod[2][k];
				}
			}
			else {//Bundary
				int boundaryIndex = _currentBoundaryAndDeepIndex - 1;
				for (int l = 0; l < 4; l++) {
					int neighbourNode = indicesGlobal[l];
					int _neighbourBoundaryAndDeepIndex = boundaryAndDeepIndex[neighbourNode];
					if (_neighbourBoundaryAndDeepIndex < 0) {//Deep
						boundaryDeepMassTriplets.push_back(Eigen::Triplet< Real >(boundaryIndex, neighbourNode, cellMass[i](k, l)));
						boundaryDeepStiffnessTriplets.push_back(Eigen::Triplet< Real >(boundaryIndex, neighbourNode, cellStiffness[i](k, l)));
					}
					else {//Boundary
						boundaryBoundaryMassTriplets.push_back(Eigen::Triplet< Real >(fineBoundaryIndex[indicesInterior[k]], fineBoundaryIndex[indicesInterior[l]], cellMass[i](k, l)));
						boundaryBoundaryStiffnessTriplets.push_back(Eigen::Triplet< Real >(fineBoundaryIndex[indicesInterior[k]], fineBoundaryIndex[indicesInterior[l]], cellStiffness[i](k, l)));
					}
				}
				if (computeCellBasedStiffness) {
					//Add cell data
					int _fineBoundaryIndex = fineBoundaryIndex[indicesInterior[k]];
					for (int c = 0; c < 3; c++) {
						boundaryCellStiffnessTriplets[c].push_back(Eigen::Triplet< Real >(_fineBoundaryIndex, globalCellIndex, prod[c][k]));
					}
				}
			}

		}
	}

	for (int c = 0; c < gridChart.boundaryTriangles.size(); c++) {

		const int localCellIndex = gridChart.boundaryCellIndexToLocalCellIndex[c];
		const int globalCellIndex = localCellIndex + gridChart.globalIndexCellOffset;

		for (int b = 0; b < gridChart.boundaryTriangles[c].size(); b++) {
			const int i = gridChart.boundaryTriangles[c][b].id;
			const TriangleElementIndex & indices = gridChart.boundaryTriangles[c][b].indices;
			TriangleElementIndex fineTriangleElementIndices;
			for (int k = 0; k < 6; k++) fineTriangleElementIndices[k] = fineBoundaryIndex[indices[k]];

			//Add cell data
			
			if (computeCellBasedStiffness){
				Point< Real , 6 > prod[3];
				Point3D < Real > values[6] = { boundarySignal[fineTriangleElementIndices[0]],boundarySignal[fineTriangleElementIndices[1]],boundarySignal[fineTriangleElementIndices[2]],boundarySignal[fineTriangleElementIndices[3]],boundarySignal[fineTriangleElementIndices[4]],boundarySignal[fineTriangleElementIndices[5]] };
				for (int c = 0; c < 3; c++) {
					Point< Real , 6 > v;
					v[0] = values[0][c];
					v[1] = values[1][c];
					v[2] = values[2][c];
					v[3] = values[3][c];
					v[4] = values[4][c];
					v[5] = values[5][c];
					prod[c] = triangleElementStiffness[i] * v;
				}

				for (int k = 0; k < 6; k++) {
					for (int c = 0; c < 3; c++) {
						boundaryCellStiffnessTriplets[c].push_back(Eigen::Triplet< Real >(fineTriangleElementIndices[k], globalCellIndex, prod[c][k]));
					}
				}
			}
			for (int k = 0; k < 6; k++)for (int l = 0; l < 6; l++) {
				boundaryBoundaryMassTriplets.push_back(Eigen::Triplet< Real >(fineTriangleElementIndices[k], fineTriangleElementIndices[l], triangleElementMass[i](l, k)));
				boundaryBoundaryStiffnessTriplets.push_back(Eigen::Triplet< Real >(fineTriangleElementIndices[k], fineTriangleElementIndices[l], triangleElementStiffness[i](l, k)));
			}
		}
	}

	return 1;
}

template< typename Real >
int InitializeMassAndStiffness( const std::vector<std::vector<SquareMatrix< Real , 2 > > >& parameterMetric , const std::vector< AtlasChart >& atlasCharts , const GridAtlas& gridAtlas , const std::vector< int >& fineBoundaryIndex , const int numFineBoundaryNodes ,
	std::vector< Real>& deepMassCoefficients, std::vector< Real>& deepStiffnessCoefficients ,
	SparseMatrix< Real , int >& boundaryBoundaryMassMatrix , SparseMatrix< Real , int >& boundaryBoundaryStiffnessMatrix ,
	SparseMatrix< Real , int >& boundaryDeepMassMatrix , SparseMatrix< Real , int >& boundaryDeepStiffnessMatrix ,
	bool computeCellBasedStiffness , const std::vector< Point3D< Real > >& inputSignal , const std::vector< Point3D< Real > >& boundarySignal , std::vector< Real >& texelToCellCoeffs , SparseMatrix< Real , int > boundaryCellBasedStiffnessRHSMatrix[3] )
{

	const std::vector<GridChart> & gridCharts = gridAtlas.gridCharts;
	const std::vector<int> & boundaryAndDeepIndex = gridAtlas.boundaryAndDeepIndex;

	if(computeCellBasedStiffness) texelToCellCoeffs.resize(3 * 4 * gridAtlas.numDeepTexels);

	clock_t m_begin;
	m_begin = clock();

	std::vector< Eigen::Triplet< Real > > boundaryBoundaryMassTriplets;
	std::vector< Eigen::Triplet< Real > > boundaryBoundaryStiffnessTriplets;
	std::vector< Eigen::Triplet< Real > > boundaryDeepMassTriplets;
	std::vector< Eigen::Triplet< Real > > boundaryDeepStiffnessTriplets;
	std::vector< Eigen::Triplet< Real > > boundaryCellStiffnessTriplets[3];

#pragma omp parallel for
	for (int i = 0; i < gridCharts.size(); i++) {
		std::vector< Eigen::Triplet< Real > > chartBoundaryBoundaryMassTriplets;
		std::vector< Eigen::Triplet< Real > > chartBoundaryBoundaryStiffnessTriplets;
		std::vector< Eigen::Triplet< Real > > chartBoundaryDeepMassTriplets;
		std::vector< Eigen::Triplet< Real > > chartBoundaryDeepStiffnessTriplets;
		std::vector< Eigen::Triplet< Real > > chartBoundaryCellStiffnessTriplets[3];
		InitializeChartMassAndStiffness(parameterMetric[i], atlasCharts[i], gridCharts[i], boundaryAndDeepIndex, fineBoundaryIndex, deepMassCoefficients, deepStiffnessCoefficients, chartBoundaryBoundaryMassTriplets, chartBoundaryBoundaryStiffnessTriplets,
			chartBoundaryDeepMassTriplets, chartBoundaryDeepStiffnessTriplets,
			computeCellBasedStiffness ,inputSignal, boundarySignal, texelToCellCoeffs, chartBoundaryCellStiffnessTriplets);
#pragma omp critical
		{
			boundaryBoundaryMassTriplets.insert(boundaryBoundaryMassTriplets.end(), chartBoundaryBoundaryMassTriplets.begin(), chartBoundaryBoundaryMassTriplets.end());
			boundaryBoundaryStiffnessTriplets.insert(boundaryBoundaryStiffnessTriplets.end(), chartBoundaryBoundaryStiffnessTriplets.begin(), chartBoundaryBoundaryStiffnessTriplets.end());
			boundaryDeepMassTriplets.insert(boundaryDeepMassTriplets.end(), chartBoundaryDeepMassTriplets.begin(), chartBoundaryDeepMassTriplets.end());
			boundaryDeepStiffnessTriplets.insert(boundaryDeepStiffnessTriplets.end(), chartBoundaryDeepStiffnessTriplets.begin(), chartBoundaryDeepStiffnessTriplets.end());
			if (computeCellBasedStiffness) for (int c = 0; c < 3; c++) boundaryCellStiffnessTriplets[c].insert(boundaryCellStiffnessTriplets[c].end(), chartBoundaryCellStiffnessTriplets[c].begin(), chartBoundaryCellStiffnessTriplets[c].end());
		}
	}

	int numTexels = gridAtlas.numTexels;
	int numBoundaryTexels = gridAtlas.numBoundaryTexels;
	boundaryBoundaryMassMatrix = SetSparseMatrix( boundaryBoundaryMassTriplets , numFineBoundaryNodes , numFineBoundaryNodes , true );
	boundaryBoundaryStiffnessMatrix = SetSparseMatrix( boundaryBoundaryStiffnessTriplets , numFineBoundaryNodes , numFineBoundaryNodes , true );
	boundaryDeepMassMatrix = SetSparseMatrix( boundaryDeepMassTriplets , numBoundaryTexels , numTexels , false );
	boundaryDeepStiffnessMatrix = SetSparseMatrix( boundaryDeepStiffnessTriplets , numBoundaryTexels , numTexels , false );

	if (computeCellBasedStiffness){
		for (int c = 0; c < 3; c++) {
			boundaryCellBasedStiffnessRHSMatrix[c] = SetSparseMatrix( boundaryCellStiffnessTriplets[c] , numFineBoundaryNodes , gridAtlas.numCells , false );
		}
	}

	return 1;
}

template< typename Real >
int InitializeMassAndStiffness( std::vector< Real >& deepMassCoefficients , std::vector< Real > & deepStiffnessCoefficients ,
	SparseMatrix< Real , int >& boundaryBoundaryMassMatrix , SparseMatrix< Real , int >& boundaryBoundaryStiffnessMatrix ,
	SparseMatrix< Real , int >& boundaryDeepMassMatrix , SparseMatrix< Real , int >& boundaryDeepStiffnessMatrix ,
	const HierarchicalSystem& hierarchy , const std::vector< std::vector< SquareMatrix< Real , 2 > > >& parameterMetric , const std::vector< AtlasChart >& atlasCharts ,
	const BoundaryProlongationData& boundaryProlongation , bool computeCellBasedStiffness ,
	const std::vector< Point3D< Real > >& inputSignal , std::vector< Real >& texelToCellCoeffs , SparseMatrix< Real , int > boundaryCellBasedStiffnessRHSMatrix[3] )
{

	//(2) Initialize mass and stiffness
	deepMassCoefficients.resize(10 * hierarchy.gridAtlases[0].numDeepTexels, 0);
	deepStiffnessCoefficients.resize(10 * hierarchy.gridAtlases[0].numDeepTexels, 0);


	SparseMatrix< Real , int > fineBoundaryBoundaryMassMatrix;
	SparseMatrix< Real , int > fineBoundaryBoundaryStiffnessMatrix;

	SparseMatrix< Real , int > fineBoundaryCellStiffnessRHSMatrix[3];
	std::vector< Point3D< Real > > fineBoundarySignal;
	
	if (computeCellBasedStiffness){
		const std::vector<int> & boundaryGlobalIndex = hierarchy.gridAtlases[0].boundaryGlobalIndex;
		int numBoundaryTexels = (int)boundaryGlobalIndex.size();
		int numFineBoundaryNodes = boundaryProlongation.numFineBoundarNodes;
		std::vector< Point3D< Real > > coarseBoundarySignal;
		coarseBoundarySignal.resize(numBoundaryTexels);
		for (int i = 0; i < numBoundaryTexels; i++)  coarseBoundarySignal[i] = inputSignal[boundaryGlobalIndex[i]];
		fineBoundarySignal.resize(numFineBoundaryNodes);
		boundaryProlongation.coarseBoundaryFineBoundaryProlongation.Multiply(&coarseBoundarySignal[0], &fineBoundarySignal[0]);
	}

	//clock_t m_begin = clock();
	if( !InitializeMassAndStiffness( parameterMetric , atlasCharts , hierarchy.gridAtlases[0] , boundaryProlongation.fineBoundaryIndex , boundaryProlongation.numFineBoundarNodes , deepMassCoefficients , deepStiffnessCoefficients , fineBoundaryBoundaryMassMatrix , fineBoundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix ,
		computeCellBasedStiffness , inputSignal , fineBoundarySignal , texelToCellCoeffs , fineBoundaryCellStiffnessRHSMatrix ) )
	{
		printf("ERROR: Unable to initialize fine mass and stiffness! \n");
		return 0;
	}

	SparseMatrix< Real , int > temp = fineBoundaryBoundaryMassMatrix * boundaryProlongation.coarseBoundaryFineBoundaryProlongation;
	boundaryBoundaryMassMatrix = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction * temp;

	temp = fineBoundaryBoundaryStiffnessMatrix * boundaryProlongation.coarseBoundaryFineBoundaryProlongation;
	boundaryBoundaryStiffnessMatrix = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction * temp;

	if( 1 )
	{
		std::vector< Real > in( boundaryBoundaryMassMatrix.Rows() , 1.0 );
		std::vector< Real > out( boundaryBoundaryMassMatrix.Rows() , 0.0 );
		boundaryBoundaryMassMatrix.Multiply( GetPointer(in) , GetPointer(out) );
		for( int i=0 ; i<out.size() ; i++ ) if( out[i]==0.0 ) printf( "[WARNING] Zero row at index %d. Try running with jittering.\n" , i );
	}

	if( computeCellBasedStiffness ) for( int c=0 ; c<3 ; c++ ) boundaryCellBasedStiffnessRHSMatrix[c] = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction * fineBoundaryCellStiffnessRHSMatrix[c];
	return 1;
}

