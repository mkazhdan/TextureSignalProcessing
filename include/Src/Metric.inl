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

enum
{
	EMBEDDING_METRIC,
	UNIFORM_METRIC,
	UMBILIC_METRIC,
	EDGE_AWARE_METRIC,
};



template< typename Real >
int InitializeVectorFieldMetric( const std::vector< SquareMatrix< Real , 2 > >& embeddingMetric , const std::vector< Point2D< Real > >& vf , bool normalizeArea , std::vector< SquareMatrix< Real , 2 > >& outputMetric )
{
	int tCount = (int)embeddingMetric.size();

	//Compute vf scale
	Real vfScale = 0;
	for( int t=0 ; t<tCount ; t++ ) vfScale += Point2D< Real >::Dot( vf[t] , embeddingMetric[t] * vf[t] ) * (Real)( sqrt( embeddingMetric[t].determinant() )/2.0 );
	Real vfNormalization = (Real)( 1.0 / sqrt(vfScale) );

	outputMetric.resize( embeddingMetric.size() );
	Real totalMass = 0;

	for (int t = 0; t < tCount; t++) {
		SquareMatrix< Real , 2 > g = embeddingMetric[t];
		SquareMatrix< Real , 2 > u;
		if (Point2D< Real >::Dot( vf[t] , embeddingMetric[t] * vf[t] )>0 )
		{
			Point2D< Real > vectors[2];
			vectors[0] = vf[t];
			Point2D< Real > temp = g*vectors[0];
			vectors[1] = Point2D< Real>( -temp[1] , temp[0] );
			for (int e = 0; e < 2; e++) vectors[e] /= sqrt( Point2D< Real >::Dot(vectors[e], g*vectors[e]));

			Real principalScale = (Real)sqrt( Point2D< Real >::Dot( vf[t] , embeddingMetric[t] * vf[t] ) ) * vfNormalization;

			SquareMatrix< Real , 2 > B;
			B(0, 0) = vectors[0][0]; B(0, 1) = vectors[0][1];
			B(1, 0) = vectors[1][0]; B(1, 1) = vectors[1][1];


			SquareMatrix< Real , 2 > D;
			D(0, 0) = 0;
#if 0
			D(1, 1) = principalScale;
#else
			D(1, 1) = 1e-1;
#endif
			D(0, 1) = D(1, 0) = 0.0;

			u = g * B * D * B.inverse();
		}
		else {
			u *= 0;
		}
		outputMetric[t] = u + g*(1e-5);

		totalMass += sqrt(outputMetric[t].determinant()) / 2;
	}
	if( normalizeArea ) for( int t=0 ; t<tCount ; t++ )	outputMetric[t] /= totalMass;

	return 1;
}

#define NORMALIZE_SURFACE_EMBEDDING 1

template< typename Real >
int InitializeCurvatureBasedMetric( const TexturedMesh & mesh , bool normalizeArea , std::vector< SquareMatrix< Real , 2 > >& outputMetric , std::vector< Point3D< Real > >& vNormals , int metricMode )
{
	outputMetric.resize(mesh.triangles.size());

#if NORMALIZE_SURFACE_EMBEDDING
	Real inputMass = GetMeshArea(mesh);
	Real edgeScaling = 1.0 / sqrt(inputMass);
#endif

	Real totalMass = 0;

	for (int t = 0; t < mesh.triangles.size(); t++) {

		Point3D< Real > vPos[3];
		for (int i = 0; i < 3; i++) vPos[i] = mesh.vertices[mesh.triangles[t][i]];

		SquareMatrix< Real , 2> g;
#if NORMALIZE_SURFACE_EMBEDDING
		Point3D< Real > dv[2] = { (vPos[1] - vPos[0])*edgeScaling, (vPos[2] - vPos[0])*edgeScaling };
#else
		Point3D< Real > dv[2] = { (vPos[1] - vPos[0]), (vPos[2] - vPos[0]) };
#endif
		for (int k = 0; k < 2; k++)for (int l = 0; l < 2; l++)g(k, l) = Point3D< Real >::Dot(dv[k], dv[l]);

		Point3D< Real > vNormal[3];
		for (int i = 0; i < 3; i++) vNormal[i] = vNormals[mesh.triangles[t][i]];
		Point3D< Real > dn[2] = { vNormal[1] - vNormal[0], vNormal[2] - vNormal[0] };
		SquareMatrix< Real , 2 > gg;
		for (int k = 0; k < 2; k++)for (int l = 0; l < 2; l++) gg(k, l) = Point3D< Real >::Dot( dn[k] , dv[l] );
		gg(0, 1) = gg(1, 0) = (gg(0, 1) + gg(1, 0)) / 2.0;

		SquareMatrix< Real , 2 > S = g.inverse()*gg;

		if (metricMode == UMBILIC_METRIC) {
			Real a = (Real)1.0;
			Real b = -(S(0, 0) + S(1, 1));
			Real c = S(0, 0)*S(1, 1) - S(0, 1)*S(1, 0);
			Real discriminant = (Real)( b*b - 4.0 *a*c );
			if( discriminant<0 ){ fprintf( stderr , "[ERROR] Unexpected negative discriminant!\n" ) ; return 0; }

			discriminant = sqrt(discriminant);
			Real roots[2] = { (-b - discriminant) / (Real)(2.0*a), (-b + discriminant) / (Real)(2.0*a) };
			Point2D< Real > vectors[2] = { Point2D< Real >(1,0) , Point2D< Real >(0,1) };
			if (S(1, 0) != 0 || S(0, 1) != 0) {
				if (fabs(S(1, 0)) > fabs(S(0, 1))) {
					vectors[0] = Point2D< Real >( -S(1,0), S(0,0) - roots[0] );
				}
				else {
					vectors[0] = Point2D< Real >( -( S(1,1)-roots[0] ) , S(0,1) );
				}
			}

			Point2D< Real > temp = g*vectors[0];
			vectors[1] = Point2D< Real >(-temp[1], temp[0]);

			for (int e = 0; e < 2; e++) vectors[e] /= (Real)sqrt(Point2D< Real>::Dot(vectors[e], g*vectors[e]));

			SquareMatrix< Real , 2 > B;
			B(0, 0) = vectors[0][0]; B(0, 1) = vectors[0][1];
			B(1, 0) = vectors[1][0]; B(1, 1) = vectors[1][1];

			//Major curves
			roots[0] = fabs(roots[1] - roots[0]);
			roots[1] = 0;

			//Minor curves
			//roots[1] = fabs(roots[1] - roots[0]);
			//roots[0] = 0;

			SquareMatrix< Real , 2 > D;
			D(0,0) = roots[0]; D(1,1) = roots[1];
			D(0,1) = D(1,0) = 0;

			SquareMatrix< Real , 2 > u = g * B * D * B.inverse();

			//g = u + g*(1e-10);
			g = u + g*(1e-8); //For Camel LIC minor
		}
		else {
			Real trace = S.trace();
			Real det = S.determinant();
			Real squared_curvatures = trace*trace - 2.0*det;
			if (squared_curvatures < 0) {
				printf("Unexpected negative sum of curvatures!\n");
				return 0;
			}
			g *= (1e-10 + squared_curvatures);
		}

		totalMass += sqrt(g.determinant()) / 2;

		outputMetric[t] = g;
	}

	printf("Here Total Mass %g \n", totalMass);

	if (normalizeArea) {
		for (int t = 0; t < mesh.triangles.size(); t++) {
			outputMetric[t] /= totalMass;
		}
	}
	return 1;
}

template< typename Real >
int InitializePrincipalCurvatureDirection( const TexturedMesh& mesh , std::vector< Point2D< Real > >& principalDirection , std::vector< Point3D< Real > >& vNormals , bool useSmallestCurvarture=true )
{
	principalDirection.resize(mesh.triangles.size());

#if NORMALIZE_SURFACE_EMBEDDING
	Real inputMass = GetMeshArea(mesh);
	Real edgeScaling = (Real)( 1.0 / sqrt(inputMass) );
#endif
	Real totalMass = 0;
	for (int t = 0; t < mesh.triangles.size(); t++) {

		Point3D< Real > vPos[3];
		for (int i = 0; i < 3; i++) vPos[i] = mesh.vertices[mesh.triangles[t][i]];

		SquareMatrix< Real , 2 > g;
#if NORMALIZE_SURFACE_EMBEDDING
		Point3D< Real > dv[2] = { (vPos[1] - vPos[0])*edgeScaling, (vPos[2] - vPos[0])*edgeScaling };
#else
		Point3D< Real > dv[2] = { (vPos[1] - vPos[0]), (vPos[2] - vPos[0]) };
#endif
		for (int k = 0; k < 2; k++)for (int l = 0; l < 2; l++)g(k, l) = Point3D< Real >::Dot(dv[k], dv[l]);

		Point3D< Real > vNormal[3];
		for (int i = 0; i < 3; i++) vNormal[i] = vNormals[mesh.triangles[t][i]];
		Point3D< Real > dn[2] = { vNormal[1] - vNormal[0], vNormal[2] - vNormal[0] };
		SquareMatrix< Real , 2 > gg;
		for (int k = 0; k < 2; k++)for (int l = 0; l < 2; l++) gg(k, l) = Point3D< Real >::Dot(dn[k], dv[l]);
		gg(0, 1) = gg(1, 0) = (gg(0, 1) + gg(1, 0)) / 2.0;

		SquareMatrix< Real , 2 > S = g.inverse()*gg;

		Real a = (Real)1.0;
		Real b = -(S(0, 0) + S(1, 1));
		Real c = S(0, 0)*S(1, 1) - S(0, 1)*S(1, 0);
		Real discriminant = (Real)( b*b - 4.0 *a*c );
		if (discriminant < 0) {
			printf("Unexpected negative discriminant! \n");
			return 0;
		}

		discriminant = sqrt(discriminant);
		Real root = (-b - discriminant) / (2.0*a);
		Point2D< Real > vector =  Point2D< Real >(1,0);
		if (S(1, 0) != 0 || S(0, 1) != 0) {
			if (fabs(S(1, 0)) > fabs(S(0, 1))) {
				vector = Point2D<Real>(-S(1, 0), S(0, 0) - root);
			}
			else {
				vector = Point2D<Real>(-(S(1, 1) - root), S(0, 1));
			}
		}
		if (useSmallestCurvarture) {
			principalDirection[t] = vector;
		}
		else {
			Point2D<Real> temp = g*vector;
			principalDirection[t] = Point2D<Real>(-temp[1], temp[0]);
		}
		principalDirection[t] /= (Real)sqrt( Point2D<Real>::Dot( principalDirection[t] , g*principalDirection[t] ) );
	}
	return 1;
}

template< typename Real >
int InitializeEmbeddingMetric( const TexturedMesh& mesh , bool normalizeArea , std::vector< SquareMatrix< Real , 2 > >& embeddingMetric )
{
	embeddingMetric.resize(mesh.triangles.size());

	Real totalMass = 0;

	for (int t = 0; t < mesh.triangles.size(); t++) {

		Point3D< Real > vPos[3];
		for (int i = 0; i < 3; i++) vPos[i] = mesh.vertices[mesh.triangles[t][i]];

		SquareMatrix< Real , 2> g;
		Point3D< Real > dv[2] = { vPos[1] - vPos[0], vPos[2] - vPos[0] };
		for (int k = 0; k < 2; k++)for (int l = 0; l < 2; l++)g(k, l) = Point3D< Real >::Dot(dv[k], dv[l]);

		totalMass += sqrt(g.determinant()) / 2;

		embeddingMetric[t] = g;
	}

	if (normalizeArea) {
		for (int t = 0; t < mesh.triangles.size(); t++) {
			embeddingMetric[t] /= totalMass;
		}
	}

	return 1;
}

template< typename Real >
int InitializeUniformMetric( const TexturedMesh& mesh , bool normalizeArea , std::vector< SquareMatrix< Real , 2 > >& embeddingMetric )
{
	embeddingMetric.resize(mesh.triangles.size());

	Real totalMass = 0;

	for (int t = 0; t < mesh.triangles.size(); t++) {

		Point2D< Real > tPos[3];
		for (int i = 0; i < 3; i++) tPos[i] = mesh.textureCoordinates[3 * t + i];

		SquareMatrix< Real , 2 > g;
		Point2D< Real > dt[2] = { tPos[1] - tPos[0], tPos[2] - tPos[0] };
		for (int k = 0; k < 2; k++)for (int l = 0; l < 2; l++)g(k, l) = Point2D< Real >::Dot(dt[k], dt[l]);

		totalMass += sqrt(g.determinant()) / 2;

		embeddingMetric[t] = g;
	}

	if (normalizeArea) {
		for (int t = 0; t < mesh.triangles.size(); t++) {
			embeddingMetric[t] /= totalMass;
		}
	}

	return 1;
}

template< typename Real >
int InitializeParameterMetric( const TexturedMesh& mesh , const std::vector< SquareMatrix< Real , 2 > >& embeddingMetric , const std::vector< AtlasChart >& atlasCharts , std::vector< std::vector< SquareMatrix< Real , 2 > > >& parameterMetric )
{
	parameterMetric.resize(atlasCharts.size());
	for (int i = 0; i < atlasCharts.size(); i++) {
		parameterMetric[i].resize(atlasCharts[i].meshTriangleIndices.size());
		for (int k = 0; k < atlasCharts[i].meshTriangleIndices.size(); k++) {
			int t = atlasCharts[i].meshTriangleIndices[k];
			Point3D< Real > vPos[3];
			for (int i = 0; i < 3; i++) vPos[i] = mesh.vertices[mesh.triangles[t][i]];

			SquareMatrix< Real , 2 > embedding_metric = embeddingMetric[t];

			Point2D< Real > tPos[3];
			for (int i = 0; i < 3; i++) tPos[i] = mesh.textureCoordinates[3 * t + i];

			Point2D< Real > dp[2] = { tPos[1] - tPos[0], tPos[2] - tPos[0] };

			//Parametric map
			SquareMatrix<Real, 2> parametric_map_differential;
			parametric_map_differential(0, 0) = dp[0][0];
			parametric_map_differential(0, 1) = dp[0][1];

			parametric_map_differential(1, 0) = dp[1][0];
			parametric_map_differential(1, 1) = dp[1][1];

			SquareMatrix< Real , 2 > inverse_parametric_map_differential = parametric_map_differential.inverse();
			SquareMatrix< Real , 2 > parameter_metric = inverse_parametric_map_differential.transpose() * embedding_metric * inverse_parametric_map_differential;
			parameterMetric[i][k] = parameter_metric;
		}
	}
	return 1;
}

template< typename Real >
int InitializeMetric( TexturedMesh& mesh , const int metricMode , const std::vector< AtlasChart >& atlasCharts , std::vector< std::vector< SquareMatrix< Real , 2 > > >& parameterMetric )
{
	std::vector< SquareMatrix< Real , 2 > > surfaceMetric;
	std::vector< SquareMatrix< Real , 2 > > embeddingMetric;

	InitializeEmbeddingMetric( mesh , true , embeddingMetric );

	if     ( metricMode==EMBEDDING_METRIC ) surfaceMetric = embeddingMetric;
	else if( metricMode==UNIFORM_METRIC ) InitializeUniformMetric( mesh , true , surfaceMetric );
	else{ fprintf( stderr , "[ERROR] Unrecognized  metric!\n") ; return 0; }
	InitializeParameterMetric( mesh , surfaceMetric , atlasCharts , parameterMetric );
	return 1;
}

template< typename Real >
int InitializeAnisotropicMetric( TexturedMesh & mesh , const std::vector< AtlasChart > & atlasCharts , std::vector< std::vector< SquareMatrix< Real , 2 > > >& parameterMetric , std::vector< Point2D< Real > >& vf ,  bool useSmallestCurvarture )
{
	std::vector< SquareMatrix< Real , 2 > > surfaceMetric;
	std::vector< SquareMatrix< Real , 2 > > embeddingMetric;

	InitializeEmbeddingMetric( mesh , true , embeddingMetric );

	if( !vf.size() )
	{
		UpdateNormals( mesh );
		Eigen::SparseMatrix< Real > meshMassMatrix;
		Eigen::SparseMatrix< Real > meshStiffnessMatrix;
		InitializeMeshMatrices( mesh , meshMassMatrix , meshStiffnessMatrix );
		Eigen::SimplicialLDLT< Eigen::SparseMatrix< Real > > meshSolver( meshMassMatrix + meshStiffnessMatrix*1e-4 );
		SmoothSignal( meshMassMatrix , meshSolver , mesh.normals , true );
		InitializePrincipalCurvatureDirection( mesh , vf , mesh.normals , useSmallestCurvarture );
	}
	InitializeVectorFieldMetric( embeddingMetric , vf , true , surfaceMetric );
	InitializeParameterMetric( mesh , surfaceMetric , atlasCharts , parameterMetric );
	return 1;
}