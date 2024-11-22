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

enum
{
	EMBEDDING_METRIC,
	UNIFORM_METRIC,
	UMBILIC_METRIC,
	EDGE_AWARE_METRIC,
};



template< typename GeometryReal , typename LengthToAnisotropyFunctor >
void InitializeVectorFieldMetric( const std::vector< SquareMatrix< GeometryReal , 2 > >& embeddingMetric , const std::vector< Point2D< GeometryReal > >& vf , const LengthToAnisotropyFunctor &LengthToAnisotropy , bool normalizeArea , std::vector< SquareMatrix< GeometryReal , 2 > >& outputMetric )
{
	int tCount = (int)embeddingMetric.size();

	outputMetric.resize( embeddingMetric.size() );
	GeometryReal totalMass = 0;

	for( int t=0 ; t<tCount ; t++ )
	{
		SquareMatrix< GeometryReal , 2 > g = embeddingMetric[t];
		SquareMatrix< GeometryReal , 2 > gOrtho;
		GeometryReal len2 = std::max< GeometryReal >( Point2D< GeometryReal >::Dot( vf[t] , embeddingMetric[t] * vf[t] ) , 0 );
		if( len2>0 )
		{
			// Construct an orthonormal frame with basis[0] pointing in the direction of vf
			Point2D< GeometryReal > basis[2];
			{
				basis[0] = vf[t];
				Point2D< GeometryReal > temp = g*basis[0];
				basis[1] = Point2D< GeometryReal >( -temp[1] , temp[0] );
				for( int e=0 ; e<2 ; e++ ) basis[e] /= sqrt( Point2D< GeometryReal >::Dot( basis[e] , g*basis[e] ) );
			}

			// Construct the (rotation) matrix taking coordinates in the orthonormal basis to coordinates in the standard triangle basis
			SquareMatrix< GeometryReal , 2 > R;
			for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) R(i,j) = basis[i][j];

			// Construct the orthogonally projected metric tensor
			{
				// Compute the transformation projecting onto the perpendicular
				SquareMatrix< GeometryReal , 2 > P;
				P(1,1) = 1;
				P = R * P * R.inverse();
				gOrtho = P.transpose() * g * P;
			}
		}
		else gOrtho *= 0;
		// When the anisotropy is small we want to revert to the standard metric
		// When the anisotropy is large we want to scale distances along the perpendicular direction
		GeometryReal aniso = (GeometryReal)LengthToAnisotropy( sqrt(len2) );
		outputMetric[t] = gOrtho*aniso + g;

		totalMass += sqrt( outputMetric[t].determinant() ) / 2;
	}
	if( normalizeArea ) for( int t=0 ; t<tCount ; t++ )	outputMetric[t] /= totalMass;
}

#define NORMALIZE_SURFACE_EMBEDDING

template< typename GeometryReal >
struct PrincipalCurvature
{
	Point2D< GeometryReal > dirs[2];
	GeometryReal values[2];
};

template< typename GeometryReal >
void InitializePrincipalCurvatureDirection( const OrientedTexturedTriangleMesh< GeometryReal > &mesh , const std::vector< Point3D< GeometryReal > >& vNormals , std::vector< PrincipalCurvature< GeometryReal > >& principalCurvatures )
{
	principalCurvatures.resize( mesh.triangles.size() );

#ifdef NORMALIZE_SURFACE_EMBEDDING
	GeometryReal inputMass = mesh.area();
	GeometryReal edgeScaling = (GeometryReal)( 1.0 / sqrt(inputMass) );
#endif // NORMALIZE_SURFACE_EMBEDDING
	for( int t=0 ; t<mesh.triangles.size() ; t++ )
	{
		Point3D< GeometryReal > vPos[3];
		for( int i=0 ; i<3 ; i++ ) vPos[i] = mesh.vertices[ mesh.triangles[t][i] ];

		SquareMatrix< GeometryReal , 2 > g;
#ifdef NORMALIZE_SURFACE_EMBEDDING
		Point3D< GeometryReal > dv[2] = { (vPos[1] - vPos[0])*edgeScaling, (vPos[2] - vPos[0])*edgeScaling };
#else // !NORMALIZE_SURFACE_EMBEDDING
		Point3D< GeometryReal > dv[2] = { (vPos[1] - vPos[0]), (vPos[2] - vPos[0]) };
#endif // NORMALIZE_SURFACE_EMBEDDING
		for( int k=0 ; k<2 ; k++ ) for( int l=0 ; l<2 ; l++ ) g(k,l) = Point3D< GeometryReal >::Dot( dv[k] , dv[l] );

		Point3D< GeometryReal > vNormal[3];
		for( int i=0 ; i<3 ; i++ ) vNormal[i] = vNormals[ mesh.triangles[t][i] ];
		Point3D< GeometryReal > dn[2] = { vNormal[1] - vNormal[0] , vNormal[2] - vNormal[0] };
		SquareMatrix< GeometryReal , 2 > gg;
		for( int k=0 ; k<2 ; k++ ) for( int l=0 ; l<2 ; l++ ) gg(k,l) = Point3D< GeometryReal >::Dot( dn[k] , dv[l] );
		gg(0,1) = gg(1,0) = ( gg(0,1)+gg(1,0 ) )/2.0;

		SquareMatrix< GeometryReal , 2 > S = g.inverse()*gg;

		GeometryReal a = (GeometryReal)1.0;
		GeometryReal b = -S.trace();
		GeometryReal c = S.determinant();
		GeometryReal discriminant = (GeometryReal)( b*b - 4.0*a*c );

		if( discriminant<0 ) Miscellany::Throw( "Negative discriminant" );

		discriminant = (GeometryReal)sqrt(discriminant);
		GeometryReal roots[] = { (-b-discriminant) / (2.0*a) , (-b+discriminant) / (2.0*a) };
		Point2D< GeometryReal > vectors[] = { Point2D< GeometryReal >(1,0) , Point2D< GeometryReal >(0,1) };
		if( S(1,0)!=0 || S(0,1)!=0 )
		{
			if( fabs( S(1,0) )>fabs( S(0,1) ) )	vectors[0] = Point2D< GeometryReal >( -S(1,0) , S(0,0) - roots[0] );
			else                                vectors[0] = Point2D< GeometryReal >( -S(1,1) + roots[0] , S(0,1) );
		}
		{
			Point2D< GeometryReal > temp = g*vectors[0];
			vectors[1] = Point2D< GeometryReal >( -temp[1] , temp[0] );
		}
		for( int i=0 ; i<2 ; i++ )
		{
			principalCurvatures[t].values[i] = roots[i];
			GeometryReal len2 = Point2D< GeometryReal >::Dot( vectors[i] , g*vectors[i] );
			if( len2>0 ) principalCurvatures[t].dirs[i] = vectors[i] / (GeometryReal)sqrt( len2 );
			else         principalCurvatures[t].dirs[i] = vectors[i] * 0;
		}
	}
}

template< typename GeometryReal >
void InitializeEmbeddingMetric( const OrientedTexturedTriangleMesh< GeometryReal > &mesh , bool normalizeArea , std::vector< SquareMatrix< GeometryReal , 2 > > &embeddingMetric )
{
	embeddingMetric.resize( mesh.triangles.size() );

	GeometryReal totalMass = 0;

	for( int t=0 ; t<mesh.triangles.size() ; t++ )
	{
		Point3D< GeometryReal > vPos[3];
		for( int i=0 ; i<3 ; i++ ) vPos[i] = mesh.vertices[ mesh.triangles[t][i] ];

		SquareMatrix< GeometryReal , 2 > g;
		Point3D< GeometryReal > dv[2] = { vPos[1] - vPos[0], vPos[2] - vPos[0] };
		for( int k=0 ; k<2 ; k++ ) for ( int l=0 ; l<2 ; l++ ) g(k,l) = Point3D< GeometryReal >::Dot( dv[k] , dv[l] );

		totalMass += sqrt( g.determinant() ) / 2;

		embeddingMetric[t] = g;
	}

	if( normalizeArea ) for( int t=0 ; t<mesh.triangles.size() ; t++ ) embeddingMetric[t] /= totalMass;
}

template< typename GeometryReal >
void InitializeUniformMetric( const OrientedTexturedTriangleMesh< GeometryReal > &mesh , bool normalizeArea , std::vector< SquareMatrix< GeometryReal , 2 > > &embeddingMetric )
{
	embeddingMetric.resize( mesh.triangles.size() );

	GeometryReal totalMass = 0;

	for( int t=0 ; t<mesh.triangles.size() ; t++ )
	{

		Point2D< GeometryReal > tPos[3];
		for( int i=0 ; i<3 ; i++ ) tPos[i] = mesh.textureCoordinates[ 3*t+i ];

		SquareMatrix< GeometryReal , 2 > g;
		Point2D< GeometryReal > dt[2] = { tPos[1] - tPos[0], tPos[2] - tPos[0] };
		for( int k=0 ; k<2 ; k++ ) for( int l=0 ; l<2 ; l++ ) g(k,l) = Point2D< GeometryReal >::Dot( dt[k] , dt[l] );

		totalMass += sqrt(g.determinant()) / 2;

		embeddingMetric[t] = g;
	}

	if( normalizeArea ) for( int t=0 ; t<mesh.triangles.size() ; t++ ) embeddingMetric[t] /= totalMass;
}

template< typename GeometryReal >
void InitializeParameterMetric( const OrientedTexturedTriangleMesh< GeometryReal > &mesh , const std::vector< SquareMatrix< GeometryReal , 2 > > &embeddingMetric , const std::vector< AtlasChart< GeometryReal > > &atlasCharts , std::vector< std::vector< SquareMatrix< GeometryReal , 2 > > > &parameterMetric )
{
	parameterMetric.resize( atlasCharts.size() );
	for( int i=0 ; i<atlasCharts.size() ; i++ )
	{
		parameterMetric[i].resize( atlasCharts[i].meshTriangleIndices.size() );
		for( int k=0 ; k<atlasCharts[i].meshTriangleIndices.size() ; k++ )
		{
			int t = atlasCharts[i].meshTriangleIndices[k];

			SquareMatrix< GeometryReal , 2 > embedding_metric = embeddingMetric[t];

			Point2D< GeometryReal > tPos[3];
			for( int ii=0 ; ii<3 ; ii++ ) tPos[ii] = mesh.textureCoordinates[3*t+ii];

			Point2D< GeometryReal > dp[2] = { tPos[1] - tPos[0], tPos[2] - tPos[0] };

			//Parametric map
			SquareMatrix< GeometryReal , 2 > parametric_map_differential;
			parametric_map_differential(0,0) = dp[0][0];
			parametric_map_differential(0,1) = dp[0][1];
			parametric_map_differential(1,0) = dp[1][0];
			parametric_map_differential(1,1) = dp[1][1];

			SquareMatrix< GeometryReal , 2 > inverse_parametric_map_differential = parametric_map_differential.inverse();
			SquareMatrix< GeometryReal , 2 > parameter_metric = inverse_parametric_map_differential.transpose() * embedding_metric * inverse_parametric_map_differential;
			parameterMetric[i][k] = parameter_metric;
		}
	}
}

template< typename GeometryReal >
void InitializeMetric( OrientedTexturedTriangleMesh< GeometryReal > &mesh , int metricMode , const std::vector< AtlasChart< GeometryReal > > &atlasCharts , std::vector< std::vector< SquareMatrix< GeometryReal , 2 > > > &parameterMetric )
{
	std::vector< SquareMatrix< GeometryReal , 2 > > surfaceMetric;
	std::vector< SquareMatrix< GeometryReal , 2 > > embeddingMetric;

	InitializeEmbeddingMetric( mesh , true , embeddingMetric );

	if     ( metricMode==EMBEDDING_METRIC ) surfaceMetric = embeddingMetric;
	else if( metricMode==UNIFORM_METRIC   ) InitializeUniformMetric( mesh , true , surfaceMetric );
	else Miscellany::Throw( "Unrecognized  metric: %d" , metricMode );
	InitializeParameterMetric( mesh , surfaceMetric , atlasCharts , parameterMetric );
}

template< typename GeometryReal , typename LengthToAnisotropyFunctor >
void InitializeAnisotropicMetric( OrientedTexturedTriangleMesh< GeometryReal > &mesh , const std::vector< AtlasChart< GeometryReal > > &atlasCharts , const std::vector< Point2D< GeometryReal > > &vf , const LengthToAnisotropyFunctor &LengthToAnisotropy , std::vector< std::vector< SquareMatrix< GeometryReal , 2 > > > &parameterMetric )
{
	std::vector< SquareMatrix< GeometryReal , 2 > > surfaceMetric;
	std::vector< SquareMatrix< GeometryReal , 2 > > embeddingMetric;

	InitializeEmbeddingMetric( mesh , true , embeddingMetric );
	InitializeVectorFieldMetric( embeddingMetric , vf , LengthToAnisotropy , true , surfaceMetric );
	InitializeParameterMetric( mesh , surfaceMetric , atlasCharts , parameterMetric );
}