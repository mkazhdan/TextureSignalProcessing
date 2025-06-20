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



template< typename GeometryReal , typename LengthToAnisotropyFunctor >
void InitializeVectorFieldMetric
(
	const ExplicitIndexVector< AtlasMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > >& embeddingMetric ,
	const ExplicitIndexVector< AtlasMeshTriangleIndex , Point2D< GeometryReal > >& vf ,
	const LengthToAnisotropyFunctor &LengthToAnisotropy ,
	bool normalizeArea ,
	ExplicitIndexVector< AtlasMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > >& outputMetric
)
{
	unsigned int tCount = (unsigned int)embeddingMetric.size();

	outputMetric.resize( embeddingMetric.size() );
	GeometryReal totalMass = 0;

	for( unsigned int t=0 ; t<tCount ; t++ )
	{
		SquareMatrix< GeometryReal , 2 > g = embeddingMetric[ AtlasMeshTriangleIndex(t) ];
		SquareMatrix< GeometryReal , 2 > gOrtho;
		GeometryReal len2 = std::max< GeometryReal >( Point2D< GeometryReal >::Dot( vf[ AtlasMeshTriangleIndex(t) ] , embeddingMetric[ AtlasMeshTriangleIndex(t) ] * vf[ AtlasMeshTriangleIndex(t) ] ) , 0 );
		if( len2>0 )
		{
			// Construct an orthonormal frame with basis[0] pointing in the direction of vf
			Point2D< GeometryReal > basis[2];
			{
				basis[0] = vf[ AtlasMeshTriangleIndex(t) ];
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
		outputMetric[ AtlasMeshTriangleIndex(t) ] = gOrtho*aniso + g;

		totalMass += sqrt( outputMetric[ AtlasMeshTriangleIndex(t) ].determinant() ) / 2;
	}
	if( normalizeArea ) for( unsigned int t=0 ; t<tCount ; t++ ) outputMetric[ AtlasMeshTriangleIndex(t) ] /= totalMass;
}

#define NORMALIZE_SURFACE_EMBEDDING

template< typename GeometryReal >
struct PrincipalCurvature
{
	Point2D< GeometryReal > dirs[2];
	GeometryReal values[2];
};

template< typename GeometryReal >
void InitializePrincipalCurvatureDirection
(
	const TexturedTriangleMesh< GeometryReal > &mesh ,
	const std::vector< Point3D< GeometryReal > >& vNormals ,
	std::vector< PrincipalCurvature< GeometryReal > >& principalCurvatures
)
{
	principalCurvatures.resize( mesh.numTriangles() );

#ifdef NORMALIZE_SURFACE_EMBEDDING
	GeometryReal inputMass = mesh.surface.area();
	GeometryReal edgeScaling = (GeometryReal)( 1.0 / sqrt(inputMass) );
#endif // NORMALIZE_SURFACE_EMBEDDING
	for( int t=0 ; t<mesh.numTriangles() ; t++ )
	{
		Simplex< GeometryReal , 3 , 2 > s = mesh.surfaceTriangle(t);
		SquareMatrix< GeometryReal , 2 > g = s.metric();
#ifdef NORMALIZE_SURFACE_EMBEDDING
		Point3D< GeometryReal > dv[2] = { ( s[1]-s[0] )*edgeScaling , ( s[2]-s[0] )*edgeScaling };
#else // !NORMALIZE_SURFACE_EMBEDDING
		Point3D< GeometryReal > dv[2] = { ( s[1]-s[0] ), ( s[2]-s[0] ) };
#endif // NORMALIZE_SURFACE_EMBEDDING

		Point3D< GeometryReal > vNormal[3];
		for( int i=0 ; i<3 ; i++ ) vNormal[i] = vNormals[ mesh.surface.triangles[t][i] ];
		Point3D< GeometryReal > dn[2] = { vNormal[1] - vNormal[0] , vNormal[2] - vNormal[0] };
		SquareMatrix< GeometryReal , 2 > gg;
		for( int k=0 ; k<2 ; k++ ) for( int l=0 ; l<2 ; l++ ) gg(k,l) = Point3D< GeometryReal >::Dot( dn[k] , dv[l] );
		gg(0,1) = gg(1,0) = ( gg(0,1)+gg(1,0 ) )/2.0;

		SquareMatrix< GeometryReal , 2 > S = g.inverse()*gg;

		GeometryReal a = (GeometryReal)1.0;
		GeometryReal b = -S.trace();
		GeometryReal c = S.determinant();
		GeometryReal discriminant = (GeometryReal)( b*b - 4.0*a*c );

		if( discriminant<0 ) MK_THROW( "Negative discriminant: " , discriminant );

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
void InitializeEmbeddingMetric
(
	const TexturedTriangleMesh< GeometryReal > &mesh ,
	bool normalizeArea ,
	ExplicitIndexVector< AtlasMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > & embeddingMetric
)
{
	embeddingMetric.resize( mesh.numTriangles() );

	GeometryReal totalMass = 0;

	for( unsigned int t=0 ; t<mesh.numTriangles() ; t++ )
	{
		Simplex< GeometryReal , 3 , 2 > s = mesh.surfaceTriangle( t );
		totalMass += s.measure();
		embeddingMetric[ AtlasMeshTriangleIndex(t) ] = s.metric();
	}

	if( normalizeArea ) for( int t=0 ; t<mesh.numTriangles() ; t++ ) embeddingMetric[ AtlasMeshTriangleIndex(t) ] /= totalMass;
}

template< typename GeometryReal >
void InitializeUniformMetric
(
	const TexturedTriangleMesh< GeometryReal > &mesh ,
	bool normalizeArea ,
	ExplicitIndexVector< AtlasMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > &embeddingMetric
)
{
	embeddingMetric.resize( mesh.numTriangles() );

	GeometryReal totalMass = 0;

	for( unsigned int t=0 ; t<mesh.numTriangles() ; t++ )
	{
		Simplex< GeometryReal , 2 , 2 > s = mesh.textureTriangle( t );
		SquareMatrix< GeometryReal , 2 > g = s.metric();
		totalMass += s.measure();
		embeddingMetric[ AtlasMeshTriangleIndex(t) ] = g;
	}

	if( normalizeArea ) for( unsigned int t=0 ; t<mesh.numTriangles() ; t++ ) embeddingMetric[ AtlasMeshTriangleIndex(t) ] /= totalMass;
}

template< typename GeometryReal >
void InitializeParameterMetric
(
	const TexturedTriangleMesh< GeometryReal > &mesh ,
	const ExplicitIndexVector< AtlasMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > &embeddingMetric ,
	const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > &atlasCharts ,
	ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > &parameterMetrics
)
{
	parameterMetrics.resize( atlasCharts.size() );
	for( unsigned int i=0 ; i<atlasCharts.size() ; i++ )
	{
		const AtlasChart< GeometryReal > & atlasChart = atlasCharts[ ChartIndex(i) ];
		ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > & parameterMetric = parameterMetrics[ ChartIndex(i) ];
		parameterMetric.resize( atlasChart.numTriangles() );
		for( unsigned int k=0 ; k<atlasChart.numTriangles() ; k++ )
		{
			AtlasMeshTriangleIndex t = atlasChart.atlasTriangle( ChartMeshTriangleIndex(k) );

			SquareMatrix< GeometryReal , 2 > embedding_metric = embeddingMetric[t];
			Simplex< GeometryReal , 2 , 2 > simplex = mesh.textureTriangle( static_cast< unsigned int >(t) );
			Point2D< GeometryReal > dp[2] = { simplex[1]-simplex[0] , simplex[2]-simplex[0] };

			//Parametric map
			SquareMatrix< GeometryReal , 2 > parametric_map_differential;
			parametric_map_differential(0,0) = dp[0][0];
			parametric_map_differential(0,1) = dp[0][1];
			parametric_map_differential(1,0) = dp[1][0];
			parametric_map_differential(1,1) = dp[1][1];

			SquareMatrix< GeometryReal , 2 > inverse_parametric_map_differential = parametric_map_differential.inverse();
			SquareMatrix< GeometryReal , 2 > parameter_metric = inverse_parametric_map_differential.transpose() * embedding_metric * inverse_parametric_map_differential;
			parameterMetric[ ChartMeshTriangleIndex(k) ] = parameter_metric;
		}
	}
}

template< typename GeometryReal >
void InitializeMetric
(
	TexturedTriangleMesh< GeometryReal > &mesh ,
	int metricMode ,
	const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > &atlasCharts ,
	ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > &parameterMetric ,
	bool normalize = true
)
{
	ExplicitIndexVector< AtlasMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > surfaceMetric;
	ExplicitIndexVector< AtlasMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > embeddingMetric;

	InitializeEmbeddingMetric( mesh , normalize , embeddingMetric );

	if     ( metricMode==EMBEDDING_METRIC ) surfaceMetric = embeddingMetric;
	else if( metricMode==UNIFORM_METRIC   ) InitializeUniformMetric( mesh , normalize , surfaceMetric );
	else MK_THROW( "Unrecognized  metric: " , metricMode );
	InitializeParameterMetric( mesh , surfaceMetric , atlasCharts , parameterMetric );
}

template< typename GeometryReal , typename LengthToAnisotropyFunctor >
void InitializeAnisotropicMetric
(
	TexturedTriangleMesh< GeometryReal > &mesh ,
	const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > &atlasCharts ,
	const ExplicitIndexVector< AtlasMeshTriangleIndex , Point2D< GeometryReal > > &vf ,
	const LengthToAnisotropyFunctor &LengthToAnisotropy ,
	ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > &parameterMetric
)
{
	ExplicitIndexVector< AtlasMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > surfaceMetric;
	ExplicitIndexVector< AtlasMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > embeddingMetric;

	InitializeEmbeddingMetric( mesh , true , embeddingMetric );
	InitializeVectorFieldMetric( embeddingMetric , vf , LengthToAnisotropy , true , surfaceMetric );
	InitializeParameterMetric( mesh , surfaceMetric , atlasCharts , parameterMetric );
}