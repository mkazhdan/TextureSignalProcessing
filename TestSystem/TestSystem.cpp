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

#include <Src/PreProcessing.h>

#include <Eigen/Sparse>

#include <Misha/CmdLineParser.h> 
#include <Misha/Miscellany.h>
#include <Misha/Exceptions.h>
#include <Misha/FEM.h>
#include <Misha/MultiThreading.h>
#include <Src/Hierarchy.h>
#include <Src/SimpleTriangleMesh.h>
#include <Src/Basis.h>
#include <Src/Operators.h>

using namespace MishaK;
using namespace MishaK::TSP;

CmdLineParameter< std::string >
	In( "in" );

CmdLineParameter< unsigned int >
	Width( "width" , 512 ) ,
	Height( "height" , 512 ) ,
	MatrixQuadrature( "mQuadrature" , 6 ) ,
	VectorQuadrature( "vQuadrature" , 6 );

CmdLineReadable
	SanityCheck( "sanityCheck" ) ,
	Approximate( "approximate" ) ,
	Parallel( "parallel" ) ,
	Single( "single" );

CmdLineReadable* params[] =
{
	&In ,
	&Width ,
	&Height ,
	&Approximate ,
	&Parallel,
	&Single ,
	&MatrixQuadrature ,
	&VectorQuadrature ,
	&SanityCheck ,
	NULL
};

void ShowUsage( const char *ex )
{
	printf( "Usage %s:\n" , ex );

	printf( "\t --%s <input mesh and texels>\n" , In.name.c_str() );
	printf( "\t[--%s <texture width>=%d]\n" , Width.name.c_str() , Width.value );
	printf( "\t[--%s <texture height>=%d]\n" , Height.name.c_str() , Height.value );
	printf( "\t[--%s <system matrix quadrature points per triangle>=%d]\n" , MatrixQuadrature.name.c_str(), MatrixQuadrature.value );
	printf( "\t[--%s <system vector quadrature points per triangle>=%d]\n" , VectorQuadrature.name.c_str(), VectorQuadrature.value );
	printf( "\t[--%s]\n" , Single.name.c_str() );
	printf( "\t[--%s]\n" , Approximate.name.c_str() );
	printf( "\t[--%s]\n" , Parallel.name.c_str() );
	printf( "\t[--%s]\n" , SanityCheck.name.c_str() );
}

template< typename PreReal , typename Real >
Eigen::SparseMatrix< Real > Execute( unsigned int textureWidth , unsigned int textureHeight , TexturedTriangleMesh< PreReal > mesh , const RegularGrid< 2 , Real > & conf )
{
	auto Conf = [&]( typename RegularGrid< 2 >::Index I ){ return conf(I); };

	HierarchicalSystem< PreReal , Real > hierarchy;
	std::vector< TextureNodeInfo< PreReal > > textureNodes;
	MassAndStiffnessOperators< Real > massAndStiffnessOperators;
	DivergenceOperator< Real > divergenceOperator;
	ScalarIntegrator< Real > scalarIntegrator;
	unsigned int levels = 1;

	Miscellany::PerformanceMeter pMeter( '.' );
	// Initialize the system
	{
		Miscellany::PerformanceMeter pMeter( '.' );

		ExplicitIndexVector< ChartIndex , AtlasChart< PreReal > > atlasCharts;

		MultigridBlockInfo multigridBlockInfo;
		InitializeHierarchy( mesh , textureWidth , textureHeight , levels , textureNodes , hierarchy , atlasCharts , multigridBlockInfo , SanityCheck.set );
		std::cout << pMeter( "Hierarchy" ) << std::endl;

		ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< PreReal , 2 > > > parameterMetric;
		InitializeMetric( mesh , EMBEDDING_METRIC , atlasCharts , parameterMetric );

		pMeter.reset();
		OperatorInitializer::Initialize( MatrixQuadrature.value , massAndStiffnessOperators , hierarchy.gridAtlases[0] , parameterMetric , atlasCharts , divergenceOperator , scalarIntegrator , VectorQuadrature.value , Approximate.set , SanityCheck.set , Conf );
		std::cout << pMeter( "Mass and stiffness" ) << std::endl;
	}
	std::cout << pMeter( "Initialized" ) << std::endl;
	printf( "Resolution: %d / %d x %d\n" , (int)textureNodes.size() , textureWidth , textureHeight );

	Eigen::SparseMatrix< Real > M = massAndStiffnessOperators.mass();
	{
		std::vector< Real > x( textureNodes.size() ) , b1( textureNodes.size() ) , b2( textureNodes.size() );
		Eigen::Matrix< Real , Eigen::Dynamic , 1 > _x( textureNodes.size() );
		for( unsigned int i=0 ; i<textureNodes.size() ; i++ ) x[i] = _x[i] = static_cast< Real >( 1.);

		massAndStiffnessOperators.mass( x , b1 );
		auto ScalarFunction = []( Real v , SquareMatrix< Real , 2 > ){ return v; };
		scalarIntegrator( x , ScalarFunction , b2 , Conf );
		Real area1 = 0 , area2 = 0;
		for( unsigned int i=0 ; i<x.size() ; i++ ) area1 += x[i] * b1[i] , area2 += x[i] * b2[i];
		std::cout << "Area: " << mesh.surface.area() << " : " << _x.dot( M * _x ) << " : " << area1 << " : " << area2 << std::endl;
	}
	return M;
}

template< typename PreReal , typename Real >
void Execute( unsigned int textureWidth , unsigned int textureHeight )
{
	TexturedTriangleMesh< PreReal > mesh;
	unsigned int levels = 1;

	mesh.read( In.value , false , 0. );
	{
		Real scale = static_cast< Real >( 1. / sqrt( mesh.surface.area() ) );
		for( unsigned int i=0 ; i<mesh.surface.vertices.size() ; i++ ) mesh.surface.vertices[i] *= scale;
	}

	RegularGrid< 2 , Real > confidences[2] , confidence;
	for( unsigned int c=0 ; c<2 ; c++ )
	{
		confidences[c].resize( textureWidth-1 , textureHeight-1 );
		for( size_t i=0 ; i<confidences[c].size() ; i++ ) confidences[c][i] = Random< Real >();
	}
	confidence.resize( textureWidth-1 , textureHeight-1 );
	for( size_t i=0 ; i<confidence.size() ; i++ ) confidence[i] = confidences[0][i] + confidences[1][i];

#if 1
	Eigen::SparseMatrix< Real > mass[3];
	for( unsigned int c=0 ; c<2 ; c++ ) mass[c] = Execute< PreReal , Real >( textureWidth , textureHeight , mesh , confidences[c] );
	mass[2] = Execute< PreReal , Real >( textureWidth , textureHeight , mesh , confidence );
#else
	for( unsigned int c=0 ; c<2 ; c++ )
	{
		auto Conf = [&]( typename RegularGrid< 2 >::Index I ){ return confidences[c](I); };

		HierarchicalSystem< PreReal , Real > hierarchy;
		std::vector< TextureNodeInfo< PreReal > > textureNodes;
		MassAndStiffnessOperators< Real > massAndStiffnessOperators;
		DivergenceOperator< Real > divergenceOperator;
		ScalarIntegrator< Real > scalarIntegrator;


		Miscellany::PerformanceMeter pMeter( '.' );
		// Initialize the system
		{
			Miscellany::PerformanceMeter pMeter( '.' );

			ExplicitIndexVector< ChartIndex , AtlasChart< PreReal > > atlasCharts;

			MultigridBlockInfo multigridBlockInfo;
			InitializeHierarchy( mesh , textureWidth , textureHeight , levels , textureNodes , hierarchy , atlasCharts , multigridBlockInfo , SanityCheck.set );
			std::cout << pMeter( "Hierarchy" ) << std::endl;

			ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< PreReal , 2 > > > parameterMetric;
			InitializeMetric( mesh , EMBEDDING_METRIC , atlasCharts , parameterMetric );

			pMeter.reset();
			OperatorInitializer::Initialize( MatrixQuadrature.value , massAndStiffnessOperators , hierarchy.gridAtlases[0] , parameterMetric , atlasCharts , divergenceOperator , scalarIntegrator , VectorQuadrature.value , Approximate.set , SanityCheck.set , Conf );
			std::cout << pMeter( "Mass and stiffness" ) << std::endl;
		}
		std::cout << pMeter( "Initialized" ) << std::endl;
		printf( "Resolution: %d / %d x %d\n" , (int)textureNodes.size() , textureWidth , textureHeight );

		{
			std::vector< Real > x( textureNodes.size() ) , b1( textureNodes.size() ) , b2( textureNodes.size() );
			Eigen::Matrix< Real , Eigen::Dynamic , 1 > _x( textureNodes.size() );
			for( unsigned int i=0 ; i<textureNodes.size() ; i++ ) x[i] = _x[i] = static_cast< Real >( 1.);
			Eigen::SparseMatrix< Real > M = massAndStiffnessOperators.mass();

			massAndStiffnessOperators.mass( x , b1 );
			auto ScalarFunction = []( Real v , SquareMatrix< Real , 2 > ){ return v; };
			scalarIntegrator( x , ScalarFunction , b2 , Conf );
			Real area1 = 0 , area2 = 0;
			for( unsigned int i=0 ; i<x.size() ; i++ ) area1 += x[i] * b1[i] , area2 += x[i] * b2[i];
			std::cout << "Area: " << mesh.surface.area() << " : " << _x.dot( M * _x ) << " : " << area1 << " : " << area2 << std::endl;
		}
	}
	{
	}
#endif
}

int main( int argc , char* argv[] )
{
	CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}

	if( !Parallel.set ) ThreadPool::ParallelizationType = ThreadPool::ParallelType::NONE;

	try
	{
		if( Single.set ) Execute< double , float  >( Width.value , Height.value );
		else             Execute< double , double >( Width.value , Height.value );
	}
	catch( Exception &e )
	{
		printf( "%s\n" , e.what() );
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
