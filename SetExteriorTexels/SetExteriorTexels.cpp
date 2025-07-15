/*
Copyright (c) 2025, Michael Kazhdan
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
#ifdef USE_EIGEN_PARDISO
#include <Eigen/PardisoSupport>
#endif // USE_EIGEN_PARDISO

#include <Misha/CmdLineParser.h> 
#include <Misha/Miscellany.h>
#include <Misha/Exceptions.h>
#include <Src/ImageIO.h>
#include <Src/GradientDomain.h>

using namespace MishaK;
using namespace MishaK::TSP;

#ifdef USE_EIGEN_PARDISO
template< typename Real > using Solver = Eigen::PardisoLDLT< Eigen::SparseMatrix< Real > >;
#else // !USE_EIGEN_PARDISO
template< typename Real > using Solver = Eigen::SimplicialLDLT< Eigen::SparseMatrix< Real > >;
#endif // USE_EIGEN_PARDISO

CmdLineParameterArray< std::string , 2 >
	In( "in" );

CmdLineParameter< std::string >
	Out( "out" );

CmdLineParameter< unsigned int >
	QuadratureSamples( "qSamples" , 6 );

CmdLineReadable
	Verbose( "verbose" );

CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&QuadratureSamples ,
	&Verbose ,
	nullptr
};

void ShowUsage( const char *ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh and texels>\n" , In.name.c_str() );
	printf( "\t[--%s <output texels>]\n" , Out.name.c_str() );
	printf( "\t[--%s <quadrature samples per triangle>=%d]\n" , QuadratureSamples.name.c_str() , QuadratureSamples.value );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

int main( int argc , char* argv[] )
{
	static const unsigned int Channels = 3;
	static const unsigned int TextureBitDepth = 8;

	using Real = double;

	CmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}

	TexturedTriangleMesh< Real > mesh;
	RegularGrid< 2 , Point< Real , Channels > > texture;
	RegularGrid< 2 , Point< unsigned char , Channels > > mask;

	// Read in the mesh and texture
	mesh.read( In.values[0] , false , 0  , false );
	ReadImage< TextureBitDepth >( texture , In.values[1] );
	if( Verbose.set )
	{
		std::cout << "Triangles: " << mesh.numTriangles() << std::endl;
		std::cout << "Resolution: " << texture.res(0) << " x " << texture.res(1) << std::endl;
	}

	Miscellany::PerformanceMeter pMeter( '.' );

	// Construct the gradient domain object from the textured mesh
	GradientDomain< Real > gd
	(
		QuadratureSamples.value ,
		mesh.numTriangles() ,
		mesh.surface.vertices.size() ,
		mesh.texture.vertices.size() ,
		[&]( size_t t , unsigned int k ){ return mesh.surface.triangles[t][k]; } ,
		[&]( size_t v ){ return mesh.surface.vertices [v]; } ,
		[&]( size_t t , unsigned int k ){ return mesh.texture.triangles[t][k]; } ,
		[&]( size_t v ){ return mesh.texture.vertices [v]; } ,
		texture.res(0) ,
		texture.res(1)
	);
	if( Verbose.set ) std::cout << pMeter( "Gradient domain" ) << std::endl;

	// Solving E(_x) = ( x + P*_x )^t * S * ( x + P*_x )
	//               = x^t * S * x + 2 * _x^t * R * S * x + _x^t * R * S * P * _x
	// Differentiating and setting to zero:
	//		       0 = 2 * R * S * x + 2 * R * S * P * _x
	//            _x = - ( R * S * P )^{-1} * ( R * S * x )

	// Constuct the prolongation matrix from exterior texels to texels
	Eigen::SparseMatrix< Real > P , R;
	{
		size_t exteriorCount = 0;
		for( size_t n=0 ; n<gd.numNodes() ; n++ ) if( !gd.isCovered(n) ) exteriorCount++;

		std::vector< Eigen::Triplet< Real > > triplets;
		triplets.reserve( exteriorCount );

		exteriorCount = 0;
		for( size_t n=0 ; n<gd.numNodes() ; n++ ) if( !gd.isCovered(n) ) triplets.emplace_back( static_cast< unsigned int >( n ) , static_cast< unsigned int >( exteriorCount++ ) , 1 );

		P.resize( gd.numNodes() , exteriorCount );
		P.setFromTriplets( triplets.begin() , triplets.end() );

		R = P.transpose();
	}
	if( Verbose.set ) std::cout << pMeter( "Prolongation" ) << std::endl;

	std::vector< Point< Real , Channels > > x( gd.numNodes() ) , b;

	// Copy the texture values into the vector
	for( size_t n=0 ; n<gd.numNodes() ; n++ )
	{
		std::pair< unsigned int , unsigned int > coords = gd.node(n);
		x[n] = texture( coords.first , coords.second );

		// Zero out the exterior texels (not really required)
		if( !gd.isCovered(n) ) x[n] = Point< Real , Channels >();
	}

	// Get the system matrices
	Eigen::SparseMatrix< Real > RS = R * gd.stiffness();
	Eigen::SparseMatrix< Real > RSP = RS * P;
	if( Verbose.set ) std::cout << pMeter( "Matrix" ) << std::endl;

	// Construct/factor the solver
	Solver< Real > solver( RSP );
	switch( solver.info() )
	{
	case Eigen::Success: break;
	case Eigen::NumericalIssue: MK_THROW( "Failed to factor matrix (numerical issue): "            , typeid( Solver< Real > ).name() ) ; break;
	case Eigen::NoConvergence:  MK_THROW( "Failed to factor matrix (no convergence): "             , typeid( Solver< Real > ).name() ) ; break;
	case Eigen::InvalidInput:   MK_THROW( "Failed to factor matrix (invalid input): "              , typeid( Solver< Real > ).name() ) ; break;
	default:                    MK_THROW( "Failed to factor matrix (info=" , solver.info() , "): " , typeid( Solver< Real > ).name() );
	}
	if( Verbose.set ) std::cout << pMeter( "Factored" ) << std::endl;

	// Solve the system per channel
	{
		Eigen::VectorXd _x( gd.numNodes() );
		for( unsigned int c=0 ; c<Channels ; c++ )
		{
			for( size_t n=0 ; n<gd.numNodes() ; n++ ) _x[n] = x[n][c];
			Eigen::VectorXd _b = - RS * _x;
			_x += P * solver.solve( _b );
			for( size_t n=0 ; n<gd.numNodes() ; n++ ) x[n][c] = _x[n];
		}
		if( Verbose.set ) std::cout << pMeter( "Solved" ) << std::endl;
	}

	// Put the texel values back into the texture
	for( size_t n=0 ; n<gd.numNodes() ; n++ )
	{
		std::pair< unsigned int , unsigned int > coords = gd.node(n);
		texture( coords.first , coords.second ) = x[n];
	}

	if( Out.set ) WriteImage< TextureBitDepth >( texture , Out.value );

	return EXIT_SUCCESS;
}
