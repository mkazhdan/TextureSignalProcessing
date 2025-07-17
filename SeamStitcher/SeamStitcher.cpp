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

#include <set>
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
	Mask( "mask" ) ,
	Out( "out" );

CmdLineParameter< unsigned int >
	QuadratureSamples( "qSamples" , 6 );

CmdLineReadable
	ExteriorBoundary( "exterior" ) ,
	Randomize( "random" ) ,
	Verbose( "verbose" );

CmdLineReadable* params[] =
{
	&In ,
	&Out ,
	&QuadratureSamples ,
	&ExteriorBoundary ,
	&Randomize ,
	&Verbose ,
	nullptr
};

void ShowUsage( const char *ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh and texels>\n" , In.name.c_str() );
	printf( "\t[--%s <input texel mask>]\n" , Mask.name.c_str() );
	printf( "\t[--%s <output texels>]\n" , Out.name.c_str() );
	printf( "\t[--%s <quadrature samples per triangle>=%d]\n" , QuadratureSamples.name.c_str() , QuadratureSamples.value );
	printf( "\t[--%s]\n" , ExteriorBoundary.name.c_str() );
	printf( "\t[--%s]\n" , Randomize.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

int main( int argc , char* argv[] )
{
	//
	// Letting x be the input texture and y the degrees of freedom associated with the unlocked variables,
	// our goal is to minimize the energy:
	//		E(y) = ( x + P * y )^t * S * ( x + p * y )
	//           = y^t * P^t * S *  P * y + 2 * y * P^t * S * x + x^t * S * x
	// Differentiang with respect to y and setting to 0 gives:
	//		0 = 2 * ( P^t * S * P ) * y + 2 * P^t * S * x
	// <=>	( P^t * S * P ) * y = - P^t * S * x
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
	if( Mask.set )
	{
		ReadImage< TextureBitDepth >( mask , Mask.value );
		if( mask.res(0)!=texture.res(0) || mask.res(1)!=texture.res(1) )
			MK_THROW( "Texture and mask resolutions don't match: " , texture.res(0) , " x " , texture.res(1) , " != " , mask.res(0) , " x " , mask.res(1) );
	}
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

	// Compute the set of degrees of freedom
	std::set< size_t > dofSet;
	{
		if( Mask.set )
		{
			auto IsDoF = []( Point< unsigned char , Channels > pixel )
				{
					for( unsigned int c=0 ; c<Channels ; c++ ) if( pixel[c] ) return false;
					return true;
				};
			for( size_t n=0 ; n<gd.numNodes() ; n++ )
			{
				std::pair< unsigned int , unsigned int > p = gd.node(n);
				if( IsDoF( mask( p.first , p.second ) ) ) dofSet.insert( n );
			}
		}
		else
		{
			if( ExteriorBoundary.set ) for( size_t n=0 ; n<gd.numNodes() ; n++ ) if( !gd.isCovered(n) ) dofSet.insert( n );
			else
			{
				for( size_t e=0 ; e<gd.numEdges() ; e++ ) if( gd.isChartCrossing(e) )
				{
					std::pair< size_t , size_t > edge = gd.edge( e );
					dofSet.insert( edge.first ) , dofSet.insert( edge.second );
				}
			}
		}
	}
	if( Verbose.set )
	{
		std::stringstream sStream;
		sStream << "DoFs " << dofSet.size();
		std::cout << pMeter( sStream.str() ) << std::endl;
	}

	// Compute the prolongation matrix from the degrees of freedom to texels
	Eigen::SparseMatrix< double > P , Pt;
	{
		std::vector< Eigen::Triplet< double > > triplets;
		{
			triplets.reserve( dofSet.size() );
			size_t idx = 0;
			for( auto iter=dofSet.begin() ; iter!=dofSet.end() ; iter++ )
				triplets.emplace_back( static_cast< typename Eigen::SparseMatrix< double >::StorageIndex >(*iter) , static_cast< typename Eigen::SparseMatrix< double >::StorageIndex >(idx++) , 1. );
		}
		P.resize( gd.numNodes() , triplets.size() );
		P.setFromTriplets( triplets.begin() , triplets.end() );
		Pt = P.transpose();
	}
	if( Verbose.set )
	{
		std::stringstream sStream;
		sStream << "Prolongation " << P.rows() << " x " << P.cols();
		std::cout << pMeter( sStream.str() ) << std::endl;
	}

	std::vector< Point< Real , Channels > > x( gd.numNodes() ) , b( gd.numNodes() );

	// Copy the texture values into the vector
	for( size_t n=0 ; n<gd.numNodes() ; n++ )
	{
		std::pair< unsigned int , unsigned int > coords = gd.node(n);
		x[n] = texture( coords.first , coords.second );
	}

	if( Randomize.set )
		for( size_t n=0 ; n<gd.numNodes() ; n++ ) if( dofSet.find(n)!=dofSet.end() )
			for( unsigned int c=0 ; c<Channels ; c++ ) x[n][c] = Random< double >() * 255.;

	// Construct the constraints
	gd.stiffness( &x[0] , &b[0] );
	if( Verbose.set ) std::cout << pMeter( "Coefficients" ) << std::endl;

	// Compute the system matrix
	Eigen::SparseMatrix< Real > M = Pt * gd.stiffness() * P;
	if( Verbose.set ) std::cout << pMeter( "Matrix" ) << std::endl;

	// Construct/factor the solver
	Solver< Real > solver( M );
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
	for( unsigned int c=0 ; c<Channels ; c++ )
	{
		Eigen::VectorXd _b( gd.numNodes() );
		for( size_t n=0 ; n<gd.numNodes() ; n++ ) _b[n] = b[n][c];
		_b = Pt * _b;
		Eigen::VectorXd _x = P * solver.solve( _b );
		for( size_t n=0 ; n<gd.numNodes() ; n++ ) x[n][c] -= _x[n];
	}
	if( Verbose.set ) std::cout << pMeter( "Solved" ) << std::endl;

	// Put the texel values back into the texture
	for( size_t n=0 ; n<gd.numNodes() ; n++ )
	{
		std::pair< unsigned int , unsigned int > coords = gd.node(n);
		texture( coords.first , coords.second ) = x[n];
	}

	if( Out.set ) WriteImage< TextureBitDepth >( texture , Out.value );

	return EXIT_SUCCESS;
}
