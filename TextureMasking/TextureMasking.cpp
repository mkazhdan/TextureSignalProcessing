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


#include <Misha/CmdLineParser.h> 
#include <Misha/Miscellany.h>
#include <Misha/Exceptions.h>
#include <Misha/FEM.h>
#include <Misha/MultiThreading.h>
#include <Misha/Texels.h>
#include <Src/TextureIO.h>
#include <Src/MeshIO.h>

using namespace MishaK;


CmdLineParameter< std::string >
Input( "in" ) ,
Output( "out" );

CmdLineParameterArray< unsigned int , 2 >
Resolution( "res" );

CmdLineParameter< double >
CollapseEpsilon( "collapse" , 0 );

CmdLineParameter< unsigned int >
DilationRadius( "radius" , 0 );

CmdLineReadable
UseNearest( "nearest" ) ,
NodeAtCorner( "nodeAtCorner" ) ,
BoundaryOnly( "boundary" ) ,
ID( "id" ) ,
Verbose( "verbose" );

CmdLineReadable* params[] =
{
	&Input ,
	&Output ,
	&Resolution ,
	&ID ,
	&UseNearest ,
	&NodeAtCorner ,
	&CollapseEpsilon ,
	&BoundaryOnly ,
	&DilationRadius ,
	&Verbose ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n", ex );
	printf( "\t --%s <input mesh>\n" , Input.name.c_str() );
	printf( "\t --%s <texture width, texture height> \n" , Resolution.name.c_str() );
	printf( "\t[--%s <output texture>]\n" , Output.name.c_str() );
	printf( "\t[--%s <collapse epsilon>=%g]\n" , CollapseEpsilon.name.c_str() , CollapseEpsilon.value );
	printf( "\t[--%s <dilation radius>]\n" , DilationRadius.name.c_str() );
	printf( "\t[--%s]\n" , UseNearest.name.c_str() );
	printf( "\t[--%s]\n" , NodeAtCorner.name.c_str() );
	printf( "\t[--%s]\n" , BoundaryOnly.name.c_str() );
	printf( "\t[--%s]\n" , ID.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

static const unsigned int K = 2;
static const unsigned int Dim = 3;

template< bool Nearest , bool NodeAtCellCenter >
RegularGrid< K , Point< double , 3 > > Execute
( 
	const std::vector< Point< double , Dim > > & vertices ,
	const std::vector< Point< double , K > > & textureCoordinates ,
	const std::vector< SimplexIndex< K > > & simplices
)
{
	using Index = unsigned int;
	using TexelInfo = typename Texels< NodeAtCellCenter , Index >::template TexelInfo< K >;

	RegularGrid< K , Point< double , 3 > > mask( Resolution.values );

	auto TextureSimplexFunctor = [&]( size_t sIdx )
		{
			Simplex< double , K , K > simplex;
			for( unsigned int k=0 ; k<=K ; k++ ) simplex[k] = textureCoordinates[ sIdx*(K+1) + k ];
			return simplex;
		};

	Miscellany::Timer timer;
	if( ID.set )
	{
		auto RandomPoint = []( void )
			{
				Point< double , Dim > p;
				for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = Random< double >();
				return p;
			};

		for( size_t i=0 ; i<mask.size() ; i++ ) mask[i] = Point< double , 3 >( 0. , 0. , 0. );

		RegularGrid< K , TexelInfo > texelInfo = Texels< NodeAtCellCenter , Index >::GetNodeTexelInfo( simplices.size() , TextureSimplexFunctor , mask.res() , DilationRadius.value , Verbose.set );
		for( size_t i=0 ; i<texelInfo.size() ; i++ ) if( texelInfo[i].sIdx!=-1 )
		{
			srand( texelInfo[i].sIdx );
			mask[i] = RandomPoint();
		}
	}
	else
	{
		for( size_t i=0 ; i<mask.size() ; i++ ) mask[i] = Point< double , 3 >( 1. , 0. , 0. );

		if( BoundaryOnly.set )
		{
			RegularGrid< K , typename Texels< NodeAtCellCenter , Index >::template TexelInfo< 1 > > activeTexels;
			std::vector< Simplex< double , K , 1 > > edges( simplices.size()*(K+1) );
			for( unsigned int i=0 ; i<simplices.size() ; i++ ) for( unsigned int k=0 ; k<=K ; k++ )
			{
				edges[ i*(K+1) + k ][0] = textureCoordinates[ i*(K+1) + k ];
				edges[ i*(K+1) + k ][1] = textureCoordinates[ i*(K+1) + (k+1)%(K+1) ];
			}
			activeTexels = Texels< NodeAtCellCenter , Index >::template GetSupportedTexelInfo< Nearest , 1 >( edges.size() , [&]( size_t e ){ return edges[e]; } , mask.res() , DilationRadius.value , Verbose.set );
			for( size_t i=0 ; i<activeTexels.size() ; i++ ) if( activeTexels[i].sIdx!=-1 ) mask[i] = Point< double , 3 >( 0. , 0. , 1. );
		}
		else
		{
			RegularGrid< K , TexelInfo > activeTexels = Texels< NodeAtCellCenter , Index >::template GetSupportedTexelInfo< Nearest >( simplices.size() , TextureSimplexFunctor , mask.res() , DilationRadius.value , Verbose.set );
			for( size_t i=0 ; i<activeTexels.size() ; i++ ) if( activeTexels[i].sIdx!=-1 ) mask[i] = Point< double , 3 >( 0. , 0. , 1. );
		}
	}

	std::cout << "Time: " << timer() << std::endl;

	return mask;
}

int main( int argc , char* argv[] )
{
	CmdLineParse( argc-1 , argv+1 , params );
	if( !Input.set || !Resolution.set ){ ShowUsage( argv[0] ) ; return EXIT_FAILURE; }

	std::vector< Point< double , Dim > > vertices;
	std::vector< Point< double , K > > textureCoordinates;
	std::vector< SimplexIndex< K > > simplices;

	ReadTexturedMesh( Input.value , vertices , textureCoordinates , simplices );
	if( CollapseEpsilon.value>0 ) CollapseVertices( vertices , simplices , CollapseEpsilon.value );

	RegularGrid< K , Point< double , 3 > > mask;
	if( UseNearest.set )
		if( NodeAtCorner.set ) mask = Execute< true  , false >( vertices , textureCoordinates , simplices );
		else                   mask = Execute< true  , true  >( vertices , textureCoordinates , simplices );
	else
		if( NodeAtCorner.set ) mask = Execute< false , false >( vertices , textureCoordinates , simplices );
		else                   mask = Execute< false , true  >( vertices , textureCoordinates , simplices );

	if( Output.set ) WriteTexture( Output.value , mask );

	return EXIT_SUCCESS;
}
