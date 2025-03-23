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


CmdLineParameter< std::string > Input( "in" ) , Output( "out" );
CmdLineParameterArray< unsigned int , 2 > Resolution( "res" );
CmdLineReadable Cells( "cells" ) , ID( "id" );

CmdLineReadable* params[] =
{
	&Input , &Output , &Resolution , &Cells , &ID ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n", ex );
	printf( "\t --%s <input mesh>\n" , Input.name.c_str() );
	printf( "\t --%s <texture width, texture height> \n" , Resolution.name.c_str() );
	printf( "\t[--%s <output texture>]\n" , Output.name.c_str() );
	printf( "\t[--%s]\n" , Cells.name.c_str() );
	printf( "\t[--%s]\n" , ID.name.c_str() );
}

static const unsigned int K = 2;
static const unsigned int Dim = 3;

int main( int argc , char* argv[] )
{
	CmdLineParse( argc-1 , argv+1 , params );
	if( !Input.set || !Resolution.set ){ ShowUsage( argv[0] ) ; return EXIT_FAILURE; }

	std::vector< Point< double , Dim > > vertices;
	std::vector< Point< double , K > > textureCoordinates;
	std::vector< SimplexIndex< K > > simplices;

	ReadTexturedMesh( Input.value , vertices , textureCoordinates , simplices );

	RegularGrid< K , Point< double , 3 > > mask( Resolution.values );

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

		RegularGrid< K , Texels::TexelInfo > interiorTexels;
		if( !Cells.set ) interiorTexels = Texels::GetInteriorTexels( simplices.size() , [&]( size_t s  , unsigned int k ){ return textureCoordinates[ s*(K+1)+k ]; } , mask.res(0) , mask.res(1) );
		else             interiorTexels = Texels::GetActiveTexels< Dim >( simplices.size() , [&]( size_t v ){ return vertices[v]; } , [&]( size_t s ){ return simplices[s]; } , [&]( size_t s  , unsigned int k ){ return textureCoordinates[ s*(K+1)+k ]; } , mask.res(0) , mask.res(1) , 0 , false , false );
		for( size_t i=0 ; i<interiorTexels.size() ; i++ ) if( interiorTexels[i].sIdx!=-1 )
		{
			srand( interiorTexels[i].sIdx );
			mask[i] = RandomPoint();
		}
	}
	else
	{
		for( size_t i=0 ; i<mask.size() ; i++ ) mask[i] = Point< double , 3 >( 1. , 0. , 0. );

		RegularGrid< K , Texels::TexelInfo > activeTexels = Texels::GetActiveTexels< Dim >( simplices.size() , [&]( size_t v ){ return vertices[v]; } , [&]( size_t s ){ return simplices[s]; } , [&]( size_t s  , unsigned int k ){ return textureCoordinates[ s*(K+1)+k ]; } , mask.res(0) , mask.res(1) , 0 , false , false );
		for( size_t i=0 ; i<activeTexels.size() ; i++ ) if( activeTexels[i].sIdx!=-1 ) mask[i] = Point< double , 3 >( 0. , 1. , 0. );

		RegularGrid< K , Texels::TexelInfo > interiorTexels = Texels::GetInteriorTexels( simplices.size() , [&]( size_t s  , unsigned int k ){ return textureCoordinates[ s*(K+1)+k ]; } , mask.res(0) , mask.res(1) );
		for( size_t i=0 ; i<interiorTexels.size() ; i++ ) if( interiorTexels[i].sIdx!=-1 ) mask[i] = Point< double , 3 >( 0. , 0. , 1. );
	}

	std::cout << "Time: " << timer() << std::endl;

	if( Output.set ) WriteTexture( Output.value , mask );

	return EXIT_SUCCESS;
}
