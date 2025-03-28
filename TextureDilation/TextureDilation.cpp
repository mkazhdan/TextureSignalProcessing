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
#include <Src/MeshIO.h>
#include <Src/TextureIO.h>
#include "tinyexr.h"

using namespace MishaK;

CmdLineParameterArray< std::string , 2 > Input( "in" );
CmdLineParameter< std::string > Output( "out" ) , OutputTexturePositions( "outP" );
CmdLineParameter< unsigned int > DilationRadius( "radius" , 0 );
CmdLineParameter< double > CollapseEpsilon( "collapse" , 0 );
CmdLineReadable Verbose( "verbose" );

CmdLineReadable* params[] =
{
	&Input , &Output , &OutputTexturePositions , &DilationRadius , &Verbose ,
	&CollapseEpsilon ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n", ex );
	printf( "\t --%s <input mesh, input texture>\n" , Input.name.c_str() );
	printf( "\t[--%s <output texture>]\n" , Output.name.c_str() );
	printf( "\t[--%s <output texture positions>]\n" , OutputTexturePositions.name.c_str() );
	printf( "\t[--%s <dilation radius>=%d]\n" , DilationRadius.name.c_str() , DilationRadius.value );
	printf( "\t[--%s <collapse epsilon>=%g]\n" , CollapseEpsilon.name.c_str() , CollapseEpsilon.value );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}


int main( int argc , char* argv[] )
{
	static const unsigned int K = 2;
	static const unsigned int Dim = 3;

	CmdLineParse( argc-1 , argv+1 , params );
	if( !Input.set ){ ShowUsage( argv[0] ) ; return EXIT_FAILURE; }

	Miscellany::Timer timer;
	std::vector< Point< double , Dim > > vertices;
	std::vector< Point< double , K > > textureCoordinates;
	std::vector< SimplexIndex< K > > simplices;

	ReadTexturedMesh( Input.values[0] , vertices , textureCoordinates , simplices );
	if( CollapseEpsilon.value>0 ) CollapseVertices( vertices , simplices , CollapseEpsilon.value );
	RegularGrid< K , Point< double , 3 > > texture = ReadTexture( Input.values[1] );

	if( Verbose.set ) std::cout << "Vertices / texture coordinates / simplices: " << vertices.size() << " / " << textureCoordinates.size() << " / " << simplices.size() << std::endl;
	if( Verbose.set ) std::cout << "Texture resolution: " << texture.res(0) << " x " << texture.res(1) << std::endl;

	timer.reset();
	RegularGrid< 2 , Texels::TexelInfo > inputTexelInfo = Texels::GetSupportedTexelInfo< Dim , false , true >( simplices.size() , [&]( size_t v ){ return vertices[v]; } , [&]( size_t s ){ return simplices[s]; } , [&]( size_t s , unsigned int k ){ return textureCoordinates[ s*(K+1)+k ]; } , texture.res(0) , texture.res(1) , 0 , true , false );
	if( Verbose.set ) std::cout << "Got input texels: " << timer() << std::endl;

	timer.reset();
	RegularGrid< 2 , Texels::TexelInfo > dilatedTexelInfo = Texels::GetSupportedTexelInfo< Dim , false , true >( simplices.size() , [&]( size_t v ){ return vertices[v]; } , [&]( size_t s ){ return simplices[s]; } , [&]( size_t s , unsigned int k ){ return textureCoordinates[ s*(K+1)+k ]; } , texture.res(0) , texture.res(1) , DilationRadius.value , true , Verbose.set );
	if( Verbose.set ) std::cout << "Got dilated texels: " << timer() << std::endl;

	auto GetSimplex = [&]( unsigned int si )
		{
			Simplex< double , K , K > s;
			for( unsigned int k=0 ; k<=K ; k++ ) s[k] = textureCoordinates[ si*(K+1) + k ];
			for( unsigned int k=0 ; k<=K ; k++ ) for( unsigned int d=0 ; d<K ; d++ ) s[k][d] *= texture.res(d);
			return s;
		};

	auto SamplePosition = [&]( Texels::TexelInfo ti ){ return GetSimplex( ti.sIdx )( ti.bc ); };
	auto SampleValue = [&]( Texels::TexelInfo ti ){ return texture( SamplePosition( ti ) ); };

	timer.reset();
	ThreadPool::ParallelFor
		(
			0 , dilatedTexelInfo.resolution() ,
			[&]( unsigned int , size_t i ){ if( dilatedTexelInfo[i].sIdx!=-1 && inputTexelInfo[i].sIdx==-1 ) texture[i] = SampleValue( dilatedTexelInfo[i] ); }
		);
	if( Verbose.set ) std::cout << "Set dilated texels: " << timer() << std::endl;

	RegularGrid< K , Point< float , Dim > > texturePositions = Texels::GetTexelPositions< float , Dim >( simplices.size() , [&]( size_t v ){ return vertices[v]; } , [&]( size_t s ){ return simplices[s]; } , dilatedTexelInfo );

	if( Output.set ) WriteTexture( Output.value , texture );
	if( OutputTexturePositions.set )
		if( SaveEXR( &texturePositions[0][0] , (int)texture.res(0) , (int)texture.res(1) , Dim , 0, OutputTexturePositions.value.c_str() , nullptr )!=TINYEXR_SUCCESS )
			MK_ERROR_OUT( "Failed to save EXR file: " , OutputTexturePositions.value );

	return EXIT_SUCCESS;
}
