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
#include <Misha/MultiIndex.h>
#include <Src/TextureIO.h>
#include <Src/MeshIO.h>

using namespace MishaK;
using namespace MishaK::TSP;

enum Rasterization
{
	ACTIVE ,
	BOUNDARY ,
	TRIANGLE_ID ,
	UNSIGNED_NODE_INCIDENCE_COUNT ,
	SIGNED_NODE_INCIDENCE_COUNT ,
	COUNT
};

std::string RasterizationNames[] =
{
	"active" ,
	"boundary" ,
	"triangle id" ,
	"node incidence count (unsigned)" ,
	"node incidence count (signed)"
};

CmdLineParameter< std::string >
	Input( "in" ) ,
	Output( "out" );

CmdLineParameterArray< unsigned int , 2 >
	Resolution( "res" );

CmdLineParameter< double >
	CollapseEpsilon( "collapse" , 0 );

CmdLineParameter< unsigned int >
	RasterizationType( "rasterize" , Rasterization::ACTIVE ) ,
	DilationRadius( "radius" , 0 );

CmdLineReadable
	UseNearest( "nearest" ) ,
	NodeAtCorner( "nodeAtCorner" ) ,
	Verbose( "verbose" );

CmdLineReadable* params[] =
{
	&Input ,
	&Output ,
	&Resolution ,
	&UseNearest ,
	&NodeAtCorner ,
	&CollapseEpsilon ,
	&RasterizationType ,
	&DilationRadius ,
	&Verbose ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n", ex );
	printf( "\t --%s <input mesh>\n" , Input.name.c_str() );
	printf( "\t --%s <texture width, texture height> \n" , Resolution.name.c_str() );
	printf( "\t[--%s <output texture image/grid>]\n" , Output.name.c_str() );
	printf( "\t[--%s <collapse epsilon>=%g]\n" , CollapseEpsilon.name.c_str() , CollapseEpsilon.value );
	printf( "\t[--%s <dilation radius>]\n" , DilationRadius.name.c_str() );
	printf( "\t[--%s <rasterization type>=%d]\n" , RasterizationType.name.c_str() , RasterizationType.value );
	for( unsigned int i=0 ; i<Rasterization::COUNT ; i++ ) printf( "\t\t%d] %s\n" , i , RasterizationNames[i].c_str() );
	printf( "\t[--%s]\n" , UseNearest.name.c_str() );
	printf( "\t[--%s]\n" , NodeAtCorner.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

static const unsigned int K = 2;
static const unsigned int Dim = 3;

template< bool Nearest , bool NodeAtCellCenter >
RegularGrid< K , int > Execute
( 
	const std::vector< Point< double , Dim > > & vertices ,
	const std::vector< Point< double , K > > & textureCoordinates ,
	const std::vector< SimplexIndex< K > > & simplices ,
	const std::vector< SimplexIndex< K > > & textureSimplices
)
{
	using Index = unsigned int;
	using TexelInfo = typename Texels< NodeAtCellCenter , Index >::template TexelInfo< K >;

	RegularGrid< K , int > mask( Resolution.values );
	for( size_t i=0 ; i<mask.size() ; i++ ) mask[i] = 0;

	auto TextureSimplexFunctor = [&]( size_t sIdx )
		{
			Simplex< double , K , K > simplex;
			for( unsigned int k=0 ; k<=K ; k++ ) simplex[k] = textureCoordinates[ textureSimplices[sIdx][k] ];
			return simplex;
		};

	Miscellany::Timer timer;
	switch( RasterizationType.value )
	{
	case Rasterization::ACTIVE:
	{
		RegularGrid< K , TexelInfo > activeTexels = Texels< NodeAtCellCenter , Index >::template GetSupportedTexelInfo< Nearest >( simplices.size() , TextureSimplexFunctor , mask.res() , DilationRadius.value , Verbose.set );
		for( size_t i=0 ; i<activeTexels.size() ; i++ ) if( activeTexels[i].sIdx!=-1 ) mask[i] = 1;
	}
	break;
	case Rasterization::TRIANGLE_ID:
	{
		RegularGrid< K , TexelInfo > texelInfo = Texels< NodeAtCellCenter , Index >::GetNodeTexelInfo( simplices.size() , TextureSimplexFunctor , mask.res() , DilationRadius.value , Verbose.set );
		for( size_t i=0 ; i<texelInfo.size() ; i++ ) mask[i] = texelInfo[i].sIdx;
	}
	break;
	case Rasterization::BOUNDARY:
	{
		RegularGrid< K , typename Texels< NodeAtCellCenter , Index >::template TexelInfo< 1 > > activeTexels;
		std::map< MultiIndex< 2 > , unsigned int > edgeIncidenceCount;
		for( unsigned int i=0 ; i<textureSimplices.size() ; i++ )
			textureSimplices[i].template processFaces< 1 >( [&]( SimplexIndex< 1 > e ){ edgeIncidenceCount[ MultiIndex< 2 >( &e[0] ) ]++; } );

		std::vector< Simplex< double , K , 1 > > edges;
		for( unsigned int i=0 ; i<textureSimplices.size() ; i++ )
		{
			auto Kernel = [&]( SimplexIndex< 1 > e )
				{
					if( edgeIncidenceCount[ MultiIndex< 2 >( &e[0] ) ]==1 ) edges.emplace_back( textureCoordinates[ e[0] ] , textureCoordinates[ e[1] ] );
				};
			textureSimplices[i].template processFaces< 1 >( Kernel );
		}
		activeTexels = Texels< NodeAtCellCenter , Index >::template GetSupportedTexelInfo< Nearest , 1 >( edges.size() , [&]( size_t e ){ return edges[e]; } , mask.res() , DilationRadius.value , Verbose.set );
		for( size_t i=0 ; i<activeTexels.size() ; i++ ) if( activeTexels[i].sIdx!=-1 ) mask[i] = 1;

	}
	break;
	case Rasterization::UNSIGNED_NODE_INCIDENCE_COUNT:
	{
		RegularGrid< K , std::vector< Index > > texelInfo = Texels< NodeAtCellCenter , Index >::GetNodeSimplexIndices( simplices.size() , TextureSimplexFunctor , mask.res() );
		for( size_t i=0 ; i<texelInfo.size() ; i++ ) mask[i] = static_cast< unsigned int >( texelInfo[i].size() );
	}
	break;
	case Rasterization::SIGNED_NODE_INCIDENCE_COUNT:
	{
		RegularGrid< K , std::vector< Index > > texelInfo = Texels< NodeAtCellCenter , Index >::GetNodeSimplexIndices( simplices.size() , TextureSimplexFunctor , mask.res() );
		for( size_t i=0 ; i<texelInfo.size() ; i++ )
		{
			unsigned int count = 0;
			for( unsigned int j=0 ; j<texelInfo[i].size() ; j++ )
				if( TextureSimplexFunctor( texelInfo[i][j] ).volume(true)<0 ) count--;
				else                                                          count++;
			mask[i] = count;
		}
	}
	break;
	default: MK_ERROR_OUT( "Unrecognized rasterization type: " , RasterizationType.value );
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
	std::vector< SimplexIndex< K > > simplices , textureSimplices;

	ReadTexturedMesh( Input.value , vertices , textureCoordinates , simplices , textureSimplices );
	if( CollapseEpsilon.value>0 ) CollapseVertices( vertices , simplices , CollapseEpsilon.value );

	RegularGrid< K , int > mask;
	if( UseNearest.set )
		if( NodeAtCorner.set ) mask = Execute< true  , false >( vertices , textureCoordinates , simplices , textureSimplices );
		else                   mask = Execute< true  , true  >( vertices , textureCoordinates , simplices , textureSimplices );
	else
		if( NodeAtCorner.set ) mask = Execute< false , false >( vertices , textureCoordinates , simplices , textureSimplices );
		else                   mask = Execute< false , true  >( vertices , textureCoordinates , simplices , textureSimplices );

	if( Output.set )
	{
		std::string ext = ToLower( GetFileExtension( Output.value ) );
		if( ext=="grid" ) mask.write( Output.value );
		else
		{
			RegularGrid< K , Point< double , 3 > > _mask( mask.res() );
			if( RasterizationType.value==Rasterization::UNSIGNED_NODE_INCIDENCE_COUNT || RasterizationType.value==Rasterization::SIGNED_NODE_INCIDENCE_COUNT )
			{
				unsigned int mn = std::numeric_limits< unsigned int >::infinity() , mx = 0;
				for( size_t i=0 ; i<mask.size() ; i++ ) mn = std::min< unsigned int >( mn , mask[i] ) , mx = std::max< unsigned int >( mx , mask[i] );
				for( size_t i=0 ; i<mask.size() ; i++ ) _mask[i] = Point< double , 3 >(1.,1.,1.) * static_cast< double >( mask[i] - mn ) / ( mx - mn );
			}
			else if( RasterizationType.value==Rasterization::TRIANGLE_ID )
			{
				auto RandomPoint = []( void )
					{
						Point< double , Dim > p;
						for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = Random< double >();
						return p;
					};
				for( size_t i=0 ; i<mask.size() ; i++ ) if( mask[i]!=-1 )
				{
					srand( mask[i] );
					_mask[i] = RandomPoint();
				}
			}
			else
			{
				for( size_t i=0 ; i<mask.size() ; i++ ) _mask[i] = mask[i]==0 ? Point< double , 3 >(1.,0.,.0) : Point< double , 3 >(0.,0.,1.);
			}
			WriteTexture( Output.value , _mask );
		}
	}

	return EXIT_SUCCESS;
}
