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
#ifndef IMAGE_IO_INCLUDED
#define IMAGE_IO_INCLUDED

#include <Misha/Image.h>
#include <Misha/RegularGrid.h>
#include <Misha/Miscellany.h>
#include <Misha/Exceptions.h>

namespace MishaK
{
	template< typename Data > using Image = RegularGrid< 2 , Data >;

	template< unsigned int BitDepth , typename Data >
	void WriteImage( const Image< Data > &image , std::string fileName )
	{
		using CType = typename ImageChannel< BitDepth >::Type;
		static const CType Scale = ~((CType )0);
		if constexpr( std::is_same_v< Data , Point3D< double > > )
		{
			double scale = (double)Scale;
			CType * pixels = new CType[ image.res(0)*image.res(1)*3 ];
			for( unsigned int i=0 ; i<image.res(0) ; i++ ) for( unsigned int j=0 ; j<image.res(1) ; j++ ) for( int c=0 ; c<3 ; c++ ) pixels[ 3*(j*image.res(0)+i)+c ] = (CType)std::min< long long >( Scale , std::max< long long >( 0 , (long long)( image(i,j)[c] * scale + 0.5 ) ) );
			ImageWriter< BitDepth >::Write( fileName , pixels , image.res(0) , image.res(1) , 3 );
			delete[] pixels;
		}
		else if constexpr( std::is_same_v< Data , Point3D< float > > )
		{
			float scale = (float)Scale;
			CType * pixels = new CType[ image.res(0) * image.res(1) * 3 ];
			for( unsigned int i=0 ; i<image.res(0) ; i++ ) for( unsigned int j=0 ; j<image.res(1) ; j++ ) for( int c=0 ; c<3 ; c++ ) pixels[ 3*(j*image.res(0)+i)+c ] = (CType)std::min< long long >( Scale , std::max< long long >( 0 , (long long)( image(i,j)[c] * scale + 0.5f ) ) );
			ImageWriter< BitDepth >::Write( fileName , pixels , image.res(0) , image.res(1) , 3 );
			delete[] pixels;
		}
		else if constexpr( std::is_same_v< Data , Point3D< CType > > )
		{
			ImageWriter< BitDepth >::Write( fileName , (const CType*)image() , image.res(0) , image.res(1) , 3 );
		}
		else MK_THROW( "Bad data type " );
	}

	template< unsigned int BitDepth , typename Data >
	void ReadImage( Image< Data > &image , std::string fileName )
	{
		using CType = typename ImageChannel< BitDepth >::Type;
		static const CType Scale = ~((CType )0);
		unsigned int width , height;
		CType * pixels = ImageReader< BitDepth >::ReadColor( fileName , width , height );
		if( !pixels ) MK_THROW( "Failed to read image: " , fileName );
		image.resize( width , height );

		if constexpr( std::is_same_v< Data , Point3D< double > > )
		{
			double scale = (double)Scale;
			for( int i=0 ; i<(int)width ; i++ ) for( int j=0 ; j<(int)height ; j++ ) for( int c=0 ; c<3 ; c++ ) image(i,j)[c] = ( pixels[ (j*width+i)*3+c ] ) / scale;
		}
		else if constexpr( std::is_same_v< Data , Point3D< float > > )
		{
			float scale = (float)Scale;
			for( int i=0 ; i<(int)width ; i++ ) for( int j=0 ; j<(int)height ; j++ ) for( int c=0 ; c<3 ; c++ ) image(i,j)[c] = ( pixels[ (j*width+i)*3+c ] ) / scale;
		}
		else if constexpr( std::is_same_v< Data , Point3D< CType > > )
		{
			memcpy( (CType*)image() , pixels , sizeof( CType ) * image.res(0) * image.res(1) * 3 );
		}
		else MK_THROW( "Bad data type " );

		delete[] pixels;
	}
}
#endif // IMAGE_IO_INCLUDED
