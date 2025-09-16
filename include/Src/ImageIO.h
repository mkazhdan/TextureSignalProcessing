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

#include <ZLIB/zlib.h>
#define TINYEXR_USE_MINIZ 0
#define TINYEXR_IMPLEMENTATION 
#include "tinyexr.h"

#include <Misha/Image.h>
#include <Misha/RegularGrid.h>
#include <Misha/Miscellany.h>
#include <Misha/Exceptions.h>

namespace MishaK
{
	namespace TSP
	{
		template< unsigned int BitDepth , typename Data >
		void WriteImage( const RegularGrid< 2 , Data > &image , std::string fileName )
		{
			using CType = typename ImageChannel< BitDepth >::Type;
			static const CType Scale = ~((CType)0);

			std::string ext = ToLower( GetFileExtension( fileName ) );
			if( ext==std::string( "exr" ) )
			{
				static const bool IsPoint = !( std::is_same_v< Data , double > || std::is_same_v< Data , float > || std::is_same_v< Data , CType > );

				unsigned int channels;
				bool isFloat = true;
				float * values = nullptr;
				if constexpr( std::is_same_v< Data , Point< double , 4 > > || std::is_same_v< Data , Point< float , 4 > > || std::is_same_v< Data , Point< CType , 4 > > )
				{
					isFloat = !std::is_same_v< Data , Point< CType , 4 > >;
					channels = 4;
				}
				else if constexpr( std::is_same_v< Data , Point< double , 3 > > || std::is_same_v< Data , Point< float , 3 > > || std::is_same_v< Data , Point< CType , 3 > > )
				{
					isFloat = !std::is_same_v< Data , Point< CType , 3 > >;
					channels = 3;
				}
				else if constexpr( std::is_same_v< Data , Point< double , 1 > > || std::is_same_v< Data , Point< float , 1 > > || std::is_same_v< Data , Point< CType , 1 > > )
				{
					isFloat = !std::is_same_v< Data , Point< CType , 1 > >;
					channels = 1;
				}
				else if constexpr( std::is_same_v< Data , double > || std::is_same_v< Data , float > || std::is_same_v< Data , CType > )
				{
					isFloat = !std::is_same_v< Data , CType >;
					channels = 1;
				}
				else MK_THROW( "Only 1-, 3-, and, 4-channel images supported" ); 
				values = new float[ image.res(0) * image.res(1) * channels ];
				if( isFloat )
				{
					if constexpr( IsPoint ) for( unsigned int j=0 , idx=0 ; j<image.res(1) ; j++ ) for( unsigned int i=0 ; i<image.res(0) ; i++ ) for( unsigned int c=0 ; c<channels ; c++ , idx++ ) values[idx] = static_cast< float >( image(i,j)[c] );
					else                    for( unsigned int j=0 , idx=0 ; j<image.res(1) ; j++ ) for( unsigned int i=0 ; i<image.res(0) ; i++ , idx++ ) values[idx] = static_cast< float >( image(i,j) );
				}
				else
				{
					if constexpr( IsPoint ) for( unsigned int j=0 , idx=0 ; j<image.res(1) ; j++ ) for( unsigned int i=0 ; i<image.res(0) ; i++ ) for( unsigned int c=0 ; c<channels ; c++ , idx++ ) values[idx] = static_cast< float >( image(i,j)[c] ) / static_cast< float >( Scale );
					else                    for( unsigned int j=0 , idx=0 ; j<image.res(1) ; j++ ) for( unsigned int i=0 ; i<image.res(0) ; i++ , idx++ ) values[idx] = static_cast< float >( image(i,j) ) / static_cast< float >( Scale );
				}

				if( SaveEXR( values , (int)image.res(0) , (int)image.res(1) , channels , 0 , fileName.c_str() , nullptr )!=TINYEXR_SUCCESS ) MK_THROW( "Failed to save EXR file: " , fileName );

				delete[] values;

				return;
			}

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
				CType * pixels = new CType[ image.res(0)*image.res(1)*3 ];
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
		void ReadImage( RegularGrid< 2 , Data > &image , std::string fileName )
		{
			using CType = typename ImageChannel< BitDepth >::Type;
			static const CType Scale = ~((CType)0);

			std::string ext = ToLower( GetFileExtension( fileName ) );
			if( ext==std::string( "exr" ) )
			{
				const char *err = nullptr;
				int width , height;
				float * rgba;
				int ret;
				ret = LoadEXR( &rgba , &width , &height , fileName.c_str() , &err );
				if( ret!=TINYEXR_SUCCESS )
				{
					FreeEXRErrorMessage( err );
					MK_THROW( "Failed to load exr: ", ret );
				}
				image.resize( width , height );


				if constexpr( std::is_same_v< Data , Point3D< double > > || std::is_same_v< Data , Point3D< float > > )
				{
					for( int i=0 ; i<width ; i++ ) for( int j=0 ; j<height ; j++ ) for( int c=0 ; c<3 ; c++ )
						image(i,j)[c] = rgba[ (j*width+i)*4+c ];
				}
				else if constexpr( std::is_same_v< Data , Point3D< CType > > )
				{
					for( int i=0 ; i<width ; i++ ) for( int j=0 ; j<height ; j++ ) for( int c=0 ; c<3 ; c++ )
						image(i,j)[c] = static_cast< CType >( rgba[ (j*width+i)*4+c ] * static_cast< float >( Scale ) );
				}

				free( rgba );
				return;
			}

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

		void GetImageInfo( std::string fileName , unsigned int & width , unsigned int & height , unsigned int & channels , unsigned int & bitDepth )
		{
			std::string ext = ToLower( GetFileExtension( fileName ) );
			if( ext==std::string( "exr" ) )
			{
				const char *err = nullptr;
				int _width , _height;
				float * rgba;
				int ret;
				ret = LoadEXR( &rgba , &_width , &_height , fileName.c_str() , &err );
				free( rgba );
				width = static_cast< unsigned int >( _width );
				height = static_cast< unsigned int >( _height );
				channels = 4;
				bitDepth = 8;
			}
			else ImageReader< 8 >::GetInfo( fileName , width , height , channels , bitDepth );
		}
	}

}
#endif // IMAGE_IO_INCLUDED
