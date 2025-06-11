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

#ifndef TEXTURE_INCLUDED
#define TEXTURE_INCLUDED

#include <string>
#include <Misha/Image.h>
#include <Misha/RegularGrid.h>
#include <Misha/Geometry.h>

namespace MishaK
{
	namespace TSP
	{
		RegularGrid< 2 , Point< double , 3 > > ReadTexture( std::string fileName )
		{
			unsigned int res[2];
			unsigned char *pixels = ImageReader<>::ReadColor( fileName , res[0] , res[1] );

			RegularGrid< 2 , Point< double , 3 > > img;
			img.resize( res );

			for( unsigned int j=0 , idx=0 ; j<img.res(1) ; j++ ) for( unsigned int i=0 ; i<img.res(1) ; i++ ) for( unsigned int k=0 ; k<3 ; k++ , idx++ )
				img(i,res[1]-1-j)[k] = pixels[idx]/255.;

			delete[] pixels;

			return img;
		};

		void WriteTexture( std::string fileName , const RegularGrid< 2 , Point< double , 3 > > &img )
		{
			unsigned char * pixels = new unsigned char[ img.resolution() * 3 ];
			for( unsigned int j=0 , idx=0 ; j<img.res(1) ; j++ ) for( unsigned int i=0 ; i<img.res(1) ; i++ ) for( unsigned int k=0 ; k<3 ; k++ , idx++ )
				pixels[idx] = static_cast< unsigned char >( std::min< double >( 255. , std::max< double >( 0. , img(i,img.res(1)-1-j)[k]*255. ) ) );

			ImageWriter<>::Write( fileName , pixels , img.res(0) , img.res(1) , 3 );

			delete[] pixels;
		}
	}
}
#endif // TEXTURE_INCLUDED