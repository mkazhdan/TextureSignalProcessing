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
#ifndef VECTOR_IO_INCLUDED
#define VECTOR_IO_INCLUDED
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <Misha/Image.h>
#include <Misha/Miscellany.h>

template< typename T >
void ReadVector( std::vector<T> & vec , const char * fileName )
{
	FILE * file;
	file = fopen( fileName , "rb" );
	if( !file ) Miscellany::Throw( "Unable to read %s" , fileName );
	int vecSize;
	fread( &vecSize , sizeof(int) , 1 , file );
	vec.resize( vecSize );
	fread( &vec[0] , sizeof(T) , vecSize , file );
	fclose( file );
}

template<typename T>
void WriteVector(const std::vector<T> & vec, const char * fileName){

	FILE * file;
	file = fopen(fileName, "wb");
	int vecSize = (int)vec.size();
	fwrite(&vecSize, sizeof(int), 1, file);
	fwrite(&vec[0], sizeof(T), vecSize, file);
	fclose(file);
}


template <class T>
void WriteBinaryImage(const Image<T> & image, const char * fileName) {
	FILE * file;
	file = fopen(fileName, "wb");
	int width = image.width();
	int height = image.height();
	fwrite(&width, sizeof(int), 1, file);
	fwrite(&height, sizeof(int), 1, file);
	fwrite(&image[0], sizeof(T), width*height, file);
	fclose(file);
}

template< class T >
void ReadBinaryImage( Image< T > &image , const char *fileName )
{
	FILE * file;
	file = fopen(fileName, "rb");
	if( !file ) Miscellany::Throw( "Unable to read %s" , fileName );
	int width, height;
	fread(&width, sizeof(int), 1, file);
	fread(&height, sizeof(int), 1, file);
	image.resize(width, height);
	fread(&image[0], sizeof(T), width*height, file);
	fclose(file);
}

#endif //VECTOR_IO_INCLUDED