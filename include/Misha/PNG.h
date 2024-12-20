/*
Copyright (c) 2011, Michael Kazhdan and Ming Chuang
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
#ifndef PNG_INCLUDED
#define PNG_INCLUDED

#ifdef _WIN32
#include "PNG/png.h"
#else // !_WIN32
#include <png.h>
#endif // _WIN32


#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
template< unsigned int BitDepth=8 >
struct PNGReader : public ImageReader< BitDepth >
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
{
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
	using ChannelType = typename ImageChannel< BitDepth >::Type;
#endif // VARIABLE_SIZED_IMAGE_CHANNEL

	PNGReader( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
	~PNGReader( void );
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
	unsigned int nextRow( ChannelType * row );
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
	unsigned int nextRow( unsigned char* row );
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
	static bool GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
	static bool GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels , unsigned int &bitDepth );
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
protected:
	png_structp _png_ptr;
	png_infop _info_ptr;
	png_infop _end_info ;
	FILE* _fp;
	unsigned int _currentRow;
};

#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
template< unsigned int BitDepth=8 >
struct PNGWriter : public ImageWriter< BitDepth >
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
struct PNGWriter : public ImageWriter
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
{
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
	using ChannelType = typename ImageChannel< BitDepth >::Type;
#endif // VARIABLE_SIZED_IMAGE_CHANNEL

	PNGWriter( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int quality=100 );
	~PNGWriter( void );
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
	unsigned int nextRow( const ChannelType * row );
	unsigned int nextRows( const ChannelType * rows , unsigned int rowNum );
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
	unsigned int nextRow( const unsigned char* row );
	unsigned int nextRows( const unsigned char* rows , unsigned int rowNum );
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
protected:
	FILE* _fp;
	png_structp _png_ptr;
	png_infop _info_ptr;
	unsigned int _currentRow;
};

#include "PNG.inl"
#endif // PNG_INCLUDED
