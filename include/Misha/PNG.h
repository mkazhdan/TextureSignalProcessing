/*
Copyright (c) 2023, Michael Kazhdan
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

#include <stdio.h>
#include <vector>
#define NEW_ZLIB
#ifdef _WIN32
#include "PNG/png.h"
#ifdef NEW_ZLIB
#include "ZLIB/zlib.h"
#endif // NEW_ZLIB
#else // !_WIN32
#include <png.h>
#ifdef NEW_ZLIB
#include <zlib.h>
#endif // NEW_ZLIB
#endif // _WIN32

namespace MishaK
{
	template< unsigned int BitDepth=8 >
	struct PNGReader : public ImageReader< BitDepth >
	{
		using ChannelType = typename ImageChannel< BitDepth >::Type;

		PNGReader( std::string fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
		~PNGReader( void );
		unsigned int nextRow( ChannelType * row );
		static bool GetInfo( std::string fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
		static bool GetInfo( std::string fileName , unsigned int& width , unsigned int& height , unsigned int& channels , unsigned int &bitDepth );
	protected:
		png_structp _png_ptr;
		png_infop _info_ptr;
		png_infop _end_info ;
		FILE* _fp;
		unsigned int _currentRow;
	};

	template< unsigned int BitDepth=8 >
	struct PNGWriter : public ImageWriter< BitDepth >
	{
		using ChannelType = typename ImageChannel< BitDepth >::Type;

		PNGWriter( std::string fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int quality=100 );
		~PNGWriter( void );
		unsigned int nextRow( const ChannelType * row );
		unsigned int nextRows( const ChannelType * rows , unsigned int rowNum );
	protected:
		FILE* _fp;
		png_structp _png_ptr;
		png_infop _info_ptr;
		unsigned int _currentRow;
	};

#include "PNG.inl"
}
#endif //PNG_INCLUDED
