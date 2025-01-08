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
#include <stdio.h>
#include <vector>
#ifdef _WIN32
#include "PNG/png.h"
#else // !_WIN32
#include <png.h>
#endif // _WIN32
#include "Miscellany.h"
#ifdef NEW_CODE
#include "Exceptions.h"
#endif // NEW_CODE

template< unsigned int BitDepth >
inline PNGReader< BitDepth >::PNGReader( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	_currentRow = 0;

	_png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING , 0 , 0 , 0);
#ifdef NEW_CODE
	if( !_png_ptr ) ERROR_OUT( "failed to create png pointer" );
#else // !NEW_CODE
	if( !_png_ptr ) Miscellany::ErrorOut( "failed to create png pointer" );
#endif // NEW_CODE
	_info_ptr = png_create_info_struct( _png_ptr );
#ifdef NEW_CODE
	if( !_info_ptr ) ERROR_OUT( "failed to create info pointer" );
#else // 1NEW_CODE
	if( !_info_ptr ) Miscellany::ErrorOut( "failed to create info pointer" );
#endif // NEW_CODE

	_end_info = png_create_info_struct( _png_ptr );
#ifdef NEW_CODE
	if( !_end_info ) ERROR_OUT( "failed to create end pointer" );
#else // !NEW_CODE
	if( !_end_info ) Miscellany::ErrorOut( "failed to create end pointer" );
#endif // NEW_CODE


	_fp = fopen( fileName , "rb" );
#ifdef NEW_CODE
	if( !_fp ) ERROR_OUT( "Failed to open file for reading: " , std::string( fileName ) );
#else // !NEW_CODE
	if( !_fp ) Miscellany::ErrorOut( "Failed to open file for reading: %s" , fileName );
#endif // NEW_CODE
	png_init_io( _png_ptr , _fp );

	png_read_info( _png_ptr, _info_ptr );

	width = png_get_image_width( _png_ptr , _info_ptr );
	height = png_get_image_height( _png_ptr, _info_ptr );
	channels = png_get_channels( _png_ptr , _info_ptr );
	int bit_depth=png_get_bit_depth( _png_ptr , _info_ptr );
	int color_type = png_get_color_type( _png_ptr , _info_ptr );
#ifdef NEW_CODE
	if( bit_depth!=BitDepth ) ERROR_OUT( "expected " , BitDepth ,  " bits per channel: " , bit_depth );
#else // !NEW_CODE
	if( bit_depth!=BitDepth ) Miscellany::ErrorOut( "expected %d bits per channel: %d" , BitDepth , bit_depth );
#endif // NEW_CODE
	if( color_type==PNG_COLOR_TYPE_PALETTE ) png_set_expand( _png_ptr ) , printf( "Expanding PNG color pallette\n" );

	{
		long int a = 1;
		int swap = (*((unsigned char *) &a) == 1);
		if( swap ) png_set_swap( _png_ptr );
	}
}

template< unsigned int BitDepth >
inline unsigned int PNGReader< BitDepth >::nextRow( ChannelType * row )
{
	png_read_row( _png_ptr , (png_bytep)row , NULL );
	return _currentRow++;
}

template< unsigned int BitDepth >
PNGReader< BitDepth >::~PNGReader( void ){ }

template< unsigned int BitDepth >
inline bool PNGReader< BitDepth >::GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	unsigned int bitDepth;
	return GetInfo( fileName , width , height , channels , bitDepth );
}

template< unsigned int BitDepth >
inline bool PNGReader< BitDepth >::GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels , unsigned int &bitDepth )
{
	png_structp png_ptr;
	png_infop info_ptr;
	png_infop end_info ;
	FILE* fp;

	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING , 0 , 0 , 0);
#ifdef NEW_CODE
	if( !png_ptr ) ERROR_OUT( "failed to create png pointer" );
#else // !NEW_CODE
	if( !png_ptr ) Miscellany::ErrorOut( "failed to create png pointer" );
#endif // NEW_CODE
	info_ptr = png_create_info_struct( png_ptr );
#ifdef NEW_CODE
	if( !info_ptr ) ERROR_OUT( "failed to create info pointer" );
#else // !NEW_CODE
	if( !info_ptr ) Miscellany::ErrorOut( "failed to create info pointer" );
#endif // NEW_CODE
	end_info = png_create_info_struct( png_ptr );
#ifdef NEW_CODE
	if( !end_info ) ERROR_OUT( "failed to create end pointer" );
#else // !NEW_CODE
	if( !end_info ) Miscellany::ErrorOut( "failed to create end pointer" );
#endif // NEW_CODE

	fp = fopen( fileName , "rb" );
#ifdef NEW_CODE
	if( !fp ) ERROR_OUT( "Failed to open file for reading: " , std::string( fileName ) );
#else // !NEW_CODE
	if( !fp ) Miscellany::ErrorOut( "Failed to open file for reading: %s" , fileName );
#endif // NEW_CODE
	png_init_io( png_ptr , fp );

	png_read_info( png_ptr, info_ptr );

	width = png_get_image_width( png_ptr , info_ptr );
	height = png_get_image_height( png_ptr, info_ptr );
	channels = png_get_channels( png_ptr , info_ptr );
	bitDepth = png_get_bit_depth( png_ptr , info_ptr );

	png_destroy_read_struct( &png_ptr , &info_ptr , &end_info );
	fclose( fp );
	return true;
}

template< unsigned int BitDepth >
PNGWriter< BitDepth >::PNGWriter( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int /*quality*/ )
{
	_currentRow = 0;

	_png_ptr = png_create_write_struct( PNG_LIBPNG_VER_STRING , 0 , 0 , 0 );
#ifdef NEW_CODE
	if( !_png_ptr )	ERROR_OUT( "Failed to create png write struct" );
#else // !NEW_CODE
	if( !_png_ptr )	Miscellany::ErrorOut( "Failed to create png write struct" );
#endif // NEW_CODE
	_info_ptr = png_create_info_struct( _png_ptr );
#ifdef NEW_CODE
	if( !_info_ptr ) ERROR_OUT( "Failed to create png info struct" );
#else // !NEW_CODE
	if( !_info_ptr ) Miscellany::ErrorOut( "Failed to create png info struct" );
#endif // NEW_CODE

	_fp = fopen( fileName , "wb" );
#ifdef NEW_CODE
	if( !_fp ) ERROR_OUT( "Failed to open file for writing: " , std::string( fileName ) );
#else // !NEW_CODE
	if( !_fp ) Miscellany::ErrorOut( "Failed to open file for writing: %s" , fileName );
#endif // NEW_CODE
	png_init_io( _png_ptr , _fp );

	const int compression_effort = 0;  // was Z_BEST_SPEED
	png_set_compression_level( _png_ptr , compression_effort );

	int pngColorType;
	switch( channels )
	{
		case 1: pngColorType = PNG_COLOR_TYPE_GRAY ; break;
		case 3: pngColorType = PNG_COLOR_TYPE_RGB  ; break;
		case 4: pngColorType = PNG_COLOR_TYPE_RGBA ; break;
#ifdef NEW_CODE
		default: ERROR_OUT( "Only 1, 3, or 4 channel PNGs are supported" );
#else // !NEW_CODE
		default: Miscellany::ErrorOut( "Only 1, 3, or 4 channel PNGs are supported" );
#endif // NEW_CODE
	};
	png_set_IHDR( _png_ptr , _info_ptr, width , height, BitDepth , pngColorType , PNG_INTERLACE_NONE , PNG_COMPRESSION_TYPE_DEFAULT , PNG_FILTER_TYPE_DEFAULT );
	png_write_info( _png_ptr , _info_ptr );

	{
		long int a = 1;
		int swap = (*((unsigned char *) &a) == 1);
		if( swap ) png_set_swap( _png_ptr );
	}
}

template< unsigned int BitDepth >
PNGWriter< BitDepth >::~PNGWriter( void )
{
	png_write_end( _png_ptr , NULL );
	png_destroy_write_struct( &_png_ptr , &_info_ptr );
	fclose( _fp );
}

template< unsigned int BitDepth >
unsigned int PNGWriter< BitDepth >::nextRow( const ChannelType * row )
{
	png_write_row( _png_ptr , (png_bytep)row );
	return _currentRow++;
}

template< unsigned int BitDepth >
unsigned int PNGWriter< BitDepth >::nextRows( const ChannelType * rows , unsigned int rowNum )
{
	for( unsigned int r=0 ; r<rowNum ; r++ ) png_write_row( _png_ptr , (png_bytep)( rows + r * 3 * sizeof( ChannelType ) * png_get_image_width(_png_ptr, _info_ptr) ));
	return _currentRow += rowNum;
}

