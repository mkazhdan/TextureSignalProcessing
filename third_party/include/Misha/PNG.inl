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

template< unsigned int BitDepth >
PNGReader< BitDepth >::PNGReader( std::string fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	_currentRow = 0;

	_png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING , 0 , 0 , 0);
	if( !_png_ptr ) MK_ERROR_OUT( "failed to create png pointer" );
	_info_ptr = png_create_info_struct( _png_ptr );
	if( !_info_ptr ) MK_ERROR_OUT( "failed to create info pointer" );

	_end_info = png_create_info_struct( _png_ptr );
	if( !_end_info ) MK_ERROR_OUT( "failed to create end pointer" );


#if _WIN32 || _WIN64
	if( fopen_s( &_fp , fileName.c_str() , "rb" ) ) fprintf( stderr , "[ERROR] PNGReader: Failed to open file for reading: %s\n" , fileName.c_str() ) , exit( 0 );
#else // !_WIN32 && !_WIN64
	_fp = fopen( fileName.c_str() , "rb" );
	if( !_fp ) fprintf( stderr , "[ERROR] PNGReader: Failed to open file for reading: %s\n" , fileName.c_str() ) , exit( 0 );
#endif // _WIN32 || _WIN64
	png_init_io( _png_ptr , _fp );

	png_read_info( _png_ptr, _info_ptr );

	width = png_get_image_width( _png_ptr , _info_ptr );
	height = png_get_image_height( _png_ptr, _info_ptr );
	channels = png_get_channels( _png_ptr , _info_ptr );
	int bit_depth=png_get_bit_depth( _png_ptr , _info_ptr );
	int color_type = png_get_color_type( _png_ptr , _info_ptr );
	if( bit_depth!=BitDepth ) MK_ERROR_OUT( "expected " , BitDepth ,  " bits per channel: " , bit_depth );
	if( color_type==PNG_COLOR_TYPE_PALETTE ) png_set_expand( _png_ptr ) , printf( "Expanding PNG color pallette\n" );

	{
		long int a = 1;
		int swap = (*((unsigned char *) &a) == 1);
		if( swap ) png_set_swap( _png_ptr );
	}
}

template< unsigned int BitDepth >
inline unsigned int PNGReader< BitDepth >::nextRow( ChannelType* row )
{
	png_read_row( _png_ptr , (png_bytep)row , NULL );
	return _currentRow++;
}

template< unsigned int BitDepth >
PNGReader< BitDepth >::~PNGReader( void ){ }

template< unsigned int BitDepth >
inline bool PNGReader< BitDepth >::GetInfo( std::string fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	unsigned int bitDepth;
	return GetInfo( fileName , width , height , channels );
}

template< unsigned int BitDepth >
inline bool PNGReader< BitDepth >::GetInfo( std::string fileName , unsigned int& width , unsigned int& height , unsigned int& channels , unsigned int &bitDepth )
{
	png_structp png_ptr;
	png_infop info_ptr;
	png_infop end_info ;
	FILE* fp;

	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING , 0 , 0 , 0);
	if( !png_ptr ) MK_ERROR_OUT( "failed to create png pointer" );
	info_ptr = png_create_info_struct( png_ptr );
	if( !info_ptr ) MK_ERROR_OUT( "failed to create info pointer" );
	end_info = png_create_info_struct( png_ptr );
	if( !end_info ) MK_ERROR_OUT( "failed to create end pointer" );

#if _WIN32 || _WIN64
	if( fopen_s( &fp , fileName.c_str() , "rb" ) ) MK_ERROR_OUT( "Failed to open file for reading: " , fileName );
#else // !_WIN32 && !_WIN64
	fp = fopen( fileName.c_str() , "rb" );
	if( !fp ) MK_ERROR_OUT( "Failed to open file for reading: " , fileName );
#endif // _WIN32 || _WIN64
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
PNGWriter< BitDepth >::PNGWriter( std::string fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int quality )
{
	_currentRow = 0;

	_png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING , 0 , 0 , 0);
	if( !_png_ptr ) MK_ERROR_OUT( "failed to create png pointer" );
	_info_ptr = png_create_info_struct( _png_ptr );
	if( !_info_ptr ) MK_ERROR_OUT( "Failed to create png info struct" );

#if _WIN32 || _WIN64
	if( fopen_s( &_fp , fileName.c_str() , "wb" ) ) MK_ERROR_OUT( "Failed to open file for writing: " , fileName );
#else // !_WIN32 && !_WIN64
	_fp = fopen( fileName.c_str() , "wb" );
	if( !_fp ) MK_ERROR_OUT( "Failed to open file for writing: " , fileName );
#endif // _WIN32 || _WIN64
	png_init_io( _png_ptr , _fp );

	const int compression_effort = 0;  // was Z_BEST_SPEED
	png_set_compression_level( _png_ptr , compression_effort );

	int pngColorType;
	switch( channels )
	{
	case 1: pngColorType = PNG_COLOR_TYPE_GRAY ; break;
	case 3: pngColorType = PNG_COLOR_TYPE_RGB  ; break;
	case 4: pngColorType = PNG_COLOR_TYPE_RGBA ; break;
	default: MK_ERROR_OUT( "Only 1, 3, or 4 channel PNGs are supported" );
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
	for( unsigned int r=0 ; r<rowNum ; r++ ) png_write_row( _png_ptr , (png_bytep)( rows + r * 3 * sizeof( ChannelType ) * png_get_image_width( _png_ptr , _info_ptr ) ) );
	return _currentRow += rowNum;
}

#if 0
void PNGWriteColor( std::string fileName , const unsigned char* pixels , int width , int height )
{
	FILE* fp = fopen( fileName.c_str() , "wb" );
	if( !fp ) fprintf( stderr , "[ERROR] Failed to open file for writing: %s\n" , fileName.c_str() ) , exit( 0 );
	png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,0,0,0);
	if(!png_ptr)	return;
	png_infop info_ptr = png_create_info_struct(png_ptr);
	if(!info_ptr)	return;
	png_init_io(png_ptr, fp);
	// turn off compression or set another filter
	// png_set_filter(png_ptr, 0, PNG_FILTER_NONE);
	png_set_IHDR(png_ptr, info_ptr, width , height ,
		8,PNG_COLOR_TYPE_RGB,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);
	if (0) {                    // high-level write
		std::vector<unsigned char> matrix( width * height * 3 );
		std::vector<png_bytep> row_pointers( height );
		for(int y=0;y<height;y++)
		{
			row_pointers[y]=&matrix[y*width*3];
			unsigned char* buf=&matrix[y*width*3];
			for(int x=0;x<width;x++)
				for(int z=0;z<3;z++)
					*buf++ = pixels[ (y*width+x)*3 + z ];
		}
		png_set_rows(png_ptr, info_ptr, &row_pointers[0]);
		int png_transforms=0;
		png_write_png(png_ptr, info_ptr, png_transforms, NULL);
	} else {                    // low-level write
		png_write_info(png_ptr, info_ptr);
		// png_set_filler(png_ptr, 0, PNG_FILLER_AFTER);
		//  but no way to provide GRAY data with RGBA fill, so pack each row
		std::vector<unsigned char> buffer(width*3);
		for(int y=0;y<height;y++)
		{
			unsigned char* buf=&buffer[0];
			for(int x=0;x<width;x++)
				for(int z=0;z<3;z++)
					*buf++ = pixels[ (y*width+x)*3+z ];
			png_bytep row_pointer=&buffer[0];
			png_write_row(png_ptr, row_pointer);
		}
	}
}
#endif