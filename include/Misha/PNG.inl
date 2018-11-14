#include <stdio.h>
#include <vector>
#ifdef _WIN32
#include "PNG/png.h"
#else // !_WIN32
#include <png.h>
#endif // _WIN32
#include "Miscellany.h"

inline PNGReader::PNGReader( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	_currentRow = 0;

	_png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING , 0 , 0 , 0);
	if( !_png_ptr ) Miscellany::ErrorOut( "failed to create png pointer" );
	_info_ptr = png_create_info_struct( _png_ptr );
	if( !_info_ptr ) Miscellany::ErrorOut( "failed to create info pointer" );

	_end_info = png_create_info_struct( _png_ptr );
	if( !_end_info ) Miscellany::ErrorOut( "failed to create end pointer" );


	_fp = fopen( fileName , "rb" );
	if( !_fp ) Miscellany::ErrorOut( "Failed to open file for reading: %s" , fileName );
	png_init_io( _png_ptr , _fp );

	png_read_info( _png_ptr, _info_ptr );

	width = png_get_image_width( _png_ptr , _info_ptr );
	height = png_get_image_height( _png_ptr, _info_ptr );
	channels = png_get_channels( _png_ptr , _info_ptr );
	int bit_depth=png_get_bit_depth( _png_ptr , _info_ptr );
	int color_type = png_get_color_type( _png_ptr , _info_ptr );
	if( bit_depth!=8 ) Miscellany::ErrorOut( "expected 8 bits per channel" );
	if( color_type==PNG_COLOR_TYPE_PALETTE ) png_set_expand( _png_ptr ) , printf( "Expanding PNG color pallette\n" );

	{
		long int a = 1;
		int swap = (*((unsigned char *) &a) == 1);
		if( swap ) png_set_swap( _png_ptr );
	}
}
inline unsigned int PNGReader::nextRow( unsigned char* row )
{
	png_read_row( _png_ptr , (png_bytep)row , NULL );
	return _currentRow++;
}

PNGReader::~PNGReader( void ){ }

inline bool PNGReader::GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	png_structp png_ptr;
	png_infop info_ptr;
	png_infop end_info ;
	FILE* fp;

	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING , 0 , 0 , 0);
	if( !png_ptr ) Miscellany::ErrorOut( "failed to create png pointer" );
	info_ptr = png_create_info_struct( png_ptr );
	if( !info_ptr ) Miscellany::ErrorOut( "failed to create info pointer" );
	end_info = png_create_info_struct( png_ptr );
	if( !end_info ) Miscellany::ErrorOut( "failed to create end pointer" );

	fp = fopen( fileName , "rb" );
	if( !fp ) Miscellany::ErrorOut( "Failed to open file for reading: %s" , fileName );
	png_init_io( png_ptr , fp );

	png_read_info( png_ptr, info_ptr );

	width = png_get_image_width( png_ptr , info_ptr );
	height = png_get_image_height( png_ptr, info_ptr );
	channels = png_get_channels( png_ptr , info_ptr );

	png_destroy_read_struct( &png_ptr , &info_ptr , &end_info );
	fclose( fp );
	return true;
}

PNGWriter::PNGWriter( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int quality )
{
	_currentRow = 0;

	_png_ptr = png_create_write_struct( PNG_LIBPNG_VER_STRING , 0 , 0 , 0 );
	if( !_png_ptr )	Miscellany::ErrorOut( "Failed to create png write struct" );
	_info_ptr = png_create_info_struct( _png_ptr );
	if( !_info_ptr ) Miscellany::ErrorOut( "Failed to create png info struct" );

	_fp = fopen( fileName , "wb" );
	if( !_fp ) Miscellany::ErrorOut( "Failed to open file for writing: %s" , fileName );
	png_init_io( _png_ptr , _fp );

	const int compression_effort = 0;  // was Z_BEST_SPEED
	png_set_compression_level( _png_ptr , compression_effort );

	int pngColorType;
	switch( channels )
	{
	case 1: pngColorType = PNG_COLOR_TYPE_GRAY ; break;
	case 3: pngColorType = PNG_COLOR_TYPE_RGB  ; break;
	case 4: pngColorType = PNG_COLOR_TYPE_RGBA ; break;
	default: Miscellany::ErrorOut( "Only 1, 3, or 4 channel PNGs are supported" );
	};
	png_set_IHDR( _png_ptr , _info_ptr, width , height, 8 , pngColorType , PNG_INTERLACE_NONE , PNG_COMPRESSION_TYPE_DEFAULT , PNG_FILTER_TYPE_DEFAULT );
	png_write_info( _png_ptr , _info_ptr );

	{
		long int a = 1;
		int swap = (*((unsigned char *) &a) == 1);
		if( swap ) png_set_swap( _png_ptr );
	}
}
PNGWriter::~PNGWriter( void )
{
	png_write_end( _png_ptr , NULL );
	png_destroy_write_struct( &_png_ptr , &_info_ptr );
	fclose( _fp );
}
unsigned int PNGWriter::nextRow( const unsigned char* row )
{
	png_write_row( _png_ptr , (png_bytep)row );
	return _currentRow++;
}
unsigned int PNGWriter::nextRows( const unsigned char* rows , unsigned int rowNum )
{
	for( unsigned int r=0 ; r<rowNum ; r++ ) png_write_row( _png_ptr , (png_bytep)( rows + r * 3 * sizeof( unsigned char ) * png_get_image_width(_png_ptr, _info_ptr) ));
	return _currentRow += rowNum;
}

