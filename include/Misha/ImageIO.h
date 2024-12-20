#ifndef IMAGE_IO_INCLUDED
#define IMAGE_IO_INCLUDED

#if defined( NEW_CMD_LINE_PARSER ) || 1
#include <string.h>
#endif // NEW_CMD_LINE_PARSER
#include "Miscellany.h"

template< unsigned int BitDepth > struct ImageChannel;

template<> struct ImageChannel< 8>{ using Type = uint8_t; };
template<> struct ImageChannel<16>{ using Type = uint16_t; };
template<> struct ImageChannel<32>{ using Type = uint32_t; };
template<> struct ImageChannel<64>{ using Type = uint64_t; };

template< unsigned int BitDepth >
struct ImageReader
{
	using ChannelType = typename ImageChannel< BitDepth >::Type;
	virtual unsigned int nextRow( ChannelType * row ) = 0;

	static ChannelType * Read( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
	{
		ImageReader* reader = Get( fileName , width , height , channels );
		ChannelType * pixels = new ChannelType[ width*height*channels ];
		for( unsigned int j=0 ; j<height ; j++ ) reader->nextRow( pixels + j*width*channels );
		delete reader;
		return pixels;
	}

	static ChannelType * ReadColor( const char* fileName , unsigned int& width , unsigned int& height )
	{
		unsigned int channels;
		ImageReader* reader = Get( fileName , width , height , channels );
		if( channels!=1 && channels!=3 && channels!=4 ) Miscellany::ErrorOut( "Requires one-, three-, or four-channel input: %d" , channels );
		ChannelType * pixels = new ChannelType[ width*height*3 ];
		ChannelType * pixelRow = new ChannelType[ width*channels ];
		for( unsigned int j=0 ; j<height ; j++ )
		{
			reader->nextRow( pixelRow );
			if     ( channels==3 ) memcpy( pixels+j*width*3 , pixelRow , sizeof(ChannelType)*width*3 );
			else if( channels==4 ) for( unsigned int i=0 ; i<width ; i++ ) for( unsigned int c=0 ; c<3 ; c++ ) pixels[j*width*3+i*3+c] = pixelRow[i*channels+c];
			else if( channels==1 ) for( unsigned int i=0 ; i<width ; i++ ) for( unsigned int c=0 ; c<3 ; c++ ) pixels[j*width*3+i*3+c] = pixelRow[i];
		}
		delete[] pixelRow;
		delete reader;
		return pixels;
	}

	static ImageReader* Get( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
	static void GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels , unsigned int &bitDepth );
	static void GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
	virtual ~ImageReader( void ){ }
};

template< unsigned int BitDepth >
struct ImageWriter
{
	using ChannelType = typename ImageChannel< BitDepth >::Type;

	virtual unsigned int nextRow( const ChannelType * row ) = 0;
	virtual unsigned int nextRows( const ChannelType * rows , unsigned int rowNum ) = 0;
	static bool Write( const char* fileName , const ChannelType * pixels , unsigned int width , unsigned int height , int channels , int quality=100 )
	{
		ImageWriter* writer = Get( fileName , width , height , channels , quality );
		if( !writer ) return false;
		for( unsigned int j=0 ; j<height ; j++ ) writer->nextRow( pixels + j*width*channels );
		delete writer;
		return true;
	}
	static ImageWriter* Get( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int quality=100 );
	virtual ~ImageWriter( void ){ }
};

// May need to include PNG.h before JPEG.h to avoid the issue of multiple inclusions of setjmp.h under Linux
#include "PNG.h"
#include "JPEG.h"


inline char* GetImageExtension( const char* imageName )
{
	char *imageNameCopy , *temp , *ext = NULL;

	imageNameCopy = new char[ strlen(imageName)+1 ];
	strcpy( imageNameCopy , imageName );
	temp = strtok( imageNameCopy , "." );
	while( temp )
	{
		if( ext ) delete[] ext;
		ext = new char[ strlen(temp)+1 ];
		strcpy( ext , temp );
		temp = strtok( NULL , "." );
	}
	delete[] imageNameCopy;
	return ext;
}

template< unsigned int BitDepth >
inline ImageReader< BitDepth >* ImageReader< BitDepth >::Get( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	ImageReader< BitDepth > * reader = nullptr;
	char* ext = GetImageExtension( fileName );
#ifdef WIN32
	if( !_stricmp( ext , "jpeg" ) || !_stricmp( ext , "jpg" ) )
	{
		if constexpr( BitDepth==8 ) reader = new JPEGReader( fileName , width , height , channels );
	}
	else if( !_stricmp( ext , "png" ) ) reader = new PNGReader< BitDepth >( fileName , width , height , channels );
#else // !WIN32
	if( !strcasecmp( ext , "jpeg" ) || !strcasecmp( ext , "jpg" )  )
	{
		if constexpr( BitDepth==8 ) reader = new JPEGReader( fileName , width , height , channels );
	}
	else if( !strcasecmp( ext , "png" ) ) reader = new PNGReader< BitDepth >( fileName , width , height , channels );
#endif // WIN32
	delete[] ext;
	return reader;
}

template< unsigned int BitDepth >
inline void ImageReader< BitDepth >::GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	unsigned int bitDepth;
	return GetInfo( fileName , width , height , channels );
}
template< unsigned int BitDepth >
inline void ImageReader< BitDepth >::GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels , unsigned int &bitDepth )
{
#if 1
	char* ext = GetImageExtension( fileName );
#ifdef WIN32
	if( !_stricmp( ext , "jpeg" ) || !_stricmp( ext , "jpg" ) ) JPEGReader::GetInfo( fileName , width , height , channels , bitDepth );
	else if( !_stricmp( ext , "png" ) )                          PNGReader< BitDepth >::GetInfo( fileName , width , height , channels , bitDepth );
#else // !WIN32
	if( !strcasecmp( ext , "jpeg" ) || !strcasecmp( ext , "jpg" ) ) void();  //??
	// delete new JPEGReader( fileName , width , height , channels , quality );
	else if( !strcasecmp( ext , "png" ) ) void();  //??
	// delete new PNGReader( fileName , width , height , channels , quality );
#endif // WIN32
#else
	ImageReader* reader = Get( fileName , width , height , channels );
	delete reader;
#endif
}

template< unsigned int BitDepth >
inline ImageWriter< BitDepth > * ImageWriter< BitDepth >::Get( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int quality )
{
	ImageWriter< BitDepth > * writer = nullptr;
	char* ext = GetImageExtension( fileName );
#ifdef WIN32
	if( !_stricmp( ext , "jpeg" ) || !_stricmp( ext , "jpg" ) )
	{
		if constexpr( BitDepth==8 ) writer = new JPEGWriter( fileName , width , height , channels , quality );
	}
	else if( !_stricmp( ext , "png" ) ) writer = new PNGWriter< BitDepth >( fileName , width , height , channels , quality );
#else // !WIN32
	if( !strcasecmp( ext , "jpeg" ) || !strcasecmp( ext , "jpg" ) )
	{
		if constexpr( BitDepth==8 ) writer = new JPEGWriter( fileName , width , height , channels , quality );
	}
	else if( !strcasecmp( ext , "png" ) ) writer = new PNGWriter< BitDepth >( fileName , width , height , channels , quality );
#endif // WIN32

	delete[] ext;
	return writer;
}

#endif // IMAGE_INCLUDED
