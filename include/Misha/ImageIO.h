#ifndef IMAGE_IO_INCLUDED
#define IMAGE_IO_INCLUDED

#if defined( NEW_CMD_LINE_PARSER ) || 1
#include <string.h>
#endif // NEW_CMD_LINE_PARSER
#include "Miscellany.h"

#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
template< unsigned int BitDepth > struct ImageChannel;

template<> struct ImageChannel< 8>{ using Type = uint8_t; };
template<> struct ImageChannel<16>{ using Type = uint16_t; };
template<> struct ImageChannel<32>{ using Type = uint32_t; };
template<> struct ImageChannel<64>{ using Type = uint64_t; };
#endif // VARIABLE_SIZED_IMAGE_CHANNEL

#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
template< unsigned int BitDepth >
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
struct ImageReader
{
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
	using ChannelType = typename ImageChannel< BitDepth >::Type;
	virtual unsigned int nextRow( ChannelType * row ) = 0;
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
	virtual unsigned int nextRow( unsigned char* row ) = 0;
#endif // VARIABLE_SIZED_IMAGE_CHANNEL

#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
	static ChannelType * Read( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
	static unsigned char* Read( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
	{
		ImageReader* reader = Get( fileName , width , height , channels );
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
		ChannelType * pixels = new ChannelType[ width*height*channels ];
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
		unsigned char* pixels = new unsigned char[ width*height*channels ];
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
		for( unsigned int j=0 ; j<height ; j++ ) reader->nextRow( pixels + j*width*channels );
		delete reader;
		return pixels;
	}

#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
	static ChannelType * ReadColor( const char* fileName , unsigned int& width , unsigned int& height )
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
	static unsigned char* ReadColor( const char* fileName , unsigned int& width , unsigned int& height )
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
	{
		unsigned int channels;
		ImageReader* reader = Get( fileName , width , height , channels );
#ifdef RGBA_IMAGE
		if( channels!=1 && channels!=3 && channels!=4 ) Miscellany::ErrorOut( "Requires one-, three-, or four-channel input: %d" , channels );
#else // !RGBA_IMAGE
		if( channels!=1 && channels!=3 ) Miscellany::ErrorOut( "Requres one- or three-channel input: %d" , channels );
#endif // RGBA_IMAGE
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
		ChannelType * pixels = new ChannelType[ width*height*3 ];
		ChannelType * pixelRow = new ChannelType[ width*channels ];
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
		unsigned char* pixels = new unsigned char[ width*height*3 ];
		unsigned char* pixelRow = new unsigned char[ width*channels ];
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
		for( unsigned int j=0 ; j<height ; j++ )
		{
			reader->nextRow( pixelRow );
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
			if     ( channels==3 ) memcpy( pixels+j*width*3 , pixelRow , sizeof(ChannelType)*width*3 );
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
			if     ( channels==3 ) memcpy( pixels+j*width*3 , pixelRow , sizeof(unsigned char)*width*3 );
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
#ifdef RGBA_IMAGE
			else if( channels==4 ) for( unsigned int i=0 ; i<width ; i++ ) for( unsigned int c=0 ; c<3 ; c++ ) pixels[j*width*3+i*3+c] = pixelRow[i*channels+c];
#endif // RGBA_IMAGE
			else if( channels==1 ) for( unsigned int i=0 ; i<width ; i++ ) for( unsigned int c=0 ; c<3 ; c++ ) pixels[j*width*3+i*3+c] = pixelRow[i];
		}
		delete[] pixelRow;
		delete reader;
		return pixels;
	}

	static ImageReader* Get( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
	static void GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels , unsigned int &bitDepth );
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
	static void GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
	virtual ~ImageReader( void ){ }
};

#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
template< unsigned int BitDepth >
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
struct ImageWriter
{
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
	using ChannelType = typename ImageChannel< BitDepth >::Type;
#endif // VARIABLE_SIZED_IMAGE_CHANNEL

#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
	virtual unsigned int nextRow( const ChannelType * row ) = 0;
	virtual unsigned int nextRows( const ChannelType * rows , unsigned int rowNum ) = 0;
	static bool Write( const char* fileName , const ChannelType * pixels , unsigned int width , unsigned int height , int channels , int quality=100 )
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
	virtual unsigned int nextRow( const unsigned char* row ) = 0;
	virtual unsigned int nextRows( const unsigned char* rows , unsigned int rowNum ) = 0;
	static bool Write( const char* fileName , const unsigned char* pixels , unsigned int width , unsigned int height , int channels , int quality=100 )
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
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

#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
template< unsigned int BitDepth >
inline ImageReader< BitDepth >* ImageReader< BitDepth >::Get( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
inline ImageReader* ImageReader::Get( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
{
	ImageReader< BitDepth > * reader = nullptr;
	char* ext = GetImageExtension( fileName );
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
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
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
#ifdef WIN32
	if( !_stricmp( ext , "jpeg" ) || !_stricmp( ext , "jpg" ) ) reader = new JPEGReader( fileName , width , height , channels );
	else if( !_stricmp( ext , "png" ) ) reader = new PNGReader( fileName , width , height , channels );
#else // !WIN32
	if( !strcasecmp( ext , "jpeg" ) || !strcasecmp( ext , "jpg" ) ) reader = new JPEGReader( fileName , width , height , channels );
	else if( !strcasecmp( ext , "png" ) ) reader = new PNGReader( fileName , width , height , channels );
#endif // WIN32
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
	delete[] ext;
	return reader;
}

#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
template< unsigned int BitDepth >
inline void ImageReader< BitDepth >::GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	unsigned int bitDepth;
	return GetInfo( fileName , width , height , channels );
}
template< unsigned int BitDepth >
inline void ImageReader< BitDepth >::GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels , unsigned int &bitDepth )
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
inline void ImageReader::GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
{
#if 1
	char* ext = GetImageExtension( fileName );
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
#ifdef WIN32
	if( !_stricmp( ext , "jpeg" ) || !_stricmp( ext , "jpg" ) ) JPEGReader::GetInfo( fileName , width , height , channels , bitDepth );
	else if( !_stricmp( ext , "png" ) )                          PNGReader< BitDepth >::GetInfo( fileName , width , height , channels , bitDepth );
#else // !WIN32
	if( !strcasecmp( ext , "jpeg" ) || !strcasecmp( ext , "jpg" ) ) void();  //??
	// delete new JPEGReader( fileName , width , height , channels , quality );
	else if( !strcasecmp( ext , "png" ) ) void();  //??
	// delete new PNGReader( fileName , width , height , channels , quality );
#endif // WIN32
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
#ifdef WIN32
	if( !_stricmp( ext , "jpeg" ) || !_stricmp( ext , "jpg" ) ) JPEGReader::GetInfo( fileName , width , height , channels );
	else if( !_stricmp( ext , "png" ) )                          PNGReader::GetInfo( fileName , width , height , channels );
#else // !WIN32
	if( !strcasecmp( ext , "jpeg" ) || !strcasecmp( ext , "jpg" ) ) void();  //??
																			 // delete new JPEGReader( fileName , width , height , channels , quality );
	else if( !strcasecmp( ext , "png" ) ) void();  //??
												   // delete new PNGReader( fileName , width , height , channels , quality );
#endif // WIN32
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
#else
	ImageReader* reader = Get( fileName , width , height , channels );
	delete reader;
#endif
}

#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
template< unsigned int BitDepth >
inline ImageWriter< BitDepth > * ImageWriter< BitDepth >::Get( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int quality )
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
inline ImageWriter* ImageWriter::Get( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int quality )
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
{
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
	ImageWriter< BitDepth > * writer = nullptr;
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
	ImageWriter* writer = NULL;
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
	char* ext = GetImageExtension( fileName );
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
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
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
#ifdef WIN32
	if( !_stricmp( ext , "jpeg" ) || !_stricmp( ext , "jpg" ) ) writer = new JPEGWriter( fileName , width , height , channels , quality );
	else if( !_stricmp( ext , "png" ) ) writer = new PNGWriter( fileName , width , height , channels , quality );
#else // !WIN32
	if( !strcasecmp( ext , "jpeg" ) || !strcasecmp( ext , "jpg" ) ) writer = new JPEGWriter( fileName , width , height , channels , quality );
	else if( !strcasecmp( ext , "png" ) ) writer = new PNGWriter( fileName , width , height , channels , quality );
#endif // WIN32
#endif // VARIABLE_SIZED_IMAGE_CHANNEL

	delete[] ext;
	return writer;
}

#endif // IMAGE_INCLUDED
