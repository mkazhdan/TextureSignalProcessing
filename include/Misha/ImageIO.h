#ifndef IMAGE_IO_INCLUDED
#define IMAGE_IO_INCLUDED

#if defined( NEW_CMD_LINE_PARSER ) || 1
#include <string.h>
#endif // NEW_CMD_LINE_PARSER

struct ImageReader
{
	virtual unsigned int nextRow( unsigned char* row ) = 0;
	static unsigned char* Read( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
	{
		ImageReader* reader = Get( fileName , width , height , channels );
		unsigned char* pixels = new unsigned char[ width*height*channels ];
		for( unsigned int j=0 ; j<height ; j++ ) reader->nextRow( pixels + j*width*channels );
		delete reader;
		return pixels;
	}
	static unsigned char* ReadColor( const char* fileName , unsigned int& width , unsigned int& height )
	{
		unsigned int channels;
		ImageReader* reader = Get( fileName , width , height , channels );
		if( channels!=1 && channels!=3 ) fprintf( stderr , "[ERROR] ImageReader::ReadColor requres one- or three-channel input\n" ) , exit( 0 );
		unsigned char* pixels = new unsigned char[ width*height*3 ];
		unsigned char* pixelRow = new unsigned char[ width*channels];
		for( unsigned int j=0 ; j<height ; j++ )
		{
			reader->nextRow( pixelRow );
			if     ( channels==3 ) memcpy( pixels+j*width*3 , pixelRow , sizeof(unsigned char)*width*3 );
			else if( channels==1 ) for( unsigned int i=0 ; i<width ; i++ ) for( unsigned int c=0 ; c<3 ; c++ ) pixels[j*width*3+i*3+c] = pixelRow[i];
		}
		delete[] pixelRow;
		delete reader;
		return pixels;
	}

	static ImageReader* Get( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
	static void GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels );
	virtual ~ImageReader( void ){ }
};
struct ImageWriter
{
	virtual unsigned int nextRow( const unsigned char* row ) = 0;
	virtual unsigned int nextRows( const unsigned char* rows , unsigned int rowNum ) = 0;
	static bool Write( const char* fileName , const unsigned char* pixels , unsigned int width , unsigned int height , int channels , int quality=100 )
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
inline ImageReader* ImageReader::Get( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
	ImageReader* reader = NULL;
	char* ext = GetImageExtension( fileName );
#ifdef WIN32
	if( !_stricmp( ext , "jpeg" ) || !_stricmp( ext , "jpg" ) ) reader = new JPEGReader( fileName , width , height , channels );
	else if( !_stricmp( ext , "png" ) ) reader = new PNGReader( fileName , width , height , channels );
#else // !WIN32
	if( !strcasecmp( ext , "jpeg" ) || !strcasecmp( ext , "jpg" ) ) reader = new JPEGReader( fileName , width , height , channels );
	else if( !strcasecmp( ext , "png" ) ) reader = new PNGReader( fileName , width , height , channels );
#endif // WIN32

	delete[] ext;
	return reader;
}
inline void ImageReader::GetInfo( const char* fileName , unsigned int& width , unsigned int& height , unsigned int& channels )
{
#if 1
	char* ext = GetImageExtension( fileName );
#ifdef WIN32
	if( !_stricmp( ext , "jpeg" ) || !_stricmp( ext , "jpg" ) ) JPEGReader::GetInfo( fileName , width , height , channels );
	else if( !_stricmp( ext , "png" ) )                          PNGReader::GetInfo( fileName , width , height , channels );
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
inline ImageWriter* ImageWriter::Get( const char* fileName , unsigned int width , unsigned int height , unsigned int channels , unsigned int quality )
{
	ImageWriter* writer = NULL;
	char* ext = GetImageExtension( fileName );
#ifdef WIN32
	if( !_stricmp( ext , "jpeg" ) || !_stricmp( ext , "jpg" ) ) writer = new JPEGWriter( fileName , width , height , channels , quality );
	else if( !_stricmp( ext , "png" ) ) writer = new PNGWriter( fileName , width , height , channels , quality );
#else // !WIN32
	if( !strcasecmp( ext , "jpeg" ) || !strcasecmp( ext , "jpg" ) ) writer = new JPEGWriter( fileName , width , height , channels , quality );
	else if( !strcasecmp( ext , "png" ) ) writer = new PNGWriter( fileName , width , height , channels , quality );
#endif // WIN32

	delete[] ext;
	return writer;
}

#endif // IMAGE_INCLUDED
