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
#ifndef IMAGE_INCLUDED
#define IMAGE_INCLUDED

#include "Miscellany.h"
#include "Geometry.h"
#include "CmdLineParser.h"
#include "ImageIO.h"
#include "SparseMatrix.h"


template< class Data >
struct Image
{
	Image(void);
	Image(int w, int h);
	Image(const Image& img);
	~Image(void);
	void resize(int w, int h);
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
	template< unsigned int BitDepth > void read( const char* fileName );
	template< unsigned int BitDepth > void write( const char* fileName ) const;
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
	void read(const char* fileName);
	void write(const char* fileName) const;
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
	int width(void) const { return _width; }
	int height(void) const { return _height; }
	int size(void) const { return _width * _height; }
	const Data& operator()(int x, int y) const { return _pixels[y*_width + x]; }
	Data& operator()(int x, int y) { return _pixels[y*_width + x]; }
	const Data& operator[](int i) const { return _pixels[i]; }
	Data& operator[](int i) { return _pixels[i]; }
	template< class Real > Data sample(Real x, Real y) const;
	Image& operator = (const Image& img);
protected:
	int _width, _height;
	Data* _pixels;
};


template< class Data > Image< Data >::Image(void) : _width(0), _height(0), _pixels(NULL){ ; }
template< class Data > Image< Data >::Image(int w, int h) : _width(0), _height(0), _pixels(NULL){ resize(w, h); }
template< class Data > Image< Data >::Image(const Image& img) : _width(0), _height(0), _pixels(NULL)
{
	resize(img.width(), img.height());
	memcpy(_pixels, img._pixels, sizeof(Data) * img.size());
}
template< class Data > Image< Data >::~Image(void){ if (_pixels) delete[] _pixels; _pixels = NULL, _width = _height = 0; }
template< class Data >
template< class Real >
Data Image< Data >::sample(Real x, Real y) const
{
#if 1
	int ix1 = (int)floor(x), iy1 = (int)floor(y);
	Real dx = x - ix1, dy = y - iy1;
	ix1 = std::max< int >(0, std::min< int >(ix1, _width - 1));
	iy1 = std::max< int >(0, std::min< int >(iy1, _height - 1));
	int ix2 = std::min< int >(ix1 + 1, _width - 1), iy2 = std::min< int >(iy1 + 1, _height - 1);
	return
		((*this)(ix1, iy1) * (Real)(1. - dy) + (*this)(ix1, iy2) * (Real)(dy)) * (Real)(1. - dx) +
		((*this)(ix2, iy1) * (Real)(1. - dy) + (*this)(ix2, iy2) * (Real)(dy)) * (Real)(dx);
#else
	int ix1 = (int)floor(x), ix2 = ix1 + 1, iy1 = (int)floor(y), iy2 = iy1 + 1;
	float dx = x - ix1, dy = y - iy1;
	ix1 = std::max< int >(0, std::min< int >(ix1, _width - 1));
	ix2 = std::max< int >(0, std::min< int >(ix2, _width - 1));
	iy1 = std::max< int >(0, std::min< int >(iy1, _height - 1));
	iy2 = std::max< int >(0, std::min< int >(iy2, _height - 1));
	return
		(*this)(ix1, iy1) * (Real)((1. - dx) * (1. - dy)) +
		(*this)(ix1, iy2) * (Real)((1. - dx) * (dy)) +
		(*this)(ix2, iy1) * (Real)((dx)* (1. - dy)) +
		(*this)(ix2, iy2) * (Real)((dx)* (dy));
#endif
}

template< class Data > void Image< Data >::resize(int w, int h)
{
	if (_width*_height != w*h)
	{
		if (_pixels) delete[] _pixels;
		_pixels = NULL;
		if (w*h) _pixels = new Data[w*h];
		if (!_pixels) Miscellany::ErrorOut( "Failed to allocate pixels: %d x %d" , w , h );
	}
	_width = w, _height = h;
}

template< class Data > Image< Data >& Image< Data >::operator = (const Image& img)
{
	resize(img.width(), img.height());
	memcpy(_pixels, img._pixels, sizeof(Data) * img.size());
	return *this;
}

#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
template< typename Data >
template< unsigned int BitDepth >
void Image< Data >::read( const char* fileName )
{
	using CType = typename ImageChannel< BitDepth >::Type;
	static const CType Scale = ~((CType )0);
#ifdef NEW_CODE
	unsigned int width , height;
	CType * pixels = ImageReader< BitDepth >::ReadColor( fileName , width , height );
	if( !pixels ) Miscellany::Throw( "Failed to read image: %s\n" , fileName );
	resize( width , height );

	if constexpr( std::is_same_v< Data , Point3D< double > > )
	{
		double scale = (double)Scale;
		for( int i=0 ; i<(int)width ; i++ ) for( int j=0 ; j<(int)height ; j++ ) for( int c=0 ; c<3 ; c++ ) (*this)(i,j)[c] = ( pixels[ (j*width+i)*3+c ] ) / scale;
	}
	else if constexpr( std::is_same_v< Data , Point3D< float > > )
	{
		float scale = (float)Scale;
		for( int i=0 ; i<(int)width ; i++ ) for( int j=0 ; j<(int)height ; j++ ) for( int c=0 ; c<3 ; c++ ) (*this)(i,j)[c] = ( pixels[ (j*width+i)*3+c ] ) / scale;
	}
	else if constexpr( std::is_same_v< Data , Point3D< CType > > )
	{
		memcpy( _pixels , pixels , sizeof( CType ) * _width * _height * 3 );
	}
	else Miscellany::ErrorOut( "Bad data type " );

	delete[] pixels;
#else // !NEW_CODE
	unsigned int width , height , channels;
	unsigned char* pixels = ImageReader::Read( fileName , width , height , channels );
	if( !pixels ) Miscellany::Throw( "Failed to read image: %s\n" , fileName );
	if( channels!=3 ) Miscellany::Throw( "Only three channel images are supported: %d\n" , channels );
	resize( width , height );
	for( int i=0 ; i<(int)width ; i++ ) for( int j=0 ; j<(int)height ; j++ ) for( int c=0 ; c<(int)channels ; c++ ) (*this)(i,j)[c] = ( (float)pixels[ (j*width+i)*3+c ] ) / 255.f;
	delete[] pixels;
#endif // NEW_CODE
}


template< typename Data >
template< unsigned int BitDepth >
void Image< Data >::write( const char* fileName ) const
{
	using CType = typename ImageChannel< BitDepth >::Type;
	static const CType Scale = ~((CType )0);
	if constexpr( std::is_same_v< Data , Point3D< double > > )
	{
		double scale = (double)Scale;
		CType * pixels = new CType[ _width*_height*3 ];
		for( int i=0 ; i<_width ; i++ ) for( int j=0 ; j<_height ; j++ ) for( int c=0 ; c<3 ; c++ ) pixels[ 3*(j*_width+i)+c ] = (CType)std::min< unsigned long long >( Scale , std::max< unsigned long long >( 0 , (unsigned long long)( (*this)(i,j)[c] * scale + 0.5 ) ) );
		ImageWriter< BitDepth >::Write( fileName , pixels , _width , _height , 3 );
		delete[] pixels;
	}
	else if constexpr( std::is_same_v< Data , Point3D< float > > )
	{
		float scale = (float)Scale;
		CType * pixels = new CType[ _width * _height * 3 ];
		for( int i=0 ; i<_width ; i++ ) for( int j=0 ; j<_height ; j++ ) for( int c=0 ; c<3 ; c++ ) pixels[ 3*(j*_width+i)+c ] = (CType)std::min< unsigned long long >( Scale , std::max< unsigned long long >( 0 , (unsigned long long)( (*this)(i,j)[c] * scale + 0.5f ) ) );
		ImageWriter< BitDepth >::Write( fileName , pixels , _width , _height , 3 );
		delete[] pixels;
	}
	else if constexpr( std::is_same_v< Data , Point3D< CType > > )
	{
		ImageWriter< BitDepth >::Write( fileName , (CType*)_pixels , _width , _height , 3 );
	}
	else Miscellany::ErrorOut( "Bad data type " );

}

#else // !VARIABLE_SIZED_IMAGE_CHANNEL




#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
template< class Data > template < unsigned int BitDepht > void Image< Data >::read( const char* fileName ) { Miscellany::ErrorOut( "Image read not supported" ); }
template< class Data > template < unsigned int BitDepht > void Image< Data >::write( const char* fileName ) const { Miscellany::ErrorOut( "Image write not supported" ); }
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
template< class Data > void Image< Data >::read( const char* fileName ) { Miscellany::ErrorOut( "Image read not supported" ); }
template< class Data > void Image< Data >::write( const char* fileName ) const { Miscellany::ErrorOut( "Image write not supported" ); }
#endif // VARIABLE_SIZED_IMAGE_CHANNEL

template<>
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
template< unsigned int BitDepth >
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
void Image< Point3D< float > >::read( const char* fileName )
{
#ifdef NEW_CODE
	unsigned int width , height;
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
	unsigned char* pixels = ImageReader< BitDepth >::ReadColor( fileName , width , height );
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
	unsigned char* pixels = ImageReader::ReadColor( fileName , width , height );
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
	if( !pixels ) Miscellany::Throw( "Failed to read image: %s\n" , fileName );
	resize( width , height );
	for( int i=0 ; i<(int)width ; i++ ) for( int j=0 ; j<(int)height ; j++ ) for( int c=0 ; c<3 ; c++ ) (*this)(i,j)[c] = ( (float)pixels[ (j*width+i)*3+c ] ) / 255.f;
	delete[] pixels;
#else // !NEW_CODE
	unsigned int width , height , channels;
	unsigned char* pixels = ImageReader::Read( fileName , width , height , channels );
	if( !pixels ) Miscellany::Throw( "Failed to read image: %s\n" , fileName );
	if( channels!=3 ) Miscellany::Throw( "Only three channel images are supported: %d\n" , channels );
	resize( width , height );
	for( int i=0 ; i<(int)width ; i++ ) for( int j=0 ; j<(int)height ; j++ ) for( int c=0 ; c<(int)channels ; c++ ) (*this)(i,j)[c] = ( (float)pixels[ (j*width+i)*3+c ] ) / 255.f;
	delete[] pixels;
#endif // NEW_CODE
}
template<>
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
template< unsigned int BitDepth >
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
void Image< Point3D< double > >::read( const char* fileName )
{
#ifdef NEW_CODE
	unsigned int width , height;
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
	unsigned char* pixels = ImageReader< BitDepth >::ReadColor( fileName , width , height );
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
	unsigned char* pixels = ImageReader::ReadColor( fileName , width , height );
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
	if( !pixels ) Miscellany::Throw( "Failed to read image: %s\n" , fileName );
	resize( width , height );
	for( int i=0 ; i<(int)width ; i++ ) for( int j=0 ; j<(int)height ; j++ ) for( int c=0 ; c<3 ; c++) (*this)(i,j)[c] = ( (double)pixels[ (j*width+i)*3+c ] ) / 255.;
	delete[] pixels;
#else // !NEW_CODE
	unsigned int width , height , channels;
	unsigned char* pixels = ImageReader::Read( fileName , width , height , channels );
	if( !pixels ) Miscellany::Throw( "Failed to read image: %s\n" , fileName );
	if( channels!=3 ) Miscellany::Throw( "Only three channel images are supported: %d\n" , channels );
	resize( width , height );
	for( int i=0 ; i<(int)width ; i++ ) for( int j=0 ; j<(int)height ; j++ ) for( int c=0 ; c<(int)channels ; c++) (*this)(i,j)[c] = ( (double)pixels[ (j*width+i)*3+c ] ) / 255.;
	delete[] pixels;
#endif // NEW_CODE
}

template<>
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
template< unsigned int BitDepth >
void Image< Point3D< typename ImageChannel< BitDepth >::Type > >::read( const char* fileName )
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
void Image< Point3D< unsigned char > >::read( const char* fileName )
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
{
	unsigned int w , h , c;
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
	typename ImageChannel< BitDepth >::Type * pixels = ImageReader< BitDepth >::Read( fileName , w , h , c );
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
	unsigned char* pixels = ImageReader::Read( fileName , w , h , c );
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
	if( !pixels ) Miscellany::Throw( "Failed to read image: %s\n" , fileName );
	if( c!=3 ) Miscellany::Throw( "Only three channel images are supported: %d\n" , c );
	resize( w , h );
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
	memcpy( _pixels , pixels , sizeof( typename ImageChannel< BitDepth >::Type ) * _width * _height*3 );
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
	memcpy( _pixels , pixels , sizeof(unsigned char) * _width * _height*3 );
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
	delete[] pixels;
}

template<>
void Image< Point3D< float > >::write( const char* fileName ) const
{
	unsigned char* pixels = new unsigned char[ _width*_height*3 ];
	for( int i=0 ; i<_width ; i++ ) for( int j=0 ; j<_height ; j++ ) for( int c=0 ; c<3 ; c++ ) pixels[ 3*(j*_width+i)+c ] = (unsigned char)std::min< int >( 255 , std::max< int >( 0 , (int)( (*this)(i,j)[c] * 255.f + 0.5f ) ) );
	ImageWriter::Write( fileName , pixels , _width , _height , 3 );
	delete[] pixels;
}
template<>
void Image< Point3D< double > >::write( const char* fileName ) const
{
	unsigned char* pixels = new unsigned char[ _width*_height*3 ];
	for( int i=0 ; i<_width ; i++ ) for( int j=0 ; j<_height ; j++ ) for( int c=0 ; c<3 ; c++ ) pixels[ 3*(j*_width+i)+c ] = (unsigned char)std::min< int >( 255 , std::max< int >( 0 , (int)( (*this)(i,j)[c] * 255.0 + 0.5 ) ) );
	ImageWriter::Write( fileName , pixels , _width , _height , 3 );
	delete[] pixels;
}
template<>
void Image< Point3D< unsigned char > >::write( const char *fileName ) const { ImageWriter::Write( fileName , (const unsigned char*)_pixels , _width , _height , 3 ); }
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
#endif //IMAGE_INCLUDED