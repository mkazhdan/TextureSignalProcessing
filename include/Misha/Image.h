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
	void read(const char* fileName);
	void write(const char* fileName) const;
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
template< class Data > void Image< Data >::read( const char* fileName ) { Miscellany::ErrorOut( "Image read not supported" ); }
template< class Data > void Image< Data >::write( const char* fileName ) const { Miscellany::ErrorOut( "Image write not supported" ); }

template<>
void Image< Point3D< float > >::read( const char* fileName )
{
	unsigned int width , height , channels;
	unsigned char* pixels = ImageReader::Read( fileName , width , height , channels );
	if( !pixels ) Miscellany::Throw( "Failed to read image: %s\n" , fileName );
	if( channels!=3 ) Miscellany::Throw( "Only three channel images are supported: %d\n" , channels );
	resize( width , height );
	for( int i=0 ; i<(int)width ; i++ ) for( int j=0 ; j<(int)height ; j++ ) for( int c=0 ; c<(int)channels ; c++ ) (*this)(i,j)[c] = ( (float)pixels[ (j*width+i)*3+c ] ) / 255.f;
	delete[] pixels;
}
template<>
void Image< Point3D< double > >::read( const char* fileName )
{
	unsigned int width , height , channels;
	unsigned char* pixels = ImageReader::Read( fileName , width , height , channels );
	if( !pixels ) Miscellany::Throw( "Failed to read image: %s\n" , fileName );
	if( channels!=3 ) Miscellany::Throw( "Only three channel images are supported: %d\n" , channels );
	resize( width , height );
	for( int i=0 ; i<(int)width ; i++ ) for( int j=0 ; j<(int)height ; j++ ) for( int c=0 ; c<(int)channels ; c++) (*this)(i,j)[c] = ( (double)pixels[ (j*width+i)*3+c ] ) / 255.;
	delete[] pixels;
}
template<>
void Image< Point3D< unsigned char > >::read( const char* fileName )
{
	unsigned int w , h , c;
	unsigned char* pixels = ImageReader::Read( fileName , w , h , c );
	if( !pixels ) Miscellany::Throw( "Failed to read image: %s\n" , fileName );
	if( c!=3 ) Miscellany::Throw( "Only three channel images are supported: %d\n" , c );
	resize( w , h );
	memcpy( _pixels , pixels , sizeof(unsigned char) * _width * _height*3 );
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

#endif //IMAGE_INCLUDED