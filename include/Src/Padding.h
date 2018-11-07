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
#pragma once

class Padding {
public:
	Padding() {
		left = bottom = right = top = 0;
		nonTrivial = false;
	}
	int left;
	int bottom;
	int right;
	int top;
	bool nonTrivial;
};

void ComputePadding( Padding& padding , const int& width , const int& height , std::vector< Point2D< double > >& textureCoordinates , bool verbose )
{
	Point2D<double> pixMinCorner(0.5 / double(width), 0.5 / double(height));
	Point2D<double> pixMaxCorner((double(width) - 0.5) / double(width), (double(height) - 0.5) / double(height));

	Point2D<double> texMinCorner = textureCoordinates[0];
	Point2D<double> texMaxCorner = textureCoordinates[0];
	for (int i = 0; i < textureCoordinates.size(); i++) {
		for (int c = 0; c < 2;c++) texMinCorner[c] = std::min<double>(texMinCorner[c], textureCoordinates[i][c]);
		for (int c = 0; c < 2;c++) texMaxCorner[c] = std::max<double>(texMaxCorner[c], textureCoordinates[i][c]);
	}

	padding.left = texMinCorner[0] < pixMinCorner[0] ? (int)ceil((pixMinCorner[0] - texMinCorner[0]) * double(width)) : 0;
	padding.bottom = texMinCorner[1] < pixMinCorner[1] ? (int)ceil((pixMinCorner[1] - texMinCorner[1]) * double(height)) : 0;

	padding.right = texMaxCorner[0] > pixMaxCorner[0] ? (int)ceil((texMaxCorner[0] - pixMaxCorner[0]) * double(width)) : 0;
	padding.top = texMaxCorner[1] > pixMaxCorner[1] ? (int)ceil((texMaxCorner[1] - pixMaxCorner[1]) * double(height)) : 0;

	//Make image dimensions are multiples of 8 (Hardware texture mapping seems to fail if not)
	int newWidth = width + padding.left + padding.right;
	int newHeight = height + padding.bottom + padding.top;

	int paddedWidth = 8 * (((newWidth - 1) / 8) + 1);
	int paddedHeight = 8 * (((newHeight - 1) / 8) + 1);
	padding.left += (paddedWidth - newWidth);
	padding.bottom += (paddedHeight - newHeight);

	if( padding.left || padding.right || padding.bottom || padding.top ) padding.nonTrivial = true;
	if( padding.nonTrivial )
	{
		if( verbose ) printf( "Padding applied : Left %d. Right %d. Bottom %d. Top %d.\n" , padding.left , padding.right , padding.bottom , padding.top );
	}
	else
	{
		if( verbose ) printf( "No padding required!\n" );
	}

}

template< class DataType >
void PadImage( Padding &padding , Image< DataType > &im )
{
	int newWidth = im.width() + padding.left + padding.right;
	int newHeight = im.height() + padding.bottom + padding.top;

	Image< DataType > newIm;
	newIm.resize( newWidth , newHeight );
	for( int i=0 ; i<newWidth ; i++ ) for( int j=0 ; j<newHeight ; j++ )
	{
		int ni = std::min< int >( std::max< int >( 0 , i-padding.left   ) , im.width()  - 1 );
		int nj = std::min< int >( std::max< int >( 0 , j-padding.bottom ) , im.height() - 1 );
		newIm(i,j) = im(ni,nj);
	}
	im = newIm;
}

template<class DataType>
void UnpadImage( Padding & padding , Image<DataType> & im )
{
	int outputWidth = im.width() - padding.left - padding.right;
	int outputHeight = im.height() - padding.bottom - padding.top;
	Image<Point3D<float>> newIm;
	newIm.resize(outputWidth, outputHeight);
	for (int i = 0; i < outputWidth; i++)for (int j = 0; j < outputHeight; j++)newIm(i, j) = im(padding.left + i, padding.bottom + j);
	im = newIm;
}


void PadTextureCoordinates(Padding & padding, int width, int height, std::vector<Point2D<double>>  & textureCoordinates){
	int newWidth = width + padding.left + padding.right;
	int newHeight = height + padding.bottom + padding.top;

	for (int i = 0; i < textureCoordinates.size(); i++) {
		textureCoordinates[i][0] = (textureCoordinates[i][0] * double(width) + double(padding.left)) / double(newWidth);
		textureCoordinates[i][1] = (textureCoordinates[i][1] * double(height) + double(padding.bottom)) / double(newHeight);
	}
}

void UnpadTextureCoordinates( Padding& padding , int width , int height , std::vector< Point2D< double > >& textureCoordinates )
{
	int newWidth = width + padding.left + padding.right;
	int newHeight = height + padding.bottom + padding.top;

	for (int i = 0; i < textureCoordinates.size(); i++) {
		textureCoordinates[i][0] = (textureCoordinates[i][0] * double(newWidth) - double(padding.left)) / double(width);
		textureCoordinates[i][1] = (textureCoordinates[i][1] * double(newHeight) - double(padding.bottom)) / double(height);
	}
}