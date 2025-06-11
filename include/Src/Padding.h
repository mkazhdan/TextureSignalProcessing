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

namespace MishaK
{
	namespace TSP
	{
		// A structure for computing and (un)setting the padding needed to ensure that texture coordinates fall within the rectangle defined by the _centers_ of the corner texels.
		// Given a texture coordinate (s) indexing a texture map of width W, we have:
		//		s -> W*s
		// Offsetting W -> W+D we want the associated texture coordinate t to satisfy:
		//		(W+D)*t = D + W*s
		//		t = ( W*s + D ) / (W+D)
		struct Padding
		{
			unsigned int left , right , bottom , top;

			Padding( void ) : left(0) , bottom(0) , right(0) , top(0){}

			unsigned int width( void ) const { return left+right; }
			unsigned int height( void ) const { return bottom+top; }

			template< typename GeometryReal >
			static Padding Init( unsigned int width , unsigned int height , const std::vector< Point2D< GeometryReal > > &textureCoordinates , bool verbose=false )
			{
				Padding padding;
				Point2D< GeometryReal > pixMinCorner( (GeometryReal)0.5/width , (GeometryReal)0.5/height );
				Point2D< GeometryReal > pixMaxCorner( (GeometryReal)( width-0.5 )/ width , (GeometryReal)( height-0.5 )/height );

				Point2D< GeometryReal > texMinCorner = textureCoordinates[0];
				Point2D< GeometryReal > texMaxCorner = textureCoordinates[0];
				for( int i=0 ; i<textureCoordinates.size() ; i++ )
				{
					for( int c=0 ; c<2 ; c++ ) texMinCorner[c] = std::min< GeometryReal >( texMinCorner[c] , textureCoordinates[i][c] );
					for( int c=0 ; c<2 ; c++ ) texMaxCorner[c] = std::max< GeometryReal >( texMaxCorner[c] , textureCoordinates[i][c] );
				}

				padding.left   = texMinCorner[0] < pixMinCorner[0] ? (int)ceil( ( pixMinCorner[0]-texMinCorner[0] )*width  ) : 0;
				padding.bottom = texMinCorner[1] < pixMinCorner[1] ? (int)ceil( ( pixMinCorner[1]-texMinCorner[1] )*height ) : 0;

				padding.right = texMaxCorner[0] > pixMaxCorner[0] ? (int)ceil( ( texMaxCorner[0]-pixMaxCorner[0] )*width  ) : 0;
				padding.top   = texMaxCorner[1] > pixMaxCorner[1] ? (int)ceil( ( texMaxCorner[1]-pixMaxCorner[1] )*height ) : 0;

				// Make image dimensions multiples of 8 (Hardware texture mapping seems to fail if not)
				{
					int newWidth = width + padding.left + padding.right;
					int newHeight = height + padding.bottom + padding.top;

					int paddedWidth = 8 * (((newWidth - 1) / 8) + 1);
					int paddedHeight = 8 * (((newHeight - 1) / 8) + 1);
					padding.left += (paddedWidth - newWidth);
					padding.bottom += (paddedHeight - newHeight);
				}

				if( verbose )
				{
					if( padding.width() || padding.height() ) printf( "Padding applied : Left %d. Right %d. Bottom %d. Top %d.\n" , padding.left , padding.right , padding.bottom , padding.top );
					else                                      printf( "No padding required!\n" );
				}

				return padding;
			}

			// Add the padding to an image (set new texelv alues to closest boundary texel)
			// [WARNING] Assuming the image dimensions match those used to define the object
			template< typename DataType >
			void pad( Image< DataType > &im ) const
			{
				if( !( left || right || bottom || top ) ) return;

				unsigned int newWidth = im.res(0) + left + right;
				unsigned int newHeight = im.res(1) + bottom + top;

				Image< DataType > newIm;
				newIm.resize( newWidth , newHeight );
				for( unsigned int i=0 ; i<newWidth ; i++ ) for( unsigned int j=0 ; j<newHeight ; j++ )
				{
					unsigned int ni = std::min< int >( std::max< int >( 0 , (int)i-(int)left   ) , im.res(0) - 1 );
					unsigned int nj = std::min< int >( std::max< int >( 0 , (int)j-(int)bottom ) , im.res(1) - 1 );
					newIm(i,j) = im(ni,nj);
				}
				im = newIm;
			}

			// Remove the padding from an image
			template< class DataType >
			void unpad( Image<DataType> &im ) const
			{
				if( !( left || right || bottom || top ) ) return;

				unsigned int outputWidth = im.res(0) - left - right;
				unsigned int outputHeight = im.res(1) - bottom - top;
				Image< DataType > newIm;
				newIm.resize( outputWidth , outputHeight );
				for( unsigned int i=0 ; i<outputWidth ; i++ ) for( unsigned int j=0 ; j<outputHeight ; j++ ) newIm(i,j) = im( left+i , bottom+j );
				im = newIm;
			}

			template< typename GeometryReal >
			void pad( int width , int height , std::vector< Point2D< GeometryReal > > &textureCoordinates ) const
			{
				if( !( left || right || bottom || top ) ) return;

				int newWidth = width + left + right;
				int newHeight = height + bottom + top;

				for( int i=0 ; i<textureCoordinates.size() ; i++ )
				{
					textureCoordinates[i][0] = ( textureCoordinates[i][0]*width  + (GeometryReal)( left   ) )/newWidth;
					textureCoordinates[i][1] = ( textureCoordinates[i][1]*height + (GeometryReal)( bottom ) )/newHeight;
				}
			}

			template< typename GeometryReal >
			void unpad( int width , int height , std::vector< Point2D< GeometryReal > > &textureCoordinates ) const
			{
				if( !( left || right || bottom || top ) ) return;

				int newWidth = width + left + right;
				int newHeight = height + bottom + top;

				for( int i=0 ; i<textureCoordinates.size() ; i++ )
				{
					textureCoordinates[i][0] = ( textureCoordinates[i][0]*newWidth  - (GeometryReal)( left   ) )/width;
					textureCoordinates[i][1] = ( textureCoordinates[i][1]*newHeight - (GeometryReal)( bottom ) )/height;
				}
			}
		};
	}
}