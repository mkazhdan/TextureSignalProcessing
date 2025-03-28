/*
Copyright (c) 2025, Michael Kazhdan
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

#ifndef TEXEL_DILATION_INCLUDED
#define TEXEL_DILATION_INCLUDED

#include <vector>

#include "Misha/Exceptions.h"
#include "Misha/FEM.h"
#include "Misha/Rasterizer2D.h"
#include "Misha/RegularGrid.h"
#include "Misha/Miscellany.h"
#include "Misha/MultiThreading.h"

namespace MishaK
{
	namespace Texels
	{
		const unsigned int K = 2;

		struct TexelInfo
		{
			TexelInfo( void ) : sIdx(-1) { for( unsigned int k=0 ; k<=K ; k++ ) bc[k] = 1./(K+1); }

			unsigned int sIdx;
			Point< double , K+1 > bc;

			template< typename Real , unsigned int Dim , typename VertexEmbeddingFunctor /* = std::function< Point< double , Dim > ( size_t ) > */ , typename SimplexFunctor /* = std::function< SimplexIndex< K >( size_t ) > */ >
			Point< Real , Dim > position( const VertexEmbeddingFunctor & VF , const SimplexFunctor & SF ) const;
		};

		template< typename Real , unsigned int Dim , typename VertexEmbeddingFunctor /* = std::function< Point< double , Dim > ( size_t ) > */ , typename SimplexFunctor /* = std::function< SimplexIndex< K >( size_t ) > */ >
		RegularGrid< K , Point< Real , Dim > > GetTexelPositions( size_t simplexNum , VertexEmbeddingFunctor && VF , SimplexFunctor && SF , const RegularGrid< K , TexelInfo > &texelInfo );

#pragma message ( "[WARNING] The output should really be a grid of vectors, since multiple triangles may be supported by the same texel" )
		template< unsigned int Dim , bool Nearest , bool NodeAtCellCenter , typename VertexEmbeddingFunctor /* = std::function< Point< double , Dim > ( size_t ) > */ , typename SimplexFunctor /* = std::function< SimplexIndex< K >( size_t ) > */ , typename UVCoordinateFunctor /* = std::function< Point< double , K > ( unsigned int , unsigned int ) > */ >
		RegularGrid< K , TexelInfo > GetSupportedTexelInfo( size_t simplexNum , VertexEmbeddingFunctor && VF , SimplexFunctor && SF , UVCoordinateFunctor && UV , unsigned int width , unsigned int height , unsigned int dilationRadius , bool finalize , bool verbose );

		template< bool Nearest , bool NodeAtCellCenter , typename UVCoordinateFunctor /* = std::function< Point< double , K > ( unsigned int , unsigned int ) > */ >
		RegularGrid< K , TexelInfo > GetSupportedTexelInfo( size_t simplexNum , UVCoordinateFunctor && UV , unsigned int width , unsigned int height );

		template< bool NodeAtCellCenter , typename UVCoordinateFunctor /* = std::function< Point< double , K > ( unsigned int , unsigned int ) > */ >
		RegularGrid< K , TexelInfo > GetNodeTexelInfo( size_t simplexNum , UVCoordinateFunctor && UV , unsigned int width , unsigned int height );

		template< bool Nearest , bool NodeAtCellCenter , typename UVCoordinateFunctor /* = std::function< Point< double , K > ( unsigned int , unsigned int ) > */ >
		size_t DilateTexelInfo( UVCoordinateFunctor && UV , RegularGrid< K , TexelInfo > &texelInfo );
		template< unsigned int Dim , typename VertexEmbeddingFunctor /* = std::function< Point< double , Dim > ( size_t ) > */ , typename SimplexFunctor /* = std::function< SimplexIndex< K >( size_t ) > */ >
		void FinalizeTexelInfo( size_t simplexNum , VertexEmbeddingFunctor && VF , SimplexFunctor && SF , RegularGrid< K , TexelInfo > &texelInfo );

#include "Texels.inl"
	}
}
#endif // TEXEL_DILATION_INCLUDED