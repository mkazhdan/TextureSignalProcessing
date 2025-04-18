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

#define NEW_TEXEL_CODE

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
		const unsigned int Dim = 2;

		template< typename SimplexIndexFunctor , unsigned int K >
		constexpr bool IsValidSimplexIndexFunctor( void ){ return std::is_convertible_v< SimplexIndexFunctor , std::function< SimplexIndex< K >( size_t ) > >; }
		template< typename SimplexFunctor , unsigned int EmbeddingDim , unsigned int K >
		constexpr bool IsValidSimplexFunctor( void ){ return std::is_convertible_v< SimplexFunctor , std::function< Simplex< double , EmbeddingDim , K >( size_t ) > >; }
		template< typename VertexFunctor , unsigned int EmbeddingDim >
		constexpr bool IsValidVertexFunctor( void ){ return std::is_convertible_v< VertexFunctor , std::function< Point< double , EmbeddingDim >( size_t ) > >; }

		struct TexelInfo
		{
			TexelInfo( void ) : sIdx(-1) { for( unsigned int d=0 ; d<=Dim ; d++ ) bc[d] = 1./(Dim+1); }

			unsigned int sIdx;
			Point< double , Dim+1 > bc;

			template< typename Real , unsigned int EmbeddingDim , typename SimplexEmbeddingFunctor /* = std::function< Point< double , EmbeddingDim > ( size_t ) > */ >
			Point< Real , EmbeddingDim > position( const SimplexEmbeddingFunctor & SEF ) const;
		};

		template< bool NodeAtCellCenter >
		Point< double , Dim > TexelNodePosition( typename RegularGrid< Dim >::Index I );

		template< bool NodeAtCellCenter , unsigned int K >
		Simplex< double , Dim , K > GetTexelSpaceSimplex( Simplex< double , Dim , K > simplex , const unsigned int res[Dim] );

		template< typename Real , unsigned int EmbeddingDim , typename SimplexEmbeddingFunctor /* = std::function< Simplex< double , EmbeddingDim , Dim > > ( size_t ) > */ >
		RegularGrid< Dim , Point< Real , EmbeddingDim > > GetTexelPositions( size_t simplexNum , SimplexEmbeddingFunctor && SEF , const RegularGrid< Dim , TexelInfo > &texelInfo );

#pragma message ( "[WARNING] The output should really be a grid of vectors, since multiple triangles may be supported by the same texel" )
		template< unsigned int EmbeddingDim , bool Nearest , bool NodeAtCellCenter , typename VertexEmbeddingFunctor /* = std::function< Point< double , EmbeddingDim > ( size_t ) > */ , typename SimplexIndexFunctor /* = std::function< SimplexIndex< Dim >( size_t ) > */ , typename SimplexFunctor /* = std::function< Simplex< double , Dim , Dim > ( size_t ) > */ >
		RegularGrid< Dim , TexelInfo > GetSupportedTexelInfo( size_t simplexNum , VertexEmbeddingFunctor && VF , SimplexIndexFunctor && SIF , SimplexFunctor && SF , const unsigned int res[Dim] , unsigned int dilationRadius , bool finalize , bool verbose );

		template< bool Nearest , bool NodeAtCellCenter , typename SimplexFunctor /* = std::function< Simplex< double , Dim , Dim > ( size_t ) > */ >
		RegularGrid< Dim , TexelInfo > GetSupportedTexelInfo( size_t simplexNum , SimplexFunctor && SF , const unsigned int res[Dim] );

		template< bool Nearest , bool NodeAtCellCenter , unsigned int K=Dim , typename SimplexFunctor = std::function< Simplex< double , Dim , K > ( size_t ) > >
		RegularGrid< Dim , std::vector< size_t > > GetSupportedSimplexIndices( size_t simplexNum , SimplexFunctor && SF , unsigned int res[Dim] );

		template< bool NodeAtCellCenter , typename SimplexFunctor /* = std::function< Simplex< double , Dim , Dim > ( size_t ) > */ >
		RegularGrid< Dim , TexelInfo > GetNodeTexelInfo( size_t simplexNum , SimplexFunctor && SF , const unsigned int res[Dim] , bool forceThreadSafe=false );

		template< bool Nearest , bool NodeAtCellCenter , typename SimplexFunctor /* = std::function< Simplex< double , Dim , Dim > ( size_t ) > */ >
		size_t DilateTexelInfo( SimplexFunctor && SF , RegularGrid< Dim , TexelInfo > &texelInfo );

		template< unsigned int EmbeddingDim , typename VertexEmbeddingFunctor /* = std::function< Point< double , EmbeddingDim > ( size_t ) > */ , typename SimplexIndexFunctor /* = std::function< SimplexIndex< Dim >( size_t ) > */ >
		void FinalizeTexelInfo( size_t simplexNum , VertexEmbeddingFunctor && VF , SimplexIndexFunctor && SIF , RegularGrid< Dim , TexelInfo > &texelInfo );

#include "Texels.inl"
	}
}
#endif // TEXEL_DILATION_INCLUDED