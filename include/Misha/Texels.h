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
#include <map>
#include <set>

#include "Misha/Exceptions.h"
#include "Misha/FEM.h"
#include "Misha/Rasterizer2D.h"
#include "Misha/RegularGrid.h"
#include "Misha/Miscellany.h"
#include "Misha/MultiThreading.h"

namespace MishaK
{
	template< bool NodeAtCellCenter , typename Index=size_t , unsigned int Dim=2 >
	struct Texels
	{
		static_assert( Dim==2 , "[ERROR] Only supported for triangle meshes" );

		template< unsigned int K >
		struct TexelInfo
		{
			TexelInfo( void ) : sIdx(-1) { for( unsigned int k=0 ; k<=K ; k++ ) bc[k] = 1./(K+1); }

			Index sIdx;
			Point< double , K+1 > bc;

			template< unsigned int EmbeddingDim , typename SimplexEmbeddingFunctor /* = std::function< Simplex< double , EmbeddingDim , K > ( size_t ) > */ >
			Point< double , EmbeddingDim >
				position( const SimplexEmbeddingFunctor & SEF ) const;
		};

		static Point< double , Dim > NodePosition( typename RegularGrid< Dim >::Index I );

		template< unsigned int K >
		static Simplex< double , Dim , K > TexelSimplex( Simplex< double , Dim , K > simplex , const unsigned int res[/*Dim*/] );


		// Returns a grid with each texel storing the indices of the simplices in its support
		template< bool Nearest , unsigned int K=Dim , typename SimplexFunctor = std::function< Simplex< double , Dim , K > ( size_t ) > >
		static RegularGrid< Dim , std::vector< Index > > GetSupportedSimplexIndices( size_t simplexNum , SimplexFunctor && SF , const unsigned int res[/*Dim*/] );

		// Returns a grid with each texel storing the indices of simplices covering the node
		template< typename SimplexFunctor /* = std::function< Simplex< double , Dim , Dim > ( size_t ) > */ >
		static RegularGrid< Dim , std::vector< Index > > GetNodeSimplexIndices( size_t simplexNum , SimplexFunctor && SF , const unsigned int res[/*Dim*/] );

		// Flood-fills simplex index assignments into neighboring texels
		static size_t DilateSimplexIndices( RegularGrid< Dim , std::vector< Index > > &tIndices , unsigned int radius , bool verbose=false );

		// Converts a grid with multiple simplex indicces per texel into a grid with the closest one
		template< unsigned int K=Dim , typename SimplexFunctor = std::function< Simplex< double , Dim , K > ( size_t ) > >
		static RegularGrid< Dim , TexelInfo< K > > GetNearestTexelInfo( SimplexFunctor && SF , const RegularGrid< Dim , std::vector< Index > > & tIndices );



		// Returns a grid with each texel storing the index of a simplex closest to the node
		template< bool Nearest , unsigned int K=Dim , typename SimplexFunctor = std::function< Simplex< double , Dim , K > ( size_t ) > >
		static RegularGrid< Dim , TexelInfo< K > > GetSupportedTexelInfo( size_t simplexNum , SimplexFunctor && SF , const unsigned int res[/*Dim*/] );

		// Returns a grid with each texel storing the index of a simplex covering the node
		template< typename SimplexFunctor /* = std::function< Simplex< double , Dim , Dim > ( size_t ) > */ >
		static RegularGrid< Dim , TexelInfo< Dim > > GetNodeTexelInfo( size_t simplexNum , SimplexFunctor && SF , const unsigned int res[/*Dim*/] );

		// Flood-fills simplex index assignments into neighboring texels
		template< unsigned int K , typename SimplexFunctor /* = std::function< Simplex< double , Dim , K > ( size_t ) > */ >
		static size_t DilateTexelInfo( SimplexFunctor && SF , RegularGrid< Dim , TexelInfo< K > > &tInfo , unsigned int radius , bool verbose=false );


		// Returns a grid with each texel storing the index of a simplex closest to the node and dilates
		template< bool Nearest , unsigned int K=Dim , bool UseVector=true , typename SimplexFunctor = std::function< Simplex< double , Dim , K > ( size_t ) > >
		static RegularGrid< Dim , TexelInfo< K > > GetSupportedTexelInfo( size_t simplexNum , SimplexFunctor && SF , const unsigned int res[/*Dim*/] , unsigned int dilationRadius , bool verbose );

		// Returns a grid with each texel storing the index of a simplex covering the node and dilates
		template< bool UseVector=true , typename SimplexFunctor = std::function< Simplex< double , Dim , Dim > ( size_t ) > >
		static RegularGrid< Dim , TexelInfo< Dim > > GetNodeTexelInfo( size_t simplexNum , SimplexFunctor && SF , const unsigned int res[/*Dim*/] , unsigned int dilationRadius , bool verbose );


		// Returns a grid with each texel storing the index of a simplex closest to the node, dilates and flows to the interior
		template< unsigned int EmbeddingDim , bool Nearest , unsigned int K=Dim , bool UseVector=true , typename VertexEmbeddingFunctor = std::function< Point< double , EmbeddingDim > ( size_t ) > , typename SimplexIndexFunctor = std::function< SimplexIndex< Dim >( size_t ) > , typename SimplexFunctor = std::function< Simplex< double , Dim , K > ( size_t ) > >
		static RegularGrid< Dim , TexelInfo< K > > GetSupportedTexelInfo( size_t simplexNum , VertexEmbeddingFunctor && VF , SimplexIndexFunctor && SIF , SimplexFunctor && SF , const unsigned int res[/*Dim*/] , unsigned int dilationRadius , bool verbose );

		// Returns a grid with each texel storing the index of a simplex covering the node, dilates, and flows to the interior
		template< unsigned int EmbeddingDim , bool UseVector=true , typename VertexEmbeddingFunctor = std::function< Point< double , EmbeddingDim > ( size_t ) > , typename SimplexIndexFunctor = std::function< SimplexIndex< Dim >( size_t ) > , typename SimplexFunctor = std::function< Simplex< double , Dim , Dim > ( size_t ) > >
		static RegularGrid< Dim , TexelInfo< Dim > > GetNodeTexelInfo( size_t simplexNum , VertexEmbeddingFunctor && VF , SimplexIndexFunctor && SIF , SimplexFunctor && SF , const unsigned int res[/*Dim*/] , unsigned int dilationRadius , bool verbose );


		template< unsigned int EmbeddingDim , typename VertexEmbeddingFunctor /* = std::function< Point< double , EmbeddingDim > ( size_t ) > */ , typename SimplexIndexFunctor /* = std::function< SimplexIndex< Dim >( size_t ) > */ >
		static void FlowTexelInfoToInterior( size_t simplexNum , VertexEmbeddingFunctor && VF , SimplexIndexFunctor && SIF , RegularGrid< Dim , TexelInfo< Dim > > &texelInfo );

		// Compute embedding positions for active texels
		template< typename Real , unsigned int EmbeddingDim , unsigned int K , typename SimplexEmbeddingFunctor /* = std::function< Simplex< double , EmbeddingDim , K > > ( size_t ) > */ >
		static RegularGrid< Dim , Point< Real , EmbeddingDim > > GetTexelPositions( size_t simplexNum , SimplexEmbeddingFunctor && SEF , const RegularGrid< Dim , TexelInfo< K > > &texelInfo );

	protected:
		template< typename SimplexIndexFunctor , unsigned int K >
		static constexpr bool _IsValidSimplexIndexFunctor( void ){ return std::is_convertible_v< SimplexIndexFunctor , std::function< SimplexIndex< K >( size_t ) > >; }

		template< typename SimplexFunctor , unsigned int EmbeddingDim , unsigned int K >
		static constexpr bool _IsValidSimplexFunctor( void ){ return std::is_convertible_v< SimplexFunctor , std::function< Simplex< double , EmbeddingDim , K >( size_t ) > >; }

		template< typename VertexFunctor , unsigned int EmbeddingDim >
		static constexpr bool _IsValidVertexFunctor( void ){ return std::is_convertible_v< VertexFunctor , std::function< Point< double , EmbeddingDim >( size_t ) > >; }
	};
#include "Texels.inl"
}
#endif // TEXEL_DILATION_INCLUDED