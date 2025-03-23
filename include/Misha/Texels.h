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
		};

		template< typename Real , unsigned int Dim , typename VertexEmbeddingFunctor /* = std::function< Point< double , Dim > ( size_t ) > */ , typename SimplexFunctor /* = std::function< SimplexIndex< K >( size_t ) > */ >
		RegularGrid< K , Point< Real , Dim > > GetTexelPositions( size_t simplexNum , VertexEmbeddingFunctor && VF , SimplexFunctor && SF , const RegularGrid< K , TexelInfo > &texelInfo );

		template< unsigned int Dim , typename VertexEmbeddingFunctor /* = std::function< Point< double , Dim > ( size_t ) > */ , typename SimplexFunctor /* = std::function< SimplexIndex< K >( size_t ) > */ , typename UVCoordinateFunctor /* = std::function< Point< double , K > ( unsigned int , unsigned int ) > */ >
		RegularGrid< K , TexelInfo > GetActiveTexels( size_t simplexNum , VertexEmbeddingFunctor && VF , SimplexFunctor && SF , UVCoordinateFunctor && UV , unsigned int width , unsigned int height , unsigned int dilationRadius , bool finalize , bool verbose );

		template< typename UVCoordinateFunctor /* = std::function< Point< double , K > ( unsigned int , unsigned int ) > */ >
		RegularGrid< K , TexelInfo > GetActiveTexels( size_t simplexNum , UVCoordinateFunctor && UV , unsigned int width , unsigned int height );

		template< typename UVCoordinateFunctor /* = std::function< Point< double , K > ( unsigned int , unsigned int ) > */ >
		RegularGrid< K , TexelInfo > GetInteriorTexels( size_t simplexNum , UVCoordinateFunctor && UV , unsigned int width , unsigned int height );

		template< typename UVCoordinateFunctor /* = std::function< Point< double , K > ( unsigned int , unsigned int ) > */ >
		size_t DilateTexels( UVCoordinateFunctor && UV , RegularGrid< K , TexelInfo > &texelInfo );

		template< unsigned int Dim , typename VertexEmbeddingFunctor /* = std::function< Point< double , Dim > ( size_t ) > */ , typename SimplexFunctor /* = std::function< SimplexIndex< K >( size_t ) > */ >
		void FinalizeTexelInfo( size_t simplexNum , VertexEmbeddingFunctor && VF , SimplexFunctor && SF , RegularGrid< K , TexelInfo > &texelInfo );

#include "Texels.inl"
	}
}
#endif // TEXEL_DILATION_INCLUDED