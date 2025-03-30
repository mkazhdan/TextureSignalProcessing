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

#ifndef RASTERIZER_2D_INCLUDED
#define RASTERIZER_2D_INCLUDED

#define NEW_RASTERIZER_CODE

#include "RegularGrid.h"
#include "Geometry.h"
#include "Miscellany.h"

namespace MishaK
{
	// Nearest: Specifies whether constant or bi-linear interpolation is to be used
	// Note: If the mapping to 2D is injective, RasterizeNodes should be parallelizable (assuming simplices to not pass exactly through node centers)
	namespace Rasterizer2D
	{
		const unsigned int Dim = 2;

		using Index = typename RegularGrid< Dim >::Index;
		template< unsigned int Dim > using Range = typename RegularGrid< Dim >::Range;

		template< bool NodeAtCellCenter , typename RasterizationFunctor /* = std::function< void ( Index ) > )*/ >
		void RasterizeNodes( Simplex< double , Dim , Dim > triangle , RasterizationFunctor && F , Range< Dim > cellRange );

		template< bool Nearest , bool NodeAtCellCenter , unsigned int K , typename RasterizationFunctor /* = std::function< void ( Index ) > )*/ >
		void RasterizeSupports( Simplex< double , Dim , K > simplex , RasterizationFunctor && F , Range< Dim > cellRange );

#include "Rasterizer2D.inl"
	};
}
#endif // RASTERIZER_2D_INCLUDED