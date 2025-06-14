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
#ifndef CONSTRAINED_TRIANGULATION_INCLUDED
#define CONSTRAINED_TRIANGULATION_INCLUDED

#include <Misha/Geometry.h>

namespace MishaK
{
	namespace TSP
	{
		template< typename GeometryReal >
		void TriangulateConvexPolygon( const std::vector< Point2D< GeometryReal > > &vertices , std::vector< SimplexIndex< 2 > > &outputTriangles )
		{
			auto SquaredArea = [&]( unsigned int i0 , unsigned int i1 , unsigned int i2 )
				{
					Point< GeometryReal , 2 > d[] = { vertices[i1]-vertices[i0] , vertices[i2]-vertices[i0] };
					SquareMatrix< double , 2 > M;
					for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) M(i,j) = Point< GeometryReal , 2 >::Dot( d[i] , d[j] );
					return M.determinant();
				};
			MinimalAreaTriangulation::GetTriangulation( SquaredArea , static_cast< unsigned int >( vertices.size() ) , outputTriangles );
		}
	}
}
#endif //CONSTRAINED_TRIANGULATION_INCLUDED
