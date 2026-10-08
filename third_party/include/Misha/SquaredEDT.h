/*
Copyright (c) 2019, Michael Kazhdan
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

#ifndef SQUARED_EDT_INCLUDED
#define SQUARED_EDT_INCLUDED

#include <limits>
#include "Array.h"
#include "Geometry.h"
#include "RegularGrid.h"
#include "Rasterizer.h"
#include "MultiThreading.h"

namespace MishaK
{
	template< typename Real , unsigned int Dim >
	struct SquaredEDT
	{
		// [Saito and Toriwaki, '94]
		// New algorithms for Euclidean Distance Transformation of an n-dimensional digitized picture with applications
		template< typename BinaryType >
		static RegularGrid< Dim , unsigned int > Saito( const RegularGrid< Dim , BinaryType > &binaryGrid , bool verbose=false );
		template< typename BinaryType >
		static RegularGrid< Dim , std::pair< unsigned int , size_t > > FullSaito( const RegularGrid< Dim , BinaryType > &binaryGrid , bool verbose=false );
		template< typename IndexType , unsigned int K >
		static RegularGrid< Dim , unsigned int > Saito( const SimplicialComplex< Real , Dim , K > &simplicialComplex , unsigned int depth , unsigned int lockDepth , XForm< Real , Dim+1 > &gridToModel , Real bBoxScale , bool verbose=false );


		// [Danielsson, '90]
		// Euclidean distance mapping
		template< typename IndexType , unsigned int K >
		static RegularGrid< Dim , Real > Danielsson( const SimplicialComplex< Real , Dim , K > &simplicialComplex , unsigned int depth , unsigned int lockDepth , unsigned int radius , XForm< Real , Dim+1 > &gridToModel , Real bBoxScale , bool verbose=false );
	protected:
		template< typename IndexType , unsigned int K , unsigned int SliceDim , bool MultiThreaded >
		static typename std::enable_if< SliceDim!=0 >::type _Danielsson( const std::vector< typename Simplex< Real , Dim , K >::NearestKey > &nearestKeys , Pointer( std::pair< Real , IndexType > ) sliceNearest , const unsigned int res[SliceDim] , Point< Real , Dim > center );
		template< typename IndexType , unsigned int K , unsigned int SliceDim , bool MultiThreaded >
		static typename std::enable_if< SliceDim==0 >::type _Danielsson( const std::vector< typename Simplex< Real , Dim , K >::NearestKey > &nearestKeys , Pointer( std::pair< Real , IndexType > ) sliceNearest , const unsigned int res[] , Point< Real , Dim > center ){}
	};
#include "SquaredEDT.inl"
}
#endif // SQUARED_EDT_INCLUDED
