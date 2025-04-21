/*
Copyright (c) 2025, Fabian Prada and Michael Kazhdan
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
#ifdef DEBUG_INDEXING
	template< typename T >
	struct UnsignedIntIndex
	{
		explicit UnsignedIntIndex( unsigned int idx=-1 ) : _idx(idx){}
		explicit operator unsigned int (){ return _idx; }
	protected:
		unsigned int _idx;
	};

	struct               ChartIndex : public UnsignedIntIndex<               ChartIndex >{ using UnsignedIntIndex<               ChartIndex >::UnsignedIntIndex; };
	struct       AtlasTriangleIndex : public UnsignedIntIndex<       AtlasTriangleIndex >{ using UnsignedIntIndex<       AtlasTriangleIndex >::UnsignedIntIndex; };
	struct       ChartTriangleIndex : public UnsignedIntIndex<       ChartTriangleIndex >{ using UnsignedIntIndex<       ChartTriangleIndex >::UnsignedIntIndex; };
	struct         AtlasVertexIndex : public UnsignedIntIndex<         AtlasVertexIndex >{ using UnsignedIntIndex<         AtlasVertexIndex >::UnsignedIntIndex; };
	struct         ChartVertexIndex : public UnsignedIntIndex<         ChartVertexIndex >{ using UnsignedIntIndex<         ChartVertexIndex >::UnsignedIntIndex; };
	struct           AtlasEdgeIndex : public UnsignedIntIndex<           AtlasEdgeIndex >{ using UnsignedIntIndex<           AtlasEdgeIndex >::UnsignedIntIndex; };
	struct       AtlasHalfEdgeIndex : public UnsignedIntIndex<       AtlasHalfEdgeIndex >{ using UnsignedIntIndex<       AtlasHalfEdgeIndex >::UnsignedIntIndex; };
	struct       ChartHalfEdgeIndex : public UnsignedIntIndex<       ChartHalfEdgeIndex >{ using UnsignedIntIndex<       ChartHalfEdgeIndex >::UnsignedIntIndex; };
	struct AtlasBoundaryVertexIndex : public UnsignedIntIndex< AtlasBoundaryVertexIndex >{ using UnsignedIntIndex< AtlasBoundaryVertexIndex >::UnsignedIntIndex; };

#else // !DEBUG_INDEXING
	using               ChartIndex = unsigned int;
	using       AtlasTriangleIndex = unsigned int;
	using       ChartTriangleIndex = unsigned int;
	using         AtlasVertexIndex = unsigned int;
	using         ChartVertexIndex = unsigned int;
	using           AtlasEdgeIndex = unsigned int;
	using       ChartHalfEdgeIndex = unsigned int;
	using       AtlasHalfEdgeIndex = unsigned int;
	using AtlasBoundaryVertexIndex = unsigned int;
#endif // DEBUG_INDEXING
}