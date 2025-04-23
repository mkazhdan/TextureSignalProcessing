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

#include <optional>

namespace MishaK
{
#ifdef DEBUG_INDEXING
	template< typename T >
	struct UnsignedIntIndex
	{
		explicit UnsignedIntIndex( unsigned int idx=-1 ) : _idx(idx){}
		explicit operator unsigned int () const { return _idx; }

		UnsignedIntIndex & operator += ( unsigned int off ){ _idx += off ; return *this; }
		UnsignedIntIndex operator + ( unsigned int off ) const { return UnsignedIntIndex( _idx+off ); }

		bool operator == ( const UnsignedIntIndex &I ) const { return _idx==I._idx; }
		bool operator != ( const UnsignedIntIndex &I ) const { return _idx!=I._idx; }
		bool operator <  ( const UnsignedIntIndex &I ) const { return _idx< I._idx; }

		friend std::ostream &operator << ( std::ostream &os , const UnsignedIntIndex &I ){ return os << I._idx; }
	protected:
		unsigned int _idx;
	};

	struct               ChartIndex : public UnsignedIntIndex<               ChartIndex >{ using UnsignedIntIndex<               ChartIndex >::UnsignedIntIndex; };
	struct       AtlasTriangleIndex : public UnsignedIntIndex<       AtlasTriangleIndex >{ using UnsignedIntIndex<       AtlasTriangleIndex >::UnsignedIntIndex; };
	struct       ChartTriangleIndex : public UnsignedIntIndex<       ChartTriangleIndex >{ using UnsignedIntIndex<       ChartTriangleIndex >::UnsignedIntIndex; };
	struct         AtlasVertexIndex : public UnsignedIntIndex<         AtlasVertexIndex >{ using UnsignedIntIndex<         AtlasVertexIndex >::UnsignedIntIndex; };
	struct         ChartVertexIndex : public UnsignedIntIndex<         ChartVertexIndex >{ using UnsignedIntIndex<         ChartVertexIndex >::UnsignedIntIndex; };
	struct            GridNodeIndex : public UnsignedIntIndex<            GridNodeIndex >{ using UnsignedIntIndex<            GridNodeIndex >::UnsignedIntIndex; };
#if 0
	struct           AtlasEdgeIndex : public UnsignedIntIndex<           AtlasEdgeIndex >{ using UnsignedIntIndex<           AtlasEdgeIndex >::UnsignedIntIndex; };
	struct            GridEdgeIndex : public UnsignedIntIndex<            GridEdgeIndex >{ using UnsignedIntIndex<            GridEdgeIndex >::UnsignedIntIndex; };
//	struct     AtlasOrGridEdgeIndex : public UnsignedIntIndex<     AtlasOrGridEdgeIndex >{ using UnsignedIntIndex<     AtlasOrGridEdgeIndex >::UnsignedIntIndex; };
#endif
	struct       AtlasHalfEdgeIndex : public UnsignedIntIndex<       AtlasHalfEdgeIndex >{ using UnsignedIntIndex<       AtlasHalfEdgeIndex >::UnsignedIntIndex; };
	struct       ChartHalfEdgeIndex : public UnsignedIntIndex<       ChartHalfEdgeIndex >{ using UnsignedIntIndex<       ChartHalfEdgeIndex >::UnsignedIntIndex; };
//	struct AtlasBoundaryVertexIndex : public UnsignedIntIndex< AtlasBoundaryVertexIndex >{ using UnsignedIntIndex< AtlasBoundaryVertexIndex >::UnsignedIntIndex; };

#else // !DEBUG_INDEXING
	using               ChartIndex = unsigned int;
	using       AtlasTriangleIndex = unsigned int;
	using       ChartTriangleIndex = unsigned int;
	using         AtlasVertexIndex = unsigned int;
	using         ChartVertexIndex = unsigned int;
	using            GridNodeIndex = unsigned int;
#endif
#if 1
	using           AtlasEdgeIndex = unsigned int;
	using            GridEdgeIndex = unsigned int;
//	using     AtlasOrGridEdgeIndex = unsigned int;
#endif
	using AtlasBoundaryVertexIndex = unsigned int;
//#endif // DEBUG_INDEXING

	struct GridOrAtlasEdgeIndex
	{
		std::optional< GridEdgeIndex >  grid( void ) const { if( _type==_EdgeType::Grid  ) return  _grid ; else return std::nullopt; }
		std::optional< GridEdgeIndex > atlas( void ) const { if( _type==_EdgeType::Atlas ) return _atlas ; else return std::nullopt; }

		static GridOrAtlasEdgeIndex FromGridEdgeIndex( GridEdgeIndex idx )
		{
			GridOrAtlasEdgeIndex _idx;
			_idx._grid = idx;
			_idx._type = _EdgeType::Grid;
			return _idx;
		}
		static GridOrAtlasEdgeIndex FromAtlasEdgeIndex( AtlasEdgeIndex idx )
		{
			GridOrAtlasEdgeIndex _idx;
			_idx._atlas = idx;
			_idx._type = _EdgeType::Atlas;
			return _idx;
		}

		GridOrAtlasEdgeIndex( void ) : _type(_EdgeType::Undefined) {}

	protected:
		enum _EdgeType{ Grid , Atlas , Undefined };

		union
		{
			GridEdgeIndex _grid;
			AtlasEdgeIndex _atlas;
		};
		_EdgeType _type;
	};

#if 0
	struct AtlasOrChartVertexIndex
	{
		AtlasOrChartVertexIndex( void ) : atlas(-1) , chart(-1){}
		union
		{
			AtlasVertexIndex atlas;
			ChartVertexIndex chart;
		};
	};
#endif

	AtlasHalfEdgeIndex GetAtlasHalfEdgeIndex( AtlasTriangleIndex t , unsigned int k ){ return AtlasHalfEdgeIndex( static_cast< unsigned int >(t)*3+k ); }
	ChartHalfEdgeIndex GetChartHalfEdgeIndex( ChartTriangleIndex t , unsigned int k ){ return ChartHalfEdgeIndex( static_cast< unsigned int >(t)*3+k ); }
}