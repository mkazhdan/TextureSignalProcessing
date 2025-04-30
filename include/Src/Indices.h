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
	// A wrapper class for unsigned int that allows distinguishing between different types of indices
	template< typename T >
	struct UnsignedIntIndex
	{
		explicit UnsignedIntIndex( unsigned int idx=-1 ) : _idx(idx){}
		explicit operator unsigned int () const { return _idx; }

		UnsignedIntIndex & operator += ( unsigned int off ){ _idx += off ; return *this; }
		T operator + ( unsigned int off ) const { return T( _idx+off ); }

		bool operator == ( const UnsignedIntIndex &I ) const { return _idx==I._idx; }
		bool operator != ( const UnsignedIntIndex &I ) const { return _idx!=I._idx; }
		bool operator <  ( const UnsignedIntIndex &I ) const { return _idx< I._idx; }

		friend std::ostream &operator << ( std::ostream &os , const UnsignedIntIndex &I ){ return os << I._idx; }
	protected:
		unsigned int _idx;
	};

	struct                       ChartIndex : public UnsignedIntIndex<                       ChartIndex >{ using UnsignedIntIndex<                       ChartIndex >::UnsignedIntIndex; };
	struct           AtlasMeshTriangleIndex : public UnsignedIntIndex<           AtlasMeshTriangleIndex >{ using UnsignedIntIndex<           AtlasMeshTriangleIndex >::UnsignedIntIndex; };
	struct           ChartMeshTriangleIndex : public UnsignedIntIndex<           ChartMeshTriangleIndex >{ using UnsignedIntIndex<           ChartMeshTriangleIndex >::UnsignedIntIndex; };
	struct             AtlasMeshVertexIndex : public UnsignedIntIndex<             AtlasMeshVertexIndex >{ using UnsignedIntIndex<             AtlasMeshVertexIndex >::UnsignedIntIndex; };
	struct             ChartMeshVertexIndex : public UnsignedIntIndex<             ChartMeshVertexIndex >{ using UnsignedIntIndex<             ChartMeshVertexIndex >::UnsignedIntIndex; };
	struct             AtlasGridVertexIndex : public UnsignedIntIndex<             AtlasGridVertexIndex >{ using UnsignedIntIndex<             AtlasGridVertexIndex >::UnsignedIntIndex; };
	struct               AtlasMeshEdgeIndex : public UnsignedIntIndex<               AtlasMeshEdgeIndex >{ using UnsignedIntIndex<               AtlasMeshEdgeIndex >::UnsignedIntIndex; };
	struct               AtlasGridEdgeIndex : public UnsignedIntIndex<               AtlasGridEdgeIndex >{ using UnsignedIntIndex<               AtlasGridEdgeIndex >::UnsignedIntIndex; };
	struct           AtlasMeshHalfEdgeIndex : public UnsignedIntIndex<           AtlasMeshHalfEdgeIndex >{ using UnsignedIntIndex<           AtlasMeshHalfEdgeIndex >::UnsignedIntIndex; };
	struct           ChartMeshHalfEdgeIndex : public UnsignedIntIndex<           ChartMeshHalfEdgeIndex >{ using UnsignedIntIndex<           ChartMeshHalfEdgeIndex >::UnsignedIntIndex; };
	struct AtlasInteriorOrBoundaryNodeIndex : public UnsignedIntIndex< AtlasInteriorOrBoundaryNodeIndex >{ using UnsignedIntIndex< AtlasInteriorOrBoundaryNodeIndex >::UnsignedIntIndex; };

#else // !DEBUG_INDEXING
	using                       ChartIndex = unsigned int;
	using           AtlasMeshTriangleIndex = unsigned int;
	using           ChartMeshTriangleIndex = unsigned int;
	using             AtlasMeshVertexIndex = unsigned int;
	using             ChartMeshVertexIndex = unsigned int;
	using             AtlasGridVertexIndex = unsigned int;
	using           AtlasMeshHalfEdgeIndex = unsigned int;
	using           ChartMeshHalfEdgeIndex = unsigned int;
	using               AtlasMeshEdgeIndex = unsigned int;
	using               AtlasGridEdgeIndex = unsigned int;
	using AtlasInteriorOrBoundaryNodeIndex = unsigned int;
#endif // DEBUG_INDEXING

	struct AtlasGridOrMeshEdgeIndex
	{
		std::optional< AtlasGridEdgeIndex > grid( void ) const { if( _type==_EdgeType::Grid ) return _grid ; else return std::nullopt; }
		std::optional< AtlasMeshEdgeIndex > mesh( void ) const { if( _type==_EdgeType::Mesh ) return _mesh ; else return std::nullopt; }

		static AtlasGridOrMeshEdgeIndex FromGrid( AtlasGridEdgeIndex idx )
		{
			AtlasGridOrMeshEdgeIndex _idx;
			_idx._grid = idx;
			_idx._type = _EdgeType::Grid;
			return _idx;
		}
		static AtlasGridOrMeshEdgeIndex FromMesh( AtlasMeshEdgeIndex idx )
		{
			AtlasGridOrMeshEdgeIndex _idx;
			_idx._mesh = idx;
			_idx._type = _EdgeType::Mesh;
			return _idx;
		}

		AtlasGridOrMeshEdgeIndex( void ) : _type(_EdgeType::Undefined) {}

	protected:
		enum _EdgeType{ Grid , Mesh , Undefined };

		union
		{
			AtlasGridEdgeIndex _grid;
			AtlasMeshEdgeIndex _mesh;
		};
		_EdgeType _type;
	};

	AtlasMeshHalfEdgeIndex GetAtlasMeshHalfEdgeIndex( AtlasMeshTriangleIndex t , unsigned int k ){ return AtlasMeshHalfEdgeIndex( static_cast< unsigned int >(t)*3+k ); }
	ChartMeshHalfEdgeIndex GetChartMeshHalfEdgeIndex( ChartMeshTriangleIndex t , unsigned int k ){ return ChartMeshHalfEdgeIndex( static_cast< unsigned int >(t)*3+k ); }
	std::pair< AtlasMeshTriangleIndex , unsigned int > FactorAtlasMeshHalfEdgeIndex( AtlasMeshHalfEdgeIndex he )
	{
		unsigned int _he = static_cast< unsigned int >( he );
		return std::make_pair( AtlasMeshTriangleIndex( _he/3 ) , _he%3 );
	}
	std::pair< ChartMeshTriangleIndex , unsigned int > FactorChartMeshHalfEdgeIndex( ChartMeshHalfEdgeIndex he )
	{
		unsigned int _he = static_cast< unsigned int >( he );
		return std::make_pair( ChartMeshTriangleIndex( _he/3 ) , _he%3 );
	}
}