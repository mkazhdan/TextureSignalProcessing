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
		template< typename I=unsigned int >
		explicit UnsignedIntIndex( I idx=static_cast<I>(-1) )
			: _idx(static_cast< unsigned int >(idx) )
		{ static_assert( std::is_integral_v< I > , "[ERROR] Expected integral type" ); }

		explicit operator unsigned int () const { return _idx; }

		UnsignedIntIndex & operator += ( unsigned int off ){ _idx += off ; return *this; }
		UnsignedIntIndex & operator -= ( unsigned int off ){ _idx -= off ; return *this; }
		T operator + ( unsigned int off ) const { return T( _idx+off ); }
		T operator - ( unsigned int off ) const { return T( _idx-off ); }

		unsigned int operator - ( const UnsignedIntIndex &t ) const { return _idx - t._idx; }
		friend T & operator ++ ( T & t ){ t._idx++ ; return t; }
		friend T operator ++ ( T &t , int ){ return T( t._idx++ ); }

		bool operator == ( const UnsignedIntIndex &I ) const { return _idx==I._idx; }
		bool operator != ( const UnsignedIntIndex &I ) const { return _idx!=I._idx; }
		bool operator <  ( const UnsignedIntIndex &I ) const { return _idx< I._idx; }

		friend std::ostream &operator << ( std::ostream &os , const UnsignedIntIndex &I ){ return os << I._idx; }

		template< typename Data >
		struct IndexVector : protected std::vector< Data >
		{
			IndexVector( void ) : std::vector< Data >(){}
			IndexVector( size_t sz ) : std::vector< Data >(sz){}
			IndexVector( size_t sz , const Data & data ) : std::vector< Data >( sz , data ){}
			size_t size( void ) const { return std::vector< Data >::size(); }
			void resize( size_t sz ){ return std::vector< Data >::resize(sz); }
			void resize( size_t sz , const Data & data ){ return std::vector< Data >::resize( sz , data ); }
			void push_back( const Data & data ){ return std::vector< Data >::push_back(data); }
			Data & operator[]( const T & t ){ return std::vector< Data >::operator[]( static_cast< unsigned int >(t) ); }
			const Data & operator[]( const T &t ) const { return std::vector< Data >::operator[]( static_cast< unsigned int >(t) ); }

			typename std::vector< Data >::iterator begin( void ){ return std::vector< Data >::begin(); }
			typename std::vector< Data >::iterator end( void ){ return std::vector< Data >::end(); }
			typename std::vector< Data >::const_iterator begin( void ) const { return std::vector< Data >::begin(); }
			typename std::vector< Data >::const_iterator end( void ) const { return std::vector< Data >::end(); }

			template< class ... Args >
			void emplace_back( Args && ... args ){ std::vector< Data >::emplace_back( std::forward< Args >(args) ... ); }

			template< class InputIt >
			typename std::vector< Data >::iterator insert( typename std::vector< Data >::const_iterator pos , InputIt first , InputIt last ){ return std::vector< Data >::insert( pos , first , last ); }

			explicit operator const std::vector< Data > & () const { return *this; }
			explicit operator       std::vector< Data > & ()       { return *this; }
		};
	protected:
		unsigned int _idx;
	};

	template< typename T , typename Data >
	using IndexVector = typename UnsignedIntIndex< T >::template IndexVector< Data >;

	template< typename T , typename Data >
	Data * operator + ( Data * data , const T & t ){ return data + static_cast< unsigned int >(t); }

	template< typename T , typename Data >
	const Data * operator + ( const Data * data , const T & t ){ return data + static_cast< unsigned int >(t); }

	// Chart index
	struct                       ChartIndex : public UnsignedIntIndex<                       ChartIndex >{ using UnsignedIntIndex<                       ChartIndex >::UnsignedIntIndex; };

	// Mesh index types
	struct           AtlasMeshTriangleIndex : public UnsignedIntIndex<           AtlasMeshTriangleIndex >{ using UnsignedIntIndex<           AtlasMeshTriangleIndex >::UnsignedIntIndex; };
	struct           ChartMeshTriangleIndex : public UnsignedIntIndex<           ChartMeshTriangleIndex >{ using UnsignedIntIndex<           ChartMeshTriangleIndex >::UnsignedIntIndex; };
	struct             AtlasMeshVertexIndex : public UnsignedIntIndex<             AtlasMeshVertexIndex >{ using UnsignedIntIndex<             AtlasMeshVertexIndex >::UnsignedIntIndex; };
	struct             ChartMeshVertexIndex : public UnsignedIntIndex<             ChartMeshVertexIndex >{ using UnsignedIntIndex<             ChartMeshVertexIndex >::UnsignedIntIndex; };
	struct               AtlasMeshEdgeIndex : public UnsignedIntIndex<               AtlasMeshEdgeIndex >{ using UnsignedIntIndex<               AtlasMeshEdgeIndex >::UnsignedIntIndex; };
	struct           AtlasMeshHalfEdgeIndex : public UnsignedIntIndex<           AtlasMeshHalfEdgeIndex >{ using UnsignedIntIndex<           AtlasMeshHalfEdgeIndex >::UnsignedIntIndex; };
	struct           ChartMeshHalfEdgeIndex : public UnsignedIntIndex<           ChartMeshHalfEdgeIndex >{ using UnsignedIntIndex<           ChartMeshHalfEdgeIndex >::UnsignedIntIndex; };
	struct     AtlasMeshBoundaryVertexIndex : public UnsignedIntIndex<     AtlasMeshBoundaryVertexIndex >{ using UnsignedIntIndex<     AtlasMeshBoundaryVertexIndex >::UnsignedIntIndex; };

	// Grid index types
	struct             AtlasGridVertexIndex : public UnsignedIntIndex<             AtlasGridVertexIndex >{ using UnsignedIntIndex<             AtlasGridVertexIndex >::UnsignedIntIndex; };
	struct               AtlasGridEdgeIndex : public UnsignedIntIndex<               AtlasGridEdgeIndex >{ using UnsignedIntIndex<               AtlasGridEdgeIndex >::UnsignedIntIndex; };

	struct AtlasInteriorOrBoundaryNodeIndex : public UnsignedIntIndex< AtlasInteriorOrBoundaryNodeIndex >{ using UnsignedIntIndex< AtlasInteriorOrBoundaryNodeIndex >::UnsignedIntIndex; };

	// Cell index types
	struct           ChartInteriorCellIndex : public UnsignedIntIndex<           ChartInteriorCellIndex >{ using UnsignedIntIndex<           ChartInteriorCellIndex >::UnsignedIntIndex; };
	struct           ChartBoundaryCellIndex : public UnsignedIntIndex<           ChartBoundaryCellIndex >{ using UnsignedIntIndex<           ChartBoundaryCellIndex >::UnsignedIntIndex; };
	struct           ChartCombinedCellIndex : public UnsignedIntIndex<           ChartCombinedCellIndex >{ using UnsignedIntIndex<           ChartCombinedCellIndex >::UnsignedIntIndex; };
	struct           AtlasInteriorCellIndex : public UnsignedIntIndex<           AtlasInteriorCellIndex >{ using UnsignedIntIndex<           AtlasInteriorCellIndex >::UnsignedIntIndex; };
	struct           AtlasBoundaryCellIndex : public UnsignedIntIndex<           AtlasBoundaryCellIndex >{ using UnsignedIntIndex<           AtlasBoundaryCellIndex >::UnsignedIntIndex; };
	struct           AtlasCombinedCellIndex : public UnsignedIntIndex<           AtlasCombinedCellIndex >{ using UnsignedIntIndex<           AtlasCombinedCellIndex >::UnsignedIntIndex; };

	// Texel index types
	struct          AtlasInteriorTexelIndex : public UnsignedIntIndex<          AtlasInteriorTexelIndex >{ using UnsignedIntIndex<          AtlasInteriorTexelIndex >::UnsignedIntIndex; };
	struct           AtlasCoveredTexelIndex : public UnsignedIntIndex<           AtlasCoveredTexelIndex >{ using UnsignedIntIndex<           AtlasCoveredTexelIndex >::UnsignedIntIndex; };
	struct          AtlasBoundaryTexelIndex : public UnsignedIntIndex<          AtlasBoundaryTexelIndex >{ using UnsignedIntIndex<          AtlasBoundaryTexelIndex >::UnsignedIntIndex; };
	struct          AtlasCombinedTexelIndex : public UnsignedIntIndex<          AtlasCombinedTexelIndex >{ using UnsignedIntIndex<          AtlasCombinedTexelIndex >::UnsignedIntIndex; };

	// Boundary index types
	struct            BoundaryMidPointIndex : public UnsignedIntIndex<            BoundaryMidPointIndex >{ using UnsignedIntIndex<            BoundaryMidPointIndex >::UnsignedIntIndex; };

	// Auxilary node index type
	struct  AtlasRefinedBoundaryVertexIndex : public UnsignedIntIndex<  AtlasRefinedBoundaryVertexIndex >{ using UnsignedIntIndex<  AtlasRefinedBoundaryVertexIndex >::UnsignedIntIndex; };

#else // !DEBUG_INDEXING
	template< typename T , typename Data >
	using IndexVector = std::vector< Data >;

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

	// Cell types
	using           ChartInteriorCellIndex = unsigned int;
	using           ChartBoundaryCellIndex = unsigned int;
	using           ChartCombinedCellIndex = unsigned int;
	using           AtlasInteriorCellIndex = unsigned int;
	using           AtlasBoundaryCellIndex = unsigned int;
	using           AtlasCombinedCellIndex = unsigned int;

	// Texel types
	using          AtlasInteriorTexelIndex = unsigned int;
	using           AtlasCoveredTexelIndex = unsigned int;
	using          AtlasBoundaryTexelIndex = unsigned int;
	using          AtlasCombinedTexelIndex = unsigned int;

	// Boundary types
	using             BoundaryMiPointIndex = unsigned int;

	// Auxiliary node index type
	using AtlasRefinedMeshBoundaryVertexIndex = unsigned int;

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