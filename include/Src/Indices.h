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
#ifdef USE_UNORDERED_MAP
#include "Misha/UnorderedMapsAndSets.h"
#else // !USE_UNORDERED_MAP
#include <map>
#endif // USE_UNORDERED_MAP
#ifdef DEBUG_INDEXING
#include "Misha/ExplicitIndex.h"
#endif // DEBUG_INDEXING


namespace MishaK
{
	namespace TSP
	{
#ifdef DEBUG_INDEXING
		template< typename T , typename Data >
		using ExplicitIndexVector = typename ExplicitIndex< unsigned int , T >::template ExplicitIndexVector< Data >;

		// Chart index
		struct                       ChartIndex : public ExplicitIndex< unsigned int ,                       ChartIndex >{ using ExplicitIndex< unsigned int ,                       ChartIndex >::ExplicitIndex; };

		// Mesh index types
		struct           AtlasMeshTriangleIndex : public ExplicitIndex< unsigned int ,           AtlasMeshTriangleIndex >{ using ExplicitIndex< unsigned int ,           AtlasMeshTriangleIndex >::ExplicitIndex; };
		struct           ChartMeshTriangleIndex : public ExplicitIndex< unsigned int ,           ChartMeshTriangleIndex >{ using ExplicitIndex< unsigned int ,           ChartMeshTriangleIndex >::ExplicitIndex; };
		struct             AtlasMeshVertexIndex : public ExplicitIndex< unsigned int ,             AtlasMeshVertexIndex >{ using ExplicitIndex< unsigned int ,             AtlasMeshVertexIndex >::ExplicitIndex; };
		struct             ChartMeshVertexIndex : public ExplicitIndex< unsigned int ,             ChartMeshVertexIndex >{ using ExplicitIndex< unsigned int ,             ChartMeshVertexIndex >::ExplicitIndex; };
		struct               AtlasMeshEdgeIndex : public ExplicitIndex< unsigned int ,               AtlasMeshEdgeIndex >{ using ExplicitIndex< unsigned int ,               AtlasMeshEdgeIndex >::ExplicitIndex; };
		struct           AtlasMeshHalfEdgeIndex : public ExplicitIndex< unsigned int ,           AtlasMeshHalfEdgeIndex >{ using ExplicitIndex< unsigned int ,           AtlasMeshHalfEdgeIndex >::ExplicitIndex; };
		struct           ChartMeshHalfEdgeIndex : public ExplicitIndex< unsigned int ,           ChartMeshHalfEdgeIndex >{ using ExplicitIndex< unsigned int ,           ChartMeshHalfEdgeIndex >::ExplicitIndex; };
		struct     AtlasMeshBoundaryVertexIndex : public ExplicitIndex< unsigned int ,     AtlasMeshBoundaryVertexIndex >{ using ExplicitIndex< unsigned int ,     AtlasMeshBoundaryVertexIndex >::ExplicitIndex; };

		// Grid index types
		struct             AtlasGridVertexIndex : public ExplicitIndex< unsigned int ,             AtlasGridVertexIndex >{ using ExplicitIndex< unsigned int ,             AtlasGridVertexIndex >::ExplicitIndex; };
		struct               AtlasGridEdgeIndex : public ExplicitIndex< unsigned int ,               AtlasGridEdgeIndex >{ using ExplicitIndex< unsigned int ,               AtlasGridEdgeIndex >::ExplicitIndex; };

		struct AtlasInteriorOrBoundaryNodeIndex : public ExplicitIndex< unsigned int , AtlasInteriorOrBoundaryNodeIndex >{ using ExplicitIndex< unsigned int , AtlasInteriorOrBoundaryNodeIndex >::ExplicitIndex; };

		// Cell index types
		struct           ChartInteriorCellIndex : public ExplicitIndex< unsigned int ,           ChartInteriorCellIndex >{ using ExplicitIndex< unsigned int ,           ChartInteriorCellIndex >::ExplicitIndex; };
		struct           ChartBoundaryCellIndex : public ExplicitIndex< unsigned int ,           ChartBoundaryCellIndex >{ using ExplicitIndex< unsigned int ,           ChartBoundaryCellIndex >::ExplicitIndex; };
		struct                   ChartCellIndex : public ExplicitIndex< unsigned int ,                   ChartCellIndex >{ using ExplicitIndex< unsigned int ,                   ChartCellIndex >::ExplicitIndex; };
		struct           AtlasInteriorCellIndex : public ExplicitIndex< unsigned int ,           AtlasInteriorCellIndex >{ using ExplicitIndex< unsigned int ,           AtlasInteriorCellIndex >::ExplicitIndex; };
		struct           AtlasBoundaryCellIndex : public ExplicitIndex< unsigned int ,           AtlasBoundaryCellIndex >{ using ExplicitIndex< unsigned int ,           AtlasBoundaryCellIndex >::ExplicitIndex; };
		struct                   AtlasCellIndex : public ExplicitIndex< unsigned int ,                   AtlasCellIndex >{ using ExplicitIndex< unsigned int ,                   AtlasCellIndex >::ExplicitIndex; };

		// Texel index types
		struct          AtlasInteriorTexelIndex : public ExplicitIndex< unsigned int ,          AtlasInteriorTexelIndex >{ using ExplicitIndex< unsigned int ,          AtlasInteriorTexelIndex >::ExplicitIndex; };
		struct           AtlasCoveredTexelIndex : public ExplicitIndex< unsigned int ,           AtlasCoveredTexelIndex >{ using ExplicitIndex< unsigned int ,           AtlasCoveredTexelIndex >::ExplicitIndex; };
		struct          AtlasBoundaryTexelIndex : public ExplicitIndex< unsigned int ,          AtlasBoundaryTexelIndex >{ using ExplicitIndex< unsigned int ,          AtlasBoundaryTexelIndex >::ExplicitIndex; };
		struct                  AtlasTexelIndex : public ExplicitIndex< unsigned int ,                  AtlasTexelIndex >{ using ExplicitIndex< unsigned int ,                  AtlasTexelIndex >::ExplicitIndex; };

		// Boundary index types
		struct            BoundaryMidPointIndex : public ExplicitIndex< unsigned int ,            BoundaryMidPointIndex >{ using ExplicitIndex< unsigned int ,            BoundaryMidPointIndex >::ExplicitIndex; };

		// Auxilary node index type
		struct       ChartBoundaryTriangleIndex : public ExplicitIndex< unsigned int ,       ChartBoundaryTriangleIndex >{ using ExplicitIndex< unsigned int ,       ChartBoundaryTriangleIndex >::ExplicitIndex; };
		struct  AtlasRefinedBoundaryVertexIndex : public ExplicitIndex< unsigned int ,  AtlasRefinedBoundaryVertexIndex >{ using ExplicitIndex< unsigned int ,  AtlasRefinedBoundaryVertexIndex >::ExplicitIndex; };
		struct    AtlasRefinedBoundaryEdgeIndex : public ExplicitIndex< unsigned int ,    AtlasRefinedBoundaryEdgeIndex >{ using ExplicitIndex< unsigned int ,    AtlasRefinedBoundaryEdgeIndex >::ExplicitIndex; };

#else // !DEBUG_INDEXING
		template< typename T , typename Data >
		using ExplicitIndexVector = std::vector< Data >;

		// Chart index
		using                       ChartIndex = unsigned int;

		using           AtlasMeshTriangleIndex = unsigned int;
		using           ChartMeshTriangleIndex = unsigned int;
		using             AtlasMeshVertexIndex = unsigned int;
		using             ChartMeshVertexIndex = unsigned int;
		using           AtlasMeshHalfEdgeIndex = unsigned int;
		using           ChartMeshHalfEdgeIndex = unsigned int;
		using               AtlasMeshEdgeIndex = unsigned int;
		using     AtlasMeshBoundaryVertexIndex = unsigned int;
		using       ChartBoundaryTriangleIndex = unsigned int;

		using             AtlasGridVertexIndex = unsigned int;
		using               AtlasGridEdgeIndex = unsigned int;

		using AtlasInteriorOrBoundaryNodeIndex = unsigned int;

		// Cell types
		using           ChartInteriorCellIndex = unsigned int;
		using           ChartBoundaryCellIndex = unsigned int;
		using                   ChartCellIndex = unsigned int;
		using           AtlasInteriorCellIndex = unsigned int;
		using           AtlasBoundaryCellIndex = unsigned int;
		using                   AtlasCellIndex = unsigned int;

		// Texel types
		using          AtlasInteriorTexelIndex = unsigned int;
		using           AtlasCoveredTexelIndex = unsigned int;
		using          AtlasBoundaryTexelIndex = unsigned int;
		using                  AtlasTexelIndex = unsigned int;

		// Boundary types
		using            BoundaryMidPointIndex = unsigned int;

		// Auxiliary node index type
		using       ChartBoundaryTriangleIndex = unsigned int;
		using  AtlasRefinedBoundaryVertexIndex = unsigned int;
		using    AtlasRefinedBoundaryEdgeIndex = unsigned int;

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

		inline AtlasMeshHalfEdgeIndex GetAtlasMeshHalfEdgeIndex( AtlasMeshTriangleIndex t , unsigned int k ){ return AtlasMeshHalfEdgeIndex( static_cast< unsigned int >(t)*3+k ); }
		inline ChartMeshHalfEdgeIndex GetChartMeshHalfEdgeIndex( ChartMeshTriangleIndex t , unsigned int k ){ return ChartMeshHalfEdgeIndex( static_cast< unsigned int >(t)*3+k ); }
		inline std::pair< AtlasMeshTriangleIndex , unsigned int > FactorAtlasMeshHalfEdgeIndex( AtlasMeshHalfEdgeIndex he )
		{
			unsigned int _he = static_cast< unsigned int >( he );
			return std::make_pair( AtlasMeshTriangleIndex( _he/3 ) , _he%3 );
		}
		inline std::pair< ChartMeshTriangleIndex , unsigned int > FactorChartMeshHalfEdgeIndex( ChartMeshHalfEdgeIndex he )
		{
			unsigned int _he = static_cast< unsigned int >( he );
			return std::make_pair( ChartMeshTriangleIndex( _he/3 ) , _he%3 );
		}
	}

#ifdef USE_UNORDERED_MAP
	template< typename Key , typename Value >
	using Map = UnorderedMap< Key , Value >;
#else // !USE_UNORDERED_MAP
	template< typename Key , typename Value >
	using Map = std::map< Key , Value >;
#endif // USE_UNORDERED_MAP
}