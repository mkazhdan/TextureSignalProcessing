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
#ifndef HIERARCHICAL_SYSTEM_INCLUDED
#define HIERARCHICAL_SYSTEM_INCLUDED

#include <algorithm>
#include <unordered_set>
#include <Eigen/Sparse>
#include <Misha/LinearSolvers.h>
#include <Misha/RegularGrid.h>
#include <Misha/Miscellany.h>
#include <Misha/Exceptions.h>
#include <Misha/Texels.h>
#include <Misha/Rasterizer2D.h>
#include <Misha/Geometry.h>
#include <Misha/MultiThreading.h>
#include <Misha/Atomic.h>
#include <Misha/Atomic.Geometry.h>
#include <Misha/SimplexBasis.h>
#include "IndexedPolygon.h"
//#include "ImageIO.h"

namespace MishaK
{
	namespace TSP
	{
		template< typename Real >
		struct BilinearElementScalarSample
		{
			struct SampleData
			{
				Real dualValues[4];				// The integrated values of the four incident bilinear basis functions, dualized
				SampleData( void ){ _weights[0] = _weights[1] = _weights[2] = _weights[3] = 0; dualValues[0] = dualValues[1] = dualValues[2] = dualValues[3] = 0; }
				void init( Point2D< Real > p )
				{
					_weights[0] = ( 1 - p[0] ) * ( 1 - p[1] );
					_weights[1] = (     p[0] ) * ( 1 - p[1] );
					_weights[2] = (     p[0] ) * (     p[1] );
					_weights[3] = ( 1 - p[0] ) * (     p[1] );
				}
				template< typename T >
				T operator()( const T values[4] ) const
				{
					return values[0] * _weights[0] + values[1] * _weights[1] + values[2] * _weights[2] + values[3] * _weights[3];
				}
			protected:
				Real _weights[4];
			};
			SquareMatrix< Real , 2 > tensor;	// The inverse metric tensor defined by the intersecting triangle
			unsigned int cellOffset;
			bool operator < ( const BilinearElementScalarSample& sample ) const { return cellOffset<sample.cellOffset; }

			BilinearElementScalarSample( void ) : _sampleNum(0) , _samples(nullptr) {}
			BilinearElementScalarSample( unsigned int sz ) : _sampleNum(0) , _samples(nullptr) { resize(sz); }
			BilinearElementScalarSample( const BilinearElementScalarSample& bilinearElementScalarSample ) : _sampleNum(0) , _samples(nullptr)
			{
				resize( bilinearElementScalarSample._sampleNum );
				memcpy( _samples , bilinearElementScalarSample._samples , sizeof( SampleData ) * _sampleNum );
				tensor = bilinearElementScalarSample.tensor;
				cellOffset = bilinearElementScalarSample.cellOffset;
			}

			BilinearElementScalarSample& operator = ( const BilinearElementScalarSample& bilinearElementScalarSample )
			{
				resize( bilinearElementScalarSample._sampleNum );
				memcpy( _samples , bilinearElementScalarSample._samples , sizeof( SampleData ) * _sampleNum );
				tensor = bilinearElementScalarSample.tensor;
				cellOffset = bilinearElementScalarSample.cellOffset;
				return *this;
			}

			~BilinearElementScalarSample( void ){ resize(0); }

			void resize( unsigned int sz )
			{
				if( _samples ){ delete[] _samples ; _samples = nullptr; }
				_sampleNum = 0;
				if( sz ){ _samples = new SampleData[sz] ; _sampleNum = sz; }
			}
			unsigned int size( void ) const { return _sampleNum; }
			SampleData& operator[]( unsigned int idx ){ return _samples[idx]; }
			const SampleData& operator[]( unsigned int idx ) const { return _samples[idx]; }
		protected:
			unsigned int _sampleNum;
			SampleData* _samples;
		};

		template< typename Real >
		struct QuadraticElementScalarSample
		{
			struct SampleData
			{
				Real dualValues[6];				// The integrated values of the six incident quadratic basis functions, dualized
				SampleData( void )
				{
					_weights[0] = _weights[1] = _weights[2] = _weights[3] = _weights[4] = _weights[5] = 0;
					dualValues[0] = dualValues[1] = dualValues[2] = dualValues[3] = dualValues[4] = dualValues[5] = 0;
				}
				void init( Point2D< Real > p )
				{
					for( unsigned int i=0 ; i<6 ; i++ ) _weights[i] = static_cast< Real >( SimplexElements< 2 , 2 >::Element( _idx[i] )( Point2D< double >(p) ) );
				}
				template< typename T >
				T operator()( const T values[6] ) const
				{
					return values[0]*_weights[0] + values[1]*_weights[1] + values[2]*_weights[2] + values[3]*_weights[3] + values[4]*_weights[4] + values[5]*_weights[5];
				}
			protected:
				Real _weights[6];
			};
			SquareMatrix< Real , 2 > tensor;	// The inverse metric tensor defined by the intersecting triangle
			AtlasInteriorOrBoundaryNodeIndex fineNodes[6];

			QuadraticElementScalarSample( void ) : _sampleNum(0) , _samples(nullptr) {}
			QuadraticElementScalarSample( unsigned int sz ) : _sampleNum(0) , _samples(nullptr) { resize(sz); }
			QuadraticElementScalarSample( const QuadraticElementScalarSample& quadraticElementScalarSample ) : _sampleNum(0) , _samples(nullptr)
			{
				resize( quadraticElementScalarSample._sampleNum );
				memcpy( _samples , quadraticElementScalarSample._samples , sizeof( SampleData ) * _sampleNum );
				tensor = quadraticElementScalarSample.tensor;
				memcpy( fineNodes , quadraticElementScalarSample.fineNodes , sizeof(int)*6 );
			}
			QuadraticElementScalarSample& operator = ( const QuadraticElementScalarSample& quadraticElementScalarSample )
			{
				resize( quadraticElementScalarSample._sampleNum );
				memcpy( _samples , quadraticElementScalarSample._samples , sizeof( SampleData ) * _sampleNum );
				tensor = quadraticElementScalarSample.tensor;
				memcpy( fineNodes , quadraticElementScalarSample.fineNodes , sizeof(int)*6 );
				return *this;
			}
			~QuadraticElementScalarSample( void ){ resize(0); }
			void resize( unsigned int sz )
			{
				if( _samples ){ delete[] _samples ; _samples = nullptr; }
				_sampleNum = 0;
				if( sz ){ _samples = new SampleData[sz] ; _sampleNum = sz; }
			}
			unsigned int size( void ) const { return _sampleNum; }
			SampleData& operator[]( unsigned int idx )       { return _samples[idx]; }
			const SampleData& operator[]( unsigned int idx ) const { return _samples[idx]; }
		protected:
			unsigned int _sampleNum;
			SampleData* _samples;
			static const unsigned int _idx[];
		};

		template< typename Real >
		const unsigned int QuadraticElementScalarSample< Real >::_idx[] = 
		{
			SimplexElements< 2 , 2 >::NodeIndex( 0 , 0 ) ,
			SimplexElements< 2 , 2 >::NodeIndex( 1 , 1 ) ,
			SimplexElements< 2 , 2 >::NodeIndex( 2 , 2 ) ,
			SimplexElements< 2 , 2 >::NodeIndex( 1 , 2 ) ,
			SimplexElements< 2 , 2 >::NodeIndex( 2 , 0 ) ,
			SimplexElements< 2 , 2 >::NodeIndex( 0 , 1 ),
		};

		template< typename Real >
		struct BilinearElementGradientSample
		{
			struct SampleData
			{
				Point2D< Real > dualGradients[4];	// The integrated gradients of the four incident bilinear basis functions, dualized
				void init( Point2D< Real > p )
				{
					_weights[0] = Point2D< Real >(   p[1] - 1 ,    p[0] - 1 );
					_weights[1] = Point2D< Real >( - p[1] + 1 ,  - p[0]     );
					_weights[2] = Point2D< Real >(   p[1]     ,    p[0]     );
					_weights[3] = Point2D< Real >( - p[1]     ,  - p[0] + 1 );
				}
				template< typename T >
				Point2D< T > operator()( const T values[4] ) const
				{
					return Point2D< T >
						(
							values[0] * _weights[0][0] + values[1] * _weights[1][0] + values[2] * _weights[2][0] + values[3] * _weights[3][0] ,
							values[0] * _weights[0][1] + values[1] * _weights[1][1] + values[2] * _weights[2][1] + values[3] * _weights[3][1]
						);
				}
			protected:
				Point2D< Real > _weights[4];
			};
			SquareMatrix< Real , 2 > tensor;		// The inverse metric tensor defined by the intersecting triangle
			unsigned int cellOffset;
			bool operator < ( const BilinearElementGradientSample& sample ) const { return cellOffset < sample.cellOffset; }

			BilinearElementGradientSample( void ) : _sampleNum(0) , _samples(nullptr) {}
			BilinearElementGradientSample( unsigned int sz ) : _sampleNum(0) , _samples(nullptr) { resize(sz); }
			BilinearElementGradientSample( const BilinearElementGradientSample& bilinearElementGradientSample ) : _sampleNum(0) , _samples(nullptr)
			{
				resize( bilinearElementGradientSample._sampleNum );
				memcpy( _samples , bilinearElementGradientSample._samples , sizeof( SampleData ) * _sampleNum );
				tensor = bilinearElementGradientSample.tensor;
				cellOffset = bilinearElementGradientSample.cellOffset;
			}
			BilinearElementGradientSample& operator = ( const BilinearElementGradientSample& bilinearElementGradientSample )
			{
				resize( bilinearElementGradientSample._sampleNum );
				memcpy( _samples , bilinearElementGradientSample._samples , sizeof( SampleData ) * _sampleNum );
				tensor = bilinearElementGradientSample.tensor;
				cellOffset = bilinearElementGradientSample.cellOffset;
				return *this;
			}
			~BilinearElementGradientSample( void ){ resize(0); }
			void resize( unsigned int sz )
			{
				if( _samples ){ delete[] _samples ; _samples = nullptr; }
				_sampleNum = 0;
				if( sz ){ _samples = new SampleData[sz] ; _sampleNum = sz; }
			}
			unsigned int size( void ) const { return _sampleNum; }
			SampleData& operator[]( unsigned int idx )       { return _samples[idx]; }
			const SampleData& operator[]( unsigned int idx ) const { return _samples[idx]; }
		protected:
			unsigned int _sampleNum;
			SampleData* _samples;
		};

		template< typename Real >
		struct QuadraticElementGradientSample
		{
			struct SampleData
			{
				Point2D< Real > dualGradients[6];	// The integrated gradients of the six incident quadratic basis functions, dualized
				void init( Point2D< Real > p )
				{
					for( unsigned int i=0 ; i<6 ; i++ )
					{
						const Point< Polynomial::Polynomial< 2 , 1 , double > , 2 > & d = SimplexElements< 2 , 2 >::Differential( _idx[i] );
						_weights[i] = Point2D< Real >( d[0]( Point2D< double >(p) ) , d[1]( Point2D< double >(p) ) );
					}
				}
				template< typename T >
				Point2D< T > operator()( const T values[6] ) const
				{
					return Point2D< T >
						(
							values[0] * _weights[0][0] + values[1] * _weights[1][0] + values[2] * _weights[2][0] + values[3] * _weights[3][0] + values[4] * _weights[4][0] + values[5] * _weights[5][0] ,
							values[0] * _weights[0][1] + values[1] * _weights[1][1] + values[2] * _weights[2][1] + values[3] * _weights[3][1] + values[4] * _weights[4][1] + values[5] * _weights[5][1]
						);
				}
			protected:
				Point2D< Real > _weights[6];
			};
			SquareMatrix< Real , 2 > tensor;		// The inverse metric tensor defined by the intersecting triangle
			AtlasInteriorOrBoundaryNodeIndex fineNodes[6];

			QuadraticElementGradientSample( void ) : _sampleNum(0) , _samples(nullptr) {}
			QuadraticElementGradientSample( unsigned int sz ) : _sampleNum(0) , _samples(nullptr) { resize(sz); }
			QuadraticElementGradientSample( const QuadraticElementGradientSample& quadraticElementGradientSample ) : _sampleNum(0) , _samples(nullptr)
			{
				resize( quadraticElementGradientSample._sampleNum );
				memcpy( _samples , quadraticElementGradientSample._samples , sizeof( SampleData ) * _sampleNum );
				tensor = quadraticElementGradientSample.tensor;
				memcpy( fineNodes , quadraticElementGradientSample.fineNodes , sizeof(int)*6 );
			}
			QuadraticElementGradientSample& operator = ( const QuadraticElementGradientSample& quadraticElementGradientSample )
			{
				resize( quadraticElementGradientSample._sampleNum );
				memcpy( _samples , quadraticElementGradientSample._samples , sizeof( SampleData ) * _sampleNum );
				tensor = quadraticElementGradientSample.tensor;
				memcpy( fineNodes , quadraticElementGradientSample.fineNodes , sizeof(int)*6 );
				return *this;
			}
			~QuadraticElementGradientSample( void ){ resize(0); }
			void resize( unsigned int sz )
			{
				if( _samples ){ delete[] _samples ; _samples = nullptr; }
				_sampleNum = 0;
				if( sz ){ _samples = new SampleData[sz] ; _sampleNum = sz; }
			}
			unsigned int size( void ) const { return _sampleNum; }
			SampleData& operator[]( unsigned int idx )       { return _samples[idx]; }
			const SampleData& operator[]( unsigned int idx ) const { return _samples[idx]; }
		protected:
			unsigned int _sampleNum;
			SampleData* _samples;
			static const unsigned int _idx[];
		};

		template< typename Real >
		const unsigned int QuadraticElementGradientSample< Real >::_idx[] = 
		{
			SimplexElements< 2 , 2 >::NodeIndex( 0 , 0 ) ,
			SimplexElements< 2 , 2 >::NodeIndex( 1 , 1 ) ,
			SimplexElements< 2 , 2 >::NodeIndex( 2 , 2 ) ,
			SimplexElements< 2 , 2 >::NodeIndex( 1 , 2 ) ,
			SimplexElements< 2 , 2 >::NodeIndex( 2 , 0 ) ,
			SimplexElements< 2 , 2 >::NodeIndex( 0 , 1 ),
		};

		template< typename _Real >
		struct ScalarElementSamples
		{
			typedef _Real Real;
			typedef BilinearElementScalarSample< Real > Bilinear;
			typedef QuadraticElementScalarSample< Real > Quadratic;
			std::vector< std::vector< Bilinear > > bilinear;
			std::vector< Quadratic > quadratic;
			void resize( size_t sz ){ bilinear.resize( sz ); }
			void sort( void ){ for( int i=0 ; i<bilinear.size() ; i++ ) std::sort( bilinear[i].begin() , bilinear[i].end() ); }
		};

		template< typename _Real >
		struct GradientElementSamples
		{
			typedef _Real Real;
			typedef BilinearElementGradientSample< Real > Bilinear;
			typedef QuadraticElementGradientSample< Real > Quadratic;
			std::vector< std::vector< Bilinear > > bilinear;
			std::vector< Quadratic > quadratic;
			void resize( size_t sz ){ bilinear.resize( sz ); }
			void sort( void ){ for( int i=0 ; i<bilinear.size() ; i++ ) std::sort( bilinear[i].begin() , bilinear[i].end() ); }
		};

		struct RasterLine
		{
			AtlasTexelIndex lineStartIndex;
			AtlasTexelIndex lineEndIndex;
			AtlasTexelIndex prevLineIndex;
			AtlasTexelIndex nextLineIndex;
			AtlasInteriorTexelIndex coeffStartIndex;
		};

		struct DeepLine
		{
			AtlasInteriorTexelIndex coarseLineStartIndex;
			AtlasInteriorTexelIndex coarseLineEndIndex;
			AtlasInteriorTexelIndex finePrevLineIndex;
			AtlasInteriorTexelIndex fineCurrentLineIndex;
			AtlasInteriorTexelIndex fineNextLineIndex;
		};

		struct SegmentedRasterLine
		{
			std::vector< RasterLine > segments;
		};

		struct MultigridBlockInfo
		{
			MultigridBlockInfo( unsigned int p_blockWidth=128 , unsigned int p_blockHeight=16 , unsigned int p_paddingWidth=2 , unsigned int p_paddingHeight=2 )
			{
				blockWidth = p_blockWidth;
				blockHeight = p_blockHeight;
				paddingWidth = p_paddingWidth;
				paddingHeight = p_paddingHeight;
			}
			unsigned int blockWidth;
			unsigned int blockHeight;
			unsigned int paddingWidth;
			unsigned int paddingHeight;
		};

		struct BlockDeepSegment
		{
			AtlasTexelIndex currentStart;
			AtlasTexelIndex currentEnd;
			AtlasTexelIndex previousStart;
			AtlasTexelIndex nextStart;
			AtlasInteriorTexelIndex deepStart;
		};

		struct BlockDeepSegmentedLine
		{
			std::vector< BlockDeepSegment > blockDeepSegments;
		};

		struct BlockTask
		{
			std::vector< BlockDeepSegmentedLine > blockPaddedSegmentedLines;
			std::vector< BlockDeepSegmentedLine > blockDeepSegmentedLines;
		};

		struct ThreadTask
		{
			int taskDeepTexels;
			std::vector< BlockTask > blockTasks;
		};

		bool threadTaskComparison( const ThreadTask & task1 , const ThreadTask & task2 ){ return task1.taskDeepTexels < task2.taskDeepTexels; }

		struct InteriorTexelToCellLine
		{
			AtlasTexelIndex texelStartIndex;
			AtlasTexelIndex texelEndIndex;
			AtlasInteriorTexelIndex coeffOffset;
			int length;
			AtlasCellIndex previousCellStartIndex;
			AtlasCellIndex     nextCellStartIndex;
		};


		struct ProlongationLine
		{
			AtlasTexelIndex startIndex;
			AtlasTexelIndex centerLineIndex;
			AtlasTexelIndex prevLineIndex;
			AtlasTexelIndex nextLineIndex;
			int length;
			bool alignedStart;
		};

		template< typename GeometryReal >
		using AuxiliaryNode = NodeInfo< GeometryReal , AtlasInteriorOrBoundaryNodeIndex >;

		enum struct CellType
		{
			Exterior ,		// No geometry passes through the cell
			Boundary ,		// A boundary edge passes through the cell
			Interior		// A triangle passes through the cell, but not a boundary edge
		};

		enum struct TexelType
		{
			Unsupported ,					// No geometry passes through the support of the node
			BoundarySupportedAndUncovered ,	// The support of the texel includes boundary nodes but the node is not covered by a triangle
			BoundarySupportedAndCovered ,	// The support of the texel includes boundary nodes and the node is covered by a triangle
			InteriorSupported				// The support of the texel consists of interior nodes
		};

		struct CellIndex
		{
			CellIndex( void ) : combined(-1) , interior(-1) , boundary(-1){}
			ChartCellIndex combined;
			ChartInteriorCellIndex interior;
			ChartBoundaryCellIndex boundary;
		};

		struct TexelIndex
		{
			TexelIndex( void ) : combined(-1) , interior(-1) , covered(-1){}
			AtlasTexelIndex combined;
			AtlasCoveredTexelIndex covered;
			AtlasInteriorTexelIndex interior;
		};

		bool IsCovered( TexelType texelType ){ return texelType==TexelType::BoundarySupportedAndCovered || texelType==TexelType::InteriorSupported; }
		bool IsBoundarySupported( TexelType texelType ){ return texelType==TexelType::BoundarySupportedAndCovered || texelType==TexelType::BoundarySupportedAndUncovered; }

		std::string CellTypeName( CellType cellType )
		{
			if     ( cellType==CellType::Exterior ) return std::string( "Exterior" );
			else if( cellType==CellType::Boundary ) return std::string( "Boundary" );
			else if( cellType==CellType::Interior ) return std::string( "Interior" );
			else
			{
				MK_THROW( "Unrecognized cell type" );
				return std::string();
			}
		}

		std::string TexelTypeName( TexelType texelType )
		{
			if     ( texelType==TexelType::Unsupported ) return std::string( "Unsupported" );
			else if( texelType==TexelType::BoundarySupportedAndCovered ) return std::string( "Boundary supported and covered" );
			else if( texelType==TexelType::BoundarySupportedAndUncovered ) return std::string( "Boundary supported and uncovered" );
			else if( texelType==TexelType::InteriorSupported ) return std::string( "Interior supported" );
			else
			{
				MK_THROW( "Unrecognized texel type" );
				return std::string();
			}
		}

		template< typename GeometryReal >
		struct GridChart
		{
			Point2D< GeometryReal > corner;
			unsigned int cornerCoords[2];
			int centerOffset[2];
			GeometryReal cellSizeW;
			GeometryReal cellSizeH;
			unsigned int width;
			unsigned int height;
			unsigned int atlasWidth;
			unsigned int atlasHeight;

			Point2D< GeometryReal > nodePosition( unsigned int i , unsigned int j ) const { return Point2D< GeometryReal >( i*cellSizeW , j*cellSizeH ); }

			AtlasGridVertexIndex nodeIndex( unsigned int i , unsigned int j ) const { return AtlasGridVertexIndex( atlasWidth*(j+cornerCoords[1]) + (i+cornerCoords[0]) ); }
			bool factorNodeIndex( AtlasGridVertexIndex g , unsigned int &i , unsigned int &j ) const
			{
				unsigned int idx = static_cast< unsigned int >( g );
				if( idx<atlasWidth*atlasHeight )
				{
					i = idx % atlasWidth , j = idx / atlasWidth;
					if( i<cornerCoords[0] || j<cornerCoords[1] || i>=cornerCoords[0]+width || j>=cornerCoords[1]+height ) return false;
					i -= cornerCoords[0] , j -= cornerCoords[1];
					return true;
				}
				else return false;
			}

			AtlasGridEdgeIndex edgeIndex( unsigned int i , unsigned int j , unsigned int dir ) const { return AtlasGridEdgeIndex( atlasWidth*atlasHeight*dir + (j+cornerCoords[1])*atlasWidth + (i+cornerCoords[0]) ); }

			bool factorEdgeIndex( AtlasGridEdgeIndex e , unsigned int &i , unsigned int &j , unsigned int &dir ) const
			{
				unsigned int idx = static_cast< unsigned int >( e );
				dir = 0;
				if( idx<atlasWidth*atlasHeight )
				{
					i = idx % atlasWidth , j = idx / atlasWidth;
					if( i<cornerCoords[0] || j<cornerCoords[1] || i>=cornerCoords[0]+width || j>=cornerCoords[1]+height ) return false;
					i -= cornerCoords[0] , j -= cornerCoords[1];
					return true;
				}
				else
				{
					idx -= atlasWidth*atlasHeight;
					dir++;
					if( idx<atlasWidth*atlasHeight )
					{
						i = idx % atlasWidth , j = idx / atlasWidth;
						if( i<cornerCoords[0] || j<cornerCoords[1] || i>=cornerCoords[0]+width || j>=cornerCoords[1]+height ) return false;
						i -= cornerCoords[0] , j -= cornerCoords[1];
						return true;
					}
					else return false;
				}
			}

			// Texel indices (atlas)
			RegularGrid< 2 , TexelIndex > texelIndices;

			// Cell indices (chart)
			RegularGrid< 2 , CellIndex > cellIndices;

			// The indices of the incident nodes
			ExplicitIndexVector< ChartCellIndex , BilinearElementIndex< AtlasTexelIndex > > combinedCellCombinedTexelBilinearElementIndices;

			// For interior cells, the indices of the incident interior nodes
			ExplicitIndexVector< ChartInteriorCellIndex , BilinearElementIndex< AtlasCoveredTexelIndex > > interiorCellCoveredTexelBilinearElementIndices;

			// For interior cells, the indices of the incident nodes
			ExplicitIndexVector< ChartInteriorCellIndex , BilinearElementIndex< AtlasTexelIndex > > interiorCellCombinedTexelBilinearElementIndices;

			// Maps converting boundary/interiorl cell indices to combined cell indices
			std::vector< ChartCellIndex > interiorCellIndexToCombinedCellIndex;
			std::vector< ChartCellIndex > boundaryCellIndexToCombinedCellIndex;

			const size_t numInteriorCells( void ) const { return interiorCellIndexToCombinedCellIndex.size(); }
			const size_t numBoundaryCells( void ) const { return boundaryCellIndexToCombinedCellIndex.size(); }

			ExplicitIndexVector< ChartBoundaryCellIndex , std::vector< BoundaryIndexedTriangle< GeometryReal > > > boundaryTriangles;

			ChartBoundaryTriangleIndex endBoundaryTriangleIndex;

			std::vector< AuxiliaryNode< GeometryReal > > auxiliaryNodes;

			RegularGrid< 2 , CellType > cellType;
			RegularGrid< 2 , TexelType > texelType;
			RegularGrid< 2 , AtlasMeshTriangleIndex > triangleID;
			RegularGrid< 2 , Point2D< GeometryReal > > barycentricCoords;

			AtlasInteriorCellIndex chartToAtlasInteriorCellIndex( ChartInteriorCellIndex idx ) const { return static_cast< AtlasInteriorCellIndex >( static_cast< unsigned int >(idx) + _interiorCellOffset ); }
			ChartInteriorCellIndex atlasToChartInteriorCellIndex( AtlasInteriorCellIndex idx ) const { return static_cast< ChartInteriorCellIndex >( static_cast< unsigned int >(idx) - _interiorCellOffset ); }
			AtlasCellIndex chartToAtlasCombinedCellIndex( ChartCellIndex idx ) const { return static_cast< AtlasCellIndex >( static_cast< unsigned int >(idx) + _combinedCellOffset ); }
			ChartCellIndex atlasToChartCombinedCellIndex( AtlasCellIndex idx ) const { return static_cast< ChartCellIndex >( static_cast< unsigned int >(idx) - _combinedCellOffset ); }

			// [WARNING] This should really be done through friendship
			void setInteriorCellOffset( unsigned int interiorCellOffset ){ _interiorCellOffset = interiorCellOffset; }
			void setCombinedCellOffset( unsigned int combinedCellOffset ){ _combinedCellOffset = combinedCellOffset; }
		protected:
			unsigned int _interiorCellOffset;
			unsigned int _combinedCellOffset;
		};

		struct TexelInfo
		{
			ChartIndex chartID;
			unsigned int ci;
			unsigned int cj;
			TexelType texelType;
		};

		struct TexelLineInfo
		{
			TexelLineInfo() : lineIndex(-1) , offset(-1) {}
			int lineIndex;
			int offset;
		};

		template< typename ... T > struct GridAtlas;

		template<>
		struct GridAtlas<>
		{
			struct IndexConverter
			{
				AtlasTexelIndex boundaryToCombined( AtlasBoundaryTexelIndex idx ) const { return _boundaryToCombined[idx]; }
				AtlasTexelIndex interiorToCombined( AtlasInteriorTexelIndex idx ) const { return _interiorToCombined[idx]; }
				AtlasBoundaryTexelIndex combinedToBoundary( AtlasTexelIndex idx ) const { return AtlasBoundaryTexelIndex( _combinedToBoundaryOrInterior[idx].first ? _combinedToBoundaryOrInterior[idx].second : static_cast< unsigned int >(-1) ); }
				AtlasInteriorTexelIndex combinedToInterior( AtlasTexelIndex idx ) const { return AtlasInteriorTexelIndex( _combinedToBoundaryOrInterior[idx].first ? static_cast< unsigned int >(-1) : _combinedToBoundaryOrInterior[idx].second ); }

				const ExplicitIndexVector< AtlasBoundaryTexelIndex , AtlasTexelIndex > &boundaryToCombined( void ) const { return _boundaryToCombined; }
				const ExplicitIndexVector< AtlasInteriorTexelIndex , AtlasTexelIndex > &interiorToCombined( void ) const { return _interiorToCombined; }

				size_t numCombined( void ) const { return _combinedToBoundaryOrInterior.size(); }
				size_t numBoundary( void ) const { return _boundaryToCombined.size(); }
				size_t numInterior( void ) const { return _interiorToCombined.size(); }
			protected:
				template< typename GeometryReal >
				friend void InitializeIndexConverter( const ExplicitIndexVector< ChartIndex , GridChart< GeometryReal > > & , AtlasTexelIndex , IndexConverter & );

				ExplicitIndexVector< AtlasTexelIndex , std::pair< bool , unsigned int > > _combinedToBoundaryOrInterior;
				ExplicitIndexVector< AtlasBoundaryTexelIndex , AtlasTexelIndex > _boundaryToCombined;
				ExplicitIndexVector< AtlasInteriorTexelIndex , AtlasTexelIndex > _interiorToCombined;
			};
		};

		template< typename GeometryReal , typename MatrixReal >
		struct GridAtlas< GeometryReal , MatrixReal >
		{
			std::vector< ThreadTask > threadTasks;
			typename GridAtlas<>::IndexConverter indexConverter;
			ExplicitIndexVector< AtlasTexelIndex , TexelInfo > texelInfo;
			ExplicitIndexVector< ChartIndex , GridChart< GeometryReal > > gridCharts;
			std::vector< SegmentedRasterLine > segmentedLines;
			std::vector< RasterLine > rasterLines;
			std::vector< RasterLine > restrictionLines;
			std::vector< DeepLine > deepLines;
			std::vector< ProlongationLine > prolongationLines;
			Eigen::SparseMatrix< MatrixReal > coarseToFineNodeProlongation;

			AtlasTexelIndex endCombinedTexelIndex;
			AtlasCoveredTexelIndex endCoveredTexelIndex;
			AtlasInteriorTexelIndex endInteriorTexelIndex;
			AtlasBoundaryTexelIndex endBoundaryTexelIndex;
			AtlasCellIndex endCombinedCellIndex;
			AtlasBoundaryCellIndex endBoundaryCellIndex;
			AtlasInteriorCellIndex endInteriorCellIndex;

			AtlasRefinedBoundaryVertexIndex endBoundaryVertexIndex;
			BoundaryMidPointIndex endMidPointIndex;
			unsigned int numFineNodes;
		};

		struct BoundaryDeepIndex
		{
			AtlasBoundaryTexelIndex boundaryIndex;
			AtlasTexelIndex combinedIndex;
			AtlasInteriorTexelIndex interiorIndex;
			unsigned int offset;
		};

		template< typename Real >
		struct BoundaryBoundaryIndex
		{
			AtlasBoundaryTexelIndex coarsePrincipalBoundaryIndex;
			AtlasBoundaryTexelIndex coarseSecondaryBoundaryIndex;
			AtlasInteriorTexelIndex fineInteriorIndex;
			unsigned int offset;
			Real weight;
		};

		template< typename GeometryReal , typename MatrixReal >
		struct HierarchicalSystem
		{
			std::vector< GridAtlas< GeometryReal , MatrixReal > > gridAtlases;

			std::vector< SparseMatrix< MatrixReal , int > > boundaryRestriction;
			std::vector< SparseMatrix< MatrixReal , int > > prolongation;
			std::vector< SparseMatrix< MatrixReal , int > > boundaryCoarseFineProlongation;
			std::vector< SparseMatrix< MatrixReal , int > > boundaryFineCoarseRestriction;
			std::vector< ExplicitIndexVector< AtlasBoundaryTexelIndex , BoundaryBoundaryIndex< MatrixReal > > > boundaryBoundaryIndices;
			std::vector< ExplicitIndexVector< AtlasBoundaryTexelIndex , BoundaryDeepIndex > > boundaryDeepIndices;
		};


		template< typename Real >
		struct SystemCoefficients
		{
			std::vector< Real > deepCoefficients;
			SparseMatrix< Real , int > boundaryDeepMatrix;
			SparseMatrix< Real , int > boundaryBoundaryMatrix;
		};

		template< typename DataType >
		struct MultigridLevelVariables
		{
			std::vector< DataType > x;
			std::vector< DataType > rhs;
			std::vector< DataType > residual;
			std::vector< DataType > boundary_rhs;
			std::vector< DataType > boundary_value;
			std::vector< DataType > variable_boundary_value;
		};

		template< typename Real >
		struct MultigridLevelIndices
		{
			std::vector< ThreadTask > threadTasks;
			ExplicitIndexVector< AtlasBoundaryTexelIndex , AtlasTexelIndex > boundaryToCombined;
			std::vector< SegmentedRasterLine > segmentedLines;
			std::vector< RasterLine > rasterLines;
			std::vector< RasterLine > restrictionLines;
			std::vector< ProlongationLine > prolongationLines;
			SparseMatrix< Real , int > boundaryRestriction;
		};

		template< typename DirectSolver >
		struct VCycleSolvers
		{
			std::vector< DirectSolver > boundary;
			DirectSolver coarse;
		};
	}
}
#include "SparseMatrixParser.h"
#include "Atlas.h"
#include "ConstrainedTriangulation.h"
#include "PolygonClipping.h"


namespace MishaK
{
	namespace TSP
	{
#include "Metric.inl"
#include "BoundaryQuadraticElements.inl"
#include "ProlongationAndRestriction.inl"
#include "HierarchyConstruction.inl"
#include "IterativeSolvers.inl"
	}
}

#endif// HIERARCHICAL_SYSTEM_INCLUDED
