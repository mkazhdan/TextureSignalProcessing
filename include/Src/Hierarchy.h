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
#include <Misha/Image.h>
#include <Misha/LinearSolvers.h>
#include <Misha/RegularGrid.h>
#include <Misha/Miscellany.h>
#include <Misha/Exceptions.h>
#ifdef USE_RASTERIZER
#include <Misha/Texels.h>
#include <Misha/Rasterizer2D.h>
#include <Misha/Geometry.h>
#include <Misha/MultiThreading.h>
#include <Misha/Atomic.h>
#include <Misha/Atomic.Geometry.h>
#endif // USE_RASTERIZER
#include "IndexedPolygon.h"
#ifdef PRE_CLIP_TRIANGLES
#include "PolygonClipping.h"
#endif // PRE_CLIP_TRIANGLES

namespace MishaK
{
	template< typename Data > using Image = RegularGrid< 2 , Data >;

	template< unsigned int BitDepth , typename Data >
	void WriteImage( const Image< Data > &image , std::string fileName )
	{
		using CType = typename ImageChannel< BitDepth >::Type;
		static const CType Scale = ~((CType )0);
		if constexpr( std::is_same_v< Data , Point3D< double > > )
		{
			double scale = (double)Scale;
			CType * pixels = new CType[ image.res(0)*image.res(1)*3 ];
			for( unsigned int i=0 ; i<image.res(0) ; i++ ) for( unsigned int j=0 ; j<image.res(1) ; j++ ) for( int c=0 ; c<3 ; c++ ) pixels[ 3*(j*image.res(0)+i)+c ] = (CType)std::min< long long >( Scale , std::max< long long >( 0 , (long long)( image(i,j)[c] * scale + 0.5 ) ) );
			ImageWriter< BitDepth >::Write( fileName , pixels , image.res(0) , image.res(1) , 3 );
			delete[] pixels;
		}
		else if constexpr( std::is_same_v< Data , Point3D< float > > )
		{
			float scale = (float)Scale;
			CType * pixels = new CType[ image.res(0) * image.res(1) * 3 ];
			for( unsigned int i=0 ; i<image.res(0) ; i++ ) for( unsigned int j=0 ; j<image.res(1) ; j++ ) for( int c=0 ; c<3 ; c++ ) pixels[ 3*(j*image.res(0)+i)+c ] = (CType)std::min< long long >( Scale , std::max< long long >( 0 , (long long)( image(i,j)[c] * scale + 0.5f ) ) );
			ImageWriter< BitDepth >::Write( fileName , pixels , image.res(0) , image.res(1) , 3 );
			delete[] pixels;
		}
		else if constexpr( std::is_same_v< Data , Point3D< CType > > )
		{
			ImageWriter< BitDepth >::Write( fileName , (const CType*)image() , image.res(0) , image.res(1) , 3 );
		}
		else MK_THROW( "Bad data type " );
	}

	template< unsigned int BitDepth , typename Data >
	void ReadImage( Image< Data > &image , std::string fileName )
	{
		using CType = typename ImageChannel< BitDepth >::Type;
		static const CType Scale = ~((CType )0);
		unsigned int width , height;
		CType * pixels = ImageReader< BitDepth >::ReadColor( fileName , width , height );
		if( !pixels ) MK_THROW( "Failed to read image: " , fileName );
		image.resize( width , height );

		if constexpr( std::is_same_v< Data , Point3D< double > > )
		{
			double scale = (double)Scale;
			for( int i=0 ; i<(int)width ; i++ ) for( int j=0 ; j<(int)height ; j++ ) for( int c=0 ; c<3 ; c++ ) image(i,j)[c] = ( pixels[ (j*width+i)*3+c ] ) / scale;
		}
		else if constexpr( std::is_same_v< Data , Point3D< float > > )
		{
			float scale = (float)Scale;
			for( int i=0 ; i<(int)width ; i++ ) for( int j=0 ; j<(int)height ; j++ ) for( int c=0 ; c<3 ; c++ ) image(i,j)[c] = ( pixels[ (j*width+i)*3+c ] ) / scale;
		}
		else if constexpr( std::is_same_v< Data , Point3D< CType > > )
		{
			memcpy( (CType*)image() , pixels , sizeof( CType ) * image.res(0) * image.res(1) * 3 );
		}
		else MK_THROW( "Bad data type " );

		delete[] pixels;
	}

	template< typename Real >
	struct BilinearElementScalarSample
	{
		struct SampleData
		{
			SampleData( void ){ memset( dualValues , 0 , sizeof(Real) * 4 ); }
			Point2D< Real > pos;			// The sample position inside the element
			Real dualValues[4];				// The integrated values of the four incident bilinear basis functions, dualized
		};
		SquareMatrix< Real , 2 > tensor;	// The inverse metric tensor defined by the intersecting triangle
		unsigned int cellOffset;
		bool operator < ( const BilinearElementScalarSample& sample ) const { return cellOffset<sample.cellOffset; }

		BilinearElementScalarSample( void ) : _sampleNum(0) , _samples(NULL) {}
		BilinearElementScalarSample( unsigned int sz ) : _sampleNum(0) , _samples(NULL) { resize(sz); }
		BilinearElementScalarSample( const BilinearElementScalarSample& bilinearElementScalarSample ) : _sampleNum(0) , _samples(NULL)
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
			if( _samples ){ delete[] _samples ; _samples = NULL; }
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
	struct QuadraticElementScalarSample
	{
		struct SampleData
		{
			SampleData( void ){ memset( dualValues , 0 , sizeof(Real) * 6 ); }
			Point2D< Real > pos;			// The sample position inside the element
			Real dualValues[6];				// The integrated values of the six incident quadratic basis functions, dualized
		};
		SquareMatrix< Real , 2 > tensor;	// The inverse metric tensor defined by the intersecting triangle
		AtlasInteriorOrBoundaryNodeIndex fineNodes[6];

		QuadraticElementScalarSample( void ) : _sampleNum(0) , _samples(NULL) {}
		QuadraticElementScalarSample( unsigned int sz ) : _sampleNum(0) , _samples(NULL) { resize(sz); }
		QuadraticElementScalarSample( const QuadraticElementScalarSample& quadraticElementScalarSample ) : _sampleNum(0) , _samples(NULL)
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
			if( _samples ){ delete[] _samples ; _samples = NULL; }
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
	struct BilinearElementGradientSample
	{
		struct SampleData
		{
			Point2D< Real > pos;				// The sample position inside the element
			Point2D< Real > dualGradients[4];	// The integrated gradients of the four incident bilinear basis functions, dualized
		};
		SquareMatrix< Real , 2 > tensor;		// The inverse metric tensor defined by the intersecting triangle
		unsigned int cellOffset;
		bool operator < ( const BilinearElementGradientSample& sample ) const { return cellOffset < sample.cellOffset; }

		BilinearElementGradientSample( void ) : _sampleNum(0) , _samples(NULL) {}
		BilinearElementGradientSample( unsigned int sz ) : _sampleNum(0) , _samples(NULL) { resize(sz); }
		BilinearElementGradientSample( const BilinearElementGradientSample& bilinearElementGradientSample ) : _sampleNum(0) , _samples(NULL)
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
			if( _samples ){ delete[] _samples ; _samples = NULL; }
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
			Point2D< Real > pos;				// The sample position inside the element
			Point2D< Real > dualGradients[6];	// The integrated gradients of the six incident quadratic basis functions, dualized
		};
		SquareMatrix< Real , 2 > tensor;		// The inverse metric tensor defined by the intersecting triangle
		AtlasInteriorOrBoundaryNodeIndex fineNodes[6];

		QuadraticElementGradientSample( void ) : _sampleNum(0) , _samples(NULL) {}
		QuadraticElementGradientSample( unsigned int sz ) : _sampleNum(0) , _samples(NULL) { resize(sz); }
		QuadraticElementGradientSample( const QuadraticElementGradientSample& quadraticElementGradientSample ) : _sampleNum(0) , _samples(NULL)
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
			if( _samples ){ delete[] _samples ; _samples = NULL; }
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
		int lineStartIndex;
		int lineEndIndex;
		int prevLineIndex;
		int nextLineIndex;
		int coeffStartIndex;
	};

	struct DeepLine
	{
		int coarseLineStartIndex;
		int coarseLineEndIndex;
		int finePrevLineIndex;
		int fineCurrentLineIndex;
		int fineNextLineIndex;
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
		int currentStart;
		int currentEnd;
		int previousStart;
		int nextStart;
		int deepStart;
	};

	struct BlockDeepSegmentedLine
	{
		std::vector<BlockDeepSegment> blockDeepSegments;
	};

	struct BlockTask
	{
		std::vector<BlockDeepSegmentedLine> blockPaddedSegmentedLines;
		std::vector<BlockDeepSegmentedLine> blockDeepSegmentedLines;
	};

	struct ThreadTask
	{
		int taskDeepTexels;
		std::vector< BlockTask > blockTasks;
	};

	bool threadTaskComparison( const ThreadTask & task1 , const ThreadTask & task2 ){ return task1.taskDeepTexels < task2.taskDeepTexels; }

	struct InteriorTexelToCellLine
	{
		int texelStartIndex;
		int texelEndIndex;
		int coeffOffset;
		int length;
		AtlasCombinedCellIndex previousCellStartIndex;
		AtlasCombinedCellIndex     nextCellStartIndex;
	};


	struct ProlongationLine
	{
		int startIndex;
		int length;
		int centerLineIndex;
		int prevLineIndex;
		int nextLineIndex;
		bool alignedStart;
	};

	struct RestrictionLine
	{
		RestrictionLine( void ) { startIndex = length = centerLineIndex = prevLineIndex = nextLineIndex = -1; }
		int startIndex;
		int length;
		int centerLineIndex;
		int prevLineIndex;
		int nextLineIndex;
		int outputStart;
	};

	template< typename GeometryReal >
	struct AuxiliaryNode
	{
		Point2D< GeometryReal > position;
		AtlasInteriorOrBoundaryNodeIndex index;

		AuxiliaryNode( Point2D< GeometryReal > position=Point2D< GeometryReal >() , AtlasInteriorOrBoundaryNodeIndex index=static_cast< AtlasInteriorOrBoundaryNodeIndex >(-1) ) : position(position) , index(index) {}
	};

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
		ChartCombinedCellIndex combined;
		ChartInteriorCellIndex interior;
		ChartBoundaryCellIndex boundary;
	};

	struct TexelIndex
	{
		TexelIndex( void ) : combined(-1) , interior(-1) , covered(-1){}
		unsigned int combined , interior , covered;
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
		int cornerCoords[2];
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

		// Texel indices (global)
		Image< TexelIndex > texelIndices;

		// Cell indices (local)
		Image< CellIndex > cellIndices;

		// The indices of the incident nodes
		std::vector< BilinearElementIndex > combinedCellCombinedBilinearElementIndices;

		// For interior cells, the indices of the incident interior nodes
		std::vector< BilinearElementIndex > interiorCellInteriorBilinearElementIndices;

		// For interior cells, the indices of the incident nodes
		std::vector< BilinearElementIndex > interiorCellCombinedBilinearElementIndices;

		// Maps converting boundary/interiorl cell indices to combined cell indices
		std::vector< ChartCombinedCellIndex > interiorCellIndexToCombinedCellIndex;
		std::vector< ChartCombinedCellIndex > boundaryCellIndexToCombinedCellIndex;

		const size_t numInteriorCells( void ) const { return interiorCellIndexToCombinedCellIndex.size(); }
		const size_t numBoundaryCells( void ) const { return boundaryCellIndexToCombinedCellIndex.size(); }

		IndexVector< ChartBoundaryCellIndex , std::vector< IndexedPolygon< GeometryReal > > > boundaryPolygons;
		std::vector< std::vector< BoundaryIndexedTriangle< GeometryReal > > > boundaryTriangles;

		unsigned int numBoundaryTriangles;

		std::vector< AuxiliaryNode< GeometryReal > > auxiliaryNodes;

		Image< CellType > cellType;
		Image< TexelType > texelType;
		Image< AtlasMeshTriangleIndex > triangleID;
		Image< Point2D< GeometryReal > > barycentricCoords;
#ifdef PRE_CLIP_TRIANGLES
		Image< std::vector< std::pair< ChartMeshTriangleIndex , CellClippedTriangle< GeometryReal > > > > clippedTriangles;
#endif // PRE_CLIP_TRIANGLES

		AtlasInteriorCellIndex chartToAtlasInteriorCellIndex( ChartInteriorCellIndex idx ) const { return static_cast< AtlasInteriorCellIndex >( static_cast< unsigned int >(idx) + _interiorCellOffset ); }
		ChartInteriorCellIndex atlasToChartInteriorCellIndex( AtlasInteriorCellIndex idx ) const { return static_cast< ChartInteriorCellIndex >( static_cast< unsigned int >(idx) - _interiorCellOffset ); }
		AtlasCombinedCellIndex chartToAtlasCombinedCellIndex( ChartCombinedCellIndex idx ) const { return static_cast< AtlasCombinedCellIndex >( static_cast< unsigned int >(idx) + _combinedCellOffset ); }
		ChartCombinedCellIndex atlasToChartCombinedCellIndex( AtlasCombinedCellIndex idx ) const { return static_cast< ChartCombinedCellIndex >( static_cast< unsigned int >(idx) - _combinedCellOffset ); }

		// [WARNING] This should really be done through friendship
		void setInteriorCellOffset( unsigned int interiorCellOffset ){ _interiorCellOffset = interiorCellOffset; }
		void setCombinedCellOffset( unsigned int combinedCellOffset ){ _combinedCellOffset = combinedCellOffset; }
	protected:
		unsigned int _interiorCellOffset;
		unsigned int _combinedCellOffset;
	};

	struct GridNodeInfo
	{
		int chartID;
		int ci;
		int cj;
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
			unsigned int boundaryToSupported( unsigned int idx ) const { return _boundaryToSupported[idx]; }
			unsigned int     deepToSupported( unsigned int idx ) const { return     _deepToSupported[idx]; }
			unsigned int supportedToBoundary( unsigned int idx ) const { return _supportedToBoundaryOrDeep[idx].first ? _supportedToBoundaryOrDeep[idx].second : static_cast< unsigned int >(-1); }
			unsigned int supportedToDeep    ( unsigned int idx ) const { return _supportedToBoundaryOrDeep[idx].first ? static_cast< unsigned int >(-1) : _supportedToBoundaryOrDeep[idx].second; }

			const std::vector< unsigned int > &boundaryToSupported( void ) const { return _boundaryToSupported; }
			const std::vector< unsigned int > &    deepToSupported( void ) const { return     _deepToSupported; }

			size_t numSupported( void ) const { return _supportedToBoundaryOrDeep.size(); }
			size_t numBoundary( void ) const { return _boundaryToSupported.size(); }
			size_t numDeep( void ) const { return _deepToSupported.size(); }
		protected:
			template< typename GeometryReal >
			friend void InitializeIndexConverter( const std::vector< GridChart< GeometryReal > > & , unsigned int , IndexConverter & );

			std::vector< std::pair< bool , unsigned int > > _supportedToBoundaryOrDeep;
			std::vector< unsigned int > _boundaryToSupported;
			std::vector< unsigned int >     _deepToSupported;
		};
	};

	template< typename GeometryReal , typename MatrixReal >
	struct GridAtlas< GeometryReal , MatrixReal >
	{
		std::vector< ThreadTask > threadTasks;
		typename GridAtlas<>::IndexConverter indexConverter;
		std::vector< GridNodeInfo > nodeInfo;
		std::vector< GridChart< GeometryReal > > gridCharts;
		std::vector< SegmentedRasterLine > segmentedLines;
		std::vector< RasterLine > rasterLines;
		std::vector< RasterLine > restrictionLines;
		std::vector< DeepLine > deepLines;
		std::vector< ProlongationLine > prolongationLines;
		Eigen::SparseMatrix< MatrixReal > coarseToFineNodeProlongation;

		unsigned int numTexels;
		unsigned int numInteriorTexels;
		unsigned int numDeepTexels;
		unsigned int numBoundaryTexels;
		unsigned int numCells;
		unsigned int numBoundaryCells;
		unsigned int numInteriorCells;
		unsigned int numBoundaryNodes;
		unsigned int numMidPoints;
		unsigned int numFineNodes;
	};

	struct BoundaryDeepIndex
	{
		int boundaryIndex;
		int deepGlobalIndex;
		int deepIndex;
		int offset;
	};

	template< typename Real >
	struct BoundaryBoundaryIndex
	{
		int coarsePrincipalBoundaryIndex;
		int coarseSecondaryBoundaryIndex;
		int fineDeepIndex;
		int offset;
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
		std::vector< std::vector< BoundaryBoundaryIndex< MatrixReal > > > boundaryBoundaryIndices;
		std::vector< std::vector< BoundaryDeepIndex > > boundaryDeepIndices;
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
		std::vector<ThreadTask> threadTasks;
		std::vector< unsigned int > boundaryToSupported;
		std::vector<SegmentedRasterLine> segmentedLines;
		std::vector<RasterLine> rasterLines;
		std::vector<RasterLine> restrictionLines;
		std::vector<ProlongationLine> prolongationLines;
		SparseMatrix<Real, int> boundaryRestriction;
	};

	template< typename DirectSolver >
	struct VCycleSolvers
	{
		std::vector< DirectSolver > boundary;
		DirectSolver coarse;
	};
}

#include "SparseMatrixParser.h"
#include "Atlas.h"
#include "ConstrainedTriangulation.h"
#include "MapLoop.h"

namespace MishaK
{
#include "Metric.inl"
#include "BoundaryQuadraticElements.inl"
#include "ProlongationAndRestriction.inl"
#include "HierarchyConstruction.inl"
#include "UpdateCoefficients.inl"
#include "IterativeSolvers.inl"

}

#endif//HIERARCHICAL_SYSTEM_INCLUDED
