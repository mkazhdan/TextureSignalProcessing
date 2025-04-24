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
#pragma once

template< typename GeometryReal >
void InitializeIndexConverter
(
	const std::vector< GridChart< GeometryReal > > &gridCharts ,
	unsigned int numSupportedTexels ,
	typename GridAtlas<>::IndexConverter & indexConverter
) 
{
	indexConverter._supportedToBoundaryOrDeep.resize( numSupportedTexels );
	unsigned int globalIndex = 0;
	unsigned int boundaryIndex = 0;
	unsigned int deepIndex = 0;

	for( const GridChart< GeometryReal > &gridChart : gridCharts ) for( size_t i=0 ; i<gridChart.texelType.size() ; i++ ) if( gridChart.texelType[i]!=TexelType::Unsupported )
	{
		if( IsBoundarySupported( gridChart.texelType[i] ) )
		{
			indexConverter._boundaryToSupported.push_back( globalIndex );
			indexConverter._supportedToBoundaryOrDeep[ globalIndex ] = std::pair< bool , unsigned int >( true , boundaryIndex++ );
		}
		else if( gridChart.texelType[i]==TexelType::InteriorSupported )
		{
			indexConverter._deepToSupported.push_back( globalIndex );
			indexConverter._supportedToBoundaryOrDeep[ globalIndex ] = std::pair< bool , unsigned int >( false , deepIndex++ );
		}
		if( gridChart.texelIndices[i].combined!=globalIndex ) MK_THROW( "Unexpected global index: actual " , gridChart.texelIndices[i].combined , " , expected " , globalIndex );
		globalIndex++;
	}
}

//Node type : inactive(-1) , exterior (0), interior boundary (1), interior deep (2) hybryd (both deep and boundary for the solver)(3).
//Cell type : inactive(-1) , boundary (0), interior (1).
template< typename GeometryReal >
void InitializeGridChartsActiveNodes
(
	unsigned int chartID ,
	const AtlasChart< GeometryReal > & atlasChart ,
	GridChart< GeometryReal > & gridChart ,
	std::vector< GridNodeInfo > & nodeInfo ,
	std::vector< RasterLine > & rasterLines ,
	std::vector< SegmentedRasterLine > & segmentedLines ,
	std::vector< ThreadTask > & threadTasks ,
	unsigned int & combinedTexelIndex ,
	unsigned int & interiorTexelIndex ,
	unsigned int &     deepTexelIndex ,
	unsigned int & boundaryTexelIndex ,
	unsigned int & combinedCellIndex ,
	unsigned int & boundaryCellIndex ,
	unsigned int & interiorCellIndex ,
	const MultigridBlockInfo & multigridBlockInfo
)
{
	unsigned int width = gridChart.width;
	unsigned int height = gridChart.height;
	GeometryReal cellSizeW = gridChart.cellSizeW;
	GeometryReal cellSizeH = gridChart.cellSizeH;
	Image< TexelType > &texelType = gridChart.texelType;
	texelType.resize( width , height );
	for( int i=0 ; i<texelType.size() ; i++ ) texelType[i] = TexelType::Unsupported;
#ifdef DEBUG_ATLAS
	Image< int > nodeOwner;
	nodeOwner.resize( width , height );
#endif // DEBUG_ATLAS

	Image< CellType > &cellType = gridChart.cellType;
	cellType.resize( width-1 , height-1 );
	for( int i=0 ; i<cellType.size() ; i++ ) cellType[i] = CellType::Exterior;

#ifdef NEW_INDEXING
	Image< ChartMeshTriangleIndex > chartTriangleID;
	chartTriangleID.resize( width , height );
	for( unsigned int i=0 ; i<chartTriangleID.size() ; i++ ) chartTriangleID[i] = ChartMeshTriangleIndex( -1 );
#else // !NEW_INDEXING
	Image< unsigned int > &triangleID = gridChart.triangleID;
	triangleID.resize( width , height );
	for ( int i=0 ; i<triangleID.size() ; i++ ) triangleID[i] = -1;
#endif // NEW_INDEXING

	Image< Point2D< GeometryReal > > &barycentricCoords = gridChart.barycentricCoords;
	barycentricCoords.resize( width , height );

#ifdef PRE_CLIP_TRIANGLES
#ifdef NEW_INDEXING
	Image< std::vector< std::pair< ChartMeshTriangleIndex , CellClippedTriangle< GeometryReal > > > > &clippedTriangles = gridChart.clippedTriangles;
#else // !NEW_INDEXING
	Image< std::vector< std::pair< unsigned int , CellClippedTriangle< GeometryReal > > > > &clippedTriangles = gridChart.clippedTriangles;
#endif // NEW_INDEXING
	clippedTriangles.resize( width-1 , height-1 );
#endif // PRE_CLIP_TRIANGLES

#ifdef USE_RASTERIZER
	using Range = RegularGrid< 2 >::Range;
	using Index = RegularGrid< 2 >::Index;
	Range nodeRange , cellRange;
	nodeRange.second[0] = width , cellRange.second[0] = width-1;
	nodeRange.second[1] = height , cellRange.second[1] = height-1;

#ifdef NEW_INDEXING
	auto GetTriangle = [&]( ChartMeshTriangleIndex t , bool textureSpace=true )
		{
			Simplex< double , 2 , 2 > simplex;
			SimplexIndex< 2 , ChartMeshVertexIndex > tri = atlasChart.triangleIndex(t);
			for( unsigned int k=0 ; k<=2 ; k++ )
			{
				simplex[k] = atlasChart.vertex( tri[k] ) - gridChart.corner;
				if( textureSpace ) simplex[k][0] /= cellSizeW , simplex[k][1] /= cellSizeH;
			}
			return simplex;
		};
#else // !NEW_INDEXING
	auto GetTriangle = [&]( unsigned int t , bool textureSpace=true )
		{
			Simplex< double , 2 , 2 > simplex;
			for( unsigned int i=0 ; i<=2 ; i++ )
			{
				simplex[i] = atlasChart.vertices[ atlasChart.triangles[t][i] ] - gridChart.corner;
				if( textureSpace ) simplex[i][0] /= cellSizeW , simplex[i][1] /= cellSizeH;
			}
			return simplex;
		};
#endif // NEW_INDEXING

#ifdef NEW_INDEXING
	auto GetEdge = [&]( ChartMeshHalfEdgeIndex he , bool textureSpace=true )
#else // !NEW_INDEXING
	auto GetEdge = [&]( unsigned int he , bool textureSpace=true )
#endif // NEW_INDEXING
		{
			Simplex< double , 2 , 1 > simplex;
#ifdef NEW_INDEXING
			SimplexIndex< 1 , ChartMeshVertexIndex > eIndex = atlasChart.edgeIndex( he );
#else // !NEW_INDEXING
			SimplexIndex< 1 > eIndex = atlasChart.edgeIndex( he );
#endif // NEW_INDEXING
			for( unsigned int i=0 ; i<2 ; i++ )
			{
#ifdef NEW_INDEXING
				simplex[i] = atlasChart.vertex( eIndex[i] ) - gridChart.corner;
#else // !NEW_INDEXING
				simplex[i] = atlasChart.vertices[ eIndex[i] ] - gridChart.corner;
#endif // NEW_INDEXING
				if( textureSpace ) simplex[i][0] /= cellSizeW , simplex[i][1] /= cellSizeH;
			}
			return simplex;
		};

	//(1) Add nodes covered by the triangles
#else // !USE_RASTERIZER
	//(1) Add interior texels
#endif // USE_RASTERIZER

#ifdef USE_RASTERIZER
	// Rasterize triangles into nodes/cells
#ifdef NEW_INDEXING
	for( unsigned int t=0 ; t<atlasChart.numTriangles() ; t++ )
#else // !NEW_INDEXING
	for( unsigned int t=0 ; t<atlasChart.triangles.size() ; t++ )
#endif // NEW_INDEXING
	{
		// Compute the associated triangle in (shifted) texel coordinates
#ifdef NEW_INDEXING
		Simplex< double , 2 , 2 > simplex = GetTriangle( ChartMeshTriangleIndex(t) );
#else // !NEW_INDEXING
		Simplex< double , 2 , 2 > simplex = GetTriangle( t );
#endif // NEW_INDEXING

		auto Kernel = [&]( Index I )
			{
				Point3D< double > bc = simplex.barycentricCoordinates( Point2D< double >( I[0] , I[1] ) );
				Point2D< GeometryReal > newBC = Point2D< GeometryReal >( Point2D< double >( bc[1] , bc[2] ) );

				texelType(I) = TexelType::BoundarySupported;

#ifdef NEW_INDEXING
				if( chartTriangleID(I)==ChartMeshTriangleIndex(-1) )
				{
					chartTriangleID(I) = ChartMeshTriangleIndex(t);
#else // !NEW_INDEXING
				if( triangleID(I)==-1 )
				{
					triangleID(I) = t;
#endif // NEW_INDEXING
					barycentricCoords(I) = newBC;
				}
				else
				{
					Point2D< GeometryReal > oldBC = barycentricCoords(I);
#ifdef NEW_INDEXING
					Simplex< double , 2 , 2 > oldSimplex = GetTriangle( chartTriangleID(I) );
#else // !NEW_INDEXING
					Simplex< double , 2 , 2 > oldSimplex = GetTriangle( triangleID(I) );
#endif // NEW_INDEXING
					Point< double , 2 > oldP = oldSimplex( oldBC ) , newP = simplex( newBC );
					double oldD2 = Point< double , 2 >::SquareDistance( oldP , oldSimplex.nearest(oldP) );
					double newD2 = Point< double , 2 >::SquareDistance( newP ,    simplex.nearest(newP) );
					if( newD2<oldD2 )
					{
#ifdef NEW_INDEXING
						chartTriangleID(I) = ChartMeshTriangleIndex(t);
#else // !NEW_INDEXING
						triangleID(I) = t;
#endif // NEW_INDEXING
						barycentricCoords(I) = newBC;
					}
				}
			};
		// Process texels whose support overlaps a triangle
		Rasterizer2D::RasterizeSupports< false , false >( simplex , Kernel , nodeRange );

#ifdef PRE_CLIP_TRIANGLES
		{
#ifdef NEW_INDEXING
			Simplex< double , 2 , 2 > normalizedSimplex = GetTriangle( ChartMeshTriangleIndex(t) , false );
#else // !NEW_INDEXING
			Simplex< double , 2 , 2 > normalizedSimplex = GetTriangle( t , false );
#endif // NEW_INDEXING

			auto Kernel = [&]( Index I )
				{
					cellType(I) = CellType::Interior;
					CellClippedTriangle poly( normalizedSimplex );
#ifdef NEW_INDEXING
					if( ClipTriangleToPrimalCell( poly , I[0] , I[1] , gridChart.cellSizeW , gridChart.cellSizeH ) ) clippedTriangles(I).push_back( std::make_pair( ChartMeshTriangleIndex(t) , poly ) );
#else // !NEW_INDEXING
					if( ClipTriangleToPrimalCell( poly , I[0] , I[1] , gridChart.cellSizeW , gridChart.cellSizeH ) ) clippedTriangles(I).push_back( std::make_pair( t, poly ) );
#endif // NEW_INDEXING
					else MK_THROW( "Expected triangle to intersect cell" );
				};
			// Process cells
			Rasterizer2D::RasterizeSupports< true , true >( simplex , Kernel , cellRange );
		}
#else // !PRE_CLIP_TRIANGLES
		Rasterizer2D::RasterizeSupports< true , true >( simplex , [&]( Index I ){ cellType(I)=CellType::Interior; } , cellRange );
#endif // PRE_CLIP_TRIANGLES
	}

	// Over-write the node designation for nodes covered by a triangle
#ifdef NEW_INDEXING
	for( int t=0 ; t<atlasChart.numTriangles() ; t++ )
		Rasterizer2D::RasterizeNodes< false >( GetTriangle( ChartMeshTriangleIndex(t) ) , [&]( Index I ){ texelType(I) = TexelType::BoundarySupportedAndCovered; } , nodeRange );
#else // !NEW_INDEXING
	for( int t=0 ; t<atlasChart.triangles.size() ; t++ )
		Rasterizer2D::RasterizeNodes< false >( GetTriangle( t ) , [&]( Index I ){ texelType(I) = TexelType::BoundarySupportedAndCovered; } , nodeRange );
#endif // NEW_INDEXING

	// Rasterize boundary edges into cells
	for( unsigned int e=0 ; e<atlasChart.boundaryHalfEdges.size() ; e++ )
	{
#ifdef NEW_INDEXING
		unsigned int tIndex = static_cast< unsigned int >( atlasChart.boundaryHalfEdges[e] ) / 3;
		unsigned int kIndex = static_cast< unsigned int >( atlasChart.boundaryHalfEdges[e] ) % 3;

		Simplex< double , 2 , 2 > simplex = GetTriangle( ChartMeshTriangleIndex(tIndex) );
#else // !NEW_INDEXING
		int tIndex = atlasChart.boundaryHalfEdges[e] / 3;
		int kIndex = atlasChart.boundaryHalfEdges[e] % 3;

		Simplex< double , 2 , 2 > simplex = GetTriangle( tIndex );
#endif // NEW_INDEXING
		Simplex< double , 2 , 1 > subSimplex = GetEdge( atlasChart.boundaryHalfEdges[e] );
		Rasterizer2D::RasterizeSupports< true , true >( subSimplex , [&]( Index I ){ cellType(I)=CellType::Boundary; } , cellRange );
	}

	// Over-write the node designation for nodes whose support is entirely interior
	{
		auto Kernel = [&]( Index I )
			{
#ifdef NEW_INDEXING
				if( chartTriangleID(I)!=ChartMeshTriangleIndex(-1) )
#else // !NEW_INDEXING
				if( triangleID(I)!=-1 )
#endif // NEW_INDEXING
				{
					unsigned int bCount = 0;
					Range::Intersect( cellRange , Range::CellsSupportedOnNode(I) ).process( [&]( Index I ){ if( cellType(I)==CellType::Boundary ) bCount++; } );
					if( !bCount ) texelType(I) = TexelType::InteriorSupported;
				}
			};
		nodeRange.process( Kernel );
	}
#else // !USE_RASTERIZER

	for( int t=0 ; t<atlasChart.triangles.size() ; t++ )
	{
#ifdef USE_RASTERIZER
		// Compute the associated triangle in (shifted) texel coordinates
		Simplex< double , 2 , 2 > simplex = GetTriangle( t );
		auto Kernel = [&]( Index I )
			{
#ifdef DEBUG_ATLAS
				if( texelType(I)!=TexelType::Unassigned ) MK_WARN_ONCE( "Texel ( " , I[0]+gridChart.cornerCoords[0] , " , " , I[1]+gridChart.cornerCoords[1] , " ) owned by multiple triangles (mapping likely not injective): " , atlasChart.triangles[t]() , " " , nodeOwner(I) );
#else // !DEBUG_ATLAS
				if( texelType(I)!=TexelType::Unassigned ) MK_WARN( "Node ( " , i , " , " , j , " ) in chart " , chartID , " already covered [" , t , "]" );
#endif // DEBUG_ATLAS
				texelType(I) = TexelType::BoundarySupportedAndCovered;
#ifdef DEBUG_ATLAS
				nodeOwner(I) = atlasChart.triangles[t]();
#endif // DEBUG_ATLAS
				triangleID(I) = t;
				Point3D< double > bc = simplex.barycentricCoordinates( Point2D< double >( I[0] , I[1] ) );
				barycentricCoords(I) = Point2D< double >( bc[1] , bc[2] );
			};
		Rasterizer2D::RasterizeNodes< false >( simplex , Kernel , nodeRange );
		Rasterizer2D::RasterizeSupports< true , true >( simplex , [&]( Index I ){ cellType(I)=CellType::Interior; } , cellRange );
#else // !USE_RASTERIZER
		Point2D< GeometryReal > tPos[3];
		for( int i=0 ; i<3 ; i++ ) tPos[i] = atlasChart.vertices[ atlasChart.triangles[t][i] ] - gridChart.corner;
		int minCorner[2] , maxCorner[2];
		GetTriangleIntegerBBox( tPos , (GeometryReal)1./cellSizeW , (GeometryReal)1./cellSizeH , minCorner , maxCorner );

		SquareMatrix< GeometryReal , 2 > barycentricMap = GetBarycentricMap(tPos);

		for( int j=minCorner[1] ; j<=maxCorner[1] ; j++ ) for( int i=minCorner[0] ; i<=maxCorner[0] ; i++ )
		{
			Point2D< GeometryReal > texel_pos = Point2D< GeometryReal >( (GeometryReal)i*cellSizeW , (GeometryReal)j*cellSizeH ) - tPos[0];
			Point2D< GeometryReal > barycentricCoord = barycentricMap * texel_pos;
			if( barycentricCoord[0]>=0 && barycentricCoord[1]>=0 && ( barycentricCoord[0]+barycentricCoord[1] )<=1 )
			{
#ifdef DEBUG_ATLAS
				if( texelType(i,j)!=TexelType::Unsupported ) MK_WARN( "Texel ( " , i+gridChart.cornerCoords[0] , " , " , j+gridChart.cornerCoords[1] , " ) covered by multiple triangles: " , atlasChart.triangles[t]() , " " , nodeOwner(i,j) );
#else // !DEBUG_ATLAS
				if( texelType(i,j)!=TexelType::Unsupported ) MK_WARN( "Node ( " , i , " , " , j , " ) in chart " , chartID , " already covered [" , t , "]" );
#endif // DEBUG_ATLAS
				texelType(i,j) = TexelType::BoundarySupportedAndCovered;
#ifdef DEBUG_ATLAS
				nodeOwner(i,j) = atlasChart.triangles[t]();
#endif // DEBUG_ATLAS
				triangleID(i,j) = t;
				barycentricCoords(i,j) = barycentricCoord;
			}
		}
#endif // USE_RASTERIZER
	}

	//(2) Add texels adjacent to boundary cells
	for( int e=0 ; e<atlasChart.boundaryHalfEdges.size() ; e++ )
	{
		int tIndex = atlasChart.boundaryHalfEdges[e] / 3;
		int kIndex = atlasChart.boundaryHalfEdges[e] % 3;

#ifdef USE_RASTERIZER
		Simplex< double , 2 , 2 > simplex = GetTriangle( tIndex );
		Simplex< double , 2 , 1 > subSimplex;
		subSimplex[0] = simplex[kIndex] , subSimplex[1] = simplex[(kIndex+1)%3];
		auto Kernel = [&]( Index I )
			{
				if( texelType(I)!=TexelType::BoundarySupportedAndCovered )
				{
					texelType(I) = TexelType::BoundarySupported;
					Point3D< double > bc = simplex.barycentricCoordinates( Point2D< double >( I[0] , I[1] ) );
					Point2D< GeometryReal > newBC = Point2D< GeometryReal >( Point2D< double >( bc[1] , bc[2] ) );

					// If no triangle has been assigned, assign this one
					if( triangleID(I)==-1 )
					{
						triangleID(I) = tIndex;
						barycentricCoords(I) = newBC;
					}
					// Otherwise, update if closer
					else
					{
						Point2D< GeometryReal > oldBC = barycentricCoords(I);
						Simplex< double , 2 , 2 > oldSimplex = GetTriangle( triangleID(I) );
						Point< double , 2 > oldP = oldSimplex( oldBC ) , newP = simplex( newBC );
						double oldD2 = Point< double , 2 >::SquareDistance( oldP , oldSimplex.nearest(oldP) );
						double newD2 = Point< double , 2 >::SquareDistance( newP ,    simplex.nearest(newP) );
						if( newD2<oldD2 )
						{
							triangleID(I) = tIndex;
							barycentricCoords(I) = newBC;
						}
					}
				}
			};
		Rasterizer2D::RasterizeSupports< false , false >( subSimplex , Kernel , nodeRange );
		Rasterizer2D::RasterizeSupports< true , true >( subSimplex , [&]( Index I ){ cellType(I)=CellType::Boundary; } , cellRange );
#else // !USE_RASTERIZER
		Point2D< GeometryReal > ePos[2];
		SimplexIndex< 1 > eIndex = atlasChart.edgeIndex( atlasChart.boundaryHalfEdges[e] );
		ePos[0] = atlasChart.vertices[ eIndex[0] ] - gridChart.corner;
		ePos[1] = atlasChart.vertices[ eIndex[1] ] - gridChart.corner;

		int minCorner[2] , maxCorner[2];
		GetEdgeIntegerBBox( ePos , (GeometryReal)1./cellSizeW , (GeometryReal)1./cellSizeH , minCorner , maxCorner);
		Point2D< GeometryReal > tPos[3];
		for (int k = 0; k < 3; k++) tPos[k] = atlasChart.vertices[atlasChart.triangles[tIndex][k]] - gridChart.corner;

		SquareMatrix< GeometryReal , 2 > barycentricMap = GetBarycentricMap(tPos);

		Point2D< GeometryReal > edgeNormal;
		GeometryReal edgeLevel;
		Point2D< GeometryReal > edgeDirection = ePos[1] - ePos[0];
		edgeNormal = Point2D< GeometryReal >( edgeDirection[1] , -edgeDirection[0] );
		edgeNormal /= Point2D< GeometryReal >::Length( edgeNormal );
		edgeLevel = ( Point2D< GeometryReal >::Dot( edgeNormal , ePos[0] ) + Point2D< GeometryReal >::Dot( edgeNormal , ePos[1] ) ) / 2;

		//(2.1) Add texels adjacent to cell intersecting boundary edges

		for( int c=0 ; c<2 ; c++ ) for( int j=minCorner[1] ; j<=maxCorner[1] ; j++ ) for( int i=minCorner[0] ; i<=maxCorner[0] ; i++ )
		{
			Point2D< GeometryReal > cellNode[2] =
			{
				Point2D< GeometryReal >( (GeometryReal)i*cellSizeW , (GeometryReal)j*cellSizeH ),
				Point2D< GeometryReal >( (GeometryReal)i*cellSizeW , (GeometryReal)j*cellSizeH )
			};
			if( c==0 ) cellNode[1][c] += cellSizeW;
			else       cellNode[1][c] += cellSizeH;
			Point2D< GeometryReal > cellSide = cellNode[1] - cellNode[0];
			Point2D< GeometryReal > cellSideNormal = Point2D< GeometryReal >(cellSide[1], -cellSide[0]);
			cellSideNormal /= Point2D< GeometryReal >::Length(cellSideNormal);
			GeometryReal cellLevel = ( Point2D< GeometryReal >::Dot( cellSideNormal , cellNode[0] ) + Point2D< GeometryReal >::Dot( cellSideNormal , cellNode[1] ) ) / 2;

			// Are the nodes on opposite sides of the edge
			bool oppositeEdgeSide = ( Point2D< GeometryReal >::Dot( edgeNormal , cellNode[0] ) - edgeLevel ) * ( Point2D< GeometryReal >::Dot( edgeNormal , cellNode[1] ) - edgeLevel )<0;
			// Are the edge end-points on opposite sides of the node
			bool oppositeCellSide = ( Point2D< GeometryReal >::Dot( cellSideNormal , ePos[1] ) - cellLevel ) * ( Point2D< GeometryReal >::Dot( cellSideNormal , ePos[0] ) - cellLevel )<0;

			// Do the edge and the side of the node cross
			if( oppositeEdgeSide && oppositeCellSide )
			{
				if( c==0 ) cellType(i,j-1) = cellType(i,j) = CellType::Boundary;
				if( c==1 ) cellType(i-1,j) = cellType(i,j) = CellType::Boundary;


				for( int dn=-1 ; dn<=1 ; dn++ ) for( int dc=0 ; dc<2 ; dc++ ) //Update nodes on adjacent cells
				{
					int nIndices[2] = { i, j };
					nIndices[(1 - c)] += dn;
					nIndices[c] += dc;
					nIndices[0] = std::min< int >( std::max< int >( 0 , nIndices[0] ) , width  - 1 );
					nIndices[1] = std::min< int >( std::max< int >( 0 , nIndices[1] ) , height - 1 );
					if( texelType( nIndices[0] , nIndices[1] )!=TexelType::BoundarySupportedAndCovered )
					{
						texelType(nIndices[0],nIndices[1]) = TexelType::BoundarySupported;

						Point2D< GeometryReal > texel_pos = Point2D< GeometryReal >( (GeometryReal)nIndices[0]*cellSizeW , (GeometryReal)nIndices[1]*cellSizeH ) - tPos[0];
						Point2D< GeometryReal > barycentricCoord = barycentricMap*texel_pos;

						if( triangleID( nIndices[0] , nIndices[1] )==-1 )
						{
							triangleID(nIndices[0], nIndices[1]) = tIndex;
							barycentricCoords(nIndices[0], nIndices[1]) = barycentricCoord;
						}
						else //Update the position to the closest triangle
						{
							Point2D< GeometryReal > oldBarycentricCoord = barycentricCoords( nIndices[0] , nIndices[1] );
							Point3D< GeometryReal > oldBarycentricCoord3( (GeometryReal)1. - oldBarycentricCoord[0] - oldBarycentricCoord[1] , oldBarycentricCoord[0] , oldBarycentricCoord[1] );
							Point3D< GeometryReal > newBarycentricCoord3( (GeometryReal)1. -    barycentricCoord[0] -    barycentricCoord[1] ,    barycentricCoord[0] ,    barycentricCoord[1] );
							GeometryReal minOld = std::min< GeometryReal >( std::min< GeometryReal >( oldBarycentricCoord3[0] , oldBarycentricCoord3[1] ) , oldBarycentricCoord3[2] );
							GeometryReal minNew = std::min< GeometryReal >( std::min< GeometryReal >( newBarycentricCoord3[0] , newBarycentricCoord3[1] ) , newBarycentricCoord3[2] );
							if( minNew>minOld )
							{
								triangleID( nIndices[0] , nIndices[1] ) = tIndex;
								barycentricCoords( nIndices[0] , nIndices[1] ) = barycentricCoord;
							}
						}
					}
				}
			}
		}

		//(2.2) Add texels adjacent to cells that contain triangles
		if ((minCorner[0] + 1 == maxCorner[0]) && (minCorner[1] + 1 == maxCorner[1]))
		{
			cellType(minCorner[0], minCorner[1]) = CellType::Boundary;
			for (int di = 0; di < 2; di++) for (int dj = 0; dj < 2; dj++) {
				int nIndices[2] = { minCorner[0] + di, minCorner[1] + dj };
				nIndices[0] = std::min<int>(std::max<int>(0, nIndices[0]), width - 1);
				nIndices[1] = std::min<int>(std::max<int>(0, nIndices[1]), height - 1);
				if( texelType( nIndices[0], nIndices[1] )!=TexelType::BoundarySupportedAndCovered )
				{
					texelType( nIndices[0] , nIndices[1] ) = TexelType::BoundarySupported;

					Point2D< GeometryReal > texel_pos = Point2D< GeometryReal >( (GeometryReal)nIndices[0]*cellSizeW , (GeometryReal)nIndices[1]*cellSizeH ) - tPos[0];
					Point2D< GeometryReal > barycentricCoord = barycentricMap*texel_pos;

					if (triangleID(nIndices[0], nIndices[1]) == -1) {
						triangleID(nIndices[0], nIndices[1]) = tIndex;
						barycentricCoords(nIndices[0], nIndices[1]) = barycentricCoord;
					}
					else //Update the position to the closest triangle
					{
						Point2D< GeometryReal > oldBarycentricCoord = barycentricCoords( nIndices[0] , nIndices[1] );
						Point3D< GeometryReal > oldBarycentricCoord3( (GeometryReal)1. - oldBarycentricCoord[0] - oldBarycentricCoord[1] , oldBarycentricCoord[0] , oldBarycentricCoord[1] );
						Point3D< GeometryReal > newBarycentricCoord3( (GeometryReal)1. -    barycentricCoord[0] -    barycentricCoord[1] ,    barycentricCoord[0] ,    barycentricCoord[1] );
						GeometryReal minOld = std::min< GeometryReal >( std::min< GeometryReal >( oldBarycentricCoord3[0] , oldBarycentricCoord3[1] ) , oldBarycentricCoord3[2] );
						GeometryReal minNew = std::min< GeometryReal >( std::min< GeometryReal >( newBarycentricCoord3[0] , newBarycentricCoord3[1] ) , newBarycentricCoord3[2] );
						if( minNew>minOld )
						{
							triangleID( nIndices[0] , nIndices[1] ) = tIndex;
							barycentricCoords( nIndices[0] , nIndices[1] ) = barycentricCoord;
						}
					}

				}
			}
		}
#endif // USE_RASTERIZER
	}

#ifdef USE_RASTERIZER
#else // !USE_RASTERIZER
	//(3) Add interior cells
	for( int j=0 ; j<height-1 ; j++ ) for( int i=0 ; i<width-1 ; i++ )
		if( texelType(i,j)==TexelType::BoundarySupportedAndCovered && texelType(i+1,j)==TexelType::BoundarySupportedAndCovered && texelType(i,j+1)==TexelType::BoundarySupportedAndCovered && texelType(i+1,j+1)==TexelType::BoundarySupportedAndCovered && cellType(i,j)==CellType::Exterior ) cellType(i,j) = CellType::Interior;
#endif // USE_RASTERIZER

#ifdef USE_RASTERIZER
	//(3) Mark deep nodes
	{
		auto Kernel = [&]( Index I )
			{
				bool validDeepNode = true;
				Range::Intersect( cellRange , Range::CellsSupportedOnNode(I) ).process( [&]( Index I ){ validDeepNode &= ( cellType(I)==CellType::Interior ); } );
				if( validDeepNode ) texelType(I) = TexelType::InteriorSupported;
			};
		nodeRange.process( Kernel );
	}
#else // !USE_RASTERIZER
	for( int j=1 ; j<height-1 ; j++ ) for( int i=1 ; i<width-1 ; i++ )
	{
		bool validDeepNode = true;
		for( int di=-1 ; di<1 ; di++ ) for( int dj=-1 ; dj<1 ; dj++ ) if( cellType( i+di , j+dj )!=CellType::Interior ) validDeepNode = false;
		if( validDeepNode ) texelType(i,j) = TexelType::InteriorSupported;
	}
#endif // USE_RASTERIZER
#endif // USE_RASTERIZER

#ifdef NEW_INDEXING
	// Transform from chart triangle indices to atlas triangle indices
	{
		Image< AtlasMeshTriangleIndex > &atlasTriangleID = gridChart.triangleID;
		atlasTriangleID.resize( width , height );
		for( unsigned int i=0 ; i<atlasTriangleID.size() ; i++ ) atlasTriangleID[i] = AtlasMeshTriangleIndex( -1 );

		for( size_t i=0 ; i<chartTriangleID.size() ; i++ ) if( chartTriangleID[i]!=ChartMeshTriangleIndex(-1) ) atlasTriangleID[i] = atlasChart.atlasTriangle( chartTriangleID[i] );
	}
#else // !NEW_INDEXING
	// Make the triangle IDs chart-local
	for( size_t i=0 ; i<triangleID.size() ; i++ ) if( triangleID[i]!=-1 ) triangleID[i] = atlasChart.atlasTriangle( triangleID[i] );
#endif // NEW_INDEXING

	// (5) Enumerate variables in raster order
	gridChart.combinedCellOffset = combinedCellIndex;
	gridChart.interiorCellOffset = interiorCellIndex;

	Image< TexelIndex > & texelIndices = gridChart.texelIndices;
	texelIndices.resize( width , height );

	// Set the texel indices
	for( unsigned int j=0 ; j<height ; j++ ) for( unsigned int i=0 ; i<width ; i++ ) if( gridChart.texelType(i,j)!=TexelType::Unsupported )
	{
		texelIndices(i,j).combined = combinedTexelIndex++;		

		GridNodeInfo currentNodeInfo;
		currentNodeInfo.ci = i;
		currentNodeInfo.cj = j;
		currentNodeInfo.chartID = chartID;
		currentNodeInfo.texelType = gridChart.texelType(i,j);
		nodeInfo.push_back( currentNodeInfo );
		if( IsCovered( gridChart.texelType(i,j) ) ) texelIndices(i,j).interiorOrCovered = interiorTexelIndex++;
		if( IsBoundarySupported( gridChart.texelType(i,j) ) ) boundaryTexelIndex++;
		if( gridChart.texelType(i,j)==TexelType::InteriorSupported ) texelIndices(i,j).interior = deepTexelIndex++;
	}

#ifdef USE_RASTERIZER
	// Sanity check: Confirm that all active/interior cells are incident on active/non-boundary nodes
	{
		auto Kernel = [&]( Index I )
			{
				if( gridChart.cellType(I)!=CellType::Exterior )
				{
					unsigned int count = 0;
					auto SubKernel = [&]( Index I ){ if( texelIndices(I).combined!=-1 ){ count++; } };
					Range::NodesSupportedOnCell( I ).process( SubKernel );
					if( count!=4 ) MK_THROW( "Active cell adjacent to inactive node" );

				}
				if( gridChart.cellType(I)==CellType::Interior )
				{
					unsigned int count = 0;
					auto SubKernel = [&]( Index I ){ if( texelIndices(I).interiorOrCovered!=-1 ){ count++; } };
					Range::NodesSupportedOnCell( I ).process( SubKernel );
					if( count!=4 ) MK_THROW( "Interior cell adjacent to non interior node" );
				}
			};
		cellRange.process( Kernel );
	}
#endif // USE_RASTERIZER


	Image< CellIndex > & cellIndices = gridChart.cellIndices;
	cellIndices.resize( width-1 , height-1 );

	std::vector< BilinearElementIndex > & combinedCellCombinedBilinearElementIndices = gridChart.combinedCellCombinedBilinearElementIndices;
	std::vector< BilinearElementIndex > & interiorCellInteriorBilinearElementIndices = gridChart.interiorCellInteriorBilinearElementIndices;
	std::vector< BilinearElementIndex > & interiorCellCombinedBilinearElementIndices = gridChart.interiorCellCombinedBilinearElementIndices;

	std::vector< unsigned int > & interiorCellIndexToCombinedCellIndex = gridChart.interiorCellIndexToCombinedCellIndex;
	std::vector< unsigned int > & boundaryCellIndexToCombinedCellIndex = gridChart.boundaryCellIndexToCombinedCellIndex;

	unsigned int _combinedCellIndex = 0;
	unsigned int _boundaryCellIndex = 0;
	unsigned int _interiorCellIndex = 0;

	// For all cells supporting the triangle:
	// -- set the (combined) indices of the corners
	// -- if the cell insterior
	// ---- set the (interior/covered) indices of the corners
	for( unsigned int j=0 ; j<height-1 ; j++ ) for( unsigned int i=0 ; i<width-1 ; i++ ) if( gridChart.cellType(i,j)!=CellType::Exterior )
	{
#ifdef USE_RASTERIZER
		combinedCellCombinedBilinearElementIndices.emplace_back( texelIndices(i,j).combined , texelIndices(i+1,j).combined , texelIndices(i+1,j+1).combined , texelIndices(i,j+1).combined );

		if( gridChart.cellType(i,j)==CellType::Boundary )
		{
			cellIndices(i,j).boundary = _boundaryCellIndex++;
			boundaryCellIndexToCombinedCellIndex.push_back( _combinedCellIndex );
		}
		else
		{
			cellIndices(i,j).interior = _interiorCellIndex++;
			interiorCellIndexToCombinedCellIndex.push_back( _combinedCellIndex );

			interiorCellInteriorBilinearElementIndices.emplace_back( texelIndices(i,j).interiorOrCovered , texelIndices(i+1,j).interiorOrCovered , texelIndices(i+1,j+1).interiorOrCovered , texelIndices(i,j+1).interiorOrCovered );
			interiorCellCombinedBilinearElementIndices.emplace_back( texelIndices(i,j).combined , texelIndices(i+1,j).combined , texelIndices(i+1,j+1).combined , texelIndices(i,j+1).combined );
		}
		cellIndices(i,j).combined = _combinedCellIndex++;
#else // !USE_RASTERIZER
		int globalTexelIndices[4] = { texelIndices(i,j).combined , texelIndices(i+1,j).combined , texelIndices(i+1,j+1).combined , texelIndices(i,j+1).combined };
		if( globalTexelIndices[0]!=-1 && globalTexelIndices[1]!=-1 && globalTexelIndices[2]!=-1 && globalTexelIndices[3] != -1 )
			bilinearElementIndices.push_back( BilinearElementIndex( globalTexelIndices[0] , globalTexelIndices[1] , globalTexelIndices[2] , globalTexelIndices[3] ) );
		else MK_THROW( "Active cell adjacent to unactive node" );

		if( gridChart.cellType(i,j)==CellType::Boundary )
		{
			cellIndices(i,j).boundary = _boundaryCellIndex;
			boundaryCellIndexToCombinedCellIndex.push_back( _combinedCellIndex );
			_boundaryCellIndex++;
		}
		else
		{
			cellIndices(i,j).interior = _interiorCellIndex;
			interiorCellIndexToCombinedCellIndex.push_back( _combinedCellIndex );

			int globalTexelInteriorIndices[4] = { texelIndices(i,j).interiorOrCovered , texelIndices(i+1,j).interiorOrCovered , texelIndices(i+1,j+1).interiorOrCovered , texelIndices(i,j+1).interiorOrCovered };
			if( globalTexelInteriorIndices[0]!=-1 && globalTexelInteriorIndices[1]!=-1 && globalTexelInteriorIndices[2]!=-1 && globalTexelInteriorIndices[3]!=-1)
			{
				interiorCellInteriorBilinearElementIndices.push_back( BilinearElementIndex( globalTexelInteriorIndices[0] , globalTexelInteriorIndices[1] , globalTexelInteriorIndices[2] , globalTexelInteriorIndices[3] ) );
				interiorCellCombinedBilinearElementIndices.push_back( BilinearElementIndex( globalTexelIndices[0] , globalTexelIndices[1] , globalTexelIndices[2] , globalTexelIndices[3] ) );
			}
			else MK_THROW( "Interior cell adjacent to non interior node" );
			_interiorCellIndex++;
		}
		cellIndices(i,j).combined = _combinedCellIndex++;
#endif // USE_RASTERIZER
	}

	combinedCellIndex += _combinedCellIndex;
	boundaryCellIndex += _boundaryCellIndex;
	interiorCellIndex += _interiorCellIndex;

	// (6) Construct raster lines
	for( unsigned int j=0 ; j<height ; j++ )
	{
		bool firstSegment = true;
		bool previousDeep = false;
		int rasterStart = -1;
		for( unsigned int i=0 ; i<width ; i++ )
		{
			bool currentIsDeep = texelType(i,j)==TexelType::InteriorSupported;
			if( currentIsDeep && !previousDeep )
			{
				rasterStart = i; // Start of raster line
				if( firstSegment )
				{
					firstSegment = false;
					SegmentedRasterLine newSegmentLine;
					segmentedLines.push_back( newSegmentLine );
				}
			}
			if( !currentIsDeep && previousDeep ) // End of raster line
			{ 
				RasterLine newLine;
				newLine.lineStartIndex  = texelIndices( rasterStart , j   ).combined;
				newLine.lineEndIndex    = texelIndices(         i-1 , j   ).combined;
				newLine.prevLineIndex   = texelIndices( rasterStart , j-1 ).combined;
				newLine.nextLineIndex   = texelIndices( rasterStart , j+1 ).combined;
				newLine.coeffStartIndex = texelIndices( rasterStart , j   ).interior;

				if( newLine.lineStartIndex==-1 || newLine.lineEndIndex==-1 || newLine.prevLineIndex==-1 || newLine.nextLineIndex==-1 ) MK_THROW( "Inavlid Indexing" );
				rasterLines.push_back( newLine );

				SegmentedRasterLine & newSegmentLine = segmentedLines.back();
				newSegmentLine.segments.push_back( newLine );
			}
			previousDeep = currentIsDeep;
		}
	}

	//Initialize thread tasks
	int blockHorizontalOffset = multigridBlockInfo.blockWidth - multigridBlockInfo.paddingWidth;
	int blockVerticalOffset = multigridBlockInfo.blockHeight - multigridBlockInfo.paddingHeight;
	int numHorizontalBlocks = ((width - multigridBlockInfo.paddingWidth - 1) / blockHorizontalOffset) + 1;
	int numVerticalBlocks = ((height - multigridBlockInfo.paddingHeight - 1) / blockVerticalOffset) + 1;

	for( unsigned int bj=0 ; bj<numVerticalBlocks ; bj++ )
	{
		ThreadTask threadTask;
		int taskDeepTexels = 0;

		int blockVerticalStart = bj*blockVerticalOffset;
		int blockVerticalEnd = std::min<int>((bj + 1)*blockVerticalOffset + multigridBlockInfo.paddingHeight + 2 - 1, height - 1);

		for( unsigned int bi=0 ; bi<numHorizontalBlocks ; bi++ )
		{
			BlockTask blockTask;

			int blockHorizontalStart = bi*blockHorizontalOffset;
			int blockHorizontalEnd = std::min<int>((bi + 1)*blockHorizontalOffset + multigridBlockInfo.paddingWidth + 2 - 1, width - 1);

			//Deep texel within rows[blockVerticalStart + 1,blockVerticalEnd - 1]column [blockHorizontalStart + 1,blockHorizontalEnd - 1] 
			std::vector<BlockDeepSegmentedLine>  &  blockDeepSegmentedLines = blockTask.blockDeepSegmentedLines;
			for (int j = blockVerticalStart + 1; j <= blockVerticalEnd - 1; j++) {
				BlockDeepSegmentedLine segmentedLine;
				int offset = blockHorizontalStart + 1;
				int segmentStart = -1;
				bool previousDeep = false;
				while (offset <= blockHorizontalEnd - 1)
				{
					bool currentDeep = (texelIndices(offset, j).interior != -1);
					if (currentDeep && !previousDeep) {//Start segment
						segmentStart = offset;
					}

					if ((previousDeep && !currentDeep) || (currentDeep && offset == blockHorizontalEnd - 1)) {//Terminate segment
						BlockDeepSegment deepSegment;
						deepSegment.currentStart               = texelIndices(segmentStart, j).combined;
						if (currentDeep)deepSegment.currentEnd = texelIndices(offset, j).combined;
						else deepSegment.currentEnd            = texelIndices(offset - 1, j).combined;
						deepSegment.previousStart              = texelIndices(segmentStart, j - 1).combined;
						deepSegment.nextStart                  = texelIndices(segmentStart, j + 1).combined;

						deepSegment.deepStart = texelIndices(segmentStart, j).interior;

						segmentedLine.blockDeepSegments.push_back(deepSegment);
					}
					previousDeep = currentDeep;
					offset++;
				}
				if (segmentedLine.blockDeepSegments.size()) {
					blockDeepSegmentedLines.push_back(segmentedLine);
				}
			}

			if (blockDeepSegmentedLines.size()) {

				//Deep texel within rows[blockVerticalStart,blockVerticalEnd]column [blockHorizontalStart,blockHorizontalEnd] 
				std::vector<BlockDeepSegmentedLine>  &  blockPaddedSegmentedLines = blockTask.blockPaddedSegmentedLines;
				for (int j = blockVerticalStart; j <= blockVerticalEnd; j++) {
					BlockDeepSegmentedLine segmentedLine;
					int offset = blockHorizontalStart;
					int segmentStart = -1;
					bool previousDeep = false;
					while (offset <= blockHorizontalEnd)
					{
						bool currentDeep = (texelIndices(offset, j).interior != -1);
						if (currentDeep && !previousDeep) {//Start segment
							segmentStart = offset;
						}

						if ((previousDeep && !currentDeep) || (currentDeep && offset == blockHorizontalEnd)) {//Terminate segment
							BlockDeepSegment deepSegment;
							deepSegment.currentStart               = texelIndices(segmentStart, j).combined;
							if (currentDeep)deepSegment.currentEnd = texelIndices(offset, j).combined;
							else deepSegment.currentEnd            = texelIndices(offset - 1, j).combined;
							deepSegment.previousStart              = texelIndices(segmentStart, j - 1).combined;
							deepSegment.nextStart                  = texelIndices(segmentStart, j + 1).combined;

							deepSegment.deepStart = texelIndices(segmentStart, j).interior;

							segmentedLine.blockDeepSegments.push_back(deepSegment);
						}
						previousDeep = currentDeep;
						offset++;
					}
					if (segmentedLine.blockDeepSegments.size()) {
						blockPaddedSegmentedLines.push_back(segmentedLine);
					}
				}
				threadTask.blockTasks.push_back(blockTask);
			}
		}
		if( threadTask.blockTasks.size() )
		{
			threadTask.taskDeepTexels = taskDeepTexels;
			threadTasks.push_back( threadTask );
		}
	}
}

template< typename GeometryReal >
void InitializeGridCharts
(
	const std::vector< AtlasChart< GeometryReal > > &atlasCharts ,
	unsigned int width , 
	unsigned int height ,
	unsigned int level ,
	std::vector< GridNodeInfo > &nodeInfo ,
	std::vector< GridChart< GeometryReal > > &gridCharts ,
	std::vector< RasterLine > &rasterLines ,
	std::vector< SegmentedRasterLine > &segmentedLines ,
	std::vector< ThreadTask > &threadTasks ,
	unsigned int &numTexels ,
	unsigned int &numInteriorTexels ,
	unsigned int &numDeepTexels ,
	unsigned int &numBoundaryTexels ,
	unsigned int &numCells ,
	unsigned int &numBoundaryCells ,
	unsigned int &numInteriorCells ,
	const MultigridBlockInfo &multigridBlockInfo
)
{
	gridCharts.resize( atlasCharts.size() );

	numTexels = 0;
	numInteriorTexels = 0;
	numDeepTexels = 0;
	numBoundaryTexels = 0;
	numCells = 0;
	numBoundaryCells = 0;
	numInteriorCells = 0;

	GeometryReal cellSize[] = { (GeometryReal)(1<<level) / width , (GeometryReal)(1<<level) / height };

	unsigned int _width  = (  width + ( (1<<level)-1 ) )>>level;
	unsigned int _height = ( height + ( (1<<level)-1 ) )>>level;
	unsigned int approxGridSize = _width * _height;

	for( unsigned int i=0 ; i<atlasCharts.size() ; i++ )
	{
		{
			int halfSize[2][2];
			for( unsigned int c=0 ; c<2 ; c++ )
			{
				halfSize[c][0] = (int)ceil( ( atlasCharts[i].gridOrigin[c] - atlasCharts[i].minCorner[c] ) / cellSize[c] );
				halfSize[c][1] = (int)ceil( ( atlasCharts[i].maxCorner[c] - atlasCharts[i].gridOrigin[c] ) / cellSize[c] );
				gridCharts[i].corner[c] = atlasCharts[i].gridOrigin[c] - cellSize[c] * (GeometryReal)halfSize[c][0];
				gridCharts[i].centerOffset[c] = halfSize[c][0];
				gridCharts[i].cornerCoords[c] = atlasCharts[i].originCoords[c] - halfSize[c][0];
			}

			gridCharts[i].width  = halfSize[0][0] + halfSize[0][1] + 1;
			gridCharts[i].height = halfSize[1][0] + halfSize[1][1] + 1;
		}
		gridCharts[i].cellSizeW = cellSize[0];
		gridCharts[i].cellSizeH = cellSize[1];
		gridCharts[i].atlasWidth = _width;
		gridCharts[i].atlasHeight = _height;
		InitializeGridChartsActiveNodes( i , atlasCharts[i] , gridCharts[i] , nodeInfo , rasterLines , segmentedLines , threadTasks , numTexels , numInteriorTexels , numDeepTexels , numBoundaryTexels , numCells , numBoundaryCells , numInteriorCells , multigridBlockInfo );
	}
}

template< typename GeometryReal >
void InitializeTextureNodes( const std::vector< GridChart< GeometryReal > > &gridCharts , std::vector< TextureNodeInfo< GeometryReal > > &textureNodes )
{
	for( unsigned int c=0 ; c<gridCharts.size() ; c++ )
	{
		const GridChart< GeometryReal > &gridChart = gridCharts[c];
		const Image< TexelIndex > & texelIndices = gridChart.texelIndices;
		for( unsigned int j=0 ; j<gridChart.height ; j++ ) for( unsigned int i=0 ; i<gridChart.width ; i++ ) if( texelIndices(i,j).combined!=-1 )
		{
			TextureNodeInfo< GeometryReal > textureNode;
			textureNode.ci = gridChart.cornerCoords[0] + i;
			textureNode.cj = gridChart.cornerCoords[1] + j;
			textureNode.chartID = c;
			textureNode.tID = gridChart.triangleID(i,j);
			textureNode.barycentricCoords = gridChart.barycentricCoords(i, j);
			textureNode.isInterior = IsCovered( gridChart.texelType(i,j) );
			textureNodes.push_back( textureNode );
		}
	}
}

template< typename GeometryReal >
void InitializeCellNodes( const std::vector< GridChart< GeometryReal > > &gridCharts , std::vector< BilinearElementIndex > &combinedCellCombinedBilinearElementIndices )
{
	for( int c=0 ; c<gridCharts.size() ; c++ )
	{
		const std::vector< BilinearElementIndex >& _combinedCellCombinedBilinearElementIndices = gridCharts[c].combinedCellCombinedBilinearElementIndices;
		combinedCellCombinedBilinearElementIndices.insert( combinedCellCombinedBilinearElementIndices.end() , _combinedCellCombinedBilinearElementIndices.begin() , _combinedCellCombinedBilinearElementIndices.end() );
	}
}

template< typename GeometryReal , typename MatrixReal >
void InitializeAtlasHierachicalBoundaryCoefficients
(
	const GridAtlas< GeometryReal , MatrixReal > &fineAtlas ,
	GridAtlas< GeometryReal , MatrixReal > &coarseAtlas ,
	SparseMatrix< MatrixReal , int > &boundaryCoarseFineProlongation ,
	std::vector< BoundaryDeepIndex > &boundaryDeepIndices ,
	std::vector< BoundaryBoundaryIndex< MatrixReal > > &boundaryBoundaryIndices
)
{
	std::vector< Eigen::Triplet< MatrixReal > > prolongationTriplets;

	const std::vector<GridNodeInfo> & coarseNodeInfo = coarseAtlas.nodeInfo;
	const typename GridAtlas<>::IndexConverter & coarseIndexConverter = coarseAtlas.indexConverter;
	const typename GridAtlas<>::IndexConverter & fineIndexConverter = fineAtlas.indexConverter;
	for( unsigned int i=0 ; i<coarseIndexConverter.numBoundary() ; i++ )
	{
		const GridNodeInfo & currentNode = coarseNodeInfo[ coarseIndexConverter.boundaryToSupported(i) ];
		const GridChart< GeometryReal > &coarseChart = coarseAtlas.gridCharts[ currentNode.chartID ];
		const GridChart< GeometryReal > &fineChart = fineAtlas.gridCharts[ currentNode.chartID ];

		unsigned int coarseChartWidth = coarseChart.texelType.res(0);
		unsigned int coarseChartHeight = coarseChart.texelType.res(1);

		unsigned int fineChartWidth = fineChart.texelType.res(0);
		unsigned int fineChartHeight = fineChart.texelType.res(1);

		int ci = currentNode.ci;
		int cj = currentNode.cj;

		int numBoundaryNeighbours = 0;// -1: not defined, 0: boundary, 1: interior 
		int boundary_id[9];
		int boundary_offset_i[9];
		int boundary_offset_j[9];
		for( int li=-1 ; li<=1 ; li++ ) for( int lj=-1 ; lj<=1 ; lj++ )
		{
			int pi = ci + li;
			int pj = cj + lj;
			if( pi>=0 && pi<coarseChartWidth && pj>=0 && pj<coarseChartHeight )
			{
				unsigned int coarseGlobalIndex = coarseChart.texelIndices(pi, pj).combined;
				if( coarseGlobalIndex!=-1 )
				{
					unsigned int coarseBoundaryIndex = coarseIndexConverter.supportedToBoundary( coarseGlobalIndex );
					unsigned int coarseDeepIndex = coarseIndexConverter.supportedToDeep( coarseGlobalIndex );
					if( coarseDeepIndex!=-1 )
					{
						//Deep
						BoundaryDeepIndex bdIndex;
						bdIndex.boundaryIndex = i;
						bdIndex.deepGlobalIndex = coarseGlobalIndex;
						bdIndex.deepIndex = coarseDeepIndex;
						bdIndex.offset = (1 - lj) * 3 + (1 - li);
						boundaryDeepIndices.push_back(bdIndex);
					}
					else if( coarseBoundaryIndex!=-1 )
					{
						//Boundary
						boundary_id[ numBoundaryNeighbours ] = coarseBoundaryIndex;
						boundary_offset_i[numBoundaryNeighbours] = 2 * li;
						boundary_offset_j[numBoundaryNeighbours] = 2 * lj;
						numBoundaryNeighbours++;
					}
					else MK_ERROR_OUT( "Expected a supported index" );
				}
			}
		}

		int ri = ci - coarseChart.centerOffset[0];
		int rj = cj - coarseChart.centerOffset[1];

		int fi = fineChart.centerOffset[0] + 2 * ri;
		int fj = fineChart.centerOffset[1] + 2 * rj;


		for( int di=-1 ; di<=1 ; di++ ) for( int dj=-1 ; dj<=1 ; dj++ )
		{
			int pi = fi + di;
			int pj = fj + dj;
			if (pi >= 0 && pi < fineChartWidth && pj >= 0 && pj < fineChartHeight) {
				int fineGlobalIndex = fineChart.texelIndices(pi, pj).combined;
				if( fineGlobalIndex!=-1 )
				{
					unsigned int fineBoundaryIndex = fineIndexConverter.supportedToBoundary( fineGlobalIndex );
					unsigned int fineDeepIndex = fineIndexConverter.supportedToDeep( fineGlobalIndex );
					MatrixReal weight = (MatrixReal)( 1.0 / (1.0 + std::abs(di) ) / ( 1.0 + std::abs(dj) ) );
					if( fineBoundaryIndex!=-1 )
					{
						// Boundary
						prolongationTriplets.emplace_back( fineBoundaryIndex , i , weight );

						for( int ki=-1; ki<=1 ; ki++ ) for( int kj=-1 ; kj<=1 ; kj++ )
						{
							int qi = pi + ki;
							int qj = pj + kj;
							if (qi >= 0 && qi < fineChartWidth && qj >= 0 && qj < fineChartHeight) {
								int neighboutFineGlobalIndex = fineChart.texelIndices(qi, qj).combined;
								if( neighboutFineGlobalIndex!=-1 )
								{
									unsigned int _fineDeepIndex = fineIndexConverter.supportedToDeep( neighboutFineGlobalIndex );
									if( _fineDeepIndex!=-1 )
									{
										//Deep
										int oi = di + ki;
										int oj = dj + kj;
										for (int n = 0; n < numBoundaryNeighbours; n++) {
											MatrixReal diff_i = (MatrixReal)std::abs( oi - boundary_offset_i[n] );
											MatrixReal diff_j = (MatrixReal)std::abs( oj - boundary_offset_j[n] );
											if (diff_i < 1.5 && diff_j < 1.5) {
												MatrixReal weight2 = (MatrixReal)( 1.0 / (1.0 + diff_i ) / ( 1.0 + diff_j ) );
												BoundaryBoundaryIndex< MatrixReal > bbIndex;
												bbIndex.coarsePrincipalBoundaryIndex = i;
												bbIndex.coarseSecondaryBoundaryIndex = boundary_id[n];
												bbIndex.fineDeepIndex = _fineDeepIndex;
												bbIndex.offset = (1 - kj) * 3 + (1 - ki);
												bbIndex.weight = weight * weight2;
												boundaryBoundaryIndices.push_back(bbIndex);
											}
										}
									}
								}
							}
						}
					}
					else if( fineDeepIndex!=-1 )
					{
						// Deep
						for( int ki=-1 ; ki<=1 ; ki++ ) for( int kj=-1 ; kj<=1 ; kj++ )
						{
							int oi = di + ki;
							int oj = dj + kj;

							for (int n = 0; n < numBoundaryNeighbours; n++)
							{
								MatrixReal diff_i = (MatrixReal)std::abs( oi - boundary_offset_i[n] );
								MatrixReal diff_j = (MatrixReal)std::abs( oj - boundary_offset_j[n] );
								if (diff_i < 1.5 && diff_j < 1.5)
								{
									MatrixReal weight2 = (MatrixReal)( 1.0 / ( 1.0 + diff_i ) / ( 1.0 + diff_j ) );
									BoundaryBoundaryIndex< MatrixReal > bbIndex;
									bbIndex.coarsePrincipalBoundaryIndex = i;
									bbIndex.coarseSecondaryBoundaryIndex = boundary_id[n];
									bbIndex.fineDeepIndex = fineDeepIndex;
									bbIndex.offset = (1 + kj) * 3 + (1 + ki);
									bbIndex.weight = weight * weight2;
									boundaryBoundaryIndices.push_back(bbIndex);
								}
							}
						}
					}
					else MK_ERROR_OUT( "Expected supported index" );
				}
			}
		}
	}

	int numCoarseBoundaryTexels = coarseAtlas.numBoundaryTexels;
	int numFineBoundaryTexels = fineAtlas.numBoundaryTexels;
	boundaryCoarseFineProlongation = SetSparseMatrix( prolongationTriplets , numFineBoundaryTexels , numCoarseBoundaryTexels , false );
}


template< typename GeometryReal , typename MatrixReal >
void InitializeHierarchy
(
	unsigned int width ,
	unsigned int height ,
	HierarchicalSystem< GeometryReal , MatrixReal > &hierarchy ,
	std::vector< AtlasChart< GeometryReal > > &atlasCharts ,
	unsigned int levels ,
	const MultigridBlockInfo &multigridBlockInfo
)
{
	std::vector< GridAtlas< GeometryReal , MatrixReal > > &gridAtlases = hierarchy.gridAtlases;
	gridAtlases.resize( levels );

	for( unsigned int l=0 ; l<levels ; l++ )
	{
		InitializeGridCharts( atlasCharts , width , height , l , gridAtlases[l].nodeInfo , gridAtlases[l].gridCharts , gridAtlases[l].rasterLines , gridAtlases[l].segmentedLines , gridAtlases[l].threadTasks , gridAtlases[l].numTexels , gridAtlases[l].numInteriorTexels , gridAtlases[l].numDeepTexels , gridAtlases[l].numBoundaryTexels , gridAtlases[l].numCells , gridAtlases[l].numBoundaryCells , gridAtlases[l].numInteriorCells , multigridBlockInfo );

		if( gridAtlases[l].numTexels!=gridAtlases[l].numBoundaryTexels+gridAtlases[l].numDeepTexels )
			MK_THROW( "Boundary and deep texels does not form a partition: " , gridAtlases[l].numTexels , " != " , gridAtlases[l].numBoundaryTexels , " + " , gridAtlases[l].numDeepTexels );

		InitializeIndexConverter( gridAtlases[l].gridCharts , gridAtlases[l].numTexels , gridAtlases[l].indexConverter );
	}

	hierarchy.boundaryRestriction.resize( levels-1 );
	for( unsigned int i=0 ; i<levels-1 ; i++ ) InitializeAtlasHierachicalRestriction( gridAtlases[i] , gridAtlases[i+1] , hierarchy.boundaryRestriction[i] );
	for( unsigned int i=0 ; i<levels-1 ; i++ ) InitializeAtlasHierachicalProlongation( gridAtlases[i] , gridAtlases[i+1] );


	hierarchy.boundaryCoarseFineProlongation.resize(levels);
	hierarchy.boundaryFineCoarseRestriction.resize(levels);
	hierarchy.boundaryDeepIndices.resize(levels);
	hierarchy.boundaryBoundaryIndices.resize(levels);

	for( unsigned int i=1 ; i<levels ; i++ )
	{
		InitializeAtlasHierachicalBoundaryCoefficients( hierarchy.gridAtlases[i-1] , hierarchy.gridAtlases[i] , hierarchy.boundaryCoarseFineProlongation[i] , hierarchy.boundaryDeepIndices[i] , hierarchy.boundaryBoundaryIndices[i] );
		hierarchy.boundaryFineCoarseRestriction[i-1] = hierarchy.boundaryCoarseFineProlongation[i].transpose();
	}
}

template< typename GeometryReal , typename MatrixReal >
void InitializeHierarchy
(
	const TexturedTriangleMesh< GeometryReal > &mesh ,
	unsigned int width ,
	unsigned int height ,
	unsigned int levels ,
	std::vector< TextureNodeInfo< GeometryReal > > &textureNodes ,
	std::vector< BilinearElementIndex > &cellNodes ,
	HierarchicalSystem< GeometryReal , MatrixReal > &hierarchy ,
	std::vector< AtlasChart< GeometryReal > > &atlasCharts ,
	const MultigridBlockInfo &multigridBlockInfo ,
	bool computeProlongation
)
{
	typename AtlasChart< GeometryReal >::AtlasInfo atlasInfo;

	//(1) Initialize atlas charts
	atlasCharts = AtlasChart< GeometryReal >::GetCharts( mesh , width , height , atlasInfo );

	//(2) Initialize hierarchy
	InitializeHierarchy( width , height , hierarchy , atlasCharts , levels , multigridBlockInfo );

	//(3) Initialize fine level texture nodes and cells
	InitializeTextureNodes( hierarchy.gridAtlases[0].gridCharts , textureNodes );
	InitializeCellNodes( hierarchy.gridAtlases[0].gridCharts , cellNodes );

	//(4) Initialize boundary triangulation
	InitializeBoundaryTriangulation( hierarchy.gridAtlases[0] , atlasCharts , atlasInfo );

	//(5) Unnnecesary
	if( computeProlongation )
	{
		Eigen::SparseMatrix< MatrixReal > __prolongation;
		InitializeProlongation( hierarchy.gridAtlases[0].numInteriorTexels , hierarchy.gridAtlases[0].numFineNodes , hierarchy.gridAtlases[0].numTexels , hierarchy.gridAtlases[0].gridCharts , hierarchy.gridAtlases[0].nodeInfo , __prolongation );
		hierarchy.gridAtlases[0].coarseToFineNodeProlongation = __prolongation;
	}
}