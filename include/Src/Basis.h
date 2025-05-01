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

#include <Misha/Miscellany.h>
#include <Misha/Exceptions.h>
#include <Misha/RightTriangleQuadrature.h>
#include "SimpleTriangleMesh.h"
#include "Indices.h"

namespace MishaK
{
	template< typename GeometryReal >
	struct TextureNodeInfo
	{
		TextureNodeInfo( void )
			: tID(-1) , ci( static_cast< unsigned int >(-1) ) , cj( static_cast< unsigned int >(-1) ) , chartID( static_cast< unsigned int >(-1) ) , isInterior(false) {}
		TextureNodeInfo( AtlasMeshTriangleIndex tID , Point2D< GeometryReal > barycentricCoords , unsigned int ci , unsigned int cj , ChartIndex chartID , bool isInterior )
			: tID(tID) , barycentricCoords(barycentricCoords) , ci(ci) , cj(cj) , chartID(chartID) , isInterior(isInterior) {}
		AtlasMeshTriangleIndex tID;
		Point< GeometryReal , 2 > barycentricCoords;
		unsigned int ci , cj;
		ChartIndex chartID;
		bool isInterior;

		operator MeshSample< GeometryReal >() const { return MeshSample< GeometryReal >( static_cast< unsigned int >(tID) , barycentricCoords ); }
	};

	struct BilinearElementIndex
	{
	protected:
		unsigned int v[4];
	public:
		BilinearElementIndex( void ) { v[0] = v[1] = v[2] = v[3] = 0; }
		BilinearElementIndex( unsigned int v0 , unsigned int v1 , unsigned int v2 , unsigned int v3 ) { v[0] = v0 , v[1] = v1 , v[2] = v2 , v[3] = v3; }
		unsigned int &operator[]( unsigned int idx )       { return v[idx]; }
		unsigned int  operator[]( unsigned int idx ) const { return v[idx]; }
	};

	struct QuadraticElementIndex
	{
	protected:
		AtlasInteriorOrBoundaryNodeIndex v[6];
	public:
		QuadraticElementIndex( void ) { v[0] = v[1] = v[2] = v[3] = v[4] = v[5] = AtlasInteriorOrBoundaryNodeIndex(0); }
		QuadraticElementIndex( AtlasInteriorOrBoundaryNodeIndex v0 , AtlasInteriorOrBoundaryNodeIndex v1 , AtlasInteriorOrBoundaryNodeIndex v2 , AtlasInteriorOrBoundaryNodeIndex v3 , AtlasInteriorOrBoundaryNodeIndex v4 , AtlasInteriorOrBoundaryNodeIndex v5 ){ v[0] = v0 , v[1] = v1 , v[2] = v2 , v[3] = v3 , v[4] = v4 , v[5] = v5; }
		AtlasInteriorOrBoundaryNodeIndex &operator[]( unsigned int idx )       { return v[idx]; }
		AtlasInteriorOrBoundaryNodeIndex  operator[]( unsigned int idx ) const { return v[idx]; }
	};

	template< typename GeometryReal >
	GeometryReal QuadraticElementValue( int elementIndex , Point2D< GeometryReal > pos )
	{
		switch( elementIndex )
		{
		case 0:	return 2 * pos[0] * pos[0] + 4 * pos[0] * pos[1] + 2 * pos[1] * pos[1] - 3 * pos[0] - 3 * pos[1] + 1;
		case 1: return 2 * pos[0] * pos[0] - 1 * pos[0];
		case 2: return 2 * pos[1] * pos[1] - 1 * pos[1];
		case 3: return 4 * pos[0] * pos[1];
		case 4: return -4 * pos[0] * pos[1] - 4 * pos[1] * pos[1] + 4 * pos[1];
		case 5: return -4 * pos[0] * pos[0] - 4 * pos[0] * pos[1] + 4 * pos[0];
		default: MK_THROW( "Element out of bounds" );
		}
		return (GeometryReal)0;
	}

	template< typename GeometryReal >
	Point2D< GeometryReal > QuadraticElementGradient( int elementIndex , Point2D< GeometryReal > pos )
	{
		switch( elementIndex )
		{
		case 0: return Point2D< GeometryReal >( 4 * pos[0] + 4 * pos[1] - 3 , 4 * pos[0] + 4 * pos[1] - 3 );
		case 1: return Point2D< GeometryReal >( 4 * pos[0] - 1 , 0 );
		case 2: return Point2D< GeometryReal >( 0 , 4 * pos[1] - 1 );
		case 3: return Point2D< GeometryReal >( 4 * pos[1] , 4 * pos[0]);
		case 4: return Point2D< GeometryReal >( -4 * pos[1] , -4 * pos[0] - 8 * pos[1] + 4 );
		case 5: return Point2D< GeometryReal >( -8 * pos[0] - 4 * pos[1] + 4 , -4 * pos[0] );
		default: MK_THROW( "Element out of bounds" );
		}
		return Point2D< GeometryReal >();
	}
	template< typename GeometryReal >
	void QuadraticElementValuesAndGradients( Point2D< GeometryReal > pos , GeometryReal values[] , Point2D< GeometryReal > gradients[] )
	{
		GeometryReal xx = pos[0]*pos[0] , xy = pos[0]*pos[1] , yy = pos[1]*pos[1] , x=pos[0] , y=pos[1];
		values[0] =  2*xx + 4*xy + 2*yy - 3*x - 3*y + 1 , gradients[0] = Point2D< GeometryReal >(  4*x + 4*y - 3 ,  4*x + 4*y - 3 );
		values[1] =  2*xx               - 1*x           , gradients[1] = Point2D< GeometryReal >(  4*x       - 1 ,              0 );
		values[2] =                2*yy       - 1*y     , gradients[2] = Point2D< GeometryReal >(              0 ,  4*y       - 1 );
		values[3] =         4*xy                        , gradients[3] = Point2D< GeometryReal >(        4*y     ,  4*x           );
		values[4] =        -4*xy - 4*yy       + 4*y     , gradients[4] = Point2D< GeometryReal >(      - 4*y     , -4*x - 8*y + 4 );
		values[5] = -4*xx - 4*xy        + 4*x           , gradients[5] = Point2D< GeometryReal >( -8*x - 4*y + 4 , -4*x           );
	}

	template< class Real , typename T >
	T BilinearValue( const T values[4] , Point2D< Real > pos )
	{
		Real x = pos[0] , y = pos[1];
		return
			values[0] * ( 1 - x ) * ( 1 - y ) +
			values[1] * (     x ) * ( 1 - y ) +
			values[2] * (     x ) * (     y ) +
			values[3] * ( 1 - x ) * (     y ) ;
	}
	template< class Real , typename T >
	Point2D< T > BilinearGradient( const T values[4] , Point2D< Real > pos )
	{
		Real x = pos[0] , y = pos[1];
		return
			Point2D< T >( values[0] * (   y - 1 ) , values[0] * (   x - 1 ) ) +
			Point2D< T >( values[1] * ( - y + 1 ) , values[1] * ( - x     ) ) +
			Point2D< T >( values[2] * (   y     ) , values[2] * (   x     ) ) +
			Point2D< T >( values[3] * ( - y     ) , values[3] * ( - x + 1 ) ) ;
	}
	// For the node at (0.0,0.0):
	//		F(x,y) = a x^2 + b y^2 + c xy + d x + e y + f
	//		F(0.0,0.0) = 1 => f = 1                 => F(x,y) = a x^2 + b y^2 + c xy + d x + e y + 1
	//		F(1.0,0.0) = 0 => a + d + 1 = 0         => F(x,y) = a x^2 + b y^2 + c xy - (a+1) x + e y + 1
	//		F(0.0,1.0) = 0 => b + e + 1 = 0         => F(x,y) = a x^2 + b y^2 + c xy - (a+1) x - (b+1) y + 1
	//		F(0.5,0.0) = 0 => a/4 - (a+1)/2 + 1 = 0 => F(x,y) = 2 x^2 + b y^2 + c xy - 3 x - (b+1) y + 1
	//		F(0.0,0.5) = 0 => b/4 - (b+1)/2 + 1 = 0 => F(x,y) = 2 x^2 + 2 y^2 + c xy - 3 x - 3 y + 1
	//		F(0.5,0.5) = 0 => c/4 - 1 = 0           => F(x,y) = 2 x^2 + 2 y^2 + 4 xy - 3 x - 3 y + 1
	// For the node at (1.0,0.0):
	//		F(x,y) = a x^2 + b y^2 + c xy + d x + e y + f
	//		F(0.0,0.0) = 0 => f = 0             => F(x,y) = a x^2 + b y^2 + c xy + d x + e y
	//		F(1.0,0.0) = 1 => a + d = 1         => F(x,y) = a x^2 + b y^2 + c xy + (1-a) x + e y
	//		F(0.0,1.0) = 0 => b + e = 0         => F(x,y) = a x^2 + b y^2 + c xy + (1-a) x - b y
	//		F(0.5,0.0) = 0 => a/4 + (1-a)/2 = 0 => F(x,y) = 2 x^2 + b y^2 + c xy - x - b y
	//		F(0.0,0.5) = 0 => b/4 - b/2 = 0     => F(x,y) = 2 x^2         + c xy - x
	//		F(0.5,0.5) = 0 => c/4 = 0           => F(x,y) = 2 x^2                - x
	// For the node at (0.0,1.0):
	//		F(x,y) = 2 y^2 - y
	// For the node at (0.5,0.0):
	//		F(x,y) = a x^2 + c xy + d x
	//		F(1.0,0.0) = 0 => a + d = 0     => F(x,y) =  a x^2 + c xy - a x
	//		F(0.5,0.0) = 1 => a/4 - a/2 = 1 => F(x,y) = -4 x^2 + c xy + 4 x
	//		F(0.5,0.5) = 0 => c/4 + 1 = 0   => F(x,y) = -4 x^2 - 4 xy + 4 x
	// For the node at (0.0,0.5):
	//		F(x,y) = -4 y^2 - 4 xy + 4y
	// For the node at (0.5,0.5):
	//		F(x,y) = c xy
	//		F(0.5,0.5) = 1 => c/4 = 1 => F(x,y) = 4 xy

	// 0 -> (0.0,0.0)
	// 1 -> (1.0,0.0)
	// 2 -> (0.0,1.0)
	// 3 -> (0.5,0.5)
	// 4 -> (0.0,0.5)
	// 5 -> (0.5,0.0)
	/////////////////////
	// Function values //
	//////////////////////////////////////////////////////////////////
	// (0.0,0.0) -> F(x,y) =   2 x^2 + 2 y^2 + 4 xy - 3 x - 3 y + 1 //
	// (1.0,0.0) -> F(x,y) =   2 x^2                - 1 x           //
	// (0.0,1.0) -> F(x,y) =           2 y^2              - 1 y     //
	// (0.5,0.5) -> F(x,y) =                   4 xy                 //
	// (0.0,0.5) -> F(x,y) =         - 4 y^2 - 4 xy       + 4 y     //
	// (0.5,0.0) -> F(x,y) = - 4 x^2         - 4 xy + 4 x           //
	//////////////////////////////////////////////////////////////////

	////////////////////////
	// Function gradients //
	/////////////////////////////////////////////////////////////////
	// (0.0,0.0) -> G(x,y) = (   4 x + 4 y - 3 ,   4 y + 4 x - 3 ) //
	// (1.0,0.0) -> G(x,y) = (   4 x       - 1 ,               0 ) //
	// (0.0,1.0) -> G(x,y) = (               0 ,   4 y       - 1 ) //
	// (0.5,0.5) -> G(x,y) = (         4 y     ,         4 x     ) //
	// (0.0,0.5) -> G(x,y) = (       - 4 y     , - 8 y - 4 x + 4 ) //
	// (0.5,0.0) -> G(x,y) = ( - 8 x - 4 y + 4 ,       - 4 x     ) //
	/////////////////////////////////////////////////////////////////
	template< class Real , typename T >
	T QuadraticValue( const T values[6] , Point2D< Real > pos )
	{
		Real xx = pos[0] * pos[0] , xy = pos[0]*pos[1] , yy = pos[1] * pos[1] , x = pos[0] , y = pos[1];
		return
			values[0] * (   2 * xx + 2 * yy + 4 * xy - 3 * x - 3 * y + 1 ) +
			values[1] * (   2 * xx                   - 1 * x             ) +
			values[2] * (            2 * yy                  - 1 * y     ) +
			values[3] * (                     4 * xy                     ) +
			values[4] * (          - 4 * yy - 4 * xy         + 4 * y     ) +
			values[5] * ( - 4 * xx          - 4 * xy + 4 * x             ) ;
	}
	template< class Real , typename T >
	Point2D< T > QuadraticGradient( const T values[6] , Point2D< Real > pos )
	{
		Real x = pos[0] , y = pos[1];
		return 
			Point2D< T >( values[0] * (   4 * x + 4 * y - 3 ) , values[0] * (   4 * y + 4 * x - 3 ) ) +
			Point2D< T >( values[1] * (   4 * x         - 1 ) , values[1] * (                   0 ) ) +
			Point2D< T >( values[2] * (                   0 ) , values[2] * (   4 * y         - 1 ) ) +
			Point2D< T >( values[3] * (           4 * y     ) , values[3] * (           4 * x     ) ) +
			Point2D< T >( values[4] * (         - 4 * y     ) , values[4] * ( - 8 * y - 4 * x + 4 ) ) +
			Point2D< T >( values[5] * ( - 8 * x - 4 * y + 4 ) , values[5] * (         - 4 * x     ) ) ;
	}

	template< typename GeometryReal >
	GeometryReal LinearElementValue( int elementIndex , Point2D< GeometryReal > pos )
	{
		switch( elementIndex )
		{
		case 0: return (GeometryReal)1. - pos[0] - pos[1];
		case 1: return pos[0];
		case 2: return pos[1];
		default: MK_THROW( "Element out of bounds" );
		}
		return (GeometryReal)0;
	}

	template< typename GeometryReal >
	Point2D< GeometryReal > LinearElementGradient( int elementIndex )
	{
		switch ( elementIndex )
		{
		case 0: return Point2D< GeometryReal >(-1,-1);
		case 1: return Point2D< GeometryReal >( 1, 0);
		case 2: return Point2D< GeometryReal >( 0, 1);
		default: MK_THROW( "Element out of bounds" );
		}
		return Point2D< GeometryReal >();
	}

	template< typename GeometryReal >
	GeometryReal BilinearElementValue( int elementIndex , Point2D<GeometryReal> pos )
	{
		switch( elementIndex )
		{
		case 0: return ( (GeometryReal)1.-pos[0]) * ( (GeometryReal)1.-pos[1]);
		case 1: return pos[0] * ( (GeometryReal)1.-pos[1]);
		case 2: return pos[0] * pos[1];
		case 3: return ( (GeometryReal)1.-pos[0]) * pos[1];
		default: MK_THROW( "Element out of bounds" );
		}
		return (GeometryReal)0;
	}

	template< typename GeometryReal >
	Point2D< GeometryReal > BilinearElementGradient( int elementIndex , Point2D< GeometryReal > pos)
	{
		switch( elementIndex )
		{
		case 0: return Point2D< GeometryReal >(pos[1]-(GeometryReal)1. , pos[0]-(GeometryReal)1. );
		case 1: return Point2D< GeometryReal >( (GeometryReal)1.-pos[1] , -pos[0] );
		case 2: return Point2D< GeometryReal >( pos[1] , pos[0] );
		case 3: return Point2D< GeometryReal >( -pos[1] , (GeometryReal)1.-pos[0] );
		default: MK_THROW( "Element out of bounds" );
		}
		return Point2D< GeometryReal >();
	}

	template< typename GeometryReal >
	void BilinearElementValuesAndGradients( Point2D< GeometryReal > pos , GeometryReal values[] , Point2D< GeometryReal > gradients[] )
	{
		GeometryReal x1 = (GeometryReal)1.-pos[0] , y1 = (GeometryReal)1.-pos[1] , x2 = pos[0] , y2 = pos[1];
		values[0] = x1*y1 , gradients[0] = Point2D< GeometryReal >(-y1,-x1);
		values[1] = x2*y1 , gradients[1] = Point2D< GeometryReal >( y1,-x2);
		values[2] = x2*y2 , gradients[2] = Point2D< GeometryReal >( y2, x2);
		values[3] = x1*y2 , gradients[3] = Point2D< GeometryReal >(-y2, x1);
	}

	template< typename GeometryReal >
	void ReducedVectorFieldBasis( Point2D< GeometryReal > pos , Point2D< GeometryReal > direction[] )
	{
		GeometryReal x = pos[0], y = pos[1];
		direction[0][0] = (GeometryReal)1.-y , direction[0][1] = (GeometryReal)0.  ;
		direction[1][0] = (GeometryReal)0.   , direction[1][1] = (GeometryReal)1.-x;
		direction[2][0] =                  y , direction[2][1] = (GeometryReal)0.  ;
		direction[3][0] = (GeometryReal)0.   , direction[3][1] =	              x;
	}
}