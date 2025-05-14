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
#include <Misha/SimplexBasis.h>
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

	template< typename IndexType >
	struct BilinearElementIndex
	{
		BilinearElementIndex( void ) { _v[0] = _v[1] = _v[2] = _v[3] = IndexType(-1); }
		BilinearElementIndex( IndexType v0 , IndexType v1 , IndexType v2 , IndexType v3 ) { _v[0] = v0 , _v[1] = v1 , _v[2] = v2 , _v[3] = v3; }
		IndexType &operator[]( unsigned int idx )       { return _v[idx]; }
		IndexType  operator[]( unsigned int idx ) const { return _v[idx]; }
	protected:
		IndexType _v[4];
	};

	struct QuadraticElement
	{
		struct Index
		{
			Index( void ) { _v[0] = _v[1] = _v[2] = _v[3] = _v[4] = _v[5] = AtlasInteriorOrBoundaryNodeIndex(0); }
			Index( AtlasInteriorOrBoundaryNodeIndex v0 , AtlasInteriorOrBoundaryNodeIndex v1 , AtlasInteriorOrBoundaryNodeIndex v2 , AtlasInteriorOrBoundaryNodeIndex v3 , AtlasInteriorOrBoundaryNodeIndex v4 , AtlasInteriorOrBoundaryNodeIndex v5 ){ _v[0] = v0 , _v[1] = v1 , _v[2] = v2 , _v[3] = v3 , _v[4] = v4 , _v[5] = v5; }
			AtlasInteriorOrBoundaryNodeIndex &operator[]( unsigned int idx )       { return _v[idx]; }
			AtlasInteriorOrBoundaryNodeIndex  operator[]( unsigned int idx ) const { return _v[idx]; }
		protected:
			AtlasInteriorOrBoundaryNodeIndex _v[6];
		};

		template< typename Real >
		static Real Value( unsigned int elementIndex , Point2D< Real > pos )
		{
			return static_cast< Real >( SimplexElements< 2 , 2 >::Element( _idx[elementIndex] )( Point2D< double >(pos) ) );
		}

		template< typename Real >
		static Point2D< Real > Differential( unsigned int elementIndex , Point2D< Real > pos )
		{
			const Point< Polynomial::Polynomial< 2 , 1 , double > , 2 > & d = SimplexElements< 2 , 2 >::Differential( _idx[elementIndex] );
			return Point2D< Real >( d[0]( Point2D< double >(pos) ) , d[1]( Point2D< double >(pos) ) );
		}

		template< typename Real >
		static void ValuesAndDifferentials( Point2D< Real > pos , Real values[6] , Point2D< Real > differentials[6] )
		{
			Point< double , 6 > _values = SimplexElements< 2 , 2 >::Elements()( pos );
			Point< Point< double , 2 > , 6 , double > _differentials = SimplexElements< 2 , 2 >::Differentials()( pos );
			for( unsigned int i=0 ; i<6 ; i++ ) values[i] = _values[ _idx[i] ] , differentials[i] = _differentials[ _idx[i] ];
		}

	protected:
		static const unsigned int _idx[];
	};

	inline const unsigned int QuadraticElement::_idx[] = 
	{
		SimplexElements< 2 , 2 >::NodeIndex( 0 , 0 ) ,
		SimplexElements< 2 , 2 >::NodeIndex( 1 , 1 ) ,
		SimplexElements< 2 , 2 >::NodeIndex( 2 , 2 ) ,
		SimplexElements< 2 , 2 >::NodeIndex( 1 , 2 ) ,
		SimplexElements< 2 , 2 >::NodeIndex( 2 , 0 ) ,
		SimplexElements< 2 , 2 >::NodeIndex( 0 , 1 ),
	};

	struct BilinearElement
	{
		template< class Real , typename T >
		static T Value( const T values[4] , Point2D< Real > pos )
		{
			Real x = pos[0] , y = pos[1];
			return
				values[0] * ( 1 - x ) * ( 1 - y ) +
				values[1] * (     x ) * ( 1 - y ) +
				values[2] * (     x ) * (     y ) +
				values[3] * ( 1 - x ) * (     y ) ;
		}

		template< class Real , typename T >
		static Point2D< T > Differential( const T values[4] , Point2D< Real > pos )
		{
			Real x = pos[0] , y = pos[1];
			return
				Point2D< T >( values[0] * (   y - 1 ) , values[0] * (   x - 1 ) ) +
				Point2D< T >( values[1] * ( - y + 1 ) , values[1] * ( - x     ) ) +
				Point2D< T >( values[2] * (   y     ) , values[2] * (   x     ) ) +
				Point2D< T >( values[3] * ( - y     ) , values[3] * ( - x + 1 ) ) ;
		}

	};

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