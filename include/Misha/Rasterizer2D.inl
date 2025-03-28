/*
Copyright (c) 2025, Michael Kazhdan
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


template< bool NodeAtCellCenter , unsigned int Degree , typename RasterizationFunctor /* = std::function< void ( Index ) > )*/ >
void _Rasterize( Triangle triangle , RasterizationFunctor && F , Range< 2 > cellRange );

template< bool NodeAtCellCenter , typename RasterizationFunctor /* = std::function< void ( Index ) > )*/ >
void RasterizeNodes( Triangle triangle , RasterizationFunctor && F , Range< 2 > cellRange )
{
	return _Rasterize< NodeAtCellCenter , (unsigned int)-1 >( triangle , std::forward< RasterizationFunctor >( F ) , cellRange );
}

template< bool Nearest , bool NodeAtCellCenter , typename RasterizationFunctor /* = std::function< void ( Index ) > )*/ >
void RasterizeSupports( Triangle triangle , RasterizationFunctor && F , Range< 2 > cellRange )
{
	if constexpr( Nearest ) return _Rasterize< NodeAtCellCenter , 0 >( triangle , std::forward< RasterizationFunctor >( F ) , cellRange );
	else                    return _Rasterize< NodeAtCellCenter , 1 >( triangle , std::forward< RasterizationFunctor >( F ) , cellRange );
}

template< bool NodeAtCellCenter >
constexpr double _Offset( void )
{
	if constexpr( NodeAtCellCenter ) return 0.5;
	else                              return 0.0;
}

template< unsigned int Degree >
constexpr double _SupportRadius( void )
{
	if constexpr( Degree==-1 ) return 0.;
	else					   return ( Degree + 1. ) / 2.;
}

// Get the 1D range of indices whose support overlaps the segment
template< bool NodeAtCellCenter , unsigned int Degree >
static Range< 1 > _GetCellRange( double s1 , double s2 )
{
	Range< 1 > range;
	// Solve for the smallest integer I s.t.:
	//	I + Offset + SupportRadius >= s1
	//	I >= s1 - Offset - SupportRadius
	range.first[0] = (int)std::ceil( s1 - _Offset< NodeAtCellCenter >() - _SupportRadius< Degree >() );

	// Solve for the largest integer I s.t.:
	//	I + Offset - SupportRadius < s2
	//	I < s2 - Offset + SupportRadius
	range.second[0] = (int)std::floor( s2 - _Offset< NodeAtCellCenter >() + _SupportRadius< Degree >() );
	if( ( s2 - _Offset< NodeAtCellCenter >() + _SupportRadius< Degree >() )==range.second[0] ) range.second[0]--;

	return range;
}

// Rasterizes the triangle obtained by connecting a horizontal line segment to a point
template< bool NodeAtCellCenter , unsigned int Degree , typename RasterizationFunctor /* = std::function< void ( Index ) > )*/ >
void _Rasterize( double y , double x0 , double x1 , Point< double , 2 > tip , RasterizationFunctor &&  F , const Range< 1 > cellRanges[2] )
{
	double y0 = y , y1 = tip[1];
	if( y0>y1 ) std::swap( y0 , y1 );
	if( x0>x1 ) std::swap( x0 , x1 );

	// For a given height, gives the span of the intersection of the horizontal line with the triangle
	auto Intersection = [&]( double _y )
		{
			_y = std::max< double >( y0 , std::min< double >( y1 , _y ) );
			// Solve for s s.t.:
			//	y*(1-s) + tip[1]*s = _y
			//  (tip[1]-y)*s = _y - y
			//	s = (_y-y) / (tip[1]-y)
			double s = (_y-y) / (tip[1]-y);
			return std::pair< double , double >( x0*(1-s) + tip[0]*s , x1*(1-s) + tip[0]*s );
		};

	// Returns the horizontal span of indices intersecting the triangle, for a fixed height index
	auto HorizontalCellRange = [&]( int iy )
		{
			if constexpr( Degree==-1 )
			{
				std::pair< double , double > r = Intersection( iy + _Offset< NodeAtCellCenter >() );
				return Range< 1 >::Intersect( cellRanges[0] , _GetCellRange< NodeAtCellCenter , Degree >( r.first , r.second ) );
			}
			else
			{
				std::pair< double , double > r1 = Intersection( iy + _Offset< NodeAtCellCenter >() - _SupportRadius< Degree >() );
				std::pair< double , double > r2 = Intersection( iy + _Offset< NodeAtCellCenter >() + _SupportRadius< Degree >() );
				return Range< 1 >::Intersect( cellRanges[0] , _GetCellRange< NodeAtCellCenter , Degree >( std::min< double >( r1.first , r2.first ) , std::max< double >( r1.second , r2.second ) ) );
			}
		};

	Range< 1 > iyRange = Range< 1 >::Intersect( cellRanges[1] , _GetCellRange< NodeAtCellCenter , Degree >( y0 , y1 ) );

	for( int iy=iyRange.first[0] ; iy<=iyRange.second[0] ; iy++ )
	{
		Range< 1 > ixRange = HorizontalCellRange( iy );
		for( int ix=ixRange.first[0] ; ix<=ixRange.second[0] ; ix++ ) F( Index(ix,iy) );
	}
}


template< bool NodeAtCellCenter , unsigned int Degree , typename RasterizationFunctor /* = std::function< void ( Index ) > )*/ >
void _Rasterize( Triangle triangle , RasterizationFunctor && F , Range< 2 > cellRange )
{
	static_assert( std::is_convertible_v< RasterizationFunctor , std::function< void ( Index ) > > , "[ERROR] RasterizationFunctor poorly formed" );

	// Split the 2D range into two 1D ranges
	Range< 1 > cellRanges[2];
	for( unsigned int d=0 ; d<2 ; d++ ) cellRanges[d].first[0] = cellRange.first[d] , cellRanges[d].second[0] = cellRange.second[d];

	// Find the index of the middle vertex
	unsigned int i1 = -1;
	for( unsigned int d=0 ; d<=2 ; d++ )
		if( ( triangle[d][1]<=triangle[(d+1)%3][1] && triangle[d][1]>=triangle[(d+2)%3][1] ) || ( triangle[d][1]>=triangle[(d+1)%3][1] && triangle[d][1]<=triangle[(d+2)%3][1] ) )
			i1 = d;
	if( i1==-1 ) MK_ERROR_OUT( "Could not find middle vertex: " , triangle );
	unsigned int i0 = (i1+2)%3 , i2 = (i1+1)%3;

	double x1 = triangle[i1][0] , y = triangle[i1][1];

	// All three vertices at the same height
	if     ( y==triangle[i0][1] && y==triangle[i2][1] ) return;
	// First and second vertices at the same height
	else if( y==triangle[i0][1]                       ) _Rasterize< NodeAtCellCenter , Degree >( y , x1 , triangle[i0][0] , triangle[i2] , F , cellRanges );
	// Second and third vertices at the same height
	else if(                       y==triangle[i2][1] ) _Rasterize< NodeAtCellCenter , Degree >( y , x1 , triangle[i2][0] , triangle[i0] , F , cellRanges );
	// All vertices at different heights
	else
	{
		// Solve for s s.t.:
		//	triangle[i0][1]*(1-s) + triangle[i2][1]*s = triangle[i1][1]
		//	(triangle[i2][1]-triangle[i0][1]) * s = triangle[i1][1] - triangle[i0][1]
		//	s = (triangle[i1][1] - triangle[i0][1]) / (triangle[i2][1]-triangle[i0][1])
		double s = (triangle[i1][1] - triangle[i0][1]) / (triangle[i2][1]-triangle[i0][1]);
		double x2 = triangle[i0][0]*(1-s) + triangle[i2][0]*s;
		_Rasterize< NodeAtCellCenter , Degree >( y , x1 , x2 , triangle[i0] , F , cellRanges );
		_Rasterize< NodeAtCellCenter , Degree >( y , x1 , x2 , triangle[i2] , F , cellRanges );
	}
}
