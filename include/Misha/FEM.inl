/*
Copyright (c) 2015, Michael Kazhdan and Fabian Prada
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

/////////
// FEM //
/////////

template< class Real >
inline Real FEM::Area( const SquareMatrix< Real , 2 >& tensor , const Point2D< Real > t[3] )
{
	return Dot( tensor , TangentVector< Real >( t[2] - t[0] ) , Rotate90( tensor , TangentVector< Real >( t[1] - t[0] ) ) ) / (Real)2;
}
template< class Real >
inline Real FEM::Area( const SquareMatrix< Real , 2 >& tensor , Point2D< Real > v1 , Point2D< Real > v2 , Point2D< Real > v3 )
{
	return Dot( tensor , TangentVector< Real >( v3 - v1 ) , Rotate90( tensor , TangentVector< Real >( v2 - v1 ) ) ) / (Real)2;
}
template< class Real >
inline Real FEM::Angle( const SquareMatrix< Real , 2 >& tensor , FEM::TangentVector< Real > v1 , FEM::TangentVector< Real > v2 )
{
	return (Real)acos( Dot( tensor , v1 , v2 ) / sqrt( SquareLength( tensor , v1 ) * SquareLength( tensor , v2 ) ) );
}
template< class Real >
inline Real FEM::Angle( const SquareMatrix< Real , 2 >& tensor , FEM::CotangentVector< Real > v1 , FEM::CotangentVector< Real > v2 )
{
	return (Real)acos( Dot( tensor , v1 , v2 ) / sqrt( SquareLength( tensor , v1 ) * SquareLength( tensor , v2 ) ) );
}
template< class Real >
inline SquareMatrix< Real , 2 > FEM::CoJ( const SquareMatrix< Real , 2 >& tensor )
{
	SquareMatrix< Real , 2 > J;
	J(0,1) = (Real)1 , J(1,0) = (Real)-1;
	return tensor * J / (Real)sqrt( tensor.determinant() );
}
template< class Real >
inline SquareMatrix< Real , 2 > FEM::J( const SquareMatrix< Real , 2 >& tensor )
{
	SquareMatrix< Real , 2 > J;
	J(0,1) = (Real)1 , J(1,0) = (Real)-1;
	return tensor.inverse() * J * (Real)sqrt( tensor.determinant() );
}
template< class Real > inline FEM::CotangentVector< Real > FEM::Rotate90( const SquareMatrix< Real , 2 >& tensor , FEM::CotangentVector< Real > v ){ return CoJ( tensor ) * v; }
template< class Real > inline FEM::  TangentVector< Real > FEM::Rotate90( const SquareMatrix< Real , 2 >& tensor , FEM::  TangentVector< Real > v ){ return   J( tensor ) * v; }

template< class Real > SquareMatrix< Real , 2 > FEM::MakeConformal( const SquareMatrix< Real , 2 >& sourceTensor , const SquareMatrix< Real , 2 >& targetTensor ){ return targetTensor * (Real)sqrt( sourceTensor.determinant() / targetTensor.determinant() ); }
template< class Real > SquareMatrix< Real , 2 > FEM::MakeAuthalic ( const SquareMatrix< Real , 2 >& sourceTensor , const SquareMatrix< Real , 2 >& targetTensor ){ return sourceTensor * (Real)sqrt( targetTensor.determinant() / sourceTensor.determinant() ); }
template< class Real >
SquareMatrix< Real , 2 > FEM::TensorRoot( const SquareMatrix< Real , 2 >& tensor )
{
	// Code borrowed from: https://en.wikipedia.org/wiki/Square_root_of_a_2_by_2_matrix
	SquareMatrix< Real , 2 > root = tensor;
	Real det = tensor.determinant();
	if( det<0 ) fprintf( stderr , "[ERROR] FEM::TensorRoot: Negative determinant: %g\n" , det ) , exit( 0 );
	Real s = (Real)sqrt( det );
	Real disc = (Real)( tensor.trace() + 2. * s );
	if( disc<0 ) fprintf( stderr , "[ERROR] FEM::TensorRoot: Negative discriminant: %g\n" , disc ) , exit( 0 );
	root(0,0) += s , root(1,1) += s;
	return root / (Real)sqrt( disc );
}
////////////////////////
// FEM::RightTriangle //
////////////////////////

template< class Real >
inline Point2D< Real > FEM::RightTriangle< Real >::Center( const SquareMatrix< Real , 2 >& tensor , int centerType )
{
	static const double SQRT_THREE_QUARTERS = sqrt( 3./4 );
	switch( centerType )
	{
	case CENTER_BARYCENTRIC:
		return Point2D< Real >( (Real)1./3 , (Real)1./3 );
	case CENTER_INCENTRIC:
		{
			Real lengths[] =
			{
				(Real)sqrt( Point2D< Real >::Dot( RightTriangle< Real >::EdgeDirections[0] , tensor * RightTriangle< Real >::EdgeDirections[0] ) ) ,
				(Real)sqrt( Point2D< Real >::Dot( RightTriangle< Real >::EdgeDirections[1] , tensor * RightTriangle< Real >::EdgeDirections[1] ) ) ,
				(Real)sqrt( Point2D< Real >::Dot( RightTriangle< Real >::EdgeDirections[2] , tensor * RightTriangle< Real >::EdgeDirections[2] ) )
			};
			Real lSum = lengths[0] + lengths[1] + lengths[2];
			return Point2D< Real >( lengths[1] / lSum , lengths[2] / lSum );
		}
	case CENTER_CIRCUMCENTRIC:
		{
			Real maxDet = 0;
			Point2D< Real > c;
			for( int j=0 ; j<3 ; j++ )
			{
				Point2D< Real > c1 = RightTriangle< Real >::EdgeMidpoints[(j+1)%3];
				Point2D< Real > c2 = RightTriangle< Real >::EdgeMidpoints[(j+2)%3];
				TangentVector< Real > v1 = Rotate90( tensor , RightTriangle< Real >::EdgeDirections[(j+1)%3] );
				TangentVector< Real > v2 = Rotate90( tensor , RightTriangle< Real >::EdgeDirections[(j+2)%3] );
				// Solve for s and t such that:
				//		c1 + s * v1 = c2 + t * v2
				// =>	c1 - c2 = - s * v1 + t * v2
				// =>	c1 - c2 = ( -v1 | v2 ) * ( s , t )^t
				SquareMatrix< Real , 2 > M;
				M( 0 , 0 ) = -v1[0] , M( 0 , 1 ) = -v1[1];
				M( 1 , 0 ) =  v2[0] , M( 1 , 1 ) =  v2[1];
				Real det = (Real)fabs( M.determinant() );
				if( det>maxDet )
				{
					Point2D< Real > x = M.inverse() * ( c1 - c2 );
					c = ( c1 + v1 * x[0] + c2 + v2 * x[1] ) / static_cast< Real >( 2 );
					maxDet = det;
				}
			}
			return c;
		}
	case CENTER_ISOGONIC:
		{
			Point2D< Real > eVertices[3];
			// Solve for the positive value of t such that:
			//		|| Corners[j+1] - ( midPoint + t * perpDir ) ||^2 = || Corners[j+1] - Corners[j+2] ||^2
			// <=>	t^2 || perpDir ||^2 = || Corners[j+1] - Corners[j+2] ||^2 - || Corners[j+1] - midPoint ||^2
			// <=>	t^2 || Corners[j+1] - Corners[j+2] ||^2 = || Corners[j+1] - Corners[j+2] ||^2 - || ( Corners[j+1] - Corners[j+2] )/2 ||^2
			// <=>	t^2 || Corners[j+1] - Corners[j+2] ||^2 = 3/4 || Corners[j+1] - Corners[j+2] ||^2
			// <=>	t = sqrt( 3/4 )
			for( int j=0 ; j<3 ; j++ ) eVertices[j] = RightTriangle< Real >::EdgeMidpoints[j] - Rotate90( tensor ,  RightTriangle< Real >::EdgeDirections[j] ) * (Real)SQRT_THREE_QUARTERS;
			Real maxDet = 0;
			Point2D< Real > c;
			for( int j=0 ; j<3 ; j++ )
			{
				Point2D< Real > c1 = eVertices[(j+1)%3] , c2 = eVertices[(j+2)%3];
				Point2D< Real > v1 = RightTriangle< Real >::Corners[(j+1)%3] - c1 , v2 = RightTriangle< Real >::Corners[(j+2)%3] - c2;

				// Solve for s and t such that:
				//		c1 + s * v1 = c2 + t * v2
				// =>	c1 - c2 = - s * v1 + t * v2
				// =>	c1 - c2 = ( -v1 | v2 ) * ( s , t )^t
				SquareMatrix< Real , 2 > M;
				M( 0 , 0 ) = -v1[0] , M( 0 , 1 ) = -v1[1];
				M( 1 , 0 ) =  v2[0] , M( 1 , 1 ) =  v2[1];
				Real det = (Real)fabs( M.determinant() );
				if( det>maxDet )
				{
					Point2D< Real > x = M.inverse() * ( c1 - c2 );
					c = ( c1 + v1 * x[0] + c2 + v2 * x[1] ) / static_cast< Real >( 2 );
					maxDet = det;
				}
			}
			return c;
		}
	default:
		fprintf( stderr , "[ERROR] FEM::RightTriangle::Center: Unrecognized center type: %d\n" , centerType ) , exit( 0 );
	}
}
template< class Real >
inline Point3D< Real > FEM::RightTriangle< Real >::SubTriangleAreas( const SquareMatrix< Real , 2 >& tensor , Point2D< Real > center )
{
	Point2D< Real > triangle[3];
	triangle[2] = center;
	Point3D< Real > areas;
	for( int i=0 ; i<3 ; i++ )
	{
		triangle[0] = RightTriangle< Real >::Corners[(i+1)%3] , triangle[1] = RightTriangle< Real >::Corners[(i+2)%3];
		areas[i] = Area( tensor , triangle );
	}
	return areas;
}
template< class Real >
inline Point3D< Real > FEM::RightTriangle< Real >::CenterAreas( const SquareMatrix< Real , 2 >& tensor , int centerType ){ return SubTriangleAreas( tensor , Center( tensor , centerType ) ); }

template< class Real >
inline Point2D< Real > FEM::RightTriangle< Real >::EdgeReflect( const SquareMatrix< Real , 2 >& tensor , int e , Point2D< Real > p )
{
	Point2D< Real > c = Corners[(e+1)%3] , v = p - c ,  perp = Rotate90( tensor , EdgeDirections[e] );
	return c + v - ( 2 * Point2D< Real >::Dot( perp , tensor * v ) / Point2D< Real >::Dot( perp , tensor * perp ) ) * perp;
}

template< class Real >
template< unsigned int BasisType , class V >
inline V FEM::RightTriangle< Real >::EvaluateScalarField( ConstPointer( V ) c , const Point2D< Real >& p )
{
	switch( BasisType )
	{
		case BASIS_0_WHITNEY: return c[0] * (Real)( 1. - p[0] - p[1] ) + c[1] * p[0] + c[2] * p[1];
		default: fprintf( stderr , "[ERROR] FEM::RightTriangle::EvaluateScalarField: unrecognized basis type: %d\n" , BasisType ) , exit( 0 );
	}
}
template< class Real >
template< unsigned int BasisType , class V >
FEM::CotangentVector< V > FEM::RightTriangle< Real >::EvaluateScalarFieldDifferential( ConstPointer( V ) c , const Point2D< Real >& p )
{
	switch( BasisType )
	{
		case BASIS_0_WHITNEY: return CotangentVector< V >( c[0]*CornerDifferentials[0][0] + c[1]*CornerDifferentials[1][0] + c[2]*CornerDifferentials[2][0] , c[0]*CornerDifferentials[0][1] + c[1]*CornerDifferentials[1][1] + c[2]*CornerDifferentials[2][1] );
		default: fprintf( stderr , "[ERROR] FEM::RightTriangle::EvaluateScalarFieldDifferential: unrecognized basis type: %d\n" , BasisType ) , exit( 0 );
	}
}

template< class Real >
template< unsigned int BasisType , class V >
FEM::TangentVector< V > FEM::RightTriangle< Real >::EvaluateScalarFieldGradient( const SquareMatrix< Real , 2 >& tensor , ConstPointer( V ) c , const Point2D< Real >& p )
{
	return Dual( tensor , EvaluateScalarFieldDifferential< BasisType >( c , p ) );
}

template< class Real >
template< unsigned int BasisType , class V >
inline FEM::CotangentVector< V > FEM::RightTriangle< Real >::EvaluateCovectorField( const SquareMatrix< Real , 2 >& tensor , ConstPointer( V ) c , const Point2D< Real >& p )
{
	switch( BasisType )
	{
		case BASIS_1_CONFORMING:
		{
			SquareMatrix< Real , 2 > J = CoJ( tensor );
			CotangentVector< Real > g[] = { CornerDifferentials[0] , J *CornerDifferentials[0] , CornerDifferentials[1] , J * CornerDifferentials[1] , CornerDifferentials[2] , J * CornerDifferentials[2] };
			return CotangentVector< V >( c[0]*g[0][0] + c[1]*g[1][0] + c[2]*g[2][0] + c[3]*g[3][0] + c[4]*g[4][0] + c[5]*g[5][0] , c[0]*g[0][1] + c[1]*g[1][1] + c[2]*g[2][1] + c[3]*g[3][1] + c[4]*g[4][1] + c[5]*g[5][1] );
		}
		case BASIS_1_WHITNEY:
		{
			Point3D< Real > v( (Real)( 1. - p[0] - p[1] ) , p[0] , p[1] );
			V n[] = { c[1]*v[2] - c[2]*v[1] , c[2]*v[0] - c[0]*v[2] , c[0]*v[1] - c[1]*v[0] };
#if 1
			return CotangentVector< V >( n[0]*CornerDifferentials[0][0] + n[1]*CornerDifferentials[1][0] + n[2]*CornerDifferentials[2][0] , n[0]*CornerDifferentials[0][1] + n[1]*CornerDifferentials[1][1] + n[2]*CornerDifferentials[2][1] );
#else
			return CotangentVector< V >( n[0]*CornerDifferentials[0][0] + n[1]*CornerDifferentials[1][0] + n[2]*CornerDifferentials[2][0] , n[0]*CornerDifferentials[0][0] + n[1]*CornerDifferentials[1][1] + n[2]*CornerDifferentials[2][1] );
#endif
		}
		case BASIS_1_TRIANGLE_CONSTANT: return CotangentVector< V >( c[0] , c[1] );
		default: fprintf( stderr , "[ERROR] FEM::RightTriangle::EvaluateCovectorField: unrecognized basis type: %d\n" , BasisType ) , exit( 0 );
	}
}
template< class Real >
template< unsigned int BasisType , class V >
inline V FEM::RightTriangle< Real >::EvaluateDensityField( const SquareMatrix< Real , 2 >& tensor , ConstPointer( V ) c , const Point2D< Real >& p )
{
	switch( BasisType )
	{
		case BASIS_2_WHITNEY: return c[0] / Area( tensor );
		default: fprintf( stderr , "[ERROR] FEM::RightTriangle::EvaluateDensityField: unrecognized basis type: %d\n" , BasisType ) , exit( 0 );
	}
}
template< class Real >
template< unsigned int InBasisType , unsigned int OutBasisType >
typename FEM::BasisInfoSystem2< Real , InBasisType , OutBasisType >::MaskMatrix FEM::RightTriangle< Real >::GetDMask( bool& redundant )
{
	auto Fail = [&] ( void ){ fprintf( stderr , "[ERROR] FEM::RightTriangle::GetDMask: %s -> %s\n" , BasisNames[ InBasisType ] , BasisNames[ OutBasisType ] ) , exit( 0 ); };

	TestBasisType(  InBasisType , "FEM::RightTriangle::GetDMask" , false );
	TestBasisType( OutBasisType , "FEM::RightTriangle::GetDMask" , false );

	Matrix< unsigned char , FEM::BasisInfo< InBasisType >::Coefficients , FEM::BasisInfo< OutBasisType >::Coefficients > D;
	D *= 0;

	if( InBasisType==BASIS_0_WHITNEY )
	{
		if     ( OutBasisType==BASIS_1_WHITNEY           ) for( int e=0 ; e<3 ; e++ ) D( (e+1)%3 , e ) = D( (e+2)%3 , e ) = 1;
		else if( OutBasisType==BASIS_1_CONFORMING        ) for( int v=0 ; v<3 ; v++ ) D( v , v*2 ) = 1;
		else if( OutBasisType==BASIS_1_TRIANGLE_CONSTANT ) for( int i=0 ; i<2 ; i++ ) D( 0 , i ) = D( 1+i , i ) = 1;
		else Fail();
	}
	else if( OutBasisType==BASIS_2_VERTEX_CONSTANT )
	{
		if     ( InBasisType==BASIS_1_CONFORMING                                        ) for( int v=0 ; v<3 ; v++ ) for( int i=0 ; i<3 ; i++ ) D( i*2+1 , v ) = 1;
		else if( InBasisType==BASIS_1_WHITNEY || InBasisType==BASIS_1_TRIANGLE_CONSTANT ) for( int v=0 ; v<3 ; v++ ) for( int i=0 ; i<BasisInfo< InBasisType >::Coefficients ; i++ ) D(i,v) = 1;
		else Fail();
	}
	else if( InBasisType==BASIS_1_WHITNEY )
	{
		if( OutBasisType==BASIS_2_WHITNEY ) D(0,0) = D(1,0) = D(2,0) = 1;
		else Fail();
	}
	else Fail();
	redundant = ( OutBasisType == BASIS_1_WHITNEY || OutBasisType==BASIS_1_CONFORMING );
	return D;
}
template< class Real >
template< unsigned int InBasisType , unsigned int OutBasisType >
typename FEM::BasisInfoSystem2< Real , InBasisType , OutBasisType >::Matrix FEM::RightTriangle< Real >::GetDMatrix( const SquareMatrix< Real , 2 >& tensor )
{
	auto Fail = [&] ( void ){ fprintf( stderr , "[ERROR] FEM::RightTriangle::GetDMatrix: %s -> %s\n" , BasisNames[ InBasisType ] , BasisNames[ OutBasisType ] ) , exit( 0 ); };

	TestBasisType(  InBasisType , "FEM::RightTriangle::GetDMatrix" , false );
	TestBasisType( OutBasisType , "FEM::RightTriangle::GetDMatrix" , false );

	Matrix< Real , FEM::BasisInfo< InBasisType >::Coefficients , FEM::BasisInfo< OutBasisType >::Coefficients > D;

	if( InBasisType==BASIS_0_WHITNEY )
	{
		if     ( OutBasisType==BASIS_1_WHITNEY           ) for( int e=0 ; e<3 ; e++ ) D( (e+1)%3 , e ) = -1 , D( (e+2)%3 , e ) = 1;
		else if( OutBasisType==BASIS_1_CONFORMING        ) for( int v=0 ; v<3 ; v++ ) D( v , v*2 ) = 1;
		else if( OutBasisType==BASIS_1_TRIANGLE_CONSTANT ) for( int i=0 ; i<2 ; i++ ) D( 0 , i ) = -1 , D( 1+i , i ) = 1;
		else Fail();
	}
	else if( OutBasisType==BASIS_2_VERTEX_CONSTANT )
	{
		if( InBasisType==BASIS_1_WHITNEY || InBasisType==BASIS_1_CONFORMING || InBasisType==BASIS_1_TRIANGLE_CONSTANT )
		{
			Point2D< Real > center = Center( tensor , CENTER_CIRCUMCENTRIC );
			for( int v=0 ; v<3 ; v++ )
			{
				Point< Real , BasisInfo< InBasisType >::Coefficients > d = IntegrationDual< InBasisType >( tensor , EdgeMidpoints[(v+2)%3] , center ) + IntegrationDual< InBasisType >( tensor , center , EdgeMidpoints[(v+1)%3] );
				for( int i=0 ; i<BasisInfo< InBasisType >::Coefficients ; i++ ) D(i,v) = d[i];
			}

		}
		else Fail();
	}
	else if( InBasisType==BASIS_1_WHITNEY )
	{
		if( OutBasisType==BASIS_2_WHITNEY ) D(0,0) = D(1,0) = D(2,0) = 1;
		else Fail();
	}
	else Fail();
	return D;
}
template< class Real >
template< unsigned int BasisType >
typename FEM::BasisInfoSystem< Real , BasisType >::Matrix FEM::RightTriangle< Real >::GetMassMatrix( const SquareMatrix< Real , 2 >& tensor )
{
	SquareMatrix< Real , BasisInfo< BasisType >::Coefficients > mass;

	switch( BasisType )
	{
		case BASIS_0_WHITNEY:
		{
			{
				mass(0,0) = mass(1,1) = mass(2,2) = (Real)1./6;
				mass(1,0) = mass(0,1) = mass(1,2) = mass(2,1) = mass(0,2) = mass(2,0) = (Real)1./12;	
				mass *= (Real)( sqrt( tensor.determinant() ) / 2. );
			}
			break;
		}
		case BASIS_1_CONFORMING:
		{
			{
				SquareMatrix< Real , 2 > iTensor = tensor.inverse() * (Real)( sqrt( tensor.determinant() ) / 2. );
				for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) mass( 2*i , 2*j ) = mass( 2*i+1 , 2*j+1 ) = Point2D< Real >::Dot( CornerDifferentials[i] , iTensor * CornerDifferentials[j] );
			}
			break;
		}
		case BASIS_1_WHITNEY:
		{
			{
				SquareMatrix< Real , 3 > M = GetMassMatrix< BASIS_0_WHITNEY >( tensor ) , S;
				SquareMatrix< Real , 2 > iTensor = tensor.inverse();
				for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) S(i,j) = Point2D< Real >::Dot( CornerDifferentials[i] , iTensor * CornerDifferentials[j] );
		
				const int I[][2] = { {1,2} , {2,0} , {0,1} };
				for( int i=0 ; i<3;  i++ ) for( int j=0 ; j<3 ; j++ )
					mass(i,j) = M( I[i][0] , I[j][0] ) * S( I[i][1] , I[j][1] ) + M( I[i][1] , I[j][1] ) * S( I[i][0] , I[j][0] ) - M( I[i][0] , I[j][1] ) * S( I[i][1] , I[j][0] ) - M( I[i][1] , I[j][0] ) * S( I[i][0] , I[j][1] );
			}
			break;
		}
		case BASIS_1_TRIANGLE_CONSTANT:
		{
			SquareMatrix< Real , 2 > iTensor = tensor.inverse();
			for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) mass(i,j) = iTensor(i,j);
			mass *= (Real)( sqrt( tensor.determinant() ) / 2. );
			break;
		}
		case BASIS_2_WHITNEY:
		{
			mass(0,0) = (Real)( 1./(Real)( sqrt( tensor.determinant() ) / 2. ) );
		}
		break;
		default: TestBasisType( BasisType , "FEM::RightTriangle::GetMassMatrix" , true );
	}

	return mass;
}
template< class Real >
template< unsigned int BasisType >
typename FEM::BasisInfoSystem< Real , BasisType >::Matrix FEM::RightTriangle< Real >::GetMassMatrix( const SquareMatrix< Real , 2 >& tensor , const SquareMatrix< Real , 2 >& newTensor )
{
	SquareMatrix< Real , BasisInfo< BasisType >::Coefficients > mass;
	SquareMatrix< Real , 2 > tensor_i = tensor.inverse() , newCotensor = tensor_i.transpose() * newTensor * tensor_i;

	switch( BasisType )
	{
		case BASIS_1_CONFORMING:
		{
			CotangentVector< Real > grads[] = { CornerDifferentials[0] , Rotate90( tensor , CornerDifferentials[0] ) , CornerDifferentials[1] , Rotate90( tensor , CornerDifferentials[1] ) , CornerDifferentials[2] , Rotate90( tensor , CornerDifferentials[2] ) };
			for( int i=0 ; i<6 ; i++ ) for( int j=0 ; j<6 ; j++ ) mass(i,j) = Point2D< Real >::Dot( grads[i]  , newCotensor * grads[j] );
			mass *= (Real)( sqrt( tensor.determinant() ) / 2. );
			break;
		}
		case BASIS_1_TRIANGLE_CONSTANT:
		{
			for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) mass(i,j) = newCotensor(i,j);
			mass *= (Real)( sqrt( tensor.determinant() ) / 2. );
			break;
		}
		case BASIS_1_WHITNEY:
		{
			SquareMatrix< Real , 3 > M = GetMassMatrix< BASIS_0_WHITNEY >( tensor ) , S;
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) S(i,j) = Point2D< Real >::Dot( CornerDifferentials[i] , newCotensor * CornerDifferentials[j] );
		
			const int I[][2] = { {1,2} , {2,0} , {0,1} };
			for( int i=0 ; i<3;  i++ ) for( int j=0 ; j<3 ; j++ )
				mass(i,j) = M( I[i][0] , I[j][0] ) * S( I[i][1] , I[j][1] ) + M( I[i][1] , I[j][1] ) * S( I[i][0] , I[j][0] ) - M( I[i][0] , I[j][1] ) * S( I[i][1] , I[j][0] ) - M( I[i][1] , I[j][0] ) * S( I[i][0] , I[j][1] );
			break;
		}
		default: TestBasisType( BasisType , "FEM::RightTriangle::GetCovectorMassMatrix" , true );
	}

	return mass;
}
template< class Real >
template< unsigned int BasisType >
typename FEM::BasisInfoSystem< Real , BasisType >::Point FEM::RightTriangle< Real >::GetDiagonalMassMatrix( const SquareMatrix< Real , 2 >& tensor )
{
	Point< Real , BasisInfo< BasisType >::Coefficients > mass;

	switch( BasisType )
	{
		case BASIS_0_WHITNEY:
			mass[0] = mass[1] = mass[2] = (Real)( sqrt( tensor.determinant() ) / 6. );
			break;
		case BASIS_1_CONFORMING:
		{
			Point3D< Real > areas = CenterAreas( tensor , CENTER_CIRCUMCENTRIC ) * static_cast< Real >( 2 );
			mass[0] = mass[1] = areas[0] / SquareLength( tensor , EdgeDirections[0] );
			mass[2] = mass[3] = areas[1] / SquareLength( tensor , EdgeDirections[1] );
			mass[4] = mass[5] = areas[2] / SquareLength( tensor , EdgeDirections[2] );
			break;
		}
		case BASIS_1_WHITNEY:
		{
			Point3D< Real > areas = CenterAreas( tensor , CENTER_CIRCUMCENTRIC ) * static_cast< Real >( 2 );
			mass[0] = areas[0] / SquareLength( tensor , EdgeDirections[0] );
			mass[1] = areas[1] / SquareLength( tensor , EdgeDirections[1] );
			mass[2] = areas[2] / SquareLength( tensor , EdgeDirections[2] );
			break;
		}
		case BASIS_2_WHITNEY:
			mass[0] = (Real)( 1./( sqrt( tensor.determinant() ) / 2. ) );
			break;
		default: TestBasisType( BasisType , "FEM::RightTriangle::GetDiagonalMassMatrix" , true );
	}

	return mass;
}
template< class Real >
template< unsigned int BasisType >
typename FEM::BasisInfoSystem< Real , BasisType >::MaskMatrix FEM::RightTriangle< Real >::GetMassMask( bool useTensor )
{
	SquareMatrix< unsigned char , BasisInfo< BasisType >::Coefficients > M;
	M *= 0;

	if( useTensor )
		switch( BasisType )
		{
			case BASIS_0_WHITNEY:           for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) M(i,j) = 1 ; break;
			case BASIS_1_CONFORMING:        for( int i=0 ; i<6 ; i++ ) for( int j=0 ; j<6 ; j++ ) M(i,j) = 1 ; break;
			case BASIS_1_WHITNEY:           for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) M(i,j) = 1 ; break;
			case BASIS_1_TRIANGLE_CONSTANT: for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) M(i,j) = 1 ; break;
			case BASIS_2_WHITNEY:                                                                 M(0,0) = 1 ; break;
			default: TestBasisType( BasisType , "FEM::RightTriangle::GetMassMask" , true );
		}
	else
		switch( BasisType )
		{
			case BASIS_0_WHITNEY:           for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) M(i,j) = 1 ; break;
			case BASIS_1_CONFORMING:        for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) M(2*i,2*j) = M(2*i+1,2*j+1) = 1 ; break;
			case BASIS_1_WHITNEY:           for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) M(i,j) = 1 ; break;
			case BASIS_1_TRIANGLE_CONSTANT: for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) M(i,j) = 1 ; break;
			case BASIS_2_WHITNEY:                                                                 M(0,0) = 1 ; break;
			default: TestBasisType( BasisType , "FEM::RightTriangle::GetMassMask" , true );
		}
	return M;
}

template< class Real >
template< unsigned int BasisType >
Real FEM::RightTriangle< Real >::Integrate( const SquareMatrix< Real , 2 >& tensor , ConstPointer( Real ) linear )
{
	Real integral = 0;
	switch( BasisType )
	{
		case BASIS_0_WHITNEY:
		{
			SquareMatrix< Real , 3 > mass = GetMassMatrix< BASIS_0_WHITNEY >( tensor );
			for( int i=0 ; i<3 ; i++ ) integral += linear[i] * ( mass(i,0) + mass(i,1) + mass(i,2) );
			break;
		}
		default: TestBasisType( BasisType , "FEM::RightTriangle::Integrate" , true );
	}
	return integral;
}

template< typename Real >
template< unsigned int Degree >
Real FEM::RightTriangle< Real >::Integrate( const SquareMatrix< Real , 2 > &tensor , const Polynomial::Polynomial< 2 , Degree , Real > &p )
{
	return p.integrateUnitRightSimplex() * sqrt( tensor.determinant() );
}


template< class Real >
template< unsigned int BasisType >
typename FEM::BasisInfoSystem< Real , BasisType >::Point FEM::RightTriangle< Real >::IntegrationDual( const SquareMatrix< Real , 2 >& tensor , ConstPointer( Real ) linear )
{
	Point< Real , BasisInfo< BasisType >::Coefficients > integralDual;
	switch( BasisType )
	{
		case BASIS_0_WHITNEY:
		{
			SquareMatrix< Real , 3 > mass = GetMassMatrix< BASIS_0_WHITNEY >( tensor );
			for( int i=0 ; i<3 ; i++ ) for( int j=0 ; j<3 ; j++ ) integralDual[i] += linear[j] * mass(i,j);
			break;
		}
		default: TestBasisType( BasisType , "FEM::RightTriangle::IntegrationDual" , true );
	}
	return integralDual;
}
template< class Real >
template< unsigned int BasisType >
typename FEM::BasisInfoSystem< Real , BasisType >::Point FEM::RightTriangle< Real >::IntegrationDual( const SquareMatrix< Real , 2 >& tensor , ConstPointer( CotangentVector< Real > ) linear )
{
	Point< Real , BasisInfo< BasisType >::Coefficients > integralD;
	switch( BasisType )
	{
		case BASIS_1_TRIANGLE_CONSTANT:
		{
			SquareMatrix< Real , 2 > tensor_i = tensor.inverse();
			Point2D< Real > dualBasis[] = { tensor_i * Point2D< Real >( 1 , 0 ) , tensor_i * Point2D< Real >( 0 , 1 ) };
			Real _linear[3];
			for( int i=0 ; i<2 ; i++ )
			{
				for( int j=0 ; j<3 ; j++ ) _linear[j] = Point2D< Real >::Dot( dualBasis[i] , linear[j] );
				integralD[i] = Integrate< BASIS_0_WHITNEY >( tensor , ( ConstPointer( Real ) )GetPointer( _linear , 3 ) );
			}
			break;
		}
		case BASIS_1_CONFORMING:
		{
			SquareMatrix< Real , 2 > tensor_i = tensor.inverse();
			Point2D< Real > dualBasis[] =
			{
				tensor_i * CornerDifferentials[0] , tensor_i * Rotate90( tensor , CornerDifferentials[0] ) ,
				tensor_i * CornerDifferentials[1] , tensor_i * Rotate90( tensor , CornerDifferentials[1] ) ,
				tensor_i * CornerDifferentials[2] , tensor_i * Rotate90( tensor , CornerDifferentials[2] )
			};
			Real _linear[3];
			for( int i=0 ; i<6 ; i++ )
			{
				for( int j=0 ; j<3 ; j++ ) _linear[j] = Point2D< Real >::Dot( dualBasis[i] , linear[j] );
				integralD[i] = Integrate< BASIS_0_WHITNEY >( tensor , ( ConstPointer( Real ) )GetPointer( _linear , 3 ) );
			}
			break;
		}
		case BASIS_1_WHITNEY:
		{
			// \phi_{ij}(p) = \phi_i(p) \nabla \phi_j - \phi_j(p) \nabla \phi_i
			// phi[i](p) = hat[i+1](p) * CornerDifferentials[i+2] - hat[i+2](p) * CornerDifferentials[i+1]
			// < phi[i](p) , L(p) > = \sum_j < phi[i](p) , linear[j] * hat[j](p) >
			//                      = \sum_j < hat[i+1](p) * CornerDifferentials[i+2] - hat[i+2](p) * CornerDifferentials[i+1] , linear[j] * hat[j](p) >
			//                      = \sum_j < CornerDifferentials[i+2] , linear[j] > * Integral{ hat[i+1](p) * hat[j](p) } - < CornerDifferentials[i+1] , linear[j] >  * Integral{ hat[i+2](p) * hat[j](p) }
			SquareMatrix< Real , 3 > mass = GetMassMatrix< BASIS_0_WHITNEY >( tensor );
			SquareMatrix< Real , 2 > tensor_i = tensor.inverse();
			for( int j=0 ; j<3 ; j++ )
			{
				CotangentVector< Real > linearD = tensor_i * linear[j];
				for( int i=0 ; i<3 ; i++ ) 
				{
					integralD[i] += mass( (i+1)%3 , j ) * Point2D< Real >::Dot( CornerDifferentials[(i+2)%3] , linearD );
					integralD[i] -= mass( (i+2)%3 , j ) * Point2D< Real >::Dot( CornerDifferentials[(i+1)%3] , linearD );
				}
			}
			break;
		}
		default: TestBasisType( BasisType , "FEM::RightTriangle::IntegrationDual" , true );
	}
	return integralD;
}
template< class Real >
template< unsigned int BasisType >
typename FEM::BasisInfoSystem< Real , BasisType >::Point FEM::RightTriangle< Real >::IntegrationDual( const SquareMatrix< Real , 2 >& tensor , Point2D< Real > p , Point2D< Real > q )
{
	Point< Real , BasisInfo< BasisType >::Coefficients > integralDual , coefficients;
	switch( BasisType )
	{
		case BASIS_1_WHITNEY:
		case BASIS_1_CONFORMING:
		case BASIS_1_TRIANGLE_CONSTANT:
		{
			TangentVector< Real >dir = q - p;
			for( int i=0 ; i<BasisInfo< BasisType >::Coefficients ; i++ )
			{
				coefficients[i] = 1;
				CotangentVector< Real > _p = EvaluateCovectorField< BasisType >( tensor , ( ConstPointer( Real ) )GetPointer( &coefficients[0] , BasisInfo< BasisType >::Coefficients ) , p );
				CotangentVector< Real > _q = EvaluateCovectorField< BasisType >( tensor , ( ConstPointer( Real ) )GetPointer( &coefficients[0] , BasisInfo< BasisType >::Coefficients ) , q );
				integralDual[i] = ( Point2D< Real >::Dot( dir , _p ) + Point2D< Real >::Dot( dir , _q ) ) / 2;
				coefficients[i] = 0;
			}
			break;
		}
		default: TestBasisType( BasisType , "FEM::RightTriangle::IntegrationDual" , true );
	}
	return integralDual;
}

////////////////////////////////////////////////////
// FEM::RiemannianMesh::TangentVectorFieldWrapper //
////////////////////////////////////////////////////
template< class Real , unsigned int BasisType >
template< typename Index >
FEM::TangentVectorFieldWrapper< Real , BasisType >::TangentVectorFieldWrapper( const RiemannianMesh< Real , Index >* mesh , ConstPointer( Real ) coefficients , bool precomputeInverses ) : _mesh(mesh) , _coefficients(coefficients)
{
	if( precomputeInverses )
	{
		_gInverse = AllocPointer< SquareMatrix< Real , 2 > >( _mesh->tCount() );
		ThreadPool::ParallelFor( 0 , _mesh->tCount() , [&]( unsigned int , size_t t ){ _gInverse[t] = _mesh->g(t).inverse(); } );
	}
	else _gInverse = NullPointer< SquareMatrix< Real , 2 > >();
}

template< class Real , unsigned int BasisType >
FEM::TangentVectorFieldWrapper< Real , BasisType >::~TangentVectorFieldWrapper( void ){ FreePointer( _gInverse ); }

template< class Real , unsigned int BasisType >
FEM::TangentVector< Real > FEM::TangentVectorFieldWrapper< Real , BasisType >::operator() ( const FEM::SamplePoint< Real >& p ) const
{
	bool isAligned;
	Real c[BasisInfo< BasisType >::Coefficients];
	for( int i=0 ; i<BasisInfo< BasisType >::Coefficients ; i++ )
	{
		int idx =_mesh->template index< BasisType >( p.tIdx , i , isAligned );
		c[i] = isAligned ? _coefficients[ idx ] : -_coefficients[ idx ] ;
	}
	if( _gInverse ) return _gInverse[ p.tIdx ]          * RightTriangle< Real >::template EvaluateCovectorField< BasisType >( _mesh->g(p.tIdx) , ( ConstPointer( Real ) )GetPointer( &c[0] , BasisInfo< BasisType >::Coefficients ) , p.p );
	else            return _mesh->g( p.tIdx ).inverse() * RightTriangle< Real >::template EvaluateCovectorField< BasisType >( _mesh->g(p.tIdx) , ( ConstPointer( Real ) )GetPointer( &c[0] , BasisInfo< BasisType >::Coefficients ) , p.p );
}

/////////////////////////
// FEM::RiemannianMesh //
/////////////////////////
template< class Real , typename Index >
FEM::RiemannianMesh< Real , Index >::RiemannianMesh( Pointer( TriIndex ) t , size_t tC ) : _triangles(t) , _tCount(tC) , _edgeMap( t , tC )
{
	_g = AllocPointer< SquareMatrix< Real , 2 > >( _tCount );
	for( int t=0 ; t<_tCount ; t++ ) _g[t] = SquareMatrix< Real , 2 >::Identity();
	_vCount = 0;
	for( size_t i=0 ; i<_tCount ; i++ ) for( int j=0 ; j<3 ; j++ ) _vCount = std::max< size_t >( _vCount , _triangles[i][j]+1 );
}

template< class Real , typename Index >
FEM::RiemannianMesh< Real , Index >::~RiemannianMesh( void ){ FreePointer( _g ); }

template< class Real , typename Index >
template< unsigned int BasisType >
size_t FEM::RiemannianMesh< Real , Index >::dimension( void ) const
{
	size_t dim;
	switch( BasisInfo< BasisType >::ElementType )
	{
		case ELEMENT_VERTEX:   dim = vCount() * BasisInfo< BasisType >::CoefficientsPerElement ; break;
		case ELEMENT_EDGE:     dim = eCount() * BasisInfo< BasisType >::CoefficientsPerElement ; break;
		case ELEMENT_TRIANGLE: dim = tCount() * BasisInfo< BasisType >::CoefficientsPerElement ; break;
		default:
			TestBasisType( BasisType , "FEM::RiemannianMesh::dimension" , false );
			TestElementType( BasisInfo< BasisType >::ElementType , "FEM::RiemannianMesh::dimension" , true );
	}
	return dim;
}

template< class Real , typename Index >
template< unsigned int BasisType >
int FEM::RiemannianMesh< Real , Index >::index( int t , int idx ) const
{
	bool isAligned;
	return index< BasisType >( t , idx , isAligned );
}

template< class Real , typename Index >
template< unsigned int BasisType >
int FEM::RiemannianMesh< Real , Index >::index( int t , int idx , bool& isAligned ) const
{
	isAligned = true;
	int i , e = idx / BasisInfo< BasisType >::CoefficientsPerElement , c = idx % BasisInfo< BasisType >::CoefficientsPerElement;
	switch( BasisInfo< BasisType >::ElementType )
	{
		case ELEMENT_VERTEX:   i = _triangles[t][e]                   ; break;
		case ELEMENT_EDGE:     i = _edgeMap.edge( t*3+e , isAligned ) ; break;
		case ELEMENT_TRIANGLE: i = t                                  ; break;
		default: fprintf( stderr , "[ERROR] FEM::RiemannianMesh::index: unrecognized element type: %d\n" , BasisInfo< BasisType >::ElementType ) , exit( 0 );
	}
	return i * BasisInfo< BasisType >::CoefficientsPerElement + c;
}

template< class Real , typename Index >
Pointer( FEM::CoordinateXForm< Real > ) FEM::RiemannianMesh< Real , Index >::getCoordinateXForms( void ) const
{
	Pointer( CoordinateXForm< Real > ) xForms = NewPointer< CoordinateXForm< Real > >( _tCount*3 );
	setCoordinateXForms( xForms );
	return xForms;
}

template< class Real , typename Index >
FEM::CoordinateXForm< Real > FEM::RiemannianMesh< Real , Index >::xForm( int he ) const
{
	CoordinateXForm< Real > xForm;
	int ohe = _edgeMap.opposite(he);
	if( ohe==-1 ) MK_ERROR_OUT( "FEM::RiemannianMesh::xForm: Boundary edge" );

	// The two triangles on this edge
	int tIdx[] = { he/3 , ohe/3 };

	// The end-points of the edge
	int  _v =  he%3 ,  v[] = { ( he+1)%3 , ( he+2)%3 };
	int _ov = ohe%3 , ov[] = { (ohe+1)%3 , (ohe+2)%3 };

	// The direction of the edge
	TangentVector< Real > edgeDir = RightTriangle< Real >::EdgeDirections[_v] , oEdgeDir = -RightTriangle< Real >::EdgeDirections[_ov];
	// The perpendicular direction to the edge
	TangentVector< Real > edgePerpDir = Rotate90( _g[ tIdx[0] ] , edgeDir ) , oEdgePerpDir = Rotate90( _g[ tIdx[1] ] , oEdgeDir );
		
	// The linear part of the transformation should map ( edgeDir , edgePerpDir ) -> ( oEdgeDir , oEdgePerpDir )
	{
		SquareMatrix< Real , 2 > M , oM;
		M( 0 , 0 ) =     edgeDir[0] , M( 0 , 1 ) =     edgeDir[1];
		M( 1 , 0 ) = edgePerpDir[0] , M( 1 , 1 ) = edgePerpDir[1];
		oM( 0 , 0 ) =     oEdgeDir[0] , oM( 0 , 1 ) =     oEdgeDir[1];
		oM( 1 , 0 ) = oEdgePerpDir[0] , oM( 1 , 1 ) = oEdgePerpDir[1];
		xForm.linear = oM * M.inverse();
	}

	// The transformation should also take ( v[0] + v[1] ) / 2 to ( ov[0] + ov[1] ) / 2
	xForm.constant = RightTriangle< Real >::EdgeMidpoints[_ov] - xForm.linear * ( RightTriangle< Real >::EdgeMidpoints[_v] );
	return xForm;
}

template< class Real , typename Index >
void FEM::RiemannianMesh< Real , Index >::setCoordinateXForms( Pointer( CoordinateXForm< Real > ) xForms ) const
{
	ThreadPool::ParallelFor
		(
			0 , _edgeMap.size() ,
			[&]( unsigned int , size_t e )
			{
				const int* he = _edgeMap[(int)e];
				if( he[1]!=-1 )
				{
					xForms[ he[0] ] = xForm( he[0] );
					xForms[ he[1] ] = xForms[ he[0] ].inverse();
				}
				else xForms[ he[0] ] = CoordinateXForm< Real >();
			}
		);
}

template< class Real , typename Index >
void FEM::RiemannianMesh< Real , Index >::edgeVertices( int edge , int& v1 , int& v2 ) const
{
	int he = _edgeMap[edge][0];
	int  t =  he / 3 ,  v =  he % 3;
	v1 = _triangles[t][(v+1)%3] , v2 = _triangles[t][(v+2)%3];
}

template< class Real , typename Index >
bool FEM::RiemannianMesh< Real , Index >::oppositeEdgeVertices( int edge , int& v1 , int& v2 ) const
{
	int he = _edgeMap[edge][0] , ohe = _edgeMap[edge][1];
	int  t =  he / 3 ,  v =  he % 3;
	int ot = ohe / 3 , ov = ohe % 3;
	if( ohe==-1 ) return false;
	v1 = _triangles[t][v] , v2 = _triangles[ot][ov];
	return true;
}

template< class Real , typename Index >
bool FEM::RiemannianMesh< Real , Index >::edgeFlip( int edge , Real eps )
{
	int he = _edgeMap[edge][0] , ohe = _edgeMap[edge][1];
	int  t =  he / 3 ,  v =  he % 3;
	int ot = ohe / 3 , ov = ohe % 3;

	// First test that the edge is not on the boundary
	if( ohe==-1 ) return false;

	// Check that the two vertices we would like to link are not already linked
	// [WARNING] The call to getVertexCorners function will fail on a mesh with boundary
	{
		int neighbor = _triangles[ot][ov];
		std::vector< int > corners = getVertexCorners( t , v );
		for( int i=0 ; i<corners.size() ; i++ )
		{
			int _t = corners[i]/3;
			for( int j=0 ; j<3 ; j++ ) if( _triangles[_t][j]==neighbor ) return false;
		}
	}

	// Get the coordinates of the old and new edges and test that if we can flip
	Point2D< Real > newEdge[] = { RightTriangle< Real >::Corners[v] , xForm( ohe )( RightTriangle< Real >::Corners[ ov ] ) };
	Point2D< Real > oldEdge[] = { RightTriangle< Real >::Corners[(v+1)%3] , RightTriangle< Real >::Corners[(v+2)%3] };
	{

		// Solve for the point of intersection between the two corresponding line segments
		//  => Solve for s and t such that:
		//		newEdge[0] + s * ( newEdge[1]-newEdge[0]) = oldEdge[0] + t * ( oldEdge[1] - oldEdge[0] )
		// <=>	( newEdge[1]-newEdge[0] | -oldEdge[1]+oldEdge[0] ) * (s,t)^t = oldEdge[0] - newEdge[0]
		SquareMatrix< Real , 2 > M;
		for( int i=0 ; i<2 ; i++ ) M(0,i) = newEdge[1][i]-newEdge[0][i] , M(1,i) = -oldEdge[1][i]+oldEdge[0][i];
		Point2D< Real > st = M.inverse() * ( oldEdge[0] - newEdge[0] );

		// Test that the intersection point is in bounds
		if( st[0]<=eps || st[0]>=1-eps || st[1]<=eps || st[1]>=1-eps ) return false;
	}

	//     o              2,1
	//    / \             /|\
	//   /   \           / | \
	//  /     \         /  |  \
	// |-------|  =>  0| o | t |0
	//  \  e  /         \  |  /
	//   \   /           \ | /
	//    \ /             \|/
	//     v              1,2

	//  t*3 + ( v+1) -> ot*3+2
	//  t*3 + ( v+2) ->  t*3+1
	// ot*3 + (ov+1) ->  t*3+2
	// ot*3 + (ov+2) -> ot*3+1

	// The new triangles
	TriIndex tris[] = { TriIndex( _triangles[t][(v+1)%3] , _triangles[ot][ov] , _triangles[t][v] ) , TriIndex( _triangles[t][(v+2)%3] , _triangles[t][v] , _triangles[ot][ov] ) };

	// The new metrics
	SquareMatrix< Real , 2 > tensors[2];
	tensors[0](0,0) = Point2D< Real >::Dot( RightTriangle< Real >::Corners[ov] - RightTriangle< Real >::Corners[(ov+2)%3] , _g[ot] * ( RightTriangle< Real >::Corners[ov] - RightTriangle< Real >::Corners[(ov+2)%3] ) );
	tensors[0](1,1) = Point2D< Real >::Dot( RightTriangle< Real >::Corners[ v] - RightTriangle< Real >::Corners[( v+1)%3] , _g[ t] * ( RightTriangle< Real >::Corners[ v] - RightTriangle< Real >::Corners[( v+1)%3] ) );
	tensors[0](1,0) = tensors[0](0,1) = ( tensors[0](0,0) + tensors[0](1,1) - Point2D< Real >::Dot( newEdge[1]-newEdge[0] , _g[t] * ( newEdge[1]-newEdge[0] ) ) ) / (Real)2.;

	tensors[1](0,0) = Point2D< Real >::Dot( RightTriangle< Real >::Corners[ v] - RightTriangle< Real >::Corners[( v+2)%3] , _g[ t] * ( RightTriangle< Real >::Corners[ v] - RightTriangle< Real >::Corners[( v+2)%3] ) );
	tensors[1](1,1) = Point2D< Real >::Dot( RightTriangle< Real >::Corners[ov] - RightTriangle< Real >::Corners[(ov+1)%3] , _g[ot] * ( RightTriangle< Real >::Corners[ov] - RightTriangle< Real >::Corners[(ov+1)%3] ) );
	tensors[1](1,0) = tensors[1](0,1) = ( tensors[1](0,0) + tensors[1](1,1) - Point2D< Real >::Dot( newEdge[1]-newEdge[0] , _g[t] * ( newEdge[1]-newEdge[0] ) ) ) / (Real)2.;

	int old_he[] = { t*3+(v+1)%3 , ot*3+(ov+1)%3 , t*3+(v+2)%3 , ot*3+(ov+2)%3 } , new_he[] = { 3*ot+2 , 3*t+2 , 3*t+1 , 3*ot+1 };
	int edges[4];
	for( int i=0 ; i<4 ; i++ ) edges[i] = _edgeMap.edge( old_he[i] );
	for( int i=0 ; i<4 ; i++ )
	{
		if     ( _edgeMap._e2he[2*edges[i]+0]==old_he[i] ) _edgeMap._e2he[2*edges[i]+0] = new_he[i] , _edgeMap._he2e[ new_he[i] ] =  edges[i]+1;
		else if( _edgeMap._e2he[2*edges[i]+1]==old_he[i] ) _edgeMap._e2he[2*edges[i]+1] = new_he[i] , _edgeMap._he2e[ new_he[i] ] = -edges[i]-1;
		else fprintf( stderr , "[ERROR] FEM::RiemannianMesh::edgeFlip: Unmatched half-edge\n" ) , exit( 0 );
	}
	_edgeMap._e2he[2*edge+0] = 3* t , _edgeMap._he2e[3* t] =  edge+1;
	_edgeMap._e2he[2*edge+1] = 3*ot , _edgeMap._he2e[3*ot] = -edge-1;

	// Set the triangles and metrics
	_triangles[t] = tris[0] , _triangles[ot] = tris[1];
	_g[t] = tensors[0] , _g[ot] = tensors[1];
	return true;
}

template< class Real , typename Index >
void FEM::RiemannianMesh< Real , Index >::sanityCheck( ConstPointer( CoordinateXForm< Real > ) xForms , Real eps ) const
{
#if 1
	MK_WARN_ONCE( "Method unsupported" );
#else
	static const Point2D< Real > Corners[] = { Point2D< Real >(0,0) , Point2D< Real >(1,0) , Point2D< Real >(0,1) };
	for( int t=0 ; t<_tCount ; t++ )
	{
		bool success = true;

		// Check that the metric tensor is symmetric
		if( fabs( _g[t](1,0) - _g[t](0,1) )>eps ) fprintf( stderr , "[ERROR] Metric tensor is not symetric: %g != %g\n" , _g[t](1,0) , _g[t](0,1) ) , success = false;

		// The characteristic polynomial of the metric tensor is:
		// P(x) = x^2 - Tr * x + det
		double a = 1. , b = -_g[t].trace() , c = _g[t].determinant();
		double disc = b*b - 4.*a*c;
		if( disc<0 ) fprintf( stderr , "[ERROR] Vanishing discriminant: %g\n" , disc ) , success = false;
		else
		{
			double x[] = { ( -b - sqrt( disc ) ) / (2.*a) , ( -b + sqrt( disc ) ) / (2.*a) };
			if( x[0]<=0 ) fprintf( stderr , "[ERROR] Metric tensor is not positive definite: %g %g\n" , x[0] , x[1] ) , success = false;
		}

		if( !success ) fprintf( stderr , "[ERROR] Failure on triangle: %d\n" , t ) , exit( 0 );
	}
	for( int e=0 ; e<_tCount*3 ; e++ ) if( _edgeMap.opposite( e )!=-1 )
	{
		bool success = true;
		int oe = _edgeMap.opposite( e );
		int t = e/3 , v = e%3 , ot = oe/3 , ov = oe%3;

		// Test the opposite of the opposite is itself
		if( e!=_edgeMap.opposite( oe ) ) fprintf( stderr , "[ERROR] edge is not the opposite of its opposite: %d != %d\n" , e , _edgeMap.opposite( oe ) ) , success = false;

		// Test that the shared vertices are the same (assumes oriented)
		if( _triangles[t][(v+1)%3]!=_triangles[ot][(ov+2)%3] || _triangles[t][(v+2)%3]!=_triangles[ot][(ov+1)%3] )
			fprintf( stderr , "[ERROR] Vertices don't match: %d %d != %d %d\n" , _triangles[t][(v+1)%3] , _triangles[t][(v+2)%3] , _triangles[ot][(ov+2)%3] , _triangles[ot][(ov+1)%3] ) , success = false;

		// Test that transforming across the edge and then back gives the identity
		{
			FEM::CoordinateXForm< Real > xForm = edges[oe] * edges[e];
			Real cError = ( xForm.constant ).squareNorm() , lError = ( xForm.linear - SquareMatrix< Real , 2 >::Identity() ).squareNorm();
			if( cError > eps*eps || lError > eps*eps ) fprintf( stderr , "[ERROR] edge transformations are not inverses: %g %g\n" , cError , lError ) , success = false;
		}

		// Test that the lengths of the shared edge agree
		{
			Real l1 = (Real)sqrt( Point2D< Real >::Dot( Corners[( v+1)%3] - Corners[( v+2)%3] , _g[ t] * ( Corners[( v+1)%3] - Corners[( v+2)%3] ) ) );
			Real l2 = (Real)sqrt( Point2D< Real >::Dot( Corners[(ov+1)%3] - Corners[(ov+2)%3] , _g[ot] * ( Corners[(ov+1)%3] - Corners[(ov+2)%3] ) ) );
			if( fabs( l1-l2 )>eps ) fprintf( stderr , "[ERROR] edge lengths don't agree: %g != %g\n" , l1 , l2 ) , success = false;
		}

		// Test that the edges map into each other
		{
			Point2D< Real > e1 = Corners[(v+1)%3] - Corners[(v+2)%3];
			Point2D< Real > e2 = -edges[oe].linear * ( Corners[(ov+1)%3] - Corners[(ov+2)%3]  );
			if( Point2D< Real >::Dot( e1-e2 , _g[t]*(e1-e2) )>eps*eps ) fprintf( stderr , "[ERROR] edges don't match: %g %g %g\n" , e1[0]-e2[0] , e1[1]-e2[1] , e1[2]-e2[2] ) , success = false;
		}

		// Test that the perpendiculars map into each other
		{
			Point2D< Real > p1 = Rotate90( _g[t] , Corners[(v+1)%3] - Corners[(v+2)%3] );
			Point2D< Real > p2 = -edges[oe].linear * Rotate90( _g[ot] , Corners[(ov+1)%3] - Corners[(ov+2)%3] );
			if( Point2D< Real >::Dot( p1-p2 , _g[t]*(p1-p2) )>eps*eps ) fprintf( stderr , "[ERROR] edge perps don't match: %g %g %g\n" , p1[0]-p2[0] , p1[1]-p2[1] , p1[2]-p2[2] ) , success = false;
		}

		if( !success ) fprintf( stderr , "[ERROR] Failure on edge between triangles: %d %d\n" , t , ot ) , exit( 0 );
	}
#endif
}

template< class Real , typename Index >
bool FEM::RiemannianMesh< Real , Index >::isVoronoiEdge( int e , Real eps ) const
{
	int he = _edgeMap[e][0] , ohe = _edgeMap[e][1];
	int  t =  he / 3 ,  v =  he % 3;
	int ot = ohe / 3 , ov = ohe % 3;
	if( ohe==-1 ) return true;
	Point2D< Real > center = RightTriangle< Real >::Center( _g[t] , RightTriangle< Real >::CENTER_CIRCUMCENTRIC ) , oVertex = xForm( ohe )( RightTriangle< Real >::Corners[ov] );
	return Point2D< Real >::Dot( center - oVertex , _g[t] * ( center - oVertex ) )+eps > Point2D< Real >::Dot( center - RightTriangle< Real >::Corners[0] , _g[t] * ( center - RightTriangle< Real >::Corners[0] ) );
}

template< class Real , typename Index >
std::vector< FEM::SamplePoint< Real > > FEM::RiemannianMesh< Real , Index >::randomSamples( unsigned int count ) const
{
	std::vector< FEM::SamplePoint< Real > > samples( count );
	Real* cumAreas = new Real[ tCount() ];
	cumAreas[0] = area(0);
	for( int i=1 ; i<tCount() ; i++ ) cumAreas[i] = cumAreas[i-1] + area(i);

	ThreadPool::ParallelFor
		(
			0 , count ,
			[&]( unsigned int  , size_t i )
			{
				Real r1 = Random< Real >() , r2 = Random< Real >() , r3 = Random< Real >();

				// Choose a random triangle
				{
					Real r = r1 * cumAreas[ tCount()-1 ];
					Real *a , *b;
					a = cumAreas - 1;
					b = cumAreas + tCount() - 1;

					int d;
					while( (d=(int)(b-a))>1 )
					{
						Real *i = a + d/2;
						if( r<*i ) b=i;
						else       a=i;
					}
					samples[i].tIdx = (int)( b - cumAreas );
				}
				// Choose a random point in the triangle
				{
					if( r2+r3>1 ) r2 = 1-r2 , r3 = 1-r3;
					samples[i].p = Point2D< Real >( (Real)r2 , (Real) r3 );
				}
			}
		);
	delete[] cumAreas;

	return samples;
}

template< class Real , typename Index >
FEM::CoordinateXForm< Real > FEM::RiemannianMesh< Real , Index >::getVertexCoordinateXForm( ConstPointer( CoordinateXForm< Real > ) xForms , int t , int v ) const
{
	const int VertexToEdgeMap[] = { 1 , 2 , 0 };
	const int EdgeToVertexMap[] = { 1 , 2 , 0 };
	// Assume that the mesh is oriented
	FEM::CoordinateXForm< Real > xForm;
	int currentT = t , currentV = v;
	do
	{
		int edge = currentT*3 + VertexToEdgeMap[ currentV ] , oEdge = _edgeMap.opposite( edge );
		xForm = xForms[ edge ] * xForm;
		if( oEdge==-1 ) fprintf( stderr , "[ERROR] Boundary vertex\n" ) , exit( 0 );
		currentT = oEdge / 3;
		currentV = EdgeToVertexMap[ oEdge%3 ];
	}
	while( currentT!=t );
	return xForm;
}

template< class Real , typename Index >
std::vector< int > FEM::RiemannianMesh< Real , Index >::getVertexCorners( int t , int v ) const
{
#if 1
	auto AddVertexCorners = [&]( int t , int v , const int *Map , std::vector< int > &neighbors )
	{
		int currentT = t , currentV = v;
		do
		{
			neighbors.push_back( currentT*3 + currentV );

			int he = currentT*3 + Map[ currentV ] , ohe = _edgeMap.opposite( he );
			if( ohe==-1 ) return false;
			currentT = ohe / 3;
			currentV = Map[ ohe%3 ];
		}
		while( currentT!=t );
		return true;
	};

	const int CCWMap[] = { 1 , 2 , 0 }; // Vertex to incoming edge (CCW)
	const int  CWMap[] = { 2 , 0 , 1 }; // Vertex to outgoing edge (CW)
	std::vector< int > neighbors;
	if( !AddVertexCorners( t , v , CCWMap , neighbors ) )
	{
		std::vector< int > _neighbors , __neighbors;
		AddVertexCorners( t , v , CWMap , _neighbors );
		__neighbors.reserve( neighbors.size() + _neighbors.size()-1 );
		for( int i=(int)_neighbors.size()-1 ; i!=0 ; i-- ) __neighbors.push_back( _neighbors[i] );
		for( int i=0 ; i<(int)neighbors.size() ; i++ ) __neighbors.push_back( neighbors[i] );
		return __neighbors;
	}
	else return neighbors;
#else
	// Circulate CCW
	const int VertexToEdgeMap[] = { 1 , 2 , 0 };
	const int EdgeToVertexMap[] = { 1 , 2 , 0 };
	// Assume that the mesh is oriented
	std::vector< int > neighbors;
	int currentT = t , currentV = v;
	do
	{
		int he = currentT*3 + VertexToEdgeMap[ currentV ] , ohe = _edgeMap.opposite( he );
		neighbors.push_back( currentT*3 + currentV );
		if( ohe==-1 ) fprintf( stderr , "[ERROR] Boundary vertex\n" ) , exit( 0 );
		currentT = ohe / 3;
		currentV = EdgeToVertexMap[ ohe%3 ];
	}
	while( currentT!=t );
	return neighbors;
#endif
}

template< class Real , typename Index >
template< typename CornerFunctor >
void FEM::RiemannianMesh< Real , Index >::processCorners( int t , int v , CornerFunctor &F ) const
{
	// Assume that the mesh is oriented and that we are not at a boundary vertex
	const int VertexToIncomingEdge[] = { 1 , 2 , 0 };
	const int EdgeToStartVertex[] = { 1 , 2 , 0 };
	int currentT = t , currentV = v;
	do
	{
		F( currentT , currentV );
		int he = currentT*3 + VertexToIncomingEdge[ currentV ] , ohe = _edgeMap.opposite( he );
		if( ohe==-1 ) fprintf( stderr , "[ERROR] Boundary vertex\n" ) , exit( 0 );
		currentT = ohe / 3;
		currentV = EdgeToStartVertex[ ohe%3 ];
	}
	while( currentT!=t );
}

template< class Real , typename Index >
Real FEM::RiemannianMesh< Real , Index >::getVertexConeAngle( int t , int v ) const
{
	Real angle = (Real)0;

	auto CornerFunctor = [&]( int _t , int _v ){ angle += RightTriangle< Real >::Angle( _g[_t] , _v ); };
	processOneRing( t , v , CornerFunctor );
	return angle;
}

template< class Real , typename Index >
FEM::CoordinateXForm< Real > FEM::RiemannianMesh< Real , Index >::exp( ConstPointer( CoordinateXForm< Real > ) xForms , HermiteSamplePoint< Real >& p , Real eps , bool noWarning ) const
{
	HermiteSamplePoint< Real > startP = p;
	CoordinateXForm< Real > xForm = CoordinateXForm< Real >::Identity();
	if( !Point2D< Real >::SquareNorm( p.v ) ) return xForm;
	const int MAX_ITERS = 10000;
	int count = 0;
	int inEdge = -1;
#if 1
	// If the starting point happens to be on an edge
	{
		int idx = -1;
		if     ( p.p[0]<=0 && p.v[0]<0 ) idx = 1;
		else if( p.p[1]<=0 && p.v[1]<0 ) idx = 2;
		else if( p.p[0]+p.p[1]>=1 && p.v[0]+p.v[1]>0 ) idx = 0;
		if( idx!=-1 )
		{
			int he = p.tIdx * 3 + idx , ohe = _edgeMap.opposite( he );
			if( ohe==-1 )
			{
				MK_WARN_ONCE( "Hit boundary" );
				return xForm;
			}

			const CoordinateXForm< Real >& edge = xForms[he];
			p.tIdx = ohe/3;
			p.p = edge( p.p ) ; p.v = edge.linear * p.v;
			inEdge = ohe%3;
			xForm = edge * xForm;
		}
	}
#endif
	while( count<MAX_ITERS )
	{
		if( !Point2D< Real >::SquareNorm( p.v ) ) return xForm;

		// Intersect the ray p + s*v with each of the three edges
		// Bottom edge:   p[1] + s * v[1] = 0                         => s = -p[1]/v[1]
		// Left edge:     p[0] + s * v[0] = 0                         => s = -p[0]/v[0]
		// Diagonal edge: p[1] + s * v[1] = 1 - ( p[0] + s * v[0] )   => s = ( 1 - p[0]  - p[1] ) / ( v[1] + v[0] )
		Real s = 0;

		int idx = -1;
		{

			Real _s[] = { -p.p[1] / p.v[1]  , -p.p[0] / p.v[0] , ( Real(1.) - p.p[0]  - p.p[1] ) / ( p.v[1] + p.v[0] ) };
			s = std::numeric_limits< Real >::infinity();
			if( inEdge!=2 && _s[0]>0 && _s[0]<s ) idx = 2 , s = _s[0];
			if( inEdge!=1 && _s[1]>0 && _s[1]<s ) idx = 1 , s = _s[1];
			if( inEdge!=0 && _s[2]>0 && _s[2]<s ) idx = 0 , s = _s[2];
		}
		if( idx==-1 ) return xForm;

		if( s>1 ) // The end-point is within the triangle
		{
			p.p += p.v , p.v -= p.v;
			return xForm;
		}
		else // The end-point is outside the triangle
		{
			// Advance along the flow until you hit the edge	
			int he = p.tIdx*3 + idx , ohe = _edgeMap.opposite( he );
			p.p += p.v*s ; p.v -= p.v*s;
			if( opposite( he )==-1 )
			{
				MK_WARN_ONCE( "Hit boundary" );
				return xForm;
			}

			const CoordinateXForm< Real >& edge = xForms[he];
			p.tIdx = ohe/3 ; p.p = edge( p.p ) ; p.v = edge.linear * p.v;
			inEdge = ohe%3;
			xForm = edge * xForm;
		}
		count++;
	}
	if( !noWarning ) MK_WARN( "Failed to converge exp after " , MAX_ITERS , " iterations" );
	return xForm;
}

template< class Real , typename Index >
FEM::CoordinateXForm< Real > FEM::RiemannianMesh< Real , Index >::flow( ConstPointer( CoordinateXForm< Real > ) xForms , const TangentVectorField< Real >& vf , Real flowTime , SamplePoint< Real >& p , Real minStepSize , Real eps , std::vector< SamplePoint< Real > >* path , bool noWarning ) const
{
	CoordinateXForm< Real > xForm = CoordinateXForm< Real >::Identity();
	int MAX_ITERS = 1000000;
	int count = 0;
	int inEdge = -1;
	Real direction = (flowTime<0) ? (Real)-1. : (Real)1.;
	Real stepSizeLeft = minStepSize;
	TangentVector< Real > v = vf( p ) * direction;
	Real vLength = (Real)sqrt( Point2D< Real >::Dot( v , _g[p.tIdx]*v ) );
	flowTime *= direction;
	if( path ) path->push_back( p );
	while( count<MAX_ITERS )
	{
		if( !Point2D< Real >::SquareNorm( v ) ) return xForm;
		// Intersect the ray p + s * v with each of the three edges
		// Bottom edge:   p[1] + s * v[1] = 0                         => s = -p[1]/v[1]
		// Left edge:     p[0] + s * v[0] = 0                         => s = -p[0]/v[0]
		// Diagonal edge: p[1] + s * v[1] = 1 - ( p[0] + s * v[0] )   => s = ( 1 - p[0]  - p[1] ) / ( v[1] + v[0] )
		Real s = 0;
		int idx = -1;
		{
			Real _s[] = { -p.p[1] / v[1]  , -p.p[0] / v[0] , ( Real(1.) - p.p[0]  - p.p[1] ) / ( v[1] + v[0] ) };
			s = 1e10;
			if( inEdge!=2 && _s[0]>0 && _s[0]<s ) idx = 2 , s = _s[0];
			if( inEdge!=1 && _s[1]>0 && _s[1]<s ) idx = 1 , s = _s[1];
			if( inEdge!=0 && _s[2]>0 && _s[2]<s ) idx = 0 , s = _s[2];
		}
		if( idx==-1 ) return xForm;
		CotangentVector< Real > gv = _g[p.tIdx] * v;
		Real stepSize = vLength * s;
		bool updateVector = false;
		if( minStepSize>0 && stepSize>stepSizeLeft )
		{
			s = stepSizeLeft / vLength;
			updateVector = true;
		}

		// If we can finish the flow
		if( flowTime<s )
		{
			p.p += v * flowTime;
			if( path ) path->push_back( p );
			return xForm;
		}
		else
		{
			p.p += v * s , flowTime -= s;

			// If we do not cross a boundary, change direction
			if( updateVector )
			{
				// If the the vectors are oppositely oriented, terminate the flow
				if( Point2D< Real >::Dot( gv , vf( p ) )*direction < 0 ) return xForm;

				v = vf( p ) * direction;
				vLength = (Real)sqrt( Point2D< Real >::Dot( v , _g[p.tIdx]*v ) );
				stepSizeLeft = minStepSize;
				inEdge = -1;
			}
			// If we cross the boundary, transport the direction of the flow
			else // The end-point is outside the triangle
			{
				// Advance along the flow until you hit the edge	
				int he = p.tIdx*3 + idx , ohe = _edgeMap.opposite( he );
				if( opposite( he )==-1 )
				{
					MK_WARN_ONCE( "Hit boundary" );
					return xForm;
				}

				const CoordinateXForm< Real >& edge = xForms[he];
				// Switch into the next triangle
				p.tIdx = ohe/3 ; p.p = edge( p.p ) ; v = edge.linear * v;

				// Mark the edge we came in on
				inEdge = ohe%3;

				// Accumulate the transformations
				xForm = edge * xForm;

				stepSizeLeft -= stepSize;
			}

			if( path ) path->push_back( p );
		}
		count++;
	}
	if( !noWarning ) fprintf( stderr , "[WARNING] Failed to converge flow after %d iterations\n" , MAX_ITERS );
	return xForm;
}

/////////////////////////
// FEM::RiemannianMesh //
/////////////////////////
template< class Real , typename Index >
inline void FEM::RiemannianMesh< Real , Index >::makeArea( Real area )
{
	double scale = 0;
	std::vector< double > _scales( ThreadPool::NumThreads() , 0 );
	ThreadPool::ParallelFor( 0 , _tCount , [&]( unsigned int t , size_t i ){ _scales[t] += sqrt( _g[i].determinant() );} );
	for( unsigned int i=0 ; i<_scales.size() ; i++ ) scale += _scales[i];
	scale = 2. / scale;
	scale *= area;
	ThreadPool::ParallelFor( 0 , _tCount , [&]( unsigned int , size_t i ){ _g[i] *= (Real)scale; } );
}

template< class Real , typename Index >
inline Real FEM::RiemannianMesh< Real , Index >::area( void ) const
{
	Real area = 0;
	std::vector< Real > _areas( ThreadPool::NumThreads() , 0 );
	ThreadPool::ParallelFor( 0 , _tCount , [&]( unsigned int t , size_t i ){ _areas[t] += (Real)sqrt( _g[i].determinant() ); } );
	for( unsigned int t=0 ; t<_areas.size() ; t++ ) area += _areas[t];
	return area / (Real)2.;
}

template< class Real , typename Index >
inline Real FEM::RiemannianMesh< Real , Index >::area( int idx ) const { return (Real)sqrt( _g[idx].determinant() ) / (Real)2.; }

template< class Real , typename Index >
inline Real FEM::RiemannianMesh< Real , Index >::squareEdgeLength( int heIdx ) const { return SquareLength( _g[heIdx/3] , RightTriangle< Real >::EdgeDirections[heIdx%3] ); }

#if 1
template< class Real , typename Index >
template< unsigned int Dim >
unsigned int FEM::RiemannianMesh< Real , Index >::setMetricFromEmbedding( ConstPointer( Point< Real , Dim > ) pointList , bool warn ){ return setMetricFromEmbedding< Dim >( [&]( unsigned int i ){ return pointList[i]; } , warn ); }

template< class Real , typename Index >
template< unsigned int Dim >
unsigned int FEM::RiemannianMesh< Real , Index >::setMetricFromEmbedding( std::function< Point< Real , Dim > (unsigned int) > PointList , bool warn )
{
	std::mutex mut;
	unsigned int badDeterminantCount = 0;
	ThreadPool::ParallelFor
		(
			0 , _tCount ,
			[&]( unsigned int , size_t i )
			{
				Point< Real , Dim > e[] = { PointList( _triangles[i][1] ) - PointList( _triangles[i][0] ) , PointList( _triangles[i][2] ) - PointList( _triangles[i][0] ) };
				for( int j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ ) _g[i](j,k) = Point< Real , Dim >::Dot( e[j] , e[k] );
				_g[i](0,1) = _g[i](1,0) = (Real)( ( _g[i](0,1) + _g[i](1,0) )/2 );

				if( !_g[i].determinant() )
				{
					std::lock_guard< std::mutex > lock( mut );
					if( warn )
					{
						fprintf( stderr , "[WARNING] Vanishing metric tensor determinant[%d]\n" , (int)i );
						for( int j=0 ; j<3 ; j++ ) std::cerr << "\t[" << _triangles[i][j] << "] = " << PointList( _triangles[i][j] ) << std::endl;
						std::cerr << "\t" << e[0] << "\t" << e[1] << std::endl;
					}
					badDeterminantCount++;
				}
			}
		);
	return badDeterminantCount;
}
#else
template< class Real , typename Index >
template< class Vertex >
void FEM::RiemannianMesh< Real , Index >::setMetricFromEmbedding( ConstPointer( Vertex ) vertices )
{

	std::mutex mut;
	ThreadPool::ParallelFor
		(
			0 , _tCount ,
			[&]( unsigned int , size_t i )
			{
				Point3D< Real > e[] = { Point3D< Real >( vertices[ _triangles[i][1] ].data< 0 >() ) - Point3D< Real >( vertices[ _triangles[i][0] ].data< 0 >() ) , Point3D< Real >( vertices[ _triangles[i][2] ].data< 0 >() ) - Point3D< Real >( vertices[ _triangles[i][0] ].data< 0 >() ) };
				for( int j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ ) _g[i](j,k) = Point3D< Real >::Dot( e[j] , e[k] );
				_g[i](0,1) = _g[i](1,0) = (Real)( ( _g[i](0,1) + _g[i](1,0) )/2 );

				if( !_g[i].determinant() )
				{
					std::lock_guard< std::mutex > lock( mut );
					fprintf( stderr , "[WARNING] Vanishing metric tensor determinant[%d]\n" , i );
					for( int j=0 ; j<3 ; j++ ) fprintf( stderr , "\tv[%d] = %.6f %.6f %.6f\n" , _triangles[i][j] , vertices[ _triangles[i][j] ].data<0>()[0] , vertices[ _triangles[i][j] ].data<0>()[1] , vertices[ _triangles[i][j] ].data<0>()[2] );
					printf( "\t%g %g %g\t%g %g %g\n" , e[0][0] , e[0][1] , e[0][2] , e[1][0] , e[1][1] , e[1][2] );
				}
			}
		);
}
#endif

template< class Real , typename Index >
void FEM::RiemannianMesh< Real , Index >::setMetricFromEdgeLengths( ConstPointer( Real ) edgeLengths )
{
	// The conditions the matrix needs to satisfy are:
	// -- < g[i] * (1,0) , (1,0) > = edgeLengths[i*3+2]^2
	//		=> g[i](0,0) = edgeLengths[i*3+2]^2
	// -- < g[i] * (0,1) , (0,1) > = edgeLengths[i*3+1]^2
	//		=> g[i](1,1) = edgeLengths[i*3+1]^2
	// -- < g[i] * (-1,1) , (-1,1) > = edgeLengths[i*3+0]^2
	//		=> g[i](0,0) + g[i](1,1) - g[i](0,1) - g[i](1,0) = edgeLengths[i*3+0]^2
	//		=> - g[i](0,1) - g[i](1,0) = edgeLengths[i*3+0]^2 - edgeLengths[i*3+2]^2 - edgeLengths[i*3+1]^2
	//		=>  g[i](0,1) = g[i](1,0) = ( edgeLengths[i*3+2]^2 + edgeLengths[i*3+1]^2 - edgeLengths[i*3+0]^2 ) / 2

	ThreadPool::ParallelFor
		(
			0 , _tCount ,
			[&]( unsigned int , size_t i )
			{
				_g[i](0,0) = edgeLengths[i*3+2] * edgeLengths[i*3+2];
				_g[i](1,1) = edgeLengths[i*3+1] * edgeLengths[i*3+1];
				_g[i](0,1) = _g[i](1,0) = ( _g[i](0,0) + _g[i](1,1) - edgeLengths[i*3] * edgeLengths[i*3] ) / (Real)2.;
			}
		);
}

template< class Real , typename Index >
void FEM::RiemannianMesh< Real , Index >::setMetricFromSquareEdgeLengths( ConstPointer( Real ) squareEdgeLengths )
{
	ThreadPool::ParallelFor
		(
			0 , _tCount ,
			[&]( unsigned int , size_t i )
			{
				_g[i](0,0) = squareEdgeLengths[i*3+2];
				_g[i](1,1) = squareEdgeLengths[i*3+1];
				_g[i](0,1) = _g[i](1,0) = ( _g[i](0,0) + _g[i](1,1) - squareEdgeLengths[i*3] ) / (Real)2.;
			}
		);
}

template< class Real , typename Index >
template< unsigned int BasisType , class V >
V FEM::RiemannianMesh< Real , Index >::evaluateScalarField( ConstPointer( V ) coefficients , const SamplePoint< Real >& p ) const
{
	V v;
	TestBasisType( BasisType , "FEM::RiemannianMesh::evaluateScalarField" , false );
	if( BasisInfo< BasisType >::Dimension==0 )
	{
		V _coefficients[ BasisInfo< BasisType >::Coefficients ];
		for( int i=0 ; i<BasisInfo< BasisType >::Coefficients ; i++ )
		{
			bool isAligned;
			int ii = index< BasisType >( p.tIdx , i , isAligned );
			_coefficients[i] = isAligned ? coefficients[ii] : -coefficients[ii];
		}
		v = RightTriangle< Real >::template EvaluateScalarField< BasisType >( ( ConstPointer( V ) )GetPointer( _coefficients , BasisInfo< BasisType >::Coefficients ) , p.p );
	}
	else TestBasisType( BasisType , "FEM::RiemannianMesh::evaluateScalarField" , false );
	return v;
}

template< class Real , typename Index >
template< unsigned int BasisType , class V >
FEM::CotangentVector< V > FEM::RiemannianMesh< Real , Index >::evaluateScalarFieldDifferential( ConstPointer( V ) coefficients , const SamplePoint< Real >& p ) const
{
	CotangentVector< V > v;
	TestBasisType( BasisType , "FEM::RiemannianMesh::evaluateScalarField" , false );
	if( BasisInfo< BasisType >::Dimension==0 )
	{
		Real _coefficients[ BasisInfo< BasisType >::Coefficients ];
		for( int i=0 ; i<BasisInfo< BasisType >::Coefficients ; i++ )
		{
			bool isAligned;
			int ii = index< BasisType >( p.tIdx , i , isAligned );
			_coefficients[i] = isAligned ? coefficients[ii] : -coefficients[ii];
		}
		v = RightTriangle< Real >::template EvaluateScalarFieldDifferential< BasisType >( ( ConstPointer( Real ) )GetPointer( _coefficients , BasisInfo< BasisType >::Coefficients ) , p.p );
	}
	else TestBasisType( BasisType , "FEM::RiemannianMesh::evaluateScalarFieldDifferential" , false );
	return v;
}

template< class Real , typename Index >
template< unsigned int BasisType , class V >
FEM::TangentVector< V > FEM::RiemannianMesh< Real , Index >::evaluateScalarFieldGradient( ConstPointer( V ) coefficients , const SamplePoint< Real >& p ) const
{
	return Dual( _g[p.tIdx] , evaluateScalarFieldDifferential( coefficients , p ) );
}

template< class Real , typename Index >
template< unsigned int BasisType , class V >
FEM::CotangentVector< V > FEM::RiemannianMesh< Real , Index >::evaluateCovectorField( ConstPointer( V ) coefficients , const SamplePoint< Real >& p ) const
{
	CotangentVector< V > v;
	TestBasisType( BasisType , "FEM::RiemannianMesh::evaluateCovectorField" , false );
	if( BasisInfo< BasisType >::Dimension==1 )
	{
		Real _coefficients[ BasisInfo< BasisType >::Coefficients ];
		for( int i=0 ; i<BasisInfo< BasisType >::Coefficients ; i++ )
		{
			bool isAligned;
			int ii = index< BasisType >( p.tIdx , i , isAligned );
			_coefficients[i] = isAligned ? coefficients[ii] : -coefficients[ii];
		}
		v = RightTriangle< Real >::template EvaluateCovectorField< BasisType >( GetPointer( _coefficients , BasisInfo< BasisType >::Coefficients ) , p.p );
	}
	else TestBasisType( BasisType , "FEM::RiemannianMesh::evaluateScalarField" , false );
	return v;
}

template< class Real , typename Index >
#ifdef EIGEN_WORLD_VERSION 
template< unsigned int BasisType , bool UseEigen >
std::conditional_t< UseEigen , Eigen::SparseMatrix< Real > , SparseMatrix< Real , int > > FEM::RiemannianMesh< Real , Index >::massMatrix( bool lump , ConstPointer( SquareMatrix< Real , 2 > ) newTensors ) const
#else // !EIGEN_WORLD_VERSION 
template< unsigned int BasisType >
SparseMatrix< Real , int > FEM::RiemannianMesh< Real , Index >::massMatrix( bool lump , ConstPointer( SquareMatrix< Real , 2 > ) newTensors ) const
#endif // EIGEN_WORLD_VERSION 
{
#ifdef EIGEN_WORLD_VERSION 
	std::conditional_t< UseEigen , Eigen::SparseMatrix< Real > , SparseMatrix< Real , int > > M;
	std::vector< Eigen::Triplet< Real > > triplets;
#else // !EIGEN_WORLD_VERSION 
	SparseMatrix< Real , int > M;
#endif // EIGEN_WORLD_VERSION
	if( BasisType==BASIS_2_VERTEX_CONSTANT )
	{
#ifdef EIGEN_WORLD_VERSION
		if constexpr( UseEigen )
		{
			M.resize( _vCount , _vCount );
			triplets.resize( _tCount*3 );

			ThreadPool::ParallelFor
			(
				0 , _tCount ,
				[&]( unsigned int , size_t i )
				{
					Real a = area( (int)i ) / (Real)3.;
					if( newTensors ) a *= newTensors[i].determinant() / _g[i].determinant() / _g[i].determinant();
					for( int j=0 ; j<3 ; j++ ) triplets[i*3+j] = Eigen::Triplet< Real >( _triangles[i][j] , _triangles[i][j] , a );
				}
			);
			M.setFromTriplets( triplets.begin() , triplets.end() );
		}
		else
#endif // EIGEN_WORLD_VERSION
		{
			M.resize( _vCount );
			ThreadPool::ParallelFor( 0 , _vCount , [&]( unsigned int , size_t i ){ M.SetRowSize( i , 1 ) , M[i][0] = MatrixEntry< Real , int >( (int)i , 0 ); } );
			ThreadPool::ParallelFor
			(
				0 , _tCount ,
				[&]( unsigned int , size_t i )
				{
					Real a = area( (int)i ) / (Real)3.;
					if( newTensors ) a *= newTensors[i].determinant() / _g[i].determinant() / _g[i].determinant();
					for( int j=0 ; j<3 ; j++ ) AddAtomic( M[ _triangles[i][j] ][0].Value , a );
				}
			);
			ThreadPool::ParallelFor( 0 , M.rows , [&]( unsigned int , size_t i ){ M[i][0].Value = (Real)( 1. / M[i][0].Value ); } );
		}
		return M;
	}
	auto mask = RightTriangle< Real >::template GetMassMask< BasisType >( newTensors!=( ConstPointer( SquareMatrix< Real , 2 > ) )NullPointer< SquareMatrix< Real , 2 > >() );
	Point< int , BasisInfo< BasisType >::Coefficients > nonZeroCount;
	for( int i=0 ; i<BasisInfo< BasisType >::Coefficients ; i++ ) for( int j=0 ; j<BasisInfo< BasisType >::Coefficients ; j++ ) if( mask(i,j) ) nonZeroCount[i]++;

	Pointer( std::atomic< int > ) rowSizes = NullPointer< std::atomic< int > >();
#ifdef EIGEN_WORLD_VERSION
	if constexpr( UseEigen )
	{
		M.resize( dimension< BasisType >() , dimension< BasisType >() );
		if( lump ) triplets.resize( _tCount * BasisInfo< BasisType >::Coefficients );
		else       triplets.resize( _tCount * BasisInfo< BasisType >::Coefficients * BasisInfo< BasisType >::Coefficients );
		ThreadPool::ParallelFor( 0 , triplets.size() , [&]( unsigned int , size_t i ){ triplets[i] = Eigen::Triplet< Real >(0,0,(Real)0); } );
	}
	else
#endif // EIGEN_WORLD_VERSION
	{
		M.resize( dimension< BasisType >() );
		rowSizes = NewPointer< std::atomic< int > >( M.rows ); // need to support atomic increment and set, which is not supported with OpenMP
		ThreadPool::ParallelFor( 0 , M.rows , [&]( unsigned int , size_t i ){ rowSizes[i] = 0; } );


		// First, set the row sizes
		ThreadPool::ParallelFor
		(
			0 , _tCount ,
			[&]( unsigned int , size_t t )
			{
				if( lump ) for( int i=0 ; i<BasisInfo< BasisType >::Coefficients ; i++ ) rowSizes[ index< BasisType >( (int)t , i ) ]++;
				else       for( int i=0 ; i<BasisInfo< BasisType >::Coefficients ; i++ ) rowSizes[ index< BasisType >( (int)t , i ) ] += nonZeroCount[i];
			}
		);
		ThreadPool::ParallelFor( 0 , M.rows , [&]( unsigned int , size_t i ){ M.SetRowSize( i , rowSizes[i] ) , rowSizes[i] = 0; } );
	};

	// Next, set the entries
	ThreadPool::ParallelFor
		(
			0 , _tCount ,
			[&]( unsigned int , size_t t )
			{
				if( lump )
				{
					auto m = RightTriangle< Real >::template GetDiagonalMassMatrix< BasisType >( _g[t] );
					for( int i=0 ; i<BasisInfo< BasisType >::Coefficients ; i++ )
					{
						int ii = index< BasisType >( (int)t , i );
#ifdef EIGEN_WORLD_VERSION
						if constexpr( UseEigen )
							triplets[t*BasisInfo< BasisType >::Coefficients+i] = Eigen::Triplet< Real >( ii , ii , m[i] );
						else
#endif // EIGEN_WORLD_VERSION
							M[ ii ][ rowSizes[ii]++ ] = MatrixEntry< Real , int >( ii , m[i] );
					}
				}
				else
				{
					auto m = newTensors ? RightTriangle< Real >::template GetMassMatrix< BasisType >( _g[t] , newTensors[t] ) : RightTriangle< Real >::template GetMassMatrix< BasisType >( _g[t] );
					for( int i=0 ; i<BasisInfo< BasisType >::Coefficients ; i++ ) 
					{
						bool iAligned;
						int ii = index< BasisType >( (int)t , i , iAligned );

						int idx = 0;
						for( int j=0 ; j<BasisInfo< BasisType >::Coefficients ; j++ ) if( mask(i,j) )
						{
							bool jAligned;
							int jj = index< BasisType >( (int)t , j , jAligned );
#ifdef EIGEN_WORLD_VERSION
							if constexpr( UseEigen )
								triplets[ t * BasisInfo< BasisType >::Coefficients * BasisInfo< BasisType >::Coefficients + BasisInfo< BasisType >::Coefficients * i + j ] = Eigen::Triplet( ii , jj , iAligned==jAligned ? m(i,j) : - m(i,j) );
							else
#endif // EIGEN_WORLD_VERSION
								M[ ii ][ rowSizes[ii]++ ] = MatrixEntry< Real , int >( jj , iAligned==jAligned ? m(i,j) : - m(i,j) );
						}
					}
				}
			}
		);
#ifdef EIGEN_WORLD_VERSION
	if constexpr( !UseEigen ) DeletePointer( rowSizes );
#else // !EIGEN_WORLD_VERSION
	DeletePointer( rowSizes );
#endif // EIGEN_WORLD_VERSION

#ifdef EIGEN_WORLD_VERSION
	if constexpr( UseEigen ) M.setFromTriplets( triplets.begin() , triplets.end() );
	else
#endif // EIGEN_WORLD_VERSION
	{
		// Collapse the duplicate entries (and sort)
		ThreadPool::ParallelFor
		(
			0 , M.rows ,
			[&]( unsigned int , size_t i )
			{
				std::sort( ( MatrixEntry< Real , int > * )GetAddress( M[i] ) , GetAddress( M[i] ) + M.rowSizes[i] , []( MatrixEntry< Real , int > e1 , MatrixEntry< Real , int > e2 ){ return e1.N<e2.N; } );
				int idx = 0;
				for( int j=1 ; j<M.rowSizes[i] ; j++ )
					if( M[i][j].N==M[i][idx].N ) M[i][idx].Value += M[i][j].Value;
					else M[i][++idx] = M[i][j];
				M.ResetRowSize( i , idx+1 );
				for( int j=0 ; j<M.rowSizes[i] ; j++ ) if( M[i][j].N==i ) std::swap( M[i][j] , M[i][0] );
			}
		);
	}
	return M;
}

template< class Real , typename Index >
#ifdef EIGEN_WORLD_VERSION 
template< unsigned int BasisType , unsigned int Degree , bool UseEigen , typename CotangentVectorFieldFunctor /* = std::function< RightTriangle< Real >::CotangentVectorField< Degree > ( unsigned int tIdx ) > */ >
std::conditional_t< UseEigen , Eigen::SparseMatrix< Real > , SparseMatrix< Real , int > > FEM::RiemannianMesh< Real , Index >::derivation( CotangentVectorFieldFunctor v ) const
#else // !EIGEN_WORLD_VERSION 
template< unsigned int BasisType , unsigned int Degree , typename CotangentVectorFieldFunctor /* = std::function< RightTriangle< Real >::CotangentVectorField< Degree > ( unsigned int tIdx ) > */ >
SparseMatrix< Real , int > FEM::RiemannianMesh< Real , Index >::derivation( CotangentVectorFieldFunctor v ) const
#endif // EIGEN_WORLD_VERSION 

{
	// [WARNING] Hard-coded for Hat basis functions

#ifdef EIGEN_WORLD_VERSION
	std::conditional_t< UseEigen , Eigen::SparseMatrix< Real > , SparseMatrix< Real , int > > M;
	std::vector< Eigen::Triplet< Real > > triplets;
#else // !EIGEN_WORLD_VERSION
	SparseMatrix< Real , int > M;
#endif // EIGEN_WORLD_VERSION

	if( BasisType!=BASIS_0_WHITNEY ) MK_ERROR_OUT( "Expected Whitney 0-form basis" );
	typename RightTriangle< Real >::template ScalarField< 1 > hat[3];
	typename RightTriangle< Real >::template VectorField< 0 > dHat[3];

	hat[0].coefficient(0,0) =  1;
	hat[0].coefficient(1,0) = -1;
	hat[0].coefficient(0,1) = -1;
	hat[1].coefficient(1,0) = 1;
	hat[2].coefficient(0,1) = 1;
	for( unsigned int i=0 ; i<3 ; i++ ) dHat[i] = hat[i].differential();

	Pointer( std::atomic< int > ) rowSizes = NullPointer< std::atomic< int > >();

#ifdef EIGEN_WORLD_VERSION
	if constexpr( UseEigen )
	{
		M.resize( dimension< BasisType >() , dimension< BasisType >() );
		triplets.resize( _tCount*3*3 );
	}
#endif // EIGEN_WORLD_VERSION
	{
		M.resize( dimension< BasisType >() );
		rowSizes = NewPointer< std::atomic< int > >( M.rows ); // need to support atomic increment and set, which is not supported with OpenMP
		ThreadPool::ParallelFor( 0 , M.rows , [&]( unsigned int , size_t i ){ rowSizes[i] = 0; } );

		// First, set the row sizes
		ThreadPool::ParallelFor( 0 , _tCount , [&]( unsigned int , size_t t ){ for( int i=0 ; i<3 ; i++ ) rowSizes[ index< BasisType >( t , i ) ] += 3; } );
		ThreadPool::ParallelFor( 0 , M.rows , [&]( unsigned int , size_t i ){ M.SetRowSize( i , rowSizes[i] ) ; rowSizes[i] = 0; } );
	}


	// Next, set the entries
	ThreadPool::ParallelFor
		(
			0 , _tCount ,
			[&]( unsigned int , size_t t )
			{
				typename RightTriangle< Real >::template CotangentVectorField< Degree > vf = _g[t].inverse() * v( t ) * (Real)sqrt( _g[t].determinant() );

				for( unsigned int i=0 ; i<3 ; i++ )
				{
					typename RightTriangle< Real >::template ScalarField< Degree > f = vf.first * dHat[i].first + vf.second * dHat[i].second;
					unsigned int ii = index< BasisType >( t , i );
					for( unsigned int j=0 ; j<3 ; j++ )
					{
						unsigned int jj = index< BasisType >( t , j );
#ifdef EIGEN_WORLD_VERSION
						if constexpr( UseEigen )
							triplets[t*3*3 + i*3 + j ] = Eigen::Triplet< Real >( jj , jj , ( f * hat[j] ).integrateUnitRightSimplex() );
						else
#endif // EIGEN_WORLD_VERSION
							M[ jj ][ rowSizes[jj]++ ] = MatrixEntry< Real , int >( jj , ( f * hat[j] ).integrateUnitRightSimplex() );
					}
				}
			}
		);
#ifdef EIGEN_WORLD_VERSION
	if constexpr( UseEigen ) M.setFromTriplets( triplets.begin() , triplets.end() );
	else
#endif // EIGEN_WORLD_VERSION
	{
		DeletePointer( rowSizes );

		// Collapse the duplicate entries (and sort)
		ThreadPool::ParallelFor
		(
			0 , M.rows ,
			[&]( unsigned int , size_t i )
			{
				std::sort( ( MatrixEntry< Real , int > * )GetAddress( M[i] ) , GetAddress( M[i] ) + M.rowSizes[i] , []( MatrixEntry< Real , int > e1 , MatrixEntry< Real , int > e2 ){ return e1.N<e2.N; } );
				int idx = 0;
				for( int j=1 ; j<M.rowSizes[i] ; j++ )
					if( M[i][j].N==M[i][idx].N ) M[i][idx].Value += M[i][j].Value;
					else M[i][++idx] = M[i][j];
				M.ResetRowSize( i , idx+1 );
				for( int j=0 ; j<M.rowSizes[i] ; j++ ) if( M[i][j].N==i ) std::swap( M[i][j] , M[i][0] );
			}
		);
	}

	return M;
}

template< class Real , typename Index >
#ifdef EIGEN_WORLD_VERSION
template< unsigned int InBasisType , unsigned int OutBasisType , bool UseEigen >
std::conditional_t< UseEigen , Eigen::SparseMatrix< Real > , SparseMatrix< Real , int > > FEM::RiemannianMesh< Real , Index >::dMatrix( void ) const
#else // !EIGEN_WORLD_VERSION
template< unsigned int InBasisType , unsigned int OutBasisType >
SparseMatrix< Real , int > FEM::RiemannianMesh< Real , Index >::dMatrix( void ) const
#endif // EIGEN_WORLD_VERSION
{
	TestBasisType(  InBasisType , "FEM::RiemannianMesh::dMatrix" , false );
	TestBasisType( OutBasisType , "FEM::RiemannianMesh::dMatrix" , false );

	auto Fail = [&] ( void ){ fprintf( stderr , "[ERROR] FEM::RiemannianMesh::dMatrix: %s -> %s\n" , BasisNames[ InBasisType ] , BasisNames[ OutBasisType ] ) , exit( 0 ); };

	bool redundant;
	Matrix< unsigned char , BasisInfo< InBasisType >::Coefficients , BasisInfo< OutBasisType >::Coefficients > mask = RightTriangle< Real >::template GetDMask< InBasisType , OutBasisType >( redundant );
	Point< int , BasisInfo< OutBasisType >::Coefficients > nonZeroCount;
	for( int i=0 ; i<BasisInfo< InBasisType >::Coefficients ; i++ ) for( int j=0 ; j<BasisInfo< OutBasisType >::Coefficients ; j++ ) if( mask(i,j) ) nonZeroCount[j]++;

	Pointer( std::atomic< int > ) rowSizes = NullPointer< std::atomic< int > >();
#ifdef EIGEN_WORLD_VERSION
	std::conditional_t< UseEigen , Eigen::SparseMatrix< Real > , SparseMatrix< Real , int > > D;
	std::vector< Eigen::Triplet< Real > > triplets;
#else // !EIGEN_WORLD_VERSION
	SparseMatrix< Real , int > D;
#endif // EIGEN_WORLD_VERSION

#ifdef EIGEN_WORLD_VERSION
	if constexpr( UseEigen )
	{
		D.resize( dimension< OutBasisType >() , dimension< InBasisType >() );
		triplets.resize( _tCount * BasisInfo< OutBasisType >::Coefficients * BasisInfo< InBasisType >::Coefficients );
	}
	else
#endif // EIGEN_WORLD_VERSION
	{
		D.resize( dimension< OutBasisType >() );
		rowSizes = NewPointer< std::atomic< int > >( D.rows ); // need to support atomic increment and set, which is not supported with OpenMP
		ThreadPool::ParallelFor( 0 , D.rows , [&]( unsigned int , size_t i ){ rowSizes[i] = 0; } );
	}

#ifdef EIGEN_WORLD_VERSION
	if constexpr( !UseEigen )
#endif // EIGEN_WORLD_VERSION
	{
		// First, set the row sizes
		ThreadPool::ParallelFor( 0 , _tCount , [&]( unsigned int , size_t t ){ for( int j=0 ; j<BasisInfo< OutBasisType >::Coefficients ; j++ ) rowSizes[ index< OutBasisType >( (int)t , j ) ] += nonZeroCount[j]; } );
		ThreadPool::ParallelFor( 0 , D.rows , [&]( unsigned int , size_t i ){ D.SetRowSize( i , rowSizes[i] ) , rowSizes[i] = 0; } );
	}

	// Next, set the entries
	ThreadPool::ParallelFor
		(
			0 , _tCount ,
			[&]( unsigned int , size_t t )
			{
				auto d = RightTriangle< Real >::template GetDMatrix< InBasisType , OutBasisType >( _g[t] );
				for( int j=0 ; j<BasisInfo< OutBasisType >::Coefficients ; j++ ) 
				{
					bool jAligned ; int jj = index< OutBasisType >( (int)t , j , jAligned );
					for( int i=0 ; i<BasisInfo< InBasisType >::Coefficients ; i++ ) if( mask(i,j) )
					{
						bool iAligned ; int ii = index< InBasisType >( (int)t , i , iAligned );
#ifdef EIGEN_WORLD_VERSION
						if constexpr( UseEigen )
							triplets[ t * BasisInfo< OutBasisType >::Coefficients * BasisInfo< InBasisType >::Coefficients + j * BasisInfo< InBasisType >::Coefficients + i ]
							= Eigen::Triplet< Real >( jj , ii , iAligned==jAligned ? d(i,j) : -d(i,j) );
						else
#endif // EIGEN_WORLD_VERSION
							D[ jj ][ rowSizes[jj]++ ] = MatrixEntry< Real , int >( ii , iAligned==jAligned ? d(i,j) : -d(i,j) );
					}
				}
			}
		);
#ifdef EIGEN_WORLD_VERSION
	if constexpr( UseEigen ) D.setFromTriplets( triplets.begin() , triplets.end() );
	else
#endif // EIGEN_WORLD_VERSION
	{
		DeletePointer( rowSizes );

		// Collapse the duplicate entries
		ThreadPool::ParallelFor
			(
				0 , D.rows ,
				[&]( unsigned int , size_t i )
				{
					std::sort( D[i] , D[i] + D.rowSizes[i] , []( MatrixEntry< Real , int > e1 , MatrixEntry< Real , int > e2 ){ return e1.N<e2.N; } );
					int idx = 0;
					if( D.rowSizes[i] )
					{	
						for( int j=1 ; j<D.rowSizes[i] ; j++ )
							if( D[i][j].N==D[i][idx].N ) D[i][idx].Value += D[i][j].Value;
							else D[i][++idx] = D[i][j];
						D.ResetRowSize( i , idx+1 );
					}
				}
			);
	}

	if( redundant )
	{
		std::vector< int > valence;
		if( BasisInfo< OutBasisType >::ElementType==ELEMENT_VERTEX )
		{
			valence.resize( _vCount , 0 );
			ThreadPool::ParallelFor( 0 , _tCount , [&]( unsigned int , size_t t ){ for( int v=0 ; v<3 ; v++ ) AddAtomic( valence[ _triangles[t][v] ] , 1 ); } );
		}
		else if( BasisInfo< OutBasisType >::ElementType==ELEMENT_EDGE )
		{
			valence.resize( _edgeMap.size() , 0 );
			ThreadPool::ParallelFor( 0 , _tCount , [&]( unsigned int , size_t t ){ for( int e=0 ; e<3 ; e++ ) AddAtomic( valence[ _edgeMap.edge( (int)t*3+e ) ] , 1 ); } );
		}
		else TestElementType( BasisInfo< OutBasisType >::ElementType , "FEM::RiemannianMesh::dMatrix" , true );

#ifdef EIGEN_WORLD_VERSION
		if constexpr( UseEigen )
		{
			ThreadPool::ParallelFor
			(
				0 , D.outerSize() ,
				[&]( unsigned int , size_t i )
				{
					for( typename Eigen::SparseMatrix< Real >::InnerIterator it(D,i) ; it ; ++it )
						it.valueRef() /= valence[it.row()/BasisInfo< OutBasisType >::CoefficientsPerElement];
				}
			);
		}
		else
#endif // EIGEN_WORLD_VERSION
		{
			ThreadPool::ParallelFor
				(
					0 , D.rows ,
					[&]( unsigned int , size_t i )
					{
						Real scale = (Real)( 1. / valence[i/BasisInfo< OutBasisType >::CoefficientsPerElement] );
						for( int j=0 ; j<D.rowSizes[i] ; j++ ) D[i][j].Value *= scale;
					}
				);
		}
	}
	return D;
}

template< class Real , typename Index >
#ifdef EIGEN_WORLD_VERSION
template< unsigned int BasisType , unsigned int PreBasisType , unsigned int PostBasisType , bool UseEigen >
std::conditional_t< UseEigen , Eigen::SparseMatrix< Real > , SparseMatrix< Real , int > > FEM::RiemannianMesh< Real , Index >::stiffnessMatrix( ConstPointer( SquareMatrix< Real , 2 > ) newTensors ) const
#else // !EIGEN_WORLD_VERSION
template< unsigned int BasisType , unsigned int PreBasisType , unsigned int PostBasisType >
SparseMatrix< Real , int > FEM::RiemannianMesh< Real , Index >::stiffnessMatrix( ConstPointer( SquareMatrix< Real , 2 > ) newTensors ) const
#endif // EIGEN_WORLD_VERSION
{
#ifdef EIGEN_WORLD_VERSION
	using Matrix = std::conditional_t< UseEigen , Eigen::SparseMatrix< Real > , SparseMatrix< Real , int > >;
	auto MassMatrixInverse = [] ( const Matrix& M )
#else // !EIGEN_WORLD_VERSION
	auto MassMatrixInverse = [] ( const SparseMatrix< Real , int >& M )
#endif // EIGEN_WORLD_VERSION
	{
#ifdef EIGEN_WORLD_VERSION
		if constexpr( UseEigen )
		{
			Matrix M_inverse;
			M_inverse.resize( M.rows() , M.cols() );
			std::vector< Eigen::Triplet< Real > > triplets( M.outerSize() );
			ThreadPool::ParallelFor
			(
				0 , M.outerSize() ,
				[&]( unsigned int , size_t i )
				{
					Real sum = 0;
					for( typename Matrix::InnerIterator it( M , i ) ; it ; ++it ) sum += it.value();
					triplets[i] = Eigen::Triplet< Real >( (unsigned int)i , (unsigned int)i , (Real)1./sum );
				}
			);
			return M_inverse;
		}
		else
#endif // EIGEN_WORLD_VERSION
		{
			SparseMatrix< Real , int > M_inverse;
			M_inverse.resize( M.rows );
			ThreadPool::ParallelFor
			(
				0 , M.rows ,
				[&]( unsigned int , size_t i )
				{
					M_inverse.SetRowSize( i , 1 );
					Real sum = 0;
					for( int j=0 ; j<M.rowSizes[i] ; j++ ) sum += M[i][j].Value;
					M_inverse[i][0] = MatrixEntry< Real , int >( (int)i , (Real)(1./sum) );
				}
			);
			return M_inverse;
		}
	};

#ifdef EIGEN_WORLD_VERSION
	Matrix S;
#else // !EIGEN_WORLD_VERSION
	SparseMatrix< Real , int > S;
#endif // EIGEN_WORLD_VERSION
	TestBasisType( BasisType , "FEM::RiemannianMesh::stiffnessMatrix" , false );
	switch( BasisInfo< BasisType >::Dimension )
	{
		case 0:
		{	
			TestBasisType( PostBasisType , "FEM::RiemannianMesh::stiffnessMatrix" , false );
			if( BasisInfo< PostBasisType >::Dimension!=BasisInfo< BasisType >::Dimension+1 ) fprintf( stderr , "[ERROR] FEM::RiemannianMesh::stiffnessMatrix: Incompatible basis type: %s -> %s\n" , BasisNames[BasisType] , BasisNames[PostBasisType] ) , exit( 0 );
#ifdef EIGEN_WORLD_VERSION
			Matrix M2 = massMatrix< PostBasisType , UseEigen >( BasisInfo< BasisType >::Lumpable , newTensors );
			Matrix D1 = dMatrix< BasisType , PostBasisType , UseEigen >( );
#else // !EIGEN_WORLD_VERSION
			SparseMatrix< Real , int > M2 = massMatrix< PostBasisType >( BasisInfo< BasisType >::Lumpable , newTensors );
			SparseMatrix< Real , int > D1 = dMatrix< BasisType , PostBasisType >( );
#endif // EIGEN_WORLD_VERSION
			S = D1.transpose() * M2 * D1;

			break;
		}
		case 1:
		{
			TestBasisType(  PreBasisType , "FEM::RiemannianMesh::stiffnessMatrix" , false );
			TestBasisType( PostBasisType , "FEM::RiemannianMesh::stiffnessMatrix" , false );
			if( BasisInfo< PreBasisType >::Dimension!=BasisInfo< BasisType >::Dimension-1 ) MK_ERROR_OUT( "Incompatible basis type: " , BasisNames[PreBasisType] , " - > " , BasisNames[BasisType] );
			if( BasisInfo< PostBasisType >::Dimension!=BasisInfo< BasisType >::Dimension+1 ) MK_ERROR_OUT( "Incompatible basis type: " , BasisNames[BasisType] , " -> " , BasisNames[PostBasisType] );
#ifdef EIGEN_WORLD_VERSION
			Matrix M0 = massMatrix<  PreBasisType , UseEigen >( true , newTensors );
			Matrix M1 = massMatrix<     BasisType , UseEigen >( BasisInfo< BasisType >::Lumpable , newTensors );
			Matrix M2 = massMatrix< PostBasisType , UseEigen >( BasisInfo< BasisType >::Lumpable , newTensors );
			Matrix D0 = dMatrix< PreBasisType ,     BasisType , UseEigen >( );
			Matrix D1 = dMatrix<    BasisType , PostBasisType , UseEigen >( );
#else // !EIGEN_WORLD_VERSION
			SparseMatrix< Real , int > M0 = massMatrix<  PreBasisType >( true , newTensors );
			SparseMatrix< Real , int > M1 = massMatrix<     BasisType >( BasisInfo< BasisType >::Lumpable , newTensors );
			SparseMatrix< Real , int > M2 = massMatrix< PostBasisType >( BasisInfo< BasisType >::Lumpable , newTensors );
			SparseMatrix< Real , int > D0 = dMatrix< PreBasisType ,     BasisType >( );
			SparseMatrix< Real , int > D1 = dMatrix<    BasisType , PostBasisType >( );
#endif // EIGEN_WORLD_VERSION
#if 1
			S = M1 * D0 * MassMatrixInverse( M0 ) * D0.transpose() * M1 + D1.transpose() * M2 * D1;
#else
			// Just the part that doesn't require lumping/inversion
			S = D1.transpose() * M2 * D1;
#endif
			break;
		}
		case 2:
		{
			TestBasisType(  PreBasisType , "FEM::RiemannianMesh::stiffnessMatrix" , false );
			if( BasisInfo< PreBasisType >::Dimension!=BasisInfo< BasisType >::Dimension-1 ) MK_ERROR_OUT( "Incompatible basis type: " , BasisNames[PreBasisType] , " -> " , BasisNames[BasisType] );
#ifdef EIGEN_WORLD_VERSION
			Matrix M0 = massMatrix< PreBasisType , UseEigen >( true , newTensors );
			Matrix M1 = massMatrix<    BasisType , UseEigen >( BasisInfo< BasisType >::Lumpable , newTensors );
			Matrix D0 = dMatrix< PreBasisType , BasisType , UseEigen >( );
#else // !EIGEN_WORLD_VERSION
			SparseMatrix< Real , int > M0 = massMatrix< PreBasisType >( true , newTensors );
			SparseMatrix< Real , int > M1 = massMatrix<    BasisType >( BasisInfo< BasisType >::Lumpable , newTensors );
			SparseMatrix< Real , int > D0 = dMatrix< PreBasisType , BasisType >( );
#endif // EIGEN_WORLD_VERSION
			S = M1 * D0 * MassMatrixInverse( M0 ) * D0.transpose() * M1;
			break;
		}
		default: TestBasisType( BasisType , "FEM::RiemannianMesh::stiffnessMatrix" , true );
	}

#ifdef EIGEN_WORLD_VERSION
	if constexpr( !UseEigen )
#endif // EIGEN_WORLD_VERSION
	{
		ThreadPool::ParallelFor
		(
			0 , S.rows ,
			[&]( unsigned int , size_t i )
			{
				std::sort( S[i] , S[i] + S.rowSizes[i] , []( MatrixEntry< Real , int > e1 , MatrixEntry< Real , int > e2 ){ return e1.N<e2.N; } );
				int idx=0;
				for( int j=1 ; j<S.rowSizes[i] ; j++ )
					if( S[i][j].N==S[i][idx].N ) S[i][idx].Value += S[i][j].Value;
					else idx++ , S[i][idx] = S[i][j];
				S.ResetRowSize( i , idx+1 );
				for( int j=0 ; j<S.rowSizes[i] ; j++ ) if( S[i][j].N==i ) std::swap( S[i][j] , S[i][0] );
			}
		);
	}
	return S;
}

template< class Real , typename Index >
#ifdef EIGEN_WORLD_VERSION
template< unsigned int BasisType , bool UseEigen >
std::conditional_t< UseEigen , Eigen::SparseMatrix< Real > , SparseMatrix< Real , int > > FEM::RiemannianMesh< Real , Index >::stiffnessMatrix( void ) const
#else // !EIGEN_WORLD_VERSION
template< unsigned int BasisType >
SparseMatrix< Real , int > FEM::RiemannianMesh< Real , Index >::stiffnessMatrix( void ) const
#endif // EIGEN_WORLD_VERSION
{
	TestBasisType( BasisType , "FEM::RiemannianMesh::stiffnessMatrix" , false );
#ifdef EIGEN_WORLD_VERSION
	std::conditional_t< UseEigen , Eigen::SparseMatrix< Real > , SparseMatrix< Real , int > > S;
	switch( BasisInfo< BasisType >::Dimension )
	{
	case 0: S = stiffnessMatrix< BasisType , BasisType , BASIS_1_WHITNEY , UseEigen >( ) ; break;
	case 1:
		if( BasisType==BASIS_1_WHITNEY ) S = stiffnessMatrix< BasisType , BASIS_0_WHITNEY , BASIS_2_WHITNEY , UseEigen >( );
		else                             S = stiffnessMatrix< BasisType , BASIS_0_WHITNEY , BASIS_2_VERTEX_CONSTANT , UseEigen >( );
		break;
	case 2: S = stiffnessMatrix< BasisType , BASIS_1_WHITNEY , BasisType , UseEigen>( ) ; break;
	default: TestBasisType( BasisType , "FEM::RiemannianMesh::stiffnessMatrix" , true );
}
#else // !EIGEN_WORLD_VERSION
	SparseMatrix< Real , int > S;
	switch( BasisInfo< BasisType >::Dimension )
	{
	case 0: S = stiffnessMatrix< BasisType , BasisType       , BASIS_1_WHITNEY >( ) ; break;
	case 1:
		if( BasisType==BASIS_1_WHITNEY ) S = stiffnessMatrix< BasisType , BASIS_0_WHITNEY , BASIS_2_WHITNEY  >( );
		else                             S = stiffnessMatrix< BasisType , BASIS_0_WHITNEY , BASIS_2_VERTEX_CONSTANT >( );
		break;
	case 2: S = stiffnessMatrix< BasisType , BASIS_1_WHITNEY , BasisType       >( ) ; break;
	default: TestBasisType( BasisType , "FEM::RiemannianMesh::stiffnessMatrix" , true );
	}
#endif // EIGEN_WORLD_VERSION
	return S;
}

template< class Real , typename Index >
inline Real FEM::RiemannianMesh< Real , Index >::getIntegral( ConstPointer( Real ) coefficients ) const
{
	Real integral = (Real)0;
	std::vector< Real > _integrals( ThreadPool::NumThreads() , 0 );
	ThreadPool::ParallelFor
		(
			0 , _tCount ,
			[&]( unsigned int t , size_t i )
			{
				SquareMatrix< Real , 3 > mass = FEM::RightTriangle< Real >::template GetMassMatrix< FEM::BASIS_0_WHITNEY >( _g[i] );
				for( int j=0 ; j<3 ; j++ )
				{
					Real sum = (Real)0;
					for( int k=0 ; k<3 ; k++ ) sum += mass(j,k);
					_integrals[t] += coefficients[ _triangles[i][j] ] * sum;
				}
			}
		);
	for( unsigned int t=0 ; t<_integrals.size() ; t++ ) integral += _integrals[t];
	return integral;
}

template< class Real , typename Index >
inline Real FEM::RiemannianMesh< Real , Index >::getDotProduct( ConstPointer( Real ) coefficients1 , ConstPointer( Real ) coefficients2 , bool lump ) const
{
	Real dotProduct = (Real)0;
	std::vector< Real > _dotProducts( ThreadPool::NumThreads() , 0 );
	ThreadPool::ParallelFor
		(
			0 , _tCount ,
			[&]( unsigned int t , size_t i )
			{
				if( lump )
				{
					Point< Real , 3 > mass = RightTriangle< Real >::GetDiagonalMassMatrix( _g[i] );
					for( int j=0 ; j<3 ; j++ ) _dotProducts[t] += mass[j] * coefficients1[ _triangles[i][j] ] * coefficients2[ _triangles[i][j] ];
				}
				else
				{
					SquareMatrix< Real , 3 > mass = RightTriangle< Real >::GetMassMatrix( _g[i] );
					for( int j=0 ; j<3 ; j++ ) for( int k=0 ; k<3 ; k++ ) _dotProducts[t] += mass(j,k) * coefficients1[ _triangles[i][j] ] * coefficients2[ _triangles[i][k] ];
				}
			}
		);
	for( unsigned int t=0 ; t<_dotProducts.size() ; t++ ) dotProduct += _dotProducts[t];
	return dotProduct;
}

template< typename Index >
FEM::EdgeMap::EdgeMap( ConstPointer( SimplexIndex< 2 , Index > ) triangles , size_t tCount )
{
	_eCount = 0;
	std::unordered_map< long long , int > edgeMap;
	for( int t=0 ; t<tCount ; t++ ) for( int v=0 ; v<3 ; v++ )
	{
		long long key = EdgeKey( triangles[t][(v+1)%3] , triangles[t][(v+2)%3] );
		if( edgeMap.find(key)==edgeMap.end() ) edgeMap[key] = (int)_eCount++;
	}
	_he2e = AllocPointer< int >( tCount*3 );
	_e2he = AllocPointer< int >( 2*_eCount );
	memset( _e2he , -1 , sizeof(int) * 2 * _eCount );

	for( int t=0 ; t<tCount ; t++ ) for( int v=0 ; v<3 ; v++ )
	{
		long long key = EdgeKey( triangles[t][(v+1)%3] , triangles[t][(v+2)%3] );
		int heIdx = t*3+v , eIdx = edgeMap[key];
		if( _e2he[2*eIdx]==-1 ) _e2he[eIdx*2+0] = heIdx , _he2e[heIdx] =  (eIdx+1);
		else                    _e2he[eIdx*2+1] = heIdx , _he2e[heIdx] = -(eIdx+1);
	}
}

template< typename Real >
template< unsigned int Degree >
typename FEM::RightTriangle< Real >::template CotangentVectorField< Degree-1 > FEM::RightTriangle< Real >::ScalarField< Degree >::differential( void ) const { return std::make_pair( d(0) , d(1) ); }

template< typename Real , unsigned int Degree1 , unsigned int Degree2 >
typename FEM::RightTriangle< Real >::template TangentVectorField< Degree1 + Degree2 > operator * ( const typename FEM::RightTriangle< Real >::template ScalarField< Degree1 > &s , const typename FEM::RightTriangle< Real >::template TangentVectorField< Degree2 > & v ){ return std::make_pair( v.first*s , v.second*s ); }

template< typename Real , unsigned int Degree1 , unsigned int Degree2 >
typename FEM::RightTriangle< Real >::template TangentVectorField< Degree1 + Degree2 > operator * ( const typename FEM::RightTriangle< Real >::template TangentVectorField< Degree1 > & v , const typename FEM::RightTriangle< Real >::template ScalarField< Degree2 > &s ){ return std::make_pair( v.first*s , v.second*s ); }

template< typename Real , unsigned int Degree >
typename FEM::RightTriangle< Real >::template TangentVectorField< Degree > operator * ( const SquareMatrix< double , 2 > &m , const typename FEM::RightTriangle< Real >::template TangentVectorField< Degree > &v )
{
	return std::make_pair( m(0,0) * v.first + m(1,0) * v.second , m(0,1) * v.first + m(1,1) * v.second );
}

template< typename Real , unsigned int Degree1 , unsigned int Degree2 >
typename FEM::RightTriangle< Real >::template CotangentVectorField< Degree1 + Degree2 > operator * ( const typename FEM::RightTriangle< Real >::template ScalarField< Degree1 > &s , const typename FEM::RightTriangle< Real >::template CotangentVectorField< Degree2 > & v ){ return std::make_pair( v.first*s , v.second*s ); }

template< typename Real , unsigned int Degree1 , unsigned int Degree2 >
typename FEM::RightTriangle< Real >::template CotangentVectorField< Degree1 + Degree2 > operator * ( const typename FEM::RightTriangle< Real >::template CotangentVectorField< Degree1 > & v , const typename FEM::RightTriangle< Real >::template ScalarField< Degree2 > &s ){ return std::make_pair( v.first*s , v.second*s ); }

template< typename Real , unsigned int Degree >
typename FEM::RightTriangle< Real >::template CotangentVectorField< Degree > operator * ( const SquareMatrix< double , 2 > &m , const typename FEM::RightTriangle< Real >::template CotangentVectorField< Degree > &v )
{
	return std::make_pair( m(0,0) * v.first + m(1,0) * v.second , m(0,1) * v.first + m(1,1) * v.second );
}
