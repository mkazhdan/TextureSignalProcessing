/*
Copyright (c) 2019, Michael Kazhdan
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

///////////////////
// Polynomial 0D //
///////////////////
template< unsigned int Degree , typename Real >
Polynomial< 0 , Degree , Real >::Polynomial( void ){ memset( _coefficients , 0 , sizeof( _coefficients ) ); }

template< unsigned int Degree , typename Real >
Polynomial< 0 , Degree , Real >::Polynomial( Real c ) : Polynomial() { _coefficients[0] = c; }

template< unsigned int Degree , typename Real >
unsigned int Polynomial< 0 , Degree , Real >::_setCoefficients( const Real *coefficients , unsigned int maxDegree )
{
	_coefficients[0] = coefficients[0];
	return 1;
}
template< unsigned int Degree , typename Real >
unsigned int Polynomial< 0 , Degree , Real >::_getCoefficients( Real *coefficients , unsigned int maxDegree ) const
{
	coefficients[0] = _coefficients[0];
	return 1;
}
template< unsigned int Degree , typename Real >
Point< Real , Polynomial< 0 , Degree , Real >::NumCoefficients > Polynomial< 0 , Degree , Real >::coefficients( void ) const
{
	return Point< Real , 1 >( _coefficients[0] );
}
template< unsigned int Degree , typename Real >
Polynomial< 0 , Degree , Real >::Polynomial( Point< Real , NumCoefficients > coefficients ){ _setCoefficients( coefficients , 0 ); }

template<unsigned int Degree , typename Real >
void Polynomial< 0 , Degree , Real >::SetDegrees( unsigned int coefficientIndex , unsigned int degrees[/*0*/] ){}

template< unsigned int Degree , typename Real >
template< unsigned int _Degree >
Polynomial< 0 , Degree , Real >::Polynomial( const Polynomial< 0 , _Degree , Real > &p )
{
	_coefficients[0] = p._coefficients[0];
}

template< unsigned int Degree , typename Real >
template< unsigned int _Degree >
Polynomial< 0 , Degree , Real > &Polynomial< 0 , Degree , Real >::operator= ( const Polynomial< 0 , _Degree , Real > &p )
{
	_coefficients[0] = p._coefficients[0];
	return *this;
}

template< unsigned int Degree , typename Real >
const Real &Polynomial< 0 , Degree , Real >::_coefficient( const unsigned int indices[] , unsigned int maxDegree ) const
{
	return _coefficients[0];
}

template< unsigned int Degree , typename Real >
Real &Polynomial< 0 , Degree , Real >::_coefficient( const unsigned int indices[] , unsigned int maxDegree )
{
	return _coefficients[0];
}

template< unsigned int Degree , typename Real >
Real Polynomial< 0 , Degree , Real >::_evaluate( const Real coordinates[] , unsigned int maxDegree ) const { return _coefficients[0]; }

template< unsigned int Degree , typename Real >
template< unsigned int _Dim >
Polynomial< _Dim , Degree , Real > Polynomial< 0 , Degree , Real >::_pullBack( const Matrix< Real , _Dim+1 , 0 > &A , unsigned int maxDegree ) const
{
	return Polynomial< _Dim , Degree , Real >( _coefficients[0] );
}

template< unsigned int Degree , typename Real >
bool Polynomial< 0 , Degree , Real >::_isZero( unsigned int maxDegree ) const { return _coefficients[0]==0; }

template< unsigned int Degree , typename Real >
bool Polynomial< 0 , Degree , Real >::_isConstant( unsigned int maxDegree ) const { return true; }

template< unsigned int Degree , typename Real >
const Real &Polynomial< 0 , Degree , Real >::coefficient( void ) const { return _coefficients[0]; }

template< unsigned int Degree , typename Real >
Real &Polynomial< 0 , Degree , Real >::coefficient( void ) { return _coefficients[0]; }

template< unsigned int Degree , typename Real >
Real &Polynomial< 0 , Degree , Real >::coefficient( const unsigned int *d ) { return _coefficients[0]; }

template< unsigned int Degree , typename Real >
const Real &Polynomial< 0 , Degree , Real >::coefficient( const unsigned int *d ) const { return _coefficients[0]; }

template< unsigned int Degree , typename Real >
Real &Polynomial< 0 , Degree , Real >::coefficient( unsigned int *d ) { return _coefficients[0]; }

template< unsigned int Degree , typename Real >
const Real &Polynomial< 0 , Degree , Real >::coefficient( unsigned int *d ) const { return _coefficients[0]; }

template< unsigned int Degree , typename Real >
Real Polynomial< 0 , Degree , Real >::operator()( void ) const { return _coefficients[0]; }

template< unsigned int Degree , typename Real >
Real Polynomial< 0 , Degree , Real >::operator()( Point< Real , 0 > p ) const { return _coefficients[0]; }

template< unsigned int Degree , typename Real >
Polynomial< 0 , (Degree>1) ? Degree-1 : 0 , Real > Polynomial< 0 , Degree , Real >::d( unsigned int ) const
{
	return Polynomial< 0 , (Degree>1) ? Degree-1 : 0 , Real >((Real)0);
}

template< unsigned int Degree , typename Real >
template< unsigned int _Dim >
Polynomial< _Dim-1 , Degree , Real > Polynomial< 0 , Degree , Real >::operator()( const Matrix< Real , _Dim , 0 > &A ) const
{
	static_assert( _Dim!=0 , "Dimension cannot be negative" );
	return Polynomial< _Dim-1 , Degree , Real >( _coefficients[0] );
}

template< unsigned int Degree , typename Real >
Real Polynomial< 0 , Degree , Real >::integrateUnitCube( void ) const { return _coefficients[0]; }

template< unsigned int Degree , typename Real >
Real Polynomial< 0 , Degree , Real >::integrateUnitRightSimplex( void ) const { return _coefficients[0]; }

template< unsigned int Degree , typename Real >
void Polynomial< 0 , Degree , Real >::Scale( Real s ){ _coefficients[0] *= s; }

template< unsigned int Degree , typename Real >
void Polynomial< 0 , Degree , Real >::Add( const Polynomial< 0 , Degree , Real > &p ){ _coefficients[0] += p._coefficients[0]; }

template< unsigned int Degree1 , unsigned int Degree2 , typename Real >
Polynomial< 0 , Degree1 + Degree2 , Real > operator * ( const Polynomial< 0 , Degree1 , Real > &p1 , const Polynomial< 0 , Degree2 , Real > &p2 )
{
	return Polynomial< 0 , Degree1 + Degree2 , Real >( p1._coefficients[0] * p2._coefficients[0] );
}

template< unsigned int Degree1 , unsigned int Degree2 , typename Real >
Polynomial< 0 , Max< Degree1 , Degree2 >::Value , Real > operator + ( const Polynomial< 0 , Degree1 , Real > &p1 , const Polynomial< 0 , Degree2 , Real > &p2 )
{
	return Polynomial< 0 , Max< Degree1 , Degree2 >::Value , Real >( p1._coefficients[0] + p2._coefficients[0] );
}

template< unsigned int Degree1 , unsigned int Degree2 , typename Real >
Polynomial< 0 , Max< Degree1 , Degree2 >::Value , Real > operator - ( const Polynomial< 0 , Degree1 , Real > &p1 , const Polynomial< 0 , Degree2 , Real > &p2 ){ return p1 + (-p2); }

template< unsigned int Degree , typename Real >
bool Polynomial< 0 , Degree , Real >::_print( std::ostream &ostream , const std::string varNames[] , bool first ) const
{
	if( _coefficients[0] )
	{
		if( first )
		{
			if( _coefficients[0]>0 ) ostream << " " << _coefficients[0];
			else                     ostream <<        _coefficients[0];
		}
		else
		{
			if( _coefficients[0]<0 ) ostream << " - " << -_coefficients[0];
			else                     ostream << " + " <<  _coefficients[0];
		}
		first = false;
	}
	return first;
}

template< unsigned int Degree , typename Real >
bool Polynomial< 0 , Degree , Real >::_print( std::ostream &ostream , const std::string varNames[] , std::string suffix , bool first ) const
{
	if( _coefficients[0] )
	{
		if( first )
		{
			if     ( _coefficients[0]== 1 ) ostream << " " << suffix;
			else if( _coefficients[0]==-1 ) ostream << "-" << suffix; 
			else
			{
				if( _coefficients[0]>0 ) ostream << " " << _coefficients[0] << "*" << suffix;
				else                     ostream <<        _coefficients[0] << "*" << suffix;
			}
		}
		else
		{
			if( _coefficients[0]<0 )
			{
				if( _coefficients[0]==-1 ) ostream << " - " << suffix;
				else                       ostream << " - " << -_coefficients[0] << "*" << suffix;
			}
			else
			{
				if( _coefficients[0]==1 ) ostream << " + " << suffix;
				else                      ostream << " + " <<  _coefficients[0] << "*" << suffix;
			}
		}
		first = false;
	}
	return first;
}

template< unsigned int Degree , typename Real >
Matrix< Real , Polynomial< 0 , Degree , Real >::NumCoefficients , Polynomial< 0 , Degree , Real >::NumCoefficients > Polynomial< 0 , Degree , Real >::EvaluationMatrix( const Point< Real , 0 > positions[NumCoefficients] )
{
	SquareMatrix< Real , NumCoefficients > E;
	E(0,0) = 1;
	return E;
}

///////////////////////////////
// Polynomial Dim-dimensions //
///////////////////////////////
template< unsigned int Dim , unsigned int Degree , typename Real > Polynomial< Dim , Degree , Real >::Polynomial( void ){}

template< unsigned int Dim , unsigned int Degree , typename Real > Polynomial< Dim , Degree , Real >::Polynomial( Real c )
{
	_polynomials[0] = Polynomial< Dim-1 , Degree , Real >( c );
}

template< unsigned int Dim , unsigned int Degree , typename Real >
unsigned int Polynomial< Dim , Degree , Real >::_getCoefficients( Real *coefficients , unsigned int maxDegree ) const
{
	unsigned int offset = 0;
	for( unsigned int d=0 ; d<=maxDegree ; d++ ) offset += _polynomials[d]._getCoefficients( coefficients + offset , maxDegree - d );
	return offset;
}
template< unsigned int Dim , unsigned int Degree , typename Real >
unsigned int Polynomial< Dim , Degree , Real >::_setCoefficients( const Real *coefficients , unsigned int maxDegree )
{
	unsigned int offset = 0;
	for( unsigned int d=0 ; d<=maxDegree ; d++ ) offset += _polynomials[d]._setCoefficients( coefficients + offset , maxDegree - d );
	return offset;
}
template< unsigned int Dim , unsigned int Degree , typename Real >
Polynomial< Dim , Degree , Real >::Polynomial( Point< Real , NumCoefficients > coefficients ){ _setCoefficients( &coefficients[0] , Degree ); }

template< unsigned int Dim , unsigned int Degree , typename Real >
Point< Real , Polynomial< Dim , Degree , Real >::NumCoefficients > Polynomial< Dim , Degree , Real >::coefficients( void ) const
{
	Point< Real , NumCoefficients > c;
	_getCoefficients( &c[0] , Degree );
	return c;
}

template< unsigned int Dim , unsigned int Degree , typename Real >
template< unsigned int _D >
typename std::enable_if< (_D<Degree ) >::type Polynomial< Dim , Degree , Real >::_SetDegrees( unsigned int coefficientIndex , unsigned int degrees[Dim] )
{
	if( coefficientIndex < Polynomial< Dim-1 , Degree-_D , Real >::NumCoefficients )
	{
		degrees[0] = _D;
		Polynomial< Dim-1 , Degree-_D , Real >::SetDegrees( coefficientIndex , degrees+1 );
	}
	else _SetDegrees< _D+1 >( coefficientIndex - Polynomial< Dim-1 , Degree-_D , Real >::NumCoefficients , degrees );
}

template< unsigned int Dim , unsigned int Degree , typename Real >
template< unsigned int _D >
typename std::enable_if< (_D==Degree ) >::type Polynomial< Dim , Degree , Real >::_SetDegrees( unsigned int coefficientIndex , unsigned int degrees[Dim] )
{
	if( coefficientIndex < Polynomial< Dim-1 , 0 , Real >::NumCoefficients )
	{
		degrees[0] = _D;
		Polynomial< Dim-1 , Degree-_D , Real >::SetDegrees( coefficientIndex , degrees+1 );
	}
	else MK_ERROR_OUT( "coefficient index too big" );
}

template< unsigned int Dim , unsigned int Degree , typename Real >
void Polynomial< Dim , Degree , Real >::SetDegrees( unsigned int coefficientIndex , unsigned int degrees[Dim] )
{
	return _SetDegrees< 0 >( coefficientIndex , degrees );
}

template< unsigned int Dim , unsigned int Degree , typename Real >
template< unsigned int _Degree >
Polynomial< Dim , Degree , Real >::Polynomial( const Polynomial< Dim , _Degree , Real > &p )
{
	for( int d=0 ; d<=Degree && d<=_Degree ; d++ ) _polynomials[d] = p._polynomials[d];
	for( int d=_Degree+1 ; d<=Degree ; d++ ) _polynomials[d] = Polynomial< Dim-1 , Degree , Real >();
}

template< unsigned int Dim , unsigned int Degree , typename Real >
template< unsigned int _Degree >
Polynomial< Dim , Degree , Real > &Polynomial< Dim , Degree , Real >::operator= ( const Polynomial< Dim , _Degree , Real > &p )
{
	for( int d=0 ; d<=Degree && d<=_Degree ; d++ ) _polynomials[d] = p._polynomials[d];
	for( int d=_Degree+1 ; d<=Degree ; d++ ) _polynomials[d] = Polynomial< Dim-1 , Degree , Real >();
	return *this;
}

template< unsigned int Dim , unsigned int Degree , typename Real >
const Real &Polynomial< Dim , Degree , Real >::_coefficient( const unsigned int indices[] , unsigned int maxDegree ) const
{
	if( indices[0]>maxDegree ) MK_ERROR_OUT( "degree out of bounds: " , indices[0] , " > " , maxDegree );
	return _polynomials[ indices[0] ]._coefficient( indices+1 , maxDegree-indices[0] );
}

template< unsigned int Dim , unsigned int Degree , typename Real >
Real& Polynomial< Dim , Degree , Real >::_coefficient( const unsigned int indices[] , unsigned int maxDegree )
{
	if( indices[0]>maxDegree ) MK_ERROR_OUT( "degree out of bounds: " , indices[0] , " > " , maxDegree );
	return _polynomials[ indices[0] ]._coefficient( indices+1 , maxDegree-indices[0] );
}

template< unsigned int Dim , unsigned int Degree , typename Real >
Real Polynomial< Dim , Degree , Real >::_evaluate( const Real coordinates[] , unsigned int maxDegree ) const
{
	Real sum = 0 , tmp = 1;
	for( unsigned int d=0 ; d<=maxDegree ; d++ )
	{
		sum += _polynomials[d]._evaluate( coordinates+1 , maxDegree-d ) * tmp;
		tmp *= coordinates[0];
	}
	return sum;
}

template< unsigned int Dim , unsigned int Degree , typename Real >
template< unsigned int _Dim >
Polynomial< _Dim , Degree , Real > Polynomial< Dim , Degree , Real >::_pullBack( const Matrix< Real , _Dim+1 , Dim > &A , unsigned int maxDegree ) const
{
	unsigned int indices[ _Dim>0 ? _Dim : 1 ]; 
	Polynomial< _Dim , 1 , Real > _p;
	Matrix< Real , _Dim+1 , Dim-1 > _A;

	for( int i=0 ; i<_Dim+1 ; i++ ) for( int j=1 ; j<Dim ; j++ ) _A(i,j-1) = A(i,j) ;

	for( int i=0 ; i<_Dim ; i++ ) indices[i] = 0;
	_p.coefficient( indices ) = A( _Dim , 0 );
	for( int i=0 ; i<_Dim ; i++ )
	{
		indices[i] = 1;
		_p.coefficient( indices ) = A( i , 0 );
		indices[i] = 0;
	}

	Polynomial< _Dim , Degree , Real > p( 0. ) , __p( 1. );
	for( unsigned int d=0 ; d<=maxDegree ; d++ )
	{
		p += _polynomials[d].template _pullBack< _Dim >( _A , maxDegree-d ) * __p;
		__p = __p * _p;
	}
	return p;
}

template< unsigned int Dim , unsigned int Degree , typename Real >
bool Polynomial< Dim , Degree , Real >::_isZero( unsigned int maxDegree ) const
{
	for( unsigned int d=0 ; d<=maxDegree ; d++ ) if( !_polynomials[d]._isZero( maxDegree-d ) ) return false;
	return true;
}

template< unsigned int Dim , unsigned int Degree , typename Real >
bool Polynomial< Dim , Degree , Real >::_isConstant( unsigned int maxDegree ) const
{
	if( !_polynomials[0]._isConstant( Degree ) ) return false;
	for( unsigned int d=1 ; d<=maxDegree ; d++ ) if( !_polynomials[d]._isZero( maxDegree-d ) ) return false;
	return true;
}

template< unsigned int Dim , unsigned int Degree , typename Real >
template< typename ... UnsignedInts >
const Real &Polynomial< Dim , Degree , Real >::coefficient( UnsignedInts ... indices ) const
{
	static_assert( sizeof...(indices)==Dim  , "[ERROR] Polynomial< Dim , Degree , Real >::coefficient: Invalid number of indices" );
	unsigned int _indices[] = { indices ... };
	return _coefficient( _indices , Degree );
}

template< unsigned int Dim , unsigned int Degree , typename Real >
template< typename ... UnsignedInts >
Real &Polynomial< Dim , Degree , Real >::coefficient( UnsignedInts ... indices )
{
	static_assert( sizeof...(indices)==Dim , "[ERROR] Polynomial< Dim , Degree , Real >::coefficient: Invalid number of indices" );
	unsigned int _indices[] = { (unsigned int)indices ... };
	return _coefficient( _indices , Degree );
}

template< unsigned int Dim , unsigned int Degree , typename Real >
const Real &Polynomial< Dim , Degree , Real >::coefficient( const unsigned int indices[Dim] ) const
{
	return _coefficient( indices , Degree );
}

template< unsigned int Dim , unsigned int Degree , typename Real >
Real &Polynomial< Dim , Degree , Real >::coefficient( const unsigned int indices[Dim] )
{
	return _coefficient( indices , Degree );
}

template< unsigned int Dim , unsigned int Degree , typename Real >
const Real &Polynomial< Dim , Degree , Real >::coefficient( unsigned int indices[Dim] ) const
{
	return _coefficient( indices , Degree );
}

template< unsigned int Dim , unsigned int Degree , typename Real >
Real &Polynomial< Dim , Degree , Real >::coefficient( unsigned int indices[Dim] )
{
	return _coefficient( indices , Degree );
}

template< unsigned int Dim , unsigned int Degree , typename Real >
template< typename ... Reals >
Real Polynomial< Dim , Degree , Real >::operator()( Reals ... coordinates ) const
{
	static_assert( sizeof...(coordinates)==Dim , "[ERROR] Polynomial< Dim , Degree , Real >::operator(): Invalid number of coordinates" );
	Real _coordinates[] = { coordinates... };
	return _evaluate( _coordinates , Degree );
}

template< unsigned int Dim , unsigned int Degree , typename Real >
Real Polynomial< Dim , Degree , Real >::operator()( Point< Real , Dim > p ) const { return _evaluate( &p[0] , Degree ); }

template< unsigned int Dim , unsigned int Degree , typename Real >
template< typename ... Reals >
Point< Real , Dim > Polynomial< Dim , Degree , Real >::gradient( Reals ... coordinates ) const
{
	static_assert( sizeof...(coordinates)==Dim , "[ERROR] Polynomial< Dim , Degree , Real >::operator(): Invalid number of coordinates" );
	Real _coordinates[] = { coordinates... };
	Point< double , Dim > p( _coordinates );
	return gradient( p );
}

template< unsigned int Dim , unsigned int Degree , typename Real >
Point< Real , Dim > Polynomial< Dim , Degree , Real >::gradient( Point< Real , Dim > p ) const
{
	Point< Real , Dim > grad;
	if constexpr( Dim>0 ) for( unsigned int d=0 ; d<Dim ; d++ ) grad[d] = this->d( d )( p );
	return grad;
}

template< unsigned int Dim , unsigned int Degree , typename Real >
template< typename ... Reals >
SquareMatrix< Real , Dim > Polynomial< Dim , Degree , Real >::hessian( Reals ... coordinates ) const
{
	static_assert( sizeof...(coordinates)==Dim , "[ERROR] Polynomial< Dim , Degree , Real >::operator(): Invalid number of coordinates" );
	Real _coordinates[] = { coordinates... };
	Point< double , Dim > p( _coordinates );
	return hessian( p );
}

template< unsigned int Dim , unsigned int Degree , typename Real >
SquareMatrix< Real , Dim > Polynomial< Dim , Degree , Real >::hessian( Point< Real , Dim > p ) const
{
	SquareMatrix< Real , Dim > hessian;
	if constexpr( Dim>0 ) 
	{
		for( unsigned int d1=0 ; d1<Dim ; d1++ )
		{
			Point< double , Dim > grad = d(d1).gradient( p );
			for( unsigned int d2=0 ; d2<Dim ; d2++ ) hessian( d1 , d2 ) = grad[d2];
		}
		hessian = ( hessian + hessian.transpose() ) / (Real)2;
	}
	return hessian;
}

/** This method returns the partial derivative with respect to the prescribed dimension.*/
template< unsigned int Dim , unsigned int Degree , typename Real >
Polynomial< Dim , (Degree>1) ? Degree-1 : 0 , Real > Polynomial< Dim , Degree , Real >::d( int dim ) const
{
	Polynomial< Dim , (Degree>1) ? Degree-1 : 0 , Real > derivative;
	if( dim==0 ) 
	{
		for( int d=0 ; d<Degree ; d++ )
		{
			Real scale = (Real)(d+1);
			derivative._polynomials[d] = _polynomials[d+1];// * scale;
			derivative._polynomials[d] *= scale;
		}
	}
	else for( int d=0 ; d<Degree ; d++ ) derivative._polynomials[d] = _polynomials[d].d( dim-1 );
	return derivative;
}

template< unsigned int Dim , unsigned int Degree , typename Real >
template< unsigned int _Dim >
Polynomial< _Dim-1 , Degree , Real > Polynomial< Dim , Degree , Real >::operator()( Matrix< Real , _Dim , Dim > A ) const
{
	static_assert( _Dim!=0 , "Dimension cannot be negative" );
	return _pullBack< _Dim-1 >( A , Degree );
}

template< unsigned int Dim , unsigned int Degree , typename Real >
template< unsigned int _Dim >
Polynomial< _Dim-1 , Degree , Real > Polynomial< Dim , Degree , Real >::pullBack( Matrix< Real , _Dim , Dim > A ) const
{
	static_assert( _Dim!=0 , "Dimension cannot be negative" );
	return _pullBack< _Dim-1 >( A , Degree );
}

template< unsigned int Dim , unsigned int Degree , typename Real >
Polynomial< 1 , Degree , Real > Polynomial< Dim , Degree , Real >::operator()( const Ray< Real , Dim > & ray ) const
{
	Matrix< Real , 2 , Dim > A;
	for( unsigned int i=0 ; i<Dim ; i++ ) A(0,i) = ray.direction[i] , A(1,i) = ray.position[i];
	return _pullBack< 1 >( A , Degree );
}

template< unsigned int Dim , unsigned int Degree , typename Real >
Real Polynomial< Dim , Degree , Real >::integrateUnitCube( void ) const
{
	// I_d = \int_0^1 ... \int_0^1 x_n^d * P_d(x_1,...,x_{n-1}) dx_n ... dx_1
	//     = 1/(d+1) * \int_0^1 ... \int_0^1 P_d(x_1,...,x_{n-1}) dx_{n-1} ... dx_1
	Real integral = 0;
	for( int d=0 ; d<=Degree ; d++ ) integral += 1./(d+1) * _polynomials[d].integrateUnitCube();
	return integral;
}

template< unsigned int Dim , unsigned int Degree , typename Real >
Real Polynomial< Dim , Degree , Real >::integrateUnitRightSimplex( void ) const
{
	// I_d = \int_0^1 \int_0^{1-x_1} ... \int_0^{1-x_1-x_2...-x_{n-1}} x_n^d * P_d(x_1,...,x_{n-1}) dx_n ... dx_1
	//     = 1/(d+1) * \int_0^1 ... \int_0^1 P_d(x_1,...,x_{n-1}) * (1 - x_1 - x_2 - ... - x_{n-1} )^{d+1} dx_{n-1} ... dx_1
	Real integral = 0;
	Polynomial< Dim-1 , Degree+1 , Real > p;
	Polynomial< Dim-1 , 1 , Real > _p;
	{
		unsigned int indices[ Dim>1 ? Dim-1 : 1 ];
		for( int d=0 ; d<Dim-1 ; d++ ) indices[d] = 0;
		_p.coefficient( indices ) = 1;
		for( int d=0 ; d<Dim-1 ; d++ )
		{
			indices[d] = 1;
			_p.coefficient( indices ) = -1;
			indices[d] = 0;
		}
	}
	 p = _p;
	for( int d=0 ; d<=Degree ; d++ )
	{
		integral += ( _polynomials[d] * p ).integrateUnitRightSimplex() / (Real)( d+1 );
		p = p * _p;
	}
	return integral;
}

template< unsigned int Dim , unsigned int Degree , typename Real >
Matrix< Real , Polynomial< Dim , Degree , Real >::NumCoefficients , Polynomial< Dim , Degree , Real >::NumCoefficients > Polynomial< Dim , Degree , Real >::EvaluationMatrix( const Point< Real , Dim > positions[NumCoefficients] )
{
	SquareMatrix< Real , NumCoefficients > E;
	unsigned int degrees[ Dim ];
	for( unsigned int i=0 ; i<NumCoefficients ; i++ ) for( unsigned int j=0 ; j<NumCoefficients ; j++ )
	{
		SetDegrees( i , degrees );
		Real value = 1;
		for( int d=0 ; d<Dim ; d++ ) value *= pow( positions[j][d] , (Real)(double)degrees[d] );
		E(i,j) = value;
	}
	return E;
}

template< unsigned int Dim , unsigned int Degree , typename Real >
bool Polynomial< Dim , Degree , Real >::_print( std::ostream &ostream , const std::string varNames[] , bool first ) const
{
	first = _polynomials[0]._print( ostream , varNames , first );
	for( int d=1 ; d<=Degree ; d++ )
		if( d==1 ) first = _polynomials[d]._print( ostream , varNames , varNames[Dim-1] , first );
		else       first = _polynomials[d]._print( ostream , varNames , varNames[Dim-1] + "^" + std::to_string( d ) , first );
	return first;
}
template< unsigned int Dim , unsigned int Degree , typename Real >
bool Polynomial< Dim , Degree , Real >::_print( std::ostream &ostream , const std::string varNames[] , std::string suffix , bool first ) const
{
	first = _polynomials[0]._print( ostream , varNames , suffix , first );
	for( int d=1 ; d<=Degree ; d++ )
		if( d==1 ) first = _polynomials[d]._print( ostream , varNames , varNames[Dim-1] + "*" + suffix , first );
		else       first = _polynomials[d]._print( ostream , varNames , varNames[Dim-1] + "^" + std::to_string( d ) + "*" + suffix , first );
	return first;
}

template< unsigned int Dim , unsigned int Degree , typename Real >
std::ostream &operator << ( std::ostream &stream , const Polynomial< Dim , Degree , Real > &poly )
{
	std::string varNames[Dim];
	if      constexpr ( Dim==1 ) varNames[0] = "x";
	else if constexpr ( Dim==2 ) varNames[0] = "y" , varNames[1] = "x";
	else if constexpr ( Dim==3 ) varNames[0] = "z" , varNames[1] = "y" , varNames[2] = "x";
	else for( int i=0 ; i<Dim ; i++ ) varNames[i] = std::string( "x" ) + std::string( "_" ) + std::to_string( Dim-i-1 );
	if( poly._print( stream , varNames , true ) ) stream << "0";
	return stream;
}

template< unsigned int Degree , typename Real >
std::ostream &operator << ( std::ostream &stream , const Polynomial< 0 , Degree , Real > &poly )
{
	std::string varNames[1] = { "" };
	if( poly._print( stream , varNames , true ) ) stream << "0";
	return stream;
}

template< unsigned int Dim , unsigned int Degree , typename Real >
void Polynomial< Dim , Degree , Real >::Scale( Real s )
{
	for( int d=0 ; d<=Degree ; d++ ) _polynomials[d] *= s;
}

template< unsigned int Dim , unsigned int Degree , typename Real >
void Polynomial< Dim , Degree , Real >::Add( const Polynomial< Dim , Degree , Real > &p )
{
	for( int d=0 ; d<=Degree ; d++ ) _polynomials[d] += p._polynomials[d];
}

template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 , typename Real >
Polynomial< Dim , Degree1 + Degree2 , Real > operator * ( const Polynomial< Dim , Degree1 , Real > &p1 , const Polynomial< Dim , Degree2 , Real > &p2 )
{
	Polynomial< Dim , Degree1 + Degree2 , Real > p;
	for( int d1=0 ; d1<=Degree1 ; d1++ ) for( int d2=0 ; d2<=Degree2 ; d2++ ) p._polynomials[ d1+d2 ] += p1._polynomials[d1] * p2._polynomials[d2];
	return p;
}

template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 , typename Real >
Polynomial< Dim , Max< Degree1 , Degree2 >::Value , Real > operator + ( const Polynomial< Dim , Degree1 , Real > &p1 , const Polynomial< Dim , Degree2 , Real > &p2 )
{
	Polynomial< Dim , Max< Degree1 , Degree2 >::Value , Real > p;
	for( int d=0 ; d<=Degree1 ; d++ ) p._polynomials[d] += p1._polynomials[d];
	for( int d=0 ; d<=Degree2 ; d++ ) p._polynomials[d] += p2._polynomials[d];
	return p;
}

template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 , typename Real >
Polynomial< Dim , Max< Degree1 , Degree2 >::Value , Real > operator - ( const Polynomial< Dim , Degree1 , Real > &p1 , const Polynomial< Dim , Degree2 , Real > &p2 ){ return p1 + (-p2); }

template< unsigned int Degree , typename Real >
unsigned int Roots( const Polynomial< 1 , Degree , Real > &p , Real *r , double eps )
{
	MK_ERROR_OUT( "Root functionality not supported for polynomial of degree = %d" , Degree );
	return 0;
}

template<>
inline unsigned int Roots( const Polynomial< 1 , 1 , double > &p , double *r , double eps )
{
	if( fabs( p.coefficient(0u) )<eps )
	{
		r[0] = 0;
		return 1;
	}
	else if( p.coefficient(1u)==0 ) return 0;
	else
	{
		r[0] = -p.coefficient(0u) / p.coefficient(1u);
		return 1;
	}
}

template<>
inline unsigned int Roots( const Polynomial< 1 , 2 , double > &p , double *r , double eps )
{
	if( fabs( p.coefficient(0u) )<eps )
	{
		Polynomial< 1 , 1 , double > _p;
		for( unsigned int i=0 ; i<=1 ; i++ ) _p.coefficient(i) = p.coefficient(i+1);
		r[0] = 0;
		return Roots( _p , r+1 , eps )+1;
	}
	else if( !p.coefficient(2u) ) return Roots( Polynomial< 1 , 1 , double >( p ) , r , eps );
	double disc = p.coefficient(1u)*p.coefficient(1u) - 4. * p.coefficient(0u) * p.coefficient(2u);
	if( disc<0 ) return 0;
	else if( disc==0 )
	{
		r[0] = - p.coefficient(1u) / ( 2 * p.coefficient(2u) );
		return 1;
	}
	else
	{
		disc = sqrt(disc);
		r[0] = ( -p.coefficient(1u) - disc ) / (2 * p.coefficient(2u) );
		r[1] = ( -p.coefficient(1u) + disc ) / (2 * p.coefficient(2u) );
		return 2;
	}
}

template<>
inline unsigned int Roots( const Polynomial< 1 , 3 , double > &p , double *r , double eps )
{
	if( fabs( p.coefficient(0u) )<eps )
	{
		Polynomial< 1 , 2 , double > _p;
		for( unsigned int i=0 ; i<=2 ; i++ ) _p.coefficient(i) = p.coefficient(i+1);
		r[0] = 0;
		return Roots( _p , r+1 , eps )+1;
	}
	else if( !p.coefficient(3u) ) return Roots( Polynomial< 1 , 2 , double >( p ) , r , eps );
	return SergeyKhashin::Poly34::SolveP3( r , p.coefficient(2u)/p.coefficient(3u) , p.coefficient(1u)/p.coefficient(3u) , p.coefficient(0u)/p.coefficient(3u) , eps );
}

template<>
inline unsigned int Roots( const Polynomial< 1 , 4 , double > &p , double *r , double eps )
{
	if( fabs( p.coefficient(0u) )<eps )
	{
		Polynomial< 1 , 3 , double > _p;
		for( unsigned int i=0 ; i<=3 ; i++ ) _p.coefficient(i) = p.coefficient(i+1);
		r[0] = 0;
		return Roots( _p , r+1 , eps )+1;
	}
	else if( !p.coefficient(4u) ) return Roots( Polynomial< 1 , 3 , double >( p ) , r , eps );
	return SergeyKhashin::Poly34::SolveP4( r , p.coefficient(3u)/p.coefficient(4u) , p.coefficient(2u)/p.coefficient(4u) , p.coefficient(1u)/p.coefficient(4u) , p.coefficient(0u)/p.coefficient(4u) , eps );
}

template<>
inline unsigned int Roots( const Polynomial< 1 , 5 , double > &p , double *r , double eps )
{
	if( fabs( p.coefficient(0u) )<eps )
	{
		Polynomial< 1 , 4 , double > _p;
		for( unsigned int i=0 ; i<=4 ; i++ ) _p.coefficient(i) = p.coefficient(i+1);
		r[0] = 0;
		return Roots( _p , r+1 , eps )+1;
	}
	else if( !p.coefficient(5u) ) return Roots( Polynomial< 1 , 4 , double >( p ) , r , eps );
	return SergeyKhashin::Poly34::SolveP5( r , p.coefficient(4u)/p.coefficient(5u) , p.coefficient(3u)/p.coefficient(5u) , p.coefficient(2u)/p.coefficient(5u) , p.coefficient(1u)/p.coefficient(5u) , p.coefficient(0u)/p.coefficient(5u) , eps );
}
