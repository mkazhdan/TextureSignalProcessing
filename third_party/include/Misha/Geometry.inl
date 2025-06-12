/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
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

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif // M_PI


inline long long HalfEdgeKey( int i1 , int i2 )
{
	return ( ( (long long) i1 )<<32 ) | ( (long long) i2 );
}
inline long long EdgeKey( int i1 , int i2 )
{
	if( i1>i2 ) return HalfEdgeKey( i1 , i2 );
	else		return HalfEdgeKey( i2 , i1 );
}
inline void FactorEdgeKey( long long key , int& idx1 , int& idx2 )
{
	long long i1 , i2;
	i1 = key>>32;
	i2 = (key<<32)>>32;
	idx1 = int( i1 );
	idx2 = int( i2 );
}
inline long long OppositeHalfEdgeKey( long long key )
{
	int i1 , i2;
	FactorEdgeKey( key , i1 , i2 );
	return HalfEdgeKey( i2 , i1 );
}

#if 0
#if !FAST_POINT
/////////////
// Point3D //
/////////////
template< class Real >
Point3D<Real> Point3D<Real>::CrossProduct( const Point3D<Real>& p1 , const Point3D<Real> & p2 )
{
	Point3D<Real> p;
	p.coords[0]= p1.coords[1]*p2.coords[2]-p1.coords[2]*p2.coords[1];
	p.coords[1]=-p1.coords[0]*p2.coords[2]+p1.coords[2]*p2.coords[0];
	p.coords[2]= p1.coords[0]*p2.coords[1]-p1.coords[1]*p2.coords[0];
	return p;
}
#endif // FAST_POINT
#endif

////////////
// Matrix //
////////////
template<class Real,int Cols,int Rows>
void Matrix<Real,Cols,Rows>::Add(const Matrix<Real,Cols,Rows>& m)
{
	for(int i=0;i<Cols;i++)	for(int j=0;j<Rows;j++)	coords[i][j]+=m.coords[i][j];
}
template<class Real,int Cols,int Rows>
void Matrix<Real,Cols,Rows>::Scale(Real s)
{
	for(int i=0;i<Cols;i++)	for(int j=0;j<Rows;j++)	coords[i][j]*=s;
}
template<class Real,int Cols,int Rows>
Real Matrix<Real,Cols,Rows>::InnerProduct(const Matrix<Real,Cols,Rows>& m) const
{
	Real dot=0;
	for(int i=0;i<Cols;i++)
		for(int j=0;j<Rows;j++)
			dot+=m.coords[i][j]*coords[i][j];
	return dot;
}
template<class Real,int Cols,int Rows>
template<int Cols1>
Matrix<Real,Cols1,Rows> Matrix<Real,Cols,Rows>::operator * (const Matrix<Real,Cols1,Cols>& m) const
{
	Matrix<Real,Cols1,Rows> n;
	for(int i=0;i<Cols1;i++)
		for(int j=0;j<Rows;j++)
			for(int k=0;k<Cols;k++)
				n.coords[i][j]+=m.coords[i][k]*coords[k][j];
	return n;
}

template< typename Real , int Cols , int Rows >
template< typename T >
Point< T , Rows , Real > Matrix< Real , Cols , Rows >::operator * ( const Point< T , Cols , Real >& v ) const
{
	Point< T , Rows , Real > out;
	for( int j=0 ; j<Cols ; j++ )
	{
		const Real* _coords = coords[j];
		T _v = v.coords[j];
		for( int i=0 ; i<Rows ; i++ ) out.coords[i] += _v * _coords[i];
	}
	return out;
}

template< typename Real , int Cols , int Rows >
template< typename T >
Point< T , Rows , Real > Matrix< Real , Cols , Rows >::operator () ( const Point< T , Cols , Real >& v) const { return (*this)*v; }

template< class Real , int Cols , int Rows >
template< typename T >
Real Matrix< Real , Cols , Rows >::operator () ( const Point< T , Rows , Real >& v1 , const Point< T , Cols , Real > &v2 ) const { return Point< T , Rows , Real >::Dot( v1 , (*this)*v2 ); }

template<class Real,int Cols,int Rows>
Matrix<Real,Rows,Cols> Matrix<Real,Cols,Rows>::transpose(void) const
{
	Matrix<Real,Rows,Cols> out;
	for(int i=0;i<Cols;i++)
		for(int j=0;j<Rows;j++)
			out.coords[j][i]=coords[i][j];
	return out;
}

//////////////////
// SquareMatrix //
//////////////////
template<> inline double Matrix< double , 1 , 1 >::determinant( void ) const { return coords[0][0];}
template<> inline double Matrix< double , 2 , 2 >::determinant( void ) const { return coords[0][0]*coords[1][1] - coords[0][1]*coords[1][0]; }
template<> inline double Matrix< double , 3 , 3 >::determinant( void ) const
{
	return
		coords[0][0]*( coords[1][1]*coords[2][2] - coords[2][1]*coords[1][2] ) +
		coords[1][0]*( coords[2][1]*coords[0][2] - coords[0][1]*coords[2][2] ) +
		coords[2][0]*( coords[0][1]*coords[1][2] - coords[0][2]*coords[1][1] ) ;
}
template< class Real , int Dim >
Real Matrix< Real , Dim , Dim >::subDeterminant( int c , int r ) const
{
	Matrix< Real , Dim-1 , Dim-1 > temp;
	for( int i=0 , ii=0 ; i<Dim ; i++ )
	{
		if( i==c ) continue;
		for( int j=0 , jj=0 ; j<Dim ; j++ )
		{
			if( j==r ) continue;
			temp.coords[ii][jj] = coords[i][j];
			jj++;
		}
		ii++;
	}
	return Real( temp.determinant() );
}

template< class Real , int Dim >
template< typename T >
Point< T , Dim , Real > Matrix< Real , Dim , Dim >::operator * ( const Point< T , Dim , Real >& v ) const
{
	Point< T , Dim , Real > out;
	for( int j=0 ; j<Cols ; j++ )
	{
		const Real* _coords = coords[j];
		T _v = v.coords[j];
		for( int i=0 ; i<Dim ; i++ ) out.coords[i] += _v * _coords[i];
	}
	return out;
}

template< class Real , int Dim >
template< typename T >
Point< T , Dim , Real > Matrix< Real , Dim , Dim >::operator () ( const Point< T , Dim , Real >& v ) const { return (*this)*v; }

template< class Real , int Dim >
template< typename T >
Real Matrix< Real , Dim , Dim >::operator () ( const Point< T , Dim , Real >& v1 , const Point< T , Dim , Real > &v2 ) const { return Point< T , Dim , Real >::Dot( v1 , (*this)*v2 ); }

template< class Real , int Dim >
void Matrix< Real , Dim , Dim >::Add(const Matrix< Real , Dim , Dim >& m)
{
	for(int i=0;i<Dim;i++)	for(int j=0;j<Dim;j++)	coords[i][j]+=m.coords[i][j];
}

template< class Real , int Dim >
void Matrix< Real , Dim , Dim >::Scale(Real s)
{
	for(int i=0;i<Dim;i++)	for(int j=0;j<Dim;j++)	coords[i][j]*=s;
}

template< class Real , int Dim >
Real Matrix< Real , Dim , Dim >::InnerProduct(const Matrix< Real , Dim , Dim >& m) const
{
	Real dot=0;
	for(int i=0;i<Dim;i++)
		for(int j=0;j<Dim;j++)
			dot+=m.coords[i][j]*coords[i][j];
	return dot;
}

template< class Real , int Dim >
template< int Cols1 >
Matrix< Real , Cols1 , Dim > Matrix< Real , Dim , Dim >::operator * ( const Matrix< Real , Cols1 , Dim > &m ) const
{
	Matrix< Real , Cols1 , Dim > n;
	for(int i=0;i<Cols1;i++)
		for(int j=0;j<Dim;j++)
			for(int k=0;k<Dim;k++)
				n.coords[i][j]+=m.coords[i][k]*coords[k][j];
	return n;
}

template< class Real , int Dim >
Matrix< Real , Dim , Dim > Matrix< Real , Dim , Dim >::transpose(void) const
{
	Matrix< Real , Dim , Dim > out;
	for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) out.coords[j][i] = coords[i][j];
	return out;
}

template< class Real , int Dim >
Real Matrix< Real , Dim , Dim >::determinant( void ) const
{
	Real det = Real(1);
	// Gaussian Elimination
	Matrix xForm , temp;
	xForm = (*this);
	for( int i=0 ; i<Dim ; i++ )
	{
		int p = i ; Real v = (Real)fabs( xForm(i,i) );
		for( int j=i+1 ; j<Dim ; j++ ) if( fabs( xForm(i,j) )>v ) p = j , v = (Real)fabs( xForm(i,j) );
		if( !v ) return Real(0);
		if( i!=p ) det *= -xForm(i,p);
		else       det *=  xForm(i,p);
		temp.SetIdentity();
		Real scl = Real(1)/xForm(i,p);
		for( int j=0 ; j<Dim ; j++ ) temp(p,j) = - xForm(i,j) * scl;
		temp(i,p) = Real(1);
		temp(p,p) = -xForm(i,i) * scl;
		temp(i,i) = Real(0);
		temp(p,i) = scl; // Note that this is last so that if p=i the value is the right one
		xForm = temp * xForm;
	}
	return det;
}
template< class Real , int Dim >
Real Matrix< Real , Dim , Dim >::trace( void ) const
{
	Real tr = (Real)0;
	for( int i=0 ; i<Dim ; i++ ) tr += coords[i][i];
	return tr;
}
template< class Real , int Dim >
Matrix< Real , Dim , Dim > Matrix< Real , Dim , Dim >::inverse( void ) const
{
	bool success;
	Matrix inv = inverse( success );
	if( !success ) fprintf( stderr , "[WARNING] Failed to invert matrix\n" );
	return inv;
}

template< class Real , int Dim >
Matrix< Real , Dim , Dim > Matrix< Real , Dim , Dim >::inverse( bool& success ) const
{
	// Gaussian Elimination
	Matrix xForm , iXForm , temp;
	iXForm.SetIdentity() , xForm = (*this);
	for( int i=0 ; i<Dim ; i++ )
	{
		int p = i ; Real v = (Real)fabs( xForm(i,i) );
		for( int j=i+1 ; j<Dim ; j++ ) if( fabs( xForm(i,j) )>v ) p = j , v = (Real)fabs( xForm(i,j) );
		if( v==(Real)0. )
		{
			//			fprintf( stderr , "[WARNING] Failed to invert matrix\n" );
			success = false;
			return Matrix();
		}
		// temp(i,j): mapping of the i-th row to the j-th row
		temp.SetIdentity();
		Real scl = Real(1)/xForm(i,p);
		for( int j=0 ; j<Dim ; j++ ) temp(p,j) = - xForm(i,j) * scl;
		temp(i,p) = Real(1);
		temp(p,p) = -xForm(i,i) * scl;
		temp(i,i) = Real(0);
		temp(p,i) = scl; // Note that this is last so that if p=i the value is the right one
		xForm = temp * xForm , iXForm = temp * iXForm;
	}
	success = true;
	return iXForm;
}
template< >
inline Matrix< float , 2 , 2 > Matrix< float , 2 , 2 >::inverse( bool& success ) const
{
	Matrix iXForm;
	float det = ( coords[0][0]*coords[1][1]-coords[0][1]*coords[1][0] );
	if( !det ) success = false;
	float d = 1.f / det;
	iXForm.coords[0][0] =  coords[1][1] * d;
	iXForm.coords[1][1] =  coords[0][0] * d;
	iXForm.coords[0][1] = -coords[0][1] * d;
	iXForm.coords[1][0] = -coords[1][0] * d;
	success = true;
	return iXForm;
}
template< >
inline Matrix< double , 2 , 2 > Matrix< double , 2 , 2 >::inverse( bool& success ) const
{
	Matrix iXForm;
	double det = ( coords[0][0]*coords[1][1]-coords[0][1]*coords[1][0] );
	if( !det ) success = false;
	double d = 1. / det;
	iXForm.coords[0][0] =  coords[1][1] * d;
	iXForm.coords[1][1] =  coords[0][0] * d;
	iXForm.coords[0][1] = -coords[0][1] * d;
	iXForm.coords[1][0] = -coords[1][0] * d;
	success = true;
	return iXForm;
}

template< class Real , int Dim >
Polynomial::Polynomial< 1 , Dim , Real , Real > Matrix< Real , Dim , Dim >::_characteristicPolynomial( Matrix< char , Dim , Dim > mask ) const
{
	if constexpr( Dim==1 )
	{
		if( mask.coords[0][0] ) return Polynomial::Polynomial< 1 , 1 , Real , Real >( Point< Real , 2 >( coords[0][0] , (Real)-1 ) );
		else                    return Polynomial::Polynomial< 1 , 1 , Real , Real >( Point< Real , 2 >( coords[0][0] , (Real) 0 ) );
	}
	else if constexpr( Dim==2 )
	{
		Polynomial::Polynomial< 1 , 1 , Real , Real > c[2][2];
		for( unsigned int i=0 ; i<2 ; i++ ) for( unsigned int j=0 ; j<2 ; j++ )
			if( mask.coords[i][j] ) c[i][j] = Polynomial::Polynomial< 1 , 1 , Real , Real >( Point< Real , 2 >( coords[i][j] , (Real)-1 ) );
			else                    c[i][j] = Polynomial::Polynomial< 1 , 1 , Real , Real >( Point< Real , 2 >( coords[i][j] , (Real) 0 ) );
		return c[0][0] * c[1][1] - c[0][1] * c[1][0];
	}
	else if constexpr( Dim==3 )
	{
		Polynomial::Polynomial< 1 , 1 , Real , Real > c[3][3];
		for( unsigned int i=0 ; i<3 ; i++ ) for( unsigned int j=0 ; j<3 ; j++ )
			if( mask.coords[i][j] ) c[i][j] = Polynomial::Polynomial< 1 , 1 , Real , Real >( Point< Real , 2 >( coords[i][j] , (Real)-1 ) );
			else                    c[i][j] = Polynomial::Polynomial< 1 , 1 , Real , Real >( Point< Real , 2 >( coords[i][j] , (Real) 0 ) );
		return
			c[0][0]*( c[1][1]*c[2][2] - c[2][1]*c[1][2] ) +
			c[1][0]*( c[2][1]*c[0][2] - c[0][1]*c[2][2] ) +
			c[2][0]*( c[0][1]*c[1][2] - c[0][2]*c[1][1] ) ;
	}
	else
	{
		Polynomial::Polynomial< 1 , Dim , Real , Real > cPoly;

		Matrix< Real , Dim-1 , Dim-1 > temp;
		Matrix< char , Dim-1 , Dim-1 > _mask;
		for( int c=0 ; c<Dim ; c++ )
		{
			for( int i=0 , ii=0 ; i<Dim ; i++ )
			{
				if( i==c ) continue;
				for( int j=1 ; j<Dim ; j++ )
				{
					temp.coords[ii][j-1] = coords[i][j];
					_mask.coords[ii][j-1] = mask.coords[i][j];
				}
				ii++;
			}
			Real sign = c&1 ? (Real)-1 : (Real)1;
			if( mask.coords[c][0] ) cPoly += temp._characteristicPolynomial( _mask ) * Polynomial::Polynomial< 1 , 1 , Real , Real >( Point< Real , 2 >( coords[c][0] , (Real)-1) ) * sign;
			else                    cPoly += temp._characteristicPolynomial( _mask ) * Polynomial::Polynomial< 1 , 1 , Real , Real >( Point< Real , 2 >( coords[c][0] , (Real) 0) ) * sign;
		}
		return cPoly;
	}
}

template< class Real , int Dim >
Polynomial::Polynomial< 1 , Dim , Real , Real > Matrix< Real , Dim , Dim >::characteristicPolynomial( void ) const
{
	Matrix< char , Dim , Dim > mask;
	for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) mask.coords[i][j] = i==j ? 1 : 0;
	return _characteristicPolynomial( mask );
}


template<class Real,int Dim>
void Matrix< Real , Dim , Dim >::Multiply( const Matrix< Real , Dim , Dim > &m )
{
	Matrix temp=*this;
	for(int i=0;i<Dim;i++)
		for(int j=0;j<Dim;j++)
		{
			this->coords[i][j]=0;
			for(int k=0;k<Dim;k++)	this->coords[i][j]+=temp.coords[k][j]*m.coords[i][k];
		}
}
template<class Real,int Dim>
void Matrix< Real , Dim , Dim >::SetIdentity(void)
{
	memset(this->coords,0,sizeof(Real)*Dim*Dim);
	for(int i=0;i<Dim;i++)	this->coords[i][i]=1;
}
template<class Real,int Dim>
template<class Real2>
Point<Real2,Dim-1> Matrix< Real , Dim , Dim >::operator () (const Point<Real2,Dim-1>& v) const
{
	Real2 scale=1;
	Point<Real2,Dim-1> out;
	for(int i=0;i<Dim-1;i++)
	{
		for( int j=0 ; j<Dim-1 ; j++ ) out.coords[i] += v.coords[j]*Real2( this->coords[j][i] );
		out.coords[i] += Real2( this->coords[Dim-1][i] );
		scale += Real2( this->coords[i][Dim-1] );
	}
	for(int i=0;i<Dim-1;i++)	out.coords[i]/=scale;
	return out;
}

template<class Real>
//Real Random(void){return Real(rand())/RAND_MAX;}
Real Random( void )
{
	static const unsigned long long ULL_RAND_MAX = ( unsigned long long )( RAND_MAX ) + 1;
	static const double RAND_MAX_SQUARED = double( ULL_RAND_MAX * ULL_RAND_MAX - 1 );
	unsigned long long r1 = rand() , r2 = rand();
	long long foo = r1 * ULL_RAND_MAX + r2;
	return Real( double(foo) / RAND_MAX_SQUARED );
}

template<class Real>
Real Random2( void )
{
	long long temp= (long long) ( rand() )*RAND_MAX+rand();
	return Real( (double(temp)/((size_t)RAND_MAX+1))/((size_t)RAND_MAX+1) );
}

template< typename Real , unsigned int Dim >
Point< Real , Dim > RandomBallPoint( void )
{
	Point< Real , Dim > p;
	while( true )
	{
		for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = (Real)( 1.0 - 2.0*Random2< Real >() );
		double l = p.squareNorm();
		if( l<=1 ) return p;
	}
}

template< typename Real , unsigned int Dim >
Point< Real , Dim > RandomSimplexPoint( void )
{
	auto InSimplex = []( Point< double , Dim > p )
		{
			double sum = 0;
			for( unsigned int d=0 ; d<Dim ; d++ ) sum += p[d];
			return sum<1.;
		};
	Point< Real , Dim > p;
	while( true )
	{
		for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = Random< double >();
		if( InSimplex( p ) ) return p;
	}
}

template< typename Real , unsigned int K >
Point< Real , K+1 > RandomBarycentricCoordinates( void )
{
	Point< Real , K+1 > bc;
	Point< Real , K > s = RandomSimplexPoint< Real , K >();
	bc[0] = (Real)1.;
	for( unsigned int k=0 ; k<K ; k++ ) bc[0] -= s[k] , bc[k+1] = s[k];
	return bc;
}


template< typename Real , unsigned int Dim >
Point< Real , Dim > RandomSpherePoint( void )
{
	Point< Real , Dim > p = RandomBallPoint< Real , Dim >();
	Real l = (Real)sqrt( p.squareNorm() );
	return p / l;
}

template< typename Real , unsigned int Dim >
SquareMatrix< double , Dim > RandomRotationMatrix( void )
{
	Point< double , Dim > frame[Dim];
	for( unsigned int d=0 ; d<Dim-1 ; d++ )
	{
		frame[d] = RandomSpherePoint< double , Dim >();
		while( true )
		{
			for( unsigned int dd=0 ; dd<d ; dd++ ) frame[d] -= Point< double , Dim >::Dot( frame[dd] , frame[d] ) * frame[dd];
			if( frame[d].squareNorm()>1e-10 )
			{
				frame[d] /= sqrt( frame[d].squareNorm() );
				break;
			}
		}
	}
	frame[Dim-1] = Point< double , Dim >::CrossProduct( frame );
	SquareMatrix< double , Dim > R;
	for( unsigned int i=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j<Dim ; j++ ) R(i,j) = frame[i][j];
	return R;
}


template<class Real>
XForm3x3<Real> RotationMatrix( const Point3D<Real>& axis , const Real& angle )
{
	double a = cos( angle / 2 );
	double b , c , d;
	Point3D< Real > ax = axis * Real(sin( angle / 2 ) / Point3D< Real >::Length( axis ));
	b = ax[0] , c = ax[1] , d = ax[2];
	return RotationMatrix< Real >( Real( a ) , Real( b ) , Real( c ) , Real( d ) );
}
template<class Real>
XForm3x3<Real> RotationMatrix( Real a , Real b , Real c , Real d )
{
	XForm3x3< Real > rot;
	rot( 0 , 0 ) = 1 - 2*c*c - 2*d*d;
	rot( 1 , 0 ) = 2*b*c - 2*a*d;
	rot( 2 , 0 ) = 2*b*d + 2*a*c;
	rot( 0 , 1 ) = 2*b*c + 2*a*d;
	rot( 1 , 1 ) = 1 - 2*b*b - 2*d*d;
	rot( 2 , 1 ) = 2*c*d - 2*a*b;
	rot( 0 , 2 ) = 2*b*d - 2*a*c;
	rot( 1 , 2 ) = 2*c*d + 2*a*b;
	rot( 2 , 2 ) = 1 - 2*b*b - 2*c*c;
	return rot;
}

template< class Real >
Point3D< Real > RandomTrianglePoint( const Point3D< Real >& v1 , const Point3D< Real >& v2 , const Point3D< Real >& v3 )
{
	// pick random point in triangle
	Real r1 = (Real)sqrt( Random< Real >() ) , r2 = Random< Real >();
	return v1 * (1-r1) + v2 * r1 * (1-r2) + v3 * r1 * r2;
}

////////////////////////////////////////////////////////////////////////////////
/*! Uses a barycentric coordinate vector to interpolate three data values
//  @param[in]  coords      bary centric coordinates
//  @param[in]  d0, d1, d2  data values to interpolate
//  @return     interpolated data value
*///////////////////////////////////////////////////////////////////////////////
template<typename BaryCoords, typename DataType>
inline DataType BarycentricInterpolate(const BaryCoords &coords
	, const DataType &d0, const DataType &d1, const DataType &d2)
{
	// Use barycentric coordinates normalized w/ L1 norm
	return (coords[0] * d0  + coords[1] * d1 + coords[2] * d2)
		/ (coords[0] + coords[1] + coords[2]);
}

////////////////////////////////////////////////////////////////////////////////
/*! Computes a triangle's circumscribed circle
//  http://en.wikipedia.org/wiki/Circumscribed_circle
//  @param[in]  p0, p1, p2      triangle vertex positions
//  @param[in]  tri             triangle to process
//  @param[out] center          incircle center
*///////////////////////////////////////////////////////////////////////////////
template<typename PointType, typename Real>
inline void Circumcircle(const PointType &p0, const PointType &p1
	, const PointType &p2, PointType &center, Real &radius)
{
	Point3D<Real> e[3];
	e[0] = Point3D<Real>(p2 - p1);
	e[1] = Point3D<Real>(p0 - p2);
	e[2] = Point3D<Real>(p1 - p0);
	Real a2 = SquareLength(e[0]);
	Real b2 = SquareLength(e[1]);
	Real c2 = SquareLength(e[2]);
	Real a = sqrt(a2);
	Real b = sqrt(b2);
	Real c = sqrt(c2);
	Real doubleA = Length(Point3D<Real>::CrossProduct(e[0], e[1]));
	// Radius =  (a * b * c) / (4A)
	// (a, b, and c are edge lengths, A is area)
	radius = (a * b * c) / (2 * doubleA);
	// Circumcenter Barycentric Coordinates:
	//  (a^2 (b^2 + c^2 - a^2), b^2 (c^2 + a^2 - b^2), c^2 (a^2 + b^2 - c^2))
	Point3D<Real> centerBaryCoords(a2 * (b2 + c2 - a2), b2 * (c2 + a2 - b2)
		, c2 * (a2 + b2 - c2));
	center = BarycentricInterpolate(centerBaryCoords, p0, p1, p2);
}

////////////////////////////////////////////////////////////////////////////////
/*! Computes a triangle's inscribed circle
//  http://en.wikipedia.org/wiki/Incircle
//  @param[in]  p0, p1, p2      triangle vertex positions
//  @param[out] center          incircle center
//  @param[out] radius          incircle radius
*///////////////////////////////////////////////////////////////////////////////
template<typename PointType, typename Real>
inline void Incircle(const PointType &p0, const PointType &p1
	, const PointType &p2, PointType &center, Real &radius)
{
	Point3D<Real> e[3];
	e[0] = Point3D<Real>(p2 - p1);
	e[1] = Point3D<Real>(p0 - p2);
	e[2] = Point3D<Real>(p1 - p0);
	Real a = Length(e[0]);
	Real b = Length(e[1]);
	Real c = Length(e[2]);
	Real doubleA = Length(Point3D<Real>::CrossProduct(e[0], e[1]));
	// Radius =  (2A) / (a + b + c)
	// (a, b, and c are edge lengths, A is area)
	radius = doubleA / (a + b + c);
	// Incenter Barycentric Coordinates: (a, b, c)
	Point3D<Real> centerBaryCoords(a, b, c);
	center = BarycentricInterpolate(centerBaryCoords, p0, p1, p2);
}

template< class Real >
Point3D< Real > NearestPointOnEdge( Point3D< Real > point , const Point3D< Real > edge[2] , Real& b0 , Real& b1 )
{
	Point3D< Real > d = edge[1] - edge[0] , p = point - edge[0];
	Real dot = Point3D< Real >::Dot( p , d );
	if( dot<0 ) 
	{
		b0 = 1.;
		return edge[0];
	}
	else if( dot>Point3D< Real >::SquareNorm( d ) ) { 
		b1 = 1.; 
		return edge[1];
	}
	else
	{
		// Solve for the minimizer of:
		//                            E(t) = || p - t*d ||^2
		//                                 = || p ||^2 - 2 t < p , d > + t^2 || d ||^2
		//            =>             0 = -< p , d > + t || d ||^2
		//            <=>    t = < p , d > / || d ||^2
		Real t = dot / Point3D< Real >::SquareNorm( d );
		b0 = (Real)( 1.-t ) , b1 = t;
		return edge[0] + d * t;
	}
}

template< class Real >
Point3D< Real > NearestPointOnTriangle( Point3D< Real > point , const Point3D< Real > triangle[3] , Real* b )
{

	b[0] = b[1] = b[2] = 0;
	Point3D< Real > d[] = { triangle[1]-triangle[0] , triangle[2]-triangle[0] } , p = point - triangle[0] , n = CrossProduct( d[0] , d[1] );

	if( !Length(n) ) return ( triangle[0] + triangle[1] + triangle[2] ) / (Real)3.;


	if     ( Point3D< Real >::Dot( point-triangle[0] , CrossProduct( n , triangle[1]-triangle[0] ) )<0 ){ Point3D< Real > edge[] = { triangle[0] , triangle[1] } ; return NearestPointOnEdge( point , edge , b[0] , b[1] ); }
	else if( Point3D< Real >::Dot( point-triangle[1] , CrossProduct( n , triangle[2]-triangle[1] ) )<0 ){ Point3D< Real > edge[] = { triangle[1] , triangle[2] } ; return NearestPointOnEdge( point , edge , b[1] , b[2] ); }
	else if( Point3D< Real >::Dot( point-triangle[2] , CrossProduct( n , triangle[0]-triangle[2] ) )<0 ){ Point3D< Real > edge[] = { triangle[2] , triangle[0] } ; return NearestPointOnEdge( point , edge , b[2] , b[0] ); }
	else
	{
		// Solve for the minimizer of:
		//                            E(s,t) = || p - s*d[0]-t*d[1] ||^2
		//                                   = || p ||^2 - 2 s < p , d[0] > - 2 t < p , d[1] > + 2 s t < d[0] , d[1] > + s^2 || d[0] ||^2 + t^2 || d[1] ||^2
		//   =>  (0,0) = ( -< p , d[0] > + t < d[0] , d[1] > + s || d[0] ||^2 , -< p , d[1] > + s < d[0] , d[1] > + t || d[1] ||^2
		//            <=> | < p , d[0] > | = | < d[0] , d[0] >   < d[0] , d[1] > | | s |
		//                | < p , d[1] > |   | < d[0] , d[1] >   < d[1] , d[1] > | | t |
		SquareMatrix< Real , 2 > M , M_inverse;
		M(0,0) = Point3D< Real >::SquareNorm( d[0] ) , M(1,0) = M(0,1) = Point3D< Real >::Dot( d[0] , d[1] ) , M(1,1) = Point3D< Real >::SquareNorm( d[1] );
		Real det = M(0,0)*M(1,1) - M(0,1)*M(1,0);
		M_inverse(0,0) = M(1,1) , M_inverse(0,1) = -M(0,1) , M_inverse(1,0) = -M(1,0) , M_inverse(1,1) = M(0,0);
		M_inverse /= det;
		Point2D< Real > st = M_inverse * Point2D< Real >( Point3D< Real >::Dot( p , d[0] ) , Point3D< Real >::Dot( p , d[1] ) );
		b[0] = (Real)( 1. - st[0] - st[1] ) , b[1] = (Real)( st[0] , b[2] = st[1] );
		Point3D< Real > ret = triangle[0] * ( Real )( 1. - st[0] - st[1] ) + d[0] * st[0] + d[1] * st[1];

		return triangle[0] * ( Real )( 1. - st[0] - st[1] ) + d[0] * st[0] + d[1] * st[1];
	}
}

/////////////
// Simplex //
/////////////
template< class Real , unsigned int Dim , unsigned int K >
void Simplex< Real , Dim , K >::split( Point< Real , Dim > pNormal , Real pOffset , std::vector< Simplex >& back , std::vector< Simplex >& front ) const
{
	Real values[K+1];
	for( unsigned int k=0 ; k<=K ; k++ ) values[k] = Point< Real , Dim >::Dot( p[k] , pNormal ) - pOffset;
	return split( values , back , front );
}
template< class Real , unsigned int Dim , unsigned int K >
void Simplex< Real , Dim , K >::split( const Real values[K+1] , std::vector< Simplex >& back , std::vector< Simplex >& front ) const
{
	bool frontSet = false , backSet = false;

	// Evaluate the hyper-plane's function at the vertices and mark if strictly front/back vertices have been found
	for( unsigned int k=0 ; k<=K ; k++ ) backSet |= ( values[k]<0 ) , frontSet |= ( values[k]>0 );

	// If all the vertices are behind or on, or all the vertices are in front or on, we are done.
	if( !frontSet ){ back.push_back( *this ) ; return; }
	if( !backSet ){ front.push_back( *this ) ; return; }

	// Pick some intersection of the hyper-plane with a simplex edge
	unsigned int v1 , v2;
	Point< Real , Dim > midPoint;
	{
		for( unsigned int i=0 ; i<K ; i++ ) for( unsigned int j=i+1 ; j<=K ; j++ ) if( values[i]*values[j]<0 )
		{
			v1 = i , v2 = j;
			Real t1 = values[i] / ( values[i] - values[j] ) , t2 = (Real)( 1. - t1 );
			midPoint = p[j]*t1 + p[i]*t2;
		}
	}
	// Iterate over each face of the simplex, split it with the hyper-plane and connect the sub-simplices to the mid-point
	for( unsigned int i=0 ; i<=K ; i++ )
	{
		if( i!=v1 && i!=v2 ) continue;
		Simplex< Real , Dim , K-1 > f;		// The face
		Simplex< Real , Dim , K > s;		// The sub-simplex
		Real v[K];
		for( unsigned int j=0 , idx=0 ; j<=K ; j++ ) if( j!=i ){ f[idx] = p[j] , v[idx] = values[j] ; idx++ ; }
		std::vector< Simplex< Real , Dim , K-1 > > _back , _front;
		f.split( v , _back , _front );
		s[i] = midPoint;

		for( unsigned int j=0 ; j<_back.size() ; j++ ) 
		{
			for( unsigned int k=0 ; k<K ; k++ ) s[ k<i ? k : k+1 ] = _back[j][k];
			back.push_back( s );
		}

		for( unsigned int j=0 ; j<_front.size() ; j++ ) 
		{
			for( unsigned int k=0 ; k<K ; k++ ) s[ k<i ? k : k+1 ] = _front[j][k];
			front.push_back( s );
		}
	}
}

//////////////////
// SimplexIndex //
//////////////////

template< unsigned int K , typename Index >
template< unsigned int _K , typename FaceFunctor /* = std::function< void ( SimplexIndex< _K , Index > )*/ >
void SimplexIndex< K , Index >::ProcessFaces( FaceFunctor F )
{
	SimplexIndex< K , Index > si;
	for( unsigned int k=0 ; k<=K ; k++ ) si[k] = k;
	si.template processFaces< _K >( F );
}

template< unsigned int K , typename Index >
template< unsigned int _K , typename FaceFunctor /* = std::function< void ( SimplexIndex< _K , Index > )*/ >
void SimplexIndex< K , Index >::processFaces( FaceFunctor F ) const
{
	static_assert( _K<=K , "[ERROR] Face dimension too high" );
	if constexpr( K==_K ) F( *this );
	else for( unsigned int k=0 ; k<=K ; k++ ) _processFaces< _K >( F , k );
}

template< unsigned int K , typename Index >
template< unsigned int _K , typename ... UInts , typename FaceFunctor /* = std::function< void ( SimplexIndex< _K , Index > )*/ >
void SimplexIndex< K , Index >::_processFaces( FaceFunctor F , unsigned int faceIndex , UInts ... faceIndices ) const
{
	if constexpr( K-_K==sizeof...(UInts)+1 ) F( face( faceIndex , faceIndices... ) );
	else for( unsigned int f=0 ; f<faceIndex ; f++ ) _processFaces< _K >( F , f , faceIndex , faceIndices ... );
}

template< unsigned int K , typename Index >
template< typename ... UInts >
SimplexIndex< K - (unsigned int)sizeof...( UInts ) - 1 , Index > SimplexIndex< K , Index >::Face( unsigned int faceIndex , UInts ... faceIndices )
{
	SimplexIndex< K , Index > si;
	for( unsigned int k=0 ; k<=K ; k++ ) si[k] = k;
	return si.face( faceIndex , faceIndices ...  );
}

template< unsigned int K , typename Index >
template< typename ... UInts >
SimplexIndex< K - (unsigned int)sizeof...( UInts ) - 1 > SimplexIndex< K , Index >::face( unsigned int faceIndex , UInts ... faceIndices ) const
{
	static_assert( sizeof...(UInts)<K , "[ERROR] Too many indices" );
	bool flagged[K+1];
	{
		const unsigned int idx[] = { faceIndex , faceIndices ... };
		for( unsigned int k=0 ; k<=K ; k++ ) flagged[k] = false;
		for( unsigned int i=0 ; i<=sizeof...(faceIndices) ; i++ ) flagged[ idx[i] ] = true;
	}
	SimplexIndex< K - sizeof...( UInts ) - 1 > si;

	unsigned int idx=0;
	for( unsigned int k=0 ; k<=K ; k++ ) if( !flagged[k] ) si[idx++] = operator[]( k );

	return si;
}
template< unsigned int K , typename Index >
SimplexIndex< K-1 , Index > SimplexIndex< K , Index >::_Face( bool &oriented , unsigned int f )
{
	SimplexIndex< K-1 , Index > fi;
	unsigned int i=0;
	// Generate the face:
	//		(-1)^f * { 0 , 1 , ... , f-1 , f+1 , ... , K }
	for( unsigned int k=0 ; k<=K ; k++ ) if( k!=f ) fi[i++] = k;
	oriented = (f%2)==0;
	return fi;
}
template< unsigned int K , typename Index >
SimplexIndex< K-1 , Index > SimplexIndex< K , Index >::_face( bool &oriented , unsigned int f ) const
{
	SimplexIndex< K-1 , Index > s;
	unsigned int i=0;
	// Generate the face:
	//		(-1)^f * { 0 , 1 , ... , f-1 , f+1 , ... , K }
	for( unsigned int k=0 ; k<=K ; k++ ) if( k!=f ) s[i++] = idx[k];
	oriented = (f%2)==0;
	return s;
}

template< unsigned int K , typename Index >
template< typename Real , typename Vertex >
void SimplexIndex< K , Index >::split( const Real values[K+1] , std::vector< Vertex > &vertices , EdgeTable< Index > &edgeTable , std::vector< SimplexIndex >& back , std::vector< SimplexIndex >& front ) const
{
	bool frontSet = false , backSet = false;

	// Evaluate the hyper-plane's function at the vertices and mark if strictly front/back vertices have been found
	for( unsigned int k=0 ; k<=K ; k++ ) backSet |= ( values[k]<0 ) , frontSet |= ( values[k]>0 );
	// If all the vertices are behind or on, or all the vertices are in front or on, we are done.
	if( !frontSet ){ back.push_back( *this ) ; return; }
	if( !backSet ){ front.push_back( *this ) ; return; }

	// Pick some intersection of the hyper-plane with a simplex edge
	unsigned int v1=-1 , v2=-1;
	auto InitializationFunction = [&]( void )
		{
			Real t1 = values[v1] / ( values[v1] - values[v2] ) , t2 = (Real)( 1. - t1 );
			vertices.push_back( vertices[ idx[v2] ]*t1 + vertices[ idx[v1] ]*t2 );
			return (Index)( vertices.size()-1 );
		};
	for( unsigned int i=0 ; i<K ; i++ ) for( unsigned int j=i+1 ; j<=K ; j++ ) if( values[i]*values[j]<0 ) v1 = i , v2 = j;
	Index midPointIndex = edgeTable( idx[v1] , idx[v2] , InitializationFunction ); 

	// Iterate over each face of the simplex, split it with the hyper-plane and connect the sub-simplices to the mid-point
	for( unsigned int i=0 ; i<=K ; i++ )
	{
		if( i!=v1 && i!=v2 ) continue;
		Real _values[K];
		SimplexIndex< K-1 , Index > f;		// The face
		SimplexIndex< K , Index > s;		// The sub-simplex
		for( unsigned int j=0 , _idx=0 ; j<=K ; j++ ) if( j!=i ){ f[_idx] = idx[j] ; _values[_idx] = values[j] ; _idx++; }
		std::vector< SimplexIndex< K-1 , Index > > _back , _front;
		f.split( _values , vertices , edgeTable , _back , _front );
		s[i] = midPointIndex;

		for( unsigned int j=0 ; j<_back.size() ; j++ ) 
		{
			for( unsigned int k=0 ; k<K ; k++ ) s[ k<i ? k : k+1 ] = _back[j][k];
			back.push_back( s );
		}

		for( unsigned int j=0 ; j<_front.size() ; j++ ) 
		{
			for( unsigned int k=0 ; k<K ; k++ ) s[ k<i ? k : k+1 ] = _front[j][k];
			front.push_back( s );
		}
	}
}

template< unsigned int K , typename Index >
bool SimplexIndex< K , Index >::sort( void )
{
	unsigned int indices[K+1];
	for( unsigned int k=0 ; k<=K ; k++ ) indices[k] = k;
	return sort( indices );
}

template< unsigned int K , typename Index >
bool SimplexIndex< K , Index >::sort( const Index indices[] )
{
	// Find the permutation that orders the face indices from smallest to largest
	unsigned int permutation[ K+1 ] , temp[ K+1 ];
	for( unsigned int k=0 ; k<=K ; k++ ) permutation[k] = k;
	std::sort( permutation , permutation+K+1 , [&]( unsigned int i1 , unsigned int i2 ){ return indices[ idx[i1] ]<indices[ idx[i2] ]; } );

	for( int k=0 ; k<=K ; k++ ) temp[k] = idx[ permutation[k] ];
	for( int k=0 ; k<=K ; k++ ) idx[k] = temp[k];

	unsigned int count = 0;
	for( unsigned int i=0 ; i<=K ; i++ ) for( unsigned int j=0 ; j<i ; j++ ) if( permutation[i]>permutation[j] ) count++;

	return (count&1)==0;
}

//////////////////////////////
// MinimalAreaTriangulation //
//////////////////////////////
#ifdef NEW_MAT_CODE
template< typename Real , unsigned int Dim >
double MinimalAreaTriangulation::_Area( Point< Real , Dim > v0 , Point< Real , Dim > v1 , Point< Real , Dim > v2 )
{
	Point< Real , Dim > d[] = { v1-v0 , v2-v0 };
	SquareMatrix< double , 2 > M;
	for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) M(i,j) = Point< Real , Dim >::Dot( d[i] , d[j] );
	return sqrt( std::max< double >( 0 , M.determinant() ) ) / 2.;
}

template< typename Index , typename AreaFunctor /*=std::function< double (unsigned int , unsigned int ,  unsigned int ) > */ >
void MinimalAreaTriangulation::GetTriangulation( const AreaFunctor & AF , unsigned int vNum , std::vector< SimplexIndex< 2 , Index > >& triangles )
{
	static_assert( std::is_convertible_v< AreaFunctor , std::function< double (unsigned int,unsigned int,unsigned int) > > , "[ERROR] AreaFunctor poorly formed" );
	if( vNum<3 ) MK_ERROR_OUT( "Expected at least three vertices: " , vNum );

	triangles.resize( vNum-2 );

	if( vNum==3 )
	{
		triangles[0][0]=0;
		triangles[0][1]=1;
		triangles[0][2]=2;
		return;
	}
	else if( vNum==4 )
	{
		SimplexIndex< 2 , Index > tIndex[2][2];
		double area[] = { 0 , 0 };

		tIndex[0][0][0]=0;
		tIndex[0][0][1]=1;
		tIndex[0][0][2]=2;
		tIndex[0][1][0]=2;
		tIndex[0][1][1]=3;
		tIndex[0][1][2]=0;

		tIndex[1][0][0]=0;
		tIndex[1][0][1]=1;
		tIndex[1][0][2]=3;
		tIndex[1][1][0]=3;
		tIndex[1][1][1]=1;
		tIndex[1][1][2]=2;

		for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) area[i] += AF( tIndex[i][j][0] , tIndex[i][j][1] , tIndex[i][j][2] );

		if( area[0]>area[1] ) triangles[0] = tIndex[1][0] , triangles[1] = tIndex[1][1];
		else                  triangles[0] = tIndex[0][0] , triangles[1] = tIndex[0][1];
		return;
	}

	MinimalAreaTriangulation mat;
	unsigned int eCount = vNum;
	mat._bestTriangulation = new double[eCount*eCount];
	mat._midPoint = new unsigned int[eCount*eCount];
	for( unsigned int i=0 ; i<eCount*eCount ; i++ ) mat._bestTriangulation[i] = -1 , mat._midPoint[i] = -1;
	mat._getArea( 0 , 1 , AF , vNum );
	unsigned int idx = 0;
	mat._getTriangulation( 0 , 1 , AF , vNum , triangles , idx );
}

template< class Real , unsigned int Dim , typename Index >
void MinimalAreaTriangulation::GetTriangulation( const std::vector< Point< Real , Dim > >& vertices , std::vector< SimplexIndex< 2 , Index > >& triangles )
{
	return GetTriangulation( [&]( unsigned int i , unsigned int j , unsigned int k ){ return _Area( vertices[i] , vertices[j] , vertices[k] ); } , static_cast< unsigned int >( vertices.size() ) , triangles );
}

template< typename AreaFunctor /*=std::function< double (unsigned int , unsigned int ,  unsigned int ) > */ >
double MinimalAreaTriangulation::GetArea( const AreaFunctor & AF , unsigned int vNum )
{
	static_assert( std::is_convertible_v< AreaFunctor , std::function< double (unsigned int,unsigned int,unsigned int) > > , "[ERROR] AreaFunctor poorly formed" );
	if( vNum<3 ) MK_ERROR_OUT( "Expected at least three vertices: " , vNum );

	MinimalAreaTriangulation mat;
	unsigned int eCount = vNum;
	mat._bestTriangulation = new double[eCount*eCount];
	mat._midPoint = new unsigned int[eCount*eCount];
	for( int i=0 ; i<eCount*eCount ; i++ ) mat._bestTriangulation[i]=-1 , mat._midPoint[i] = -1;
	return mat._getArea( 0 , 1 , AF , vNum );
}

template< class Real , unsigned int Dim >
double MinimalAreaTriangulation::GetArea( const std::vector< Point< Real , Dim > > &vertices )
{
	return GetArea( [&]( unsigned int i , unsigned int j , unsigned int k ){ return _Area( vertices[i] , vertices[j] , vertices[k] ); } , vertices.size() );
}

template< typename Index , typename AreaFunctor /*=std::function< double (unsigned int , unsigned int ,  unsigned int ) > */ >
void MinimalAreaTriangulation::_getTriangulation( unsigned int i , unsigned int j , const AreaFunctor & AF , unsigned int vNum , std::vector< SimplexIndex< 2 , Index > > &triangles , unsigned int &idx )
{
	SimplexIndex< 2 , Index > tIndex;
	unsigned int eCount = vNum;
	unsigned int ii = i;
	if( i<j ) ii += eCount;
	if( j+1>=ii ) return;
	ii = _midPoint[i*eCount+j];
	if( ii>=0 )
	{
		tIndex[0] = (Index)i;
		tIndex[1] = (Index)j;
		tIndex[2] = (Index)ii;
		triangles[idx++] = tIndex;
		_getTriangulation( i , ii , AF , vNum , triangles , idx );
		_getTriangulation( ii , j , AF , vNum , triangles , idx );
	}
}

template< typename AreaFunctor /*=std::function< double (unsigned int , unsigned int ,  unsigned int ) > */ >
double MinimalAreaTriangulation::_getArea( unsigned int i , unsigned int j , const AreaFunctor & AF , unsigned int vNum )
{
	double a = std::numeric_limits< double >::infinity() , temp;
	unsigned int eCount = vNum;
	unsigned int idx = i*eCount+j;
	unsigned int ii = i;
	if( i<j ) ii += eCount;
	if( j+1>=ii)
	{
		_bestTriangulation[idx]=0;
		return 0;
	}
	unsigned int mid=-1;
	for( unsigned int r=j+1 ; r<ii ; r++ )
	{
		unsigned int rr=r%eCount;
		unsigned int idx1=i*eCount+rr , idx2=rr*eCount+j;

		temp = AF( rr , i , j );

		if( _bestTriangulation[idx1]>0 )
		{
			temp += _bestTriangulation[idx1];
			if( temp>a ) continue;
			if( _bestTriangulation[idx2]>0 ) temp += _bestTriangulation[idx2];
			else temp += _getArea( rr , j , AF , vNum );
		}
		else
		{
			if( _bestTriangulation[idx2]>0 ) temp += _bestTriangulation[idx2];
			else temp += _getArea( rr , j , AF , vNum );
			if( temp>a ) continue;
			temp += _getArea( i , rr , AF , vNum );
		}

		if( temp<a ) a = temp , mid = rr;
	}
	_bestTriangulation[idx]=a;
	_midPoint[idx]=mid;

	return a;
}

#else // !NEW_MAT_CODE
template< typename Real , unsigned int Dim >
double MinimalAreaTriangulation< Real , Dim >::_Area( Point< Real , Dim > v0 , Point< Real , Dim > v1 , Point< Real , Dim > v2 )
{
	Point< Real , Dim > d[] = { v1-v0 , v2-v0 };
	SquareMatrix< double , 2 > M;
	for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) M(i,j) = Point< Real , Dim >::Dot( d[i] , d[j] );
	return sqrt( std::max< double >( 0 , M.determinant() ) ) / 2.;
}

template< class Real , unsigned int Dim >
MinimalAreaTriangulation< Real , Dim >::MinimalAreaTriangulation( void ) : _bestTriangulation(nullptr) , _midPoint(nullptr) {}

template< class Real , unsigned int Dim >
MinimalAreaTriangulation< Real , Dim >::~MinimalAreaTriangulation( void )
{
	delete[] _bestTriangulation;
	delete[] _midPoint;
	_bestTriangulation = nullptr;
	_midPoint = nullptr;
}

template< class Real , unsigned int Dim >
template< typename Index >
void MinimalAreaTriangulation< Real , Dim >::GetTriangulation( const std::vector< Point< Real , Dim > >& vertices , std::vector< SimplexIndex< 2 , Index > >& triangles )
{
	triangles.resize( vertices.size() - 2 );
	if( vertices.size()==3 )
	{
		triangles[0][0]=0;
		triangles[0][1]=1;
		triangles[0][2]=2;
		return;
	}
	else if( vertices.size()==4 )
	{
		SimplexIndex< 2 , Index > tIndex[2][2];
		double area[] = { 0 , 0 };

		tIndex[0][0][0]=0;
		tIndex[0][0][1]=1;
		tIndex[0][0][2]=2;
		tIndex[0][1][0]=2;
		tIndex[0][1][1]=3;
		tIndex[0][1][2]=0;

		tIndex[1][0][0]=0;
		tIndex[1][0][1]=1;
		tIndex[1][0][2]=3;
		tIndex[1][1][0]=3;
		tIndex[1][1][1]=1;
		tIndex[1][1][2]=2;

		for( int i=0 ; i<2 ; i++ ) for( int j=0 ; j<2 ; j++ ) area[i] += _Area( vertices[ tIndex[i][j][0] ] , vertices[ tIndex[i][j][1] ] , vertices[ tIndex[i][j][2] ] );

		if( area[0]>area[1] ) triangles[0] = tIndex[1][0] , triangles[1] = tIndex[1][1];
		else                  triangles[0] = tIndex[0][0] , triangles[1] = tIndex[0][1];
		return;
	}

	MinimalAreaTriangulation mat;
	size_t eCount=vertices.size();
	mat._bestTriangulation = new double[eCount*eCount];
	mat._midPoint = new size_t[eCount*eCount];
	for( unsigned int i=0 ; i<eCount*eCount ; i++ ) mat._bestTriangulation[i] = -1 , mat._midPoint[i] = -1;
	mat._getArea( 0 , 1 , vertices );
	//	triangles.clear();
	size_t idx = 0;
	//	GetTriangulation(0,1,vertices,triangles);
	mat._getTriangulation( 0 , 1 , vertices , triangles , idx );
}

template< class Real , unsigned int Dim >
double MinimalAreaTriangulation< Real , Dim >::GetArea( const std::vector< Point< Real , Dim > > &vertices )
{
	MinimalAreaTriangulation mat;
	size_t eCount = vertices.size();
	mat._bestTriangulation = new double[eCount*eCount];
	mat._midPoint = new size_t[eCount*eCount];
	for( int i=0 ; i<eCount*eCount ; i++ ) mat._bestTriangulation[i]=-1 , mat._midPoint[i] = -1;
	return mat._getArea( 0 , 1 , vertices );
}

template< typename Real , unsigned int Dim >
template< typename Index >
void MinimalAreaTriangulation< Real , Dim >::_getTriangulation( size_t i , size_t j , const std::vector< Point< Real , Dim > > &vertices , std::vector< SimplexIndex< 2 , Index > > &triangles , size_t &idx )
{
	SimplexIndex< 2 , Index > tIndex;
	size_t eCount = vertices.size();
	size_t ii = i;
	if( i<j ) ii += eCount;
	if( j+1>=ii ) return;
	ii = _midPoint[i*eCount+j];
	if( ii>=0 )
	{
		tIndex[0] = (Index)i;
		tIndex[1] = (Index)j;
		tIndex[2] = (Index)ii;
		triangles[idx++] = tIndex;
		_getTriangulation( i , ii , vertices , triangles , idx );
		_getTriangulation( ii , j , vertices , triangles , idx );
	}
}

template< class Real , unsigned int Dim >
double MinimalAreaTriangulation< Real , Dim >::_getArea( size_t i , size_t j , const std::vector< Point< Real , Dim > > &vertices )
{
	double a = std::numeric_limits< double >::infinity() , temp;
	size_t eCount = vertices.size();
	size_t idx = i*eCount+j;
	size_t ii = i;
	if( i<j ) ii += eCount;
	if( j+1>=ii)
	{
		_bestTriangulation[idx]=0;
		return 0;
	}
	size_t mid=-1;
	for( size_t r=j+1 ; r<ii ; r++ )
	{
		size_t rr=r%eCount;
		size_t idx1=i*eCount+rr,idx2=rr*eCount+j;
		temp = _Area( vertices[rr] , vertices[i] , vertices[j] );

		if( _bestTriangulation[idx1]>0 )
		{
			temp += _bestTriangulation[idx1];
			if( temp>a ) continue;
			if( _bestTriangulation[idx2]>0 ) temp += _bestTriangulation[idx2];
			else temp += _getArea( rr , j , vertices );
		}
		else
		{
			if( _bestTriangulation[idx2]>0 ) temp += _bestTriangulation[idx2];
			else temp += _getArea( rr , j , vertices );
			if( temp>a ) continue;
			temp += _getArea( i , rr , vertices );
		}

		if( temp<a ) a = temp , mid = rr;
	}
	_bestTriangulation[idx]=a;
	_midPoint[idx]=mid;

	return a;
}
#endif // NEW_MAT_CODE

//////////////////////
// EarTriangulation //
//////////////////////

template< typename Index , typename Real >
void EarTriangulation::GetTriangulation( const std::vector< Point< Real , 2 > > &vertices , std::vector< SimplexIndex< 2 , Index > > &triangles )
{
	struct PolygonVertex
	{
		unsigned int idx;
		bool isEar;
		PolygonVertex *prev , *next;

		PolygonVertex( void ) : idx(-1) , isEar(false) , prev(NULL) , next(NULL) {}
	};

	std::vector< PolygonVertex > polygonVertices( vertices.size() );
	for( unsigned int i=0 ; i<vertices.size() ; i++ )
	{
		polygonVertices[i].idx = i;
		polygonVertices[i].prev = &polygonVertices[( i+vertices.size()-1 ) % vertices.size() ];
		polygonVertices[i].next = &polygonVertices[( i+vertices.size()+1 ) % vertices.size() ];
	}
	PolygonVertex *polygon = &polygonVertices[0];

	auto ProcessPolygon =[&]< typename F /*=std::function< void ( PolygonVertex * ) >*/ >( F f )
	{
		PolygonVertex *v = polygon;
		do
		{
			f(v);
			v = v->next;
		}
		while( v!=polygon );
	};

	auto PolygonSize = [&]( void )
		{
			unsigned int sz = 0;
			ProcessPolygon( [&]( PolygonVertex * ){ sz++; } );
			return sz;
		};

	// From [O'Rourke]

	auto RotateCCW90 = []( Point< double , 2 > v )
		{
			return Point< Real , 2 >( -v[1] , v[0] );
		};

	auto TotalArea2 = [&]( void )
		{
			Point< Real , 2 > center;
			ProcessPolygon( [&]( PolygonVertex *v ){ center += vertices[v->idx]; } );
			center /= PolygonSize();

			Real area = 0;
			auto A = [&]( PolygonVertex *v )
				{
					Point< Real , 2 > v1 = vertices[v->idx] - center , v2 = vertices[v->next->idx ] - center;
					area += Point< Real , 2 >::Dot( RotateCCW90( v1 ) , v2 );
				};
			ProcessPolygon( A );
			return area;
		};

	auto Area2 = [&]( unsigned int i0 , unsigned int i1 , unsigned int i2 )
		{
			Point< Real , 2 > v1 = vertices[i1] - vertices[i0];
			Point< Real , 2 > v2 = vertices[i2] - vertices[i0];
			return Point< Real , 2 >::Dot( RotateCCW90(v1) , v2 );
		};


	auto Left = [&]( unsigned int i0 , unsigned int i1 , unsigned int i2 )
		{
			return Area2(i0,i1,i2)>0;
		};

	auto LeftOn = [&]( unsigned int i0 , unsigned int i1 , unsigned int i2 )
		{
			return Area2(i0,i1,i2)>=0;
		};

	auto Collinear = [&]( unsigned int i0 , unsigned int i1 , unsigned int i2 )
		{
			return Area2(i0,i1,i2)==0;
		};

	auto IntersectProp = [&]( unsigned int a , unsigned int b , unsigned int c , unsigned int d )
		{
			if( Collinear(a,b,c) || Collinear(a,b,d) || Collinear(c,d,a) || Collinear(c,d,b) ) return false;
			else return ( Left(a,b,c) ^ Left(a,b,d) ) && ( Left(c,d,a) ^ Left(c,d,b) );
		};

	auto Between = [&]( unsigned int a , unsigned int b , unsigned int c )
		{
			if( !Collinear(a,b,c) ) return false;
			Point< double , 2 > _a = vertices[a] , _b = vertices[b] , _c = vertices[c];
			if( _a[0]!=_b[0] ) return ( _a[0]<=_c[0] && _c[0]<=_b[0] ) || ( _a[0]>=_c[0] && _c[0]>=_b[0] );
			else               return ( _a[1]<=_c[1] && _c[1]<=_b[1] ) || ( _a[1]>=_c[1] && _c[1]>=_b[1] );
		};

	auto Intersect = [&]( unsigned int a , unsigned int b , unsigned int c , unsigned int d )
		{
			return IntersectProp(a,b,c,d) || Between(a,b,c) || Between(a,b,d) || Between(c,d,a) || Between(c,d,b);
		};

	auto InCone = [&]( const PolygonVertex *a , const PolygonVertex *b )
		{
			const PolygonVertex *a0 = a->prev , *a2 = a->next;
			if( LeftOn( a->idx , a2->idx , a0->idx ) ) return    Left  ( a->idx , b->idx , a0->idx ) && Left  ( b->idx , a->idx , a2->idx );
			else                                       return !( LeftOn( a->idx , b->idx , a2->idx ) && LeftOn( b->idx , a->idx , a0->idx ) );
		};

	auto IsDiagonal = [&]( const PolygonVertex *a , const PolygonVertex *b )
		{
			PolygonVertex *c = polygon;
			do
			{
				PolygonVertex *n = c->next;
				if( c!=a && n!=a && c!=b && n!=b && Intersect( a->idx , b->idx , c->idx , n->idx ) ) return false;
				c = n;
			}
			while( c!=polygon );
			return true;
		};

	auto IsEar = [&]( const PolygonVertex *v )
		{
			return InCone( v->prev , v->next ) && InCone( v->next , v->prev ) && IsDiagonal( v->prev , v->next );
		};

	if( TotalArea2()<0 ) for( unsigned int i=0 ; i<polygonVertices.size() ; i++ ) std::swap( polygonVertices[i].prev , polygonVertices[i].next );

	std::vector< PolygonVertex * > earVertices;
	// Mark the ear vertices
	ProcessPolygon( [&]( PolygonVertex *v ){ if( ( v->isEar = IsEar( v ) ) ) earVertices.push_back(v); } );

	while( triangles.size()<vertices.size()-2 )
	{
		if( !earVertices.size() )
		{
			//			MK_WARN( "Expected an ear vertex: " , PolygonSize() );
			ProcessPolygon( [&]( PolygonVertex *v ){ if( ( v->isEar = IsEar( v ) ) ) earVertices.push_back(v); } );
			if( !earVertices.size() )
			{
				//				ProcessPolygon( [&]( PolygonVertex *v ){ std::cout << v->idx << " : " << vertices[ v->idx ] << std::endl; }  );
				//				ProcessPolygon( [&]( PolygonVertex *v ){ poly.push_back( v->idx ); }  );
				//				MK_ERROR_OUT( "Could not find ears" );
				//std::cout << polys.size() << std::endl;
				return;
			}
		}
#if 1
		auto EarQuality = [&]( PolygonVertex *v )
			{
				Point< Real , 2 > v1 = vertices[ v->next->idx ] - vertices[ v->idx ];
				Point< Real , 2 > v2 = vertices[ v->prev->idx ] - vertices[ v->idx ];
				v1 /= (Real)sqrt( Point< Real , 2 >::SquareNorm( v1 ) );
				v2 /= (Real)sqrt( Point< Real , 2 >::SquareNorm( v2 ) );
				return Point< Real , 2 >::Dot( v1 , v2 );
			};

		unsigned int idx = 0;
		for( unsigned int i=1 ; i<earVertices.size() ; i++ ) if( EarQuality( earVertices[i] ) > EarQuality( earVertices[idx] ) ) idx = i;
		PolygonVertex *ear = earVertices[idx];
		earVertices[idx] = earVertices.back();
		earVertices.pop_back();
#else
		PolygonVertex *ear = earVertices.back();
		earVertices.pop_back();
#endif

		{
			SimplexIndex< 2 , Index > si;
			si[0] = ear->idx;
			si[1] = ear->next->idx;
			si[2] = ear->prev->idx;
			triangles.push_back( si );
		}

		PolygonVertex *prev = ear->prev , *next = ear->next;

		prev->next = next;
		next->prev = prev;
		if( ear==polygon ) polygon = next;

		if( !prev->isEar ) if( ( prev->isEar=IsEar( prev ) ) ) earVertices.push_back( prev );
		if( !next->isEar ) if( ( next->isEar=IsEar( next ) ) ) earVertices.push_back( next );
	}
}