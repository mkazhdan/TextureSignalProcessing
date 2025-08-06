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

#ifndef GEOMETRY_INCLUDED
#define GEOMETRY_INCLUDED

#define NEW_MAT_CODE
#define NEW_GEOMETRY_CODE

//#include <concepts>
#include <cmath>
#include <cassert>
#include <complex>
#include <vector>
#include <unordered_map>
#include <initializer_list>
#include <cstring>
#include <algorithm>
#include <functional>
#ifndef _WIN32
#include <cstdlib>
#endif
#include <float.h>
#include <unordered_map>
#include "Algebra.h"
#include "Exceptions.h"
#include "Array.h"

namespace MishaK
{
	template< typename V > struct FieldOf{ using F = typename V::R; };
	template<> struct FieldOf< float >{ using F = float; };
	template<> struct FieldOf< double >{ using F = double; };


	// An empty type
	template< typename Real >
	struct EmptyVectorType
	{
		EmptyVectorType& operator += ( const EmptyVectorType& p ){ return *this; }
		EmptyVectorType& operator -= ( const EmptyVectorType& p ){ return *this; }
		EmptyVectorType& operator *= ( Real s )                  { return *this; }
		EmptyVectorType& operator /= ( Real s )                  { return *this; }
		EmptyVectorType  operator +  ( const EmptyVectorType& p ) const { EmptyVectorType _p = *this ; _p += p ; return _p; }
		EmptyVectorType  operator -  ( const EmptyVectorType& p ) const { EmptyVectorType _p = *this ; _p -= p ; return _p; }
		EmptyVectorType  operator *  ( Real s )                   const { EmptyVectorType _p = *this ; _p *= s ; return _p; }
		EmptyVectorType  operator /  ( Real s )                   const { EmptyVectorType _p = *this ; _p /= s ; return _p; }

		friend std::ostream &operator << ( std::ostream &os , const EmptyVectorType &v ){ return os; }
	};
	template< typename Real > EmptyVectorType< Real > operator * ( Real s , EmptyVectorType< Real > v ){ return v*s; }

	template< typename _Real , typename ... VectorTypes >
	struct VectorTypeUnion
	{
	protected:
		typedef std::tuple< VectorTypes ... > _VectorTuple;
	public:
		typedef _Real Real;
		template< unsigned int I > using VectorType = typename std::tuple_element< I , _VectorTuple >::type;
		template< unsigned int I >       VectorType< I >& get( void )       { return std::get< I >( _vectorTypeTuple ); }
		template< unsigned int I > const VectorType< I >& get( void ) const { return std::get< I >( _vectorTypeTuple ); }

		VectorTypeUnion& operator += ( const VectorTypeUnion& p ){ _add<0>( p ) ; return *this; }
		VectorTypeUnion& operator -= ( const VectorTypeUnion& p ){ _sub<0>( p ) ; return *this; }
		VectorTypeUnion& operator *= ( Real s )                  { _mul<0>( s ) ; return *this; }
		VectorTypeUnion& operator /= ( Real s )                  { _div<0>( s ) ; return *this; }
		VectorTypeUnion  operator +  ( const VectorTypeUnion& p ) const { VectorTypeUnion _p = *this ; _p += p ; return _p; }
		VectorTypeUnion  operator -  ( const VectorTypeUnion& p ) const { VectorTypeUnion _p = *this ; _p -= p ; return _p; }
		VectorTypeUnion  operator *  ( Real s )                   const { VectorTypeUnion _p = *this ; _p *= s ; return _p; }
		VectorTypeUnion  operator /  ( Real s )                   const { VectorTypeUnion _p = *this ; _p /= s ; return _p; }

		VectorTypeUnion( void ){}
		VectorTypeUnion( const VectorTypes & ... vectors ){ _set< 0 >( vectors ... ); }

		friend std::ostream &operator << ( std::ostream &os , const VectorTypeUnion &v )
		{
			os << "{ ";
			v._streamOut< 0 >( os );
			os << " }";
			return os;
		}
	protected:
		std::tuple< VectorTypes ... > _vectorTypeTuple;
		template< unsigned int I , typename _Vector , typename ... _Vectors > void _set( const _Vector &vector , const _Vectors & ... vectors ){ get< I >() = vector ; _set< I+1 >( vectors ... ); }
		template< unsigned int I , typename _Vector                         > void _set( const _Vector &vector                                ){ get< I >() = vector ;                             }
		template< unsigned int I > typename std::enable_if< I!=sizeof...(VectorTypes) >::type _add( const VectorTypeUnion& p ){ get<I>() += p.get<I>() ; _add< I+1 >( p ); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(VectorTypes) >::type _add( const VectorTypeUnion& p ){ }
		template< unsigned int I > typename std::enable_if< I!=sizeof...(VectorTypes) >::type _sub( const VectorTypeUnion& p ){ get<I>() -= p.get<I>() ; _sub< I+1 >( p ); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(VectorTypes) >::type _sub( const VectorTypeUnion& p ){ }
		template< unsigned int I > typename std::enable_if< I!=sizeof...(VectorTypes) >::type _mul( Real s ){ get<I>() *= s ; _mul< I+1 >( s ); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(VectorTypes) >::type _mul( Real s ){ }
		template< unsigned int I > typename std::enable_if< I!=sizeof...(VectorTypes) >::type _div( Real s ){ get<I>() /= s ; _div< I+1 >( s ); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(VectorTypes) >::type _div( Real s ){ }
		template< unsigned int I > typename std::enable_if< I!=sizeof...(VectorTypes) >::type _streamOut( std::ostream &os ) const { os << get<I>() ; if( I!=sizeof...(VectorTypes)-1 ) os << " , "; _streamOut< I+1 >( os ); }
		template< unsigned int I > typename std::enable_if< I==sizeof...(VectorTypes) >::type _streamOut( std::ostream &os ) const { }
	};
	template< typename Real , typename ... Vectors >
	VectorTypeUnion< Real , Vectors ... > operator * ( Real s , VectorTypeUnion< Real , Vectors ... > vu ){ return vu * s; }

	template< class Real > Real Random( void );

	template< typename Real , int Rows , int Cols > class Matrix;
	template< typename Real , int Dim > using SquareMatrix = Matrix< Real , Dim , Dim >;

	template< typename Real , unsigned int Dim > class XForm;


	template< typename T , unsigned int Dim , typename Real=T > struct Point;

	template< typename T , unsigned int Dim , typename Real >
	struct Point : public InnerProductSpace< Real , Point< T , Dim , Real > >
	{
		typedef InnerProductSpace< T , Point > IPS;

		template< class ... Points >
		static void _AddColumnVector( SquareMatrix< T , Dim >& x , int c , Point point , Points ... points )
		{
			for( int r=0 ; r<Dim ; r++ ) x( c , r ) = point[r];
			if constexpr( sizeof...(Points) ) _AddColumnVector( x , c+1 , points ... );
		}

		template< typename ... Ts >
		void _set( T t , Ts ... ts )
		{
			coords[ Dim-1-sizeof...(Ts) ] = t;
			if constexpr( sizeof...(Ts) ) _set( ts... );
		}
	public:
		/////////////////////////////////
		// Inner product space methods //
		void Add( const Point &p ){ for( int d=0 ; d<Dim ; d++ ) coords[d] += p.coords[d]; }
		void Scale( Real s ){ for( int d=0 ; d<Dim ; d++ ) coords[d] *= s; }
		Real InnerProduct( const Point &p )	const
		{
			Real dot={};
#if 0
			if constexpr( std::derived_from< T , InnerProductSpace< Real , T > > )
#else // __cplusplus<202002l
			if constexpr( std::is_base_of_v< InnerProductSpace< Real , T > , T > )
#endif // __cplusplus
				for( int i=0 ; i<Dim ; i++ ) dot += T::Dot( p.coords[i] , coords[i] );
			else
				for( int i=0 ; i<Dim ; i++ ) dot += p.coords[i] * coords[i];
			return dot;
		}
		/////////////////////////////////

		T coords[Dim];

		template< typename ... Ts >
		Point( Ts ... ts )
		{
			static_assert( sizeof...(Ts)==0 || sizeof...(Ts)==Dim , "[ERROR] Invalid number of coefficients" );

			if constexpr( sizeof...(Ts)==0 ) for( unsigned int d=0 ; d<Dim ; d++ ) coords[d] = T{};
			else _set( ts... );
		}
		Point( T *c ){ for( unsigned int d=0 ; d<Dim ; d++ ) coords[d] = c[d]; }
		Point( const T *c ){ for( unsigned int d=0 ; d<Dim ; d++ ) coords[d] = c[d]; }

#if 0 // def NEW_CODE
		template< typename _T , typename _Real >
		explicit operator Point< _T , Dim , _Real >() const { Point< _T , Dim , _Real > p ; for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = static_cast< _T >( coords[d] ) ; return p; }
#else // !NEW_CODE
		template< typename _T , typename _Real >
		Point( const Point< _T , Dim , _Real > &p ){ for( unsigned int d=0 ; d<Dim ; d++ ) coords[d] = static_cast< T >( p[d] ); }

		template< typename _T , typename _Real >
		Point( Point< _T , Dim , _Real > &p ){ for( unsigned int d=0 ; d<Dim ; d++ ) coords[d] = static_cast< T >( p[d] ); }
#endif // NEW_CODE

		T& operator [] ( int idx ) { return coords[idx]; }
		const T& operator [] ( int idx ) const { return coords[idx]; }

		volatile T& operator [] ( int idx ) volatile { return coords[idx]; }
		const volatile T& operator [] ( int idx ) const volatile { return coords[idx]; }

		template< class ... Points > static Point CrossProduct( Points ... points )
		{
			static_assert( sizeof ... ( points )==Dim-1 , "Number of points in cross-product must be one less than the dimension" );
			SquareMatrix< Real , Dim > x;
			_AddColumnVector( x , 0 , points ... );
			Point p;
			for( int d=0 ; d<Dim ; d++ ) p[d] = ( d&1 ) ? -x.subDeterminant( Dim-1 , d ) : x.subDeterminant( Dim-1 , d );
			return p;
		}

		static Point CrossProduct( const Point *points )
		{
			XForm< Real , Dim > x;
			for( unsigned int d=0 ; d<Dim-1 ; d++ ) for( unsigned int c=0 ; c<Dim ; c++ ) x(d,c) = points[d][c];
			Point p;
			for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = ( d&1 ) ? -x.subDeterminant( Dim-1 , d ) : x.subDeterminant( Dim-1 , d );
			return p;
		}

		static Point CrossProduct( Point *points ){ return CrossProduct( ( const Point * )points ); }


		static Point< Real , Dim > Min( Point< Real , Dim > p , Point< Real , Dim > q ){ Point< Real , Dim > m ; for( int d=0 ; d<Dim ; d++ ) m[d] = std::min< Real >( p[d] , q[d] ) ; return m; }
		static Point< Real , Dim > Max( Point< Real , Dim > p , Point< Real , Dim > q ){ Point< Real , Dim > m ; for( int d=0 ; d<Dim ; d++ ) m[d] = std::max< Real >( p[d] , q[d] ) ; return m; }

		friend std::ostream &operator << ( std::ostream &os , const Point &p )
		{
			os << "( ";
			for( int d=0 ; d<Dim ; d++ )
			{
				if( d ) os << " , ";
				os << p[d];
			}
			return os << " )";
		}

	};

	template< typename T , typename Real >
	struct Point< T , (unsigned int)-1 , Real > : public InnerProductSpace< Real , Point< T , (unsigned int)-1 , Real > >
	{
	public:
		/////////////////////////////////
		// Inner product space methods //
		void Add( const Point& p )
		{
			if( !_dim ){ _resize( p._dim ) ; for( unsigned int i=0 ; i<_dim ; i++ ) _coords[i] = p._coords[i]; }
			else if( _dim==p._dim ) for( unsigned int i=0 ; i<_dim ; i++ ) _coords[i] += p._coords[i];
			else MK_ERROR_OUT( "Dimensions don't match: " , _dim , " != " , p._dim );
		}
		void Scale( Real s ){ for( unsigned int i=0 ; i<_dim ; i++ ) (*this)[i] *= s; }
		Real InnerProduct( const Point &p ) const
		{
			if( _dim!=p._dim ) MK_ERROR_OUT( "Dimensions differ: " , _dim , " != " , p._dim );

			Real dot={};
#if 0
			if constexpr( std::derived_from< T , InnerProductSpace< Real , T > > )
#else // __cplusplus<202002l
			if constexpr( std::is_base_of_v< InnerProductSpace< Real , T > , T > )
#endif // __cplusplus
				for( int i=0 ; i<_dim ; i++ ) dot += T::Dot( p._coords[i] , _coords[i] );
			else
				for( int i=0 ; i<_dim ; i++ ) dot += p._coords[i] * _coords[i];
			return dot;
		}
		/////////////////////////////////

		Point( void ) : _coords(nullptr) , _dim(0){}
		Point( size_t dim ) : _coords(nullptr) , _dim(0) { if( dim ){ _resize( (unsigned int)dim ) ; for( unsigned int d=0 ; d<dim ; d++ ) _coords[d] = T{}; } }
		Point( const Point &p ) : _coords(nullptr) , _dim(0) { if( p._dim ){ _resize( p._dim ) ; for( unsigned int d=0 ; d<_dim ; d++ ) _coords[d] = p._coords[d]; } }
		template< typename _T , typename _Real >
		Point( const Point< _T , (unsigned int)-1 , _Real > &p ) : Point( p.dim() ) { for( unsigned int d=0 ; d<_dim ; d++ ) _coords[d] = static_cast< T >( p[d] ); }
		template< typename _T , typename _Real >
		Point( Point< _T , (unsigned int)-1 , _Real > &p ) : Point( p.dim() ) { for( unsigned int d=0 ; d<_dim ; d++ ) _coords[d] = static_cast< T >( p[d] ); }

		~Point( void ){ delete[] _coords ; _coords = nullptr; }



		Point &operator = ( const Point &p )
		{
			if( _dim!=p._dim && _dim!=0 ) MK_ERROR_OUT( "Dimensions don't match: " , _dim , " != " , p._dim );

			if( !_dim ) _resize( p._dim );
			for( unsigned int d=0 ; d<_dim ; d++ ) _coords[d] = p._coords[d];

			return *this;
		}

		unsigned int dim( void ) const { return _dim; }
		T &operator[]( size_t idx ){ return _coords[idx]; }
		const T &operator[]( size_t idx ) const { return _coords[idx]; }

		friend std::ostream &operator << ( std::ostream &os , const Point &p )
		{
			os << "( ";
			for( size_t d=0 ; d<p.dim() ; d++ )
			{
				if( d ) os << " , ";
				os << p[d];
			}
			return os << " )";
		}
	protected:
		T *_coords;
		unsigned int _dim;
		void _resize( unsigned int dim )
		{
			delete[] _coords;
			_coords = new Real[dim];
			_dim = dim;
		}
	};

	template< typename T , unsigned int Dim , typename Real > Point< T , Dim , Real > operator * ( Real s , Point< T , Dim , Real > p ){ return p*s; }

	/** This templated class represents a Ray.*/
	template< typename T , unsigned int Dim , typename Real=T >
	class Ray
	{
	public:
		/** The starting point of the ray */
		Point< T , Dim , Real> position;

		/** The direction of the ray */
		Point< T , Dim , Real> direction;

		/** The default constructor */
		Ray( void ){}

		/** The constructor settign the the position and direction of the ray */
		Ray( const Point< T , Dim , Real> &position , const Point< T , Dim , Real> &direction ) : position(position) , direction(direction) {}

		/** This method computes the translation of the ray by p and returns the translated ray.*/
		Ray  operator +  ( const Point< T , Dim , Real> &p ) const { return Ray( position + p , direction );}

		/** This method translates the current ray by p.*/
		Ray &operator += ( const Point< T , Dim , Real> &p ){ position += p ; return *this; }

		/** This method computes the translation of the ray by -p and returns the translated ray.*/
		Ray  operator -  ( const Point< T , Dim , Real> &p ) const { return Ray( position - p , direction );}

		/** This method translates the current ray by -p.*/
		Ray &operator -= ( const Point< T , Dim , Real> &p ){ position -= p ; return *this; }

		/** This method returns the point at a distance of t along the ray. */
		Point< T , Dim , Real> operator() ( Real t ) const { return position + direction*t; }
	};


	template< class Real , int Cols , int Rows >
	class Matrix : public InnerProductSpace< Real , Matrix< Real , Cols , Rows > >
	{
	public:
		//////////////////////////
		// Vector space methods //
		void Add            ( const Matrix& m );
		void Scale          ( Real s );
		Real InnerProduct   ( const Matrix& m ) const;
		//////////////////////////

		Real coords[Cols][Rows];
		Matrix ( void ) { memset( coords , 0 , sizeof( Real ) * Cols * Rows ); }
		template<class Real2>
		operator Matrix< Real2 , Cols , Rows > ( void ) const
		{
			Matrix< Real2, Cols , Rows > m;
			for( int c=0 ; c<Cols ; c++ ) for ( int r=0 ; r<Rows ; r++ ) m.coords[c][r] = Real2( coords[c][r] ); 
			return m;
		}
		template<int C,int R>
		Matrix( const Matrix< Real , C , R> &m )
		{
			for(int i=0;i<Cols && i<C;i++)
				for(int j=0;j<Rows && j<R;j++)
					coords[i][j]=m.coords[i][j];
		}
		Real& operator () ( unsigned int c , unsigned int r ) { return coords[c][r]; }
		const Real& operator () ( unsigned int c , unsigned int r ) const { return coords[c][r]; }

		volatile Real& operator () ( unsigned int c , unsigned int r ) volatile { return coords[c][r]; }
		volatile const Real& operator () ( unsigned int c , unsigned int r ) volatile const { return coords[c][r]; }

		template<int Cols1>
		Matrix<Real,Cols1,Rows> operator * ( const Matrix< Real , Cols1 , Cols >& m ) const;

		Matrix< Real , Rows , Cols > transpose( void ) const;

		template< typename T >
		Point< T , Rows , Real > operator * ( const Point<  T , Cols , Real > &v ) const;

		template< typename T >
		Point< T , Rows , Real > operator() ( const Point<  T , Cols , Real > &v ) const;

		template< typename T >
		Real operator () ( const Point< T , Rows , Real >& v1 , const Point< T , Cols , Real >& v2 ) const;

		friend std::ostream &operator << ( std::ostream &os , const Matrix &m )
		{
			os << "{ ";
			for( int r=0 ; r<Rows ; r++ )
			{
				if( r ) os << " , ";
				os << "{ ";
				for( int c=0 ; c<Cols ; c++ )
				{
					if( c ) os << " , ";
					os << m.coords[c][r];
				}
				os << " }";
			}
			return os << " }";
		}
	};

	template< class Real , int Rows >
	class Matrix< Real , 0 , Rows > : public InnerProductSpace< Real , Matrix< Real , 0 , Rows > >
	{
		static const unsigned int Cols = 0;
	public:
		//////////////////////////
		// Vector space methods //
		void Add            ( const Matrix& m ){}
		void Scale          ( Real s ){}
		Real InnerProduct   ( const Matrix& m ) const { return (Real)0; }
		//////////////////////////

		Matrix( void ){}

		template< class Real2 >
		operator Matrix< Real2 , Cols , Rows > ( void ) const{}

		template< int C , int R >
		Matrix( const Matrix< Real , C , R > &m ){}

		Real& operator () ( unsigned int c , unsigned int r ) { MK_ERROR_OUT( "Should not be accessing the entries of this matrix" ) ; Real v=0 ; return v; }
		const Real& operator () ( unsigned int c , unsigned int r ) const { MK_ERROR_OUT( "Should not be accessing the entries of this matrix" ) ; Real v=0 ; return v; }

		volatile Real& operator () ( unsigned int c , unsigned int r ) volatile { MK_ERROR_OUT( "Should not be accessing the entries of this matrix" ) ; Real v=0 ; return v; }
		volatile const Real& operator () ( unsigned int c , unsigned int r ) volatile const { MK_ERROR_OUT( "Should not be accessing the entries of this matrix" ) ; Real v=0 ; return v; }

		template< int Cols1 >
		Matrix< Real , Cols1 , Rows > operator * ( const Matrix< Real , Cols1 , Cols >& m ) const { return Matrix< Real , Cols1 , Rows >(); }

		Matrix< Real , Rows , Cols > transpose( void ) const{ return Matrix< Real , Rows , Cols >(); }

		template< class Real2 >
		Point< Real2 , Rows > operator * ( const Point< Real2 , Cols >& v ) const { return Point< Real2 , Rows >(); }

		template< class Real2 >
		Point< Real2 , Rows > operator () ( const Point< Real2 , Cols >& v ) const { return Point< Real2 , Rows >(); }

		friend std::ostream &operator << ( std::ostream &os , const Matrix &m ){ return os <<  "{ }"; }
	};

	template< class Real , int Cols >
	class Matrix< Real , Cols , 0 > : public InnerProductSpace< Real , Matrix< Real , Cols , 0 > >
	{
		static const unsigned int Rows = 0;
	public:
		//////////////////////////
		// Vector space methods //
		void Add            ( const Matrix& m ){}
		void Scale          ( Real s ){}
		Real InnerProduct   ( const Matrix& m ) const { return (Real)0; }
		//////////////////////////

		Matrix( void ){}

		template< class Real2 >
		operator Matrix< Real2 , Cols , Rows > ( void ) const{}

		template< int C , int R >
		Matrix( const Matrix< Real , C , R > &m ){}

		Real& operator () ( unsigned int c , unsigned int r ) { MK_ERROR_OUT( "Should not be accessing the entries of this matrix" ) ; return (Real&)*(Real*)NULL; }
		const Real& operator () ( unsigned int c , unsigned int r ) const { MK_ERROR_OUT( "Should not be accessing the entries of this matrix" ) ; return (Real&)*(Real*)NULL; }

		volatile Real& operator () ( unsigned int c , unsigned int r ) volatile { MK_ERROR_OUT( "Should not be accessing the entries of this matrix" ) ; return (Real&)*(Real*)NULL; }
		volatile const Real& operator () ( unsigned int c , unsigned int r ) volatile const { MK_ERROR_OUT( "Should not be accessing the entries of this matrix" ) ; return (Real&)*(Real*)NULL; }

		template< int Cols1 >
		Matrix< Real , Cols1 , Rows > operator * ( const Matrix< Real , Cols1 , Cols >& m ) const { return Matrix< Real , Cols1 , Rows >(); }

		Matrix< Real , Rows , Cols > transpose( void ) const{ return Matrix< Real , Rows , Cols >(); }

		template< class Real2 >
		Point< Real2 , Rows > operator * ( const Point< Real2 , Cols >& v ) const { return Point< Real2 , Rows >(); }

		template< class Real2 >
		Point< Real2 , Rows > operator () ( const Point< Real2 , Cols >& v ) const { return Point< Real2 , Rows >(); }

		friend std::ostream &operator << ( std::ostream &os , const Matrix &m ){ return os <<  "{ }"; }
	};

	// Need forward declaration to support the characteristic polynomial
	//	namespace Polynomial{ template< unsigned int Dim , unsigned int Degree , typename T , typename Real=T > class Polynomial; }
	namespace Polynomial{ template< unsigned int Dim , unsigned int Degree , typename T , typename Real > class Polynomial; }

	template< class Real , int Dim >
	class Matrix< Real , Dim , Dim > : public Algebra< Real , Matrix< Real , Dim , Dim > > , public InnerProductSpace< Real , Matrix< Real , Dim , Dim > >
	{
	public:
		static const int Cols = Dim;
		static const int Rows = Dim;
		////////////////////////////////
		// Vector space methods       //
		void Add            ( const Matrix& m );
		void Scale          ( Real s );
		Real InnerProduct   ( const Matrix& m ) const;
		////////////////////////////////
		// Additional algebra methods //
		void Multiply       ( const Matrix& m );
		void SetIdentity    ( void );
		////////////////////////////////

		Real coords[Dim][Dim];
		Matrix ( void ) { memset( coords , 0 , sizeof( Real ) * Cols * Rows ); }
		template< class Real2 >
		operator Matrix< Real2 , Cols , Rows > ( void ) const
		{
			Matrix< Real2, Cols , Rows > m;
			for( int c=0 ; c<Cols ; c++ ) for ( int r=0 ; r<Rows ; r++ ) m.coords[c][r] = Real2( coords[c][r] ); 
			return m;
		}
		template<int C,int R>
		Matrix( const Matrix< Real , C , R> &m )
		{
			for(int i=0;i<Cols && i<C;i++)
				for(int j=0;j<Rows && j<R;j++)
					coords[i][j]=m.coords[i][j];
		}
		Real& operator () ( unsigned int c , unsigned int r ){ return coords[c][r]; }
		const Real& operator () ( unsigned int c , unsigned int r ) const { return coords[c][r]; }

		volatile Real& operator () ( unsigned int c , unsigned int r ) volatile { return coords[c][r]; }
		volatile const Real& operator () ( unsigned int c , unsigned int r ) volatile const { return coords[c][r]; }

		template<int Cols1>
		Matrix< Real , Cols1 , Dim > operator * ( const Matrix< Real , Cols1 , Dim >& m ) const;

		Matrix< Real , Dim , Dim > transpose( void ) const;

		template< typename T >
		Point< T , Dim , Real > operator * ( const Point< T , Dim , Real > &v ) const;

		template< typename T >
		Point< T , Dim , Real > operator() ( const Point< T , Dim , Real > &v ) const;

		template< typename T >
		Real operator () ( const Point< T , Dim , Real >& v1 , const Point< T , Dim , Real >& v2 ) const;

		friend std::ostream &operator << ( std::ostream &os , const Matrix &m )
		{
			os << "{ ";
			for( int r=0 ; r<Rows ; r++ )
			{
				if( r ) os << " , ";
				os << "{ ";
				for( int c=0 ; c<Cols ; c++ )
				{
					if( c ) os << " , ";
					os << m.coords[c][r];
				}
				os << " }";
			}
			return os << " }";
		}

		static Matrix Identity( void )
		{
			Matrix M;
			M.SetIdentity();
			return M;
		}
		Real subDeterminant( int c , int r ) const;
		Real determinant( void ) const;
		Real trace( void ) const;
		Matrix inverse( bool& success ) const;
		Matrix inverse( void ) const;
		class Polynomial::Polynomial< 1 , Dim , Real , Real > characteristicPolynomial( void ) const;

		template< class Real2 > Point< Real2 , Dim-1 > operator () ( const Point< Real2 , Dim-1 >& v ) const;
	protected:
		friend Matrix< double , Dim+1 , Dim+1 >;
		class Polynomial::Polynomial< 1 , Dim , Real , Real > _characteristicPolynomial( Matrix< char , Dim , Dim > mask ) const;
	};

#if 0
	template< class Real , int Dim1 , int Dim2 > Matrix< Real , Dim2 , Dim1 > operator * ( const SquareMatrix< Real , Dim1 >& m1 , const Matrix< Real , Dim2 , Dim1 >& m2 ){ return ( Matrix< Real , Dim1 , Dim1 > )m1 * m2; }
	template< class Real , int Dim1 , int Dim2 > Matrix< Real , Dim1 , Dim2 > operator * ( const Matrix< Real , Dim1 , Dim2 >& m1 , const SquareMatrix< Real , Dim1 >& m2 ){ return m1 * ( Matrix< Real , Dim1 , Dim1 > )m2; }
#endif

	template< class Real >
	class Matrix< Real , 0 , 0 > : public Algebra< Real , Matrix< Real , 0 , 0 > > , public InnerProductSpace< Real , Matrix< Real , 0 , 0 > >
	{
	public:
		static const unsigned int Dim = 0;
		static const int Cols = 0;
		static const int Rows = 0;

		////////////////////////////////
		// Vector space methods       //
		void Add            ( const Matrix& m ){}
		void Scale          ( Real s ){}
		Real InnerProduct   ( const Matrix& m ) const { return (Real)0; }
		////////////////////////////////
		// Additional algebra methods //
		void Multiply ( const Matrix& m ){;}
		void SetIdentity( void ){;}
		////////////////////////////////

		Matrix( void ){}

		template< class Real2 >
		operator Matrix< Real2 , Cols , Rows > ( void ) const{}

		template< int C , int R >
		Matrix( const Matrix< Real , C , R > &m ){}

		Real& operator () ( unsigned int c , unsigned int r ){ MK_ERROR_OUT( "Should not be accessing the entries of this matrix" ) ; Real v=0 ; return v; }
		const Real& operator () ( unsigned int c , unsigned int r ) const { MK_ERROR_OUT( "Should not be accessing the entries of this matrix" ) ; Real v=0 ; return v; }

		volatile Real& operator () ( unsigned int c , unsigned int r ) volatile { MK_ERROR_OUT( "Should not be accessing the entries of this matrix" ) ; Real v=0 ; return v; }
		volatile const Real& operator () ( unsigned int c , unsigned int r ) volatile const { MK_ERROR_OUT( "Should not be accessing the entries of this matrix" ) ; Real v=0 ; return v; }

		template< int Cols1 >
		Matrix< Real , Cols1 , Rows > operator * ( const Matrix< Real , Cols1 , Cols >& m ) const { return Matrix< Real , Cols1 , Rows >(); }

		template< class Real2 >
		Point< Real2 , Rows > operator * ( const Point< Real2 , Cols >& v ) const { return Point< Real2 , Rows >(); }

		template< class Real2 >
		Point< Real2 , Rows > operator () ( const Point< Real2 , Cols >& v ) const { return Point< Real2 , Rows >(); }

		friend std::ostream &operator << ( std::ostream &os , const Matrix &m ){ return os <<  "{ }"; }

		Real determinant( void ) const { return 0; }
		Real trace( void ) const { return 0; }
		Matrix transpose( void ) const { return Matrix(); }
	};

	template< typename T , unsigned int Dim1 , unsigned int Dim2 , typename Real >
	Matrix< Real , Dim2 , Dim1 > OuterProduct( Point< T , Dim1 , Real > p1 , Point< T , Dim2 , Real > p2 )
	{
		Matrix< Real , Dim2 , Dim1 > op;
#if 0
		if constexpr( std::derived_from< T , InnerProductSpace< Real , T > > )
#else // __cplusplus<202002l
		if constexpr( std::is_base_of_v< InnerProductSpace< Real , T > , T > )
#endif // __cplusplus
			for( unsigned int i=0 ; i<Dim1 ; i++ ) for( unsigned int j=0 ; j<Dim2 ; j++ ) op(j,i) = T::Dot( p1[i] , p2[j] );
		else
			for( unsigned int i=0 ; i<Dim1 ; i++ ) for( unsigned int j=0 ; j<Dim2 ; j++ ) op(j,i) = p1[i] * p2[j];
		return op;
	}
	template< typename T , unsigned int Dim , typename Real >
	SquareMatrix< Real , Dim > OuterProduct( Point< T , Dim , Real > p1 , Point< T , Dim , Real > p2 )
	{
		SquareMatrix< Real , Dim > op;
#if 0
		if constexpr( std::derived_from< T , InnerProductSpace< Real , T > > )
#else // __cplusplus<202002l
		if constexpr( std::is_base_of_v< InnerProductSpace< Real , T > , T > )
#endif // __cplusplus
			for( unsigned int i=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j<Dim ; j++ ) op(j,i) = T::Dot( p1[i] , p2[j] );
		else
			for( unsigned int i=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j<Dim ; j++ ) op(j,i) = p1[i] * p2[j];
		return op;
	}

	template< class V , int Dim , class _R = typename V::R >
	class Gradient : public VectorSpace< _R , Gradient< V , Dim , _R > >
	{
	public:
		//////////////////////////
		// Vector space methods //
		void Add            ( const Gradient& g ) { for( int c=0  ; c<Dim ; c++ ) grad[c] += g.grad[c]; }
		void Scale          ( _R s ) { for( int c=0 ; c<Dim ; c++ ) grad[c] *= s; }
		//                      //
		//////////////////////////

		V grad[Dim];
#if 1
		Gradient( void ) { for( int d=0 ; d<Dim ;  d++ ) grad[d] = V{}; }
#else
		Gradient( void ) { for( int d=0 ; d<Dim ;  d++ ) grad[d] *= 0; }
#endif
		V& operator[] ( int idx ) { return grad[idx]; }
		const V& operator[] ( int idx ) const { return grad[idx]; }

		template< class V2 , class _R2>
		operator Gradient< V2, Dim , _R2 > ( void ) const
		{
			Gradient< V2 , Dim , _R2 > g;
			for( int d=0 ; d<Dim ; d++ ) g.grad[d] = V2( grad[d] ); 
			return g;
		}

		template< class Real >
		Gradient Project( const Point< Real , Dim >& dir ) const
		{
			V dot;
			Gradient g;
#if 1
			g = V{};
			dot = V{};
#else
			g *= 0;
			dot *= 0;
#endif
			Real len = Real( sqrt( Point< Real , Dim >::SquareNorm( dir ) ) );
			if( !len ) return g;
			Point< Real , Dim > _dir = dir / len;
			for( int d=0 ; d<Dim ; d++ ) dot += grad[d] * _dir[d];
			for( int d=0 ; d<Dim ; d++ ) g.grad[d] = dot * _dir[d];
			return g;
		}
	};

	template< class V , int Dim , class _R = typename V::R >
	class ConstantFunction : public VectorSpace< _R , ConstantFunction< V , Dim , _R > >
	{
	public:
		V value;
		Gradient< V , Dim , _R > grad;
		ConstantFunction( void ) { value *= 0 , grad *= 0;}

		template< class Real > V operator( ) ( const Point< Real , Dim >& p ) const { return value; }
		template< class Real > Gradient< V , Dim , _R > gradient( const Point< Real , Dim >& p ) const { return grad; }

		//////////////////////////
		// Vector space methods //
		void Add            ( const ConstantFunction& cf ) { value += cf.value; }
		void Scale          ( _R s ) { value *= s , this->offset *= s; }
		//////////////////////////
	};

	//template< class V , int Dim , class _R = typename V::R >
	template< class V , int Dim , class _R = typename FieldOf< V >::F >
	class LinearFunction : public VectorSpace< _R , LinearFunction< V , Dim , _R > >
	{
	public:
#if 1
		std::conditional_t< std::is_same_v< V , _R > , Point< _R , Dim > , Gradient< V , Dim , _R > > grad;
#else
		Gradient< V , Dim , _R > grad;
#endif
		V offset;
		LinearFunction( void ) : offset(V{}) {}
		template< class Real >
		V operator( ) ( const Point< Real , Dim >& p ) const
		{
			V v{};
			for( int d=0 ; d<Dim ; d++ ) v += grad[d] * p[d];
			v -= offset;
			return v;
		}
		template< class Real >
		LinearFunction fitToHyperplane( const Point< Real , Dim >& p , const Point< Real , Dim >& n ) const
		{
			LinearFunction f;
			Real len = Point< Real , Dim >::SquareNorm( n );
			if( !len )
			{
				f.grad *= 0;
				f.offset = -(*this)( p );
			}
			else
			{
				Point< Real , Dim > normal = n / Real( sqrt( double( len ) ) );
				V dot;
				dot *= 0;
				for( int d=0 ; d<Dim ; d++ ) dot += grad[d] * normal[d];
				for( int d=0 ; d<Dim ; d++ ) f.grad[d] = grad[d] - dot * normal[d];
				f.offset *= 0;
				f.offset = -(*this)( p ) + f( p );
			}
			return f;
		}
		template< class V2 , class _R2 >
		operator LinearFunction< V2 , Dim , _R2 > ( void ) const
		{
			LinearFunction< V2 , Dim , _R2 > lf;
			lf.offset = V2 ( offset );
			lf.grad = Gradient< V2 , Dim , _R2 >( grad );
			return lf;
		}
		template< class Real >
		Gradient< V , Dim , _R > gradient( const Point< Real , Dim >& p ) const { return grad; }

		// Warning, this function requires the taking of an inverse, which may fail...
		template< class Real >
		static LinearFunction BestFit( ConstPointer( Point< Real , Dim > ) points , ConstPointer( V ) values , int count )
		{
			LinearFunction lf;
			V constraint[Dim];
			SquareMatrix< Real , Dim > M , Minv;
#if 1
			for( int d=0 ; d<Dim ; d++ ) constraint[d] = V{};
#else
			M *= 0;
			for( int d=0 ; d<Dim ; d++ ) constraint[d] *= 0;
#endif
			for( int i=0 ; i<count ; i++ )
			{
				for( int k=0 ; k<Dim ; k++ ) for( int l=0 ; l<Dim ; l++ ) M( k , l ) += points[i][k] * points[i][l];
				for( int j=0 ; j<count ; j++ ) for( int k=0 ; k<Dim ; k++ ) for( int l=0 ; l<Dim ; l++ ) M( k , l ) -= points[i][k] * points[j][l] / Real( count ); 

				for( int d=0 ; d<Dim ; d++ ) constraint[d] += values[i] * points[i][d];
				for( int j=0 ; j<count ; j++ ) for( int d=0 ; d<Dim ; d++ ) constraint[d] -= values[j] * points[i][d] / Real( count );
			}
			bool success;
			Minv = M.inverse( success );
			if( !success ){ MK_THROW( "Could not invert matrix" ); }

#if 1
			// [WARNING] Apparently "offset" could have been initialized to something awful so that multiplication by zero is still awful
			lf.grad = std::conditional_t< std::is_same_v< V , _R > , Point< _R , Dim > , Gradient< V , Dim , _R > >{};
			lf.offset = V{};
#else
			lf *= 0;
#endif
			for( int c=0 ; c<Dim ; c++ ) for( int r=0 ; r<Dim ; r++ ) lf.grad[r] += constraint[c] * Minv( c , r );
			for( int i=0 ; i<count ; i++ )
			{
				for( int d=0 ; d<Dim ; d++ ) lf.offset += lf.grad[d] * points[i][d];
				lf.offset -= values[i];
			}
			lf.offset /= Real( count );
			return lf;
		}

		friend std::ostream &operator << ( std::ostream &stream , const LinearFunction &lf )
		{
			stream << " < " << lf.grad << " , . >";
			if     ( lf.offset<0 ) return stream << " + " << - lf.offset;
			else if( lf.offset>0 ) return stream << " - " <<   lf.offset;
			else                   return stream;
		};


		//////////////////////////
		// Vector space methods //
		void Add            ( const LinearFunction& lf ) { grad += lf.grad , offset += lf.offset; }
		void Scale          ( _R s ) { grad *= s , offset *= s; }
		//////////////////////////
	};

	template< class Real , int Dim >
	struct OrientedPoint
	{
		Point< Real , Dim > position , normal;
		template< class Real2 > operator Point< Real2, Dim > ( void ) const { return Point< Real2 , Dim >( position ); }
	};

	template< typename Real , typename Data >
	struct ProjectiveData
	{
		Data data;
		Real weight;
		ProjectiveData( Data d=Data() , Real w=(Real)0 ) : data(d) , weight(w) { ; }
		operator Data (){ return weight!=0 ? data/weight : data*weight; }
		Data value( void ) const { return weight!=0 ? data/weight : data*weight; }
		ProjectiveData& operator += ( const ProjectiveData& p ){ data += p.data , weight += p.weight ; return *this; }
		ProjectiveData& operator -= ( const ProjectiveData& p ){ data -= p.data , weight += p.weight ; return *this; }
		ProjectiveData& operator *= ( Real s ){ data *= s , weight *= std::abs(s) ; return *this; }
		ProjectiveData& operator /= ( Real s ){ data /= s , weight /= std::abs(s) ; return *this; }
		ProjectiveData  operator +  ( const ProjectiveData& p ) const { return ProjectiveData( data+p.data , weight+p.weight ); }
		ProjectiveData  operator -  ( const ProjectiveData& p ) const { return ProjectiveData( data-p.data , weight+p.weight ); }
		ProjectiveData  operator *  ( Real s ) const { return ProjectiveData( data*s , weight*std::abs(s) ); }
		ProjectiveData  operator /  ( Real s ) const { return ProjectiveData( data/s , weight/std::abs(s) ); }
		ProjectiveData  operator -  ( void ) const{ return ProjectiveData( -data , weight ); }
	};

	template< class Real > using Point2D = Point< Real , 2 >;
	template< class Real > using Point3D = Point< Real , 3 >;

	template<class Real>
	class OrientedPoint2D : public OrientedPoint<Real,2>{;};
	template<class Real>
	class OrientedPoint3D : public OrientedPoint<Real,3>{;};

	template< typename Real , unsigned int Dim >
	class XForm : public SquareMatrix< Real , Dim >
	{
		using SquareMatrix< Real , Dim >::coords;
	public:
		using SquareMatrix< Real , Dim >::Identity;
		using SquareMatrix< Real , Dim >::operator *;

		XForm( void ) : SquareMatrix< Real , Dim >(){;}
		XForm( const SquareMatrix< Real , Dim > &xForm ){ memcpy( coords , xForm.coords , sizeof(Real)*Dim*Dim ); };
		XForm( const SquareMatrix< Real , Dim-1 > &L , Point< Real , Dim-1 > t=Point< Real , Dim >() ) : XForm()
		{
			for( unsigned int i=0 ; i<Dim-1 ; i++ ) for( unsigned int j=0 ; j<Dim-1 ; j++ ) coords[i][j] = L.coords[i][j];
			for( unsigned int j=0 ; j<Dim-1 ; j++ ) coords[Dim-1][j] = t[j];
			coords[Dim-1][Dim-1] = (Real)1.;
		}
		Point< Real , Dim > operator * ( Point< Real , Dim > p ) const { return SquareMatrix< Real , Dim >::operator * ( p ); }
		Point< Real , Dim-1 > operator * ( Point< Real , Dim-1 > p ) const
		{
			Point< Real , Dim > _p;
			for( unsigned int d=0 ; d<Dim-1 ; d++ ) _p[d] = p[d];
			_p[Dim-1] = 1;
			_p = SquareMatrix< Real , Dim >::operator * ( _p );
			_p /= _p[Dim-1];
			for( int d=0 ; d<Dim-1 ; d++ ) p[d] = _p[d];
			return p;
		}

		static XForm Scale( Real s )
		{
			XForm xForm = Identity();
			for( int d=0 ; d<Dim-1 ; d++ ) xForm(d,d) = s;
			return xForm;
		}

		static XForm Scale( Point< Real , Dim-1 > s )
		{
			XForm xForm = Identity();
			for( int d=0 ; d<Dim-1 ; d++ ) xForm(d,d) = s[d];
			return xForm;
		}

		static XForm Translate( Point< Real , Dim-1 > t )
		{
			XForm xForm = Identity();
			for( int d=0 ; d<Dim-1 ; d++ ) xForm(Dim-1,d) = t[d];
			return xForm;
		}
	};
	template< typename Real > using XForm4x4 = XForm< Real , 4 >;
	template< typename Real > using XForm3x3 = XForm< Real , 3 >;
	template< typename Real > using XForm2x2 = XForm< Real , 2 >;

	///////////////
	// Simplices //
	///////////////
	template< unsigned int K , typename Index=unsigned int > struct SimplexIndex;

	template< unsigned int K > struct Factorial{ static const unsigned long long Value = Factorial< K-1 >::Value * K; };
	template<> struct Factorial< 0 >{ static const unsigned long long Value = 1; };

	template< class Real , unsigned int Dim , unsigned int K >
	struct Simplex
	{	
		Point< Real , Dim > p[K+1];
		Simplex( void ){ static_assert( K<=Dim , "[ERROR] Bad simplex dimension" ); }

		Simplex( const Point< Real , Dim > p[K+1] ) : Simplex() { _init( p ); }

		template< typename ... Points >
		Simplex( Point< Real , Dim > p , Points ... ps ) : Simplex()
		{
			static_assert( sizeof...(Points)==K , "[ERROR] Wrong number of points" );
			const Point< Real , Dim > _p[] = { p , ps... };
			_init( _p );
		}

		Point< Real , Dim >& operator[]( unsigned int k ){ return p[k]; }
		const Point< Real , Dim >& operator[]( unsigned int k ) const { return p[k]; }
		Real measure( void ) const { return (Real)sqrt( squareMeasure() ); }

		template< unsigned int _K=K >
		std::enable_if_t< _K==Dim , Real > volume( bool signedVolume=false ) const
		{
			SquareMatrix< double , K > M;
			for( unsigned int k=0 ; k<K ; k++ )
			{
				Point< double , Dim > d = p[k+1] - p[0];
				for( unsigned int j=0 ; j<K ; j++ ) M(k,j) = d[j];
			}
			return ( signedVolume ? -M.determinant() : fabs( -M.determinant() ) ) / Factorial< K >::Value;
		}

		Real squareMeasure( void ) const { return metric().determinant() / ( Factorial< K >::Value * Factorial< K >::Value ); }
		SquareMatrix< Real , K > metric( void ) const
		{
			SquareMatrix< Real , K > m;
			for( unsigned int i=1 ; i<=K ; i++ ) for( unsigned int j=1 ; j<=K ; j++ ) m(i-1,j-1) = Point< Real , Dim >::Dot( p[i]-p[0] , p[j]-p[0] );
			return m;
		}

		Point< Real , Dim > center( void ) const
		{
			Point< Real , Dim > c;
			for( unsigned int k=0 ; k<=K ; k++ ) c += p[k];
			return c / (K+1);
		}

		Point< Real , Dim > randomSample( void ) const
		{
			while( true )
			{
				Point< Real , K > bc;
				Real sum = 0;
				for( unsigned int d=0 ; d<K ; d++ ) sum += ( bc[d] = Random< Real >() );
				if( sum<=1 )
				{
					Point< Real , Dim > q = p[0];
					for( unsigned int d=0 ; d<K ; d++ ) q += ( p[d+1] - p[0] ) * bc[d];
					return q;
				}
			}
		}

		void split( const Real values[K+1] , std::vector< Simplex >& back , std::vector< Simplex >& front ) const;
		void split( Point< Real , Dim > pNormal , Real pOffset , std::vector< Simplex >& back , std::vector< Simplex >& front ) const;
		Point< Real , Dim > operator()( const Real weights[K+1] ) const
		{
			Point< Real , Dim > q;
			for( unsigned int k=0 ; k<=K ; k++ ) q += p[k] * weights[k];
			return q;
		}

		Point< Real , Dim > operator()( Point< Real , K+1 > bc ) const
		{
			Point< Real , Dim > q;
			for( unsigned int k=0 ; k<=K ; k++ ) q += p[k] * bc[k];
			return q;
		}

		Point< Real , Dim > operator()( Point< Real , K > x ) const
		{
			Point< Real , Dim > q = p[0];
			for( unsigned int k=0 ; k<K ; k++ ) q += ( p[k+1] - p[0] ) * x[k];
			return q;
		}


		template< unsigned int _K=K >
		typename std::enable_if< _K==Dim-1 , Point< Real , Dim > >::type normal( void ) const
		{
			Point< Real , Dim > d[Dim-1];
			for( int k=1 ; k<Dim ; k++ ) d[k-1] = p[k] - p[0];
			return Point< Real , Dim >::CrossProduct( d );
		}

		template< unsigned int _K=K >
		typename std::enable_if< _K==Dim , bool >::type isInterior( Point< Real , Dim > p ) const
		{
			struct SimplexIndex< Dim > si;
			for( unsigned int d=0 ; d<=Dim ; d++ ) si[d] = d;
			unsigned int count = 0;
			for( unsigned int d=0 ; d<=Dim ; d++ )
			{
				bool oriented = d%2;
				SimplexIndex< Dim-1 > fi = si.face( d );
				Simplex< Real , Dim , Dim-1 > f;
				for( unsigned int d=0 ; d<Dim ; d++ ) f[d] = this->p[ fi[d] ];
				Point< Real , Dim > c = f.center();
				Point< Real , Dim > n = oriented ? f.normal() : -f.normal();
				if( Point< Real , Dim >::Dot( p-c , n )<0 ) count++;
			}
			return count==0 || count==Dim+1;
		}

		Point< Real , K+1 > barycentricCoordinates( Point< Real , Dim > q ) const
		{
			SquareMatrix< Real , K > M;
			Point< Real , K > dot;
			for( unsigned int i=0 ; i<K ; i++ )
			{
				dot[i] = Point< Real , Dim >::Dot( q-p[0] , p[i+1]-p[0] );
				for( unsigned int j=0 ; j<K ; j++ ) M(i,j) = Point< Real , Dim >::Dot( p[i+1]-p[0] , p[j+1]-p[0] );
			}
			Point< Real , K > _b = M.inverse() * dot;
			Point< Real , K+1 > b;
			b[0] = 1.;
			for( unsigned int k=0 ; k<K ; k++ ) b[0] -= _b[k] , b[k+1] = _b[k];
			return b;
		}

		std::pair< Real , Point< Real , K+1 > > barycentricCoordinates( Ray< Real , Dim > r ) const
		{
			// Solve for (t,a_1,..,a_K) minimizing:
			//	E(t,a_1,...,a_K) = || r(t) - (1-a_1-...-a_k)p[0] - a_1*p[1] - ... a_K*p[K] ||^2
			//                   = || r.position-p[0] + r.direction * t + a_1*(p[0]-p[1]) + ... + a_K*(p[0]-p[K]) ||^2
			//                   = || v + A * a ||^2
			//                   = || v ||^2 + 2 * a^t * A^t * v + a^t * A^t * A * a
			Matrix< Real , K+1 , Dim > A;
			for( unsigned int k=0 ; k<=K ; k++ ) for( unsigned int d=0 ; d<Dim ; d++ )
				if( k==0 ) A(k,d) = r.direction[d];
				else       A(k,d) = p[0][d] - p[k][d];
			Point< Real , K+1 > a = - ( A.transpose() * A ).inverse() * ( A.transpose() * ( r.position - p[0] ) );
			Real sum = 0;
			Real t = a[0];
			for( unsigned int k=1 ; k<=K ; k++ ) sum += a[k];
			a[0] = (Real)1 - sum;
			return std::pair< Real , Point< Real , K+1 > >( t , a );
		}

		Point< Real , K+1 > nearestBC( Point< Real , Dim > p ) const
		{
			Point< Real , K+1 > bc = barycentricCoordinates( p );
			if constexpr( K==0 ) return bc;
#ifdef NEW_GEOMETRY_CODE
			else if constexpr( K==1 )
			{
				if( bc[0]<0 ) return Point< Real , 2 >( 0 , 1 );
				else if( bc[1]<0 ) return Point< Real , 2 >( 1 , 0 );
				else return bc;
			}
#endif // NEW_GEOMETRY_CODE
			else
			{
//#ifdef NEW_GEOMETRY_CODE
				// I think this is right, but it needs double-checking
#if 0
				for( unsigned int k=0 ; k<=K ; k++ ) if( bc[k]<0 )
				{
					Simplex< double , Dim , K-1 > _s;
					for( unsigned int i=0 ; i<K ; i++ ) _s[i] = this->p[(k+1+i)%(K+1)];
					Point< Real , K > _bc = _s.nearestBC( p );
					bc[k] = 0;
					for( unsigned int i=0 ; i<K ; i++ ) bc[(k+1+i)%(K+1)] = _bc[k];
					return bc;
				}
				return bc;
#else // !NEW_GEOMETRY_CODE
				unsigned int count = 0;
				for( unsigned int k=0 ; k<=K ; k++ ) if( bc[k]<0 ) count++;
				if( !count ) return bc;
				else
				{
					Real dist2 = std::numeric_limits< Real >::infinity();
					Point< Real , K+1 > bc;
					for( unsigned int k=0 ; k<=K ; k++ )
					{
						Simplex< Real , Dim , K-1 > _simplex;
						for( unsigned int _k=0 , idx=0 ; _k<=K ; _k++ ) if( k!=_k ) _simplex[idx++] = operator[](_k);
						Point< Real , K > _bc = _simplex.nearestBC( p );
						double _dist2 = Point< double , Dim >::SquareNorm( _simplex(_bc) - p );
						if( _dist2<dist2 )
						{
							dist2 = _dist2;
							for( unsigned int _k=0 , idx=0 ; _k<=K ; _k++ )
								if( k!=_k ) bc[_k] = _bc[idx++];
								else        bc[_k] = 0;
						}
					}
					return bc;
				}
#endif // NEW_GEOMETRY_CODE
			}
		}
		Point< Real , Dim > nearest( Point< Real , Dim > p ) const { return operator()( nearestBC(p) ); }


		friend std::ostream &operator << ( std::ostream &os , const Simplex &s )
		{
			os << "{ ";
			for( int k=0 ; k<=K ; k++ )
			{
				if( k ) os << " , ";
				os << s[k];
			}
			return os << " }";
		}

	protected:
		void _init( const Point< Real , Dim > p[K+1] ){ for( unsigned int k=0 ; k<=K ; k++ ) this->p[k] = p[k]; }

	};

	template< class Real , unsigned int Dim >	
	struct Simplex< Real , Dim , 0 >
	{
		Point< Real , Dim > p[1];
		Point< Real , Dim >& operator[]( unsigned int k ){ return p[k]; }
		const Point< Real , Dim >& operator[]( unsigned int k ) const { return p[k]; }
		Real squareMeasure( void ) const { return (Real)1.; }
		Real measure( void ) const { return (Real)1.; }
		Point< Real , Dim > center( void ) const { return p[0]; }
		Point< Real , Dim > operator()( const Real weights[1] ) const { return p[0] * weights[0]; }
		Point< Real , Dim > randomSample( void ) const { return p[0]; }
		void split( const Real values[1] , std::vector< Simplex >& back , std::vector< Simplex >& front ) const
		{
			if( values[0] ) back.push_back( *this );
			else            front.push_back( *this );
		}
		void split( Point< Real , Dim > pNormal , Real pOffset , std::vector< Simplex >& back , std::vector< Simplex >& front ) const
		{
			const Real values[] = { Point< Real , Dim >::Dot( p[0] , pNormal ) - pOffset };
			return split( values , back , front );
		}

		Point< Real , 1 > nearestBC( Point< Real , Dim > p ) const { return Point< Real , 1 >( (Real)1 ); }
		Point< Real , Dim > nearest( Point< Real , Dim > p ) const { return operator()( nearestBC(p) ); }
		Point< Real , 1 > barycentricCoordinates( Point< Real , Dim > q ) const { return Point< Real , 1 >( 1 ); };
		Point< Real , Dim > operator()( Point< Real , 1 > bc ) const { return p[0] * bc[0]; }

		template< unsigned int _K=0 >
		std::enable_if_t< _K==Dim , double > volume( bool ) const { return 1.; }

		friend std::ostream &operator << ( std::ostream &os , const Simplex &s )
		{
			return os << "{ " << s[0] << " }";
		}
	};

	template< class Real , unsigned int Dim , unsigned int K >
	Simplex< Real , Dim , K > operator * ( XForm< Real , Dim+1 > xForm , Simplex< Real , Dim , K > simplex )
	{
		for( unsigned int k=0 ; k<=K ; k++ ) simplex[k] = xForm * simplex[k];
		return simplex;
	}

	template< typename Index >
	struct EdgeTable
	{
	protected:
		struct _EdgeKey
		{
			Index key1 , key2;
			_EdgeKey( Index k1=0 , Index k2=0 ) : key1(k1) , key2(k2) {}
			bool operator == ( const _EdgeKey &key ) const  { return key1==key.key1 && key2==key.key2; }
			struct Hasher{ size_t operator()( const _EdgeKey &key ) const { return (size_t)( key.key1 * key.key2 ); } };
		};
		std::unordered_map< _EdgeKey , Index , typename _EdgeKey::Hasher > _edgeTable;
	public:
		template< typename InitializationFunction >
		Index &operator()( Index v1 , Index v2 , InitializationFunction &initializationFunction )
		{
			auto iter = _edgeTable.find( _EdgeKey(v1,v2) );
			if( iter==_edgeTable.end() )
			{
				Index idx = initializationFunction();
				_edgeTable[ EdgeKey(v1,v2) ] = idx;
				return idx;
			}
			else return iter->second;
		};
	};

	/////////////////
	// Permutation //
	/////////////////
	// This class represents a permutation, which can be thought of as a map from the range [0,N) to itself
	template< unsigned int N >
	struct Permutation
	{
		// The identity permutation
		Permutation( void ){ for( unsigned int n=0 ; n<N ; n++ ) _indices[n] = n; }

		// The permutation that describes how a list of N elements should be sorted
		// Specifically, P[i] gives the index of the i-th element in the sorted list
		Permutation( std::function< bool ( unsigned int , unsigned int ) > SortFunction )
		{
			unsigned int indices[N];
			for( unsigned int n=0 ; n<N ; n++ ) indices[n] = n;
			std::sort( indices , indices+N , SortFunction );
			for( unsigned int n=0 ; n<N ; n++ ) _indices[ indices[n] ] = n;
		}

		// Returns the index of the i-th element after applying the permtuation
		unsigned int operator[]( unsigned int i ) const { return _indices[i]; }

		// Returns the inverse of a permutation
		Permutation inverse( void ) const
		{
			Permutation p;
			for( unsigned int n=0 ; n<N ; n++ ) p._indices[ _indices[n] ] = n;
			return p;
		}

		// Returns the composition of two permutations
		Permutation operator * ( Permutation p ) const
		{
			Permutation _p;
			for( unsigned int n=0 ; n<N ; n++ ) _p._indices[n] = _indices[ p[n] ];
			return _p;
		}

		// Returns the matrix form of a permutation
		template< typename Real >
		SquareMatrix< Real , N > toMatrix( void ) const
		{
			SquareMatrix< Real , N > p;
			for( unsigned int n=0 ; n<N ; n++ ) p(n,_indices[n]) = 1;
			return p;
		};

		// Transposes the i-th and j-th elements
		void transpose( unsigned int i , unsigned int j ){ std::swap< unsigned int >( _indices[i] , _indices[j] ); }

		// Returns 0/1 depending on whether the parity of the permutation is even/odd
		unsigned int parity( void ) const
		{
			unsigned int count = 0;
			for( unsigned int i=0 ; i<N ; i++ ) for( unsigned int j=0 ; j<i ; j++ ) if( _indices[i]>_indices[j] ) count++;
			return count&1;
		}

		friend std::ostream &operator << ( std::ostream &os , const Permutation &p )
		{
			os << "{ ";
			for( int n=0 ; n<N ; n++ )
			{
				if( n ) os << " , ";
				os << p._indices[n];
			}
			return os << " }";
		}

	protected:
		unsigned int _indices[N];
	};

	//////////////////
	// SimplexIndex //
	//////////////////
	//template< unsigned int K , typename Index=unsigned int >
	template< unsigned int K , typename Index >
	struct SimplexIndex
	{
		Index idx[K+1];
		template< class ... Ints >
		SimplexIndex( Ints ... values ){ static_assert( sizeof...(values)==K+1 || sizeof...(values)==0 , "[ERROR] Invalid number of coefficients" ) ; _init( 0 , (Index)values ... ); }
		Index &operator[] ( unsigned int i ) { return idx[i] ;}
		const Index &operator[] ( unsigned int i ) const { return idx[i]; }
		template< typename Real , typename Vertex >
		void split( const Real values[K+1] , std::vector< Vertex > &vertices , EdgeTable< Index > &edgeTable , std::vector< SimplexIndex >& back , std::vector< SimplexIndex >& front ) const;
		template< typename ... UInts >
		SimplexIndex< K - (unsigned int)sizeof...( UInts ) - 1 > face( unsigned int faceIndex , UInts ... faceIndices ) const;
		//	template< typename ... UInts >
		//	SimplexIndex< K - (unsigned int)sizeof...( UInts ) - 1 > face( unsigned int faceIndex , UInts ... faceIndices ) const { bool oriented ; return face( oriented , faceIndex , faceIndices... ); }

		// Invokes the function on each of the _K-dimensional faces
		template< unsigned int _K , typename FaceFunctor /* = std::function< void ( SimplexIndex< _K , Index > )*/ >
		void processFaces( FaceFunctor F ) const;

		// Invokes the function on each of the _K-dimensional faces
		template< unsigned int _K , typename FaceFunctor /* = std::function< void ( SimplexIndex< _K , Index > )*/ >
		static void ProcessFaces( FaceFunctor F );

		// Sorts the indices and returns a boolean indicating if the permutation is even
		bool sort( void );
		bool sort( const Index indices[] );

		Permutation< K+1 > getPermutation( const Index indices[] ) const { return Permutation< K+1 >( [&]( unsigned int i1 , unsigned int i2 ){ return indices[ idx[i1] ]<indices[ idx[i2] ]; } ); }

		static SimplexIndex Canonical( void )
		{
			SimplexIndex< K , Index > si;
			for( unsigned int k=0 ; k<=K ; k++ ) si[k] = k;
			return si;
		}
		template< typename ... UInts >
		static SimplexIndex< K - (unsigned int)sizeof...( UInts ) - 1 , Index > Face( unsigned int faceIndex , UInts ... faceIndices );
		//	template< typename ... UInts >
		//	static SimplexIndex< K - (unsigned int)sizeof...( UInts ) - 1 , Index > Face( unsigned int faceIndex , UInts ... faceIndices ){ bool oriented ; return Face( oriented , faceIndex , faceIndices... ); }

		bool operator < ( const SimplexIndex &si ) const
		{
			for( unsigned int k=0 ; k<=K ; k++ ) if( idx[k]!=si.idx[k] ) return idx[k]<si.idx[k];
			return false;
		}
		bool operator == ( const SimplexIndex &si ) const
		{
			for( unsigned int k=0 ; k<=K ; k++ ) if( idx[k]!=si.idx[k] ) return false;
			return true;
		}
		bool operator != ( const SimplexIndex &si ) const { return !( operator==(si) ); }

		struct Hasher{ std::size_t operator()( const SimplexIndex & si ) const { return static_cast< std::size_t >( si[0] ); } };

	protected:
		SimplexIndex< K-1 , Index > _face( bool &oriented , unsigned int faceIndex ) const;

		static SimplexIndex< K-1 , Index > _Face( bool &oriented , unsigned int k );

		template< unsigned int _K , typename ... UInts , typename FaceFunctor /* = std::function< void ( SimplexIndex< _K , Index > )*/ >
		void _processFaces( FaceFunctor F , unsigned int faceIndex , UInts ... faceIndices ) const;
		void _init( unsigned int k )
		{
			if( !k ) for( unsigned int k=0 ; k<=K ; k++ ) idx[k] = static_cast< Index >(k);
			else MK_ERROR_OUT( "Should never be called" );
		}
		template< class ... Ints > void _init( unsigned int k , Index v , Ints ... values )
		{
			idx[k] = v;
			if( k+1<=K ) _init( k+1 , values ... );
		}

		friend std::ostream &operator << ( std::ostream &os , const SimplexIndex &s )
		{
			os << "{ ";
			for( int d=0 ; d<=K ; d++ )
			{
				if( d ) os << " , ";
				os << s[d];
			}
			return os << " }";
		}
	};

	template< typename Index >
	struct SimplexIndex< 0 , Index >
	{
		Index idx[1];
		SimplexIndex( Index i=0 ){ idx[0]=i; }
		Index &operator[] ( unsigned int i ) { return idx[i] ;}
		const Index &operator[] ( unsigned int i ) const { return idx[i]; }
		template< typename Real , typename Vertex >
		void split( const Real values[1] , std::vector< Vertex > &vertices , EdgeTable< Index > &edgeTable , std::vector< SimplexIndex >& back , std::vector< SimplexIndex >& front ) const
		{
			if( values[0]<0 ) back.push_back( *this );
			else              front.push_back( *this );
		}
		bool sort( void ){ return true; }
		bool sort( const Index indices[] ){ return true; }

		Permutation< 1 > getPermutation( const Index indices[] ) const { return Permutation<1>(); }

		// Invokes the function on each of the _K-dimensional faces
		template< unsigned int _K , typename FaceFunctor /* = std::function< void ( SimplexIndex< _K , Index > )*/ >
		void processFaces( FaceFunctor F ) const
		{
			static_assert( _K<=0 , "[ERROR] Sub-simplex dimension larger than simplex dimension" );
			F( *this );
		}

		// Invokes the function on each of the _K-dimensional faces
		template< unsigned int _K , typename FaceFunctor /* = std::function< void ( SimplexIndex< _K , Index > )*/ >
		static void ProcessFaces( FaceFunctor F )
		{
			SimplexIndex< 0 , Index > si;
			for( unsigned int k=0 ; k<=0 ; k++ ) si[k] = k;
			si.template processFaces< _K >( F );
		}

		bool operator <  ( const SimplexIndex &si ) const { return idx[0]< si.idx[0]; }
		bool operator == ( const SimplexIndex &si ) const { return idx[0]==si.idx[0]; }
		bool operator != ( const SimplexIndex &si ) const { return idx[0]!=si.idx[0]; }

		struct Hasher{ std::size_t operator()( const SimplexIndex & si ) const { return static_cast< std::size_t >( si[0] ); } };

	protected:
		friend std::ostream &operator << ( std::ostream &os , const SimplexIndex &s )
		{
			os << "{ ";
			for( int d=0 ; d<=0 ; d++ )
			{
				if( d ) os << " , ";
				os << s[d];
			}
			return os << " }";
		}
	};

	template< typename Real , unsigned int Dim , unsigned int K >
	struct SimplicialComplex
	{
		SimplicialComplex( const std::vector< Simplex< Real , Dim , K > > &simplices ) : _simplices( simplices ){}
		virtual size_t size( void ) const { return _simplices.size(); }
		virtual Simplex< Real , Dim , K > operator[]( size_t idx ) const { return _simplices[idx]; }
	protected:
		SimplicialComplex( void ) : _simplices(__simplices) {}
		const std::vector< Simplex< Real , Dim , K > > &_simplices;
		const std::vector< Simplex< Real , Dim , K > > __simplices;
	};

	template< typename Real , unsigned int Dim , unsigned int K , typename IndexType >
	struct IndexedSimplicialComplex : public SimplicialComplex< Real , Dim , K >
	{
		IndexedSimplicialComplex( const std::vector< Point< Real , Dim > > &vertices , const std::vector< SimplexIndex< K , IndexType > > &simplices ) : _vertices(vertices) , _simplices(simplices){}
		IndexedSimplicialComplex( IndexedSimplicialComplex && isc )
		{
			std::swap( _vertices , isc._vertices );
			std::swap( _simplices , isc._simplices );
		}

		size_t size( void ) const { return _simplices.size(); }
		Simplex< Real , Dim , K > operator[]( size_t idx ) const
		{
			Simplex< Real , Dim , K > s;
			for( unsigned int k=0 ; k<=K ; k++ ) s[k] = _vertices[ _simplices[idx][k] ];
			return s;
		}
	protected:
		const std::vector< Point< Real , Dim > > &_vertices;
		const std::vector< SimplexIndex< K , IndexType > > &_simplices;
	};

	template< typename Real , unsigned int Dim , unsigned int K >
	struct TransformedSimplicialComplex : public SimplicialComplex< Real , Dim , K >
	{
		TransformedSimplicialComplex( const SimplicialComplex< Real , Dim , K > &simplicialComplex , const XForm< Real , Dim+1 > &xForm ) : _simplicialComplex(simplicialComplex) , _xForm(xForm){}
		size_t size( void ) const { return _simplicialComplex.size(); }
		Simplex< Real , Dim , K > operator[]( size_t idx ) const { return _xForm * _simplicialComplex[idx]; }
	protected:
		const SimplicialComplex< Real , Dim , K > &_simplicialComplex;
		XForm< Real , Dim+1 > _xForm;
	};


	template< typename Real , unsigned int Dim > Point< Real , Dim > RandomBallPoint( void );
	template< typename Real , unsigned int Dim > Point< Real , Dim > RandomSpherePoint( void );
	template< typename Real , unsigned int Dim > Point< Real , Dim > RandomSimplexPoint( void );
	template< typename Real , unsigned int K   > Point< Real , K+1 > RandomBarycentricCoordinates( void );
	template< typename Real , unsigned int Dim > SquareMatrix< double , Dim > RandomRotationMatrix( void );

	template<class Real>
	XForm3x3<Real> RotationMatrix( Real a , Real b , Real c , Real d );

	template<class Real>
	XForm3x3<Real> RotationMatrix( const Point3D<Real>& axis , const Real& angle );

	typedef SimplexIndex< 1 , int > EdgeIndex;
	typedef SimplexIndex< 2 , int > TriangleIndex;
	typedef SimplexIndex< 3 , int > TetrahedronIndex;

	template< class Real > Point3D< Real > NearestPointOnTriangle( Point3D< Real > point , const Point3D< Real > triangle[3] , Real* b );
	template< class Real > Point3D< Real > NearestPointOnEdge( Point3D< Real > point , const Point3D< Real > edge[2] , Real& b0 , Real& b1 );

#ifdef NEW_MAT_CODE
	struct MinimalAreaTriangulation
	{
		template< class Real , unsigned int Dim >
		static double GetArea( const std::vector< Point< Real , Dim > > &vertices );

		template< typename AreaFunctor /*=std::function< double (unsigned int , unsigned int ,  unsigned int ) > */ >
		static double GetArea( const AreaFunctor & AF , unsigned int vNum );

		template< class Real , unsigned int Dim , typename Index >
		static void GetTriangulation( const std::vector< Point< Real , Dim > > &vertices , std::vector< SimplexIndex< 2 , Index > > &triangles );

		template< typename Index , typename AreaFunctor /*=std::function< double (unsigned int , unsigned int ,  unsigned int ) > */ >
		static void GetTriangulation( const AreaFunctor & AF , unsigned int vNum , std::vector< SimplexIndex< 2 , Index > > &triangles );

	protected:
		template< class Real , unsigned int Dim >
		static double _Area( Point< Real , Dim > v0 , Point< Real , Dim > v1 , Point< Real , Dim > v2 );

		MinimalAreaTriangulation( void ) : _bestTriangulation(nullptr) , _midPoint(nullptr){}
		~MinimalAreaTriangulation( void ) { delete[] _bestTriangulation ; delete[] _midPoint; }

		double *_bestTriangulation;
		unsigned int *_midPoint;

		template< typename AreaFunctor /*=std::function< double (unsigned int , unsigned int ,  unsigned int ) > */ >
		double _getArea( unsigned int i , unsigned int j , const AreaFunctor & AF , unsigned int vNum );

		template< typename Index , typename AreaFunctor /*=std::function< double (unsigned int , unsigned int ,  unsigned int ) > */ >
		void _getTriangulation( unsigned int i , unsigned int j , const AreaFunctor & AF , unsigned int vNum , std::vector< SimplexIndex< 2 , Index > > &triangles , unsigned int &idx );
	};
#else // !NEW_MAT_CODE
	template< class Real , unsigned int Dim >
	struct MinimalAreaTriangulation
	{
		static double GetArea( const std::vector< Point< Real , Dim > > &vertices );

		template< typename AreaFunctor /*=std::function< double (unsigned int , unsigned int ,  unsigned int ) > */ >
		static double GetArea( const AreaFunctor & AF , unsigned int vNum );

		template< typename Index >
		static void GetTriangulation( const std::vector< Point< Real , Dim > > &vertices , std::vector< SimplexIndex< 2 , Index > > &triangles );

		template< typename Index , typename AreaFunctor /*=std::function< double (unsigned int , unsigned int ,  unsigned int ) > */ >
		static void GetTriangulation( const AreaFunctor & AF , unsigned int vNum , std::vector< SimplexIndex< 2 , Index > > &triangles );

	protected:
		static double _Area( Point< Real , Dim > v0 , Point< Real , Dim > v1 , Point< Real , Dim > v2 );

		MinimalAreaTriangulation( void );
		~MinimalAreaTriangulation( void );

		double *_bestTriangulation;
		size_t *_midPoint;

		template< typename AreaFunctor /*=std::function< double (unsigned int , unsigned int ,  unsigned int ) > */ >
		double _getArea( size_t i , size_t j , const AreaFunctor & AF , unsigned int vNum );

		double _getArea( size_t i , size_t j , const std::vector< Point< Real , Dim > > &vertices );


		template< typename Index , typename AreaFunctor /*=std::function< double (unsigned int , unsigned int ,  unsigned int ) > */ >
		void _getTriangulation( size_t i , size_t j , const AreaFunctor & AF , unsigned int vNum , std::vector< SimplexIndex< 2 , Index > > &triangles , size_t &idx );

		template< typename Index >
		void _getTriangulation( size_t i , size_t j , const std::vector< Point< Real , Dim > > &vertices , std::vector< SimplexIndex< 2 , Index > > &triangles , size_t &idx );
	};
#endif // NEW_MAT_CODE

	struct EarTriangulation
	{
		template< typename Index , typename Real >
		static void GetTriangulation( const std::vector< Point< Real , 2 > > &vertices , std::vector< SimplexIndex< 2 , Index > > &triangles );
	};
#include "Geometry.inl"
}

#endif // GEOMETRY_INCLUDED
