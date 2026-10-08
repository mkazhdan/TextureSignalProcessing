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

#ifndef POLYNOMIAL_INCLUDED
#define POLYNOMIAL_INCLUDED
#include <iostream>
#include <array>
#include <type_traits>

#include "Geometry.h"
#include "Algebra.h"
#include "Poly34.h"
#include "Exceptions.h"

#define NEW_POLYNOMIAL // Consolidates the general and base classes

namespace MishaK
{
	namespace Polynomial
	{
#ifdef NEW_POLYNOMIAL
		/** Helper functionality for computing the minimum of two integers.*/
		template< unsigned int D1 , unsigned int D2 > constexpr unsigned int Max( void ){ if constexpr( D1>D2 ) return D1 ; else return D2; };

		template< unsigned int Dim , unsigned int Degree >
		static constexpr unsigned int NumCoefficients( void )
		{
			if constexpr( Dim==0 || Degree==0 ) return 1;
			else return NumCoefficients< Dim-1 , Degree >() + NumCoefficients< Dim , Degree-1 >();
		}

		/** The generic, recursively defined, Polynomial class of total degree Degree. */
		template< unsigned int Dim , unsigned int Degree , typename T , typename Real=T >
		class Polynomial : public VectorSpace< Real , Polynomial< Dim , Degree , T , Real > >
		{
			template< unsigned int _Dim , unsigned int _Degree , typename _T , typename _Real > friend class Polynomial;
			template< unsigned int _Dim , unsigned int _Degree , typename _T , typename _Real > friend std::ostream &operator << ( std::ostream & , const Polynomial< _Dim , _Degree , _T , _Real > & );

			template< unsigned int _Dim , unsigned int Degree1 , unsigned int Degree2 , typename T1 , typename T2 , typename _Real > friend Polynomial< _Dim , Degree1 + Degree2 , decltype( std::declval<T1>() * std::declval<T2>() ) , _Real > operator * ( const Polynomial< _Dim , Degree1 , T1 , _Real > & , const Polynomial< _Dim , Degree2 , T2 , _Real > & );
			template< unsigned int _Dim , unsigned int Degree1 , unsigned int Degree2 , typename _T , typename _Real > friend Polynomial< _Dim , Max< Degree1 , Degree2 >() , _Real > operator + ( const Polynomial< _Dim , Degree1 , _T , _Real > & , const Polynomial< _Dim , Degree2 , _T , _Real > & );
			template< unsigned int _Dim , unsigned int Degree1 , unsigned int Degree2 , typename _T , typename _Real > friend Polynomial< _Dim , Max< Degree1 , Degree2 >() , _Real > operator - ( const Polynomial< _Dim , Degree1 , _T , _Real > & , const Polynomial< _Dim , Degree2 , _T , _Real > & );


			/** The polynomials in Dim-1 dimensions.
			*** The total polynomial is assumed to be _polynomials[0] + _polynomials[1] * (x_Dim) + _polynomials[2] * (x_Dim)^2 + ... */
			std::conditional_t< Dim!=0 , std::array< Polynomial< Dim-1 , Degree , T , Real > , Degree+1 > , std::array< T , 1 > > _coefficients;

			/** This method returns the specified coefficient of the polynomial.*/
			const T & _coefficient( const unsigned int indices[] , unsigned int maxDegree ) const;

			/** This method returns the specified coefficient of the polynomial.*/
			T & _coefficient( const unsigned int indices[] , unsigned int maxDegree );

			/** This method evaluates the polynomial at the specified set of coordinates.*/
			T _evaluate( const Real coordinates[] , unsigned int maxDegree ) const;

			/** This method computes the pull-back of the polynomial given the map from a _Dim-dimensional space to the Dim-dimensional space given by x -> A({x,1}). */	
			template< unsigned int _Dim >
			Polynomial< _Dim , Degree , T , Real > _pullBack( const Matrix< Real , _Dim+1 , Dim > &A , unsigned int maxDegree ) const;

			/** This method returns true if the polynomial is zero. */
			bool _isZero( unsigned int maxDegree ) const;

			/** This method returns true if the polynomial is a constant. */
			bool _isConstant( unsigned int maxDegree ) const;

			bool _print( std::ostream &ostream , const std::string varNames[] , bool first ) const;
			bool _print( std::ostream &ostream , const std::string varNames[] , std::string suffix , bool first ) const;

			unsigned int _setCoefficients( const T *coefficients , unsigned int maxDegree );
			unsigned int _getCoefficients(       T *coefficients , unsigned int maxDegree ) const;

			template< unsigned int _D=0 >
			static void _SetDegrees( unsigned int coefficientIndex , unsigned int degrees[/*Dim*/] );

			T _integrateUnitRightSimplex( void ) const;

		public:
			/** The number of coefficients in the polynomial. */
			static const unsigned int NumCoefficients = NumCoefficients< Dim , Degree >();

			/** The partial derivative type. */
			using DerivativeType = Polynomial< Dim , (Degree>1) ? Degree-1 : 0 , Point< T , Dim , Real > , Real >;

			/** The partial derivative type. */
			using PartialDerivativeType = Polynomial< Dim , (Degree>1) ? Degree-1 : 0 , T , Real >;

			/** The partial derivative type. */
			using PartialIntegralType = Polynomial< Dim , Degree+1 , T , Real >;

			/** The default constructor initializes the coefficients to zero.*/
			Polynomial( void );

			/** This constructor creates a constant polynomial */
			Polynomial( T c );

			/** This constructor initializes the coefficients. */
			Polynomial( const T coefficients[NumCoefficients] );

			/** The constructor copies over as much of the polynomial as will fit.*/
			template< unsigned int _Degree >
			Polynomial( const Polynomial< Dim , _Degree , T , Real > &p );

			/** The equality operator copies over as much of the polynomial as will fit.*/
			template< unsigned int _Degree >
			Polynomial & operator = ( const Polynomial< Dim , _Degree , T , Real > &p );

			/** This method returns the associated coefficient of the polynomial */
			template< typename ... UnsignedInts > requires( std::is_integral_v< UnsignedInts > && ... )
			const T & coefficient( UnsignedInts ... indices ) const;

			/** This method returns the associated coefficient of the polynomial */
			template< typename ... UnsignedInts > requires( std::is_integral_v< UnsignedInts > && ... )
			T & coefficient( UnsignedInts ... indices );

			/** This method returns the associated coefficient of the polynomial */
			T &coefficient( const unsigned int indices[ /*Dim*/ ] );

			/** This method returns the associated coefficient of the polynomial */
			const T & coefficient( const unsigned int indices[ /*Dim*/ ] ) const;

			/** This method returns the associated coefficient of the polynomial */
			T & operator[]( unsigned int idx );

			/** This method returns the associated coefficient of the polynomial */
			const T & operator[]( unsigned int idx ) const;

			/** This static method sets the degrees of the associated coefficient */
			static void SetDegrees( unsigned int coefficientIndex , unsigned int degrees[ /*Dim*/ ] );

			/** This method returns a vector of the coefficients of the polynomial */
			Point< T , Polynomial< Dim , Degree , T , Real >::NumCoefficients , Real > coefficients( void ) const;

			/** This method evaluates the polynomial at the prescribed point.*/
			template< typename ... Reals > requires( std::is_convertible_v< Reals , Real > && ... )
			T operator()( Reals ... coordinates ) const;

			/** This method evaluates the polynomial at the prescribed point.*/
			T operator()( Point< Real , Dim > p ) const;

			/** This method evaluates the gradient of the polynomial at the prescribed point.*/
			template< typename ... Reals > requires( std::is_convertible_v< Reals , Real > && ... )
			Point< T , Dim , Real > gradient( Reals ... coordinates ) const;

			/** This method evaluates the gradient of the  polynomial at the prescribed point.*/
			Point< T , Dim , Real > gradient( Point< Real , Dim > p ) const;

			/** This method evaluates the Hessian of the polynomial at the prescribed point.*/
			template< typename ... Reals > requires( std::is_convertible_v< Reals , Real > && ... )
			SquareMatrix< T , Dim > hessian( Reals ... coordinates ) const;

			/** This method evaluates the Hessian of the polynomial at the prescribed point.*/
			SquareMatrix< T , Dim > hessian( Point< Real , Dim > p ) const;

			/** This method returns the partial derivative with respect to the prescribed dimension.*/
			PartialDerivativeType d( unsigned int dim ) const;

			/** This method returns the gradient/differential of the polynomial.*/
			DerivativeType d( void ) const;

			/** This method returns the integral of the polynomial w.r.t. to a single parameter.*/
			PartialIntegralType integral( unsigned int dim ) const;

			/** This method computes the pull-back of the polynomial given the map from a (_Dim-1)-dimensional space to the Dim-dimensional space given by x -> A({x,1}). */	
			template< unsigned int _Dim >
			Polynomial< _Dim-1 , Degree , T , Real > operator()( Matrix< Real , _Dim , Dim > A ) const;
			template< unsigned int _Dim >
			Polynomial< _Dim-1 , Degree , T , Real > pullBack( Matrix< Real , _Dim , Dim > A ) const;
			Polynomial< 1 , Degree , T , Real > operator()( const Ray< Real , Dim > &ray ) const;

			/** Integrate the polynomial over a unit cube */
			T integrateUnitCube( void ) const;

			/** Integrate the polynomial over a unit right simplex */
			T integrateUnitRightSimplex( void ) const;

			/** Returns the matrix taking in the ccoefficients of the polynomial and returning the values at the prescribed points */
			static Matrix< Real , Polynomial< Dim , Degree , T , Real >::NumCoefficients , Polynomial< Dim , Degree , T , Real >::NumCoefficients > EvaluationMatrix( const Point< Real , Dim > positions[NumCoefficients] );

			/** Shifts the polynomial */
			Polynomial shift( Point< Real , Dim > s ) const;

			/////////////////////////
			// VectorSpace methods //
			/////////////////////////
			/** This method scales the polynomial */
			void Scale( Real s );

			/** This method adds in the polynomial */
			void Add( const Polynomial & p );
		};

		/** This function returns the sum of two polynomials. */
		template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 , typename T , typename Real >
		Polynomial< Dim , Max< Degree1 , Degree2 >() , T , Real > operator + ( const Polynomial< Dim , Degree1 , T , Real > &p1 , const Polynomial< Dim , Degree2 , T , Real > &p2 );

		/** This function returns the difference of two polynomials. */
		template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 , typename T , typename Real >
		Polynomial< Dim , Max< Degree1 , Degree2 >() , T , Real > operator - ( const Polynomial< Dim , Degree1 , T , Real > &p1 , const Polynomial< Dim , Degree2 , T , Real > &p2 );

		/** This function prints out the polynomial.*/
		template< unsigned int Dim , unsigned int Degree , typename T , typename Real >
		std::ostream &operator << ( std::ostream &stream , const Polynomial< Dim , Degree , T , Real > &poly );

		/** This function returns the product of two polynomials. */
		template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 , typename T1 , typename T2 , typename Real >
		Polynomial< Dim , Degree1 + Degree2 , decltype( std::declval< T1 >() * std::declval< T2 >() ) , Real > operator * ( const Polynomial< Dim , Degree1 , T1 , Real > & , const Polynomial< Dim , Degree2 , T2 , Real > & );

#pragma message( "[WARNING] When the first argument is also a polynomial, multiplication is ambiguous. Depending on the compiler to choose the more specialized interpretation of mutiplying two polynomials." )
		/** This function applies a transformation to a polynomial by applying the transformation to the coefficients. */
		template< typename XForm , unsigned int Dim , unsigned int Degree , typename T , typename Real > 
		Polynomial< Dim , Degree , decltype( std::declval< XForm >() * std::declval< T >() ) , Real > operator * ( const XForm & xForm , const Polynomial< Dim , Degree , T , Real > & p );

		/** This function applies a polynomial of transformations to an argument by applying the transformation to the coefficients. */
		template< unsigned int Dim , unsigned int Degree , typename T , typename Real , typename Arg > 
		Polynomial< Dim , Degree , decltype( std::declval< T >() * std::declval< Arg >() ) , Real > operator * ( const Polynomial< Dim , Degree , T , Real > & p , const Arg & arg );

#else // !NEW_POLYNOMIAL
		/** Helper functionality for computing the minimum of two integers.*/
		template< unsigned int D1 , unsigned int D2 > struct Min{ static const unsigned int Value = D1<D2 ? D1 : D2; };
		/** Helper functionality for computing the maximum of two integers.*/
		template< unsigned int D1 , unsigned int D2 > struct Max{ static const unsigned int Value = D1>D2 ? D1 : D2; };

		template< unsigned int K >
		constexpr unsigned int Factorial( void )
		{
			if constexpr( K==0 ) return 1;
			else return Factorial< K-1 >() * K;
		}

/** The generic, recursively defined, Polynomial class of total degree Degree. */
		template< unsigned int Dim , unsigned int Degree , typename T , typename Real=T >
		class Polynomial : public VectorSpace< Real , Polynomial< Dim , Degree , T , Real > >
		{
			template< unsigned int _Dim , unsigned int _Degree , typename _T , typename _Real > friend class Polynomial;
			template< unsigned int _Dim , unsigned int _Degree , typename _T , typename _Real > friend std::ostream &operator << ( std::ostream & , const Polynomial< _Dim , _Degree , _T , _Real > & );

			template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 , typename T1 , typename T2 , typename Real > friend Polynomial< Dim , Degree1 + Degree2 , decltype( std::declval<T1>() * std::declval<T2>() ) , Real > operator * ( const Polynomial< Dim , Degree1 , T1 , Real > & , const Polynomial< Dim , Degree2 , T2 , Real > & );
			template< unsigned int _Dim , unsigned int Degree1 , unsigned int Degree2 , typename _T , typename _Real > friend Polynomial< _Dim , Max< Degree1 , Degree2 >::Value , _Real > operator + ( const Polynomial< _Dim , Degree1 , _T , _Real > & , const Polynomial< _Dim , Degree2 , _T , _Real > & );
			template< unsigned int _Dim , unsigned int Degree1 , unsigned int Degree2 , typename _T , typename _Real > friend Polynomial< _Dim , Max< Degree1 , Degree2 >::Value , _Real > operator - ( const Polynomial< _Dim , Degree1 , _T , _Real > & , const Polynomial< _Dim , Degree2 , _T , _Real > & );


			/** The polynomials in Dim-1 dimensions.
			*** The total polynomial is assumed to be _polynomials[0] + _polynomials[1] * (x_Dim) + _polynomials[2] * (x_Dim)^2 + ... */
			std::array< Polynomial< Dim-1 , Degree , T , Real > , Degree+1 > _polynomials;

			/** This method returns the specified coefficient of the polynomial.*/
			const T &_coefficient( const unsigned int indices[] , unsigned int maxDegree ) const;

			/** This method returns the specified coefficient of the polynomial.*/
			T &_coefficient( const unsigned int indices[] , unsigned int maxDegree );

			/** This method evaluates the polynomial at the specified set of coordinates.*/
			T _evaluate( const Real coordinates[] , unsigned int maxDegree ) const;

			/** This method computes the pull-back of the polynomial given the map from a _Dim-dimensional space to the Dim-dimensional space given by x -> A({x,1}). */	
			template< unsigned int _Dim >
			Polynomial< _Dim , Degree , T , Real > _pullBack( const Matrix< Real , _Dim+1 , Dim > &A , unsigned int maxDegree ) const;

			/** This method returns true if the polynomial is zero. */
			bool _isZero( unsigned int maxDegree ) const;

			/** This method returns true if the polynomial is a constant. */
			bool _isConstant( unsigned int maxDegree ) const;

			bool _print( std::ostream &ostream , const std::string varNames[] , bool first ) const;
			bool _print( std::ostream &ostream , const std::string varNames[] , std::string suffix , bool first ) const;

			unsigned int _setCoefficients( const T *coefficients , unsigned int maxDegree );
			unsigned int _getCoefficients(       T *coefficients , unsigned int maxDegree ) const;

			template< unsigned int _D=0 >
			static void _SetDegrees( unsigned int coefficientIndex , unsigned int degrees[Dim] );

			T _integrateUnitRightSimplex( void ) const;

			static constexpr unsigned int _NumCoefficients( void )
			{
				if constexpr( Degree ) return Polynomial< Dim-1 , Degree , T , Real >::NumCoefficients + Polynomial< Dim , Degree-1 , Real >::NumCoefficients;
				else                   return Polynomial< Dim-1 , Degree , T , Real >::NumCoefficients;
			}
		public:
			/** The number of coefficients in the polynomial. */
			static const unsigned int NumCoefficients = _NumCoefficients();

			/** The partial derivative type. */
			using DerivativeType = Polynomial< Dim , (Degree>1) ? Degree-1 : 0 , Point< T , Dim , Real > , Real >;
				
			/** The partial derivative type. */
			using PartialDerivativeType = Polynomial< Dim , (Degree>1) ? Degree-1 : 0 , T , Real >;

			/** The partial derivative type. */
			using PartialIntegralType = Polynomial< Dim , Degree+1 , T , Real >;

			/** The default constructor initializes the coefficients to zero.*/
			Polynomial( void );

			/** This constructor creates a constant polynomial */
			Polynomial( T c );

			/** This constructor initializes the coefficients. */
			Polynomial( const T coefficients[NumCoefficients] );

			/** The constructor copies over as much of the polynomial as will fit.*/
			template< unsigned int _Degree >
			Polynomial( const Polynomial< Dim , _Degree , T , Real > &p );

			/** The equality operator copies over as much of the polynomial as will fit.*/
			template< unsigned int _Degree >
			Polynomial &operator = ( const Polynomial< Dim , _Degree , T , Real > &p );

			/** This method returns the associated coefficient of the polynomial */
			template< typename ... UnsignedInts >
			const T &coefficient( unsigned int idx , UnsignedInts ... indices ) const;

			/** This method returns the associated coefficient of the polynomial */
			template< typename ... UnsignedInts >
			T &coefficient( unsigned int idx , UnsignedInts ... indices );

			/** This method returns the associated coefficient of the polynomial */
			template< typename ... Ints >
			const T &coefficient( int idx , Ints ... indices ) const;

			/** This method returns the associated coefficient of the polynomial */
			template< typename ... Ints >
			T &coefficient( int idx , Ints ... indices );

			/** This method returns the associated coefficient of the polynomial */
			T &coefficient( const unsigned int indices[ Dim ] );

			/** This method returns the associated coefficient of the polynomial */
			const T &coefficient( const unsigned int indices[ Dim ] ) const;

			/** This method returns the associated coefficient of the polynomial */
			T &operator[]( unsigned int idx );

			/** This method returns the associated coefficient of the polynomial */
			const T &operator[]( unsigned int idx ) const;

			/** This static method sets the degrees of the associated coefficient */
			static void SetDegrees( unsigned int coefficientIndex , unsigned int degrees[Dim] );

			/** This method returns a vector of the coefficients of the polynomial */
			Point< T , Polynomial< Dim , Degree , T , Real >::NumCoefficients , Real > coefficients( void ) const;

			/** This method evaluates the polynomial at the prescribed point.*/
			template< typename ... Reals >
			T operator()( Real coord , Reals ... coordinates ) const;

			/** This method evaluates the polynomial at the prescribed point.*/
			T operator()( Point< Real , Dim > p ) const;

			/** This method evaluates the gradient of the polynomial at the prescribed point.*/
			template< typename ... Reals >
			Point< T , Dim , Real > gradient( Reals ... coordinates ) const;

			/** This method evaluates the gradient of the  polynomial at the prescribed point.*/
			Point< T , Dim , Real > gradient( Point< Real , Dim > p ) const;

			/** This method evaluates the Hessian of the polynomial at the prescribed point.*/
			template< typename ... Reals >
			SquareMatrix< T , Dim > hessian( Reals ... coordinates ) const;

			/** This method evaluates the Hessian of the polynomial at the prescribed point.*/
			SquareMatrix< T , Dim > hessian( Point< Real , Dim > p ) const;

			/** This method returns the partial derivative with respect to the prescribed dimension.*/
			PartialDerivativeType d( unsigned int dim ) const;

			/** This method returns the gradient/differential of the polynomial.*/
			DerivativeType d( void ) const;

			/** This method returns the integral of the polynomial w.r.t. to a single parameter.*/
			PartialIntegralType integral( unsigned int dim ) const;

			/** This method computes the pull-back of the polynomial given the map from a (_Dim-1)-dimensional space to the Dim-dimensional space given by x -> A({x,1}). */	
			template< unsigned int _Dim >
			Polynomial< _Dim-1 , Degree , T , Real > operator()( Matrix< Real , _Dim , Dim > A ) const;
			template< unsigned int _Dim >
			Polynomial< _Dim-1 , Degree , T , Real > pullBack( Matrix< Real , _Dim , Dim > A ) const;
			Polynomial< 1 , Degree , T , Real > operator()( const Ray< Real , Dim > &ray ) const;

			/** Integrate the polynomial over a unit cube */
			T integrateUnitCube( void ) const;

			/** Integrate the polynomial over a unit right simplex */
			T integrateUnitRightSimplex( void ) const;

			/** Returns the matrix taking in the ccoefficients of the polynomial and returning the values at the prescribed points */
			static Matrix< Real , Polynomial< Dim , Degree , T , Real >::NumCoefficients , Polynomial< Dim , Degree , T , Real >::NumCoefficients > EvaluationMatrix( const Point< Real , Dim > positions[NumCoefficients] );

			/** Shifts the polynomial */
			Polynomial shift( Point< Real , Dim > s ) const;

			/////////////////////////
			// VectorSpace methods //
			/////////////////////////
			/** This method scales the polynomial */
			void Scale( Real s );

			/** This method adds in the polynomial */
			void Add( const Polynomial & p );
		};

		/** This function returns the sum of two polynomials. */
		template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 , typename T , typename Real >
		Polynomial< Dim , Max< Degree1 , Degree2 >::Value , T , Real > operator + ( const Polynomial< Dim , Degree1 , T , Real > &p1 , const Polynomial< Dim , Degree2 , T , Real > &p2 );

		/** This function returns the difference of two polynomials. */
		template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 , typename T , typename Real >
		Polynomial< Dim , Max< Degree1 , Degree2 >::Value , T , Real > operator - ( const Polynomial< Dim , Degree1 , T , Real > &p1 , const Polynomial< Dim , Degree2 , T , Real > &p2 );

		/** This function prints out the polynomial.*/
		template< unsigned int Dim , unsigned int Degree , typename T , typename Real >
		std::ostream &operator << ( std::ostream &stream , const Polynomial< Dim , Degree , T , Real > &poly );


		/** A specialized instance of the Polynomial class in one variable */
		template< unsigned int Degree , typename T , typename Real >
		class Polynomial< 0 , Degree , T , Real > : public VectorSpace< Real , Polynomial< 0 , Degree , T , Real > >
		{
#pragma message( "[WARNING] Should merge this into the generic polynomial class" )
		public:
			/** The number of variables.*/
			static const unsigned int Dim = 0;
		protected:
			template< unsigned int Dim , unsigned int _Degree , typename _T , typename _Real > friend class Polynomial;
			template<                    unsigned int _Degree , typename _T , typename _Real > friend std::ostream &operator << ( std::ostream & , const Polynomial< Dim , _Degree , _Real > & );

			template< unsigned int Degree1 , unsigned int Degree2 , typename _T , typename _Real > friend Polynomial< 0 , Max< Degree1 , Degree2 >::Value , _T , _Real > operator + ( const Polynomial< 0 , Degree1 , _T , _Real > & , const Polynomial< 0 , Degree2 , _T , _Real > & );
			template< unsigned int Degree1 , unsigned int Degree2 , typename _T , typename _Real > friend Polynomial< 0 , Max< Degree1 , Degree2 >::Value , _T , _Real > operator - ( const Polynomial< 0 , Degree1 , _T , _Real > & , const Polynomial< 0 , Degree2 , _T , _Real > & );

		protected:
			/** The coefficients of the polynomial. */
			T _coefficients[1];

			/** This method returns the specified coefficient of the polynomial.*/
			const T &_coefficient( const unsigned int indices[] , unsigned int maxDegree ) const;

			/** This method returns the specified coefficient of the polynomial.*/
			T &_coefficient( const unsigned int indices[] , unsigned int maxDegree );

			/** This method evaluates the polynomial at the specified set of coordinates.*/
			T _evaluate( const Real coordinates[] , unsigned int maxDegree ) const;

			/** This method computes the pull-back of the polynomial given the map from a _Dim-dimensional space to the 1-dimensional space given by x -> A({x,1}). */	
			template< unsigned int _Dim >
			Polynomial< _Dim , Degree , T , Real > _pullBack( const Matrix< Real , _Dim+1 , 0 > &A , unsigned int maxDegree ) const;

			/** This method returns true if the polynomial is zero. */
			bool _isZero( unsigned int maxDegree ) const;

			/** This method returns true if the polynomial is a constant. */
			bool _isConstant( unsigned int maxDegree ) const;

			bool _print( std::ostream &ostream , const std::string varNames[] , bool first ) const;
			bool _print( std::ostream &ostream , const std::string varNames[] , std::string suffix , bool first ) const;

			unsigned int _setCoefficients( const T *coefficients , unsigned int maxDegree );
			unsigned int _getCoefficients(       T *coefficients , unsigned int maxDegree ) const;

			T _integrateUnitRightSimplex( void ) const;
		public:

			/** The number of coefficients in the polynomial. */
			static const unsigned int NumCoefficients = 1;

			/** The partial derivative type. */
			using PartialDerivativeType = Polynomial< Dim , (Degree>1) ? Degree-1 : 0 , T , Real >;

			/** The first derivative type. */
			using DerivativeType = Polynomial< Dim , (Degree>1) ? Degree-1 : 0 , Point< T , Dim , Real > , Real >;

			/** The integral type. */
			using PartialIntegralType = Polynomial< Dim , Degree+1 , T , Real >;

			/** The default constructor initializes the coefficients to zero.*/
			Polynomial( void );

			/** This constructor creates a constant polynomial */
			Polynomial( T c );

			/** This constructor initializes the coefficients (starting with lower degrees). */
			Polynomial( const T coefficients[NumCoefficients] );

			/** The constructor copies over as much of the polynomial as will fit.*/
			template< unsigned int _Degree >
			Polynomial( const Polynomial< 0 , _Degree , T , Real > &p );

			/** The equality operator copies over as much of the polynomial as will fit.*/
			template< unsigned int _Degree >
			Polynomial &operator= ( const Polynomial< 0 , _Degree , T , Real > &p );

			/** This method returns the coefficient.*/
			const T &coefficient( void ) const;

			/** This method returns the coefficient.*/
			T &coefficient( void );

			/** This method returns the associated coefficient of the polynomial */
			const T &coefficient( const unsigned int indices[0] ) const;

			/** This method returns the associated coefficient of the polynomial */
			T &coefficient( const unsigned int indices[0] );

			/** This method returns the associated coefficient of the polynomial */
			T &operator[]( unsigned int idx );

			/** This method returns the associated coefficient of the polynomial */
			const T &operator[]( unsigned int idx ) const;

			/** This static method sets the degrees of the associated coefficient */
			static void SetDegrees( unsigned int coefficientIndex , unsigned int degrees[0] );

			/** This method returns a vector of the coefficients of the polynomial */
			Point< T , Polynomial< 0 , Degree , T , Real >::NumCoefficients , Real > coefficients( void ) const;

			/** This method evaluates the polynomial at a given value.*/
			T operator()( void ) const;

			/** This method evaluates the polynomial at the prescribed point.*/
			T operator()( Point< Real , 0 > p ) const;

			/** This method returns the derivative of the polynomial.*/
			PartialDerivativeType d( unsigned int d=0 ) const;

			/** This method returns the gradient/differential of the polynomial.*/
			DerivativeType d( void ) const;

			/** This method returns the integral of the polynomial w.r.t. to a single parameter.*/
			PartialIntegralType integral( unsigned int dim ) const;

			/** This method computes the pull-back of the polynomial given the map from a (_Dim-1)-dimensional space to the 1-dimensional space given by x -> A({x,1}). */	
			template< unsigned int _Dim >
			Polynomial< _Dim-1 , Degree , T , Real > operator()( const Matrix< Real , _Dim , 0 > &A ) const;

			/** Integrate the polynomial over a unit interval */
			T integrateUnitCube( void ) const;

			/** Integrate the polynomial over a unit interval */
			T integrateUnitRightSimplex( void ) const;

			/** Returns the matrix taking in the ccoefficients of the polynomial and returning the values at the prescribed points */
			static Matrix< Real , Polynomial< 0 , Degree , T , Real >::NumCoefficients  , Polynomial< 0 , Degree , T , Real >::NumCoefficients> EvaluationMatrix( const Point< Real , 0 > positions[NumCoefficients] );

			/** Shifts the polynomial */
			Polynomial shift( Point< Real , Dim > s ) const;

			/////////////////////////
			// VectorSpace methods //
			/////////////////////////
			/** This method scales the polynomial*/
			void Scale( Real s );

			/** This method adds in the polynomial */
			void Add( const Polynomial & p );
		};

		/** This function returns the sum of two polynomials. */
		template< unsigned int Degree1 , unsigned int Degree2 , typename T , typename Real >
		Polynomial< 0 , Max< Degree1 , Degree2 >::Value , T , Real > operator + ( const Polynomial< 0 , Degree1 , T , Real > &p1 , const Polynomial< 0 , Degree2 , T , Real > &p2 );

		/** This function returns the difference of two polynomials. */
		template< unsigned int Degree1 , unsigned int Degree2 , typename T , typename Real >
		Polynomial< 0 , Max< Degree1 , Degree2 >::Value , T , Real > operator - ( const Polynomial< 0 , Degree1 , T , Real > &p1 , const Polynomial< 0 , Degree2 , T , Real > &p2 );

		/** This function returns the product of two polynomials, where one is scalar-valued. */
		template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 , typename T1 , typename T2 , typename Real >
		Polynomial< Dim , Degree1 + Degree2 , decltype( std::declval< T1 >() * std::declval< T2 >() ) , Real > operator * ( const Polynomial< Dim , Degree1 , T1 , Real > & , const Polynomial< Dim , Degree2 , T2 , Real > & );

		/** This function applies a transformation to a polynomial by applying the transformation to the coefficients. */
		template< typename XForm , unsigned int Dim , unsigned int Degree , typename T , typename Real > 
		Polynomial< Dim , Degree , decltype( std::declval< XForm >() * std::declval< T >() ) , Real > operator * ( const XForm & xForm , const Polynomial< Dim , Degree , T , Real > & p );
#endif // NEW_POLYNOMIAL

		////////////////////////////////////////////////
		// Classes specialized for 1D, 2D, 3D, and 4D //
		////////////////////////////////////////////////
		/** A polynomial in one variable of degree Degree */
		template< unsigned int Degree >
		using Polynomial1D = Polynomial< 1 , Degree , double >;

		/** A polynomial in two variables of degree Degree */
		template< unsigned int Degree >
		using Polynomial2D = Polynomial< 2 , Degree , double >;

		/** A polynomial in three variable of degree Degree */
		template< unsigned int Degree >
		using Polynomial3D = Polynomial< 3 , Degree , double >;

		/** A polynomial in four variable of degree Degree */
		template< unsigned int Degree >
		using Polynomial4D = Polynomial< 4 , Degree , double >;

		/** Shifts the polynomial */
		template< unsigned int Degree , typename T , typename Real >
		Polynomial< 1 , Degree , T , Real > Shift( const Polynomial< 1 , Degree , T , Real > & p , Real s );

		/** Computes the indefinite integral */
		template< unsigned int Degree , typename T , typename Real >
		Polynomial< 1 , Degree+1 , T , Real > Integral( const Polynomial< 1 , Degree , T , Real > & p );

		/** Sets the roots of the 1D polynomial and returns the number of roots set.
		*** The method is only specialized for degrees 1, 2, 3, and 4. */
		template< unsigned int Degree , typename Real >
		unsigned int Roots( const Polynomial< 1 , Degree , Real > &p , Real *r , double eps=1e-14 );

#include "Polynomial.inl"
	}
}
#endif // POLYNOMIAL_INCLUDED
