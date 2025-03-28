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
#include "Geometry.h"
#include "Algebra.h"
#include "Poly34.h"
#include "Exceptions.h"

namespace MishaK
{
	namespace Polynomial
	{
		/** Helper functionality for computing the minimum of two integers.*/
		template< unsigned int D1 , unsigned int D2 > struct Min{ static const unsigned int Value = D1<D2 ? D1 : D2; };
		/** Helper functionality for computing the maximum of two integers.*/
		template< unsigned int D1 , unsigned int D2 > struct Max{ static const unsigned int Value = D1>D2 ? D1 : D2; };

		/** The generic, recursively defined, Polynomial class of total degree Degree. */
		template< unsigned int Dim , unsigned int Degree , typename Real >
		class Polynomial : public VectorSpace< Real , Polynomial< Dim , Degree , Real > >
		{
			template< unsigned int _Dim , unsigned int _Degree , typename _Real > friend class Polynomial;
			template< unsigned int _Dim , unsigned int Degree1 , unsigned int Degree2 , typename _Real > friend Polynomial< _Dim , Degree1 + Degree2 , _Real > operator * ( const Polynomial< _Dim , Degree1 , _Real > & , const Polynomial< _Dim , Degree2 , _Real > & );
			template< unsigned int _Dim , unsigned int Degree1 , unsigned int Degree2 , typename _Real > friend Polynomial< _Dim , Max< Degree1 , Degree2 >::Value , _Real > operator + ( const Polynomial< _Dim , Degree1 , _Real > & , const Polynomial< _Dim , Degree2 , _Real > & );
			template< unsigned int _Dim , unsigned int Degree1 , unsigned int Degree2 , typename _Real > friend Polynomial< _Dim , Max< Degree1 , Degree2 >::Value , _Real > operator - ( const Polynomial< _Dim , Degree1 , _Real > & , const Polynomial< _Dim , Degree2 , _Real > & );
			template< unsigned int _Dim , unsigned int _Degree , typename _Real > friend std::ostream &operator << ( std::ostream & , const Polynomial< _Dim , _Degree , _Real > & );

			/** The polynomials in Dim-1 dimensions.
			*** The total polynomial is assumed to be _polynomials[0] + _polynomials[1] * (x_Dim) + _polynomials[2] * (x_Dim)^2 + ... */
			Polynomial< Dim-1 , Degree , Real > _polynomials[Degree+1];

			/** This method returns the specified coefficient of the polynomial.*/
			const Real &_coefficient( const unsigned int indices[] , unsigned int maxDegree ) const;

			/** This method returns the specified coefficient of the polynomial.*/
			Real &_coefficient( const unsigned int indices[] , unsigned int maxDegree );

			/** This method evaluates the polynomial at the specified set of coordinates.*/
			Real _evaluate( const Real coordinates[] , unsigned int maxDegree ) const;

			/** This method computes the pull-back of the polynomial given the map from a _Dim-dimensional space to the Dim-dimensional space given by x -> A({x,1}). */	
			template< unsigned int _Dim >
			Polynomial< _Dim , Degree , Real > _pullBack( const Matrix< Real , _Dim+1 , Dim > &A , unsigned int maxDegree ) const;

			/** This method returns true if the polynomial is zero. */
			bool _isZero( unsigned int maxDegree ) const;

			/** This method returns true if the polynomial is a constant. */
			bool _isConstant( unsigned int maxDegree ) const;

			bool _print( std::ostream &ostream , const std::string varNames[] , bool first ) const;
			bool _print( std::ostream &ostream , const std::string varNames[] , std::string suffix , bool first ) const;

			unsigned int _setCoefficients( const Real *coefficients , unsigned int maxDegree );
			unsigned int _getCoefficients( Real *coefficients , unsigned int maxDegree ) const;

			template< unsigned int _D=0 >
			static typename std::enable_if< (_D<Degree ) >::type _SetDegrees( unsigned int coefficientIndex , unsigned int degrees[Dim] );
			template< unsigned int _D=0 >
			static typename std::enable_if< (_D==Degree ) >::type _SetDegrees( unsigned int coefficientIndex , unsigned int degrees[Dim] );

		public:
			/** The number of coefficients in the polynomial. */
			static const unsigned int NumCoefficients = Polynomial< Dim , Degree-1 , Real >::NumCoefficients + Polynomial< Dim-1 , Degree , Real >::NumCoefficients;

			/** The default constructor initializes the coefficients to zero.*/
			Polynomial( void );

			/** This constructor creates a constant polynomial */
			Polynomial( Real c );

			/** This constructor initializes the coefficients. */
			Polynomial( Point< Real , NumCoefficients > coefficients );

			/** The constructor copies over as much of the polynomial as will fit.*/
			template< unsigned int _Degree >
			Polynomial( const Polynomial< Dim , _Degree , Real > &p );

			/** The equality operator copies over as much of the polynomial as will fit.*/
			template< unsigned int _Degree >
			Polynomial &operator= ( const Polynomial< Dim , _Degree , Real > &p );

			/** This method returns the associated coefficient of the polynomial */
			template< typename ... UnsignedInts >
			const Real &coefficient( UnsignedInts ... indices ) const;

			/** This method returns the associated coefficient of the polynomial */
			template< typename ... UnsignedInts >
			Real &coefficient( UnsignedInts ... indices );

			/** This method returns the associated coefficient of the polynomial */
			Real &coefficient( const unsigned int indices[ Dim ] );

			/** This method returns the associated coefficient of the polynomial */
			const Real &coefficient( const unsigned int indices[ Dim ] ) const;

			/** This method returns the associated coefficient of the polynomial */
			Real &coefficient( unsigned int indices[ Dim ] );

			/** This method returns the associated coefficient of the polynomial */
			const Real &coefficient( unsigned int indices[ Dim ] ) const;

			/** This static method sets the degrees of the associated coefficient */
			static void SetDegrees( unsigned int coefficientIndex , unsigned int degrees[Dim] );

			/** This method returns a vector of the coefficients of the polynomial */
			Point< Real , Polynomial< Dim , Degree , Real >::NumCoefficients > coefficients( void ) const;

			/** This method evaluates the polynomial at the prescribed point.*/
			template< typename ... Reals >
			Real operator()( Reals ... coordinates ) const;

			/** This method evaluates the polynomial at the prescribed point.*/
			Real operator()( Point< Real , Dim > p ) const;

			/** This method evaluates the gradient of the polynomial at the prescribed point.*/
			template< typename ... Reals >
			Point< Real , Dim > gradient( Reals ... coordinates ) const;

			/** This method evaluates the gradient of the  polynomial at the prescribed point.*/
			Point< Real , Dim > gradient( Point< Real , Dim > p ) const;

			/** This method evaluates the Hessian of the polynomial at the prescribed point.*/
			template< typename ... Reals >
			SquareMatrix< Real , Dim > hessian( Reals ... coordinates ) const;

			/** This method evaluates the Hessian of the  polynomial at the prescribed point.*/
			SquareMatrix< Real , Dim > hessian( Point< Real , Dim > p ) const;

			/** This method returns the partial derivative with respect to the prescribed dimension.*/
			Polynomial< Dim , (Degree>1) ? Degree-1 : 0 , Real > d( int dim ) const;

			/** This method computes the pull-back of the polynomial given the map from a (_Dim-1)-dimensional space to the Dim-dimensional space given by x -> A({x,1}). */	
			template< unsigned int _Dim >
			Polynomial< _Dim-1 , Degree , Real > operator()( Matrix< Real , _Dim , Dim > A ) const;
			template< unsigned int _Dim >
			Polynomial< _Dim-1 , Degree , Real > pullBack( Matrix< Real , _Dim , Dim > A ) const;
			Polynomial< 1 , Degree , Real > operator()( const Ray< Real , Dim > &ray ) const;



			/** Integrate the polynomial over a unit cube */
			Real integrateUnitCube( void ) const;

			/** Integrate the polynomial over a unit right simplex */
			Real integrateUnitRightSimplex( void ) const;

			/** Returns the matrix taking in the ccoefficients of the polynomial and returning the values at the prescribed points */
			static Matrix< Real , Polynomial< Dim , Degree , Real >::NumCoefficients , Polynomial< Dim , Degree , Real >::NumCoefficients > EvaluationMatrix( const Point< Real , Dim > positions[NumCoefficients] );

			/////////////////////////
			// VectorSpace methods //
			/////////////////////////
			/** This method scales the polynomial */
			void Scale( Real s );

			/** This method adds in the polynomial */
			void Add( const Polynomial & p );
		};

		/** This function returns the product of two polynomials.*/
		template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 , typename Real >
		Polynomial< Dim , Degree1 + Degree2 , Real > operator * ( const Polynomial< Dim , Degree1 , Real > &p1 , const Polynomial< Dim , Degree2 , Real > &p2 );

		/** This function returns the sum of two polynomials. */
		template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 , typename Real >
		Polynomial< Dim , Max< Degree1 , Degree2 >::Value , Real > operator + ( const Polynomial< Dim , Degree1 , Real > &p1 , const Polynomial< Dim , Degree2 , Real > &p2 );

		/** This function returns the difference of two polynomials. */
		template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 , typename Real >
		Polynomial< Dim , Max< Degree1 , Degree2 >::Value , Real > operator - ( const Polynomial< Dim , Degree1 , Real > &p1 , const Polynomial< Dim , Degree2 , Real > &p2 );

		/** This function prints out the polynomial.*/
		template< unsigned int Dim , unsigned int Degree , typename Real >
		std::ostream &operator << ( std::ostream &stream , const Polynomial< Dim , Degree , Real > &poly );


		/** A specialized instance of the Polynomial class in one variable */
		template< unsigned int Degree , typename Real >
		class Polynomial< 0 , Degree , Real > : public VectorSpace< double , Polynomial< 0 , Degree , Real > >
		{
			template< unsigned int _Dim , unsigned int _Degree , typename _Real > friend class Polynomial;
			template< unsigned int Degree1 , unsigned int Degree2 , typename _Real > friend Polynomial< 0 , Degree1 + Degree2 , _Real > operator * ( const Polynomial< 0 , Degree1 , _Real > & , const Polynomial< 0 , Degree2 , _Real > & );
			template< unsigned int Degree1 , unsigned int Degree2 , typename _Real > friend Polynomial< 0 , Max< Degree1 , Degree2 >::Value , _Real > operator + ( const Polynomial< 0 , Degree1 , _Real > & , const Polynomial< 0 , Degree2 , _Real > & );
			template< unsigned int Degree1 , unsigned int Degree2 , typename _Real > friend Polynomial< 0 , Max< Degree1 , Degree2 >::Value , _Real > operator - ( const Polynomial< 0 , Degree1 , _Real > & , const Polynomial< 0 , Degree2 , _Real > & );
			template< unsigned int _Degree , typename _Real > friend std::ostream &operator << ( std::ostream & , const Polynomial< 0 , _Degree , _Real > & );
			template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 , typename _Real > friend Polynomial< Dim , Degree1 + Degree2 , _Real > operator * ( const Polynomial< Dim , Degree1 , _Real > & , const Polynomial< Dim , Degree2 , _Real > & );
			template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 , typename _Real > friend Polynomial< Dim , Max< Degree1 , Degree2 >::Value , _Real > operator + ( const Polynomial< Dim , Degree1 , _Real > & , const Polynomial< Dim , Degree2 , _Real > & );
			template< unsigned int Dim , unsigned int Degree1 , unsigned int Degree2 , typename _Real > friend Polynomial< Dim , Max< Degree1 , Degree2 >::Value , _Real > operator - ( const Polynomial< Dim , Degree1 , _Real > & , const Polynomial< Dim , Degree2 , _Real > & );
			template< unsigned int Dim , unsigned int _Degree , typename _Real > friend std::ostream &operator << ( std::ostream & , const Polynomial< Dim , _Degree , _Real > & );
		protected:
			/** The coefficients of the polynomial. */
			Real _coefficients[1];

			/** This method returns the specified coefficient of the polynomial.*/
			const Real &_coefficient( const unsigned int indices[] , unsigned int maxDegree ) const;

			/** This method returns the specified coefficient of the polynomial.*/
			Real &_coefficient( const unsigned int indices[] , unsigned int maxDegree );

			/** This method evaluates the polynomial at the specified set of coordinates.*/
			Real _evaluate( const Real coordinates[] , unsigned int maxDegree ) const;

			/** This method computes the pull-back of the polynomial given the map from a _Dim-dimensional space to the 1-dimensional space given by x -> A({x,1}). */	
			template< unsigned int _Dim >
			Polynomial< _Dim , Degree , Real > _pullBack( const Matrix< Real , _Dim+1 , 0 > &A , unsigned int maxDegree ) const;

			/** This method returns true if the polynomial is zero. */
			bool _isZero( unsigned int maxDegree ) const;

			/** This method returns true if the polynomial is a constant. */
			bool _isConstant( unsigned int maxDegree ) const;

			bool _print( std::ostream &ostream , const std::string varNames[] , bool first ) const;
			bool _print( std::ostream &ostream , const std::string varNames[] , std::string suffix , bool first ) const;

			unsigned int _setCoefficients( const Real *coefficients , unsigned int maxDegree );
			unsigned int _getCoefficients( Real *coefficients , unsigned int maxDegree ) const;
		public:
			/** The number of coefficients in the polynomial. */
			static const unsigned int NumCoefficients = 1;

			/** The default constructor initializes the coefficients to zero.*/
			Polynomial( void );

			/** This constructor creates a constant polynomial */
			Polynomial( Real c );

			/** This constructor initializes the coefficients (starting with lower degrees). */
			Polynomial( Point< Real , NumCoefficients > coefficients );

			/** The constructor copies over as much of the polynomial as will fit.*/
			template< unsigned int _Degree >
			Polynomial( const Polynomial< 0 , _Degree , Real > &p );

			/** The equality operator copies over as much of the polynomial as will fit.*/
			template< unsigned int _Degree >
			Polynomial &operator= ( const Polynomial< 0 , _Degree , Real > &p );

			/** This method returns the coefficient.*/
			const Real &coefficient( void ) const;

			/** This method returns the coefficient.*/
			Real &coefficient( void );

			/** This method returns the associated coefficient of the polynomial */
			Real &coefficient( const unsigned int indices[0] );

			/** This method returns the associated coefficient of the polynomial */
			const Real &coefficient( const unsigned int indices[0] ) const;

			/** This method returns the associated coefficient of the polynomial */
			Real &coefficient( unsigned int indices[0] );

			/** This method returns the associated coefficient of the polynomial */
			const Real &coefficient( unsigned int indices[0] ) const;

			/** This static method sets the degrees of the associated coefficient */
			static void SetDegrees( unsigned int coefficientIndex , unsigned int degrees[0] );

			/** This method returns a vector of the coefficients of the polynomial */
			Point< Real , Polynomial< 0 , Degree , Real >::NumCoefficients > coefficients( void ) const;

			/** This method evaluates the polynomial at a given value.*/
			Real operator()( void ) const;

			/** This method evaluates the polynomial at the prescribed point.*/
			Real operator()( Point< Real , 0 > p ) const;

			/** This method returns the derivative of the polynomial.*/
			Polynomial< 0 , (Degree>1) ? Degree-1 : 0 , Real > d( unsigned int d=0 ) const;

			/** This method computes the pull-back of the polynomial given the map from a (_Dim-1)-dimensional space to the 1-dimensional space given by x -> A({x,1}). */	
			template< unsigned int _Dim >
			Polynomial< _Dim-1 , Degree , Real > operator()( const Matrix< Real , _Dim , 0 > &A ) const;

			/** Integrate the polynomial over a unit interval */
			Real integrateUnitCube( void ) const;

			/** Integrate the polynomial over a unit interval */
			Real integrateUnitRightSimplex( void ) const;

			/** Returns the matrix taking in the ccoefficients of the polynomial and returning the values at the prescribed points */
			static Matrix< Real , Polynomial< 0 , Degree , Real >::NumCoefficients  , Polynomial< 0 , Degree , Real >::NumCoefficients> EvaluationMatrix( const Point< Real , 0 > positions[NumCoefficients] );

			/////////////////////////
			// VectorSpace methods //
			/////////////////////////
			/** This method scales the polynomial*/
			void Scale( Real s );

			/** This method adds in the polynomial */
			void Add( const Polynomial & p );
		};

		/** This function returns the product of two polynomials.*/
		template< unsigned int Degree1 , unsigned int Degree2 , typename Real >
		Polynomial< 0 , Degree1 + Degree2 , Real > operator * ( const Polynomial< 0 , Degree1 , Real > &p1 , const Polynomial< 0 , Degree2 , Real > &p2 );

		/** This function returns the sum of two polynomials. */
		template< unsigned int Degree1 , unsigned int Degree2 , typename Real >
		Polynomial< 0 , Max< Degree1 , Degree2 >::Value , Real > operator + ( const Polynomial< 0 , Degree1 , Real > &p1 , const Polynomial< 0 , Degree2 , Real > &p2 );

		/** This function returns the difference of two polynomials. */
		template< unsigned int Degree1 , unsigned int Degree2 , typename Real >
		Polynomial< 0 , Max< Degree1 , Degree2 >::Value , Real > operator - ( const Polynomial< 0 , Degree1 , Real > &p1 , const Polynomial< 0 , Degree2 , Real > &p2 );

		/** A specialization that allows us for the recursive computation of the number of coefficints, that allows us to avoid special-casing when the the Degree is zero. */
		template< unsigned int Dim , typename Real >
		class Polynomial< Dim , (unsigned int)-1 , Real >
		{
		public:
			static const unsigned int NumCoefficients = 0;
		};


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


		/** Sets the roots of the 1D polynomial and returns the number of roots set.
		*** The method is only specialized for degrees 1, 2, 3, and 4. */
		template< unsigned int Degree , typename Real >
		unsigned int Roots( const Polynomial< 1 , Degree , Real > &p , Real *r , double eps=1e-14 );

#include "Polynomial.inl"
	}
}
#endif // POLYNOMIAL_INCLUDED
