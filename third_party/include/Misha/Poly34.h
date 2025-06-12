// poly34.h : solution of cubic and quartic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com

#ifndef POLY34_INCLUDED
#define POLY34_INCLUDED

#include <math.h>

namespace SergeyKhashin
{
	namespace Poly34
	{
		// x - array of size 1
		// return 1: 1 real roots x[0]
		static unsigned int SolveP1( double *x , double a , double eps=1e-14 );												// solve equation x + a = 0

		// x - array of size 2
		// return 2: 2 real roots x[0], x[1]
		// return 0: pair of complex roots: x[0]켲*x[1]
		static unsigned int SolveP2( double *x , double a , double b , double eps=1e-14 );									// solve equation x^2 + a*x + b = 0

		// x - array of size 3
		// return 3: 3 real roots x[0], x[1], x[2]
		// return 1: 1 real root x[0] and pair of complex roots: x[1]켲*x[2]
		static unsigned int SolveP3( double *x , double a , double b , double c , double eps=1e-14 );						// solve cubic equation x^3 + a*x^2 + b*x + c = 0

		// x - array of size 4
		// return 4: 4 real roots x[0], x[1], x[2], x[3], possible multiple roots
		// return 2: 2 real roots x[0], x[1] and complex x[2]켲*x[3], 
		// return 0: two pair of complex roots: x[0]켲*x[1],  x[2]켲*x[3], 
		static unsigned int SolveP4( double *x , double a , double b , double c , double d , double eps=1e-14 );			// solve equation x^4 + a*x^3 + b*x^2 + c*x + d = 0 by Dekart-Euler method

		// x - array of size 5
		// return 5: 5 real roots x[0], x[1], x[2], x[3], x[4], possible multiple roots
		// return 3: 3 real roots x[0], x[1], x[2] and complex x[3]켲*x[4], 
		// return 1: 1 real root x[0] and two pair of complex roots: x[1]켲*x[2],  x[3]켲*x[4], 
		static unsigned int SolveP5( double *x , double a , double b , double c , double d , double e , double eps=1e-14 );	// solve equation x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0

#include "Poly34.inl"
	}
}
#endif // POLY34_INCLUDED