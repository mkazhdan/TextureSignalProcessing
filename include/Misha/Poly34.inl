// poly34.cpp : solution of cubic and quartic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com
// Thanks to Alexandr Rakhmanin <rakhmanin (at) gmail.com>
// public domain
//

#define MARC_ALEXA_FIX

//=============================================================================
// _root3, root3 from http://prografix.narod.ru
//=============================================================================
inline double _root3 ( double x )
{
	double s = 1.;
	while ( x < 1. )
	{
		x *= 8.;
		s *= 0.5;
	}
	while ( x > 8. )
	{
		x *= 0.125;
		s *= 2.;
	}
	double r = 1.5;
	r -= 1./3. * ( r - x / ( r * r ) );
	r -= 1./3. * ( r - x / ( r * r ) );
	r -= 1./3. * ( r - x / ( r * r ) );
	r -= 1./3. * ( r - x / ( r * r ) );
	r -= 1./3. * ( r - x / ( r * r ) );
	r -= 1./3. * ( r - x / ( r * r ) );
	return r * s;
}


inline double root3 ( double x )
{
	if      ( x > 0 ) return  _root3 ( x );
	else if ( x < 0 ) return -_root3 (-x );
	else              return 0.;
}

// x - array of size 1
// return 1: 1 real root x
// solve equation x + a = 0
inline unsigned int SolveP1( double *x , double a , double eps )
{
	x[0] = -a;
	return 1;
}


// x - array of size 2
// return 2: 2 real roots x[0], x[1]
// return 0: pair of complex roots: x[0]켲*x[1]
// solve equation x^2 + a*x + b = 0
inline unsigned int SolveP2( double *x , double a , double b , double eps )
{
	double D = 0.25*a*a - b;
	if( D>=0 )
	{
		D = sqrt(D);
		x[0] = -0.5*a + D;
		x[1] = -0.5*a - D;
		return 2;
	}
	x[0] = -0.5*a;
	x[1] = sqrt(-D);
	return 0;
}

// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] � i*x[2], return 1
// solve cubic equation x^3 + a*x^2 + b*x + c = 0
inline unsigned int SolveP3( double *x , double a , double b , double c , double eps )
{
	//		const double eps = 1e-14;
	const double TwoPi = 6.28318530717958648;

	double a2 = a*a;
	double q  = (a2 - 3*b)/9; 
	double r  = (a*(2*a2-9*b) + 27*c)/54;
	// equation x^3 + q*x + r = 0
	double r2 = r*r;
	double q3 = q*q*q;
	double A , B;
#ifdef MARC_ALEXA_FIX
	// Computation should be stable in the "if" part of the branch without the epsilon,
	// but this ensures that "sqrt(r2-q3)" in the "else" branch does not lead to a NaN.
	if (r2 <= q3) {//<<-- FIXED!
#else // !MARC_ALEXA_FIX
	if (r2 <= (q3 + eps)) {//<<-- FIXED!
#endif // MARC_ALEXA_FIX
		double t=r/sqrt(q3);
		if( t<-1) t=-1;
		if( t> 1) t= 1;
		t=acos(t);
#ifdef MARC_ALEXA_FIX
#pragma message( "[WARNING] Does this require a safeguard?" )
#endif // MARC_ALEXA_FIX
		a/=3; q=-2*sqrt(q);
		x[0]=q*cos(t/3)-a;
		x[1]=q*cos((t+TwoPi)/3)-a;
		x[2]=q*cos((t-TwoPi)/3)-a;
		return(3);
	} else {
		//A =-pow(fabs(r)+sqrt(r2-q3),1./3); 
		A =-root3(fabs(r)+sqrt(r2-q3)); 
		if( r<0 ) A=-A;
		B = (A==0? 0 : B=q/A);

		a/=3;
		x[0] =(A+B)-a;
		x[1] = -0.5*(A+B)-a;
		x[2] =  0.5*sqrt(3.)*(A-B);
		if(fabs(x[2])<eps) { x[2]=x[1]; return(2); }
		return(1);
	}
	}

// a>=0!
// returns:  a+i*s = sqrt(x+i*y)
inline void CSqrt( double x , double y , double &a , double &b )
{
	double r  = sqrt(x*x+y*y);
	if( y==0 ) { 
		r = sqrt(r);
		if(x>=0) { a=r; b=0; } else { a=0; b=r; }
	} else {		// y != 0
		a = sqrt(0.5*(x+r));
		b = 0.5*y/a;
	}
}

// solve equation x^4 + b*x^2 + d = 0
inline unsigned int SolveP4Bi( double *x , double b , double d )
{
	double D = b*b-4*d;
	if( D>=0 ) 
	{
		double sD = sqrt(D);
		double x1 = (-b+sD)/2;
		double x2 = (-b-sD)/2;	// x2 <= x1
		if( x2>=0 )				// 0 <= x2 <= x1, 4 real roots
		{
			double sx1 = sqrt(x1);
			double sx2 = sqrt(x2);
			x[0] = -sx1;
			x[1] =  sx1;
			x[2] = -sx2;
			x[3] =  sx2;
			return 4;
		}
		if( x1 < 0 )				// x2 <= x1 < 0, two pair of imaginary roots
		{
			double sx1 = sqrt(-x1);
			double sx2 = sqrt(-x2);
			x[0] =    0;
			x[1] =  sx1;
			x[2] =    0;
			x[3] =  sx2;
			return 0;
		}
		// now x2 < 0 <= x1 , two real roots and one pair of imginary root
		double sx1 = sqrt( x1);
		double sx2 = sqrt(-x2);
		x[0] = -sx1;
		x[1] =  sx1;
		x[2] =    0;
		x[3] =  sx2;
		return 2;
	} else { // if( D < 0 ), two pair of compex roots
		double sD2 = 0.5*sqrt(-D);
		CSqrt( -0.5*b ,  sD2, x[0] , x[1] );
		CSqrt( -0.5*b , -sD2, x[2] , x[3] );
		return 0;
	} // if( D>=0 ) 
}

inline void SWAP( double &a , double &b )
{
	double t = a;
	a = b;
	b = t;
}

// make: a <= b <= c
inline void dblSort3( double &a , double &b , double &c )
{
	if( a>b ) SWAP(a,b);	// now a<=b
	if( c<b )
	{
		SWAP(b,c);			// now a<=b, b<=c
		if( a>b ) SWAP(a,b);// now a<=b
	}
}

// solve equation x^4 + b*x^2 + c*x + d
inline unsigned int SolveP4De( double *x , double b , double c , double d , double eps )
{
	//if( c==0 ) return SolveP4Bi(x,b,d); // After that, c!=0
	if( fabs(c)<eps*(fabs(b)+fabs(d)) ) return SolveP4Bi(x,b,d); // After that, c!=0

	int res3 = SolveP3( x , 2*b , b*b-4*d , -c*c , eps );	// solve resolvent
	// by Viet theorem:  x1*x2*x3=-c*c not equals to 0, so x1!=0, x2!=0, x3!=0
	if( res3>1 )	// 3 real roots, 
	{				
		dblSort3(x[0], x[1], x[2]);	// sort roots to x[0] <= x[1] <= x[2]
		// Note: x[0]*x[1]*x[2]= c*c > 0
		if( x[0] > 0) // all roots are positive
		{
			double sz1 = sqrt(x[0]);
			double sz2 = sqrt(x[1]);
			double sz3 = sqrt(x[2]);
			// Note: sz1*sz2*sz3= -c (and not equal to 0)
			if( c>0 )
			{
				x[0] = (-sz1 -sz2 -sz3)/2;
				x[1] = (-sz1 +sz2 +sz3)/2;
				x[2] = (+sz1 -sz2 +sz3)/2;
				x[3] = (+sz1 +sz2 -sz3)/2;
				return 4;
			}
			// now: c<0
			x[0] = (-sz1 -sz2 +sz3)/2;
			x[1] = (-sz1 +sz2 -sz3)/2;
			x[2] = (+sz1 -sz2 -sz3)/2;
			x[3] = (+sz1 +sz2 +sz3)/2;
			return 4;
		} // if( x[0] > 0) // all roots are positive
		  // now x[0] <= x[1] < 0, x[2] > 0
		  // two pair of comlex roots
		double sz1 = sqrt(-x[0]);
		double sz2 = sqrt(-x[1]);
		double sz3 = sqrt( x[2]);

		if( c>0 )	// sign = -1
		{
			x[0] = -sz3/2;					
			x[1] = ( sz1 -sz2)/2;		// x[0]켲*x[1]
			x[2] =  sz3/2;
			x[3] = (-sz1 -sz2)/2;		// x[2]켲*x[3]
			return 0;
		}
		// now: c<0 , sign = +1
		x[0] =   sz3/2;
		x[1] = (-sz1 +sz2)/2;
		x[2] =  -sz3/2;
		x[3] = ( sz1 +sz2)/2;
		return 0;
	} // if( res3>1 )	// 3 real roots, 
	  // now resoventa have 1 real and pair of compex roots
	  // x[0] - real root, and x[0]>0, 
	  // x[1]켲*x[2] - complex roots, 
	  // x[0] must be >=0. But one times x[0]=~ 1e-17, so:
	if (x[0] < 0) x[0] = 0;
	double sz1 = sqrt(x[0]);
	double szr, szi;
	CSqrt(x[1], x[2], szr, szi);  // (szr+i*szi)^2 = x[1]+i*x[2]
	if( c>0 )	// sign = -1
	{
		x[0] = -sz1/2-szr;			// 1st real root
		x[1] = -sz1/2+szr;			// 2nd real root
		x[2] = sz1/2; 
		x[3] = szi;
		return 2;
	}
	// now: c<0 , sign = +1
	x[0] = sz1/2-szr;			// 1st real root
	x[1] = sz1/2+szr;			// 2nd real root
	x[2] = -sz1/2;
	x[3] = szi;
	return 2;
}

// one Newton step for x^4 + a*x^3 + b*x^2 + c*x + d
inline double N4Step( double x , double a , double b , double c , double d )
{
	double fxs= ((4*x+3*a)*x+2*b)*x+c;	// f'(x)
	if (fxs == 0) return x;	//return 1e99; <<-- FIXED!
	double fx = (((x+a)*x+b)*x+c)*x+d;	// f(x)
	return x - fx/fxs;
} 
//-----------------------------------------------------------------------------
// x - array of size 4
// return 4: 4 real roots x[0], x[1], x[2], x[3], possible multiple roots
// return 2: 2 real roots x[0], x[1] and complex x[2]켲*x[3], 
// return 0: two pair of complex roots: x[0]켲*x[1],  x[2]켲*x[3], 
inline unsigned int SolveP4( double * x , double a , double b , double c , double d , double eps ) {	// solve equation x^4 + a*x^3 + b*x^2 + c*x + d by Dekart-Euler method
	// move to a=0:
	double d1 = d + 0.25*a*( 0.25*b*a - 3./64*a*a*a - c);
	double c1 = c + 0.5*a*(0.25*a*a - b);
	double b1 = b - 0.375*a*a;
	int res = SolveP4De( x , b1 , c1 , d1 , eps );
	if( res==4) { x[0]-= a/4; x[1]-= a/4; x[2]-= a/4; x[3]-= a/4; }
	else if (res==2) { x[0]-= a/4; x[1]-= a/4; x[2]-= a/4; }
	else             { x[0]-= a/4; x[2]-= a/4; }
	// one Newton step for each real root:
	if( res>0 )
	{
		x[0] = N4Step(x[0], a,b,c,d);
		x[1] = N4Step(x[1], a,b,c,d);
	}
	if( res>2 )
	{
		x[2] = N4Step(x[2], a,b,c,d);
		x[3] = N4Step(x[3], a,b,c,d);
	}
	return res;
}


inline double F5( double t , double a , double b , double c , double d , double e )
{
	return (((((t+a)*t+b)*t+c)*t+d)*t+e);
}

// return real root of x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
inline double SolveP5_1( double a , double b , double c , double d , double e , double eps )
{
	//		const double eps = 1e-14;

	int cnt;
	if( fabs(e)<eps ) return 0;

	double brd =  fabs(a);			// brd - border of real roots
	if( fabs(b)>brd ) brd = fabs(b);
	if( fabs(c)>brd ) brd = fabs(c);
	if( fabs(d)>brd ) brd = fabs(d);
	if( fabs(e)>brd ) brd = fabs(e);
	brd++;							// brd - border of real roots

	double x0, f0;					// less than root
	double x1, f1;					// greater than root
	double x2, f2, f2s;				// next values, f(x2), f'(x2)
	double dx;

	if( e<0 ) { x0 =   0; x1 = brd; f0=e; f1=F5(x1,a,b,c,d,e); x2 = 0.01*brd; }	// positive root
	else	  { x0 =-brd; x1 =   0; f0=F5(x0,a,b,c,d,e); f1=e; x2 =-0.01*brd; }	// negative root

	if( fabs(f0)<eps ) return x0;
	if( fabs(f1)<eps ) return x1;

	// now x0<x1, f(x0)<0, f(x1)>0
	// Firstly 10 bisections
	for( cnt=0; cnt<10; cnt++)
	{
		x2 = (x0 + x1) / 2;					// next point
		//x2 = x0 - f0*(x1 - x0) / (f1 - f0);		// next point
		f2 = F5(x2,a,b,c,d,e);				// f(x2)
		if( fabs(f2)<eps ) return x2;
		if( f2>0 ) { x1=x2; f1=f2; }
		else       { x0=x2; f0=f2; }
	}

	// At each step:
	// x0<x1, f(x0)<0, f(x1)>0.
	// x2 - next value
	// we hope that x0 < x2 < x1, but not necessarily
	do {
		if(cnt++>50) break;
		if( x2<=x0 || x2>= x1 ) x2 = (x0 + x1)/2;	// now  x0 < x2 < x1
		f2 = F5(x2,a,b,c,d,e);								// f(x2)
		if( fabs(f2)<eps ) return x2;
		if( f2>0 ) { x1=x2; f1=f2; }
		else       { x0=x2; f0=f2; }
		f2s= (((5*x2+4*a)*x2+3*b)*x2+2*c)*x2+d;		// f'(x2)
		if( fabs(f2s)<eps ) { x2=1e99; continue; }
		dx = f2/f2s;
		x2 -= dx;
	} while(fabs(dx)>eps);
	return x2;
}

// solve equation x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
inline unsigned int SolveP5( double *x , double a , double b , double c , double d , double e , double eps )
{
	double r = x[0] = SolveP5_1( a , b , c , d , e , eps );
	double a1 = a+r, b1=b+r*a1, c1=c+r*b1, d1=d+r*c1;
	return 1+SolveP4( x+1 , a1 , b1 , c1 , d1 , eps );
}
