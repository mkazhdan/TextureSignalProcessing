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
#ifndef HSV_INCLUDED
#define HSV_INCLUDED


#include <Misha/Geometry.h>

namespace MishaK
{
	// h [0,1) ,s [0,1], v [0,1]
	// Hue is mapped to [0,240) deg to avoid circularity
	template< typename Real >
	Point3D< Real > HSV2RGB( Real h,  Real s , Real v )
	{
		const Real c = s*v;
		// use next line for full [0,360) hue
		const Real _h = h * (Real)6.0;
		// use next line for partial [0,240) hue
		//const Real _h = h * 4.0;
		const Real _h_mod_2 = (Real)( _h - floor( _h/2 ) * 2 );
		Real x = (Real)( c*(1 - fabs( _h_mod_2 - 1 ) ) );
		Real r, g, b;
		r = g = b = 0.0;
		if (_h <= 1){
			r = c;
			g = x;
		}
		else if (_h <= 2){
			r = x;
			g = c;
		}
		else if (_h <= 3){
			g = c;
			b = x;
		}
		else if (_h <= 4){
			g = x;
			b = c;
		}
		else if (_h <= 5){
			r = x;
			b = c;
		}
		else{
			r = c;
			b = x;
		}
		Real m = v - c;
		return Point3D< Real >( r+m , g+m , b+m );
	}
}
#endif //HSV_INCLUDED
