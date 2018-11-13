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
#ifndef CAMERA_INCLUDED
#define CAMERA_INCLUDED
#include <Misha/Geometry.h>

template< typename Real >
class Camera
{
	void _setRight( void )
	{
		right = Point3D< Real >::CrossProduct( forward , up );
		right /= Length( right );
	}
	void rotatePoint( Point3D< Real > axis , Real angle , Point3D< Real > center )
	{
		Point3D< Real > p , r , f , u;
		Point3D< Real > v[3];
		Real c , s;
		Real d[3];

		v[2] = axis/Length( axis );

		v[0] = Point3D< Real >::CrossProduct( v[2] , Point3D< Real >( 1 , 0 , 0 ) );
		if( Point3D< Real >::SquareNorm( v[0] )<.001 ) v[0] = Point3D< Real >::CrossProduct( v[2] , Point3D< Real >( 0 , 1 , 0 ) );
		v[0] /= Length( v[0] );
		v[1] = Point3D< Real >::CrossProduct( v[2] , v[0] );
		v[1] /= Length( v[1] );

		c = cos(angle);
		s = sin(angle);

		p = position-center;
		for( int j=0 ; j<3 ; j++ ) d[j] = Point3D< Real >::Dot( p , v[j] );

		position = v[2]*d[2] + v[0]*(d[0]*c+d[1]*s) + v[1]*(-d[0]*s+d[1]*c) + center;

		for( int j=0 ; j<3 ; j++ )
		{
			r[j] = Point3D< Real >::Dot(   right , v[j] );
			f[j] = Point3D< Real >::Dot( forward , v[j] );
			u[j] = Point3D< Real >::Dot(      up , v[j] );
		}

		r = v[2]*r[2]+v[0]*(r[0]*c+r[1]*s)+v[1]*(-r[0]*s+r[1]*c);
		f = v[2]*f[2]+v[0]*(f[0]*c+f[1]*s)+v[1]*(-f[0]*s+f[1]*c);
		u = v[2]*u[2]+v[0]*(u[0]*c+u[1]*s)+v[1]*(-u[0]*s+u[1]*c);

		forward	= f / Length(f);
		right	= r / Length(r);
		up		= u / Length(u);

		_setRight();
	}

public:
	Point3D< Real > position , forward , up , right;

	Camera( void )
	{
		position = Point3D< Real >( 0 , 0 , 0 );
		forward  = Point3D< Real >( 0 , 0 , 1 );
		up       = Point3D< Real >( 0 , 1 , 0 );
		_setRight();
	}
	Camera( Point3D< Real > p , Point3D< Real > f , Point3D< Real > u )
	{
		position = p , forward = f , up = u;
		_setRight();
	}
	void draw( void )
	{
		glMatrixMode( GL_MODELVIEW );        
		glLoadIdentity();
		gluLookAt(
			position[0] , position[1] , position[2] ,
			position[0]+forward[0] , position[1]+forward[1] , position[2]+forward[2] ,
			up[0] , up[1] , up[2]
		);
	}

	void translate( Point3D< Real > t ){ position += t; }
	void rotateUp     ( Real angle , Point3D< Real > p=Point3D< Real >() ){ rotatePoint( up      , angle , p ); }
	void rotateRight  ( Real angle , Point3D< Real > p=Point3D< Real >() ){ rotatePoint( right   , angle , p ); }
	void rotateForward( Real angle , Point3D< Real > p=Point3D< Real >() ){ rotatePoint( forward , angle , p ); }

	Point2D< Real > project( Point3D< Real > p , bool orthographic )
	{
		p -= position;
		Real x = Point3D< Real >::Dot( p , right ) , y = Point3D< Real >::Dot( p , up ) , z = Point3D< Real >::Dot( p , forward );
		if( orthographic ) return Point2D< Real >( x , y );
		else               return Point2D< Real >( x/z , y/1 );
	}

	void write( FILE *fp ) const
	{
		fwrite( &position , sizeof( Point3D< Real > ) , 1 , fp );
		fwrite( &forward  , sizeof( Point3D< Real > ) , 1 , fp );
		fwrite( &right    , sizeof( Point3D< Real > ) , 1 , fp );
		fwrite( &up       , sizeof( Point3D< Real > ) , 1 , fp );
	}
	void read( FILE *fp )
	{
		fread( &position , sizeof( Point3D< Real > ) , 1 , fp );
		fread( &forward  , sizeof( Point3D< Real > ) , 1 , fp );
		fread( &right    , sizeof( Point3D< Real > ) , 1 , fp );
		fread( &up       , sizeof( Point3D< Real > ) , 1 , fp );
	}
};
#endif // CAMERA_INCLUDED
