/*
Copyright (c) 2018, Michael Kazhdan
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
//#include <GL/glew.h>
//#include <GL/glut.h>

namespace MishaK
{
	class Camera
	{
		void _setRight( void )
		{
			right = Point3D< double >::CrossProduct( forward , up );
			right /= Point3D< double >::Length( right );
		}
		void rotatePoint( SquareMatrix< double , 3 > R , Point3D< double > center )
		{
			position = R * ( position - center ) + center;

			forward = R * forward;
			right = R * right;
			up = R * up;

			_setRight();
		}

		void rotatePoint( Point3D< double > axis , double angle , Point3D< double > center )
		{
			Point3D< double > p , r , f , u;
			Point3D< double > v[3];
			double c , s;
			double d[3];

			v[2] = axis/Point3D< double >::Length( axis );

			v[0] = Point3D< double >::CrossProduct( v[2] , Point3D< double >( 1 , 0 , 0 ) );
			if( Point3D< double >::SquareNorm( v[0] )<.001) v[0] = Point3D< double >::CrossProduct( v[2] , Point3D< double >( 0 , 1 , 0 ) );
			v[0] /= Point3D< double >::Length( v[0] );
			v[1] = Point3D< double >::CrossProduct( v[2] , v[0] );
			v[1] /= Point3D< double >::Length( v[1] );

			c = cos(angle);
			s = sin(angle);

			p = position-center;
			for( int j=0 ; j<3 ; j++ ) d[j] = Point3D< double >::Dot( p , v[j] );

			position = v[2]*d[2] + v[0]*(d[0]*c+d[1]*s) + v[1]*(-d[0]*s+d[1]*c) + center;

			for( int j=0 ; j<3 ; j++ )
			{
				r[j] = Point3D< double >::Dot(   right , v[j] );
				f[j] = Point3D< double >::Dot( forward , v[j] );
				u[j] = Point3D< double >::Dot(      up , v[j] );
			}

			r = v[2]*r[2]+v[0]*(r[0]*c+r[1]*s)+v[1]*(-r[0]*s+r[1]*c);
			f = v[2]*f[2]+v[0]*(f[0]*c+f[1]*s)+v[1]*(-f[0]*s+f[1]*c);
			u = v[2]*u[2]+v[0]*(u[0]*c+u[1]*s)+v[1]*(-u[0]*s+u[1]*c);

			forward	= f / Point3D< double >::Length(f);
			right	= r / Point3D< double >::Length(r);
			up		= u / Point3D< double >::Length(u);

			_setRight();
		}

	public:
		Point3D< double > position , forward , up , right;

		Camera( void )
		{
			position = Point3D< double >( 0 , 0 , 0 );
			forward  = Point3D< double >( 0 , 0 , -1 );
			up       = Point3D< double >( 0 , 1 , 0 );
			_setRight();
		}
		Camera( Point3D< double > p , Point3D< double > f , Point3D< double > u )
		{
			position = p , forward = f , up = u;
			_setRight();
		}
		void draw( void ) const
		{
			glMatrixMode( GL_MODELVIEW );        
			glLoadIdentity();
			gluLookAt(
				position[0] , position[1] , position[2] ,
				position[0]+forward[0] , position[1]+forward[1] , position[2]+forward[2] ,
				up[0] , up[1] , up[2]
			);
		}

		void translate( Point3D< double > t ){ position += t; }
		void rotate       ( SquareMatrix< double , 3  > R , Point3D< double > p = Point3D< double >() ){ rotatePoint( R , p ); }
		void rotate       ( Point3D< double > axis , double angle , Point3D< double > p = Point3D< double >() ){ rotatePoint( axis , angle , p ); }
		void rotateUp     ( double angle , Point3D< double > p=Point3D< double >() ){ rotatePoint( up      , angle , p ); }
		void rotateRight  ( double angle , Point3D< double > p=Point3D< double >() ){ rotatePoint( right   , angle , p ); }
		void rotateForward( double angle , Point3D< double > p=Point3D< double >() ){ rotatePoint( forward , angle , p ); }

		Point2D< double > project( Point3D< double > p , bool orthographic )
		{
			p -= position;
			double x = Point3D< double >::Dot( p , right ) , y = Point3D< double >::Dot( p , up ) , z = Point3D< double >::Dot( p , forward );
			if( orthographic ) return Point2D< double >( x , y );
			else               return Point2D< double >( x/z , y/1 );
		}

		bool read( FILE* fp )
		{
			Point3D< float > temp;
			if( fscanf( fp , " %f %f %f " , &(temp[0]) , &(temp[1]) , &(temp[2]) )!=3 ) return false;
#if 1
			for( unsigned int d=0 ; d<3 ; d++ ) position[d] = temp[d];
#else
			else position = Point3D< double >( temp );
#endif
			if( fscanf( fp , " %f %f %f " , &(temp[0]) , &(temp[1]) , &(temp[2]) )!=3 ) return false;
#if 1
			for( unsigned int d=0 ; d<3 ; d++ ) forward[d] = temp[d];
#else
			else forward = Point3D< double >( temp );
#endif
			if( fscanf( fp , " %f %f %f " , &(temp[0]) , &(temp[1]) , &(temp[2]) )!=3 ) return false;
#if 1
			for( unsigned int d=0 ; d<3 ; d++ ) up[d] = temp[d];
#else
			else up = Point3D< double >( temp );
#endif
			if( fscanf( fp , " %f %f %f " , &(temp[0]) , &(temp[1]) , &(temp[2]) )!=3 ) return false;
#if 1
			for( unsigned int d=0 ; d<3 ; d++ ) right[d] = temp[d];
#else
			else right = Point3D< double >( temp );
#endif
			return true;
		}
		bool write( FILE* fp ) const
		{
			fprintf( fp , "%f %f %f\n" , position[0] , position[1] , position[2] );
			fprintf( fp , "%f %f %f\n" ,  forward[0] ,  forward[1] ,  forward[2] );
			fprintf( fp , "%f %f %f\n" ,       up[0] ,       up[1] ,       up[2] );
			fprintf( fp , "%f %f %f\n" ,    right[0] ,    right[1] ,    right[2] );
			return true;
		}
	};
}
#endif // CAMERA_INCLUDED