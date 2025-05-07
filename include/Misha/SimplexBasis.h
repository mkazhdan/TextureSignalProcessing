/*
Copyright (c) 2022, Michael Kazhdan
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

#ifndef SIMPLEX_BASIS_INCLUDED
#define SIMPLEX_BASIS_INCLUDED

#include <map>
#include <set>
#include "Exceptions.h"
#include "Polynomial.h"
#include "Geometry.h"

namespace MishaK
{
	template< unsigned int Size , typename Index , bool SmallestFirst >
	struct MultiIndex
	{
		MultiIndex( void ){ for( unsigned int i=0 ; i<Size ; i++ ) _indices[i] = 0; }
		MultiIndex( const Index indices[] ){ _init( indices ); }
		MultiIndex(       Index indices[] ){ _init( indices ); }
		template< typename ... Is > MultiIndex( Is ... indices );
		bool operator < ( const MultiIndex &idx ) const;
		bool operator == ( const MultiIndex &idx ) const;
		const Index &operator[] ( unsigned int idx ) const { return _indices[idx]; }
		typedef std::          map< MultiIndex , unsigned int          > map;
		typedef std::          set< MultiIndex                         > set;

	protected:
		void _init( const Index indices[] );
		Index _indices[ Size ];
	};

	template< unsigned int Size , typename Index , bool SmallestFirst >
	std::ostream &operator << ( std::ostream &os , const MultiIndex< Size , Index , SmallestFirst > &idx );

	template< unsigned int Dim > struct RightSimplex;

	template<>
	struct RightSimplex< 0 >
	{
		static const unsigned int Dim = 0;

		template< unsigned int Degree , typename Real > using PFunction = Polynomial::Polynomial< Dim , Degree , Real >;

		static const Point< double , Dim > &Vertex( unsigned int d );

		template< unsigned int EmbeddingDimension >
		static SquareMatrix< double , Dim > Metric( const Simplex< double , EmbeddingDimension , Dim > &s );

		template< unsigned int EmbeddingDimension , typename PositionFunctor >
		static SquareMatrix< double , Dim > Metric( PositionFunctor positionFunctor );

		template< unsigned int Degree , typename Real >
		static double Integral( const PFunction< Degree , Real > &P , SquareMatrix< double , Dim > g = SquareMatrix< double , Dim >::Identity() );
	};

	template< unsigned int Dim >
	struct RightSimplex
	{
		template< unsigned int Degree , typename Real > using PFunction = Polynomial::Polynomial< Dim , Degree , Real >;
		template< unsigned int Degree , typename Real > using PVectorField = Point< Polynomial::Polynomial< Dim , Degree , Real > , Dim >;
		template< unsigned int Degree , typename Real > using PMatrixField = SquareMatrix< Polynomial::Polynomial< Dim , Degree , Real > , Dim >;

		static const unsigned int VertexNum = Dim+1;
		static const unsigned int FaceNum = Dim+1;

		static const Point< double , Dim > &Vertex( unsigned int d );
		static const SimplexIndex< Dim-1 , unsigned int > &Face( unsigned int f );

		// Returns the affine transformation associated with the permutation of indices
		static Matrix< double , Dim+1 , Dim > AffineTransform( Permutation< Dim+1 > p );

		template< unsigned int EmbeddingDimension >
		static SquareMatrix< double , Dim > Metric( const Simplex< double , EmbeddingDimension , Dim > &s );

		template< unsigned int EmbeddingDimension , typename PositionFunctor >
		static SquareMatrix< double , Dim > Metric( PositionFunctor positionFunctor );

		template< unsigned int SubDim >
		static SquareMatrix< double , SubDim > RestrictedMetric( const SquareMatrix< double , Dim > &g , SimplexIndex< SubDim , unsigned int > &subSimplex );

		static Point< double , Dim > FaceNormal( unsigned int f , SquareMatrix< double , Dim > g=SquareMatrix< double , Dim >::Identity() );

		template< unsigned int Degree , typename Real >
		static double Integral( const PFunction< Degree , Real > &P , SquareMatrix< double , Dim > g = SquareMatrix< double , Dim >::Identity() );

		template< unsigned int Degree , typename Real >
		static PVectorField< (Degree>1) ? Degree-1 : 0 , Real > Gradient( const PFunction< Degree , Real > &P , SquareMatrix< double , Dim > g=SquareMatrix< double , Dim >::Identity() );

		template< unsigned int Degree , typename Real >
		static PFunction< (Degree>1) ? Degree-1 : 0 , Real > Divergence( const PVectorField< Degree , Real > &V , SquareMatrix< double , Dim > g=SquareMatrix< double , Dim >::Identity() );

		template< unsigned int Degree , typename Real >
		static PMatrixField< (Degree>2) ? Degree-2 : 0 , Real > Hessian( const PFunction< Degree , Real > &P , SquareMatrix< double , Dim > g=SquareMatrix< double , Dim >::Identity() );

		template< unsigned int Degree , typename Real >
		static PFunction< Degree , Real > VectorFieldComponent( const PVectorField< Degree , Real > &V , Point< double , Dim > v , SquareMatrix< double , Dim > g=SquareMatrix< double , Dim >::Identity() );

		template< unsigned int Degree , typename Real >
		static PFunction< Degree , Real > PushForward( const PFunction< Degree , Real > &P , const Simplex< double , Dim , Dim > &s );

		template< unsigned int Degree , typename Real >
		static PVectorField< Degree , Real > PushForward( const PVectorField< Degree , Real > &V , const Simplex< double , Dim , Dim > &s );

		static Point< double , Dim > PushForward( Point< double , Dim > p , bool tangentVector , const Simplex< double , Dim , Dim > &s );

		static void GSOrthogonalize( Point< double , Dim > frame[Dim] , SquareMatrix< double , Dim > g=SquareMatrix< double , Dim >::Identity() );
	protected:
		struct Data
		{
			// [BEWARE SIOF]
			Data( void );
			Point< double , Dim > vertices[VertexNum];
			SimplexIndex< Dim-1 , unsigned int > faces[FaceNum];
		};
	};

	template< unsigned int D , unsigned int K >
	constexpr typename std::enable_if< K==0 , unsigned int >::type Choose( void ){ return 1; }
	template< unsigned int D , unsigned int K >
	constexpr typename std::enable_if< K!=0 , unsigned int >::type Choose( void ){ return ( Choose< D-1 , K-1 >() * D ) / K; }

	template< unsigned int Dim , unsigned int Degree >
	struct SimplexElements
	{
		static const unsigned int NodeNum = Choose< Degree+Dim , Dim >();
		using SystemMatrix = SquareMatrix< double , SimplexElements< Dim , Degree >::NodeNum >;

		static void SetElements( Polynomial::Polynomial< Dim , Degree , double > elements[ NodeNum ] );

		static double Volume( SquareMatrix< double , Dim > g );
		static SystemMatrix MassMatrix( SquareMatrix< double , Dim > g );
		static SystemMatrix GradientSquareNormMatrix( SquareMatrix< double , Dim > g );
		static SystemMatrix HessianSquareNormMatrix( SquareMatrix< double , Dim > g );
		static SystemMatrix LaplacianSquareNormMatrix( SquareMatrix< double , Dim > g );

		// For every node and every face, compute the gradient, and restrict the gradient functions to the face.
		static Matrix< Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > , SimplexElements< Dim , Degree >::NodeNum , Dim+1 > FaceGradients( SquareMatrix< double , Dim > g );
		static Matrix< Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > , SimplexElements< Dim , Degree >::NodeNum , Dim+1 > FaceGradients( SquareMatrix< double , Dim > g , const unsigned int *indices );

		// For every node and every face, compute the components of the gradient in an orthonormal frame, and restrict the component functions to the face.
		static Matrix< Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > , SimplexElements< Dim , Degree >::NodeNum , Dim+1 > FaceGradientOrthogonalComponents( SquareMatrix< double , Dim > g );
		static Matrix< Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > , SimplexElements< Dim , Degree >::NodeNum , Dim+1 > FaceGradientOrthogonalComponents( SquareMatrix< double , Dim > g , const unsigned int *indices );

		static Point< SquareMatrix< double , Dim-1 > , Dim+1 > FaceMetrics( SquareMatrix< double , Dim > g );
		static Point< SquareMatrix< double , Dim-1 > , Dim+1 > FaceMetrics( SquareMatrix< double , Dim > g , const unsigned int *indices );

		static Point< double , Dim > NodePosition( unsigned int n );
		static void FactorNodeIndex( unsigned int nodeIndex , unsigned int v[Degree] );
		static unsigned int NodeIndex( const unsigned int v[Degree] );

		static const Polynomial::Polynomial< Dim , Degree , double > & Element( unsigned int n ){ static Data data ; return data.elements[n]; }
		static const Point< Polynomial::Polynomial< Dim   , Degree-1 , double > , Dim > & Differential( unsigned int n ){ static Data data ; return data.differentials[n]; }
		static const Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > & FaceDifferential( unsigned int n , unsigned int f ){ static Data data ; return data.faceDifferentials[n][f]; }

	protected:
		struct Data
		{
			Data( void );

			unsigned int nodeEndPoints[ NodeNum ][ Degree>0 ? Degree : 1 ];
			Polynomial::Polynomial< Dim , Degree , double > elements[ NodeNum ];
			Point< Polynomial::Polynomial< Dim   , Degree-1 , double > , Dim > differentials[ NodeNum ];
			Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > faceDifferentials[ NodeNum ][ Dim+1 ];

			Point< double , Dim > nodePosition( unsigned int ) const;
		protected:
			template< unsigned int D=Degree-1 >
			void _initialize( unsigned int indices[Degree] , unsigned int max );
		};

		static unsigned int __Choose( unsigned int D , unsigned int K );
		static unsigned int  _Choose( unsigned int D , unsigned int K );
	};

#include "SimplexBasis.inl"
}
#endif // SIMPLEX_BASIS_INCLUDED