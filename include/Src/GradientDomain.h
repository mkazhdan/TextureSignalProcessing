/*
Copyright (c) 2025, Michael Kazhdan
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

#ifndef GRADIENT_DOMAIN_INCLUDED
#define GRADIENT_DOMAIN_INCLUDED

#include <Eigen/Sparse>

#include "Misha/Geometry.h"
#include "Misha/MultiThreading.h"
#include "Src/Operators.h"
#include "Src/Hierarchy.h"


namespace MishaK
{
	namespace TSP
	{
		template< typename Real >
		struct GradientDomain
		{
			// A constructor taking functors mapping:
			// 1.  Triangle index + corner index -> surface vertex index
			// 2a. Triangle index                -> metric tensor associated to the triangle
			// 2b. Surface vertex index          -> position in 3D
			// 3.  Triangle index + corner index -> texture vertex index
			// 4.  Texture vertex index          -> position in the unit square
			template
				<
				typename SurfaceCornerFunctor ,         /* = std::function< size_t ( size_t , unsigned int ) > */
				typename SurfaceVertexOrMetricFunctor , /* = std::function< Point< Real , 3 > ( size_t ) > || std::function< SquareMatrix< Real , 2 > ( size_t ) > */
				typename TextureCornerFunctor ,         /* = std::function< size_t ( size_t , unsigned int ) > */
				typename TextureVertexFunctor           /* = std::function< Point< Real , 2 > ( size_t ) > */
				>
				GradientDomain
				(
					unsigned quadraturePointsPerTriangle ,
					size_t numTriangles ,
					size_t numSurfaceVertices ,
					size_t numTextureVertices ,
					SurfaceCornerFunctor         && surfaceCornerFunctor ,
					SurfaceVertexOrMetricFunctor && surfaceVertexOrMetricFunctor ,
					TextureCornerFunctor         && textureCornerFunctor ,
					TextureVertexFunctor         && textureVertexFunctor ,
					unsigned int width ,
					unsigned int height ,
					bool normalize = true
				);

			// The number of texels
			size_t numNodes( void ) const;

			// The number of edges between texels
			size_t numEdges( void ) const;

			// The 2D index of a texel
			std::pair< unsigned int , unsigned int > node( size_t n ) const;

			// Is the texel covered by a triangle?
			bool isCovered( size_t n ) const;

			// The indices of the two texels that are the end-points of the edge
			std::pair< size_t , size_t > edge( size_t e ) const;

			// The mass matrix
			Eigen::SparseMatrix< Real > mass( void ) const;

			// The stiffness matrix
			Eigen::SparseMatrix< Real > stiffness( void ) const;

			// The divergence matrix
			Eigen::SparseMatrix< Real > divergence( void ) const;

			// The finite-differences matrix
			Eigen::SparseMatrix< Real > finiteDifferences( void ) const;

			// Computes the mass of the per-texel signal
			template< typename T > void mass( const T * in , T * out ) const;

			// Computes the stiffness of the per-texel signal
			template< typename T > void stiffness( const T * in , T * out ) const;

			// Computes the divergence of the per-edge signal
			template< typename T > void divergence( const T * in , T * out ) const;

			// Compute the finite-differences of texels along the edges
			template< typename T > void finiteDifferences( const T * in , T * out ) const;

			// Runs unit-tests to confirm the definitions of the operators are consistent
			void unitTests( unsigned int numTest , double eps , bool verbose=false ) const;

		protected:
			template< typename Functor > static constexpr bool _IsTriangleCornerFunctor( void );
			template< typename Functor > static constexpr bool _IsSurfaceMetricFunctor( void );
			template< typename Functor > static constexpr bool _IsSurfaceVertexFunctor( void );	
			template< typename Functor > static constexpr bool _IsSurfaceVertexOrMetricFunctor( void );
			template< typename Functor > static constexpr bool _IsTextureVertexFunctor( void );

			MassAndStiffnessOperators< Real > _massAndStiffnessOperators;
			DivergenceOperator< Real > _divergenceOperator;
			std::vector< TextureNodeInfo< Real > > _textureNodes;
		};
#include "GradientDomain.inl"
	}
}
#endif // GRADIENT_DOMAIN_INCLUDED