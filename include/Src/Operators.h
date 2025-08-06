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
#pragma once

#include <Misha/Miscellany.h>
#include <Misha/Exceptions.h>
#include "Hierarchy.h"
#include "EdgeIndexing.h"
#include "QuadratureIntegration.h"

namespace MishaK
{
	namespace TSP
	{
		template< typename MatrixReal >
		struct MassAndStiffnessOperators
		{
			SystemCoefficients< MatrixReal > massCoefficients , stiffnessCoefficients;
			std::vector< RasterLine > rasterLines;
			typename GridAtlas<>::IndexConverter indexConverter;

			size_t cols( void ) const { return indexConverter.numCombined(); }
			size_t rows( void ) const { return indexConverter.numCombined(); }

			template< typename Data >
			void operator()( double mWeight , double sWeight , const std::vector< Data > &in , std::vector< Data > &out , bool add=false ) const;

			template< typename Data >
			std::vector< Data > operator()( double mWeight , double sWeight , const std::vector< Data > &in ) const;

			template< typename Data >
			void operator()( double mWeight , double sWeight , ConstPointer( Data ) in , Pointer( Data ) out , bool add=false ) const;

			template< typename Data >
			void mass( const std::vector< Data > &in , std::vector< Data > &out , bool add=false ) const;

			template< typename Data >
			std::vector< Data > mass( const std::vector< Data > &in ) const;

			template< typename Data >
			void mass( ConstPointer( Data ) in , Pointer( Data ) out , bool add=false ) const;


			template< typename Data >
			void stiffness( const std::vector< Data > &in , std::vector< Data > &out , bool add=false ) const;

			template< typename Data >
			std::vector< Data > stiffness( const std::vector< Data > &in ) const;

			template< typename Data >
			void stiffness( ConstPointer( Data ) in , Pointer( Data ) out , bool add=false ) const;


			template< typename OutReal=MatrixReal >
			Eigen::SparseMatrix< OutReal > operator()( double mWeight , double sWeight ) const;

			template< typename OutReal=MatrixReal >
			Eigen::SparseMatrix< OutReal > mass( void ) const;

			template< typename OutReal=MatrixReal >
			Eigen::SparseMatrix< OutReal > stiffness( void ) const;

		protected:
			template< bool Add , bool Mass , bool Stiffness , typename Data >
			void _evaluate( double mWeight , double sWeight , ConstPointer( Data ) in , Pointer( Data ) out ) const;

			template< bool Mass , bool Stiffness , typename OutReal >
			Eigen::SparseMatrix< OutReal > _matrix( double mWeight , double sWeight ) const;
		};

		template< typename MatrixReal >
		struct DivergenceOperator
		{
			class DivergenceRasterLine
			{
			public:
				int prevEdgeRowStart;
				int currEdgeRowStart;
				int nextEdgeRowStart;
#ifdef SANITY_CHECK
#pragma message( "[WARNING] Is this the right type?" )
#endif // SANITY_CHECK
				AtlasInteriorTexelIndex deepCoefficientsStart;
				AtlasTexelIndex texelStart;
				AtlasTexelIndex texelEnd;

				static std::vector< DivergenceRasterLine > GetRasterLines
				(
					const Map< SimplexIndex< 1 , AtlasTexelIndex > , unsigned int > &coarseEdgeIndex ,
					const std::vector< RasterLine > &rasterLines
				);
			};

			std::vector< SimplexIndex< 1 , AtlasTexelIndex > > edges;
			SparseMatrix< MatrixReal , int > boundaryMatrix;
			std::vector< MatrixReal > deepCoefficients;
			std::vector< DivergenceRasterLine > rasterLines;

			size_t cols( void ) const { return edges.size(); }
			size_t rows( void ) const { return boundaryMatrix.rows; }

			template< typename Data >
			void operator()( const std::vector< Data > &edgeValues , std::vector< Data > &texelDivergence , bool add=false ) const;

			template< typename Data >
			void operator()( ConstPointer( Data ) edgeValues , Pointer( Data ) texelDivergence , bool add=false ) const;

			template< typename Data >
			std::vector< Data > operator()( const std::vector< Data > &edgeValues ) const;

			template< typename OutReal=MatrixReal >
			Eigen::SparseMatrix< OutReal > operator()( void ) const;
		};

		template< typename Real , typename SampleType >
		struct Integrator
		{
			template< typename OutData , typename InData >
			struct Scratch
			{
			protected:
				std::vector< InData > _coarseBoundaryPrimal;
				std::vector< InData > _fineBoundaryPrimal;
				std::vector< OutData > _coarseBoundaryDual;
				std::vector< OutData > _fineBoundaryDual;
				friend Integrator;
			};

			typename GridAtlas<>::IndexConverter indexConverter;
			SparseMatrix< Real , int > coarseBoundaryFineBoundaryProlongation;
			SparseMatrix< Real , int > fineBoundaryCoarseBoundaryRestriction;
			std::vector< InteriorCellLine > interiorCellLines;
			SampleType samples;

			template< typename OutData , typename InData >
			Scratch< OutData , InData > getScratch( void ) const;

			template< typename OutData , typename InData , typename SampleFunction /* = std::function< OutData ( InData , SquareMatrix< Real , 2 > ) */ >
			void operator()( const std::vector< InData > &primal , const SampleFunction & SF , std::vector< OutData > &dual , bool add=false ) const;

			template< typename OutData , typename InData , typename SampleFunction /* = std::function< OutData ( InData , SquareMatrix< Real , 2 > ) */ >
			std::vector< OutData > operator()( const std::vector< InData > &primal , const SampleFunction & SF ) const;

			template< typename OutData , typename InData , typename SampleFunction /* = std::function< OutData ( InData , SquareMatrix< Real , 2 > ) */ >
			void operator()( ConstPointer( InData ) primal , const SampleFunction & SF , Pointer( OutData ) dual , bool add=false ) const;


			template< typename OutData , typename InData , typename SampleFunction /* = std::function< OutData ( InData , SquareMatrix< Real , 2 > ) */ >
			void operator()( const std::vector< InData > &primal , const SampleFunction & SF , Scratch< OutData , InData > & scratch , std::vector< OutData > &dual , bool add=false ) const;

			template< typename OutData , typename InData , typename SampleFunction /* = std::function< OutData ( InData , SquareMatrix< Real , 2 > ) */ >
			std::vector< OutData > operator()( const std::vector< InData > &primal , const SampleFunction & SF , Scratch< OutData , InData > & scratch ) const;

			template< typename OutData , typename InData , typename SampleFunction /* = std::function< OutData ( InData , SquareMatrix< Real , 2 > ) */ >
			void operator()( ConstPointer( InData ) primal , const SampleFunction & SF , Scratch< OutData , InData > & scratch , Pointer( OutData ) dual , bool add=false ) const;
		};

		template< typename Real >
		using ScalarIntegrator = Integrator< Real , ScalarElementSamples< Real > >;

		template< typename Real >
		using GradientIntegrator = Integrator< Real , GradientElementSamples< Real > >;

		struct OperatorInitializer
		{
			template< typename GeometryReal , typename MatrixReal >
			static void Initialize
			(
				unsigned int samples ,
				MassAndStiffnessOperators< MatrixReal > & massAndStiffnessOperators ,
				const GridAtlas< GeometryReal , MatrixReal > & gridAtlas ,
				const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > & parameterMetric ,
				const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > & atlasCharts
			);

			template< typename GeometryReal , typename MatrixReal >
			static void Initialize
			(
				unsigned int samples ,
				MassAndStiffnessOperators< MatrixReal > & massAndStiffnessOperators ,
				const GridAtlas< GeometryReal , MatrixReal > & gridAtlas ,
				const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > & parameterMetric ,
				const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > & atlasCharts ,
				DivergenceOperator< MatrixReal > & divergenceOperator 
			);

			template< typename GeometryReal , typename MatrixReal >
			static void Initialize
			(
				unsigned int samples ,
				MassAndStiffnessOperators< MatrixReal > & massAndStiffnessOperators ,
				const GridAtlas< GeometryReal , MatrixReal > & gridAtlas ,
				const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > & parameterMetric ,
				const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > & atlasCharts ,
				const BoundaryProlongationData< MatrixReal > & boundaryProlongation
			);

			template< typename GeometryReal , typename MatrixReal >
			static void Initialize
			(
				unsigned int samples ,
				MassAndStiffnessOperators< MatrixReal > & massAndStiffnessOperators ,
				const GridAtlas< GeometryReal , MatrixReal > & gridAtlas ,
				const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > & parameterMetric ,
				const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > & atlasCharts ,
				const BoundaryProlongationData< MatrixReal > &boundaryProlongation ,
				DivergenceOperator< MatrixReal > & divergenceOperator 
			);

			template< typename GeometryReal , typename MatrixReal , typename SampleType >
			static void Initialize
			(
				unsigned int samples ,
				MassAndStiffnessOperators< MatrixReal > & massAndStiffnessOperators ,
				const GridAtlas< GeometryReal , MatrixReal > & gridAtlas ,
				const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > & parameterMetric ,
				const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > & atlasCharts ,
				Integrator< MatrixReal , SampleType > & integrator ,
				unsigned int numSamples ,
				bool approximate
			);

			template< typename GeometryReal , typename MatrixReal , typename SampleType >
			static void Initialize
			(
				unsigned int samples ,
				MassAndStiffnessOperators< MatrixReal > & massAndStiffnessOperators ,
				const GridAtlas< GeometryReal , MatrixReal > & gridAtlas ,
				const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > & parameterMetric ,
				const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > & atlasCharts ,
				DivergenceOperator< MatrixReal > & divergenceOperator ,
				Integrator< MatrixReal , SampleType > & integrator ,
				unsigned int numSamples ,
				bool approximate
			);

		protected:

			template< unsigned int Samples , typename GeometryReal , typename MatrixReal >
			static void _InitializeChart
			(
				const ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > &texture_metrics ,
				const AtlasChart< GeometryReal > &atlasChart ,
				const GridChart< GeometryReal > &gridChart ,
				const typename GridAtlas<>::IndexConverter & indexConverter ,
				const std::vector< AtlasInteriorOrBoundaryNodeIndex >& fineBoundaryIndex ,
				std::vector< MatrixReal >& deepMassCoefficients ,
				std::vector< MatrixReal >& deepStiffnessCoefficients ,
				std::vector< MassAndStiffnessCoefficient< MatrixReal > > & boundaryBoundaryMassAndStiffness ,
				std::vector< MassAndStiffnessCoefficient< MatrixReal > > & boundaryDeepMassAndStiffness ,
				bool computeDivergence ,
				Map< SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > , AtlasRefinedBoundaryEdgeIndex > & fineBoundaryEdgeIndex,
				Map< SimplexIndex< 1 , AtlasTexelIndex > , unsigned int > & coarseEdgeIndex,
				std::vector< Eigen::Triplet< MatrixReal > > & boundaryDeepDivergenceTriplets,
				std::vector< Eigen::Triplet< MatrixReal > > & boundaryBoundaryDivergenceTriplets,
				std::vector< MatrixReal > & deepDivergenceCoefficients
			);

			template< unsigned int Samples , typename GeometryReal , typename MatrixReal >
			static void _Initialize
			(
				const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > &parameterMetric ,
				const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > &atlasCharts ,
				const GridAtlas< GeometryReal , MatrixReal > &gridAtlas ,
				const std::vector< AtlasInteriorOrBoundaryNodeIndex > &fineBoundaryIndex ,
				int numFineBoundaryNodes ,
				std::vector< MatrixReal >& deepMassCoefficients ,
				std::vector< MatrixReal >& deepStiffnessCoefficients ,
				SparseMatrix< MatrixReal , int >& boundaryBoundaryMassMatrix ,
				SparseMatrix< MatrixReal , int >& boundaryBoundaryStiffnessMatrix ,
				SparseMatrix< MatrixReal , int >& boundaryDeepMassMatrix ,
				SparseMatrix< MatrixReal , int >& boundaryDeepStiffnessMatrix ,
				bool computeDivergence ,
				Map< SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > , AtlasRefinedBoundaryEdgeIndex > & fineBoundaryEdgeIndex,
				Map< SimplexIndex< 1 , AtlasTexelIndex > , unsigned int > & coarseEdgeIndex,
				std::vector< Eigen::Triplet< MatrixReal > > & boundaryDeepDivergenceTriplets,
				std::vector< Eigen::Triplet< MatrixReal > > & boundaryBoundaryDivergenceTriplets,
				std::vector< MatrixReal >& deepDivergenceCoefficients
			);

			template< unsigned int Samples , typename GeometryReal , typename MatrixReal >
			static void _Initialize
			(
				MassAndStiffnessOperators< MatrixReal > & massAndStiffnessOperators ,
				const GridAtlas< GeometryReal , MatrixReal > & gridAtlas ,
				const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > &parameterMetric ,
				const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > &atlasCharts ,
				const BoundaryProlongationData< MatrixReal > &boundaryProlongation ,
				bool computeDivergence ,
				DivergenceOperator< MatrixReal > & divergenceOperator
			);

			template< typename GeometryReal , typename MatrixReal >
			static void _Initialize
			(
				unsigned int samples ,
				MassAndStiffnessOperators< MatrixReal > & massAndStiffnessOperators ,
				const GridAtlas< GeometryReal , MatrixReal > & gridAtlas ,
				const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > &parameterMetric ,
				const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > &atlasCharts ,
				const BoundaryProlongationData< MatrixReal > &boundaryProlongation ,
				bool computeDivergence ,
				DivergenceOperator< MatrixReal > & divergenceOperator
			);
		};

#include "Operators.inl"
#include "UpdateCoefficients.inl"
	}
}
