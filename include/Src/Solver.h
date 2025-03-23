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

#if defined( USE_EIGEN ) || defined( USE_EIGEN_PARDISO )
#include <Eigen/Sparse>
#include <Eigen/Dense>
#endif // USE_EIGEN || USE_EIGEN_PARDISO
#ifdef USE_EIGEN_PARDISO
#include <Eigen/PardisoSupport>
#endif // USE_EIGEN_PARDISO
#if defined( USE_CHOLMOD ) || defined( USE_EIGEN )
#include <Misha/LinearSolvers.h>
#endif // USE_CHOLMOD


namespace MishaK
{
#ifdef USE_CHOLMOD

#define CHOLMOD_CHANNELS_IN_PARALLEL

#ifdef CHOLMOD_CHANNELS_IN_PARALLEL
	template< class Real , unsigned int Channels >
	class CholmodCholeskySolver
	{
	public:
		CholmodSolver<1> solver[Channels];
		std::vector< Real > out[Channels];
		std::vector< Real > in[Channels];
		template< typename _Real >
		void init( const SparseMatrix< _Real , int > & M )
		{
			ThreadPool::ParallelFor( 0 , Channels , [&]( unsigned int , size_t c ){ solver[c]._init(M); } );

			size_t numVariables = M.Rows();
			for( int c=0 ; c<Channels ; c++ ) out[c].resize( numVariables ) , in[c].resize( numVariables );
		}

		template< typename _Real >
		void update( const SparseMatrix< _Real , int >& M )
		{
			ThreadPool::ParallelFor( 0 , Channels , [&]( unsigned int , size_t c ){ solver[c]._update(M); } );
		}
	};

	template< class Real , unsigned int Channels , class DataType >
	void solve( CholmodCholeskySolver< Real , Channels >& chol , std::vector< DataType >& x0 , const std::vector< DataType >& rhs )
	{
		int numVariables = x0.size();
		ThreadPool::ParallelFor
		(
			0 , Channels ,
			[&]( unsigned int , size_t c )
			{
				for( int n=0 ; n<numVariables ; n++ ) chol.in[c][n] = rhs[n][c];
				chol.solver[c].solve( &chol.in[c][0] , &chol.out[c][0] );
			}
		);
		ThreadPool::ParallelFor( 0 , numVariables , [&]( unsigned int , size_t n ){ for( int c=0 ; c<Channels ; c++ ) x0[n][c] = chol.out[c][n]; } );
	}
	template< class Real >
	void solve( CholmodCholeskySolver< Real , 1 >& chol , std::vector< Real >& x0 , const std::vector< Real >& rhs )
	{
		size_t numVariables = x0.size();
		ThreadPool::ParallelFor( 0 , numVariables , [&]( unsigned int , size_t n ){ chol.in[0][n] = rhs[n]; } );
		chol.solver[0].solve( &chol.in[0][0] , &chol.out[0][0] );
		ThreadPool::ParallelFor( 0 , numVariables , [&]( unsigned int , size_t n ){ x0[n] = chol.out[0][n]; } );
	}
#else // !CHOLMOD_CHANNELS_IN_PARALLEL

#define CHOLMOD_CHANNELS_IN_BLOCK
#ifdef CHOLMOD_CHANNELS_IN_BLOCK
	template< class Real , unsigned int Channels >
	class CholmodCholeskySolver
	{
	public:
		CholmodSolver< Channels > solver;
		std::vector< Real > out;
		std::vector< Real > in;
		template< typename _Real >
		void init( const SparseMatrix< _Real , int >& M )
		{
			solver._init(M);
			int numVariables = M.Rows();
			out.resize( Channels*numVariables );
			in.resize( Channels*numVariables );
		}

		template< typename _Real >
		void update( const SparseMatrix< _Real , int >& M ){ solver._update(M); }
	};

	template< class Real , unsigned int Channels , class DataType >
	void solve( CholmodCholeskySolver< Real , Channels > & chol , std::vector< DataType >& x0 , const std::vector< DataType >& rhs )
	{
		int numVariables = x0.size();
		for( int n=0 ; n<numVariables ; n++ ) for( int c=0 ; c<Channels ; c++ ) chol.in[ c*numVariables+n ] = rhs[n][c];
		chol.solver.solve( &chol.in[0] , &chol.out[0] );
		ThreadPool::ParallelForm( 0 , numVariables , [&]( unsigned int , size_t n ){ x0[n][c] = chol.out[ c*numVariables+n ]; } );
	}
	template< class Real >
	void solve( CholmodCholeskySolver< Real , 1 >& chol , std::vector< Real >& x0 , const std::vector< Real >& rhs )
	{
		int numVariables = x0.size();
		ThreadPool::ParallelForm( 0 , numVariables , [&]( unsigned int , size_t n ){ chol.in[n] = rhs[n]; } );
		chol.solver.solve( &chol.in[0] , &chol.out[0] );
		ThreadPool::ParallelForm( 0 , numVariables , [&]( unsigned int , size_t n ){ x0[n] = chol.out[n]; } );
	}
#else // !CHOLMOD_CHANNELS_IN_BLOCK
	template< class Real , unsigned int Channels >
	class CholmodCholeskySolver
	{
	public:
		CholmodSolver< 1 > solver;
		std::vector< Real > out[Channels];
		std::vector< Real > in[Channels];
		template< typename _Real >
		void init( const SparseMatrix< _Real , int >& M )
		{
			solver._init(M);
			int numVariables = M.Rows();
			for( int c=0 ; c<Channels ; c++ ) out[c].resize( numVariables ) , in[c].resize( numVariables );
		}
		template< typename _Real >
		void update( const SparseMatrix< _Real , int >& M ){ solver._update(M); }
	};

	template< class Real , unsigned int Channels , class DataType >
	void solve( CholmodCholeskySolver< Real , Channels >& chol , std::vector< DataType >& x0 , const std::vector< DataType >& rhs )
	{
		int numVariables = x0.size();
		for( int c=0 ; c<Channels ; c++ )
		{
			ThreadPool::ParallelFor( 0 , numVariables , [&]( unsigned int , size_t n ){ chol.in[c][n] = rhs[n][c]; } );
			chol.solver.solve( &chol.in[c][0] , &chol.out[c][0] );
		}
		ThreadPool::ParallelFor( 0 , numVariables , [&]( unsigned int , size_t n ){ for( int c=0 ; c<Channels ; c++ ) x0[n][c] = chol.out[c][n]; } );
	}
	template< class Real >
	void solve( CholmodCholeskySolver< Real , 1 >& chol , std::vector< Real >& x0 , const std::vector< Real >& rhs )
	{
		int numVariables = x0.size();
		ThreadPool::ParallelFor( 0 , numVariables , [&]( unsigned int , size_t n ){ chol.in[0][n] = rhs[n]; } );
		chol.solver.solve( &chol.in[0][0] , &chol.out[0][0] );
		ThreadPool::ParallelFor( 0 , numVariables , [&]( unsigned int , size_t n ){ x0[n] = chol.out[0][n]; } );
	}
#endif // CHOLMOD_CHANNELS_IN_BLOCK
#endif // CHOLMOD_CHANNELS_IN_PARALLEL
#endif // USE_CHOLMOD


#ifdef USE_EIGEN
	template< typename Real , unsigned int Channels >
	class EigenCholeskySolver
	{
	public:
#if 1
		typedef Eigen::Matrix< double , Eigen::Dynamic , 1 > EigenVector;
		typedef EigenSolverCholeskyLDLt< Real , ConstPointer( MatrixEntry< Real , int > ) > EigenSolver;
#else
		typedef Eigen::Matrix< double , Eigen::Dynamic , 1 > EigenVector;
		typedef EigenSolverCholeskyLDLt< Real , ConstPointer( MatrixEntry< Real , int > ) , double > EigenSolver;
#endif
		EigenSolver* solver;
		EigenVector b[Channels] , x[Channels];

		EigenCholeskySolver( void ) : solver(NULL){ }
		~EigenCholeskySolver( void ){ if( solver ) delete solver; }
		void init( const SparseMatrix< Real , int >& M )
		{
			if( solver ) delete solver;
			solver = new EigenSolver( M , true );

			int numVariables = (int)M.Rows();
			for( int c=0 ; c<Channels ; c++ ) b[c].resize( numVariables ) , x[c].resize( numVariables );
		}
		void update( const SparseMatrix< Real , int >& M ){ solver->update( M ); }
	};

	template< class Real , unsigned int Channels , class DataType >
	void solve( EigenCholeskySolver< Real , Channels >& chol , std::vector< DataType >& x , const std::vector< DataType >& b )
	{
		int numVariables = (int)x.size();
		ThreadPool::ParallelFor( 0 , numVariables , [&]( unsigned int , size_t n ){ for( int c=0 ; c<Channels ; c++ ) chol.b[c][n] = b[n][c]; } );
#if 0
		ThreadPool::ParallelFor( 0 , Channels     , [&]( unsigned int , size_t c ){ chol.solver->solve( &chol.b[c][0] , &chol.x[c][0] ); } );
#else
		ThreadPool::ParallelFor( 0 , Channels     , [&]( unsigned int , size_t c ){ chol.solver->solve( chol.b[c] , chol.x[c] ); } );
#endif
		ThreadPool::ParallelFor( 0 , numVariables , [&]( unsigned int , size_t n ){ for( int c=0 ; c<Channels ; c++ ) x[n][c] = (Real)chol.x[c][n];} );
	}

	template< class Real >
	void solve( EigenCholeskySolver< Real , 1 >& chol , std::vector< Real >& x0 , const std::vector< Real >& rhs ){ chol.solver->solve( ( ConstPointer( Real ) )GetPointer( rhs ) , GetPointer( x0 ) ); }
#endif // USE_EIGEN


#ifdef USE_EIGEN_PARDISO
#pragma comment( lib , "mkl_intel_lp64.lib" )
#pragma comment( lib , "mkl_intel_thread.lib" )
#pragma comment( lib , "mkl_core.lib" )
#pragma comment( lib , "libiomp5md.lib" )

	//Pardiso Solver

	template< class Real , unsigned int Channels >
	class EigenPardisoSolver
	{
	public:
		Eigen::PardisoLDLT< Eigen::SparseMatrix< Real > >* solver;
		EigenVector x0_vectors[Channels];
		EigenVector rhs_vectors[Channels];
		EigenVector solution_vectors[Channels];

		template< typename _Real >
		void init( const SparseMatrix< _Real , int >& _M )
		{
			Eigen::SparseMatrix< Real > M;
			SparseMatrixParser( _M , M );

			solver = new Eigen::PardisoLDLT< Eigen::SparseMatrix< Real > >();
			solver->analyzePattern(M);
			switch( solver->info() )
			{
			case Eigen::Success: break;
			case Eigen::NumericalIssue: Miscellany::Throw( "Numerical issue" );
			case Eigen::NoConvergence:  Miscellany::Throw( "No convergence" );
			case Eigen::InvalidInput:   Miscellany::Throw( "Invalid input" );
			default:                    Miscellany::Throw( "Undetermined cause" );
			}

			int numVariables = M.rows();
			for( int c=0 ; c<Channels ; c++ )
			{
				x0_vectors[c].resize( numVariables );
				rhs_vectors[c].resize( numVariables );
				solution_vectors[c].resize( numVariables );
			}
		}
		template< typename _Real >
		void update( const SparseMatrix< _Real , int >& _M )
		{
			Eigen::SparseMatrix< Real > M;
			SparseMatrixParser( _M , M );
			solver->factorize(M);
		}
	};

	template< class Real , unsigned int Channels , class DataType >
	void solve( EigenPardisoSolver< Real , Channels >& chol , std::vector< DataType >& x0 , const std::vector< DataType >& rhs )
	{
		int numVariables = x0.size();
		ThreadPool::ParallelFor( 0 , numVariables , [&]( unsigned int , size_t n ){ for( int c=0 ; c<Channels ; c++ ) chol.x0_vectors[c][n] = x0[n][c] , chol.rhs_vectors[c][n] = rhs[n][c]; } );
		for( int c=0 ; c<3 ; c++ ) chol.solution_vectors[c] = chol.solver->solve( chol.rhs_vectors[c] );
		ThreadPool::ParallelFor( 0 , numVariables , [&]( unsigned int , size_t n ){ for( int c=0 ; c<Channels ; c++ ) x0[n][c] = chol.solution_vectors[c][n]; } );
	}
	template< class Real >
	void solve( EigenPardisoSolver< Real , 1 >& chol , std::vector< Real >& x0 , const std::vector< Real >& rhs )
	{
		int numVariables = x0.size();
		ThreadPool::ParallelFor( 0 , numVariables , [&]( unsigned int , size_t n ){ chol.x0_vectors[0][n] = x0[n] , chol.rhs_vectors[0][n] = rhs[n]; } );
		chol.solution_vectors[0] = chol.solver->solve( chol.rhs_vectors[0] );
		ThreadPool::ParallelFor( 0 , numVariables , [&]( unsigned int , size_t n ){ x0[n] = chol.solution_vectors[0][n]; } );
	}
#endif // USE_EIGEN_PARDISO
}