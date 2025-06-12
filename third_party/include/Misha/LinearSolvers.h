/*
Copyright (c) 2023, Michael Kazhdan
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
#ifndef LINEAR_SOLVERS_INCLUDE
#define LINEAR_SOLVERS_INCLUDE

#ifdef USE_PARDISO
extern "C" void pardisoinit_d(void   *, int    *,   int *, int *, double *, int * );
extern "C" void pardiso_d(void   *, int    *,   int *, int *,    int *, int *, double *, int * , int * , int * , int * , int *, int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix_d(int *, int *, double *, int *, int *, int * );
extern "C" void pardiso_chkvec_d(int *, int *, double *, int * );
extern "C" void pardiso_printstats_d(int *, int *, double *, int *, int *, int *, double *, int * );
extern "C" void pardiso_get_schur_d(void*, int*, int*, int*, double*, int*, int* );
#pragma comment( lib , "Pardiso/libpardiso500-WIN-X86-64.lib" )
#endif // USE_PARDISO

#ifdef USE_CHOLMOD
#ifdef USE_SUITESPARSE
#include <Cholmod.SuiteSparse/cholmod.h>
#if defined( WIN32 ) || defined( _WIN64 )
#pragma message( "[WARNING] Need to explicitly exclude VCOMP.lib" )
#pragma comment( lib , "SuiteSparse/libamd.lib" )
//#pragma comment( lib , "SuiteSparse/libbtf.lib" )
#pragma comment( lib , "SuiteSparse/libcamd.lib" )
#pragma comment( lib , "SuiteSparse/libccolamd.lib" )
#pragma comment( lib , "SuiteSparse/libcholmod.lib" )
#pragma comment( lib , "SuiteSparse/libcolamd.lib" )
//#pragma comment( lib , "SuiteSparse/libcxsparse.lib" )
//#pragma comment( lib , "SuiteSparse/libklu.lib" )
//#pragma comment( lib , "SuiteSparse/libldl.lib" )
//#pragma comment( lib , "SuiteSparse/libspqr.lib" )
//#pragma comment( lib , "SuiteSparse/libumfpack.lib" )
#pragma comment( lib , "SuiteSparse/metis.lib" )
#pragma comment( lib , "SuiteSparse/suitesparseconfig.lib" )
#pragma comment( lib , "SuiteSparse/libblas.lib" )
#pragma comment( lib , "SuiteSparse/liblapack.lib" )
#endif // WIN32 || _WIN64
#ifdef DLONG
typedef long long SOLVER_LONG;
#define CHOLMOD( name ) cholmod_l_ ## name
#else // !DLONG
typedef       int SOLVER_LONG;
#define CHOLMOD( name ) cholmod_ ## name
#endif // DLONG
#else // !USE_SUITESPARSE
#include <Cholmod/cholmod.h>
#if defined( WIN32 ) || defined( _WIN64 )
#pragma message( "[WARNING] Need to explicitly exclude VCOMP.lib" )
#pragma comment( lib , "CHOLMOD_FULL.lib" )
#endif // WIN32 || _WIN64
#ifdef DLONG
typedef long long SOLVER_LONG;
#define CHOLMOD( name ) cholmod_l_ ## name
#else // !DLONG
typedef       int SOLVER_LONG;
#define CHOLMOD( name ) cholmod_ ## name
#endif // DLONG
#endif // USE_SUITESPARSE
#elif defined(EIGEN_USE_MKL_ALL)
#pragma comment( lib , "mkl_core.lib" )
#pragma comment( lib , "mkl_intel_lp64.lib" )
#pragma comment( lib , "mkl_intel_thread.lib" )
#pragma comment( lib , "mkl_blas95_lp64.lib" )
#pragma comment( lib , "libiomp5md.lib" )
#endif // USE_CHOLMOD

#include "SparseMatrixInterface.h"
#include "MultiThreading.h"

namespace MishaK
{
	inline double                        SquareNorm( const double* values , int dim ){ double norm2 = 0 ; for( int i=0 ; i<dim ; i++ ) norm2 += values[i] * values[i] ; return norm2; }
	inline double                        SquareNorm( const  float* values , int dim ){ double norm2 = 0 ; for( int i=0 ; i<dim ; i++ ) norm2 += values[i] * values[i] ; return norm2; }
	template< class Type > inline double SquareNorm( const   Type* values , int dim ){ double norm2 = 0 ; for( int i=0 ; i<dim ; i++ ) norm2 += values[dim].squareNorm()  ; return norm2 ; }

	inline double                        SquareDifference( const double* values1 , const double* values2 , int dim ){ double norm2 = 0 ; for( int i=0 ; i<dim ; i++ ) norm2 += ( values1[i] - values2[i] ) * ( values1[i] - values2[i] ) ; return norm2; }
	inline double                        SquareDifference( const  float* values1 , const  float* values2 , int dim ){ double norm2 = 0 ; for( int i=0 ; i<dim ; i++ ) norm2 += ( values1[i] - values2[i] ) * ( values1[i] - values2[i] ) ; return norm2; }
	template< class Type > inline double SquareDifference( const   Type* values1 , const   Type* values2 , int dim ){ double norm2 = 0 ; for( int i=0 ; i<dim ; i++ ) norm2 += ( values1[dim] - values2[dim] ).squareNorm()  ; return norm2 ; }


	// This is the conjugate gradients solver.
	// The assumption is that the class SPDOperator defines a method operator()( const Real* , Real* ) which corresponds to applying a symmetric positive-definite operator.
	template< typename T >
	struct CGScratch
	{
		T *r , *d , *q;
		CGScratch( void ) : r(NULL) , d(NULL) , q(NULL) , _dim(0){ ; }
		CGScratch( int dim ) : r(NULL) , d(NULL) , q(NULL){ resize(dim); }
		~CGScratch( void ){ resize(0); }
		void resize( int dim )
		{
			if( dim!=_dim )
			{
				if( r ) delete[] r ; r = NULL;
				if( d ) delete[] d ; d = NULL;
				if( q ) delete[] q ; q = NULL;
				if( dim ) r = new T[dim] , d = new T[dim] , q = new T[dim];
				_dim = dim;
			}
		}
	protected:
		int _dim;
	};
	template< typename T >
	struct PreconditionedCGScratch : public CGScratch< T >
	{
		T *s;
		PreconditionedCGScratch( void ) : CGScratch< T >() , s(NULL){ ; }
		PreconditionedCGScratch( int dim ) : CGScratch< T >() { resize(dim); }
		~PreconditionedCGScratch( void ){ resize(0); }
		void resize( int dim )
		{
			if( dim!=CGScratch< T >::_dim )
			{
				if( s ) delete[] s; s = NULL;
				if( dim ) s = new T[dim];
			}
			CGScratch< T >::resize( dim );
		}
	};
	template< class Real >
	struct DiagonalPreconditioner
	{
		Real* iDiagonal;
		DiagonalPreconditioner( void ) : iDiagonal(NULL) , _dim(0){ ; }
		~DiagonalPreconditioner( void ){ if( iDiagonal ) delete[] iDiagonal ; iDiagonal = NULL; }
		template< class MatrixRowIterator >
		void set( const SparseMatrixInterface< Real , MatrixRowIterator >& M )
		{
			if( _dim!=M.Rows() )
			{
				_dim = (int)M.Rows();
				if( iDiagonal ) delete[] iDiagonal , iDiagonal = NULL;
				if( _dim>0 ) iDiagonal = new Real[_dim];
			}
			memset( iDiagonal , 0 , sizeof(Real)*_dim );
			ThreadPool::ParallelFor
			(
				0 , M.Rows() ,
				[&]( unsigned int , size_t i )
				{
					for( MatrixRowIterator iter=M.begin(i) ; iter!=M.end(i) ; iter++ ) if( iter->N==i ) iDiagonal[i] += iter->Value;
					iDiagonal[i] = (Real)1./iDiagonal[i];
				}
			);
		}
		template< typename T >
		void operator()( const T* in , T* out ) const
		{
			ThreadPool::ParallelFor( 0 , _dim , [&]( unsigned int , size_t i ){ out[i] = in[i] * iDiagonal[i]; } );
		}
	protected:
		int _dim;
	};

	template< class Real , typename T , class SPDOperator , typename TDotT >
	int SolveCG( SPDOperator& L , int iters , int dim , const T* b , T* x , TDotT dot , CGScratch< T >* scratch=NULL , double eps=1e-8 , int threads=1 , bool verbose=false )
	{
		eps *= eps;
		T *r , *d , *q;
		if( scratch ) r = scratch->r , d = scratch->d , q = scratch->q;
		else          r = new T[dim] , d = new T[dim] , q = new T[dim];
		memset( r , 0 , sizeof(T)*dim ) , memset( d , 0 , sizeof(T)*dim ) , memset( q , 0 , sizeof(T)*dim );
		double delta_new = 0 , delta_0;

		L( x , r );
		std::vector< double > _delta_news( ThreadPool::NumThreads() , 0 );
		ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int t , size_t i ){ d[i] = r[i] = b[i] - r[i] ; _delta_news[t] += dot( r[i] , r[i] ); } , threads );
		for( unsigned int t=0 ; t<_delta_news.size() ; t++ ) delta_new += _delta_news[t];

		delta_0 = delta_new;
		if( delta_new<=eps )
		{
			if( !scratch ) delete[] r , delete[] d , delete[] q;
			return 0;
		}

		int ii;
		for( ii=0 ; ii<iters && delta_new>eps*delta_0 ; ii++ )
		{
			L( d , q );
			double dDotQ = 0;
			std::vector< double > _dDotQs( ThreadPool::NumThreads() , 0 );
			ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int t , size_t i ){ _dDotQs[t] += dot( d[i] , q[i] ); } , threads );
			for( unsigned int t=0 ; t<_dDotQs.size() ; t++ ) dDotQ += _dDotQs[t];
			if( !dDotQ ) break;
			Real alpha = Real( delta_new / dDotQ );

			double delta_old = delta_new;
			delta_new = 0;

			const int RESET_COUNT = 50;
			if( (ii%RESET_COUNT)==(RESET_COUNT-1) )
			{
				ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int , size_t i ){ x[i] += d[i] * alpha; } , threads );
				L( x , r );
				std::vector< double > _delta_news( ThreadPool::NumThreads() , 0 );
				ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int t , size_t i ){ r[i] = b[i] - r[i] ; _delta_news[t] += dot( r[i] , r[i] ); } , threads );
				for( unsigned int t=0 ; t<_delta_news.size() ; t++ ) delta_new += _delta_news[t];
			}
			else
			{
				std::vector< double > _delta_news( ThreadPool::NumThreads() , 0 );
				ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int t , size_t i ){ r[i] -= q[i] * alpha ; _delta_news[t] += dot( r[i] , r[i] ) ; x[i] += d[i] * alpha; } , threads );
				for( unsigned int t=0 ; t<_delta_news.size() ; t++ ) delta_new += _delta_news[t];
			}

			Real beta = Real( delta_new / delta_old );
			ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int , size_t i ){ d[i] = r[i] + d[i] * beta; } , threads );
		}
		if( verbose )
		{
			L( x , r );
			ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int , size_t i ){ r[i] -= b[i]; } , threads );
			printf( "CG: %d %g -> %g\n" , ii , SquareNorm( b , dim ) , SquareNorm( r , dim ) );
		}
		if( !scratch ) delete[] r , delete[] d , delete[] q;
		return ii;
	}
	template< class Real , typename T , class SPDOperator , class SPDPreconditioner , typename TDotT >
	int SolvePreconditionedCG( SPDOperator& L , SPDPreconditioner& Pinverse , int iters , int dim , const T* b , T* x , TDotT dot , PreconditionedCGScratch< T >* scratch=NULL , double eps=1e-8 , int threads=1 , bool verbose=false )
	{
		eps *= eps;
		T *r , *d , *q , *s;
		if( scratch ) r = scratch->r , d = scratch->d , q = scratch->q , s = scratch->s;
		else          r = new T[dim] , d = new T[dim] , q = new T[dim] , s = new T[dim];
		memset( r , 0 , sizeof(T)*dim ) , memset( d , 0 , sizeof(T)*dim ) , memset( q , 0 , sizeof(T)*dim ) , memset( s , 0 , sizeof(T)*dim );
		double delta_new = 0 , delta_0;

		L( x , r );
		ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int , size_t i ){ r[i] = b[i] - r[i]; } , threads );
		Pinverse( r , d );
		std::vector< double > _delta_news( ThreadPool::NumThreads() , 0 );
		ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int t , size_t i ){ _delta_news[t] += dot( r[i] , d[i] ); } , threads );
		for( unsigned int t=0 ; t<_delta_news.size() ; t++ ) delta_new += _delta_news[t];

		delta_0 = delta_new;
		if( delta_new<=eps )
		{
			if( !scratch ) delete[] r , delete[] d , delete[] q;
			return 0;
		}
		int ii;
		for( ii=0 ; ii<iters && delta_new>eps*delta_0 ; ii++ )
		{
			L( d , q );
			double dDotQ = 0;
			std::vector< double > _dDotQs( ThreadPool::NumThreads() , 0 );
			ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int t , size_t i ){ _dDotQs[t] += dot( d[i] , q[i] ); } , threads );
			for( unsigned int t=0 ; t<_dDotQs.size() ; t++ ) dDotQ += _dDotQs[t];
			if( !dDotQ ) break;
			Real alpha = Real( delta_new / dDotQ );

			const int RESET_COUNT = 50;
			ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int , size_t i ){ x[i] += d[i] * alpha; } , threads );
			if( (ii%RESET_COUNT)==(RESET_COUNT-1) )
			{
				L( x , r );
				ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int , size_t i ){ r[i] = b[i] - r[i]; } , threads );
			}
			else
				ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int , size_t i ){ r[i] -= q[i] * alpha; } , threads );
			Pinverse( r , s );

			double delta_old = delta_new;
			delta_new = 0;
			std::vector< double > _delta_news( ThreadPool::NumThreads() , 0 );
			ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int t , size_t i ){ _delta_news[t] += dot( r[i] , s[i] );} , threads );
			for( unsigned int t=0 ; t<_delta_news.size() ; t++ ) delta_new += _delta_news[t];

			Real beta = Real( delta_new / delta_old );
			ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int , size_t i ){ d[i] = s[i] + d[i] * beta; } , threads );
		}
		if( verbose )
		{
			L( x , r );
			ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int , size_t i ){ r[i] -= b[i]; } , threads );
			printf( "PCCG: %d %g -> %g\n" , ii , SquareNorm( b , dim ) , SquareNorm( r , dim ) );
		}
		if( !scratch ) delete[] r , delete[] d , delete[] q , delete[] s;
		return ii;
	}


#ifdef USE_EIGEN
#include <Eigen/Sparse>
#define STORE_EIGEN_MATRIX
#ifdef EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#endif // EIGEN_USE_MKL_ALL

	template< class Real , class MatrixRowIterator >
	struct EigenSolver
	{
		virtual void update( const SparseMatrixInterface< Real , MatrixRowIterator >& M ) = 0;
		virtual void solve( ConstPointer( Real ) b , Pointer( Real ) x ) = 0;
		virtual size_t dimension( void ) const = 0;
	};

	template< class Real , class MatrixRowIterator >
	class EigenSolverCholeskyLLt : public EigenSolver< Real , MatrixRowIterator >
	{
#ifdef EIGEN_USE_MKL_ALL
		typedef Eigen::PardisoLLT< Eigen::SparseMatrix< double > > Eigen_Solver;
		typedef Eigen::VectorXd                                    Eigen_Vector;
#else // !EIGEN_USE_MKL_ALL
		typedef Eigen::SimplicialLLT< Eigen::SparseMatrix< double > > Eigen_Solver;
		typedef Eigen::VectorXd                                       Eigen_Vector;
#endif // EIGEN_USE_MKL_ALL
		Eigen_Solver _solver;
		Eigen_Vector _eigenB;
#ifdef STORE_EIGEN_MATRIX
		Eigen::SparseMatrix< double > _eigenM;
#endif // STORE_EIGEN_MATRIX
	public:
		EigenSolverCholeskyLLt( const SparseMatrixInterface< Real , MatrixRowIterator >& M , bool analyzeOnly=false )
		{
#ifdef STORE_EIGEN_MATRIX
			_eigenM.resize( int( M.Rows() ) , int( M.Rows() ) );
#else // !STORE_EIGEN_MATRIX
			Eigen::SparseMatrix< double > eigenM( int( M.Rows() ) , int( M.Rows() ) );
#endif // STORE_EIGEN_MATRIX
			std::vector< Eigen::Triplet< double > > triplets;
			triplets.reserve( M.Entries() );
			for( int i=0 ; i<M.Rows() ; i++ ) for( MatrixRowIterator iter=M.begin(i) ; iter!=M.end(i) ; iter++ ) triplets.push_back( Eigen::Triplet< double >( i , iter->N , iter->Value ) );
#ifdef STORE_EIGEN_MATRIX
			_eigenM.setFromTriplets( triplets.begin() , triplets.end() );
			_solver.analyzePattern( _eigenM );
#else // !STORE_EIGEN_MATRIX
			eigenM.setFromTriplets( triplets.begin() , triplets.end() );
			_solver.analyzePattern( eigenM );
#endif // STORE_EIGEN_MATRIX
			if( !analyzeOnly )
			{
#ifdef STORE_EIGEN_MATRIX
				_solver.factorize( _eigenM );
#else // !STORE_EIGEN_MATRIX
				_solver.factorize( eigenM );
#endif // STORE_EIGEN_MATRIX
				if( _solver.info()!=Eigen::Success ) fprintf( stderr , "[ERROR] EigenSolverCholeskyLLt::EigenSolverCholeskyLLt Failed to factorize matrix\n" ) , exit(0);
			}
			_eigenB.resize( M.Rows() );
		}
		void update( const SparseMatrixInterface< Real , MatrixRowIterator >& M )
		{
#ifdef STORE_EIGEN_MATRIX
			ThreadPool::ParallelFor( 0 , M.Rows() , [&]( unsigned int , size_t i ){ for( MatrixRowIterator iter=M.begin(i) ; iter!=M.end(i) ; iter++ ) _eigenM.coeffRef( i , iter->N ) = iter->Value; } );
			_solver.factorize( _eigenM );
#else // !STORE_EIGEN_MATRIX
			Eigen::SparseMatrix< double > eigenM( int( M.Rows() ) , int( M.Rows() ) );
			std::vector< Eigen::Triplet< double > > triplets;
			triplets.reserve( M.Entries() );
			for( int i=0 ; i<M.Rows() ; i++ ) for( MatrixRowIterator iter=M.begin(i) ; iter!=M.end(i) ; iter++ ) triplets.push_back( Eigen::Triplet< double >( i , iter->N , iter->Value ) );
			eigenM.setFromTriplets( triplets.begin() , triplets.end() );
			_solver.factorize( eigenM );
#endif // STORE_EIGEN_MATRIX
			switch( _solver.info() )
			{
			case Eigen::Success: break;
			case Eigen::NumericalIssue: fprintf( stderr , "[ERROR] EigenSolverCholeskyLLt::update Failed to factorize matrix (numerical issue)\n" ) , exit(0);
			case Eigen::NoConvergence:  fprintf( stderr , "[ERROR] EigenSolverCholeskyLLt::update Failed to factorize matrix (no convergence)\n" ) , exit(0);
			case Eigen::InvalidInput:   fprintf( stderr , "[ERROR] EigenSolverCholeskyLLt::update Failed to factorize matrix (invalid input)\n" ) , exit(0);
			default: fprintf( stderr , "[ERROR] EigenSolverCholeskyLLt::update Failed to factorize matrix\n" ) , exit(0);
			}
		}
		void solve( const Eigen_Vector &b , Eigen_Vector &x ){ x = _solver.solve( b ); }
		void solve( ConstPointer( Real ) b , Pointer( Real ) x )
		{
			ThreadPool::ParallelFor( 0 , _eigenB.size() , [&]( unsigned int , size_t i ){ _eigenB[i] = b[i]; } );
			Eigen_Vector eigenX = _solver.solve( _eigenB );
			ThreadPool::ParallelFor( 0 , eigenX.size() , [&]( unsigned int , size_t i ){ x[i] = (Real)eigenX[i]; } );
		}
		size_t dimension( void ) const { return _eigenB.size(); }
		static void Solve( const SparseMatrixInterface< Real , MatrixRowIterator >& M , ConstPointer( Real ) b , Pointer( Real ) x ){ EigenSolverCholeskyLLt solver( M ) ; solver.solve( b , x ); }
	};

	template< class Real , class MatrixRowIterator >
	class EigenSolverCholeskyLDLt : public EigenSolver< Real , MatrixRowIterator >
	{
#ifdef EIGEN_USE_MKL_ALL
		typedef Eigen::PardisoLDLT< Eigen::SparseMatrix< double > > Eigen_Solver;
		typedef Eigen::VectorXd                                     Eigen_Vector;
#else // !EIGEN_USE_MKL_ALL
		typedef Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > > Eigen_Solver;
		typedef Eigen::VectorXd                                        Eigen_Vector;
#endif // EIGEN_USE_MKL_ALL
		Eigen_Solver _solver;
		Eigen_Vector _eigenB;
	public:
		EigenSolverCholeskyLDLt( const SparseMatrixInterface< Real , MatrixRowIterator >& M , bool analyzeOnly=false )
		{
			Eigen::SparseMatrix< double > eigenM( int( M.Rows() ) , int( M.Rows() ) );
			std::vector< Eigen::Triplet<double> > triplets;
			triplets.reserve( M.Entries() );
			for( int i=0 ; i<M.Rows() ; i++ ) for( MatrixRowIterator iter=M.begin(i) ; iter!=M.end(i) ; iter++ ) triplets.push_back( Eigen::Triplet< double >( i , iter->N , iter->Value ) );
			eigenM.setFromTriplets( triplets.begin() , triplets.end() );
			_solver.analyzePattern( eigenM );
			if( !analyzeOnly )
			{
				_solver.factorize( eigenM );
				if( _solver.info()!=Eigen::Success ) fprintf( stderr , "[ERROR] EigenSolverCholeskyLDLt::EigenSolverCholeskyLDLt Failed to factorize matrix\n" ) , exit(0);
			}
			_eigenB.resize( M.Rows() );
		}
		void update( const SparseMatrixInterface< Real , MatrixRowIterator >& M )
		{
			Eigen::SparseMatrix< double > eigenM( int( M.Rows() ) , int( M.Rows() ) );
			std::vector< Eigen::Triplet<double> > triplets;
			triplets.reserve( M.Entries() );
			for( int i=0 ; i<M.Rows() ; i++ ) for( MatrixRowIterator iter=M.begin(i) ; iter!=M.end(i) ; iter++ ) triplets.push_back( Eigen::Triplet< double >( i , iter->N , iter->Value ) );
			eigenM.setFromTriplets( triplets.begin() , triplets.end() );
			_solver.factorize( eigenM );
			if( _solver.info()!=Eigen::Success ) fprintf( stderr , "[ERROR] EigenSolverCholeskyLDLt::update Failed to factorize matrix\n" ) , exit(0);
		}
		void solve( const Eigen_Vector &b , Eigen_Vector &x ){ x = _solver.solve( b ); }
		void solve( ConstPointer( Real ) b , Pointer( Real ) x )
		{
			ThreadPool::ParallelFor( 0 , _eigenB.size() , [&]( unsigned int , size_t i ){ _eigenB[i] = b[i]; } );
			Eigen_Vector eigenX = _solver.solve( _eigenB );
			ThreadPool::ParallelFor( 0 , eigenX.size() , [&]( unsigned int , size_t i ){ x[i] = (Real)eigenX[i]; } );
		}
		size_t dimension( void ) const { return _eigenB.size(); }
		static void Solve( const SparseMatrixInterface< Real , MatrixRowIterator >& M , ConstPointer( Real ) b , Pointer( Real ) x ){ EigenSolverCholeskyLDLt solver( M ) ; solver.solve( b , x ); }
	};

	template< class Real , class MatrixRowIterator >
	class EigenSolverCG : public EigenSolver< Real , MatrixRowIterator >
	{
#if 1
		//	Eigen::ConjugateGradient< Eigen::SparseMatrix< double > , Eigen::Lower , Eigen::IncompleteLUT< double > > _solver;
		Eigen::ConjugateGradient< Eigen::SparseMatrix< double > > _solver;
#else
		Eigen::BiCGSTAB< Eigen::SparseMatrix< double > > _solver;
#endif
		Eigen::VectorXd _eigenB , _eigenX;
		Eigen::SparseMatrix< double > _eigenM;
	public:
		EigenSolverCG( const SparseMatrixInterface< Real , MatrixRowIterator >& M , int iters=20 )
		{
			_eigenM.resize( (int)M.Rows() , (int)M.Rows() );
			std::vector< Eigen::Triplet< double > > triplets;
			triplets.reserve( M.Entries() );
			for( int i=0 ; i<M.Rows() ; i++ ) for( MatrixRowIterator iter=M.begin(i) ; iter!=M.end(i) ; iter++ ) triplets.push_back( Eigen::Triplet< double >( i , iter->N , iter->Value ) );
			_eigenM.setFromTriplets( triplets.begin() , triplets.end() );
			_solver.compute( _eigenM );
			_solver.analyzePattern( _eigenM );
			if( _solver.info()!=Eigen::Success ) fprintf( stderr , "[ERROR] EigenSolverCG::EigenSolverCG Failed to factorize matrix\n" ) , exit(0);
			_eigenB.resize( M.Rows() ) , _eigenX.resize( M.Rows() );
			_solver.setMaxIterations( iters );
		}
		void update( const SparseMatrixInterface< Real , MatrixRowIterator >& M )
		{
			ThreadPool::ParallelFor( 0 , M.Rows() , [&]( unsigned int , size_t i ){ for( MatrixRowIterator iter=M.begin(i) ; iter!=M.end(i) ; iter++ ) _eigenM.coeffRef( i , iter->N ) = iter->Value; } );
			_solver.compute( _eigenM );
			_solver.analyzePattern( _eigenM );
			if( _solver.info()!=Eigen::Success ) fprintf( stderr , "[ERROR] EigenSolverCG::update Failed to factorize matrix\n" ) , exit(0);
		}

		void setIters( int iters ){ _solver.setMaxIterations( iters ); }
		void solve( const Real* b , Real* x )
		{
			ThreadPool::ParallelFor( 0 , _eigenB.size() , [&]( unsigned int , size_t i ){ _eigenB[i] = b[i] , _eigenX[i] = x[i]; } );
			for( int i=0 ; i<_eigenB.size() ; i++ ) _eigenB[i] = b[i] , _eigenX[i] = x[i];
			_eigenX = _solver.solveWithGuess( _eigenB , _eigenX );
			ThreadPool::ParallelFor( 0 , _eigenX.size() , [&]( unsigned int , size_t i ){ x[i] = _eigenX[i]; } );
		}
		size_t dimension( void ) const { return _eigenB.size(); }
		static void Solve( const SparseMatrixInterface< Real , MatrixRowIterator >& M , const Real* b , Real* x , int iters ){ EigenSolverCG solver( M , iters ) ; solver.solve( b , x ); }
	};

#endif // USE_EIGEN

#ifdef USE_CHOLMOD
	class CholmodSolver
	{
		const static bool LOWER_TRIANGULAR = true;
		int dim;
		cholmod_factor* cholmod_L;
		cholmod_dense*  cholmod_b;
		cholmod_sparse* cholmod_M;
		std::vector< bool > flaggedValues;
		template< class Real , class MatrixRowIterator > void _init( const SparseMatrixInterface< Real , MatrixRowIterator >& M );
	public:
		static cholmod_common cholmod_C;
		static bool cholmod_C_set;

		template< class Real , class MatrixRowIterator >
		CholmodSolver( const SparseMatrixInterface< Real , MatrixRowIterator >& M , bool analyzeOnly=false );
		~CholmodSolver( void );

		template< class Real > void solve( ConstPointer( Real ) b , Pointer( Real ) x );
		template< class Real , class MatrixRowIterator > void update( const SparseMatrixInterface< Real , MatrixRowIterator >& M );
		int nonZeros( void ) const;

	};
	bool CholmodSolver::cholmod_C_set = false;
	cholmod_common CholmodSolver::cholmod_C;

	template< class Real , class MatrixRowIterator > CholmodSolver::CholmodSolver( const SparseMatrixInterface< Real , MatrixRowIterator >& M , bool analyzeOnly ){ _init( M ) ; if( !analyzeOnly ) update( M ); }
	template< class Real , class MatrixRowIterator >
	void CholmodSolver::_init( const SparseMatrixInterface< Real , MatrixRowIterator >& M )
	{
		{
			if( !cholmod_C_set ) CHOLMOD(start)( &cholmod_C );
			cholmod_C_set = true;
		}
		dim = (int)M.Rows();

		int maxEntries;
		if( LOWER_TRIANGULAR )
		{
			maxEntries = (int)( ( M.Entries()-M.Rows() ) / 2 + M.Rows() );
			cholmod_M = CHOLMOD(allocate_sparse)( dim , dim , maxEntries , 0 , 1 , -1 , CHOLMOD_REAL , &cholmod_C );
		}
		else
		{
			maxEntries = (int)M.Entries();
			cholmod_M = CHOLMOD(allocate_sparse)( dim , dim , maxEntries , 0 , 1 ,  0 , CHOLMOD_REAL , &cholmod_C );
		}
		cholmod_M->i = malloc( sizeof( SOLVER_LONG ) * maxEntries );
		cholmod_M->x = malloc( sizeof( double ) * maxEntries );

		SOLVER_LONG *_p = (SOLVER_LONG*)cholmod_M->p;
		SOLVER_LONG *_i = (SOLVER_LONG*)cholmod_M->i;

		int off = 0;
		dim = 0;

		for( int i=0 ; i<M.Rows() ; i++ )
		{
			_p[dim++] = off;
			for( MatrixRowIterator iter=M.begin(i) ; iter!=M.end(i) ; iter++ ) if( !LOWER_TRIANGULAR || iter->N>=i ) _i[off++] = iter->N;
		}
		_p[dim] = off;

		cholmod_L = CHOLMOD(analyze)( cholmod_M , &cholmod_C );
		cholmod_b = CHOLMOD(allocate_dense)( dim , 1 , dim , cholmod_M->xtype , &cholmod_C );
	}
	template< class Real , class MatrixRowIterator >
	void CholmodSolver::update( const SparseMatrixInterface< Real , MatrixRowIterator >& M )
	{
		double *_x = (double*)cholmod_M->x;
		int off = 0;

		SOLVER_LONG *_p = (SOLVER_LONG*)cholmod_M->p;
		ThreadPool::ParallelFor
		(
			0 , M.Rows() ,
			[&]( unsigned int , size_t i )
			{
				int off = (int)_p[i];
				for( MatrixRowIterator iter=M.begin(i) ; iter!=M.end(i) ; iter++ )if( !LOWER_TRIANGULAR || iter->N>=i ) _x[off++] = double( iter->Value );
			}
		);

		cholmod_C.print = 0;
		CHOLMOD(factorize)( cholmod_M , cholmod_L , &cholmod_C );
		if     ( cholmod_C.status==CHOLMOD_NOT_POSDEF    ) fprintf( stderr , "[ERROR] CholmodSolver::update: Matrix not positive-definite\n" )                   , exit( 0 );
		else if( cholmod_C.status==CHOLMOD_OUT_OF_MEMORY ) fprintf( stderr , "[WARNING] CholmodSolver::update: CHOLMOD ran out of memory\n" )                    , exit( 0 );
		else if( cholmod_C.status!=CHOLMOD_OK            ) fprintf( stderr , "[WARNING] CholmodSolver::update: CHOLMOD status not OK: %d\n" , cholmod_C.status ) , exit( 0 );
	}
	CholmodSolver::~CholmodSolver( void )
	{
		if( cholmod_L ) CHOLMOD(free_factor)( &cholmod_L , &cholmod_C ) , cholmod_L = NULL;
		if( cholmod_b ) CHOLMOD(free_dense )( &cholmod_b , &cholmod_C ) , cholmod_b = NULL;
		if( cholmod_M ) CHOLMOD(free_sparse)( &cholmod_M , &cholmod_C ) , cholmod_M = NULL;
	}

	template< class Real >
	void CholmodSolver::solve( ConstPointer( Real ) b , Pointer( Real ) x )
	{
		double* _b = (double*)cholmod_b->x;
		for( int i=0 ; i<dim ; i++ ) _b[i] = (double)b[i];

		cholmod_dense* cholmod_x = CHOLMOD(solve)( CHOLMOD_A , cholmod_L , cholmod_b , &cholmod_C );
		double* _x = (double*)cholmod_x->x;
		for( int i=0 ; i<dim ; i++ ) x[i] = (Real)_x[i];

		CHOLMOD(free_dense)( &cholmod_x , &cholmod_C );
	}
	int CholmodSolver::nonZeros( void ) const
	{
		long long nz = 0;
		if( cholmod_L->xtype != CHOLMOD_PATTERN && !(cholmod_L->is_super ) ) for( int i=0 ; i<cholmod_L->n ; i++ ) nz += ((SOLVER_LONG*)cholmod_L->nz)[i];
		bool examine_super = false;
		if( cholmod_L->xtype != CHOLMOD_PATTERN ) examine_super = true ;
		else                                      examine_super = ( ((int*)cholmod_L->s)[0] != (-1));
		if( examine_super )
		{
			/* check and print each supernode */
			for (int s = 0 ; s < cholmod_L->nsuper ; s++)
			{
				int k1 = ((int*)cholmod_L->super) [s] ;
				int k2 = ((int*)cholmod_L->super) [s+1] ;
				int psi = ((int*)cholmod_L->pi)[s] ;
				int psend = ((int*)cholmod_L->pi)[s+1] ;
				int nsrow = psend - psi ;
				int nscol = k2 - k1 ;
				nz += nscol * nsrow - (nscol*nscol - nscol)/2 ;
			}
		}
		return (int)nz;
	}
#endif // USE_CHOLMOD


#ifdef USE_PARDISO
	struct PardisoSolver
	{
	protected:
		int _n;
		Pointer( int ) _ia;
		Pointer( int ) _ja;
		Pointer( int ) _nnz;
		Pointer( double ) _a;
		Pointer( double ) _b;
		Pointer( double ) _x;
		int _mtype = -2;
		int _nrhs = 1;
		void *_pt[64];
		int _iparm [64];
		double _dparm [64];
		int _maxfct , _mnum , _phase , _error , _msglvl , _solver;
		int _num_procs;
		char *_var;
		double _ddum;
		int _idum;

		template< class Real , class MatrixRowIterator > void _init( const SparseMatrixInterface< Real , MatrixRowIterator >& M );
	public:

		template< class Real , class MatrixRowIterator >
		PardisoSolver( const SparseMatrixInterface< Real , MatrixRowIterator >& M , bool analyzeOnly=false );
		~PardisoSolver( void );

		template< class Real > void solve( ConstPointer( Real ) b , Pointer( Real ) x );
		template< class Real , class MatrixRowIterator > void update( const SparseMatrixInterface< Real , MatrixRowIterator >& M );
	};

	template< class Real , class MatrixRowIterator >
	PardisoSolver::PardisoSolver( const SparseMatrixInterface< Real , MatrixRowIterator >& M , bool analyzeOnly ){ _init( M ) ; if( !analyzeOnly ) update( M ); }
	PardisoSolver::~PardisoSolver( void )
	{
		FreePointer( _ia );
		FreePointer( _ja );
		FreePointer( _nnz );
		FreePointer( _a );
		FreePointer( _b );
		FreePointer( _x );
	}

	template< class Real , class MatrixRowIterator >
	void PardisoSolver::_init( const SparseMatrixInterface< Real , MatrixRowIterator >& M )
	{
		// Allocate memory
		_n = (int)M.Rows();
		int entries = (int)( ( M.Entries()-M.Rows() ) / 2 + M.Rows() );

		_ia = AllocPointer< int >( _n+1 );
		_ja = AllocPointer< int >( entries );
		_a = AllocPointer< double >( entries );
		_b = AllocPointer< double >( _n );
		_x = AllocPointer< double >( _n );
		{
			int off=0 , dim=0;
			for( int i=0 ; i<M.Rows() ; i++ )
			{
				_ia[dim++] = off+1;
				for( MatrixRowIterator iter=M.begin(i) ; iter!=M.end(i) ; iter++ ) if( iter->N>=i ) _ja[off++] = iter->N+1;
			}
			_ia[dim] = off+1;
		}

		_error = 0;
		_solver = 0;
		pardisoinit_d( _pt , &_mtype , &_solver , _iparm , _dparm , &_error );
		if( _error )
		{
			fprintf( stderr , "[ERROR] ParadisoSolver::_init\n" );
			if ( _error==-10 ) fprintf( stderr , "\tNo license file found\n" );
			if ( _error==-11 ) fprintf( stderr , "\tLicense is expired\n" );
			if ( _error==-12 ) fprintf( stderr , "\tWrong username or hostname\n" );
			exit( 0 );
		}
		else printf ( "[ PARDISO ]: License check was successful ...\n" );

		_var = getenv ( " OMP_NUM_THREADS " );
		if( _var ) sscanf( _var , " %d " , &_num_procs );
		else fprintf ( stderr , "[ERROR] PardisoSolver::_init: Set environment OMP_NUM_THREADS to 1\n" ) , exit(0);
		_iparm [2] = _num_procs ;
		_maxfct = 1;
		_mnum = 1;
		_msglvl = 1;
		_error = 0;

		pardiso_chkmatrix_d( &_mtype , &_n , _a , _ia , _ja , &_error );
		if( _error ) fprintf( stderr , "[ERROR] PardisoSolver::_init: in consistency of matrix: %d\n" , _error ) , exit( 0 );

		pardiso_printstats_d( &_mtype , &_n , _a , _ia , _ja , &_nrhs , _b , &_error );
		if( _error ) fprintf( stderr , "[ERROR] PardisoSolver::_init: stats: % d\n" , _error ) , exit( 0 );

		_phase = 11;
		_iparm [10] = 0;
		pardiso_d( _pt , &_maxfct , &_mnum , &_mtype , &_phase , &_n , _a , _ia , _ja , &_idum , &_nrhs , _iparm , &_msglvl , &_ddum , &_ddum , &_error , _dparm );
		if( _error ) fprintf( stderr , "[ERROR] PardisoSolver::_init: during symbolic factorization : %d\n" , _error ) , exit( 0 );
		printf ( "Reordering completed ...\n" );
		printf ( "Number of nonzeros in factors = %d\n" , _iparm [17]);
		printf ( "Number of factorization MFLOPS = %d\n" , _iparm [18]);

	}
	template< class Real , class MatrixRowIterator >
	void PardisoSolver::update( const SparseMatrixInterface< Real , MatrixRowIterator >& M )
	{
		ThreadPool::ParallelFor
		(
			0 , M.Rows() ,
			[&]( unsigned int , size_t i )
			{
				int off = (int)_ia[i];
				for( MatrixRowIterator iter=M.begin(i) ; iter!=M.end(i) ; iter++ )if( iter->N>=i ) _a[off++] = double( iter->Value );
			}
		);

		_phase = 22;
		_iparm[32] = 1;
		pardiso_d( _pt , &_maxfct , &_mnum , &_mtype , &_phase , &_n , _a , _ia , _ja , &_idum , &_nrhs , _iparm , &_msglvl , &_ddum , &_ddum , &_error , _dparm );
		if( _error ) fprintf( stderr , "[ERROR] PardisoSolver::update: during numerical factorization : %d\n" , _error ) , exit( 0 );
	}

	template< class Real >
	void PardisoSolver::solve( ConstPointer( Real ) b , Pointer( Real ) x )
	{
		ThreadPool::ParallelFor( 0 , _n , [&]( unsigned int , size_t i ){ _b[i] = (double)b[i]; } );

		pardiso_chkvec_d( &_n , &_nrhs , _b , &_error );
		if( _error ) fprintf( stderr , "[ERROR] PardisoSolver::solve: right hand side : %d\n" , _error ) , exit( 0 );

		_phase = 33;
		_iparm [7] = 1;
		pardiso_d( _pt , &_maxfct , &_mnum , &_mtype , &_phase , &_n , _a , _ia , _ja , &_idum , &_nrhs , _iparm , &_msglvl , _b , _x , &_error , _dparm );
		if( _error ) fprintf( stderr , "[ERROR] PardisoSolver::solve: during solution : %d\n" , _error ) , exit( 0 );

		ThreadPool::ParallelFor( 0 , _n , [&]( unsigned int , size_t i ){ x[i] = (Real)_x[i]; } );
	}
#endif // USE_PARDISO
}
#endif // LINEAR_SOLVERS_INCLUDE