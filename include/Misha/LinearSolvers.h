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
#ifndef LINEAR_SOLVERS_INCLUDE
#define LINEAR_SOLVERS_INCLUDE
#ifndef USE_EIGEN
#define USE_EIGEN 1
#define NEW_EIGEN
#endif // USE_EIGEN


#ifdef USE_CHOLMOD
#include <Cholmod/cholmod.h>
#pragma comment( lib , "CHOLMOD_FULL.lib" )
#elif defined(EIGEN_USE_MKL_ALL)
#pragma comment( lib , "mkl_core.lib" )
#pragma comment( lib , "mkl_intel_lp64.lib" )
#pragma comment( lib , "mkl_intel_thread.lib" )
#pragma comment( lib , "mkl_blas95_lp64.lib
#pragma comment( lib , "libiomp5md.lib" )
#endif // USE_CHOLMOD

#include <Misha/SparseMatrixInterface.h>

inline double                        SquareNorm(const double* values, int dim) { double norm2 = 0; for (int i = 0; i<dim; i++) norm2 += values[i] * values[i]; return norm2; }
inline double                        SquareNorm(const  float* values, int dim) { double norm2 = 0; for (int i = 0; i<dim; i++) norm2 += values[i] * values[i]; return norm2; }
template< class Type > inline double SquareNorm(const   Type* values, int dim) { double norm2 = 0; for (int i = 0; i<dim; i++) norm2 += values[dim].squareNorm(); return norm2; }

inline double                        SquareDifference(const double* values1, const double* values2, int dim) { double norm2 = 0; for (int i = 0; i<dim; i++) norm2 += (values1[i] - values2[i]) * (values1[i] - values2[i]); return norm2; }
inline double                        SquareDifference(const  float* values1, const  float* values2, int dim) { double norm2 = 0; for (int i = 0; i<dim; i++) norm2 += (values1[i] - values2[i]) * (values1[i] - values2[i]); return norm2; }
template< class Type > inline double SquareDifference(const   Type* values1, const   Type* values2, int dim) { double norm2 = 0; for (int i = 0; i<dim; i++) norm2 += (values1[dim] - values2[dim]).squareNorm(); return norm2; }


// This is the conjugate gradients solver.
// The assumption is that the class SPDOperator defines a method operator()( const Real* , Real* ) which corresponds to applying a symmetric positive-definite operator.
template< typename T >
struct CGScratch
{
	T *r, *d, *q;
	CGScratch(void) : r(NULL), d(NULL), q(NULL), _dim(0) { ; }
	CGScratch(int dim) : r(NULL), d(NULL), q(NULL) { resize(dim); }
	~CGScratch(void) { resize(0); }
	void resize(int dim)
	{
		if (dim != _dim)
		{
			if (r) delete[] r; r = NULL;
			if (d) delete[] d; d = NULL;
			if (q) delete[] q; q = NULL;
			if (dim) r = new T[dim], d = new T[dim], q = new T[dim];
			_dim = dim;
		}
	}
protected:
	int _dim;
};
template< class T >
struct PreconditionedCGScratch : public CGScratch< T >
{
	T *s;
	PreconditionedCGScratch(void) : CGScratch< T >(), s(NULL) { ; }
	PreconditionedCGScratch(int dim) : CGScratch< T >() { resize(dim); }
	~PreconditionedCGScratch(void) { resize(0); }
	void resize(int dim)
	{
		if (dim != CGScratch< T >::_dim)
		{
			if (s) delete[] s; s = NULL;
			if (dim) s = new T[dim];
		}
		CGScratch< T >::resize(dim);
	}
};
template< class Real, typename T >
struct DiagonalPreconditioner
{
	Real* iDiagonal;
	DiagonalPreconditioner(void) : iDiagonal(NULL), _dim(0) { ; }
	~DiagonalPreconditioner(void) { if (iDiagonal) delete[] iDiagonal; iDiagonal = NULL; }
	template< class MatrixRowIterator >
	void set(const SparseMatrixInterface< Real, MatrixRowIterator >& M)
	{
		if (_dim != M.Rows())
		{
			_dim = (int)M.Rows();
			if (iDiagonal) delete[] iDiagonal, iDiagonal = NULL;
			if (_dim>0) iDiagonal = new Real[_dim];
		}
		memset(iDiagonal, 0, sizeof(Real)*_dim);
#pragma omp parallel for
		for (int i = 0; i<M.Rows(); i++)
		{
			for (MatrixRowIterator iter = M.begin(i); iter != M.end(i); iter++) if (iter->N == i) iDiagonal[i] += iter->Value;
			iDiagonal[i] = (Real)1. / iDiagonal[i];
		}
	}
	void operator()(const T* in, T* out) const
	{
#pragma omp parallel for
		for (int i = 0; i<_dim; i++) out[i] = in[i] * iDiagonal[i];
	}
protected:
	int _dim;
};

template< class Real, typename T, class SPDOperator, typename TDotT >
int SolveCG(SPDOperator& L, int iters, int dim, const T* b, T* x, TDotT dot, CGScratch< T >* scratch = NULL, double eps = 1e-8, bool verbose = false)
{
	eps *= eps;
	T *r, *d, *q;
	if (scratch) r = scratch->r, d = scratch->d, q = scratch->q;
	else          r = new T[dim], d = new T[dim], q = new T[dim];
	memset(r, 0, sizeof(T)*dim), memset(d, 0, sizeof(T)*dim), memset(q, 0, sizeof(T)*dim);
	double delta_new = 0, delta_0;

	L(x, r);
#pragma omp parallel for reduction( + : delta_new )
	for( int i=0 ; i<dim ; i++ ){ d[i] = r[i] = b[i] - r[i] ; delta_new += dot(r[i], r[i]); }

	delta_0 = delta_new;
	if (delta_new<eps)
	{
		if (!scratch) delete[] r, delete[] d, delete[] q;
		return 0;
	}

	int ii;
	for (ii = 0; ii<iters && delta_new>eps*delta_0; ii++)
	{
		L(d, q);
		double dDotQ = 0;
#pragma message( "[WARNING] Why do I have to comment this out?" )
#pragma omp parallel for reduction( + : dDotQ )
		for (int i = 0; i<dim; i++) dDotQ += dot(d[i], q[i]);
		Real alpha = Real(delta_new / dDotQ);

		double delta_old = delta_new;
		delta_new = 0;

		const int RESET_COUNT = 50;
		if ((ii%RESET_COUNT) == (RESET_COUNT - 1))
		{
#pragma omp parallel for
			for (int i = 0; i<dim; i++) x[i] += d[i] * alpha;
			L(x, r);
#pragma omp parallel for reduction ( + : delta_new )
			for (int i = 0; i<dim; i++) r[i] = b[i] - r[i], delta_new += dot(r[i], r[i]);
		}
		else
#pragma omp parallel for reduction( + : delta_new )
			for (int i = 0; i<dim; i++) r[i] -= q[i] * alpha, delta_new += dot(r[i], r[i]), x[i] += d[i] * alpha;

		Real beta = Real(delta_new / delta_old);
#pragma omp parallel for
		for (int i = 0; i<dim; i++) d[i] = r[i] + d[i] * beta;
	}
	if (verbose)
	{
		L(x, r);
#pragma omp parallel for
		for (int i = 0; i<dim; i++) r[i] -= b[i];
		printf("CG: %d %g -> %g\n", ii, SquareNorm(b, dim), SquareNorm(r, dim));
	}
	if (!scratch) delete[] r, delete[] d, delete[] q;
	return ii;
}


template< class Real, typename T, class SPDOperator, class SPDPreconditioner>
int SolvePreconditionedCG(SPDOperator& L, SPDPreconditioner& Pinverse, int iters, int dim, const T* b, T* x, PreconditionedCGScratch< T >* scratch = NULL, double eps = 1e-8, bool verbose = false)
{
	int threads = omp_get_num_threads();
	eps *= eps;
	T *r, *d, *q, *s;
	if (scratch) r = scratch->r, d = scratch->d, q = scratch->q, s = scratch->s;
	else          r = new T[dim], d = new T[dim], q = new T[dim], s = new T[dim];
	memset(r, 0, sizeof(T)*dim), memset(d, 0, sizeof(T)*dim), memset(q, 0, sizeof(T)*dim), memset(s, 0, sizeof(T)*dim);
	double delta_new = 0, delta_0;

	L(x, r);
#pragma omp parallel for
	for (int i = 0; i<dim; i++) r[i] = b[i] - r[i];
	Pinverse(r, d);

#pragma omp parallel for reduction( + : delta_new )
	for (int i = 0; i<dim; i++) { double _delta_new = DotData(r[i], d[i]); delta_new += _delta_new; }

	delta_0 = delta_new;
	if (delta_new<eps)
	{
		if (!scratch) delete[] r, delete[] d, delete[] q;
		return 0;
	}
	int ii;
	for (ii = 0; ii<iters && delta_new>eps*delta_0; ii++)
	{
		L(d, q);
		double dDotQ = 0;

#pragma omp parallel for reduction( + : dDotQ )
		for (int i = 0; i<dim; i++) { double _dDotQ = DotData(d[i], q[i]); dDotQ += _dDotQ; }

		Real alpha = Real(delta_new / dDotQ);

		const int RESET_COUNT = 50;
#pragma omp parallel for
		for (int i = 0; i<dim; i++) x[i] += d[i] * alpha;
		if ((ii%RESET_COUNT) == (RESET_COUNT - 1))
		{
			L(x, r);
#pragma omp parallel for
			for (int i = 0; i<dim; i++) r[i] = b[i] - r[i];
		}
		else
#pragma omp parallel for
			for (int i = 0; i<dim; i++) r[i] -= q[i] * alpha;
		Pinverse(r, s);

		double delta_old = delta_new;
		delta_new = 0;

#if 0
		Vec4f deltas[8];
		for (int t = 0; t<8; t++) deltas[t] = Vec4f(0);
#pragma omp parallel for
		for (int i = 0; i<dim; i++) { Vec4f _delta_new = r[i] * s[i]; deltas[omp_get_thread_num()] += _delta_new; }
		for (int t = 0; t<8; t++) delta_new += horizontal_add(deltas[t]);
#else
#pragma omp parallel for reduction( + : delta_new )
		for (int i = 0; i<dim; i++) { double _delta_new = DotData(r[i], s[i]); delta_new += _delta_new; }
#endif

		Real beta = Real(delta_new / delta_old);
#pragma omp parallel for
		for (int i = 0; i<dim; i++) d[i] = s[i] + d[i] * beta;
	}
	if (!scratch) delete[] r, delete[] d, delete[] q, delete[] s;
	return ii;
}


#if USE_EIGEN
#define STORE_EIGEN_MATRIX
#ifdef EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#else // !EIGEN_USE_MKL_ALL
#include <Eigen/Sparse>
#endif // EIGEN_USE_MKL_ALL

#ifdef NEW_EIGEN
template< class Real, class MatrixRowIterator >
struct EigenSolver
{
	virtual void update(const SparseMatrixInterface< Real, MatrixRowIterator >& M) = 0;
	virtual void solve(ConstPointer(Real) b, Pointer(Real) x) = 0;
	virtual size_t dimension(void) const = 0;
};

template< typename Real, class MatrixRowIterator , typename InternalReal=Real >
class EigenSolverCholeskyLLt : public EigenSolver< Real, MatrixRowIterator >
{
#ifdef EIGEN_USE_MKL_ALL
	typedef Eigen::PardisoLLT< Eigen::SparseMatrix< InternalReal > > Eigen_Solver;
#else // !EIGEN_USE_MKL_ALL
	typedef Eigen::SimplicialLLT< Eigen::SparseMatrix< InternalReal > > Eigen_Solver;
#endif // EIGEN_USE_MKL_ALL
	typedef Eigen::Matrix< InternalReal , Eigen::Dynamic , 1 > Eigen_Vector;
	Eigen_Solver _solver;
	Eigen_Vector _eigenB;
#ifdef STORE_EIGEN_MATRIX
	Eigen::SparseMatrix< InternalReal > _eigenM;
#endif // STORE_EIGEN_MATRIX
public:
	EigenSolverCholeskyLLt( const SparseMatrixInterface< Real , MatrixRowIterator >& M , bool analyzeOnly=false )
	{
#ifdef STORE_EIGEN_MATRIX
		_eigenM.resize( (int)M.Rows() , (int)M.Rows() );
#else // !STORE_EIGEN_MATRIX
		Eigen::SparseMatrix< InternalReal > eigenM( (int)M.Rows() , (int)M.Rows() );
#endif // STORE_EIGEN_MATRIX
		std::vector< Eigen::Triplet< InternalReal > > triplets;
		triplets.reserve( M.Entries() );
		for( int i=0 ; i<M.Rows() ; i++ ) for( MatrixRowIterator iter=M.begin(i) ; iter!=M.end(i) ; iter++ ) triplets.push_back( Eigen::Triplet< InternalReal >( i , iter->N , iter->Value ) );
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
			if( _solver.info()!=Eigen::Success ) fprintf(stderr, "[ERROR] EigenSolverCholeskyLLt::EigenSolverCholeskyLLt Failed to factorize matrix\n" ) , exit(0);
		}
		_eigenB.resize( M.Rows() );
	}
	void update( const SparseMatrixInterface< Real, MatrixRowIterator >& M )
	{
#ifdef STORE_EIGEN_MATRIX
#pragma omp parallel for
		for( int i=0; i<M.Rows() ; i++) for( MatrixRowIterator iter=M.begin(i) ; iter!=M.end(i) ; iter++ ) _eigenM.coeffRef( i , iter->N ) = iter->Value;
		_solver.factorize( _eigenM );
#else // !STORE_EIGEN_MATRIX
		Eigen::SparseMatrix< InternalReal > eigenM( (int)M.Rows() , (int)M.Rows() );
		std::vector< Eigen::Triplet< InternalReal > > triplets;
		triplets.reserve( M.Entries() );
		for( int i=0 ; i<M.Rows() ; i++ ) for( MatrixRowIterator iter=M.begin(i) ; iter!=M.end(i) ; iter++ ) triplets.push_back( Eigen::Triplet< InternalReal >( i , iter->N , iter->Value ) );
		eigenM.setFromTriplets( triplets.begin() , triplets.end() );
		_solver.factorize( eigenM );
#endif // STORE_EIGEN_MATRIX
		switch( _solver.info() )
		{
		case Eigen::Success: break;
		case Eigen::NumericalIssue: fprintf( stderr , "[ERROR] EigenSolverCholeskyLLt::update Failed to factorize matrix (numerical issue)\n" ) , exit(0);
		case Eigen::NoConvergence:  fprintf( stderr , "[ERROR] EigenSolverCholeskyLLt::update Failed to factorize matrix (no convergence)\n" ) , exit(0);
		case Eigen::InvalidInput:   fprintf( stderr , "[ERROR] EigenSolverCholeskyLLt::update Failed to factorize matrix (invalid input)\n" ) , exit(0);
		default:                    fprintf( stderr , "[ERROR] EigenSolverCholeskyLLt::update Failed to factorize matrix\n" ) , exit(0);
		}
	}
	void solve( const Eigen_Vector& b , Eigen_Vector& x ){ x = _solver.solve(b); }
	void solve( ConstPointer(Real) b , Pointer(Real) x )
	{
#pragma omp parallel for
		for( int i=0 ; i<_eigenB.size() ; i++ ) _eigenB[i] = b[i];
		Eigen_Vector eigenX = _solver.solve( _eigenB );
#pragma omp parallel for
		for( int i=0 ; i<eigenX.size() ; i++ ) x[i] = (Real)eigenX[i];
	}
	size_t dimension( void ) const { return _eigenB.size(); }
	static void Solve( const SparseMatrixInterface< Real, MatrixRowIterator >& M , ConstPointer(Real) b , Pointer(Real) x ) { EigenSolverCholeskyLLt solver(M) ; solver.solve(b,x); }
};
template< typename Real, class MatrixRowIterator , typename InternalReal=Real >
class EigenSolverCholeskyLDLt : public EigenSolver< Real, MatrixRowIterator >
{
#ifdef EIGEN_USE_MKL_ALL
	typedef Eigen::PardisoLDLT< Eigen::SparseMatrix< InternalReal > > Eigen_Solver;
#else // !EIGEN_USE_MKL_ALL
	typedef Eigen::SimplicialLDLT< Eigen::SparseMatrix< InternalReal > > Eigen_Solver;
#endif // EIGEN_USE_MKL_ALL
	typedef Eigen::Matrix< InternalReal , Eigen::Dynamic , 1 > Eigen_Vector;
	Eigen_Solver _solver;
	Eigen_Vector _eigenB;
	Eigen_Vector _eigenX;
public:
	EigenSolverCholeskyLDLt( const SparseMatrixInterface< Real , MatrixRowIterator >& M , bool analyzeOnly=false )
	{
		Eigen::SparseMatrix< InternalReal > eigenM( (int)M.Rows() , (int)M.Rows() );
		std::vector< Eigen::Triplet< InternalReal > > triplets;
		triplets.reserve( M.Entries() );
		for( int i=0 ; i<M.Rows() ; i++ ) for( MatrixRowIterator iter=M.begin(i) ; iter!=M.end(i) ; iter++ ) triplets.push_back( Eigen::Triplet< InternalReal >( i , iter->N , iter->Value ) );
		eigenM.setFromTriplets( triplets.begin() , triplets.end() );
		_solver.analyzePattern( eigenM );
		if( !analyzeOnly )
		{
			_solver.factorize( eigenM );
			if( _solver.info()!=Eigen::Success ) fprintf( stderr , "[ERROR] EigenSolverCholeskyLDLt::EigenSolverCholeskyLDLt Failed to factorize matrix\n" ) , exit(0);
		}
		_eigenB.resize( M.Rows() );
		_eigenX.resize( M.Rows() );
	}
	void update( const SparseMatrixInterface< Real , MatrixRowIterator >& M )
	{
		Eigen::SparseMatrix< InternalReal > eigenM( (int)M.Rows() , (int)M.Rows() );
		std::vector< Eigen::Triplet< InternalReal > > triplets;
		triplets.reserve( M.Entries() );
		for( int i=0 ; i<M.Rows() ; i++ ) for( MatrixRowIterator iter=M.begin(i) ; iter!=M.end(i) ; iter++ ) triplets.push_back( Eigen::Triplet< InternalReal >( i , iter->N , iter->Value ) );
		eigenM.setFromTriplets( triplets.begin() , triplets.end() );
		_solver.factorize( eigenM );
		if( _solver.info()!=Eigen::Success ) fprintf(stderr, "[ERROR] EigenSolverCholeskyLDLt::update Failed to factorize matrix\n" ) , exit(0);
	}
	void solve( const Eigen_Vector& b , Eigen_Vector& x ){ x = _solver.solve(b); }
	void solve( ConstPointer(Real) b , Pointer(Real) x )
	{
#pragma omp parallel for
		for( int i=0 ; i<_eigenB.size() ; i++ ) _eigenB[i] = b[i];
		_eigenX = _solver.solve(_eigenB);
#pragma omp parallel for
		for( int i=0 ; i<_eigenX.size() ; i++ ) x[i] = (Real)_eigenX[i];
	}
	size_t dimension( void ) const { return _eigenB.size(); }
	static void Solve( const SparseMatrixInterface< Real , MatrixRowIterator >& M , ConstPointer(Real) b , Pointer(Real) x ){ EigenSolverCholeskyLDLt solver(M) ; solver.solve(b,x); }
};
template< class Real, class MatrixRowIterator >
class EigenSolverCG : public EigenSolver< Real, MatrixRowIterator >
{
#if 1
	//	Eigen::ConjugateGradient< Eigen::SparseMatrix< double > , Eigen::Lower , Eigen::IncompleteLUT< double > > _solver;
	Eigen::ConjugateGradient< Eigen::SparseMatrix< double > > _solver;
#else
	Eigen::BiCGSTAB< Eigen::SparseMatrix< double > > _solver;
#endif
	Eigen::VectorXd _eigenB, _eigenX;
	Eigen::SparseMatrix< double > _eigenM;
public:
	EigenSolverCG(const SparseMatrixInterface< Real, MatrixRowIterator >& M, int iters = 20)
	{
		_eigenM.resize((int)M.Rows(), (int)M.Rows());
		std::vector< Eigen::Triplet< double > > triplets;
		triplets.reserve(M.Entries());
		for (int i = 0; i<M.Rows(); i++) for (MatrixRowIterator iter = M.begin(i); iter != M.end(i); iter++) triplets.push_back(Eigen::Triplet< double >(i, iter->N, iter->Value));
		_eigenM.setFromTriplets(triplets.begin(), triplets.end());
		_solver.compute(_eigenM);
		_solver.analyzePattern(_eigenM);
		if (_solver.info() != Eigen::Success) fprintf(stderr, "[ERROR] EigenSolverCG::EigenSolverCG Failed to factorize matrix\n"), exit(0);
		_eigenB.resize(M.Rows()), _eigenX.resize(M.Rows());
		_solver.setMaxIterations(iters);
	}
	void update(const SparseMatrixInterface< Real, MatrixRowIterator >& M)
	{
#pragma omp parallel for
		for (int i = 0; i<M.Rows(); i++) for (MatrixRowIterator iter = M.begin(i); iter != M.end(i); iter++) _eigenM.coeffRef(i, iter->N) = iter->Value;
		_solver.compute(_eigenM);
		_solver.analyzePattern(_eigenM);
		if (_solver.info() != Eigen::Success) fprintf(stderr, "[ERROR] EigenSolverCG::update Failed to factorize matrix\n"), exit(0);
	}

	void setIters(int iters) { _solver.setMaxIterations(iters); }
	void solve(const Real* b, Real* x)
	{
#pragma omp parallel for
		for (int i = 0; i<_eigenB.size(); i++) _eigenB[i] = b[i], _eigenX[i] = x[i];
		_eigenX = _solver.solveWithGuess(_eigenB, _eigenX);
#pragma omp parallel for
		for (int i = 0; i<_eigenX.size(); i++) x[i] = _eigenX[i];
	}
	size_t dimension(void) const { return _eigenB.size(); }
	static void Solve(const SparseMatrixInterface< Real, MatrixRowIterator >& M, const Real* b, Real* x, int iters) { EigenSolverCG solver(M, iters); solver.solve(b, x); }
};

#else // !NEW_EIGEN
class EigenCholeskySolverLLt
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
	template< class Real, class MatrixRowIterator >
	EigenCholeskySolverLLt(const SparseMatrixInterface< Real, MatrixRowIterator >& M, bool analyzeOnly = false)
	{
#ifdef STORE_EIGEN_MATRIX
		_eigenM.resize(int(M.Rows()), int(M.Rows()));
#else // !STORE_EIGEN_MATRIX
		Eigen::SparseMatrix< double > eigenM(int(M.Rows()), int(M.Rows()));
#endif // STORE_EIGEN_MATRIX
		std::vector< Eigen::Triplet< double > > triplets;
		triplets.reserve(M.Entries());
		for (int i = 0; i<M.Rows(); i++) for (MatrixRowIterator iter = M.begin(i); iter != M.end(i); iter++) triplets.push_back(Eigen::Triplet< double >(i, iter->N, iter->Value));
#ifdef STORE_EIGEN_MATRIX
		_eigenM.setFromTriplets(triplets.begin(), triplets.end());
		_solver.analyzePattern(_eigenM);
#else // !STORE_EIGEN_MATRIX
		eigenM.setFromTriplets(triplets.begin(), triplets.end());
		_solver.analyzePattern(eigenM);
#endif // STORE_EIGEN_MATRIX
		if (!analyzeOnly)
		{
#ifdef STORE_EIGEN_MATRIX
			_solver.factorize(_eigenM);
#else // !STORE_EIGEN_MATRIX
			_solver.factorize(eigenM);
#endif // STORE_EIGEN_MATRIX
			if (_solver.info() != Eigen::Success) fprintf(stderr, "[ERROR] EigenCholeskySovler::EigenCholeskySolverLLt Failed to factorize matrix\n"), exit(0);
		}
		_eigenB.resize(M.Rows());
	}
	template< class Real, class MatrixRowIterator >
	void update(const SparseMatrixInterface< Real, MatrixRowIterator >& M)
	{
#ifdef STORE_EIGEN_MATRIX
#pragma omp parallel for
		for (int i = 0; i<M.Rows(); i++) for (MatrixRowIterator iter = M.begin(i); iter != M.end(i); iter++) _eigenM.coeffRef(i, iter->N) = iter->Value;
		_solver.factorize(_eigenM);
#else // !STORE_EIGEN_MATRIX
		Eigen::SparseMatrix< double > eigenM(int(M.Rows()), int(M.Rows()));
		std::vector< Eigen::Triplet< double > > triplets;
		triplets.reserve(M.Entries());
		for (int i = 0; i<M.Rows(); i++) for (MatrixRowIterator iter = M.begin(i); iter != M.end(i); iter++) triplets.push_back(Eigen::Triplet< double >(i, iter->N, iter->Value));
		eigenM.setFromTriplets(triplets.begin(), triplets.end());
		_solver.factorize(eigenM);
#endif // STORE_EIGEN_MATRIX
		switch (_solver.info())
		{
		case Eigen::Success: break;
		case Eigen::NumericalIssue: fprintf(stderr, "[ERROR] EigenCholeskySovler::EigenCholeskySolverLLt Failed to factorize matrix (numerical issue)\n"), exit(0);
		case Eigen::NoConvergence:  fprintf(stderr, "[ERROR] EigenCholeskySovler::EigenCholeskySolverLLt Failed to factorize matrix (no convergence)\n"), exit(0);
		case Eigen::InvalidInput:   fprintf(stderr, "[ERROR] EigenCholeskySovler::EigenCholeskySolverLLt Failed to factorize matrix (invalid input)\n"), exit(0);
		default: fprintf(stderr, "[ERROR] EigenCholeskySovler::EigenCholeskySolverLLt Failed to factorize matrix\n"), exit(0);
		}
	}
	template< class Real >
	void solve(ConstPointer(Real) b, Pointer(Real) x)
	{
#pragma omp parallel for
		for (int i = 0; i<_eigenB.size(); i++) _eigenB[i] = b[i];
		Eigen_Vector eigenX = _solver.solve(_eigenB);
#pragma omp parallel for
		for (int i = 0; i<eigenX.size(); i++) x[i] = (Real)eigenX[i];
	}
	size_t dimension(void) const { return _eigenB.size(); }
	template< class Real, class MatrixRowIterator >
	static void Solve(const SparseMatrixInterface< Real, MatrixRowIterator >& M, ConstPointer(Real) b, Pointer(Real) x) { EigenCholeskySolverLLt solver(M); solver.solve(b, x); }
};
class EigenCholeskySolverLDLt
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
	template< class Real, class MatrixRowIterator >
	EigenCholeskySolverLDLt(const SparseMatrixInterface< Real, MatrixRowIterator >& M, bool analyzeOnly = false)
	{
		Eigen::SparseMatrix< double > eigenM(int(M.Rows()), int(M.Rows()));
		std::vector< Eigen::Triplet<double> > triplets;
		triplets.reserve(M.Entries());
		for (int i = 0; i<M.Rows(); i++) for (MatrixRowIterator iter = M.begin(i); iter != M.end(i); iter++) triplets.push_back(Eigen::Triplet< double >(i, iter->N, iter->Value));
		eigenM.setFromTriplets(triplets.begin(), triplets.end());
		_solver.analyzePattern(eigenM);
		if (!analyzeOnly)
		{
			_solver.factorize(eigenM);
			if (_solver.info() != Eigen::Success) fprintf(stderr, "[ERROR] EigenCholeskySovler::EigenCholeskySolverLDLt Failed to factorize matrix\n"), exit(0);
		}
		_eigenB.resize(M.Rows());
	}
	template< class Real, class MatrixRowIterator >
	void update(const SparseMatrixInterface< Real, MatrixRowIterator >& M)
	{
		Eigen::SparseMatrix< double > eigenM(int(M.Rows()), int(M.Rows()));
		std::vector< Eigen::Triplet<double> > triplets;
		triplets.reserve(M.Entries());
		for (int i = 0; i<M.Rows(); i++) for (MatrixRowIterator iter = M.begin(i); iter != M.end(i); iter++) triplets.push_back(Eigen::Triplet< double >(i, iter->N, iter->Value));
		eigenM.setFromTriplets(triplets.begin(), triplets.end());
		_solver.factorize(eigenM);
		if (_solver.info() != Eigen::Success) fprintf(stderr, "[ERROR] EigenCholeskySovler::EigenCholeskySolverLDLt Failed to factorize matrix\n"), exit(0);
	}
	template< class Real >
	void solve(ConstPointer(Real) b, Pointer(Real) x)
	{
#pragma omp parallel for
		for (int i = 0; i<_eigenB.size(); i++) _eigenB[i] = b[i];
		Eigen_Vector eigenX = _solver.solve(_eigenB);
#pragma omp parallel for
		for (int i = 0; i<eigenX.size(); i++) x[i] = (Real)eigenX[i];
	}
	size_t dimension(void) const { return _eigenB.size(); }
	template< class Real, class MatrixRowIterator >
	static void Solve(const SparseMatrixInterface< Real, MatrixRowIterator >& M, ConstPointer(Real) b, Pointer(Real) x) { EigenCholeskySolverLDLt solver(M); solver.solve(b, x); }
};
class EigenCGSolver
{
#if 1
	//	Eigen::ConjugateGradient< Eigen::SparseMatrix< double > , Eigen::Lower , Eigen::IncompleteLUT< double > > _solver;
	Eigen::ConjugateGradient< Eigen::SparseMatrix< double > > _solver;
#else
	Eigen::BiCGSTAB< Eigen::SparseMatrix< double > > _solver;
#endif
	Eigen::VectorXd _eigenB, _eigenX;
	Eigen::SparseMatrix< double > _eigenM;
public:
	template< class Real, class MatrixRowIterator >
	EigenCGSolver(const SparseMatrixInterface< Real, MatrixRowIterator >& M)
	{
		_eigenM.resize((int)M.Rows(), (int)M.Rows());
		std::vector< Eigen::Triplet< double > > triplets;
		triplets.reserve(M.Entries());
		for (int i = 0; i<M.Rows(); i++) for (MatrixRowIterator iter = M.begin(i); iter != M.end(i); iter++) triplets.push_back(Eigen::Triplet< double >(i, iter->N, iter->Value));
		_eigenM.setFromTriplets(triplets.begin(), triplets.end());
		_solver.compute(_eigenM);
		_solver.analyzePattern(_eigenM);
		if (_solver.info() != Eigen::Success) fprintf(stderr, "[ERROR] EigenCGSolver::EigenCGSolver Failed to factorize matrix\n"), exit(0);
		_eigenB.resize(M.Rows()), _eigenX.resize(M.Rows());
	}
	template< class Real, class MatrixRowIterator >
	void update(const SparseMatrixInterface< Real, MatrixRowIterator >& M)
	{
#pragma omp parallel for
		for (int i = 0; i<M.Rows(); i++) for (MatrixRowIterator iter = M.begin(i); iter != M.end(i); iter++) _eigenM.coeffRef(i, iter->N) = iter->Value;
		_solver.compute(_eigenM);
		_solver.analyzePattern(_eigenM);
		if (_solver.info() != Eigen::Success) fprintf(stderr, "[ERROR] EigenCGSolver::EigenCGSolver Failed to factorize matrix\n"), exit(0);
	}

	template< class Real >
	void solve(const Real* b, Real* x, int iters = 20)
	{
		_solver.setMaxIterations(iters);
#pragma omp parallel for
		for (int i = 0; i<_eigenB.size(); i++) _eigenB[i] = b[i], _eigenX[i] = x[i];
		_eigenX = _solver.solveWithGuess(_eigenB, _eigenX);
#pragma omp parallel for
		for (int i = 0; i<_eigenX.size(); i++) x[i] = _eigenX[i];
	}
	size_t dimension(void) const { return _eigenB.size(); }
	template< class Real, class MatrixRowIterator >
	static void Solve(const SparseMatrixInterface< Real, MatrixRowIterator >& M, const Real* b, Real* x, int iters) { EigenCGSolver solver(M); solver._solver.setMaxIterations(iters); solver.solve(b, x); }
};
#endif // NEW_EIGEN
#endif // USE_EIGEN

#ifdef USE_CHOLMOD
template< int channels >
class CholmodSolver
{
public:
	const static bool LOWER_TRIANGULAR = true;
	size_t dim;
	cholmod_factor* cholmod_L;
	cholmod_dense*  cholmod_b;
	cholmod_sparse* cholmod_M;
	std::vector< bool > flaggedValues;


	template< class Real, class MatrixRowIterator > void   _init(const SparseMatrixInterface< Real, MatrixRowIterator >& M);
	template< class Real, class MatrixRowIterator > bool _update(const SparseMatrixInterface< Real, MatrixRowIterator >& M);
	CholmodSolver() { }

	cholmod_common cholmod_C;
	bool cholmod_C_set = false;


	template< class Real, class MatrixRowIterator >
	CholmodSolver(const SparseMatrixInterface< Real, MatrixRowIterator >& M);
	~CholmodSolver(void);

	template< class Real > void solve(ConstPointer(Real) b, Pointer(Real) x);
	int nonZeros(void) const;
};


template< int channels > template< class Real, class MatrixRowIterator > CholmodSolver< channels >::CholmodSolver(const SparseMatrixInterface< Real, MatrixRowIterator >& M) { _init(M), _update(M); }

template< int channels > template< class Real, class MatrixRowIterator >
void CholmodSolver< channels >::_init(const SparseMatrixInterface< Real, MatrixRowIterator >& M)
{
	if (!cholmod_C_set) cholmod_start(&cholmod_C);
	cholmod_C_set = true;

	dim = M.Rows();

	int maxEntries;
	if (LOWER_TRIANGULAR)
	{
		maxEntries = (int)((M.Entries() - M.Rows()) / 2 + M.Rows());
		cholmod_M = cholmod_allocate_sparse(dim, dim, maxEntries, 0, 1, -1, CHOLMOD_REAL, &cholmod_C);
	}
	else
	{
		maxEntries = (int)M.Entries();
		cholmod_M = cholmod_allocate_sparse(dim, dim, maxEntries, 0, 1, 0, CHOLMOD_REAL, &cholmod_C);
	}
	cholmod_M->i = malloc(sizeof(int) * maxEntries);
	cholmod_M->x = malloc(sizeof(double) * maxEntries);

	int *_p = (int*)cholmod_M->p;
	int *_i = (int*)cholmod_M->i;

	int off = 0;
	dim = 0;

	for (int i = 0; i<M.Rows(); i++)
	{
		_p[dim++] = off;
		for (MatrixRowIterator iter = M.begin(i); iter != M.end(i); iter++) if (!LOWER_TRIANGULAR || iter->N >= i) _i[off++] = iter->N;
	}
	_p[dim] = off;

	clock_t begin;
	if (0) begin = clock();
	cholmod_L = cholmod_analyze(cholmod_M, &cholmod_C);
	if (0) printf("Cholmod analyze time = %.4f \n", double(clock() - begin) / CLOCKS_PER_SEC);

	cholmod_b = cholmod_allocate_dense(dim, channels, dim, cholmod_M->xtype, &cholmod_C);
}
template< int channels > template< class Real, class MatrixRowIterator >
bool CholmodSolver< channels >::_update(const SparseMatrixInterface< Real, MatrixRowIterator >& M)
{
	double *_x = (double*)cholmod_M->x;
	int off = 0;

	int *_p = (int*)cholmod_M->p;
#pragma omp parallel for
	for (int i = 0; i<M.Rows(); i++)
	{
		int off = (int)_p[i];
		for (MatrixRowIterator iter = M.begin(i); iter != M.end(i); iter++)if (!LOWER_TRIANGULAR || iter->N >= i) _x[off++] = double(iter->Value);
	}

	cholmod_C.print = 0;

	clock_t begin;
	if (0) begin = clock();
	cholmod_factorize(cholmod_M, cholmod_L, &cholmod_C);
	if (0) printf("Cholmod factorize time = %.4f \n", double(clock() - begin) / CLOCKS_PER_SEC);

	if (cholmod_C.status == CHOLMOD_NOT_POSDEF)
	{
		fprintf(stderr, "[WARNING] Matrix not positive-definite\n");
		return false;
	}
	else if (cholmod_C.status == CHOLMOD_OUT_OF_MEMORY)
	{
		fprintf(stderr, "[WARNING] CHOLMOD ran out of memory\n");
		return false;
	}
	else if (cholmod_C.status != CHOLMOD_OK)
	{
		fprintf(stderr, "[WARNING] CHOLMOD status not OK: %d\n", cholmod_C.status);
		return false;
	}
	return true;
}

template< int channels >
CholmodSolver< channels >::~CholmodSolver(void)
{
	if (cholmod_L) cholmod_free_factor(&cholmod_L, &cholmod_C), cholmod_L = NULL;
	if (cholmod_b) cholmod_free_dense(&cholmod_b, &cholmod_C), cholmod_b = NULL;
	if (cholmod_M) cholmod_free_sparse(&cholmod_M, &cholmod_C), cholmod_M = NULL;
}

template< int channels >
template< class Real >
void CholmodSolver< channels >::solve(ConstPointer(Real) b, Pointer(Real) x)
{
	size_t numEntries = dim*channels;
	double* _b = (double*)cholmod_b->x;

#pragma omp parallel for
	for (int i = 0; i<numEntries; i++) _b[i] = (double)b[i];

	clock_t begin;
	if (0) begin = clock();
	cholmod_dense* cholmod_x = cholmod_solve(CHOLMOD_A, cholmod_L, cholmod_b, &cholmod_C);
	if (0) printf("Cholmod solve time = %.4f \n", double(clock() - begin) / CLOCKS_PER_SEC);
	double* _x = (double*)cholmod_x->x;

#pragma omp parallel for
	for (int i = 0; i<numEntries; i++) x[i] = (Real)_x[i];

	cholmod_free_dense(&cholmod_x, &cholmod_C);
}

template< int channels >
int CholmodSolver< channels >::nonZeros(void) const
{
	long long nz = 0;
	if (cholmod_L->xtype != CHOLMOD_PATTERN && !(cholmod_L->is_super)) for (int i = 0; i<cholmod_L->n; i++) nz += ((int*)cholmod_L->nz)[i];
	bool examine_super = false;
	if (cholmod_L->xtype != CHOLMOD_PATTERN) examine_super = true;
	else                                      examine_super = (((int*)cholmod_L->s)[0] != (-1));
	if (examine_super)
	{
		/* check and print each supernode */
		for (int s = 0; s < cholmod_L->nsuper; s++)
		{
			int k1 = ((int*)cholmod_L->super)[s];
			int k2 = ((int*)cholmod_L->super)[s + 1];
			int psi = ((int*)cholmod_L->pi)[s];
			int psend = ((int*)cholmod_L->pi)[s + 1];
			int nsrow = psend - psi;
			int nscol = k2 - k1;
			nz += nscol * nsrow - (nscol*nscol - nscol) / 2;
		}
	}
	return (int)nz;
}

#endif // USE_CHOLMOD
#endif // LINEAR_SOLVERS_INCLUDE