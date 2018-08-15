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

#define USE_CHOLMOD 0
#define USE_EIGEN_SIMPLICIAL 1
#define USE_EIGEN_PARDISO 0

#if USE_CHOLMOD
#include <Misha/LinearSolvers.h>

#define CHOLMOD_CHANNELS_IN_PARALLEL 1
#if CHOLMOD_CHANNELS_IN_PARALLEL
template <class Real>
class CholmodCholeskySolver3{
public:
	CholmodSolver<1> solver[3];
	std::vector<Real> out[3];
	std::vector<Real> in[3];
	void init(const SparseMatrix<double, int> & M){
#pragma omp parallel for
		for (int c = 0; c < 3; c++){
			solver[c]._init(M);
		}

		const int numVariables = M.Rows();
		for (int c = 0; c < 3; c++){
			out[c].resize(numVariables);
			in[c].resize(numVariables);
		}
	}

	void update(const SparseMatrix<double, int> & M){
#pragma omp parallel for
		for (int c = 0; c < 3; c++){
			solver[c]._update(M);
		}
	}
};

template <class Real, class DataType>
void solve( CholmodCholeskySolver3< Real >& chol , std::vector< DataType >& x0 , const std::vector< DataType >& rhs )
{
	int numVariables = x0.size();
#pragma omp parallel for
	for( int c=0 ; c<3 ; c++ )
	{
		for( int n=0 ; n<numVariables ; n++ ) chol.in[c][n] = rhs[n][c];
		chol.solver[c].solve( &chol.in[c][0] , &chol.out[c][0] );
	}
#pragma omp parallel for
	for( int n=0 ; n<numVariables ; n++ ) x0[n] = SetData< Real >( chol.out[0][n] , chol.out[1][n] , chol.out[2][n] );
}
#else 

#define CHOLMOD_CHANNELS_IN_BLOCK 1
#if CHOLMOD_CHANNELS_IN_BLOCK
template <class Real>
class CholmodCholeskySolver3{
public:
	CholmodSolver<3> solver;
	std::vector<Real> out;
	std::vector<Real> in;
	void init(const SparseMatrix<double, int> & M) {
		solver._init(M);
		const int numVariables = M.Rows();
		out.resize(3 * numVariables);
		in.resize(3 * numVariables);
	}

	void update(const SparseMatrix<double, int> & M) {
		solver._update(M);
	}
};

template <class Real, class DataType>
void solve(CholmodCholeskySolver3<Real> & chol, std::vector<DataType> & x0, const std::vector<DataType> & rhs){
	int numVariables = x0.size();
	for (int n = 0; n < numVariables; n++) for (int c = 0; c < 3; c++)chol.in[c* numVariables + n] = rhs[n][c];
	chol.solver.solve(&chol.in[0], &chol.out[0]);
	for( int n=0 ; n<numVariables ; n++ ) x0[n] = SetData< Real >( chol.out[ 0*numVariables+n ] , chol.out[ 1*numVariables+n ] , chol.out[ 2*numVariables+n ] );
}
#else
template <class Real>
class CholmodCholeskySolver3{
public:
	CholmodSolver<1> solver;
	std::vector<Real> out[3];
	std::vector<Real> in[3];
	void init(const SparseMatrix<double, int> & M) {
		solver._init(M);
		const int numVariables = M.Rows();
		for (int c = 0; c < 3; c++){
			out[c].resize(numVariables);
			in[c].resize(numVariables);
		}
	}

	void update(const SparseMatrix<double, int> & M) {
		solver._update(M);
	}
};

template <class Real, class DataType>
void solve(CholmodCholeskySolver3<Real> & chol, std::vector<DataType> & x0, const std::vector<DataType> & rhs) {
	int numVariables = x0.size();
	for (int c = 0; c < 3; c++) {
		for (int n = 0; n < numVariables; n++) chol.in[c][n] = rhs[n][c];
		chol.solver.solve(&chol.in[c][0], &chol.out[c][0]);
	}
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) x0[n] = SetData< Real >(chol.out[0][n], chol.out[1][n], chol.out[2][n]);
}
#endif
#endif

template <class Real>
class CholmodCholeskySolver1{
public:
	CholmodSolver<1> solver;
	std::vector<Real> out;
	std::vector<Real> in;
	void init(const SparseMatrix<double, int> & M) {
		solver._init(M);
		const int numVariables = M.Rows();
		out.resize(numVariables);
		in.resize(numVariables);
	}
	void update(const SparseMatrix<double, int> & M) {
		solver._update(M);
	}
};

template <class Real, class DataType>
void solve(CholmodCholeskySolver1<Real> & chol, std::vector<DataType> & x0, const std::vector<DataType> & rhs) {
	int numVariables = x0.size();
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) chol.in[n] = rhs[n];
	chol.solver.solve(&chol.in[0], &chol.out[0]);
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) x0[n] = chol.out[n];
}
#endif

#include <Eigen/Sparse>
#include <Eigen/Dense>

template <class Real>
class EigenCholeskySolver3
{
public:
	typedef typename Eigen::SimplicialLDLT< Eigen::SparseMatrix< Real > > Solver;
	//typedef typename Eigen::SimplicialLLT< Eigen::SparseMatrix< Real > > Solver;
	Solver* solver;
	typename ScalarPrecisionType< Real >::EigenVector x0_vectors[3] , rhs_vectors[3] , solution_vectors[3];

	void init(const SparseMatrix<double, int> & _M){
		Eigen::SparseMatrix< Real > M;
		SparseMatrixParser(_M, M);

		solver = new Solver();
		solver->analyzePattern(M);

		Eigen::ComputationInfo info = solver->info();
		if (info == Eigen::Success) {
		}
		else if (info == Eigen::NumericalIssue){
			printf("FAILED : Numerical issue! \n");
		}
		else if (info == Eigen::NoConvergence) {
			printf("FAILED : No convergence! \n");
		}
		else if (info == Eigen::InvalidInput) {
			printf("FAILED : Invalid input! \n");
		}
		else {
			printf("FAILED : Undetermined cause! \n");
		}

		const int numVariables = (int)M.rows();
		for (int c = 0; c < 3; c++) {
			x0_vectors[c].resize(numVariables);
			rhs_vectors[c].resize(numVariables);
			solution_vectors[c].resize(numVariables);
		}
	}
	void update(const SparseMatrix<double, int> & _M) {
		Eigen::SparseMatrix< Real > M;
		SparseMatrixParser(_M, M);
		solver->factorize(M);
	}
};

template <class Real,class DataType>
void solve(EigenCholeskySolver3<Real> & chol, std::vector<DataType> & x0, const std::vector<DataType> & rhs){
	int numVariables = (int)x0.size();
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) for (int c = 0; c < 3; c++) {
		chol.x0_vectors[c][n] = x0[n][c];
		chol.rhs_vectors[c][n] = rhs[n][c];
	}
#pragma omp parallel for
	for( int c=0 ; c<3 ; c++ ) chol.solution_vectors[c] = chol.solver->solve(chol.rhs_vectors[c]);
#pragma omp parallel for
	for( int n=0 ; n<numVariables ; n++ ) x0[n] = VectorPrecisionType< Real >::SetData( chol.solution_vectors[0][n] , chol.solution_vectors[1][n] , chol.solution_vectors[2][n] );
}

template <class Real>
class EigenCholeskySolver1{
public:
	Eigen::SimplicialLDLT< Eigen::SparseMatrix< Real > >* solver;
	typename ScalarPrecisionType< Real >::EigenVector x0_vector , rhs_vector , solution_vector;

	void init(const SparseMatrix<double, int> & _M) {
		Eigen::SparseMatrix< Real > M;
		SparseMatrixParser(_M, M);

		solver = new Eigen::SimplicialLDLT<Eigen::SparseMatrix<Real>>();
		solver->analyzePattern(M);
		Eigen::ComputationInfo info = solver->info();
		if (info == Eigen::Success) {
		}
		else if (info == Eigen::NumericalIssue) {
			printf("FAILED : Numerical issue! \n");
		}
		else if (info == Eigen::NoConvergence) {
			printf("FAILED : No convergence! \n");
		}
		else if (info == Eigen::InvalidInput) {
			printf("FAILED : Invalid input! \n");
		}
		else {
			printf("FAILED : Undetermined cause! \n");
		}

		const int numVariables = (int)M.rows();
		x0_vector.resize(numVariables);
		rhs_vector.resize(numVariables);
		solution_vector.resize(numVariables);
	}

	void update(const SparseMatrix<double, int> & _M){
		Eigen::SparseMatrix< Real > M;
		SparseMatrixParser(_M, M);
		solver->factorize(M);
	}
};

template <class Real, class DataType>
void solve(EigenCholeskySolver1<Real> & chol, std::vector<DataType> & x0, const std::vector<DataType> & rhs){
	int numVariables = (int)x0.size();
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) {
		chol.x0_vector[n] = x0[n];
		chol.rhs_vector[n] = rhs[n];
	}
	chol.solution_vector = chol.solver->solve(chol.rhs_vector);
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) x0[n] = chol.solution_vector[n];
}

#if USE_EIGEN_PARDISO
#pragma comment( lib , "mkl_intel_lp64.lib")
#pragma comment( lib , "mkl_intel_thread.lib")
#pragma comment( lib , "mkl_core.lib")
#pragma comment( lib , "libiomp5md.lib")

#include <Eigen/PardisoSupport>
//Pardiso Solver

template <class Real>
class EigenPardisoSolver3 {
public:
	Eigen::PardisoLDLT<Eigen::SparseMatrix<Real>>  * solver;
	EigenVector x0_vectors[3];
	EigenVector rhs_vectors[3];
	EigenVector solution_vectors[3];

	void init(const SparseMatrix<double, int> & _M) {
		Eigen::SparseMatrix< Real > M;
		SparseMatrixParser(_M, M);

		solver = new Eigen::PardisoLDLT<Eigen::SparseMatrix<Real>>();
		solver->analyzePattern(M);
		Eigen::ComputationInfo info = solver->info();
		if (info == Eigen::Success) {
		}
		else if (info == Eigen::NumericalIssue) {
			printf("FAILED : Numerical issue! \n");
		}
		else if (info == Eigen::NoConvergence) {
			printf("FAILED : No convergence! \n");
		}
		else if (info == Eigen::InvalidInput) {
			printf("FAILED : Invalid input! \n");
		}
		else {
			printf("FAILED : Undetermined cause! \n");
		}

		const int numVariables = M.rows();
		for (int c = 0; c < 3; c++) {
			x0_vectors[c].resize(numVariables);
			rhs_vectors[c].resize(numVariables);
			solution_vectors[c].resize(numVariables);
		}
	}
	void update(const SparseMatrix<double, int> & _M) {
		Eigen::SparseMatrix< Real > M;
		SparseMatrixParser(_M, M);
		solver->factorize(M);
	}
};

template <class Real, class DataType>
void solve(EigenPardisoSolver3<Real> & chol, std::vector<DataType> & x0, const std::vector<DataType> & rhs) {
	int numVariables = x0.size();
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) for (int c = 0; c < 3; c++) {
		chol.x0_vectors[c][n] = x0[n][c];
		chol.rhs_vectors[c][n] = rhs[n][c];
	}
	for (int c = 0; c < 3; c++) chol.solution_vectors[c] = chol.solver->solve(chol.rhs_vectors[c]);
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) x0[n] = SetData< Real >(chol.solution_vectors[0][n], chol.solution_vectors[1][n], chol.solution_vectors[2][n]);
}

template <class Real>
class EigenPardisoSolver1{
public:
	Eigen::PardisoLDLT<Eigen::SparseMatrix<Real>> * solver;
	EigenVector x0_vector;
	EigenVector rhs_vector;
	EigenVector solution_vector;

	void init(const SparseMatrix<double, int> & _M){
		Eigen::SparseMatrix< Real > M;
		SparseMatrixParser(_M, M);

		solver = new Eigen::PardisoLDLT<Eigen::SparseMatrix<Real>>();
		solver->analyzePattern(M);
		Eigen::ComputationInfo info = solver->info();
		if (info == Eigen::Success) {
		}
		else if (info == Eigen::NumericalIssue) {
			printf("FAILED : Numerical issue! \n");
		}
		else if (info == Eigen::NoConvergence) {
			printf("FAILED : No convergence! \n");
		}
		else if (info == Eigen::InvalidInput) {
			printf("FAILED : Invalid input! \n");
		}
		else {
			printf("FAILED : Undetermined cause! \n");
		}

		const int numVariables = M.rows();
		x0_vector.resize(numVariables);
		rhs_vector.resize(numVariables);
		solution_vector.resize(numVariables);
	}
	void update(const SparseMatrix<double, int> & _M) {
		Eigen::SparseMatrix< Real > M;
		SparseMatrixParser(_M, M);
		solver->factorize(M);
	}
};

template <class Real, class DataType>
void solve(EigenPardisoSolver1<Real> & chol, std::vector<DataType> & x0, const std::vector<DataType> & rhs) {
	int numVariables = x0.size();
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) {
		chol.x0_vector[n] = x0[n];
		chol.rhs_vector[n] = rhs[n];
	}
	chol.solution_vector = chol.solver->solve(chol.rhs_vector);
#pragma omp parallel for
	for (int n = 0; n < numVariables; n++) x0[n] = chol.solution_vector[n];
}
#endif
