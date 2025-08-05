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


#include <Misha/SparseMatrixInterface.h>
#include <Misha/MultiThreading.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace MishaK
{
	namespace TSP
	{
		template< typename Real , class MatrixRowIterator >
		Eigen::SparseMatrix< double > ToEigen( const SparseMatrixInterface< Real , MatrixRowIterator > &M )
		{
			Eigen::SparseMatrix< double > _M( (int)M.Rows() , (int)M.Rows() );

			std::vector< Eigen::Triplet< double > > triplets;
			triplets.reserve( M.Entries() );
			for( int i=0 ; i<M.Rows() ; i++ ) for( MatrixRowIterator iter=M.begin(i) ; iter!=M.end(i) ; iter++ ) triplets.push_back( Eigen::Triplet< double >( i , iter->N , iter->Value ) );
			_M.setFromTriplets( triplets.begin() , triplets.end() );
			return _M;
		}

		template< typename Solver >
		struct EigenSolverWrapper
		{
			EigenSolverWrapper( void ) : _solver( new Solver() ){}

			template< typename Real , class MatrixRowIterator >
			EigenSolverWrapper( const SparseMatrixInterface< Real , MatrixRowIterator > &M , bool numericFactorization=true ) : EigenSolverWrapper()
			{
				init( M , numericFactorization );
			}
			~EigenSolverWrapper( void ){ delete _solver; }

			EigenSolverWrapper( const EigenSolverWrapper & esw ) = delete;
			EigenSolverWrapper( EigenSolverWrapper && esw )
			{
				_solver = esw._solver;
				esw._solver = nullptr;
			}

			EigenSolverWrapper & operator = ( const EigenSolverWrapper & esw ) = delete;
			EigenSolverWrapper & operator = ( EigenSolverWrapper && esw ){ std::swap( esw._solver , _solver ); }

			template< typename Real , class MatrixRowIterator >
			void init( const SparseMatrixInterface< Real , MatrixRowIterator > &M , bool numericFactorization=true )
			{
				Eigen::SparseMatrix< double > _M = ToEigen( M );
				_solver->analyzePattern( _M );
				if( numericFactorization )
				{
					_solver->factorize( _M );
					if( _solver->info()!=Eigen::Success ) MK_THROW( "Failed to factorize matrix" );
				}
				_rows = _M.rows();
				_cols = _M.cols();
			}

			template< typename Real , class MatrixRowIterator >
			void update( const SparseMatrixInterface< Real , MatrixRowIterator > & M )
			{
				Eigen::SparseMatrix< double > _M = ToEigen( M );
				_solver->factorize( _M );
				if( _solver->info()!=Eigen::Success ) MK_THROW( "Failed to factorize matrix" );
				_rows = _M.rows();
				_cols = _M.cols();
			}

			template< typename Real >
			void solve( Pointer( Real ) x , ConstPointer( Real ) b )
			{
				Eigen::VectorXd B( _rows );
				ThreadPool::ParallelFor( 0 , _rows , [&]( size_t r ){ B( r ) = b[r]; } );
				Eigen::VectorXd X = _solver->solve( B );
				ThreadPool::ParallelFor( 0 , _rows , [&]( size_t r ){ x[r] = X( r ); } );
			}

			template< typename Real >
			void solve( std::vector< Real >& x , const std::vector< Real >& b )
			{
				Eigen::VectorXd B( _rows );
				ThreadPool::ParallelFor( 0 , _rows , [&]( size_t r ){ B( r ) = b[r]; } );
				Eigen::VectorXd X = _solver->solve( B );
				ThreadPool::ParallelFor( 0 , _rows , [&]( size_t r ){ x[r] = X( r ); } );
			}

			template< typename Real , unsigned int Channels >
			void solve( std::vector< Point< Real , Channels > >& x , const std::vector< Point< Real , Channels > >& b )
			{
				Eigen::MatrixXd B( _rows , Channels );
				ThreadPool::ParallelFor( 0 , _rows , [&]( size_t r ){ for( unsigned int c=0 ; c<Channels ; c++ ) B( r , c ) = b[r][c]; } );
				Eigen::MatrixXd X = _solver->solve( B );
				ThreadPool::ParallelFor( 0 , _rows , [&]( size_t r ){ for( unsigned int c=0 ; c<Channels ; c++ ) x[r][c] = X( r , c ); } );
			}
		protected:
			Solver * _solver;
			size_t _rows , _cols;
		};

	}
}
