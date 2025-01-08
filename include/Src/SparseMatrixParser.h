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
#include <set>
#include <Misha/MultiThreading.h>
#include <Misha/Atomic.h>

template< class Real , class RealOut=Real >
SparseMatrix< RealOut , int > SetSparseMatrix( const std::vector< Eigen::Triplet< Real > >& triplets , int cols , int /*rows*/ , bool rowMajor )
{
	SparseMatrix< RealOut , int > M;
	M.resize( cols );
	if( rowMajor ) ThreadPool::ParallelFor( 0 , triplets.size() , [&]( unsigned int , size_t i ){ Misha::AddAtomic< size_t >( M.rowSizes[ triplets[i].col() ] , 1 ); } );
	else           ThreadPool::ParallelFor( 0 , triplets.size() , [&]( unsigned int , size_t i ){ Misha::AddAtomic< size_t >( M.rowSizes[ triplets[i].row() ] , 1 ); } );

	ThreadPool::ParallelFor
		(
			0 , cols ,
			[&]( unsigned int , size_t i )
			{
				if( M.rowSizes[i] )
				{
					int s = (int)M.rowSizes[i];
					M.SetRowSize(i,s);
					M.rowSizes[i] = 0;
				}
			}
		);

	for( int i=0 ; i<triplets.size() ; i++ )
	{
		int col = rowMajor ? triplets[i].row() : triplets[i].col() , row = rowMajor ? triplets[i].col() : triplets[i].row();
		M[row][M.rowSizes[row]++] = MatrixEntry< RealOut , int >( col , (RealOut)triplets[i].value() );
	}

	ThreadPool::ParallelFor
	(
		0 , M.rows ,
		[&]( unsigned int , size_t i )
		{
			std::map< int , RealOut > rowMap;
			for( int j=0 ; j<M.rowSizes[i] ; j++ ) rowMap[M[i][j].N] += M[i][j].Value;
			M.SetRowSize( i , rowMap.size() );
			int count = 0;
			for( auto it=rowMap.begin() ; it!=rowMap.end() ; it++ ) M[i][count++] = MatrixEntry< RealOut , int >( it->first , it->second );
		}
	);

	return M;
}
template< typename Real >
class matrixRowEntry
{
public:
	matrixRowEntry (int p_index , Real p_value ) : index(p_index) , value(p_value) {}
	int index;
	Real value;
};

template< typename Real >
struct matrixRowEntryCompare
{
	bool operator() (const matrixRowEntry< Real > &lhs , const matrixRowEntry< Real > &rhs ) const { return lhs.index <rhs.index; }
};

template< class RealIn, class RealOut>
void SparseMatrixParser(const SparseMatrix< RealIn, int > & _M, Eigen::SparseMatrix<RealOut> & M) {
	int rows = (int)_M.rows;
	int cols = 0;
	std::vector<Eigen::Triplet<RealOut>> triplets;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < _M.rowSizes[i]; j++) {
			int col = _M[i][j].N;
			RealIn value = _M[i][j].Value;
			triplets.push_back(Eigen::Triplet<RealOut>(i, col, RealOut(value)));
			cols = std::max<int>(col, cols);
		}
	}
	M.resize(rows, cols + 1);
	M.setFromTriplets(triplets.begin(), triplets.end());
}
