/* -*- C++ -*-
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
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

#include <float.h>
#include <complex>
#include <unordered_map>
#include "Miscellany.h"

///////////////////
//  SparseMatrix //
///////////////////
///////////////////////////////////////
// SparseMatrix Methods and Memebers //
///////////////////////////////////////

template< class T , class IndexType >
SparseMatrix< T , IndexType >::SparseMatrix( void )
{
	rowSizes = NullPointer< size_t >( );
	rows = 0;
	_entries = NullPointer< Pointer( MatrixEntry< T , IndexType > ) >( );
}
////////////////////////////////////////////////////////////////////////////////
/*! Sets the column index of all allocated entries to -1 so that they are
//  treated as non-existent. This is needed because SetRowSize() uses malloc
//  instead of new and MatrixEntry's constructor is never called.
*///////////////////////////////////////////////////////////////////////////////
template <class T, class IndexType>
void SparseMatrix<T, IndexType>::invalidateEntries()
{
    IndexType numRows = Rows();
    for( IndexType r=0 ; r<numRows ; ++r )
	{
        IndexType numColumns = RowSize(r);
        MatrixEntry<T, IndexType> *row = _entries[r];
        for( IndexType c=0 ; c<numColumns ; ++c ) row[c].N = -1;
    }
}

////////////////////////////////////////////////////////////////////////////////
/*! Adds a scalar value to an element in the Matrix, using a new element if
//  necesary. If no pre-allocated space for a new element exists, false is
//  returned.
//  WARNING: no check is done to remove entries that become zero after addition
//  @param[in]  s   The scalar to add to the destination element in the matrix
//  @param[in]  i   The destination element's row
//  @param[in]  j   The destination element's column
//  @return     true if successful, false if no pre-allocated space exists for
//              the creation of a new nonzero element
*///////////////////////////////////////////////////////////////////////////////
template <class T, class IndexType>
bool SparseMatrix<T, IndexType>::addScalarToEntry(T s, IndexType i,
          IndexType j)
{
    // Don't add unset entries (possibly change this to use a tolerance)
    if ((s.real() == 0) && (s.imag() == 0))
        return true;

    MatrixEntry<T, IndexType> *row = _entries[i];
    IndexType rSize = RowSize(i);

    bool success = false;
    int availableIdx = -1;
    for (IndexType k = 0; !success && (k < rSize); ++k) {
        if (row[k].N == j) {
            row[k].Value += s;
            success = true;
        }
        if ((availableIdx == -1) && (row[k].N == (unsigned int) -1))
            availableIdx = k;
    }

    if (!success && (availableIdx != -1))   {
        row[availableIdx].Value = s;
        row[availableIdx].N = j;
        success = true;
    }

    return success;
}

template< class T , class IndexType >
SparseMatrix< T , IndexType >::SparseMatrix( size_t rows )
{
	rows = 0;
	rowSizes = NullPointer< size_t >( );
	_entries= NullPointer< Pointer( MatrixEntry< T , IndexType > ) >( );
	resize( rows );
}
template< class T , class IndexType >
SparseMatrix< T , IndexType >::SparseMatrix( const SparseMatrix& M )
{
	rowSizes = NullPointer< size_t >( );
	rows = 0;
	_entries = NullPointer< Pointer( MatrixEntry< T , IndexType > ) >( );
	resize( M.rows );
	for( int i=0 ; i<rows ; i++ )
	{
		SetRowSize( i , M.rowSizes[i] );
		for( int j=0 ; j<rowSizes[i] ; j++ ) _entries[i][j] = M._entries[i][j];
	}
}
template< class T , class IndexType >
SparseMatrix< T , IndexType >::SparseMatrix( SparseMatrix&& M )
{
	rowSizes = NullPointer< size_t >( );
	rows = 0;
	_entries = NullPointer< Pointer( MatrixEntry< T , IndexType > ) >( );

	Swap( *this , M );
}
template< class T , class IndexType >
template< class T2 , class IndexType2 >
SparseMatrix< T , IndexType >::SparseMatrix( const SparseMatrix< T2 , IndexType2 >& M )
{
	rowSizes = NullPointer< size_t >();
	rows = 0;
	_entries = NULL;
	resize( M.rows );
	for( int i=0 ; i<rows ; i++ )
	{
		SetRowSize( i , M.rowSizes[i] );
		for( int j=0 ; j<rowSizes[i] ; j++ ) _entries[i][j] = MatrixEntry< T , IndexType >( M._entries[i][j].N , T( M._entries[i][j].Value ) );
	}
}

template< class T , class IndexType >
double SparseMatrix< T , IndexType >::SquareNorm(void) const
{
	double l2=0;
	for( int i=0 ; i<rows ; i++ ) for( int j=0 ; j<rowSizes[i] ; j++ ) l2 += _entries[i][j].Value * _entries[i][j].Value;
	return l2;
}
template< class T , class IndexType >
double SparseMatrix< T , IndexType >::ASymmetricSquareNorm(void) const
{
	double l2 = 0;
	for( int i=0 ; i<rows ; i++ ) for( int j=0 ; j<rowSizes[i] ; j++ )
	{
		double t1=0 , t2=0;
		int N = _entries[i][j].N;
		if( N==i ) continue;
		// [WARNING] multi-counting 
		for( int k=0 ; k<rowSizes[i] ; k++ ) if( _entries[i][k].N==N ) t1 += _entries[i][k].Value;
		for( int k=0 ; k<rowSizes[N] ; k++ ) if( _entries[N][k].N==i ) t2 += _entries[N][k].Value;
		l2 += (t1-t2)*(t1-t2);
	}
	return l2;
}
template< class T , class IndexType >
double SparseMatrix< T , IndexType >::AHermitianSquareNorm(void) const
{
	double l2=0;
	for( int i=0 ; i<rows ; i++ ) for( int j=0 ; j<rowSizes[i] ; j++ )
	{
		std::complex< double > t , t1 = 0. , t2 = 0.;
		int N=_entries[i][j].N;
		if( N==i ) continue;
		t1 = _entries[i][j].Value;
		for( int k=0 ; k<rowSizes[N] ; k++ ) if( _entries[N][k].N==i )
		{
			t2 = _entries[N][k].Value;
			t2 = std::complex< double >( t2.real() , -t2.imag() );
			break;
		}
		t = t1-t2;
		l2 += t.real()*t.real() + t.imag()*t.imag();
	}
	return l2;
}
template< class T , class IndexType >
template< class T2 , class IndexType2 >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::copy( const SparseMatrix< T2 , IndexType2 >& M  )
{
	resize( M.rows );
	for ( int i=0 ; i<rows ; i++)
	{
		SetRowSize( i , M.rowSizes[i] );
		for( int j=0 ; j<rowSizes[i] ; j++ )
		{
			int idx = M._entries[i][j].N;
			_entries[i][j] = MatrixEntry< T , IndexType >( idx , T( M[i][j].Value ) );
		}
	}
	return *this;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::operator = ( SparseMatrix< T , IndexType >&& M )
{
	Swap( *this , M );
	return *this;
}

template< class T , class IndexType >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::operator = ( const SparseMatrix< T , IndexType >& M )
{
	resize( M.rows );
	for( int i=0 ; i<rows ; i++ )
	{
		SetRowSize( i , M.rowSizes[i] );
		for( int j=0 ; j<rowSizes[i] ; j++ ) _entries[i][j]=M._entries[i][j];
	}
	return *this;
}
template< class T , class IndexType >
template< class T2 , class IndexType2 >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::operator = (const SparseMatrix< T2 , IndexType2 >& M)
{
	resize( M.rows );
	for ( int i=0 ; i<rows ; i++ )
	{
		SetRowSize( i , M.rowSizes[i] );
		for( int j=0 ; j<rowSizes[i] ; j++ ) _entries[i][j] = MatrixEntry< T , IndexType >( M._entries[i][j].N , T( M._entries[i][j].Value ) );
	}
	return *this;
}
template< class T , class IndexType >
template< class T2 >
void SparseMatrix< T , IndexType >::operator() ( const T2* in , T2* out ) const { Interface::Multiply( in , out ); }


template< class T , class IndexType > SparseMatrix< T , IndexType >::~SparseMatrix( void ) { resize( 0 ); }

template< class T , class IndexType >
void SparseMatrix< T , IndexType >::resize( size_t r )
{
	if( rows>0 )
	{
		for( int i=0 ; i<rows ; i++ ) FreePointer( _entries[i] );
		FreePointer( _entries );
		FreePointer( rowSizes );
	}
	rows = r;
	if( r )
	{
		rowSizes = AllocPointer< size_t >( r ) , memset( rowSizes , 0 , sizeof(size_t)*r );
		_entries = AllocPointer< Pointer( MatrixEntry< T , IndexType > ) >( r );
		for( int i=0 ; i<r ; i++ ) _entries[i] = NullPointer< MatrixEntry< T , IndexType > >( );
	}
}

template< class T , class IndexType >
void SparseMatrix< T , IndexType >::SetRowSize( size_t row , size_t count )
{
	if( row>=0 && row<rows )
	{
		FreePointer( _entries[row] );
		if( count>0 )
		{
			_entries[ row ] = AllocPointer< MatrixEntry< T , IndexType > >( count );
			memset( _entries[ row ] , 0 , sizeof( MatrixEntry< T , IndexType > )*count );
		}
		rowSizes[row] = count;
	}
	else Miscellany::ErrorOut( "Row is out of bounds: 0 <= %d < %d" , (int)row , (int)rows );
}
template< class T , class IndexType >
void SparseMatrix< T , IndexType >::ResetRowSize( size_t row , size_t count )
{
	if( row>=0 && row<rows )
	{
		size_t oldCount = rowSizes[row];
		_entries[row] = ReAllocPointer< MatrixEntry< T, IndexType > >( _entries[row] , count );
		if( count>oldCount ) memset( _entries[row]+oldCount , 0 , sizeof( MatrixEntry< T , IndexType > ) * ( count - oldCount ) );
		rowSizes[row] = count;
	}
	else Miscellany::ErrorOut( "Row is out of bounds: 0 <= %d < %d" , (int)row , (int)rows );
}

/////////////////////
template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::Identity( size_t dim )
{
	SparseMatrix I;
	I.resize( dim );
	for( int i=0 ; i<dim ; i++ ) I.SetRowSize( i , 1 ) , I[i][0] = MatrixEntry< T , IndexType >( (IndexType)i , (T)1 );
	return I;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::operator *= ( T s )
{
#pragma omp parallel for
	for( int i=0 ; i<rows ; i++ ) for( int j=0 ; j<rowSizes[i] ; j++ ) _entries[i][j].Value *= s;
	return *this;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::operator /= ( T s ){ return (*this) * ( (T)1./s ); }
template< class T , class IndexType >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::operator *= ( const SparseMatrix< T , IndexType >& B )
{
	SparseMatrix foo = (*this) * B;
	(*this) = foo;
	return *this;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::operator += ( const SparseMatrix< T , IndexType >& B )
{
	SparseMatrix foo = (*this) + B;
	(*this) = foo;
	return *this;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType >& SparseMatrix< T , IndexType >::operator -= ( const SparseMatrix< T , IndexType >& B )
{
	SparseMatrix foo = (*this) - B;
	(*this) = foo;
	return *this;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::operator * ( T s ) const
{
	SparseMatrix out = (*this);
	return out *= s;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::operator / ( T s ) const { return (*this) * ( (T)1. / s ); }
template< class T , class IndexType >
Pointer( T ) SparseMatrix< T , IndexType >::operator * ( const Pointer( T ) in ) const
{
	Pointer( T ) out = AllocPointer< T >( rows );
	MultiplyParallel( in , out , omp_get_max_threads() , 0 );
	return out;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::operator * ( const SparseMatrix< T , IndexType >& B ) const
{
	SparseMatrix out;
	const SparseMatrix& A = *this;
	size_t aCols = 0 , aRows = A.rows;
	size_t bCols = 0 , bRows = B.rows;
	for( int i=0 ; i<A.rows ; i++ ) for( int j=0 ; j<A.rowSizes[i] ; j++ ) if( aCols<=A[i][j].N ) aCols = A[i][j].N+1;
	for( int i=0 ; i<B.rows ; i++ ) for( int j=0 ; j<B.rowSizes[i] ; j++ ) if( bCols<=B[i][j].N ) bCols = B[i][j].N+1;
	if( bRows<aCols ) Miscellany::ErrorOut( "Matrix sizes do not support multiplication %lld x %lld * %lld x %lld" , (unsigned long long)aRows , (unsigned long long)aCols , (unsigned long long)bRows , (unsigned long long)bCols );

	out.resize( (int)aRows );
#pragma omp parallel for
	for( int i=0 ; i<aRows ; i++ )
	{
		std::unordered_map< IndexType , T > row;
		for( int j=0 ; j<A.rowSizes[i] ; j++ )
		{
			IndexType idx1 = A[i][j].N;
			T AValue = A[i][j].Value;
			for( int k=0 ; k<B.rowSizes[idx1] ; k++ )
			{
				IndexType idx2 = B[idx1][k].N;
				T BValue = B[idx1][k].Value;
				typename std::unordered_map< IndexType , T >::iterator iter = row.find(idx2);
				if( iter==row.end() ) row[idx2] = AValue * BValue;
				else iter->second += AValue * BValue;
			}
		}
		out.SetRowSize( i , (int)row.size() );
		out.rowSizes[i] = 0;
		for( typename std::unordered_map< IndexType , T >::iterator iter=row.begin() ; iter!=row.end() ; iter++ ) out[i][ out.rowSizes[i]++ ] = MatrixEntry< T , IndexType >( iter->first , iter->second );
	}
	return out;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::operator + ( const SparseMatrix< T , IndexType >& B ) const
{
	const SparseMatrix& A = *this;
	SparseMatrix out;
	out.resize( std::max< size_t >( A.rows , B.rows ) );
#pragma omp parallel for
	for( int i=0 ; i<rows ; i++ )
	{
		std::unordered_map< IndexType , T > row;
		if( i<A.rows )
			for( int j=0 ; j<A.rowSizes[i] ; j++ )
			{
				IndexType idx = A[i][j].N;
				typename std::unordered_map< IndexType , T >::iterator iter = row.find(idx);
				if( iter==row.end() ) row[idx] = A[i][j].Value;
				else iter->second += A[i][j].Value;
			}
		if( i<B.rows )
			for( int j=0 ; j<B.rowSizes[i] ; j++ )
			{
				IndexType idx = B[i][j].N;
				typename std::unordered_map< IndexType , T >::iterator iter = row.find(idx);
				if( iter==row.end() ) row[idx] = B[i][j].Value;
				else iter->second += B[i][j].Value;
			}
		out.SetRowSize( i , row.size() );
		out.rowSizes[i] = 0;
		for( typename std::unordered_map< IndexType , T >::iterator iter=row.begin() ; iter!=row.end() ; iter++ ) out[i][ out.rowSizes[i]++ ] = MatrixEntry< T , IndexType >( iter->first , iter->second );
	}
	return out;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::operator - ( const SparseMatrix< T , IndexType >& B ) const
{
	const SparseMatrix& A = *this;
	size_t rows = std::max< size_t >( A.rows , B.rows );
	SparseMatrix out;

	out.resize( rows );
#pragma omp parallel for
	for( int i=0 ; i<rows ; i++ )
	{
		std::unordered_map< IndexType , T > row;
		if( i<A.rows )
			for( int j=0 ; j<A.rowSizes[i] ; j++ )
			{
				IndexType idx = A[i][j].N;
				typename std::unordered_map< IndexType , T >::iterator iter = row.find(idx);
				if( iter==row.end() ) row[idx] = A[i][j].Value;
				else iter->second += A[i][j].Value;
			}
		if( i<B.rows )
			for( int j=0 ; j<B.rowSizes[i] ; j++ )
			{
				IndexType idx = B[i][j].N;
				typename std::unordered_map< IndexType , T >::iterator iter = row.find(idx);
				if( iter==row.end() ) row[idx] = -B[i][j].Value;
				else iter->second -= B[i][j].Value;
			}
		out.SetRowSize( i , (int)row.size() );
		out.rowSizes[i] = 0;
		for( typename std::unordered_map< IndexType , T >::iterator iter=row.begin() ; iter!=row.end() ; iter++ ) out[i][ out.rowSizes[i]++ ] = MatrixEntry< T , IndexType >( iter->first , iter->second );
	}
	return out;
}

template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::transpose( T (*TransposeFunction)( const T& ) ) const
{
	SparseMatrix A;
	const SparseMatrix& At = *this;
	size_t aRows = 0;
	for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ ) if( aRows<=At[i][j].N ) aRows = At[i][j].N+1;

	A.resize( aRows );
	for( int i=0 ; i<aRows ; i++ ) A.rowSizes[i] = 0;
#pragma omp parallel for
	for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ )
#pragma omp atomic
		A.rowSizes[ At[i][j].N ]++;
#pragma omp parallel for
	for( int i=0 ; i<A.rows ; i++ )
	{
		size_t t = A.rowSizes[i];
		A.rowSizes[i] = 0;
		A.SetRowSize( i , t );
		A.rowSizes[i] = 0;
	}
	if( TransposeFunction ) for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ )
	{
		int ii = At[i][j].N;
		A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , IndexType >( i , TransposeFunction( At[i][j].Value ) );
	}
	else for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ )
	{
		int ii = At[i][j].N;
		A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , IndexType >( i , At[i][j].Value );
	}
	return A;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::transpose( size_t aRows , T (*TransposeFunction)( const T& ) ) const
{
	SparseMatrix A;
	const SparseMatrix& At = *this;
	size_t _aRows = 0 , aCols = At.rows;
	for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ ) if( aCols<=At[i][j].N ) _aRows = At[i][j].N+1;
	if( _aRows>aRows ) Miscellany::ErrorOut( "prescribed output dimension too low: %d < %zu" , aRows , _aRows );

	A.resize( aRows );
	for( int i=0 ; i<aRows ; i++ ) A.rowSizes[i] = 0;
#pragma omp parallel for
	for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ )
#pragma omp atomic
		A.rowSizes[ At[i][j].N ]++;
#pragma omp parallel for
	for( int i=0 ; i<A.rows ; i++ )
	{
		size_t t = A.rowSizes[i];
		A.rowSizes[i] = 0;
		A.SetRowSize( i , t );
		A.rowSizes[i] = 0;
	}
	if( TransposeFunction )
		for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ )
		{
			int ii = At[i][j].N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , IndexType >( i , TransposeFunction( At[i][j].Value ) );
		}
	else
		for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ )
		{
			int ii = At[i][j].N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , IndexType >( i , At[i][j].Value );
		}
	return A;
}
/////////////////////

// Given matrices At and B, compute A * B.
template< class T1 , class T2 , class T3 , class IndexType >
bool TransposeMultiply( const SparseMatrix< T1 , IndexType >& At , const SparseMatrix< T2 , IndexType >& B , SparseMatrix< T3 , IndexType >& out , T1 (*TransposeFunction)( const T1& ) )
{
	int aRows = 0;
	for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ ) if( aRows<=At[i][j].N ) aRows = At[i][j].N+1;
	return TransposeMultiply( At , B , out , aRows , TransposeFunction );
}
template< class T1 , class T2 , class T3 , class IndexType >
bool TransposeMultiply( const SparseMatrix< T1 , IndexType >& At , const SparseMatrix< T2 , IndexType >& B , SparseMatrix< T3 , IndexType >& out , int aRows , T1 (*TransposeFunction)( const T1& ) )
{
	int _aRows = 0     , aCols = At.rows;
	int bRows = B.rows , bCols = 0;
	for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ ) if( _aRows<=At[i][j].N ) _aRows = At[i][j].N+1;
	for( int i=0 ; i< B.rows ; i++ ) for( int j=0 ; j< B.rowSizes[i] ; j++ ) if(  bCols<= B[i][j].N )  bCols =  B[i][j].N+1;
	if( _aRows>aRows ) Miscellany::ErrorOut( "prescribed output dimension too low: %d < %d" , aRows , _aRows );
	if( bCols>At.rows ) Miscellany::ErrorOut( "Matrix sizes do not support multiplication %d x %d * %d x %d" , aRows , aCols , bRows , bCols );

	std::vector< std::unordered_map< IndexType , T3 > > rows;
	rows.resize( _aRows );
	for( int i=0 ; i<bRows ; i++ )
		for( int j=0 ; j<B.rowSizes[i] ; j++ )
		{
			IndexType idx1 = B[i][j].N;
			T2 BValue = B[i][j].Value;
			if( idx1<0 ) continue;
			for( int k=0 ; k<At.rowSizes[i] ; k++ )
			{
				T3 temp;
				IndexType idx2 = At[i][k].N;
				T1 AValue = At[i][k].Value;
				if( TransposeFunction ) temp = T3( TransposeFunction( AValue ) * BValue );
				else                    temp = T3(                    AValue   * BValue ); // temp = At( idx2 , idx1 ) * B( i , idx1 ) = A( idx1 , idx2 ) * B( i , idx1 )
				typename std::unordered_map< IndexType , T3 >::iterator iter = rows[idx2].find(idx1);
				if( iter==rows[idx2].end() ) rows[idx2][idx1] = temp;
				else iter->second += temp;
			}
		}

	out.resize( aRows );
#pragma omp parallel for
	for( int i=0 ; i<rows.size() ; i++ )
	{
		out.SetRowSize( i , rows[i].size() );
		out.rowSizes[i] = 0;
		for( typename std::unordered_map< IndexType , T3 >::iterator iter=rows[i].begin() ; iter!=rows[i].end() ; iter++ )
			out[i][ out.rowSizes[i]++ ] = MatrixEntry< T3 , IndexType >( iter->first , iter->second );
	}
	return true;
}
#if MATRIX_MULTIPLY_INTERFACE
template< class A_T , class A_const_iterator , class B_T , class B_const_iterator , class Out_T , class Out_IndexType >
bool Multiply( const SparseMatrixInterface< A_T , A_const_iterator >& A , const SparseMatrixInterface< B_T , B_const_iterator >& B , SparseMatrix< Out_T , Out_IndexType >& out , unsigned int threads )
{
	size_t aCols = 0 , aRows = A.Rows();
	size_t bCols = 0 , bRows = B.Rows();
	for( int i=0 ; i<A.Rows() ; i++ ) for( A_const_iterator iter=A.begin(i) ; iter!=A.end(i) ; iter++ ) if( aCols<=iter->N ) aCols = iter->N+1;
	for( int i=0 ; i<B.Rows() ; i++ ) for( B_const_iterator iter=B.begin(i) ; iter!=B.end(i) ; iter++ ) if( bCols<=iter->N ) bCols = iter->N+1;
	if( bRows<aCols ) Miscellany::ErrorOut( "Matrix sizes do not support multiplication %lld x %lld * %lld x %lld" , (unsigned long long)aRows , (unsigned long long)aCols , (unsigned long long)bRows , (unsigned long long)bCols );

	out.resize( (int)aRows );
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<aRows ; i++ )
	{
		std::unordered_map< Out_IndexType , Out_T > row;
		for( A_const_iterator iterA=A.begin(i) ; iterA!=A.end(i) ; iterA++ )
		{
			Out_IndexType idx1 = iterA->N;
if( idx1==-1 ) continue;
			A_T AValue = iterA->Value;
			if( idx1<0 ) continue;
			for( B_const_iterator iterB=B.begin(idx1) ; iterB!=B.end(idx1) ; iterB++ )
			{
				Out_IndexType idx2 = iterB->N;
if( idx2==-1 ) continue;
				B_T BValue = iterB->Value;
				Out_T temp = Out_T( BValue * AValue ); // temp = A( i , idx1 ) * B( idx1 , idx2 )
				typename std::unordered_map< Out_IndexType , Out_T >::iterator iter = row.find(idx2);
				if( iter==row.end() ) row[idx2] = temp;
				else iter->second += temp;
			}
		}
		out.SetRowSize( i , (int)row.size() );
		out.rowSizes[i] = 0;
		for( typename std::unordered_map< Out_IndexType , Out_T >::iterator iter=row.begin() ; iter!=row.end() ; iter++ )
			out[i][ out.rowSizes[i]++ ] = MatrixEntry< Out_T , Out_IndexType >( iter->first , iter->second );
	}
	return true;
}
template< class T , class In_const_iterator , class Out_IndexType >
bool Transpose( const SparseMatrixInterface< T , In_const_iterator >& At , SparseMatrix< T , Out_IndexType >& A , T (*TransposeFunction)( const T& ) )
{
	int aRows = 0 , aCols = (int)At.Rows();
	for( int i=0 ; i<At.Rows() ; i++ ) for( In_const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ ) if( aRows<=iter->N ) aRows = iter->N+1;

	A.resize( aRows );
	for( int i=0 ; i<aRows ; i++ ) A.rowSizes[i] = 0;
	for( int i=0 ; i<At.Rows() ; i++ ) for( In_const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ ) A.rowSizes[ iter->N ]++;
	for( int i=0 ; i<A.rows ; i++ )
	{
		int t = A.rowSizes[i];
		A.rowSizes[i] = 0;
		A.SetRowSize( i , t );
		A.rowSizes[i] = 0;
	}
	if( TransposeFunction )
		for( int i=0 ; i<At.Rows() ; i++ ) for( In_const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ )
		{
			int ii = iter->N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , Out_IndexType >( i , TransposeFunction( iter->Value ) );
		}
	else
		for( int i=0 ; i<At.Rows() ; i++ ) for( In_const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ )
		{
			int ii = iter->N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , Out_IndexType >( i , iter->Value );
		}
	return true;
}
template< class T , class In_const_iterator , class Out_IndexType >
bool Transpose( const SparseMatrixInterface< T , In_const_iterator >& At , SparseMatrix< T , Out_IndexType >& A , int aRows , T (*TransposeFunction)( const T& ) )
{
	size_t _aRows = 0 , aCols = At.Rows();
	for( int i=0 ; i<At.Rows() ; i++ ) for( In_const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ ) if( aCols<=iter->N ) _aRows = iter->N+1;
	if( _aRows>aRows ) Miscellany::ErrorOut( "prescribed output dimension too low: %d < %zu" , aRows , _aRows );

	A.resize( aRows );
	for( int i=0 ; i<aRows ; i++ ) A.rowSizes[i] = 0;
	for( int i=0 ; i<At.Rows() ; i++ ) for( In_const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ ) A.rowSizes[ iter->N ]++;
	for( int i=0 ; i<A.rows ; i++ )
	{
		int t = A.rowSizes[i];
		A.rowSizes[i] = 0;
		A.SetRowSize( i , t );
		A.rowSizes[i] = 0;
	}
	if( TransposeFunction )
		for( int i=0 ; i<At.Rows() ; i++ ) for( In_const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ )
		{
			int ii = iter->N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , Out_IndexType >( i , TransposeFunction( iter->Value ) );
		}
	else
		for( int i=0 ; i<At.Rows() ; i++ ) for( In_const_iterator iter=At.begin(i) ; iter!=At.end(i) ; iter++ )
		{
			int ii = iter->N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , Out_IndexType >( i , iter->Value );
		}
	return true;
}

#else // !MATRIX_MULTIPLY_INTERFACE
template< class T1 , class T2 , class T3 , class IndexType >
bool Multiply( const SparseMatrix< T1 , IndexType >& A , const SparseMatrix< T2 , IndexType >& B , SparseMatrix< T3 , IndexType >& out , unsigned int threads )
{
	int aCols = 0 , aRows = A.rows;
	int bCols = 0 , bRows = B.rows;
	for( int i=0 ; i<A.rows ; i++ ) for( int j=0 ; j<A.rowSizes[i] ; j++ ) if( aCols<=A[i][j].N ) aCols = A[i][j].N+1;
	for( int i=0 ; i<B.rows ; i++ ) for( int j=0 ; j<B.rowSizes[i] ; j++ ) if( bCols<=B[i][j].N ) bCols = B[i][j].N+1;
	if( bRows<aCols ) Miscellany::ErrorOut( "Matrix sizes do not support multiplication %d x %d * %d x %d" , aRows , aCols , bRows , bCols );

	out.resize( aRows );
#pragma omp parallel for num_threads( threads )
	for( int i=0 ; i<aRows ; i++ )
	{
		std::unordered_map< IndexType , T3 > row;
		for( int j=0 ; j<A.rowSizes[i] ; j++ )
		{
			IndexType idx1 = A[i][j].N;
			T1 AValue = A[i][j].Value;
			if( idx1<0 ) continue;
			for( int k=0 ; k<B.rowSizes[idx1] ; k++ )
			{
				IndexType idx2 = B[idx1][k].N;
				T2 BValue = B[idx1][k].Value;
				T3 temp = T3( BValue * AValue ); // temp = A( i , idx1 ) * B( idx1 , idx2 )
				typename std::unordered_map< IndexType , T3 >::iterator iter = row.find(idx2);
				if( iter==row.end() ) row[idx2] = temp;
				else iter->second += temp;
			}
		}
		out.SetRowSize( i , row.size() );
		out.rowSizes[i] = 0;
		for( typename std::unordered_map< IndexType , T3 >::iterator iter=row.begin() ; iter!=row.end() ; iter++ )
			out[i][ out.rowSizes[i]++ ] = MatrixEntry< T3 , IndexType >( iter->first , iter->second );
	}
	return true;
}
template< class T , class IndexType >
bool Transpose( const SparseMatrix< T , IndexType >& At , SparseMatrix< T , IndexType >& A , T (*TransposeFunction)( const T& ) )
{
	int aRows = 0 , aCols = At.rows;
	for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ ) if( aRows<=At[i][j].N ) aRows = At[i][j].N+1;


	A.resize( aRows );
	for( int i=0 ; i<aRows ; i++ ) A.rowSizes[i] = 0;
	for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ ) A.rowSizes[ At[i][j].N ]++;
	for( int i=0 ; i<A.rows ; i++ )
	{
		int t = A.rowSizes[i];
		A.rowSizes[i] = 0;
		A.SetRowSize( i , t );
		A.rowSizes[i] = 0;
	}
	if( TransposeFunction )
		for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ ) 
		{
			int ii = At[i][j].N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , IndexType >( i , TransposeFunction( At[i][j].Value ) );
		}
	else
		for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ ) 
		{
			int ii = At[i][j].N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , IndexType >( i , At[i][j].Value );
		}
	return true;
}
template< class T , class IndexType >
bool Transpose( const SparseMatrix< T , IndexType >& At , SparseMatrix< T , IndexType >& A , int aRows , T (*TransposeFunction)( const T& ) )
{
	int _aRows = 0 , aCols = At.rows;
	for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ ) if( _aRows<=At[i][j].N ) _aRows = At[i][j].N+1;
	if( _aRows>aRows ) Miscellany::ErrorOut( "prescribed output dimension too low: %d < %d" , aRows , _aRows );

	A.resize( aRows );
	for( int i=0 ; i<aRows ; i++ ) A.rowSizes[i] = 0;
	for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ ) A.rowSizes[ At[i][j].N ]++;
	for( int i=0 ; i<A.rows ; i++ )
	{
		int t = A.rowSizes[i];
		A.rowSizes[i] = 0;
		A.SetRowSize( i , t );
		A.rowSizes[i] = 0;
	}
	if( TransposeFunction )
		for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ ) 
		{
			int ii = At[i][j].N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , IndexType >( i , TransposeFunction( At[i][j].Value ) );
		}
	else
		for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ ) 
		{
			int ii = At[i][j].N;
			A[ii][ A.rowSizes[ii]++ ] = MatrixEntry< T , IndexType >( i , At[i][j].Value );
		}
	return true;
}
#endif // MATRIX_MULTIPLY_INTERFACE



#define DUMP_OUTPUT 0
inline std::complex< double > operator + ( std::complex< double > c1 , std::complex< float > c2 )
{
	return std::complex< double >( c1.real()+c2.real() , c1.imag()+c2.imag() );
}
inline std::complex< double > operator + ( std::complex< float > c1 , std::complex< double > c2 )
{
	return std::complex< double >( c1.real()+c2.real() , c1.imag()+c2.imag() );
}
inline std::complex< double > operator * ( std::complex< double > c1 , std::complex< float > c2 )
{
	return std::complex< double >( c1.real()*c2.real() - c1.imag()*c2.imag() , c1.real()*c2.imag() + c1.imag()*c2.real() );
}
inline std::complex< double > operator * ( std::complex< float > c1 , std::complex< double > c2 )
{
	return std::complex< double >( c1.real()*c2.real() - c1.imag()*c2.imag() , c1.real()*c2.imag() + c1.imag()*c2.real() );
}
inline std::complex< double > operator * ( std::complex< double > c1 , float c2 )
{
	return std::complex< double >( c1.real()*c2 , c1.imag()*c2 );
}
inline std::complex< double > operator * ( std::complex< float > c1 , double c2 )
{
	return std::complex< double >( c1.real()*c2 , c1.imag()*c2 );
}


//////////////////
// BandedMatrix //
//////////////////

template< class T , unsigned int Radius > BandedMatrix< T , Radius >::~BandedMatrix( void ){ _rows = 0 ; FreePointer( _entries ); }
template< class T , unsigned int Radius > BandedMatrix< T , Radius >::BandedMatrix( void ){ _rows = 0 , _entries = NullPointer< T >(); }
template< class T , unsigned int Radius > BandedMatrix< T , Radius >::BandedMatrix( size_t rows ){ _rows = 0 , _entries = NullPointer< T >() ; resize( rows ); }
template< class T , unsigned int Radius > BandedMatrix< T , Radius >::BandedMatrix( size_t rows , const T& clearValue ){ _rows = 0 , _entries = NullPointer< T >() ; resize( rows , clearValue ); }
template< class T , unsigned int Radius >
template< class T2 >
BandedMatrix< T , Radius >::BandedMatrix( const BandedMatrix< T2 , Radius >& M )
{
	_rows = 0 ; _entries = NullPointer< T >();
	resize( M._rows );
	for( size_t i=0 ; i<entries() ; i++ ) _entries[i] = (T)( M._entries[i] );
}
template< class T , unsigned int Radius >
template< class T2 >
BandedMatrix< T , Radius >& BandedMatrix< T , Radius >::operator = ( const BandedMatrix< T2 , Radius >& M )
{
	resize( M._rows );
	for( size_t i=0 ; i<entries() ; i++ ) _entries[i] = (T)( M._entries[i] );
	return *this;
}

template< class T , unsigned int Radius > void BandedMatrix< T , Radius >::resize( size_t rows )
{
	if( rows==_rows ) return;
	FreePointer( _entries );
	_rows = 0;
	if( rows )
	{
		_rows = rows;
		_entries = AllocPointer< T >( rows * ( 2 * Radius + 1 ) );
		if( !_entries ) Miscellany::ErrorOut( "Failed to allocate BandedMatrix::_entries ( %d x %d )" , rows , 2*Radius+1 );
	}
}
template< class T , unsigned int Radius > void BandedMatrix< T , Radius >::resize( size_t rows , const T& clearValue )
{
	resize( rows );
	for( size_t i=0 ; i<entries() ; i++ ) _entries[i] = clearValue;
}
template< class T , unsigned int Radius >
template< class T2 >
void BandedMatrix< T , Radius >::multiply( ConstPointer( T2 ) in , Pointer( T2 ) out , unsigned int threads ) const
{
	for( int i=0 ; i<Radius && i<_rows-Radius ; i++ )
	{
		T2 sum(0);
		const T* __entries = _entries + i * ( 2 * Radius + 1 );
		size_t ii = i + _rows - Radius;
		for( int j=0 ; j<=2*Radius ; j++ ) sum += (T2)( in[(ii+j)%_rows] * __entries[j] );
		out[i] = sum;
	}
	if( Radius==1 )
	{
#pragma omp parallel for num_threads( threads )
		for( int i=1 ; i<_rows-1 ; i++ )
		{
			ConstPointer( T ) __entries = _entries + i * 3;
			ConstPointer( T2 ) _in = in + i - 1;
			out[i] = (T2)( _in[0] * __entries[0] + _in[1] * __entries[1] + _in[2] * __entries[2] );
		}
	}
	else
	{
#pragma omp parallel for num_threads( threads )
		for( int i=Radius ; i<_rows-Radius ; i++ )
		{
			T2 sum(0);
			ConstPointer( T ) __entries = _entries + i * ( 2 * Radius + 1 );
			ConstPointer( T2 ) _in = in + i - Radius;
			for( int j=0 ; j<=2*Radius ; j++ ) sum += (T2)( _in[j] * __entries[j] );
			out[i] = sum;
		}
	}
	for( int i=(int)_rows-Radius ; i<_rows ; i++ )
	{
		T2 sum(0);
		const T* __entries = _entries + i * ( 2 * Radius + 1 );
		int ii = (int)( i + _rows - Radius );
		for( int j=0 ; j<=2*Radius ; j++ ) sum += (T2)( in[(ii+j)%_rows] * (T2)__entries[j] );
		out[i] = sum;
	}
}
template< class T , unsigned int Radius >
template< class T2 >
void BandedMatrix< T , Radius >::multiply2( ConstPointer( T2 ) in , Pointer( T2 ) out , unsigned int threads ) const
{
	for( int i=0 ; i<Radius && i<_rows-Radius ; i++ )
	{
		T2 sum0(0) , sum1(0);
		const T* __entries = _entries + i * ( 2 * Radius + 1 );
		size_t ii = i + _rows - Radius;
		for( int j=0 ; j<=2*Radius ; j++ )
		{
			int iii = (int)( (ii + j)%_rows );
			sum0 += (T2)( in[ iii<<1   ] * __entries[j] );
			sum1 += (T2)( in[(iii<<1)|1] * __entries[j] );
		}
		out[ i<<1   ] = sum0;
		out[(i<<1)|1] = sum1;
	}
	if( Radius==1 )
	{
#pragma omp parallel for num_threads( threads )
		for( int i=1 ; i<_rows-1 ; i++ )
		{
			ConstPointer( T ) __entries = _entries + i * 3;
			ConstPointer( T2 ) _in = in + (i-1)*2;
			out[ i<<1   ] = (T2)( _in[0] * __entries[0] + _in[2] * __entries[1] + _in[4] * __entries[2] );
			out[(i<<1)|1] = (T2)( _in[1] * __entries[0] + _in[3] * __entries[1] + _in[5] * __entries[2] );
		}
	}
	else
	{
#pragma omp parallel for num_threads( threads )
		for( int i=Radius ; i<_rows-Radius ; i++ )
		{
			T2 sum0(0) , sum1(0);
			ConstPointer( T ) __entries = _entries + i * ( 2 * Radius + 1 );
			ConstPointer( T2 ) _in = in + (i-Radius)*2;
			for( int j=0 ; j<=2*Radius ; j++ ) sum0 += (T2)( _in[j<<1] * __entries[j] ) ,  sum1 += (T2)( _in[(j<<1)|1] * __entries[j] );
			out[ i<<1   ] = sum0;
			out[(i<<1)|1] = sum1;
		}
	}
	for( int i=(int)_rows-Radius ; i<_rows ; i++ )
	{
		T2 sum0(0) , sum1(0);
		const T* __entries = _entries + i * ( 2 * Radius + 1 );
		int ii = i + (int)(_rows-Radius);
		for( int j=0 ; j<=2*Radius ; j++ )
		{
			int iii = (ii+j)%_rows;
			sum0 += (T2)( in[ iii<<1   ] * __entries[j] );
			sum1 += (T2)( in[(iii<<1)|1] * __entries[j] );
		}
		out[ i<<1   ] = sum0;
		out[(i<<1)|1] = sum1;
	}
}
template< class T , unsigned int Radius >
double BandedMatrix< T , Radius >::squareNorm( void ) const
{
	double n2 = 0;
	for( int i=0 ; i<entries() ; i++ ) n2 += _entries[i] * _entries[i];
	return n2;
}
#include <map>
template< class RealIn, class RealOut>
void CompressSparseMatrix(const  SparseMatrix< RealIn, int > & M, SparseMatrix< RealOut, int > & _M) {
	int nrows = (int)M.Rows();
	_M.resize(nrows);

	for (int i = 0; i < nrows; i++){
		//Compact Row
		std::map<int, double> rowValues;
		rowValues[i] = 0; //Ensure there is a diagonal element
		for (int j = 0; j < M.RowSize(i); j++){
			RealIn value = M[i][j].Value;
			int row = M[i][j].N;
			auto mapIter = rowValues.find(row);
			if (mapIter == rowValues.end()) {
				rowValues[row] = value;
			}
			else {
				mapIter->second += value;
			}
		}
		//Write Row
		int numEntries = (int)rowValues.size();
		_M.SetRowSize(i, numEntries);
		int offset = 0;
		for (auto iter = rowValues.begin(); iter != rowValues.end(); iter++) {
			_M[i][offset] = MatrixEntry< RealOut, int >( iter->first , (RealOut)iter->second );
			offset++;
		}
	}
}
