/*
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

#ifndef __SPARSEMATRIX_HPP
#define __SPARSEMATRIX_HPP
#define MATRIX_MULTIPLY_INTERFACE 1

#include "SparseMatrixInterface.h"
#include "Array.h"

template< class T , class IndexType > class SparseMatrix : public SparseMatrixInterface< T , ConstPointer( MatrixEntry< T , IndexType > ) >
{
	template< class T2 , class IndexType2 > friend class SparseMatrix;
	Pointer( Pointer( MatrixEntry< T , IndexType > ) ) _entries;
public:
	static void Swap( SparseMatrix& M1 , SparseMatrix& M2 )
	{
		std::swap( M1.rows , M2.rows );
		std::swap( M1.rowSizes , M2.rowSizes );
		std::swap( M1._entries , M2._entries );
	}
	typedef SparseMatrixInterface< T , ConstPointer( MatrixEntry< T , IndexType > ) > Interface;
	typedef ConstPointer( MatrixEntry< T , IndexType > ) RowIterator;

	size_t rows;
	Pointer( size_t ) rowSizes;

	SparseMatrix( void );
	SparseMatrix( const SparseMatrix& M );
	SparseMatrix( SparseMatrix&& M );
	template< class T2 , class IndexType2 >
	SparseMatrix( const SparseMatrix< T2 , IndexType2 >& M );
	~SparseMatrix();
	SparseMatrix& operator = ( SparseMatrix&& M );
	SparseMatrix< T , IndexType >& operator = ( const SparseMatrix< T , IndexType >& M );
	template< class T2 , class IndexType2 >
	SparseMatrix< T , IndexType >& operator = ( const SparseMatrix< T2 , IndexType2 >& M );

	template< class T2 > void operator()( const T2* in , T2* out ) const;

	template< class T2 , class IndexType2 >
	SparseMatrix< T , IndexType >& copy( const SparseMatrix< T2 , IndexType2 >& M );

	inline ConstPointer( MatrixEntry< T , IndexType > ) begin( size_t row ) const { return _entries[row]; }
	inline ConstPointer( MatrixEntry< T , IndexType > ) end  ( size_t row ) const { return _entries[row] + rowSizes[row]; }
	inline size_t Rows                              ( void )       const { return rows; }
	inline size_t RowSize                           ( size_t idx ) const { return rowSizes[idx]; }

	SparseMatrix( size_t rows );
	void resize	( size_t rows );
	void SetRowSize( size_t row , size_t count );
	void ResetRowSize( size_t row , size_t count );
	inline      Pointer( MatrixEntry< T , IndexType > ) operator[] ( size_t idx )       { return _entries[idx]; }
	inline ConstPointer( MatrixEntry< T , IndexType > ) operator[] ( size_t idx ) const { return _entries[idx]; }

	double SquareNorm(void) const;
	double ASymmetricSquareNorm( void ) const;
	double AHermitianSquareNorm( void ) const;

	/** Sets the column index of all allocated entries to -1 so that they are
	*  treated as non-existent. This is needed because SetRowSize() uses
	*  malloc instead of new and MatrixEntry's constructor is never called. */
	void invalidateEntries();

	/** Adds a scalar value to an element in the Matrix, using a new element if
	*  necessary. If no pre-allocated space for a new element exists, false is
	*  returned.
	*  WARNING: no check is done to remove entries that become zero after
	*  addition */
	bool addScalarToEntry( T s , IndexType i , IndexType j );
	
	// With copy move, these should be well-behaved from a memory perspective
	static SparseMatrix Identity( size_t dim );
	SparseMatrix transpose(                  T (*TransposeFunction)( const T& )=NULL ) const;
	SparseMatrix transpose( size_t outRows , T (*TransposeFunction)( const T& )=NULL ) const;
	SparseMatrix  operator *  ( T s ) const;
	SparseMatrix  operator /  ( T s ) const;
	SparseMatrix  operator *  ( const SparseMatrix& M ) const;
	SparseMatrix  operator +  ( const SparseMatrix& M ) const;
	SparseMatrix  operator -  ( const SparseMatrix& M ) const;
	SparseMatrix& operator *= ( T s );
	SparseMatrix& operator /= ( T s );
	SparseMatrix& operator *= ( const SparseMatrix& M );
	SparseMatrix& operator += ( const SparseMatrix& M );
	SparseMatrix& operator -= ( const SparseMatrix& M );

	// [WARNING] The pointer needs to be deallocated
	Pointer( T ) operator * ( const Pointer( T ) in ) const;
};
template< class T1 , class T2 , class T3 , class IndexType >
bool TransposeMultiply( const SparseMatrix< T1 , IndexType >& At , const SparseMatrix< T2 , IndexType >& B , SparseMatrix< T3 , IndexType >& out , T1 (*TransposeFunction)( const T1& )=NULL );
template< class T1 , class T2 , class T3 , class IndexType >
bool TransposeMultiply( const SparseMatrix< T1 , IndexType >& At , const SparseMatrix< T2 , IndexType >& B , SparseMatrix< T3 , IndexType >& out , size_t outRows , T1 (*TransposeFunction)( const T1& )=NULL );
#if MATRIX_MULTIPLY_INTERFACE
template< class A_T , class A_const_iterator , class B_T , class B_const_iterator , class Out_T , class Out_IndexType >
bool Multiply( const SparseMatrixInterface< A_T , A_const_iterator >& A , const SparseMatrixInterface< B_T , B_const_iterator >& B , SparseMatrix< Out_T , Out_IndexType >& out , unsigned int threads = 1 );
template< class T , class In_const_iterator , class Out_IndexType >
bool Transpose( const SparseMatrixInterface< T , In_const_iterator >& At , SparseMatrix< T , Out_IndexType >& A ,                  T (*TransposeFunction)( const T& )=NULL );
template< class T , class In_const_iterator , class Out_IndexType >
bool Transpose( const SparseMatrixInterface< T , In_const_iterator >& At , SparseMatrix< T , Out_IndexType >& A , size_t outRows , T (*TransposeFunction)( const T& )=NULL );
#else
template< class T1 , class T2 , class T3 , class IndexType >
bool Multiply(const SparseMatrix< T1, IndexType > &A, const SparseMatrix< T2, IndexType > &B, SparseMatrix< T3, IndexType > &out , unsigned int threads = 1);
template< class T , class IndexType >
bool Transpose( const SparseMatrix< T , IndexType >& At , SparseMatrix< T , IndexType >& A ,                  T (*TransposeFunction)( const T& )=NULL );
template< class T , class IndexType >
bool Transpose( const SparseMatrix< T , IndexType >& At , SparseMatrix< T , IndexType >& A , size_t outRows , T (*TransposeFunction)( const T& )=NULL );
#endif // MATRIX_MULTIPLY_INTERFACE

template< class T , unsigned int Radius > class BandedMatrixIterator
{
	size_t _rows , _row;
	int _offset;
	ConstPointer( T ) _values;
	MatrixEntry< size_t , T > _entry;
public:
	bool operator !=( const BandedMatrixIterator& i ) const { return i._row!=_row || i._offset!=_offset; }
	const MatrixEntry< T , size_t >& operator -> ( void ) const { return _entry; }
	BandedMatrixIterator& operator ++ ( void )
	{
		_offset++;
		if( _offset<=Radius ) _entry = MatrixEntry< T , size_t >( ( _rows + _row + _offset ) % _rows , _values[ Radius+_offset ] );
		return *this;
	}
	BandedMatrixIterator( size_t rows , size_t row , ConstPointer( T ) values ){ _rows = rows , _row = row , _values = values, _offset = -(int)Radius; }
	BandedMatrixIterator( size_t row ){ _row = row , _offset = Radius+1; }
};
template< class T , unsigned int Radius > class BandedMatrix : public SparseMatrixInterface< T , BandedMatrixIterator< T , Radius > >
{
	template< class T2 , unsigned int Radius2 > friend class BandedMatrix;
	size_t _rows;
	Pointer( T ) _entries;
public:
	BandedMatrix( void );
	BandedMatrix( size_t rows );
	BandedMatrix( size_t rows , const T& clearValue );
	template< class T2 > BandedMatrix( const BandedMatrix< T2 , Radius >& M );
	template< class T2 > BandedMatrix& operator = ( const BandedMatrix< T2 , Radius >& M );
	~BandedMatrix();

	BandedMatrixIterator< T , Radius > begin( size_t row ) const { return BandedMatrixIterator< T , Radius >( _rows , row , _entries + row * (2*Radius+1) ); }
	BandedMatrixIterator< T , Radius >   end( size_t row ) const { return BandedMatrixIterator< T , Radius >( row ); }
	size_t Rows( void ) const{ return _rows; }
	size_t RowSize( size_t idx ) const { return 2*Radius+1; }

	template< class T2 > void multiply ( ConstPointer( T2 ) in , Pointer( T2 ) out , unsigned int threads=1 ) const;
	template< class T2 > void multiply2( ConstPointer( T2 ) in , Pointer( T2 ) out , unsigned int threads=1 ) const;

	inline size_t rows( void ) const { return _rows; }
	inline size_t entries( void ) const { return _rows * ( 2 * Radius + 1 ); }

	void resize( size_t rows );
	void resize( size_t rows , const T& clearValue );
	inline      Pointer( T ) operator[] ( size_t idx )       { return _entries + idx * ( 2 * Radius + 1 ); }
	inline ConstPointer( T ) operator[] ( size_t idx ) const { return _entries + idx * ( 2 * Radius + 1 ); }
	double squareNorm( void ) const;
};

#include "SparseMatrix.inl"
#endif /* __SPARSEMATRIX_HPP */
