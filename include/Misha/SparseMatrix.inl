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

template< class T , class IndexType >
void SparseMatrix< T , IndexType >::write( FILE* fp ) const
{
	fwrite( &rows , sizeof(size_t) , 1 , fp );
	fwrite( rowSizes , sizeof(size_t) , rows , fp );
	for( int i=0 ; i<rows ; i++ ) fwrite( _entries[i] , sizeof( MatrixEntry< T , IndexType > ) , rowSizes[i] , fp );
}
template< class T , class IndexType >
void SparseMatrix< T , IndexType >::read( FILE* fp )
{
	size_t _rows;
	fread( &_rows , sizeof(size_t) , 1 , fp );
	resize( _rows );
	for( int i=0 ; i<rows ; i++ )
	{
		size_t _rowSize;
		fread( &_rowSize , sizeof(size_t) , 1 , fp );
		SetRowSize( i , _rowSize );
	}
	for( int i=0 ; i<rows ; i++ ) fread( _entries[i] , sizeof( MatrixEntry< T , IndexType > ) , rowSizes[i] , fp );
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
template< typename SquareNormFunctor >
double SparseMatrix< T , IndexType >::SquareNorm( SquareNormFunctor F ) const
{
	double l2=0;
	for( int i=0 ; i<rows ; i++ ) for( int j=0 ; j<rowSizes[i] ; j++ ) l2 += F( _entries[i][j].Value );
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
#if 1
		// [WARNING] multi-counting 
		for( int k=0 ; k<rowSizes[i] ; k++ ) if( _entries[i][k].N==N ) t1 += _entries[i][k].Value;
		for( int k=0 ; k<rowSizes[N] ; k++ ) if( _entries[N][k].N==i ) t2 += _entries[N][k].Value;
#else
		t1 = _entries[i][j].Value;
		for( int k=0 ; k<rowSizes[N] ; k++ ) if( _entries[N][k].N==i )
		{
			t2 = _entries[N][k].Value;
			break;
		}
#endif
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
Vector< T2 > SparseMatrix< T , IndexType >::operator * ( const Vector< T2 >& V ) const
{
	Vector< T2 > R( Rows() );
	Interface::Multiply( V , R );
	return R;
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
	else fprintf( stderr , "[ERROR] SparseMatrix::SetRowSize: Row is out of bounds: 0 <= %d < %d\n" , (int)row , (int)rows ) , exit( 0 );
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
	else fprintf( stderr , "[ERROR] SparseMatrix::ResetRowSize: Row is out of bounds: 0 <= %d < %d\n" , (int)row , (int)rows ) , exit( 0 );
}
template< class T , class IndexType >
void SparseMatrix< T , IndexType >::CollapseRow( int r )
{
	{
		std::unordered_map< IndexType , T > rowValues;
		for( int j=0 ; j<rowSizes[r] ; j++ )
		{
			int N = _entries[r][j].N;
			auto iter = rowValues.find( N );
			if( iter==rowValues.end() ) rowValues[N]  = _entries[r][j].Value;
			else                        iter->second += _entries[r][j].Value;
		}
		SetRowSize( r , (int)rowValues.size() );
		rowSizes[r] = 0;
		for( auto iter=rowValues.begin() ; iter!=rowValues.end() ; iter++ ) _entries[r][ rowSizes[r]++ ] = MatrixEntry< T , IndexType >( iter->first , iter->second );
	}
}
template< class T , class IndexType >
void SparseMatrix< T , IndexType >::CollapseRows( void )
{
	ThreadPool::ParallelFor
	(
		0 , rows ,
		[&]( unsigned int , size_t i )
		{
			std::unordered_map< IndexType , T > rowValues;
			for( int j=0 ; j<rowSizes[i] ; j++ )
			{
				int N = _entries[i][j].N;
				auto iter = rowValues.find( N );
				if( iter==rowValues.end() ) rowValues[N]  = _entries[i][j].Value;
				else                        iter->second += _entries[i][j].Value;
			}
			SetRowSize( i , (int)rowValues.size() );
			rowSizes[i] = 0;
			for( auto iter=rowValues.begin() ; iter!=rowValues.end() ; iter++ ) _entries[i][ rowSizes[i]++ ] = MatrixEntry< T , IndexType >( iter->first , iter->second );
		}
	);
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
	ThreadPool::ParallelFor( 0 , rows , [&]( unsigned int , size_t i ){ for( int j=0 ; j<rowSizes[i] ; j++ ) _entries[i][j].Value *= s; } );
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
	MultiplyParallel( in , out , ThreadPool::NumThreads() , 0 );
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
	if( bRows<aCols ) fprintf( stderr , "[Error] SparseMatrix::operator *: Matrix sizes do not support multiplication %lld x %lld * %lld x %lld\n" , (unsigned long long)aRows , (unsigned long long)aCols , (unsigned long long)bRows , (unsigned long long)bCols ) , exit( 0 );

	out.resize( (int)aRows );
	ThreadPool::ParallelFor
	(
		0 , aRows ,
		[&]( unsigned int , size_t i )
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
	);
	return out;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::operator + ( const SparseMatrix< T , IndexType >& B ) const
{
	const SparseMatrix& A = *this;
	size_t rows = std::max< size_t >( A.rows , B.rows );
	SparseMatrix out;

	out.resize( rows );
	ThreadPool::ParallelFor
	(
		0 , rows ,
		[&]( unsigned int , size_t i )
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
	);
	return out;
}
template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::operator - ( const SparseMatrix< T , IndexType >& B ) const
{
	const SparseMatrix& A = *this;
	size_t rows = std::max< size_t >( A.rows , B.rows );
	SparseMatrix out;

	out.resize( rows );
	ThreadPool::ParallelFor
	(
		0 , rows ,
		[&]( unsigned int , size_t i )
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
	);
	return out;
}

template< class T , class IndexType >
SparseMatrix< T , IndexType > SparseMatrix< T , IndexType >::transpose( T (*TransposeFunction)( const T& ) ) const
{
	SparseMatrix A;
	const SparseMatrix& At = *this;
	size_t aRows = 0 , aCols = At.rows;
	for( int i=0 ; i<At.rows ; i++ ) for( int j=0 ; j<At.rowSizes[i] ; j++ ) if( aRows<=At[i][j].N ) aRows = At[i][j].N+1;

	A.resize( aRows );
	for( int i=0 ; i<aRows ; i++ ) A.rowSizes[i] = 0;
	ThreadPool::ParallelFor( 0 , At.rows , [&]( unsigned int , size_t i ){ for( int j=0 ; j<At.rowSizes[i] ; j++ ) AddAtomic( A.rowSizes[ At[i][j].N ] , (size_t)1 ); } );

	ThreadPool::ParallelFor
	(
		0 , A.rows ,
		[&]( unsigned int , size_t i )
		{
			size_t t = A.rowSizes[i];
			A.rowSizes[i] = 0;
			A.SetRowSize( i , t );
			A.rowSizes[i] = 0;
		}
	);

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
	if( _aRows>aRows )
	{
		fprintf( stderr , "[Error] SparseMatrix::transpose: prescribed output dimension too low: %d < %d\n" , (int)aRows , (int)_aRows );
		return false;
	}

	A.resize( aRows );
	for( int i=0 ; i<aRows ; i++ ) A.rowSizes[i] = 0;

	ThreadPool::ParallelFor( 0 , At.rows , [&]( unsigned int , size_t i ){ for( int j=0 ; j<At.rowSizes[i] ; j++ ) AddAtomic( A.rowSizes[ At[i][j].N ] , 1 ); } );

	ThreadPool::ParallelFor
	(
		0 , A.rows ,
		[&]( unsigned int , size_t i )
		{
			size_t t = A.rowSizes[i];
			A.rowSizes[i] = 0;
			A.SetRowSize( i , t );
			A.rowSizes[i] = 0;
		}
	);

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
	if( _aRows>aRows )
	{
		fprintf( stderr , "[Error] TransposeMultiply: prescribed output dimension too low: %d < %d\n" , aRows , _aRows );
		return false;
	}
	if( bCols>At.rows )
	{
		fprintf( stderr , "[Error] TransposeMultiply: Matrix sizes do not support multiplication %d x %d * %d x %d\n" , aRows , aCols , bRows , bCols );
		return false;
	}

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
	ThreadPool::ParallelFor
	(
		0 , rows.size() ,
		[&]( unsigned int , size_t i )
		{
			out.SetRowSize( i , rows[i].size() );
			out.rowSizes[i] = 0;
			for( typename std::unordered_map< IndexType , T3 >::iterator iter=rows[i].begin() ; iter!=rows[i].end() ; iter++ )
				out[i][ out.rowSizes[i]++ ] = MatrixEntry< T3 , IndexType >( iter->first , iter->second );
		}
	);
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
	if( bRows<aCols )
	{
		fprintf( stderr , "[Error] Multiply: Matrix sizes do not support multiplication %lld x %lld * %lld x %lld\n" , (unsigned long long)aRows , (unsigned long long)aCols , (unsigned long long)bRows , (unsigned long long)bCols );
		return false;
	}

	out.resize( (int)aRows );
	ThreadPool::ParallelFor
	(
		0 , aRows ,
		[&]( unsigned int , size_t i )
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
		} ,
		threads
	);
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
	if( _aRows>aRows )
	{
		fprintf( stderr , "[Error] Transpose: prescribed output dimension too low: %d < %zu\n" , aRows , _aRows );
		return false;
	}

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
	if( bRows<aCols )
	{
		fprintf( stderr , "[Error] Multiply: Matrix sizes do not support multiplication %d x %d * %d x %d\n" , aRows , aCols , bRows , bCols );
		return false;
	}

	out.resize( aRows );
	ThreadPool::ParallelFor
	(
		0 , aRows ,
		[&]( unsigned int , size_t i )
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
		} ,
		threads
	);
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
	if( _aRows>aRows )
	{
		fprintf( stderr , "[Error] Transpose: prescribed output dimension too low: %d < %d\n" , aRows , _aRows );
		return false;
	}

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
template< class Data >
double Dot( const Vector< Data >& v1 , const Vector< Data >& v2 , double (*dot)( Data , Data ) )
{
	double d = 0;
	for( int i=0 ; i<v1.size() && i<v2.size() ; i++ ) d += dot( v1[i] , v2[i] );
	return d;
}
template< class Data >
std::complex< double > Dot( const Vector< Data >& v1 , const Vector< Data >& v2 , std::complex< double > (*dot)( Data , Data ) )
{
	std::complex< double > d = 0;
	for( int i=0 ; i<v1.size() && i<v2.size() ; i++ ) d += dot( v1[i] , v2[i] );
	return d;
}
template < class Matrix , class Data >
static int SolveConjugateGradient( const Matrix& SPD , const Vector< Data >& b , const int& iters , Vector< Data >& solution , double (*dot)( Data , Data ) , const double eps )
{
	Vector<Data> d,r,Md,temp;
	double alpha,beta,rDotR,oldRDotR;
	Md.resize(b.size());

	temp.resize(b.size());
	SPD.Multiply(solution,temp);
	d = r = b-temp;
	oldRDotR = rDotR = Dot( r , r , dot );
	if( Dot( b , b , dot )<=eps)
	{
		solution.SetZero();
		return 0;
	}
	int i;
	for(i=0;i<iters;i++)
	{
		double temp;
		SPD.Multiply(d,Md);
		temp = Dot( d , Md , dot );
		if( temp<=eps )
		{
			break;
		}
		alpha=rDotR/temp;
		r.SubtractScaled(Md,alpha);
		temp = Dot( r , r , dot );
		// BADNESS!!! How can the size of the residual increase?
		if(temp>2*oldRDotR)
		{
			break;
		}
		oldRDotR=rDotR;
		beta=temp/rDotR;
		solution.AddScaled(d,alpha);
		if(beta<=eps)
		{
			break;
		}
		rDotR=temp;
		Vector<Data>::Add(d,beta,r,d);
	}
	return i;
}
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
template < class Matrix , class Data >
static int SolveConjugateGradient( const Matrix& SPD , const Vector< Data >& b , const int& iters , Vector< Data >& solution , std::complex< double > (*dot)( Data , Data ) , const double eps )
{
	Vector< Data > d , r , Md,temp;
	std::complex< double > alpha;
	double rDotR , oldRDotR , beta;
	Md.resize( b.size() );

	temp.resize( b.size() );
	SPD.Multiply( solution , temp );
	d = r = b-temp;
	oldRDotR = rDotR = Dot( r , r , dot ).real();
	if( Dot( b , b , dot ).real()<=eps )
	{
		solution.SetZero();
		return 0;
	}
	int i;
	for( i=0 ; i<iters ; i++ )
	{
		std::complex< double > temp;
		SPD.Multiply( d , Md );
		temp = Dot( d , Md , dot );
		if( temp.real()*temp.real()+temp.imag()*temp.imag()<=eps*eps ) break;

		alpha = rDotR/temp;
		r.SubtractScaled( Md , alpha );
		double _temp = Dot( r , r , dot ).real();
		if( _temp>2*oldRDotR ) break;

		oldRDotR = rDotR;
		beta = _temp/rDotR;
		solution.AddScaled( d , alpha );
		if( beta<=eps ) break;
		rDotR = _temp;
		Vector< Data >::Add( d , beta , r , d );
	}
	return i;
}
template < class Matrix , class Data >
static int SolveConjugateGradient( const Matrix& SPD , const Vector< Data >& b , const int& iters , Vector< Data >& solution , const double eps )
{
	Vector<Data> d,r,Md,temp;
	double alpha,beta,rDotR,oldRDotR;
	Md.resize(b.size());

	temp.resize(b.size());
	SPD.Multiply(solution,temp);
	d=r=b-temp;
	oldRDotR=rDotR=r.Dot(r);
	if(b.Dot(b)<=eps)
	{
#if DUMP_OUTPUT
		printf("Badness0: %g %g\n",r.Dot(r),eps);
#endif // DUMP_OUTPUT
		solution.SetZero();
		return 0;
	}
	int i;
	for(i=0;i<iters;i++)
	{
		double temp;
		SPD.Multiply(d,Md);
		temp=d.Dot(Md);
		if(temp<=eps)
		{
#if DUMP_OUTPUT
			printf("Badness1: %g %g\n",temp,eps);
#endif // DUMP_OUTPUT
			break;
		}
		alpha=rDotR/temp;
		r.SubtractScaled(Md,alpha);
		temp=r.Dot(r);
		// BADNESS!!! How can the size of the residual increase?
		if(temp>2*oldRDotR)
		{
#if DUMP_OUTPUT
			printf("Badness1.5: %g %g\n",temp,oldRDotR);
#endif // DUMP_OUTPUT
			break;
		}
		oldRDotR=rDotR;
		if(temp/b.Dot(b)<=eps)
		{
#if DUMP_OUTPUT
			//			printf("Badness2: %g %g\n",temp,eps);
#endif // DUMP_OUTPUT
			//			break;
		}
		beta=temp/rDotR;
		solution.AddScaled(d,alpha);
		if(beta<=eps)
		{
#if DUMP_OUTPUT
			printf("Badness3: %g %g\n",beta,eps);
#endif // DUMP_OUTPUT
			break;
		}
		rDotR=temp;
		Vector<Data>::Add(d,beta,r,d);
	}
	return i;
}
template <class Matrix,class IPS,class Real>
static int SolveConjugateGradient(const Matrix& SPD,const Vector<IPS>& b,const int& iters,Vector<IPS>& solution,const double eps)
{
	Vector<IPS> d,r,Md,temp;
	double alpha,beta,rDotR,oldRDotR;
	Md.resize(b.size());

	temp.resize(b.size());
	SPD.Multiply(solution,temp);
	d=r=b-temp;
	oldRDotR=rDotR=r.IPSDot(r);
	if(b.IPSDot(b)<=eps)
	{
#if DUMP_OUTPUT
		printf("Badness0: %g %g\n",r.Dot(r),eps);
#endif // DUMP_OUTPUT
		solution.SetZero();
		return 0;
	}
	int i;
	for(i=0;i<iters;i++)
	{
		double temp;
		SPD.Multiply(d,Md);
		temp=d.IPSDot(Md);
		if(temp<=eps)
		{
#if DUMP_OUTPUT
			printf("Badness1: %g %g\n",temp,eps);
#endif // DUMP_OUTPUT
			break;
		}
		alpha=rDotR/temp;
		r.SubtractScaled(Md,alpha);
		temp=r.IPSDot(r);
		// BADNESS!!! How can the size of the residual increase?
		if(temp>2*oldRDotR)
		{
#if DUMP_OUTPUT
			printf("Badness1.5: %g %g\n",temp,oldRDotR);
#endif // DUMP_OUTPUT
			break;
		}
		oldRDotR=rDotR;
		if(temp/b.IPSDot(b)<=eps)
		{
#if DUMP_OUTPUT
			//			printf("Badness2: %g %g\n",temp,eps);
#endif // DUMP_OUTPUT
			//			break;
		}
		beta=temp/rDotR;
		solution.AddScaled(d,alpha);
		if(beta<=eps)
		{
#if DUMP_OUTPUT
			printf("Badness3: %g %g\n",beta,eps);
#endif // DUMP_OUTPUT
			break;
		}
		rDotR=temp;
		Vector<IPS>::Add(d,beta,r,d);
	}
	return i;
}


template <class Matrix,class IPS,class Real>
static int SolveConjugateGradient2(const Matrix& SPD,const Vector<IPS>& b,const int& iters,Vector<IPS>& x,const double eps)
{
	Vector<IPS> q,d,r;
	double delta_new,delta_old,delta_0,alpha,beta;
	q.resize(b.size());
	SPD.Multiply(x,q);
	d=r=b-q;
	delta_0=delta_new=r.IPSDot(r);
	printf("%f %f\n",x.IPSDot(x),delta_0);
	int i;
	for(i=0;i<iters && delta_new>eps*eps*delta_0;i++)
	{
		SPD.Multiply(d,q);
		alpha=delta_new/(d.IPSDot/*<double>*/(q));
		printf("\t%d] %f\n",i,d.IPSDot(q));
		printf("\t\talpha = %f\n",alpha);
		x.AddScaled(d,alpha);
		if(!(i%50))
		{
			SPD.Multiply(x,q);
			r=b-q;
		}
		else	r.SubtractScaled(q,alpha);
		delta_old=delta_new;
		delta_new=r.IPSDot/*<double>*/(r);
		beta=delta_new/delta_old;
		printf("\t\t beta = %f\n",beta);
		printf("\t\tresid = %f\n",delta_new/delta_0);
		Vector<IPS>::Add(d,beta,r,d);
	}
	printf("Iters: %d / %d\n",i,iters);
	//exit(0);
	return i;
}

// Prevent "duplicate symbol" linker errors by hiding the constants from
// external files.
namespace {
	////////////////////////////////////////////////////////////////////////////
	/*! Declares some default constants for conjugate gradient solvers
	//  @tparam Real    floating point representation
	*///////////////////////////////////////////////////////////////////////////
	template <typename Real>
	class CGConstants
	{
	public:
		/** Constant used for convergence testing */
		static const Real eps;
	};

	template<> const float CGConstants<float>::eps = float(1e-12);
	template<> const double CGConstants<double>::eps = 1e-16;
}

template< class MType , class IndexType , class VType >
int SolveConjugateGradient( const SparseMatrix< MType , IndexType >& A , const Vector<VType>& b , const int& iters , Vector<VType>& x ,
	Vector<VType> (*Multiply)(const SparseMatrix< MType , IndexType >& , const Vector<VType>& ) )
{
	VType eps = CGConstants<VType>::eps;
	Vector<VType> r = b - Multiply(A,x);
	Vector<VType> d = r;
	double delta_new = r.Dot(r);
	double delta_0 = delta_new;
	int i;
	for(i=0; i<iters && delta_new>eps*delta_0 ;i++)
	{
		Vector<VType> q = Multiply(A,d);
		double alpha = delta_new / d.Dot(q);
		x = x + d*alpha;
		if( !(i%50) )	r = b - Multiply(A,x);
		else			r = r - q*alpha;

		double delta_old = delta_new;
		delta_new = r.Dot(r);
		double beta = delta_new / delta_old;
		d = r + d*beta;
	}
	return i;
}

template< class MType , class IndexType , class VType >
int SolveConjugateGradient(const SparseMatrix<MType, IndexType> &A,
	const Vector<VType> &b, const int &iters, Vector<VType> &x,
	VType eps = CGConstants<VType>::eps)
{
	Vector<VType> r = b - A * x;
	Vector<VType> d = r;
	double delta_new = r.Dot(r);
	double delta_0 = delta_new;
	int i;
	for(i = 0; (i < iters) && (delta_new > (eps * delta_0)); i++)
	{
		Vector<VType> q = A * d;
		double alpha = delta_new / d.Dot(q);
		x = x + d*alpha;
		if( !(i%50) )	r = b - A * x;
		else			r = r - q * alpha;

		double delta_old = delta_new;
		delta_new = r.Dot(r);
		double beta = delta_new / delta_old;
		d = r + d*beta;
	}
	return i;
}


////////////////////////////////////////////////////////////////////////////////
/*! Solves a system using conjugate gradient while tracking the intermediate
//  solution at each step
//  @param[in]      A
//  @param[in]      b
//  @param[in]      iters
//  @param[inout]   x
//  @param[out]     xValues
//  @param[in]      Multiply
//  @return         Number of iterations until convergence (or iters)
*///////////////////////////////////////////////////////////////////////////////
template< class MType , class IndexType , class VType >
int SolveConjugateGradient(const SparseMatrix<MType, IndexType> &A
	, const Vector<VType> &b, const int &iters, Vector<VType> &x
	, std::vector<Vector<VType> > &xValues
	, Vector<VType> (*Multiply)(const SparseMatrix<MType, IndexType>&
		, const Vector<VType>& ) )
{
	double eps=1e-16;
	Vector<VType> r = b - Multiply(A,x);
	Vector<VType> d = r;
	double delta_new = r.Dot(r);
	double delta_0 = delta_new;
	int i;
	for(i=0; i<iters && delta_new>eps*delta_0 ;i++)
	{
		Vector<VType> q = Multiply(A,d);
		double alpha = delta_new / d.Dot(q);
		x = x + d*alpha;
		if( !(i%50) )	r = b - Multiply(A,x);
		else			r = r - q*alpha;

		xValues.push_back(x);

		double delta_old = delta_new;
		delta_new = r.Dot(r);
		double beta = delta_new / delta_old;
		d = r + d*beta;
	}
	return i;
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
		if( !_entries ) fprintf( stderr , "[ERROR] BandedMatrix::resize: Failed to allocate BandedMatrix::_entries ( %d x %d )\n" , rows , 2*Radius+1 ) , exit( 0 );
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
		ThreadPool::ParallelFor
		(
			1 , _rows-1 ,
			[&]( unsigned int , size_t i )
			{
				ConstPointer( T ) __entries = _entries + i * 3;
				ConstPointer( T2 ) _in = in + i - 1;
				out[i] = (T2)( _in[0] * __entries[0] + _in[1] * __entries[1] + _in[2] * __entries[2] );
			} ,
			threads
		);
	}
	else
	{
		ThreadPool::ParallelFor
		(
			Radius , _rows-Radius ,
			[&]( unsigned int , size_t i )
			{
				T2 sum(0);
				ConstPointer( T ) __entries = _entries + i * ( 2 * Radius + 1 );
				ConstPointer( T2 ) _in = in + i - Radius;
				for( int j=0 ; j<=2*Radius ; j++ ) sum += (T2)( _in[j] * __entries[j] );
				out[i] = sum;
			} ,
			threads
		);
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
		ThreadPool::ParallelFor
		(
			1 , _rows-1 ,
			[&]( unsigned int , size_t i )
			{
				ConstPointer( T ) __entries = _entries + i * 3;
				ConstPointer( T2 ) _in = in + (i-1)*2;
				out[ i<<1   ] = (T2)( _in[0] * __entries[0] + _in[2] * __entries[1] + _in[4] * __entries[2] );
				out[(i<<1)|1] = (T2)( _in[1] * __entries[0] + _in[3] * __entries[1] + _in[5] * __entries[2] );
			} ,
			threads
		);
	}
	else
	{
		ThreadPool::ParallelFor
		(
			Radius , _rows-Radius ,
			[&]( unsigned int , size_t i )
			{
				T2 sum0(0) , sum1(0);
				ConstPointer( T ) __entries = _entries + i * ( 2 * Radius + 1 );
				ConstPointer( T2 ) _in = in + (i-Radius)*2;
				for( int j=0 ; j<=2*Radius ; j++ ) sum0 += (T2)( _in[j<<1] * __entries[j] ) ,  sum1 += (T2)( _in[(j<<1)|1] * __entries[j] );
				out[ i<<1   ] = sum0;
				out[(i<<1)|1] = sum1;
			} ,
			threads
		);
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

template< class RealIn, class RealOut>
void CompressSparseMatrix( const SparseMatrix< RealIn , int > & M, SparseMatrix< RealOut, int > & _M )
{
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
