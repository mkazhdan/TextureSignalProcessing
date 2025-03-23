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

template< class T , class const_iterator > size_t SparseMatrixInterface< T , const_iterator >::Entries( void ) const
{
	size_t entries = 0;
	for( size_t i=0 ; i<Rows() ; i++ ) entries += RowSize( i );
	return entries;
}
template< class T , class const_iterator > double SparseMatrixInterface< T , const_iterator >::SquareNorm( void ) const
{
	double n=0;
	for( size_t i=0 ; i<Rows() ; i++ )
	{
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) n += iter->Value * iter->Value;
	}
	return n;

}
template< class T , class const_iterator > double SparseMatrixInterface< T , const_iterator >::SquareASymmetricNorm( void ) const
{
	double n=0;
	for( size_t i=0 ; i<Rows() ; i++ )
	{
		const_iterator e = end( i );
		for( const_iterator iter1 = begin( i ) ; iter1!=e ; iter1++ )
		{
			int j = iter1->N;
			const_iterator e = end( j );
			double value = 0;
			for( const_iterator iter2 = begin( j ) ; iter2!=e ; iter2++ )
			{
				int k = iter2->N;
				if( k==i ) value += iter2->Value;
			}
			n += (iter1->Value-value) * (iter1->Value-value);
		}
	}
	return n;
}
template< class T , class const_iterator > double SparseMatrixInterface< T , const_iterator >::SquareASymmetricNorm( int& idx1 , int& idx2 ) const
{
	double n=0;
	double max=0;
	for( size_t i=0 ; i<Rows() ; i++ )
	{
		const_iterator e = end( i );
		for( const_iterator iter = begin( i ) ; iter!=e ; iter++ )
		{
			int j = iter->N;
			const_iterator e = end( j );
			double value = 0;
			for( const_iterator iter2 = begin( j ) ; iter2!=e ; iter2++ )
			{
				int k = iter2->N;
				if( k==i ) value += iter2->Value;
			}
			double temp = (iter->Value-value) * (iter->Value-value);
			n += temp;
			if( temp>=max ) idx1 = i , idx2 = j , max=temp;
		}
	}
	return n;
}
template< class T , class const_iterator >
template< class T2 , class T2Real >
void SparseMatrixInterface< T , const_iterator >::Multiply( ConstPointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag ) const
{
	ConstPointer( T2 ) in = In;
	ThreadPool::ParallelFor
		(
			0 , Rows() ,
			[&]( unsigned int , size_t i )
			{
				T2 temp;
				memset( &temp , 0 , sizeof(T2) );
				ConstPointer( T2 ) _in = in;
				const_iterator e = end( i );
#if 1
				for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) temp += _in[ iter->N ] * static_cast< T2Real >( iter->Value );
#else
				for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) temp += (T2)( _in[ iter->N ] * iter->Value );
#endif
				if( multiplyFlag & MULTIPLY_NEGATE ) temp = -temp;
				if( multiplyFlag & MULTIPLY_ADD ) Out[i] += temp;
				else                              Out[i]  = temp;
			}
		);
}
template< class T , class const_iterator >
template< class T2 , class T2Real >
void SparseMatrixInterface< T , const_iterator >::MultiplyScaled( T scale , ConstPointer( T2 ) In , Pointer( T2 ) Out , int multiplyFlag ) const
{
	ConstPointer( T2 ) in = In;
	ThreadPool::ParallelFor
		(
			0 , Rows() ,
			[&]( unsigned int , size_t i )
			{
				T2 temp;
				memset( &temp , 0 , sizeof(T2) );
				ConstPointer( T2 ) _in = in;
				const_iterator e = end( i );
				for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) temp += _in[ iter->N ] * static_cast< T2Real >( iter->Value );
				temp *= scale;
				if( multiplyFlag & MULTIPLY_NEGATE ) temp = -temp;
				if( multiplyFlag & MULTIPLY_ADD ) Out[i] += temp;
				else                              Out[i]  = temp;
			}
		);
}

template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::SetDiagonal( Pointer( T2 ) diagonal ) const
{
	ThreadPool::ParallelFor
		(
			0 , Rows() ,
			[&]( unsigned int , size_t i )
			{
				diagonal[i] = (T2)0;
				const_iterator e = end( i );
				for( const_iterator iter = begin( i ) ; iter!=e ; iter++ ) if( iter->N==i ) diagonal[i] += iter->Value;
			}
		);
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::JacobiIteration( ConstPointer( T2 ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , Pointer( T2 ) Mx , T2 sor ) const
{
	Multiply( x , Mx );
	ThreadPool::ParallelFor( 0 , Rows() , [&]( unsigned int , size_t i ){ x[i] += ( b[i] - Mx[i] ) * sor / diagonal[i]; } );
}
#if 1
template< class T , class const_iterator >
template< class T2 , bool StripDiagonal >
void SparseMatrixInterface< T , const_iterator >::GSIteration( ConstPointer( T2 ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward ) const
{
	if( StripDiagonal )
	{
#define ITERATE( j )                                                                                    \
		{                                                                                               \
			T2 _b = b[j];                                                                               \
			const_iterator e = end( j );                                                                \
			for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;  \
			x[j] = _b / diagonal[j];                                                                    \
		}

		if( forward ) for( int j=0 ; j<int( Rows() ) ; j++ ){ ITERATE( j ); }
		else          for( int j=int( Rows() )-1 ; j>=0 ; j-- ){ ITERATE( j ); }
#undef ITERATE
	}
	else
	{
#define ITERATE( j )                                                                                    \
		{                                                                                               \
			T2 _b = b[j];                                                                               \
			const_iterator e = end( j );                                                                \
			for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;  \
			x[j] += _b / diagonal[j];                                                                   \
		}

		if( forward ) for( int j=0 ; j<int( Rows() ) ; j++ ){ ITERATE( j ); }
		else          for( int j=int( Rows() )-1 ; j>=0 ; j-- ){ ITERATE( j ); }
#undef ITERATE
	}
}
template< class T , class const_iterator >
template< class T2 , bool StripDiagonal >
void SparseMatrixInterface< T , const_iterator >::GSIteration( std::vector< std::vector< int > >& multiColorIndices , ConstPointer( T2 ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward ) const
{
#ifdef _WIN32
#define SetOMPParallel __pragma( omp parallel for )
#else // !_WIN32
#define SetOMPParallel _Pragma( "omp parallel for" )
#endif // _WIN32

	if( StripDiagonal )
	{
#define ITERATE( indices )                                                                                   \
		{                                                                                                    \
SetOMPParallel                                                                                               \
			for( int k=0 ; k<int( indices.size() ) ; k++ )                                                   \
			{                                                                                                \
				int jj = indices[k];                                                                         \
				T2 _b = b[jj];                                                                               \
				const_iterator e = end( jj );                                                                \
				for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;  \
				x[jj] = _b / diagonal[jj];                                                                   \
			}                                                                                                \
		}
		if( forward ) for( int j=0 ; j<multiColorIndices.size()  ; j++ ){ ITERATE( multiColorIndices[j] ); }
		else for( int j=int( multiColorIndices.size() )-1 ; j>=0 ; j-- ){ ITERATE( multiColorIndices[j] ); }
#undef ITERATE
	}
	else
	{
#define ITERATE( indices )                                                                                   \
		{                                                                                                    \
SetOMPParallel                                                                                               \
			for( int k=0 ; k<int( indices.size() ) ; k++ )                                                   \
			{                                                                                                \
				int jj = indices[k];                                                                         \
				T2 _b = b[jj];                                                                               \
				const_iterator e = end( jj );                                                                \
				for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;  \
				x[jj] += _b / diagonal[jj];                                                                  \
			}                                                                                                \
		}
		if( forward ) for( int j=0 ; j<multiColorIndices.size()  ; j++ ){ ITERATE( multiColorIndices[j] ); }
		else for( int j=int( multiColorIndices.size() )-1 ; j>=0 ; j-- ){ ITERATE( multiColorIndices[j] ); }
#undef ITERATE
	}
#undef SetOMPParallel
}
#else
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::GSIteration( ConstPointer( T2 ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward ) const
{
#define ITERATE( j )                                                                                \
	{                                                                                               \
		T2 _b = b[j];                                                                               \
		const_iterator e = end( j );                                                                \
		for( const_iterator iter = begin( j ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;  \
		x[j] += _b / diagonal[j];                                                                   \
	}

	if( forward ) for( int j=0 ; j<int( Rows() ) ; j++ ){ ITERATE( j ); }
	else          for( int j=int( Rows() )-1 ; j>=0 ; j-- ){ ITERATE( j ); }
#undef ITERATE
}
template< class T , class const_iterator >
template< class T2 >
void SparseMatrixInterface< T , const_iterator >::GSIteration( std::vector< std::vector< int > >& multiColorIndices , ConstPointer( T2 ) diagonal , ConstPointer( T2 ) b , Pointer( T2 ) x , bool forward ) const
{
#ifdef _WIN32
#define SetOMPParallel __pragma( omp parallel for )
#else // !_WIN32
#define SetOMPParallel _Pragma( "omp parallel for" )
#endif // _WIN32
#define ITERATE( indices )                                                                               \
	{                                                                                                    \
SetOMPParallel                                                                                           \
		for( int k=0 ; k<int( indices.size() ) ; k++ )                                                   \
		{                                                                                                \
			int jj = indices[k];                                                                         \
			T2 _b = b[jj];                                                                               \
			const_iterator e = end( jj );                                                                \
			for( const_iterator iter = begin( jj ) ; iter!=e ; iter++ ) _b -= x[iter->N] * iter->Value;  \
			x[jj] += _b / diagonal[jj];                                                                  \
		}                                                                                                \
	}
	if( forward ) for( int j=0 ; j<multiColorIndices.size()  ; j++ ){ ITERATE( multiColorIndices[j] ); }
	else for( int j=int( multiColorIndices.size() )-1 ; j>=0 ; j-- ){ ITERATE( multiColorIndices[j] ); }
#undef ITERATE
#undef SetOMPParallel
}
#endif
template< class SPDOperator , class T > int SolveCG( const SPDOperator& M , ConstPointer( T ) b , int iters , Pointer( T ) x , T eps , bool solveNormal )
{
	eps *= eps;
	int dim = M.Rows();
	Pointer( T ) r = AllocPointer< T >( dim );
	Pointer( T ) d = AllocPointer< T >( dim );
	Pointer( T ) q = AllocPointer< T >( dim );
	Pointer( T ) temp = NullPointer< T >( );
	memset( x , 0 , sizeof(T)* dim );
	if( solveNormal ) temp = AllocPointer< T >( dim );

	double delta_new = 0 , delta_0;
	if( solveNormal )
	{
		M.Multiply( ( ConstPointer( T ) )x , temp ) , M.Multiply( ( ConstPointer( T ) )temp , r ) , M.Multiply( ( ConstPointer( T ) )b , temp );
		std::vector< double > _delta_news( ThreadPool::NumThreads() , 0 );
		ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int t , size_t i ){ d[i] = r[i] = temp[i] - r[i] ; _delta_news[t] += r[i] * r[i]; } );
		for( unsigned int t=0 ; t<_delta_news.size() ; t++ ) delta_new += _delta_news[t];
	}
	else
	{
		M.Multiply( ( ConstPointer( T ) )x , r );
		std::vector< double > _delta_news( ThreadPool::NumThreads() , 0 );
		ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int t , size_t i ){ d[i] = r[i] = b[i] - r[i] ; _delta_news[t] += r[i] * r[i]; } );
		for( unsigned int t=0 ; t<_delta_news.size() ; t++ ) delta_new += _delta_news[t];
	}
	delta_0 = delta_new;
	if( delta_new<eps )
	{
		FreePointer( r );
		FreePointer( d );
		FreePointer( q );
		FreePointer( temp );
		return 0;
	}
	int ii;
	for( ii=0 ; ii<iters && delta_new>eps*delta_0 ; ii++ )
	{
		if( solveNormal ) M.Multiply( ( ConstPointer( T ) )d , temp ) , M.Multiply( ( ConstPointer( T ) )temp , q );
		else              M.Multiply( ( ConstPointer( T ) )d , q );
        double dDotQ = 0;
		std::vector< double > _dDotQs( ThreadPool::NumThreads() , 0 );
		ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int t , size_t i ){ _dDotQs[t] += _dDotQs[t];} );
		for( unsigned int t=0 ; t<_dDotQs.size() ; t++ ) dDotQ += _dDotQs[t];
		T alpha = T( delta_new / dDotQ );
		double delta_old = delta_new;
		delta_new = 0;
		if( (ii%50)==(50-1) )
		{
			ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int , size_t i ){ x[i] += d[i] * alpha; } );
			if( solveNormal ) M.Multiply( ( ConstPointer( T ) )x , temp ) , M.Multiply( ( ConstPointer( T ) )temp , r );
			else              M.Multiply( ( ConstPointer( T ) )x , r );
			std::vector< double > _delta_news( ThreadPool::NumThreads() , 0 );
			ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int t , size_t i ){ r[i] = b[i] - r[i] ; _delta_news[t] += r[i] * r[i] ; x[i] += d[i] * alpha; } );
			for( unsigned int t=0 ; t<_delta_news.size() ; t++ ) delta_new += _delta_news[t];
		}
		else
		{
			std::vector< double > _delta_news( ThreadPool::NumThreads() , 0 );
			ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int t , size_t i ){ r[i] -= q[i] * alpha ; _delta_news[t] += r[i] * r[i] ;  x[i] += d[i] * alpha; } );
			for( unsigned int t=0 ; t<_delta_news.size() ; t++ ) delta_new += _delta_news[t];
		}

		T beta = T( delta_new / delta_old );
		ThreadPool::ParallelFor( 0 , dim , [&]( unsigned int , size_t i ){ d[i] = r[i] + d[i] * beta; } );
	}
	FreePointer( r );
	FreePointer( d );
	FreePointer( q );
	FreePointer( temp );
	return ii;
}
