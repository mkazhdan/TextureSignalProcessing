/*
Copyright (c) 2011, Michael Kazhdan and Ming Chuang
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

#include <stdio.h>
#include <emmintrin.h>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#endif // _WIN32
#include <stddef.h>

inline bool isfinitef( float fp ){ float f=fp; return ((*(unsigned *)&f)&0x7f800000)!=0x7f800000; }


template< class C >        bool IsValid( const C& c );
#if _DEBUG
template< >         inline bool IsValid< float >( const float& f ) { return isfinitef( f ) &&  ( f==0.f || abs(f)>1e-31f ); }
#else // !_DEBUG
template< >         inline bool IsValid< float >( const float& f ) { return isfinitef( f ); }
#endif // _DEBUG
template< >         inline bool IsValid< __m128 >( const __m128& m )
{
	const __m128* addr = &m;
	if( size_t(addr) & 15 ) return false;
	else                    return true;
}
template< class C > inline bool IsValid( const C& c ){ return true; }


template< class C >
class Array
{
public:
	typedef long long difference_type;
	typedef C value_type;
	typedef C& reference;
	typedef C* pointer;
	typedef std::random_access_iterator_tag iterator_category;
protected:
	void _assertBounds( difference_type idx , const char* message=NULL ) const
	{
		if( idx<min || idx>=max )
		{
			if( message ) Miscellany::ErrorOut( "Array index out-of-bounds: %lld <= %lld < %lld [%s]" , min , idx , max , message );
			else          Miscellany::ErrorOut( "Array index out-of-bounds: %lld <= %lld < %lld" , min , idx , max );
		}
	}
	C *data , *_data;
	difference_type min , max;

public:
	difference_type minimum( void ) const { return min; }
	difference_type maximum( void ) const { return max; }

	operator C*() { return data; }
	static inline Array New( size_t size )
	{
		Array a;
		a._data = a.data = new C[size];
		a.min = 0;
#pragma message( "[WARNING] Casting unsigned to signed" )
		a.max = ( difference_type ) size;
		return a;
	}
	static inline Array Alloc( size_t size , bool clear )
	{
		Array a;
		a._data = a.data = ( C* ) malloc( size * sizeof( C ) );
		if( clear ) memset( a.data ,  0 , size * sizeof( C ) );
		//		else        memset( a.data , -1 , size * sizeof( C ) );
		a.min = 0;
#pragma message( "[WARNING] Casting unsigned to signed" )
		a.max = ( difference_type ) size;
		return a;
	}
	static inline Array AlignedAlloc( size_t size , size_t alignment , bool clear )
	{
		Array a;
		a.data = ( C* ) aligned_malloc( sizeof(C) * size , alignment );
		a._data = ( C* )( ( ( void** )a.data )[-1] );
		if( clear ) memset( a.data ,  0 , size * sizeof( C ) );
		//		else        memset( a.data , -1 , size * sizeof( C ) );
		a.min = 0;
#pragma message( "[WARNING] Casting unsigned to signed" )
		a.max = ( difference_type ) size;
		return a;
	}
	static inline Array ReAlloc( Array& a , size_t size , bool clear )
	{
		Array _a;
		_a._data = _a.data = ( C* ) realloc( a.data , size * sizeof( C ) );
		if( clear ) memset( _a.data ,  0 , size * sizeof( C ) );
		a._data = NULL;
		_a.min = 0;
#pragma message( "[WARNING] Casting unsigned to signed" )
		_a.max = ( difference_type ) size;
		return _a;
	}

	Array( void )
	{
		data = _data = NULL;
		min = max = 0;
	}
	template< class D >
	Array( Array< D >& a )
	{
		_data = NULL;
		if( !a )
		{
			data =  NULL;
			min = max = 0;
		}
		else
		{
			// [WARNING] Changing szC and szD to size_t causes some really strange behavior.
			difference_type szC = sizeof( C );
			difference_type szD = sizeof( D );
#if 1
			// Since we are not actually accessing the values, we should not be bounds testing
			data = (C*)a.ptr();
#else
			data = (C*)&a[0];
#endif
			min = ( a.minimum() * szD ) / szC;
			max = ( a.maximum() * szD ) / szC;
			if( min*szC!=a.minimum()*szD || max*szC!=a.maximum()*szD ) Miscellany::ErrorOut( "Could not convert array [ %lld , %lld ] * %lld => [ %lld , %lld ] * %lld" , a.minimum() , a.maximum() , szD , min , max , szC );
		}
	}
	static Array FromPointer( C* data , difference_type max )
	{
		Array a;
		a._data = NULL;
		a.data = data;
		a.min = 0;
		a.max = max;
		return a;
	}
	static Array FromPointer( C* data , difference_type min , difference_type max )
	{
		Array a;
		a._data = NULL;
		a.data = data;
		a.min = min;
		a.max = max;
		return a;
	}
	inline bool operator == ( const Array< C >& a ) const { return data==a.data; }
	inline bool operator != ( const Array< C >& a ) const { return data!=a.data; }
	inline bool operator == ( const C* c ) const { return data==c; }
	inline bool operator != ( const C* c ) const { return data!=c; }
	inline bool operator <  ( const Array< C >& a ) const { return data< a.data; }
	inline bool operator >  ( const Array< C >& a ) const { return data> a.data; }
	inline bool operator <= ( const Array< C >& a ) const { return data<=a.data; }
	inline bool operator >= ( const Array< C >& a ) const { return data>=a.data; }

	inline       C& operator *  ( void )       { _assertBounds( 0 , "*" ) ; return data[0]; }
	inline const C& operator *  ( void ) const { _assertBounds( 0 , "const *" ) ; return data[0]; }
	inline       C* operator -> ( void )       { _assertBounds( 0 , "->" ) ; return data; }
	inline const C* operator -> ( void ) const { _assertBounds( 0 , "const ->" ) ; return data; }
	template< class Offset > inline       C& operator [] ( Offset idx )       { _assertBounds( idx , "[]" ) ; return data[idx]; }
	template< class Offset > inline const C& operator [] ( Offset idx ) const { _assertBounds( idx , "const []" ) ; return data[idx]; }

	template< class Offset >
	inline Array operator + ( Offset idx ) const
	{
		Array a;
		a._data = _data;
		a.data = data+idx;
		a.min = min-idx;
		a.max = max-idx;
		return a;
	}
	template< class Offset >
	inline Array& operator += ( Offset idx  )
	{
		min -= idx;
		max -= idx;
		data += idx;
		return (*this);
	}
	inline Array& operator ++ ( void ) { return (*this) += 1; }
	inline Array  operator ++ ( int ) { Array tmp = (*this) ; (*this) +=1 ; return tmp; }
	template< class Offset > Array  operator -  ( Offset idx ) const { return (*this) +  (-idx); }
	template< class Offset > Array& operator -= ( Offset idx )       { return (*this) += (-idx); }
	Array& operator -- ( void ) { return (*this) -= 1; }
	Array  operator -- ( int ) { Array tmp = (*this) ; (*this) -= 1 ; return tmp; }
	difference_type operator - ( const Array& a ) const { return ( difference_type )( data - a.data ); }

	void Free( void )
	{
		if( _data )
		{
			free( _data );
		}
		(*this) = Array( );
	}
	void Delete( void )
	{
		if( _data )
		{
			delete[] _data;
		}
		(*this) = Array( );
	}
	C* ptr( void ){ return data; }
	const C* ptr( void ) const { return data; }
	bool operator !( void ) const { return data==NULL; }
	operator bool( ) const { return data!=NULL; }
};
template< class C , class Offset > Array< C > operator + ( Offset idx , const Array< C >& a ){ return (a+idx); }

template< class C >
class ConstArray
{
public:
	typedef long long difference_type;
protected:
	void _assertBounds( difference_type idx ) const
	{
		if( idx<min || idx>=max ) Miscellany::ErrorOut( "ConstArray index out-of-bounds: %lld <= %lld < %lld" , min , idx , max );
	}
protected:
	const C *data;
	difference_type min , max;
public:
	difference_type minimum( void ) const { return min; }
	difference_type maximum( void ) const { return max; }

	operator const C*() { return data; }
	inline ConstArray( void )
	{
		data = NULL;
		min = max = 0;
	}
	inline ConstArray( const Array< C >& a )
	{
		// [WARNING] Changing szC and szD to size_t causes some really strange behavior.
		data = ( const C* )a.ptr( );
		min = a.minimum();
		max = a.maximum();
	}
	template< class D >
	inline ConstArray( const Array< D >& a )
	{
		// [WARNING] Changing szC and szD to size_t causes some really strange behavior.
		difference_type szC = ( difference_type ) sizeof( C );
		difference_type szD = ( difference_type ) sizeof( D );
		data = ( const C* )a.ptr();
		min = ( a.minimum() * szD ) / szC;
		max = ( a.maximum() * szD ) / szC;
		if( min*szC!=a.minimum()*szD || max*szC!=a.maximum()*szD )
			Miscellany::ErrorOut( "Could not convert const array [ %lld , %lld ] * %lld => [ %lld , %lld ] * %lld %lld %lld %lld" , a.minimum() , a.maximum() , szD , min , max , szC , a.minimum() , a.minimum()*szD , (a.minimum()*szD)/szC );
	}
	template< class D >
	inline ConstArray( const ConstArray< D >& a )
	{
		// [WARNING] Chaning szC and szD to size_t causes some really strange behavior.
		difference_type szC = sizeof( C );
		difference_type szD = sizeof( D );
		data = ( const C*)a.ptr( );
		min = ( a.minimum() * szD ) / szC;
		max = ( a.maximum() * szD ) / szC;
		if( min*szC!=a.minimum()*szD || max*szC!=a.maximum()*szD ) Miscellany::ErrorOut( "Could not convert array [ %lld , %lld ] * %lld => [ %lld , %lld ] * %lld" , a.minimum() , a.maximum() , szD , min , max , szC );
	}
	static ConstArray FromPointer( const C* data , difference_type max )
	{
		ConstArray a;
		a.data = data;
		a.min = 0;
		a.max = max;
		return a;
	}
	static ConstArray FromPointer( const C* data , difference_type min , difference_type max )
	{
		ConstArray a;
		a.data = data;
		a.min = min;
		a.max = max;
		return a;
	}

	inline bool operator == ( const ConstArray< C >& a ) const { return data==a.data; }
	inline bool operator != ( const ConstArray< C >& a ) const { return data!=a.data; }
	inline bool operator == ( const C* c ) const { return data==c; }
	inline bool operator != ( const C* c ) const { return data!=c; }
	inline const C* operator -> ( void )
	{
		_assertBounds( 0 );
		return data;
	}
	template< class Offset > inline const C& operator[]( Offset idx ) const
	{
		_assertBounds( idx );
		return data[idx];
	}
	template< class Offset > inline ConstArray operator + ( Offset idx ) const
	{
		ConstArray a;
		a.data = data+idx;
		a.min = min-idx;
		a.max = max-idx;
		return a;
	}
	template< class Offset > inline ConstArray& operator += ( Offset idx  )
	{
		min -= idx;
		max -= idx;
		data += idx;
		return (*this);
	}
	inline ConstArray& operator ++ ( void ) { return (*this) += 1; }
	template< class Offset > ConstArray  operator -  ( Offset idx ) const { return (*this) +  (-idx); }
	template< class Offset > ConstArray& operator -= ( Offset idx )       { return (*this) += (-idx); }
	ConstArray& operator -- ( void ) { return (*this) -= 1; }
	difference_type operator - ( const ConstArray& a ) const { return ( difference_type )( data - a.data ); }
	difference_type operator - ( const Array< C >& a ) const { return ( difference_type )( data - a.ptr() ); }

	const C* ptr( void ) const { return data; }
	bool operator !( void ) { return data==NULL; }
	operator bool( ) { return data!=NULL; }
};
template< class C , class Offset > ConstArray< C > operator + ( Offset idx , const ConstArray< C >& a ){ return (a+idx); }


template< class C >
Array< C > memcpy( Array< C > destination , const void* source , size_t size )
{
	if( size>destination.maximum()*sizeof(C) ) Miscellany::ErrorOut( "Size of copy exceeds destination maximum: %lld > %lld" , ( long long )( size ) , ( long long )( destination.maximum()*sizeof( C ) ) );
	if( size ) memcpy( &destination[0] , source , size );
	return destination;
}
template< class C , class D >
Array< C > memcpy( Array< C > destination , Array< D > source , size_t size )
{
	if( size>destination.maximum()*sizeof( C ) ) Miscellany::ErrorOut( "Size of copy exceeds destination maximum: %lld > %lld" , ( long long )( size ) , ( long long )( destination.maximum()*sizeof( C ) ) );
	if( size>source.maximum()*sizeof( D ) ) Miscellany::ErrorOut( "Size of copy exceeds source maximum: %lld > %lld" , ( long long )( size ) , ( long long )( source.maximum()*sizeof( D ) ) );
	if( size ) memcpy( &destination[0] , &source[0] , size );
	return destination;
}
template< class C , class D >
Array< C > memcpy( Array< C > destination , ConstArray< D > source , size_t size )
{
	if( size>destination.maximum()*sizeof( C ) ) Miscellany::ErrorOut( "Size of copy exceeds destination maximum: %lld > %lld" , ( long long )( size ) , ( long  long )( destination.maximum()*sizeof( C ) ) );
	if( size>source.maximum()*sizeof( D ) ) Miscellany::ErrorOut( "Size of copy exceeds source maximum: %lld > %lld" , ( long long )( size ) , ( long long )( source.maximum()*sizeof( D ) ) );
	if( size ) memcpy( &destination[0] , &source[0] , size );
	return destination;
}
template< class D >
void* memcpy( void* destination , Array< D > source , size_t size )
{
	if( size>source.maximum()*sizeof( D ) ) Miscellany::ErrorOut( "Size of copy exceeds source maximum: %lld > %lld" , ( long long )( size ) , ( long long )( source.maximum()*sizeof( D ) ) );
	if( size ) memcpy( destination , &source[0] , size );
	return destination;
}
template< class D >
void* memcpy( void* destination , ConstArray< D > source , size_t size )
{
	if( size>source.maximum()*sizeof( D ) ) Miscellany::ErrorOut( "Size of copy exceeds source maximum: %lld > %lld" , ( long long )( size ) , ( long long )( source.maximum()*sizeof( D ) ) );
	if( size ) memcpy( destination , &source[0] , size );
	return destination;
}
template< class C >
Array< C > memset( Array< C > destination , int value , size_t size )
{
	if( size>destination.maximum()*sizeof( C ) ) Miscellany::ErrorOut( "Size of set exceeds destination maximum: %lld > %lld" , ( long long )( size ) , ( long long )( destination.maximum()*sizeof( C ) ) );
	if( size ) memset( &destination[0] , value , size );
	return destination;
}

template< class C >
size_t fread( Array< C > destination , size_t eSize , size_t count , FILE* fp )
{
	if( count*eSize>destination.maximum()*sizeof( C ) ) Miscellany::ErrorOut( "Size of read exceeds source maximum: %lld > %lld" , ( long long )( count*eSize ) , ( long long )( destination.maximum()*sizeof( C ) ) );
	return fread( &destination[0] , eSize , count , fp );
}
template< class C >
size_t fwrite( Array< C > source , size_t eSize , size_t count , FILE* fp )
{
	if( count*eSize>source.maximum()*sizeof( C ) ) Miscellany::ErrorOut( "Size of write exceeds source maximum: %lld > %lld" , ( long long )( count*eSize ) , ( long long )( source.maximum()*sizeof( C ) ) );
	return fwrite( &source[0] , eSize , count , fp );
}
template< class C >
size_t fwrite( ConstArray< C > source , size_t eSize , size_t count , FILE* fp )
{
	if( count*eSize>source.maximum()*sizeof( C ) ) Miscellany::ErrorOut( "Size of write exceeds source maximum: %lld > %lld" , ( long long )( count*eSize ) , ( long long )( source.maximum()*sizeof( C ) ) );
	return fwrite( &source[0] , eSize , count , fp );
}
template< class C >
void qsort( Array< C > base , size_t numElements , size_t elementSize , int (*compareFunction)( const void* , const void* ) )
{
	if( sizeof(C)!=elementSize ) Miscellany::ErrorOut( "Element sizes differ: %lld != %lld" , ( long long )( sizeof(C) ) , ( long long )( elementSize ) );
	if( base.minimum()>0 || base.maximum()<numElements ) Miscellany::ErrorOut( "Array access out of bounds: %lld <= 0 <= %lld <= %lld" , base.minimum() , base.maximum() , ( long long )( numElements ) );
	qsort( base.ptr() , numElements , elementSize , compareFunction );
}