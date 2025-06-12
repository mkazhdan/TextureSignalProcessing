/*
Copyright (c) 2025, Michael Kazhdan
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
#ifndef EXPLICIT_INDEX_INCLUDED
#define EXPLICIT_INDEX_INCLUDED

namespace MishaK
{
	// A wrapper class for unsigned int that allows distinguishing between different types of indices
	template< typename IndexType , typename Index >
	struct ExplicitIndex
	{
		template< typename I=IndexType , typename = std::enable_if_t< std::is_integral_v< I > > >
		explicit ExplicitIndex( I idx=I(-1) ) : _idx( static_cast< IndexType >(idx) ){}

		template< typename _Index >
		explicit ExplicitIndex( const ExplicitIndex< IndexType , _Index > & t ) : _idx( static_cast< IndexType >(t) ) {}

		explicit operator IndexType () const { return _idx; }

		ExplicitIndex & operator += ( IndexType off ){ _idx += off ; return *this; }
		ExplicitIndex & operator -= ( IndexType off ){ _idx -= off ; return *this; }
		Index operator + ( IndexType off ) const { return Index( _idx+off ); }
		Index operator - ( IndexType off ) const { return Index( _idx-off ); }

		IndexType operator - ( const ExplicitIndex &t ) const { return _idx - t._idx; }
		friend Index & operator ++ ( Index & t ){ t._idx++ ; return t; }
		friend Index operator ++ ( Index &t , int ){ return Index( t._idx++ ); }

		bool operator == ( const ExplicitIndex &idx ) const { return _idx==idx._idx; }
		bool operator != ( const ExplicitIndex &idx ) const { return _idx!=idx._idx; }
		bool operator <  ( const ExplicitIndex &idx ) const { return _idx< idx._idx; }
		bool operator <= ( const ExplicitIndex &idx ) const { return _idx<=idx._idx; }
		bool operator >  ( const ExplicitIndex &idx ) const { return _idx> idx._idx; }
		bool operator >= ( const ExplicitIndex &idx ) const { return _idx>=idx._idx; }

		friend std::ostream &operator << ( std::ostream &os , const ExplicitIndex &I ){ return os << I._idx; }

		template< typename Data >
		struct ExplicitIndexVector : protected std::vector< Data >
		{
			ExplicitIndexVector( void ) : std::vector< Data >(){}
			ExplicitIndexVector( size_t sz ) : std::vector< Data >(sz){}
			ExplicitIndexVector( size_t sz , const Data & data ) : std::vector< Data >( sz , data ){}
			size_t size( void ) const { return std::vector< Data >::size(); }
			void resize( size_t sz ){ return std::vector< Data >::resize(sz); }
			void resize( size_t sz , const Data & data ){ return std::vector< Data >::resize( sz , data ); }
			void push_back( const Data & data ){ return std::vector< Data >::push_back(data); }
			Data & operator[]( const Index & t ){ return std::vector< Data >::operator[]( static_cast< IndexType >(t) ); }
			const Data & operator[]( const Index &t ) const { return std::vector< Data >::operator[]( static_cast< IndexType >(t) ); }

			typename std::vector< Data >::iterator begin( void ){ return std::vector< Data >::begin(); }
			typename std::vector< Data >::iterator end( void ){ return std::vector< Data >::end(); }
			typename std::vector< Data >::const_iterator begin( void ) const { return std::vector< Data >::begin(); }
			typename std::vector< Data >::const_iterator end( void ) const { return std::vector< Data >::end(); }

			template< class ... Args >
			void emplace_back( Args && ... args ){ std::vector< Data >::emplace_back( std::forward< Args >(args) ... ); }

			template< class InputIt >
			typename std::vector< Data >::iterator insert( typename std::vector< Data >::const_iterator pos , InputIt first , InputIt last ){ return std::vector< Data >::insert( pos , first , last ); }

			explicit operator const std::vector< Data > & () const { return *this; }
			explicit operator       std::vector< Data > & ()       { return *this; }
		};
	protected:
		IndexType _idx;
	};
	template< typename T , typename Data >
	Data * operator + ( Data * data , const T & t ){ return data + static_cast< unsigned int >(t); }

	template< typename T , typename Data >
	const Data * operator + ( const Data * data , const T & t ){ return data + static_cast< unsigned int >(t); }
}
#endif // EXPLICIT_INDEX_INCLUDED