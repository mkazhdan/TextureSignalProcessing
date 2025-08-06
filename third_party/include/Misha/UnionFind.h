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
#ifndef UNION_FIND_INCLUDED
#define UNION_FIND_INCLUDED

#include <optional>

namespace MishaK
{
	// A union/find data-structure.
	// The clusters are identified by indices s.t. i == find(i)
	class UnionFind
	{
		struct Node
		{
			size_t sz;
			Node *root , *previous , *tail;
			Node( void ) { root = previous = tail = nullptr , sz=0; }
		};
		Node* _nodes;
		size_t _sz;
	public:
		UnionFind( void ) { _nodes = nullptr; }
		UnionFind( size_t sz ) : _sz(sz)
		{
			_nodes = new Node[ sz ];
			for( size_t i=0 ; i<sz ; i++ ) _nodes[i].root = _nodes[i].tail = _nodes+i , _nodes[i].previous = nullptr , _nodes[i].sz=1;
		};
		~UnionFind( void ){ if( _nodes ) delete[] _nodes; }
		bool isRoot( size_t idx ) const { return idx==find(idx); }
		size_t find( size_t idx ) const { return _nodes[idx].root - _nodes; }
		size_t size( size_t idx ) const { return _nodes[ find(idx) ].sz; }
		size_t merge( size_t idx1 , size_t idx2 )
		{
			size_t _idx1 , _idx2;
			size_t sz1 = size( idx1 ) , sz2 = size( idx2 );
			if( sz1<sz2 ) _idx1 = find( idx1 ) , _idx2 = find( idx2 );
			else          _idx1 = find( idx2 ) , _idx2 = find( idx1 );
			if( _idx1==_idx2 ) return sz1;
			// _idx1 points to the smaller of the two sets
			{
				Node* p = _nodes + _idx2;
				Node* c = _nodes + _idx1;
				p->sz = sz1 + sz2;
				_nodes[_idx2].tail->previous = c;
				_nodes[_idx2].tail = c->tail;
				while( c ) { c->root = p , c = c->previous; }
			}
			return sz1 + sz2;
		}

		std::optional< size_t > previous( size_t idx ) const
		{
			if( _nodes[idx].previous==nullptr ) return {};
			else return _nodes[idx].previous-_nodes;
		}

		size_t components( void ) const
		{
			size_t idx = 0;
			for( size_t i=0 ; i<_sz ; i++ ) if( isRoot(i) ) idx++;
			return idx;
		}

		template< typename Functor /* = std::function< void ( size_t ) */ >
		void process( size_t idx , Functor && F ) const
		{
			static_assert( std::is_convertible_v< Functor , std::function< void ( size_t ) > > , "[ERROR] Functor poorly formed" );
			F( idx );
			while( auto prev = previous(idx) )
			{
				idx = *prev;
				F( idx );
			}

		}
	};
};
#endif // UNION_FIND_INCLUDED