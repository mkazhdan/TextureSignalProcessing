#ifndef MULTI_INDEX_INCLUDED
#define MULTI_INDEX_INCLUDED

namespace MishaK
{

	template< unsigned int N , unsigned int K >
	struct ChooseTable
	{
		static constexpr unsigned int Choose( void )
		{
			static_assert( N>=K , "[ERROR] K cannot exceed N" );
			if constexpr( K==0 ) return 1;
			else return ( ChooseTable< N-1 , K-1 >::Choose() * N ) / K;
		}

		static unsigned int Choose( unsigned int n , unsigned int k )
		{
			if( k==0 ) return 1;
			else return ( Choose( n-1 , k-1 ) * n ) / k;
		}

		ChooseTable( void ){ for( unsigned int n=0 ; n<N ; n++ ) for( unsigned int k=0 ; k<K ; k++ ) _values[n][k] = Choose( n , k ); }
		const unsigned int operator()( unsigned int n , unsigned int k ) const { return _values[n][k]; }
	protected:
		unsigned int _values[N][K];
	};

	template< unsigned int K , unsigned int MaxN=(unsigned int)-1 >
	struct MultiIndex
	{
		template< unsigned int _K , unsigned int _MaxN > friend struct MultiIndex;

		MultiIndex( void ){ for( unsigned int i=0 ; i<K ; i++ ) _indices[i] = i; }
		MultiIndex( const unsigned int indices[] ){ _init( indices ); }
		MultiIndex(       unsigned int indices[] ){ _init( indices ); }
		template< typename ... UInts > MultiIndex( UInts ... indices )
		{
			static_assert( sizeof ... ( UInts )==K , "[ERROR] Wrong number of indices" );
			const unsigned int _indices[] = { (unsigned int)indices... };
			_init( _indices );
		}
		bool operator < ( const MultiIndex &idx ) const
		{
			for( unsigned int i=0 ; i<K ; i++ )
			{
				if( _indices[i]<idx._indices[i] ) return true;
				else if( _indices[i]>idx._indices[i] ) return false;
			}
			return false;
		}
		bool operator == ( const MultiIndex &idx ) const { return !( (*this)<idx ) && !( idx<(*this) ); }
		const unsigned int &operator[] ( unsigned int idx ) const { return _indices[idx]; }

#pragma message( "[WARNING] Incrementing does not match the inequality ordering" )
		MultiIndex &operator ++( void )
		{
			for( unsigned int i=0 ; i<K ; i++ )
			{
				if( i==K-1 || _indices[i]+1<_indices[i+1] )
				{
					_indices[i]++;
					for( unsigned int j=0 ; j<i ; j++ ) _indices[j] = j;
					break;
				}
			}
			return *this;
		}

		MultiIndex operator ++( int )
		{
			MultiIndex old = *this;
			operator ++ ();
			return old;
		}

		static MultiIndex Begin( unsigned int ){ return MultiIndex(); }
		static MultiIndex End( unsigned int n )
		{
			if( n<K ) return MultiIndex();
			MultiIndex mi;
			for( unsigned int k=0 ; k<K-1 ; k++ ) mi._indices[k] = k;
			mi._indices[K-1] = n;
			return mi;
		}

		static unsigned int Choose( unsigned int n , unsigned int k )
		{
			if( k==0 ) return 1;
			else return ( Choose( n-1 , k-1 ) * n ) / k;
		}

		static unsigned int Size( unsigned int N ){ return Choose( N , K ); }

		unsigned int operator()( void ) const
		{
			// If the k-th index have value x_i, then there are choose( x_i-1 , i ) possible preceeding values
			unsigned int idx = 0;
			if constexpr( MaxN==-1 ) for( int k=(int)K-1 ; k>=0 ; k-- ) idx += Choose( _indices[k] , k+1 );
			else                     for( int k=(int)K-1 ; k>=0 ; k-- ) idx += _ChooseTable( _indices[k] , k+1 );
			return idx;
		}

		template< unsigned int _K >
		bool contains( const MultiIndex< _K > &mi ) const
		{
			if constexpr( _K>K ) return false;
			unsigned int j = 0;
			for( unsigned int i=0 ; i<K && j<_K ; i++ )
			{
				if( _indices[i]>mi._indices[j] ) return false;
				else if( _indices[i]==mi._indices[j] ) j++;
			}
			return j==_K;
		}

		struct Hash
		{
#if 1
			size_t operator() ( const MultiIndex &mi ) const { return mi[0]; }
#else
			size_t operator() ( const MultiIndex &mi ) const
			{
				size_t val = 1;
				for( unsigned int s=0 ; s<K ; s++ ) val *= mi[s];
				return val;
			}
#endif
		};

		template< typename Data >
		struct Map
		{
			void resize( unsigned int N )
			{
				clear();
				_data.resize( MultiIndex::Size(N) );
				for( unsigned int i=0 ; i<_data.size() ; i++ ) _data[i].first = false;
				_indices.reserve( MultiIndex::Size(N) );
			}

			Data &operator[]( MultiIndex mi )
			{
				unsigned int idx = mi();
				if( !_data[idx].first )
				{
					_indices.push_back( idx );
					_data[idx].first = true;
					_data[idx].second = {};
				}
				return _data[idx].second;
			}

			Data &operator[]( unsigned int idx ){ return _data[ _indices[idx] ].second; }
			size_t size( void ) const { return _indices.size(); }
			void clear( void )
			{
				for( unsigned int i=0 ; i<_indices.size() ; i++ ) _data[ _indices[i] ].first = false;
				_indices.resize(0);
			}

			struct iterator
			{
				iterator( unsigned int idx , Map &map ) : _idx(idx) , _map(map) {}
				iterator &operator ++( void ){ _idx++ ; return *this; }
				iterator operator ++( int ){ return iterator( _idx++ , _map ); }
				bool operator != ( iterator iter ) const { return _idx!=iter._idx; }
				Data &operator *( void ){ return _map[_idx]; }
			protected:
				unsigned int _idx;
				Map &_map;
			};
			iterator begin( void ){ return iterator( 0 , *this ); }
			iterator end( void ){ return iterator( (unsigned int)_data.size() , *this ); }
		protected:
			std::vector< std::pair< bool , Data > > _data;
			std::vector< unsigned int > _indices;
		};

	protected:
		static std::conditional_t< MaxN==-1 , char , ChooseTable< MaxN , K+1 > > _ChooseTable;

		void _init( const unsigned int indices[] )
		{
			memcpy( _indices , indices, sizeof(unsigned int) * K );
			std::sort( _indices , _indices + K , []( unsigned int v1 , unsigned int v2 ){ return v1<v2; } );
		}
		unsigned int _indices[ K ? K : 1 ];
	};

	template< unsigned int K , unsigned int MaxN >
	std::conditional_t< MaxN==-1 , char , ChooseTable< MaxN , K+1 > > MultiIndex< K , MaxN >::_ChooseTable;

	template< unsigned int K >
	std::ostream &operator << ( std::ostream &os , const MultiIndex< K > &mi )
	{
		os << "{ " << mi[0];
		for( unsigned int i=1 ; i<K ; i++ ) os << " , " << mi[i];
		return os << " }";
	}


	template< unsigned int K >
	struct MultiIndices
	{
		struct iterator
		{
			iterator( MultiIndex< K > mi ) : mi(mi){}
			MultiIndex< K > mi;
			iterator &operator++ ( void ){ ++mi ; return *this; }
			iterator operator++ ( int ){ return iterator( mi++ ); }
			MultiIndex< K > &operator *( void ){ return mi; }
			bool operator != ( iterator iter ) const { return mi!=iter.mi; }
		};

		MultiIndices( unsigned int N ) : _N(N) { if( N<K ) MK_ERROR_OUT( "Number of indices cannot be less than size: " , K , " <= " , N ); }

		static unsigned int Choose( unsigned int n , unsigned int k )
		{
			if( k==0 ) return 1;
			else return ( Choose( n-1 , k-1 ) * n ) / k;
		}

		unsigned int size( void ) const{ return Choose( _N , K ); }

		iterator begin( void )
		{
			unsigned int idx[K];
			for( unsigned int i=0 ; i<K ; i++ ) idx[i] = i;
			return MultiIndex< K >( idx );
		}
		iterator end( void )
		{
			unsigned int idx[K];
			for( unsigned int i=0 ; i<K ; i++ ) idx[i] = _N-1-(K-1-i);
			return ++MultiIndex< K >( idx );
		}

	protected:
		unsigned int _N;

	};
}
#endif // MULTI_INDEX_INCLUDED