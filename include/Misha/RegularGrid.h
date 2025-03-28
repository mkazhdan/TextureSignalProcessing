/*
Copyright (c) 2019, Michael Kazhdan
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

#ifndef REGULAR_GRID_INCLUDED
#define REGULAR_GRID_INCLUDED

#include <type_traits>
#include "Array.h"
#include "Geometry.h"
#include "MultiThreading.h"
#include "Exceptions.h"



namespace MishaK
{
#define CLAMPED_EVALUATION
#ifdef CLAMPED_EVALUATION
#define FORCE_BILINEAR
#endif // CLAMPED_EVALUATION
	template< typename ... Type > struct RegularGridDataType{};

	template<>
	struct RegularGridDataType<>
	{
		static void Write( FILE *fp , unsigned int dim , std::string name );
		static bool Read ( FILE *fp , unsigned int dim , std::string name );
	};

	template<> struct RegularGridDataType< char >{ static const std::string Name ; static void Write( FILE *fp ){ RegularGridDataType<>::Write( fp , 1 , Name ); } ; static bool Read( FILE *fp ){ return RegularGridDataType<>::Read( fp , 1 , Name ); } };
	const std::string RegularGridDataType< char >::Name = "CHAR";
	template<> struct RegularGridDataType< unsigned char >{ static const std::string Name ; static void Write( FILE *fp ){ RegularGridDataType<>::Write( fp , 1 , Name ); } static bool Read( FILE *fp ){ return RegularGridDataType<>::Read( fp , 1 , Name ); } };
	const std::string RegularGridDataType< unsigned char >::Name = "UNSIGNED_CHAR";
	template<> struct RegularGridDataType< int >{ static const std::string Name ; static void Write( FILE *fp ){ RegularGridDataType<>::Write( fp , 1 , Name ); } static bool Read( FILE *fp ){ return RegularGridDataType<>::Read( fp , 1 , Name ); } };
	const std::string RegularGridDataType< int >::Name = "INT";
	template<> struct RegularGridDataType< unsigned int >{ static const std::string Name ; static void Write( FILE *fp ){ RegularGridDataType<>::Write( fp , 1 , Name ); } static bool Read( FILE *fp ){ return RegularGridDataType<>::Read( fp , 1 , Name ); } };
	const std::string RegularGridDataType< unsigned int >::Name = "UNSIGNED_INT";
	template<> struct RegularGridDataType< float >{ static const std::string Name ; static void Write( FILE *fp ){ RegularGridDataType<>::Write( fp , 1 , Name ); } static bool Read( FILE *fp ){ return RegularGridDataType<>::Read( fp , 1 , Name ); } };
	const std::string RegularGridDataType< float >::Name = "FLOAT";
	template<> struct RegularGridDataType< double >{ static const std::string Name ; static void Write( FILE *fp ){ RegularGridDataType<>::Write( fp , 1 , Name ); } static bool Read( FILE *fp ){ return RegularGridDataType<>::Read( fp , 1 , Name ); } };
	const std::string RegularGridDataType< double >::Name = "DOUBLE";

	template< typename Real , unsigned int Dim > struct RegularGridDataType< Point< Real , Dim > >{ static const std::string Name ; static void Write( FILE *fp ){ RegularGridDataType<>::Write( fp , Dim , Name ); } static bool Read( FILE *fp ){ return RegularGridDataType<>::Read( fp , Dim , Name ); } };
	template< typename Real , unsigned int Dim > const std::string RegularGridDataType< Point< Real , Dim > >::Name = RegularGridDataType< Real >::Name;

	template< unsigned int Dim , typename ... T > struct RegularGrid;

	template< unsigned int Dim >
	struct RegularGrid< Dim >
	{
		struct Index : public Point< int , Dim >
		{
			Index( void ){}
			Index( Point< int , Dim > p ) : Point< int , Dim >( p ){}
			template< typename ... Ints >
			Index( Ints ... is ) : Point< int , Dim >( is... ){}
			template< typename ... Indices > static Index Min( Index i , Indices ... is );
			template< typename ... Indices > static Index Max( Index i , Indices ... is );
			bool operator == ( Index i ) const;
			bool operator != ( Index i ) const;
			bool operator <  ( Index i ) const;
			bool operator <= ( Index i ) const;
			bool operator >  ( Index i ) const;
			bool operator >= ( Index i ) const;
		};

		struct Range : public std::pair< Index , Index >
		{
			using std::pair< Index , Index >::first;
			using std::pair< Index , Index >::second;

			Range( void );
			Range( Index I );
			template< typename ... Ranges > static Range Intersect( Ranges ... rs );
			Range dilate( unsigned int radius ) const;
			bool empty( void ) const;
			bool contains( Index I ) const;
			size_t size( void ) const;

			template< typename IndexFunctor /* = std::function< void ( Index ) > */ >
			void process( IndexFunctor f ) const { return this->template _process< 1 >( f ); }

			template< typename IndexFunctor /* = std::function< void ( Index ) > */  >
			void processParallel( IndexFunctor f ) const { return this->template _processParallel< 1 >( f ); }

			// IndexFunctor is a function taking in Count Index< Dim > arguments
			template< unsigned int Count , typename IndexFunctor /* = std::function< void ( Index ...  ) > */ >
			void process( IndexFunctor f ) const { return this->template _process< Count >( f ); }

			// IndexFunctor is a function taking in Count Index< Dim > arguments
			template< unsigned int Count , typename IndexFunctor /* = std::function< void ( unsigned int t , Index ...  ) > */  >
			void processParallel( IndexFunctor f ) const { return this->template _processParallel< Count >( f ); }

		protected:
			template< unsigned int Count , typename IndexFunctor /* = std::function< void ( Index ...  ) > */ , typename ... Indices /* = Index */ >
			void _process( IndexFunctor &f , Indices ... indices ) const;

			template< unsigned int Count , typename IndexFunctor /* = std::function< void ( unsigned int t , Index ...  ) > */ , typename ... Indices /* = Index */ >
			void _processParallel( IndexFunctor &f , Indices ... indices ) const;

			friend struct RegularGrid< Dim+1 >::Range;
		};

		static bool ReadDimension( std::string fileName , unsigned int &dim );
		static bool ReadHeader( std::string fileName , unsigned int &dataDim , std::string &dataName );

		static Index FactorIndex( size_t idx , const unsigned int res[Dim ] ){ Index I ; for( int d=0 ; d<Dim ; d++ ){ I[d] = idx % res[d] ; idx /= res[d]; }  return I; }
	};

	template< unsigned int Dim , typename DataType >
	struct RegularGrid< Dim , DataType >
	{
		RegularGrid( void ) : _values( NullPointer< DataType >() ){ for( unsigned int d=0 ; d<Dim ; d++ ) _res[d] = 0; }
		RegularGrid( RegularGrid &&grid ) : RegularGrid() { _Swap( *this , grid ); }
		RegularGrid( const RegularGrid &grid ) : RegularGrid()
		{
			resize( grid._res );
			for( size_t i=0 ; i<grid.resolution() ; i++ ) _values[i] = grid._values[i];
		}
		RegularGrid( const unsigned int *res ) : RegularGrid() { resize( res ); }
		RegularGrid &operator = ( RegularGrid &&grid ){ _Swap( *this , grid ) ; return *this; }
		RegularGrid &operator = ( const RegularGrid &grid )
		{
			if( this!=&grid )
			{
				resize( grid._res );
				for( size_t i=0 ; i<grid.resolution() ; i++ ) _values[i] = grid._values[i];
			}
			return *this;
		}
		~RegularGrid( void ){ DeletePointer( _values ); }

		template< typename Real > static void Read( std::string fileName , unsigned int res[Dim] , Pointer( DataType ) &values , XForm< Real , Dim+1 > &gridToModel );
		template< typename Real > static void Write( std::string fileName , const unsigned int res[Dim] , ConstPointer( DataType ) values , XForm< Real , Dim+1 > gridToModel );

		template< typename Real > void read( std::string fileName , XForm< Real , Dim+1 > &gridToModel );
		template< typename Real > void write( std::string fileName , XForm< Real , Dim+1 > gridToModel ) const;
		template< typename Real=DataType > void write( std::string fileName ) const;

		Pointer( DataType ) operator()( void ){ return _values; }
		ConstPointer( DataType ) operator()( void ) const { return _values; }

		DataType &operator[]( size_t idx ){ return _values[idx]; }
		const DataType &operator[]( size_t idx ) const { return _values[idx]; }

		template< typename Int > void factorIndex( size_t idx , Int indices[Dim] ) const { for( int d=0 ; d<Dim ; d++ ){ indices[d] = idx % _res[d] ; idx /= _res[d]; } }

		template< typename Int  >                    typename std::enable_if< std::is_integral< Int >::value , size_t >::type index( const Int coords[] ) const { return _index( coords , 1 ); }
		template< typename Int  >                    typename std::enable_if< std::is_integral< Int >::value , size_t >::type index(       Int coords[] ) const { return _index( coords , 1 ); }
		template< typename Int , typename ... Ints > typename std::enable_if< std::is_integral< Int >::value , size_t >::type index( Int coord , Ints ... coords ) const { static_assert( sizeof...(coords)+1==Dim , "[ERROR] number of coordinates does not match the number of dimensions" ) ; const Int c[] = { coord , coords ... } ; return index( c ); }

		template< typename Int >                     typename std::enable_if< std::is_integral< Int >::value >::type resize( const Int res[] );
		template< typename Int >                     typename std::enable_if< std::is_integral< Int >::value >::type resize(       Int res[] );
		template< typename Int , typename ... Ints > typename std::enable_if< std::is_integral< Int >::value >::type resize( Int res , Ints ...  ress );

		template< typename Int >                     typename std::enable_if< std::is_integral< Int >::value ,       DataType & >::type operator()( const Int coords[] )       { return operator[]( index( coords ) ); }
		template< typename Int >                     typename std::enable_if< std::is_integral< Int >::value ,       DataType & >::type operator()(       Int coords[] )       { return operator[]( index( coords ) ); }
		template< typename Int >                     typename std::enable_if< std::is_integral< Int >::value , const DataType & >::type operator()( const Int coords[] ) const { return operator[]( index( coords ) ); }
		template< typename Int >                     typename std::enable_if< std::is_integral< Int >::value , const DataType & >::type operator()(       Int coords[] ) const { return operator[]( index( coords ) ); }
		template< typename Int , typename ... Ints > typename std::enable_if< std::is_integral< Int >::value ,       DataType & >::type operator()( Int coord , Ints ... coords )       { static_assert( sizeof...(coords)+1==Dim , "[ERROR] number of coordinates does not match the number of dimensions" ) ; const Int c[] = { coord , coords ... } ; return operator()( c ); }
		template< typename Int , typename ... Ints > typename std::enable_if< std::is_integral< Int >::value , const DataType & >::type operator()( Int coord , Ints ... coords ) const { static_assert( sizeof...(coords)+1==Dim , "[ERROR] number of coordinates does not match the number of dimensions" ) ; const Int c[] = { coord , coords ... } ; return operator()( c ); }
		//	template< typename Int , typename ... Ints > typename std::enable_if< std::is_integral< Int >::value ,       DataType & >::type operator()( Point< Int , Dim > I )       { return operator()( &I[0] ); }
		//	template< typename Int , typename ... Ints > typename std::enable_if< std::is_integral< Int >::value , const DataType & >::type operator()( Point< Int , Dim > I ) const { return operator()( &I[0] ); }
		DataType &operator()( typename RegularGrid< Dim >::Index I ){ return operator()( &I[0] ); }
		const DataType &operator()( typename RegularGrid< Dim >::Index I ) const { return operator()( &I[0] ); }

#ifdef CLAMPED_EVALUATION
#ifdef FORCE_BILINEAR
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type operator()(       Real coords[] )       { return _Sample( _res , coords , _values ); }
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type operator()( const Real coords[] )       { return _Sample( _res , coords , _values ); }
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type operator()(       Real coords[] ) const { return _Sample( _res , coords , _values ); }
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type operator()( const Real coords[] ) const { return _Sample( _res , coords , _values ); }

		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type partial( unsigned int dir ,       Real coords[] )       { return _Partial( dir , _res , coords , _values ); }
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type partial( unsigned int dir , const Real coords[] )       { return _Partial( dir , _res , coords , _values ); }
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type partial( unsigned int dir ,       Real coords[] ) const { return _Partial( dir , _res , coords , _values ); }
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type partial( unsigned int dir , const Real coords[] ) const { return _Partial( dir , _res , coords , _values ); }
#else // !FORCE_BILINEAR
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type operator()(       Real coords[] , bool centered=false )       { return _Sample( _res , coords , centered , _values ); }
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type operator()( const Real coords[] , bool centered=false )       { return _Sample( _res , coords , centered , _values ); }
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type operator()(       Real coords[] , bool centered=false ) const { return _Sample( _res , coords , centered , _values ); }
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type operator()( const Real coords[] , bool centered=false ) const { return _Sample( _res , coords , centered , _values ); }

		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type partial( unsigned int dir ,       Real coords[] , bool centered=false )       { return _Partial( dir , _res , coords , centered , _values ); }
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type partial( unsigned int dir , const Real coords[] , bool centered=false )       { return _Partial( dir , _res , coords , centered , _values ); }
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type partial( unsigned int dir ,       Real coords[] , bool centered=false ) const { return _Partial( dir , _res , coords , centered , _values ); }
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type partial( unsigned int dir , const Real coords[] , bool centered=false ) const { return _Partial( dir , _res , coords , centered , _values ); }
#endif // FORCE_BILINEAR

		template< typename Real , typename ... Reals > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type operator()( Real coord , Reals ... coords  )       { static_assert( sizeof...(coords)+1==Dim , "[ERROR] number of coordinates does not match the number of dimensions" ) ; const Real c[] = { coord , coords ... } ; return operator()( c ); }
		template< typename Real , typename ... Reals > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type operator()( Real coord , Reals ... coords  ) const { static_assert( sizeof...(coords)+1==Dim , "[ERROR] number of coordinates does not match the number of dimensions" ) ; const Real c[] = { coord , coords ... } ; return operator()( c ); }

#ifdef FORCE_BILINEAR
		template< typename Real > DataType operator()( Point< Real , Dim > coords )       { return operator()( &coords[0] ); }
		template< typename Real > DataType operator()( Point< Real , Dim > coords ) const { return operator()( &coords[0] ); }

		template< typename Real > DataType partial( unsigned int dir , Point< Real , Dim > coords )       { return partial( dir , &coords[0] ); }
		template< typename Real > DataType partial( unsigned int dir , Point< Real , Dim > coords ) const { return partial( dir , &coords[0] ); }
#else // !FORCE_BILINEAR
		template< typename Real > DataType operator()( Point< Real , Dim > coords , bool centered=false )       { return operator()( &coords[0] , centered ); }
		template< typename Real > DataType operator()( Point< Real , Dim > coords , bool centered=false ) const { return operator()( &coords[0] , centered ); }

		template< typename Real > DataType partial( unsigned int dir , Point< Real , Dim > coords , bool centered=false )       { return partial( dir , &coords[0] , centered ); }
		template< typename Real > DataType partial( unsigned int dir , Point< Real , Dim > coords , bool centered=false ) const { return partial( dir , &coords[0] , centered ); }
#endif // FORCE_BILINEAR

		template< typename Real , typename ... Reals > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type partial( unsigned int dir , Real coord , Reals ... coords  )       { static_assert( sizeof...(coords)+1==Dim , "[ERROR] number of coordinates does not match the number of dimensions" ) ; const Real c[] = { coord , coords ... } ; return partial( dir , c ); }
		template< typename Real , typename ... Reals > typename std::enable_if< !std::is_integral< Real >::value , DataType >::type partial( unsigned int dir , Real coord , Reals ... coords  ) const { static_assert( sizeof...(coords)+1==Dim , "[ERROR] number of coordinates does not match the number of dimensions" ) ; const Real c[] = { coord , coords ... } ; return partial( dir , c ); }

#ifdef FORCE_BILINEAR
#else // !FORCE_BILINEAR
		template< typename Real > DataType partial( unsigned int dir , Point< Real , Dim > coords )       { return partial( dir , &coords[0] ); }
		template< typename Real > DataType partial( unsigned int dir , Point< Real , Dim > coords ) const { return partial( dir , &coords[0] ); }
#endif // FORCE_BILINEAR
#else // !CLAMPED_EVALUATION
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , ProjectiveData< Real , DataType > >::type operator()(       Real coords[] , bool centered=false )       { return _Sample( _res , coords , centered , _values ); }
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , ProjectiveData< Real , DataType > >::type operator()( const Real coords[] , bool centered=false )       { return _Sample( _res , coords , centered , _values ); }
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , ProjectiveData< Real , DataType > >::type operator()(       Real coords[] , bool centered=false ) const { return _Sample( _res , coords , centered , _values ); }
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , ProjectiveData< Real , DataType > >::type operator()( const Real coords[] , bool centered=false ) const { return _Sample( _res , coords , centered , _values ); }

		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , ProjectiveData< Real , DataType > >::type partial( unsigned int dir ,       Real coords[] , bool centered=false )       { return _Partial( dir , _res , coords , centered , _values ); }
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , ProjectiveData< Real , DataType > >::type partial( unsigned int dir , const Real coords[] , bool centered=false )       { return _Partial( dir , _res , coords , centered , _values ); }
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , ProjectiveData< Real , DataType > >::type partial( unsigned int dir ,       Real coords[] , bool centered=false ) const { return _Partial( dir , _res , coords , centered , _values ); }
		template< typename Real > typename std::enable_if< !std::is_integral< Real >::value , ProjectiveData< Real , DataType > >::type partial( unsigned int dir , const Real coords[] , bool centered=false ) const { return _Partial( dir , _res , coords , centered , _values ); }

		template< typename Real , typename ... Reals > typename std::enable_if< !std::is_integral< Real >::value , ProjectiveData< Real , DataType > >::type operator()( Real coord , Reals ... coords  )       { static_assert( sizeof...(coords)+1==Dim , "[ERROR] number of coordinates does not match the number of dimensions" ) ; const Real c[] = { coord , coords ... } ; return operator()( c ); }
		template< typename Real , typename ... Reals > typename std::enable_if< !std::is_integral< Real >::value , ProjectiveData< Real , DataType > >::type operator()( Real coord , Reals ... coords  ) const { static_assert( sizeof...(coords)+1==Dim , "[ERROR] number of coordinates does not match the number of dimensions" ) ; const Real c[] = { coord , coords ... } ; return operator()( c ); }

		template< typename Real > ProjectiveData< Real , DataType > operator()( Point< Real , Dim > coords , bool centered=false )       { return operator()( &coords[0] , centered ); }
		template< typename Real > ProjectiveData< Real , DataType > operator()( Point< Real , Dim > coords , bool centered=false ) const { return operator()( &coords[0] , centered ); }

		template< typename Real > ProjectiveData< Real , DataType > partial( unsigned int dir , Point< Real , Dim > coords , bool centered=false )       { return partial( dir , &coords[0] , centered ); }
		template< typename Real > ProjectiveData< Real , DataType > partial( unsigned int dir , Point< Real , Dim > coords , bool centered=false ) const { return partial( dir , &coords[0] , centered ); }

		template< typename Real , typename ... Reals > typename std::enable_if< !std::is_integral< Real >::value , ProjectiveData< Real , DataType > >::type partial( unsigned int dir , Real coord , Reals ... coords  )       { static_assert( sizeof...(coords)+1==Dim , "[ERROR] number of coordinates does not match the number of dimensions" ) ; const Real c[] = { coord , coords ... } ; return partial( dir , c ); }
		template< typename Real , typename ... Reals > typename std::enable_if< !std::is_integral< Real >::value , ProjectiveData< Real , DataType > >::type partial( unsigned int dir , Real coord , Reals ... coords  ) const { static_assert( sizeof...(coords)+1==Dim , "[ERROR] number of coordinates does not match the number of dimensions" ) ; const Real c[] = { coord , coords ... } ; return partial( dir , c ); }

		template< typename Real > ProjectiveData< Real , DataType > partial( unsigned int dir , Point< Real , Dim > coords )       { return partial( dir , &coords[0] ); }
		template< typename Real > ProjectiveData< Real , DataType > partial( unsigned int dir , Point< Real , Dim > coords ) const { return partial( dir , &coords[0] ); }
#endif // CLAMPED_EVALUATION

		template< typename Int >                     typename std::enable_if< std::is_integral< int >::value , bool >::type inBounds( const Int coords[] ) const { return _inBounds< Int , Dim-1 >( coords ); }
		template< typename Int >                     typename std::enable_if< std::is_integral< int >::value , bool >::type inBounds(       Int coords[] ) const { return _inBounds< Int , Dim-1 >( coords ); }
		template< typename Int , typename ... Ints > typename std::enable_if< std::is_integral< int >::value , bool >::type inBounds( Int coord , Ints ...  coords ) const { static_assert( sizeof...(coords)+1==Dim , "[ERROR] number of coordinates does not match the number of dimensions" ) ; const Int c[] = { coord , coords ... } ; return inBounds( c ); }

		const unsigned int *res( void ) const { return _res; }
		unsigned int res( unsigned int d ) const { return _res[d]; }
		size_t resolution( void ) const { return _Resolution< Dim >( _res ); }
		size_t size( void ) const { return _Resolution< Dim >( _res ); }

	protected:
		static void _Swap( RegularGrid &grid1 , RegularGrid &grid2 );

		template< unsigned int D=Dim > static typename std::enable_if< D==1 , size_t >::type _Resolution( const unsigned int res[] ) { return res[0]; }
		template< unsigned int D=Dim > static typename std::enable_if< D!=1 , size_t >::type _Resolution( const unsigned int res[] ) { return res[D-1] * _Resolution<D-1>(res); }

#ifdef CLAMPED_EVALUATION
#ifdef FORCE_BILINEAR
		template< typename Real , unsigned int D=Dim > static DataType _Sample( const unsigned int res[] , const Real coords[] , ConstPointer( DataType ) values );
		template< typename Real , unsigned int D=Dim > static DataType _Partial( unsigned int dir , const unsigned int res[] , const Real coords[] , ConstPointer( DataType ) values );
#else // !FORCE_BILINEAR
		template< typename Real , unsigned int D=Dim > static DataType _Sample( const unsigned int res[] , const Real coords[] , bool centered , ConstPointer( DataType ) values );
		template< typename Real , unsigned int D=Dim > static DataType _Partial( unsigned int dir , const unsigned int res[] , const Real coords[] , bool centered , ConstPointer( DataType ) values );
#endif // FORCE_BILINEAR
#else // !CLAMPED_EVALUATION
		template< typename Real , unsigned int D=Dim > static ProjectiveData< Real , DataType > _Sample( const unsigned int res[] , const Real coords[] , bool centered , ConstPointer( DataType ) values );
		template< typename Real , unsigned int D=Dim > static ProjectiveData< Real , DataType > _Partial( unsigned int dir , const unsigned int res[] , const Real coords[] , bool centered , ConstPointer( DataType ) values );
#endif // CLAMPED_EVALUATION

		template< typename Int , unsigned int D=0 > typename std::enable_if< D!=Dim-1 , size_t >::type _index(       Int coords[] , size_t stride ) const { return coords[0] * stride + _index< Int , D+1 >( coords+1 , stride * _res[D] ); }
		template< typename Int , unsigned int D=0 > typename std::enable_if< D==Dim-1 , size_t >::type _index(       Int coords[] , size_t stride ) const { return coords[0] * stride; }
		template< typename Int , unsigned int D=0 > typename std::enable_if< D!=Dim-1 , size_t >::type _index( const Int coords[] , size_t stride ) const { return coords[0] * stride + _index< Int , D+1 >( coords+1 , stride * _res[D] ); }
		template< typename Int , unsigned int D=0 > typename std::enable_if< D==Dim-1 , size_t >::type _index( const Int coords[] , size_t stride ) const { return coords[0] * stride; }

		template< typename Int , unsigned int D=Dim-1 > typename std::enable_if< D==0 , bool >::type _inBounds(       Int coords[] ) const { return coords[0]>=0 && coords[0]<(Int)_res[0];                                     }
		template< typename Int , unsigned int D=Dim-1 > typename std::enable_if< D!=0 , bool >::type _inBounds(       Int coords[] ) const { return coords[D]>=0 && coords[D]<(Int)_res[D] && _inBounds< Int , D-1 >( coords ); }
		template< typename Int , unsigned int D=Dim-1 > typename std::enable_if< D==0 , bool >::type _inBounds( const Int coords[] ) const { return coords[0]>=0 && coords[0]<(Int)_res[0];                                     }
		template< typename Int , unsigned int D=Dim-1 > typename std::enable_if< D!=0 , bool >::type _inBounds( const Int coords[] ) const { return coords[D]>=0 && coords[D]<(Int)_res[D] && _inBounds< Int , D-1 >( coords ); }


		unsigned int _res[Dim];
		Pointer( DataType ) _values;
	};
#include "RegularGrid.inl"
}
#endif // REGULAR_GRID_INCLUDED
