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

/////////////////////////
// RegularGridTypeData //
/////////////////////////

void RegularGridDataType<>::Write( FILE *fp , unsigned int dim , std::string name ){ fprintf( fp , "%d %s\n" , (int)dim , name.c_str() ); }
bool RegularGridDataType<>::Read( FILE *fp , unsigned int dim , std::string name )
{
	char line[1024];
	int d;
#if _WIN32 || _WIN64
	if( fscanf_s( fp , " %d %s " , &d , line , (unsigned int)sizeof(line) )!=2 ) return false;
#else // !_WIN32 && !_WIN64
	if( fscanf( fp , " %d %s " , &d , line )!=2 ) return false;
#endif // _WIN32 || _WIN64
	return d==dim && name==std::string(line);
}

/////////////////
// RegularGrid //
/////////////////
template< unsigned int Dim >
bool RegularGrid< Dim >::ReadDimension( std::string fileName , unsigned int &dim )
{
	FILE *fp = fopen( fileName.c_str() , "rb" );
	if( !fp ) return false;
	else
	{
		// Read the magic number
		int d;
		if( fscanf( fp , " G%d " , &d )!=1 ){ fclose(fp) ; return false; }
		dim = d;
		fclose( fp );
		return true;
	}
}

template< unsigned int Dim >
bool RegularGrid< Dim >::ReadHeader( std::string fileName , unsigned int &dataDim , std::string &dataName )
{
	FILE *fp = fopen( fileName.c_str() , "rb" );
	if( !fp ) return false;
	else
	{
		// Read the magic number
		int d;
		if( fscanf( fp , " G%d " , &d )!=1 || d!=Dim ){ fclose(fp) ; return false; }

		char line[1024];
		if( fscanf( fp , " %d %s " , &d , line )!=2 ){ fclose(fp) ; return false; }
		dataDim = d , dataName = std::string( line );
		fclose( fp );
	}
	return true;
}

////////////////////////
// RegularGrid::Index //
////////////////////////
template< unsigned int Dim >
template< typename ... Indices >
typename RegularGrid< Dim >::Index RegularGrid< Dim >::Index::Min( typename RegularGrid< Dim >::Index i , Indices ... is )
{
	if constexpr( sizeof...(is)==0 ) return i;
	else
	{
		Index mn = Min( is ... );
		for( unsigned int d=0 ; d<Dim ; d++ ) mn[d] = std::min< int >( i[d] , mn[d] );
		return mn;
	}
}

template< unsigned int Dim >
template< typename ... Indices >
typename RegularGrid< Dim >::Index RegularGrid< Dim >::Index::Max( typename RegularGrid< Dim >::Index i , Indices ... is )
{
	if constexpr( sizeof...(is)==0 ) return i;
	else
	{
		Index mx = Max( is ... );
		for( unsigned int d=0 ; d<Dim ; d++ ) mx[d] = std::max< int >( i[d] , mx[d] );
		return mx;
	}
}

template< unsigned int Dim >
bool RegularGrid< Dim >::Index::operator == ( typename RegularGrid< Dim >::Index i ) const
{
	for( unsigned int d=0 ; d<Dim ; d++ ) if( (*this)[d]!=i[d] ) return false;
	return true;
}

template< unsigned int Dim >
bool RegularGrid< Dim >::Index::operator != ( typename RegularGrid< Dim >::Index i ) const { return !( (*this)==i ); }

template< unsigned int Dim >
bool RegularGrid< Dim >::Index::operator < ( typename RegularGrid< Dim >::Index i ) const
{
	for( unsigned int d=0 ; d<Dim ; d++ )
		if( (*this)[d]<i[d] ) return true;
		else if( (*this)[d]>i[d] ) return false;
	return false;
}

template< unsigned int Dim >
bool RegularGrid< Dim >::Index::operator <= ( typename RegularGrid< Dim >::Index i ) const
{
	for( unsigned int d=0 ; d<Dim ; d++ )
		if( (*this)[d]<i[d] ) return true;
		else if( (*this)[d]>i[d] ) return false;
	return true;
}

template< unsigned int Dim >
bool RegularGrid< Dim >::Index::operator > ( typename RegularGrid< Dim >::Index i ) const { return !(*this<=i); }

template< unsigned int Dim >
bool RegularGrid< Dim >::Index::operator >= ( typename RegularGrid< Dim >::Index i ) const { return !(*this<i); }

////////////////////////
// RegularGrid::Range //
////////////////////////
template< unsigned int Dim >
RegularGrid< Dim >::Range::Range( void ){}

template< unsigned int Dim >
RegularGrid< Dim >::Range::Range( Index I ){ for( unsigned int d=0 ; d<Dim ; d++ ) first[d] = I[d] , second[d] = I[d]+1; }

template< unsigned int Dim >
template< typename ... Ranges > 
typename RegularGrid< Dim >::Range RegularGrid< Dim >::Range::Intersect( Ranges ... rs )
{
	Range r;
	r.first = Index::Max( rs.first ... );
	r.second = Index::Min( rs.second ... );
	return r;
}

template< unsigned int Dim >
typename RegularGrid< Dim >::Range RegularGrid< Dim >::Range::dilate( unsigned int radius ) const
{
	Range r;
	for( unsigned int d=0 ; d<Dim ; d++ ) r.first[d] = first[d]-radius , r.second[d] = second[d] + radius;
	return r;
}

template< unsigned int Dim >
bool RegularGrid< Dim >::Range::empty( void ) const
{
	for( unsigned int d=0 ; d<Dim ; d++ ) if( first[d]>=second[d] ) return true;
	return false;
}

template< unsigned int Dim >
size_t RegularGrid< Dim >::Range::size( void ) const
{
	size_t s = 1;
	for( unsigned int d=0 ; d<Dim ; d++ )
		if( first[d]>=second[d] ) return 0;
		else s *= second[d]-first[d];
	return s;
}

template< unsigned int Dim >
bool RegularGrid< Dim >::Range::contains( Index I ) const
{
	for( unsigned int d=0 ; d<Dim ; d++ ) if( I[d]<first[d] || I[d]>=second[d] ) return false;
	return true;
}

template< unsigned int Dim >
template< unsigned int Count , typename IndexFunctor /* = std::function< void ( Index ...  ) > */ , typename ... Indices /* = Index */ >
void RegularGrid< Dim >::Range::_process( IndexFunctor &f , Indices ... indices ) const
{
	if constexpr( Count==1 )
	{
		if constexpr( Dim==1 ) for( int i=first[0] ; i<second[0] ; i++ ) f( indices ... , Index(i) );
		else
		{
			typename RegularGrid< Dim-1 >::Range _r;
			for( unsigned int d=0 ; d<Dim-1 ; d++ ) _r.first[d] = first[d+1] , _r.second[d] = second[d+1];

			Index idx;
			auto _f = [&]( typename RegularGrid< Dim-1 >::Index _idx )
				{
					for( unsigned int d=0 ; d<Dim-1 ; d++ ) idx[d+1] = _idx[d];
					f( indices ... , idx );
				};
			for( int i=first[0] ; i<second[0] ; i++ )
			{
				idx[0] = i;
				_r.template _process< 1 >( _f );
			}
		}
	}
	else
	{
		auto _f = [&]( Index I ){ this->template _process< Count-1 >( f , indices ... , I ); };
		this->template _process< 1 >( _f );
	}
}

template< unsigned int Dim >
template< unsigned int Count , typename IndexFunctor /* = std::function< void ( unsigned int , Index ...  ) > */ , typename ... Indices /* = Index */ >
void RegularGrid< Dim >::Range::_processParallel( IndexFunctor &f , Indices ... indices ) const
{
	if constexpr( Count==1 )
	{
		if constexpr( Dim==1 ) ThreadPool::ParallelFor( first[0] , second[0] , [&]( unsigned int t , size_t i ){ f( t , indices ... , Index(i) ); } );
		else
		{
			typename RegularGrid< Dim-1 >::Range _r;
			for( unsigned int d=0 ; d<Dim-1 ; d++ ) _r.first[d] = first[d+1] , _r.second[d] = second[d+1];

			ThreadPool::ParallelFor
				(
					first[0] , second[0] ,
					[&]( unsigned int t , size_t i )
					{
						Index idx;
						auto _f = [&]( typename RegularGrid< Dim-1 >::Index _idx )
							{
								for( unsigned int d=0 ; d<Dim-1 ; d++ ) idx[d+1] = _idx[d];
								f( t , indices ... , idx );
							};
						idx[0] = (int)i;
						_r.template _process< 1 >( _f );
					}
				);
		}
	}
	else
	{
		MK_WARN_ONCE( "This could be parallelized" );
		auto _f = [&]( Index I ){ this->template _process< Count-1 >( 0 , f , indices ... , I ); };
		this->template _process< 1 >( _f );
	}
}


/////////////////
// RegularGrid //
/////////////////

template< unsigned int Dim , typename DataType >
void RegularGrid< Dim , DataType >::_Swap( RegularGrid &grid1 , RegularGrid &grid2 )
{
	{
		unsigned int res;
		for( int d=0 ; d<Dim ; d++ ){ res = grid1._res[d] ; grid1._res[d] = grid2._res[d] ; grid2._res[d] = res; }
	}
	{
		Pointer( DataType ) values;
		values = grid1._values ; grid1._values = grid2._values ; grid2._values = values;
	}
}

template< unsigned int Dim , typename DataType >
template< typename Int >
typename std::enable_if< std::is_integral< Int >::value >::type RegularGrid< Dim , DataType >::resize( Int res[] )
{
	if( _values ) DeletePointer( _values );
	size_t resolution = 1;
	for( int d=0 ; d<Dim ; d++ ) _res[d] = (unsigned int)res[d] , resolution *= (size_t)res[d];
	if( resolution ) _values = NewPointer< DataType >( resolution );
}

template< unsigned int Dim , typename DataType >
template< typename Int >
typename std::enable_if< std::is_integral< Int >::value >::type RegularGrid< Dim , DataType >::resize( const Int res[] )
{
	if( _values ) DeletePointer( _values );
	size_t resolution = 1;
	for( int d=0 ; d<Dim ; d++ ) _res[d] = (unsigned int)res[d] , resolution *= (size_t)res[d];
	if( resolution ) _values = NewPointer< DataType >( resolution );
}

template< unsigned int Dim , typename DataType >
template< typename Int , typename ... Ints >
typename std::enable_if< std::is_integral< Int >::value >::type RegularGrid< Dim , DataType >::resize( Int res , Ints ... ress )
{
	if constexpr( sizeof...(ress)==0 )
	{
		unsigned int r[Dim];
		for( unsigned int d=0 ; d<Dim ; d++ ) r[d] = res;
		return resize( r );
	}
	else
	{
		static_assert( sizeof...(ress)+1==Dim , "[ERROR] number of resolutions does not match the number of dimensions" );
		const Int r[] = { res , ress ... };
		return resize( r );
	}
}

template< unsigned int Dim , typename DataType >
template< typename Real , unsigned int D >
#ifdef CLAMPED_EVALUATION
#ifdef FORCE_BILINEAR
DataType RegularGrid< Dim , DataType >::_Sample( const unsigned int res[] , const Real coords[] , ConstPointer( DataType ) values )
#else // !FORCE_BILINEAR
DataType RegularGrid< Dim , DataType >::_Sample( const unsigned int res[] , const Real coords[] , bool centered , ConstPointer( DataType ) values )
#endif // FORCE_BILINEAR
#else // !CLAMPED_EVALUATION
ProjectiveData< Real , DataType > RegularGrid< Dim , DataType >::_Sample( const unsigned int res[] , const Real coords[] , bool centered , ConstPointer( DataType ) values )
#endif // CLAMPED_EVALUATION
{
#ifdef FORCE_BILINEAR
	Real coord = coords[D-1];
#else // !FORCE_BILINEAR
	Real coord = centered ? coords[D-1]-(Real)0.5 : coords[D-1];
#endif // FORCE_BILINEAR
	int iCoord1 = (int)floor(coord) , iCoord2 = (int)floor(coord)+1;
	Real dx1 = (Real)( iCoord2 - coord ) , dx2 = (Real)( coord - iCoord1 );

#ifdef CLAMPED_EVALUATION
	iCoord1 = std::max< int >( 0 , std::min< int >( res[D-1]-1 , iCoord1 ) );
	iCoord2 = std::max< int >( 0 , std::min< int >( res[D-1]-1 , iCoord2 ) );
	if constexpr( D==1 ) return values[ iCoord1 ] * dx1 + values[ iCoord2 ] * dx2;
#ifdef FORCE_BILINEAR
	else                 return _Sample< Real , D-1 >( res , coords , values + _Resolution< D-1 >(res) * iCoord1 ) * dx1 + _Sample< Real , D-1 >( res , coords , values + _Resolution< D-1 >(res) * iCoord2 ) * dx2;
#else // !FORCE_BILINEAR
	else                 return _Sample< Real , D-1 >( res , coords , centered , values + _Resolution< D-1 >(res) * iCoord1 ) * dx1 + _Sample< Real , D-1 >( res , coords , centered , values + _Resolution< D-1 >(res) * iCoord2 ) * dx2;
#endif // FORCE_BILINEAR
#else // !CLAMPED_EVALUATION
	ProjectiveData< Real , DataType > d;
	if constexpr( D==1 )
	{
		if( iCoord1>=0 && iCoord1<(int)res[0] ) d += ProjectiveData< Real , DataType >( values[ iCoord1 ] * dx1 , dx1 );
		if( iCoord2>=0 && iCoord2<(int)res[0] ) d += ProjectiveData< Real , DataType >( values[ iCoord2 ] * dx2 , dx2 );
		return d;
	}
	else
	{
		if( iCoord1>=0 && iCoord1<(int)res[D-1] ) d += _Sample< Real , D-1 >( res , coords , centered , values + _Resolution< D-1 >(res) * iCoord1 ) * dx1;
		if( iCoord2>=0 && iCoord2<(int)res[D-1] ) d += _Sample< Real , D-1 >( res , coords , centered , values + _Resolution< D-1 >(res) * iCoord2 ) * dx2;
		return d;
	}
#endif // CLAMPED_EVALUATION
}

template< unsigned int Dim , typename DataType >
template< typename Real , unsigned int D >
#ifdef CLAMPED_EVALUATION
#ifdef FORCE_BILINEAR
DataType RegularGrid< Dim , DataType >::_Partial( unsigned int dir , const unsigned int res[] , const Real coords[] , ConstPointer( DataType ) values )
#else // !FORCE_BILINEAR
DataType RegularGrid< Dim , DataType >::_Partial( unsigned int dir , const unsigned int res[] , const Real coords[] , bool centered , ConstPointer( DataType ) values )
#endif // FORCE_BILINEAR
#else // !CLAMPED_EVALUATION
ProjectiveData< Real , DataType > RegularGrid< Dim , DataType >::_Partial( unsigned int dir , const unsigned int res[] , const Real coords[] , bool centered , ConstPointer( DataType ) values )
#endif // CLAMPED_EVALUATION
{
#ifdef FORCE_BILINEAR
	Real coord = coords[D-1];
#else // !FORCE_BILINEAR
	Real coord = centered ? coords[D-1]-(Real)0.5 : coords[D-1];
#endif // FORCE_BILINEAR
	int iCoord1 = (int)floor(coord) , iCoord2 = (int)floor(coord)+1;
	Real dx1 = (Real)( iCoord2 - coord ) , dx2 = (Real)( coord - iCoord1 );

#ifdef CLAMPED_EVALUATION
	iCoord1 = std::max< int >( 0 , std::min< int >( res[D-1]-1 , iCoord1 ) );
	iCoord2 = std::max< int >( 0 , std::min< int >( res[D-1]-1 , iCoord2 ) );
	if constexpr( D==1 )
	{
		if( dir==0 ) return values[ iCoord2 ] - values[ iCoord1 ];
		else         return values[ iCoord1 ] * dx1 + values[ iCoord2 ] * dx2;
	}
	else
	{
#ifdef FORCE_BILINEAR
		if( dir==D-1 ) return _Partial< Real , D-1 >( dir , res , coords , values + _Resolution< D-1 >(res) * iCoord2 ) - _Partial< Real , D-1 >( dir , res , coords , values + _Resolution< D-1 >(res) * iCoord1 );
		else           return _Partial< Real , D-1 >( dir , res , coords , values + _Resolution< D-1 >(res) * iCoord1 ) * dx1 + _Partial< Real , D-1 >( dir , res , coords , values + _Resolution< D-1 >(res) * iCoord2 ) * dx2;
#else // !FORCE_BILINEAR
		if( dir==D-1 ) return _Partial< Real , D-1 >( dir , res , coords , centered , values + _Resolution< D-1 >(res) * iCoord2 ) - _Partial< Real , D-1 >( dir , res , coords , centered , values + _Resolution< D-1 >(res) * iCoord1 );
		else           return _Partial< Real , D-1 >( dir , res , coords , centered , values + _Resolution< D-1 >(res) * iCoord1 ) * dx1 + _Partial< Real , D-1 >( dir , res , coords , centered , values + _Resolution< D-1 >(res) * iCoord2 ) * dx2;
#endif // FORCE_BILINEAR
	}
#else // !CLAMPED_EVALUATION
	ProjectiveData< Real , DataType > data;
	if constexpr( D==1 )
	{
		if( dir==0 )
		{
			if( iCoord1>=0 && iCoord1<(int)res[0] ) data -= ProjectiveData< Real , DataType >( values[ iCoord1 ] , (Real)1 );
			if( iCoord2>=0 && iCoord2<(int)res[0] ) data += ProjectiveData< Real , DataType >( values[ iCoord2 ] , (Real)1. );
		}
		else
		{
			if( iCoord1>=0 && iCoord1<(int)res[0] ) data += ProjectiveData< Real , DataType >( values[ iCoord1 ] * dx1 , dx1 );
			if( iCoord2>=0 && iCoord2<(int)res[0] ) data += ProjectiveData< Real , DataType >( values[ iCoord2 ] * dx2 , dx2 );
		}
	}
	else
	{
		if( dir==D-1 )
		{
			if( iCoord1>=0 && iCoord1<(int)res[D-1] ) data -= _Partial< Real , D-1 >( dir , res , coords , centered , values + _Resolution< D-1 >(res) * iCoord1 );
			if( iCoord2>=0 && iCoord2<(int)res[D-1] ) data += _Partial< Real , D-1 >( dir , res , coords , centered , values + _Resolution< D-1 >(res) * iCoord2 );
		}
		else
		{
			if( iCoord1>=0 && iCoord1<(int)res[D-1] ) data += _Partial< Real , D-1 >( dir , res , coords , centered , values + _Resolution< D-1 >(res) * iCoord1 ) * dx1;
			if( iCoord2>=0 && iCoord2<(int)res[D-1] ) data += _Partial< Real , D-1 >( dir , res , coords , centered , values + _Resolution< D-1 >(res) * iCoord2 ) * dx2;
		}
	}
	return data;
#endif // CLAMPED_EVALUATION
}

template< unsigned int Dim , typename DataType >
template< typename Real >
void RegularGrid< Dim , DataType >::Write( std::string fileName , const unsigned int res[Dim] , ConstPointer( DataType ) values , XForm< Real , Dim+1 > gridToModel )
{
	FILE *fp = fopen( fileName.c_str() , "wb" );
	if( !fp ) MK_ERROR_OUT( "Failed to open grid file for writing: " , fileName );
	else
	{
		// Write the magic number
		fprintf( fp , "G%d\n" , (int)Dim );

		RegularGridDataType< DataType >::Write( fp );

		// Write the dimensions
		for( int d=0 ; d<Dim ; d++ )
		{
			fprintf( fp , "%d" , (int)res[d] );
			if( d==Dim-1 ) fprintf( fp , "\n" );
			else           fprintf( fp , " " );
		}

		// Write the transformation
		for( int j=0 ; j<Dim+1 ; j++ ) for( int i=0 ; i<Dim+1 ; i++ )
		{
			fprintf( fp , "%f" , (float)gridToModel(i,j) );
			if( i==Dim ) fprintf( fp , "\n" );
			else         fprintf( fp , " " );
		}

		// Write the grid values
		fwrite( values , sizeof(DataType) , _Resolution(res) , fp );
		fclose( fp );
	}
}

template< unsigned int Dim , typename DataType >
template< typename Real >
void RegularGrid< Dim , DataType >::write( std::string fileName ) const
{
	XForm< Real , Dim+1 > gridToModel = XForm< Real , Dim+1 >::Identity();
	Write( fileName , _res , _values , gridToModel );
}


template< unsigned int Dim , typename DataType >
template< typename Real >
void RegularGrid< Dim , DataType >::write( std::string fileName , XForm< Real , Dim+1 > gridToModel ) const
{
	Write( fileName , _res , _values , gridToModel );
}

template< unsigned int Dim , typename DataType >
template< typename Real >
void RegularGrid< Dim , DataType >::Read( std::string fileName , unsigned int res[Dim] , Pointer( DataType ) &values , XForm< Real , Dim+1 > &gridToModel )
{
	FILE *fp = fopen( fileName.c_str() , "rb" );
	if( !fp ) MK_ERROR_OUT( "Failed to open grid file for reading: " , fileName );
	else
	{
		// Read the magic number
		{
			int dim;
			if( fscanf( fp , " G%d " , &dim )!=1 ) MK_ERROR_OUT( "Failed to read magic number: " , fileName );
			if( dim!=Dim ) MK_ERROR_OUT( "Dimensions don't match: " , Dim , " != " , dim );
		}

		// Read the data type
		if( !RegularGridDataType< DataType >::Read( fp ) ) MK_ERROR_OUT( "Failed to read type" );

		// Read the dimensions
		{
			int r;
			for( int d=0 ; d<Dim ; d++ )
			{
				if( fscanf( fp , " %d " , &r )!=1 ) MK_ERROR_OUT( "Failed to read dimension[ " , d , " ]" );
				res[d] = r;
			}
		}

		// Read the transformation
		{
			float x;
			for( int j=0 ; j<Dim+1 ; j++ ) for( int i=0 ; i<Dim+1 ; i++ )
			{
				if( fscanf( fp , " %f" , &x )!=1 ) MK_ERROR_OUT( "Failed to read xForm( " , i , " , " , j , " )" );
				gridToModel(i,j) = (Real)x;
			}
		}

		// Read through the end of the line
		{
			char line[1024];
			if( !fgets( line , sizeof(line)/sizeof(char) , fp ) ) MK_ERROR_OUT( "Could not read end of line" );
		}

		values = NewPointer< DataType >( _Resolution(res) );
		// Read the grid values
		if( fread( values , sizeof(DataType) , _Resolution(res) , fp )!=_Resolution(res) ) MK_ERROR_OUT( "Failed to read values" );
		fclose( fp );
	}
}


template< unsigned int Dim , typename DataType >
template< typename Real >
void RegularGrid< Dim , DataType >::read( std::string fileName , XForm< Real , Dim+1 > &gridToModel )
{
	Read( fileName , _res , _values , gridToModel );
}
