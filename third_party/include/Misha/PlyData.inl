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

inline int ToPlyType( TypeOnDisk typeOnDisk )
{
	switch( typeOnDisk )
	{
	case TypeOnDisk::CHAR:   return PLY::Type<          char >();
	case TypeOnDisk::UCHAR:  return PLY::Type< unsigned char >();
	case TypeOnDisk::INT:    return PLY::Type<           int >();
	case TypeOnDisk::UINT:   return PLY::Type< unsigned  int >();
	case TypeOnDisk::FLOAT:  return PLY::Type<         float >();
	case TypeOnDisk::DOUBLE: return PLY::Type<        double >();
	default: MK_ERROR_OUT( "Unrecognized type: " , typeOnDisk );
	}
	return -1;
}

inline TypeOnDisk FromPlyType( int plyType )
{
	switch( plyType )
	{
	case PLY_INT:    return TypeOnDisk::INT;
	case PLY_UINT:   return TypeOnDisk::UINT;
	case PLY_CHAR:   return TypeOnDisk::CHAR;
	case PLY_UCHAR:  return TypeOnDisk::UCHAR;
	case PLY_FLOAT:  return TypeOnDisk::FLOAT;
	case PLY_DOUBLE: return TypeOnDisk::DOUBLE;
	default: MK_ERROR_OUT( "Unrecognized type: " , plyType );
	}
	return TypeOnDisk::UNKNOWN;
}

template<> inline TypeOnDisk GetTypeOnDisk<          char >( void ){ return TypeOnDisk::CHAR;  }
template<> inline TypeOnDisk GetTypeOnDisk< unsigned char >( void ){ return TypeOnDisk::UCHAR; }
template<> inline TypeOnDisk GetTypeOnDisk<           int >( void ){ return TypeOnDisk::INT;  }
template<> inline TypeOnDisk GetTypeOnDisk< unsigned  int >( void ){ return TypeOnDisk::UINT; }
template<> inline TypeOnDisk GetTypeOnDisk<         float >( void ){ return TypeOnDisk::FLOAT;  }
template<> inline TypeOnDisk GetTypeOnDisk<        double >( void ){ return TypeOnDisk::DOUBLE; }

template< typename Real >
template< typename Type >
bool DataIO< Real >::_ReadBinary( FILE *fp , Real &r )
{
	Type t;
	if( fread( &t , sizeof(Type) , 1 , fp )!=1 ) return false;
	r = (Real)t;
	return true;
}

template< typename Real >
template< typename Type >
void DataIO< Real >::_WriteBinary( FILE *fp , Real r ){ Type t = (Type)r ; fwrite(  &t , sizeof(Type) , 1 , fp ); }

template< typename Real >
bool DataIO< Real >::ReadASCII( FILE *fp , TypeOnDisk , Real &s )
{
	double d;
#if _WIN32 || _WIN64
	if( fscanf_s( fp , " %lf"  , &d ) ) return false;
#else // !_WIN32 && !_WIN64
	if( fscanf( fp , " %lf"  , &d )!=1 ) return false;
#endif // _WIN32 || _WIN64
	s = (Real)d;
	return true;
}

template< typename Real >
bool DataIO< Real >::ReadBinary( FILE *fp , TypeOnDisk typeOnDisk , Real &s )
{
	if( TypeOnDisk()==typeOnDisk ) return fread( &s , sizeof(Real) , 1 , fp )==1;
	switch( typeOnDisk )
	{
	case TypeOnDisk::CHAR:   return _ReadBinary<          char >( fp , s );
	case TypeOnDisk::UCHAR:  return _ReadBinary< unsigned char >( fp , s );
	case TypeOnDisk::INT:    return _ReadBinary<           int >( fp , s );
	case TypeOnDisk::UINT:   return _ReadBinary< unsigned  int >( fp , s );
	case TypeOnDisk::FLOAT:  return _ReadBinary<         float >( fp , s );
	case TypeOnDisk::DOUBLE: return _ReadBinary<        double >( fp , s );
	default: MK_ERROR_OUT( "Unrecognized type: " , typeOnDisk );
	}
	return true;
}
template< typename Real >
void DataIO< Real >::WriteASCII( FILE *fp , TypeOnDisk typeOnDisk , const Real &s )
{
	switch( typeOnDisk )
	{
	case TypeOnDisk::CHAR:   fprintf( fp , " %d" , (         char)s ) ; break;
	case TypeOnDisk::UCHAR:  fprintf( fp , " %u" , (unsigned char)s ) ; break;
	case TypeOnDisk::INT:    fprintf( fp , " %d" , (          int)s ) ; break;
	case TypeOnDisk::UINT:   fprintf( fp , " %u" , (unsigned  int)s ) ; break;
	case TypeOnDisk::FLOAT:  fprintf( fp , " %f" , (        float)s ) ; break;
	case TypeOnDisk::DOUBLE: fprintf( fp , " %f" , (       double)s ) ; break;
	default: MK_ERROR_OUT( "Unrecongized type: " , typeOnDisk );
	}
}

template< typename Real >
void DataIO< Real >::WriteBinary( FILE *fp , TypeOnDisk typeOnDisk , const Real &s )
{
	if( TypeOnDisk()==typeOnDisk ) fwrite( &s , sizeof(Real) , 1 , fp );
	switch( typeOnDisk )
	{
	case TypeOnDisk::CHAR:   _WriteBinary<          char >( fp , s ) ; break;
	case TypeOnDisk::UCHAR:  _WriteBinary< unsigned char >( fp , s ) ; break;
	case TypeOnDisk::INT:    _WriteBinary<           int >( fp , s ) ; break;
	case TypeOnDisk::UINT:   _WriteBinary< unsigned  int >( fp , s ) ; break;
	case TypeOnDisk::FLOAT:  _WriteBinary<         float >( fp , s ) ; break;
	case TypeOnDisk::DOUBLE: _WriteBinary<        double >( fp , s ) ; break;
	default: MK_ERROR_OUT( "Unrecongized type: " , typeOnDisk );
	}
}

template< typename Real >
bool DataIO< Real >::ReadASCII( FILE *fp , TypeOnDisk typeOnDisk , size_t sz , Real *s )
{
	for( size_t i=0 ; i<sz ; i++ ) if( !ReadASCII( fp , typeOnDisk , s[i] ) ) return false;
	return true;
}

template< typename Real >
bool DataIO< Real >::ReadBinary( FILE *fp , TypeOnDisk typeOnDisk , size_t sz , Real *s )
{
#if defined( _WIN32 ) || defined( _WIN64 )
#pragma warning( push )
#pragma warning( disable : 4068 )
#endif // _WIN32 || _WIN64
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wattribute-warning"
	if( TypeOnDisk()==typeOnDisk ) return fread( &s , sizeof(Real) , sz , fp )==sz;
	else for( size_t i=0 ; i<sz ; i++ ) if( !ReadBinary( fp , typeOnDisk , s[i] ) ) return false;
#pragma GCC diagnostic pop
#if defined( _WIN32 ) || defined( _WIN64 )
#pragma warning( pop )
#endif // _WIN32 || _WIN64
	return true;
}
template< typename Real >
void DataIO< Real >::WriteASCII( FILE *fp , TypeOnDisk typeOnDisk , size_t sz , const Real *s )
{
	for( size_t i=0 ; i<sz ; i++ ) WriteASCII( fp , typeOnDisk , s[i] );
}

template< typename Real >
void DataIO< Real >::WriteBinary( FILE *fp , TypeOnDisk typeOnDisk , size_t sz , const Real *s )
{
	if( TypeOnDisk()==typeOnDisk ) fwrite( &s , sizeof(Real) , sz , fp );
	else for( size_t i=0 ; i<sz ; i++ ) WriteBinary( fp , typeOnDisk , s[i] );
}

//////////////
// _Factory //
//////////////
template< typename _DataType , typename FactoryType >
bool _Factory< _DataType , FactoryType >::plyValidReadProperties( const std::vector< GregTurk::PlyProperty > &plyProperties ) const
{
	for( unsigned int i=0 ; i<plyReadNum() ; i++ )
	{
		std::string name = plyReadProperty(i).name;
		bool hasProperty = false;
		for( unsigned int j=0 ; j<plyProperties.size() ; j++ ) hasProperty |= name==plyProperties[j].name;
		if( !hasProperty ) return false;
	}
	return true;
}

//////////////////
// PointFactory //
//////////////////
template< typename Real , unsigned int Dim >
std::string PointFactory< Real , Dim >::_name( unsigned int idx ) const { return _header + std::string( "_" ) + std::to_string(idx); }

template< typename Real , unsigned int Dim >
GregTurk::PlyProperty PointFactory< Real , Dim >::plyReadProperty( unsigned int idx ) const
{
	if( idx>=plyReadNum() ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _name(idx) , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
}
template< typename Real , unsigned int Dim >
GregTurk::PlyProperty PointFactory< Real , Dim >::plyWriteProperty( unsigned int idx ) const
{
	if( idx>=plyWriteNum() ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _name(idx) , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
}

template< typename Real , unsigned int Dim >
GregTurk::PlyProperty PointFactory< Real , Dim >::plyStaticReadProperty( unsigned int idx ) const
{
	if( idx>=plyReadNum() ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _name(idx) , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( DataType , coords ) + sizeof(Real)*idx );
}
template< typename Real , unsigned int Dim >
GregTurk::PlyProperty PointFactory< Real , Dim >::plyStaticWriteProperty( unsigned int idx ) const
{
	if( idx>=plyWriteNum() ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _name(idx) , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( DataType , coords ) + sizeof(Real)*idx );
}

///////////////////
// MatrixFactory //
///////////////////
template< typename Real , unsigned int Cols , unsigned int Rows >
std::string MatrixFactory< Real , Cols , Rows >::_name( unsigned int idx ) const
{
	unsigned int c = idx / Rows;
	unsigned int r = idx % Rows;
	return _header + std::string( "_" ) + std::to_string(c) + std::string( "_" ) + std::to_string(r);
}
template< typename Real , unsigned int Cols , unsigned int Rows >
GregTurk::PlyProperty MatrixFactory< Real , Cols , Rows >::plyReadProperty( unsigned int idx ) const
{
	if( idx>=plyReadNum() ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _name(idx) , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
}
template< typename Real , unsigned int Cols , unsigned int Rows >
GregTurk::PlyProperty MatrixFactory< Real , Cols , Rows >::plyWriteProperty( unsigned int idx ) const
{
	if( idx>=plyWriteNum() ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _name(idx) , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
}

template< typename Real , unsigned int Cols , unsigned int Rows >
GregTurk::PlyProperty MatrixFactory< Real , Cols , Rows >::plyStaticReadProperty( unsigned int idx ) const
{
	if( idx>=plyReadNum() ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _name(idx) , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( DataType , coords ) + sizeof(Real)*idx );
}
template< typename Real , unsigned int Cols , unsigned int Rows >
GregTurk::PlyProperty MatrixFactory< Real , Cols , Rows >::plyStaticWriteProperty( unsigned int idx ) const
{
	if( idx>=plyWriteNum() ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _name(idx) , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( DataType , coords ) + sizeof(Real)*idx );
}

/////////////////////
// PositionFactory //
/////////////////////
template< typename Real , unsigned int Dim >
GregTurk::PlyProperty PositionFactory< Real , Dim >::plyReadProperty( unsigned int idx ) const
{
	if( idx>= plyReadNum() ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
}
template< typename Real , unsigned int Dim >
GregTurk::PlyProperty PositionFactory< Real , Dim >::plyWriteProperty( unsigned int idx ) const
{
	if( idx>=plyWriteNum() ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
}

template< typename Real , unsigned int Dim >
GregTurk::PlyProperty PositionFactory< Real , Dim >::plyStaticReadProperty( unsigned int idx ) const
{
	if( idx>=plyReadNum() ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( DataType , coords ) + sizeof(Real)*idx );
}
template< typename Real , unsigned int Dim >
GregTurk::PlyProperty PositionFactory< Real , Dim >::plyStaticWriteProperty( unsigned int idx ) const
{
	if( idx>=plyWriteNum() ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( DataType , coords ) + sizeof(Real)*idx );
}

template<> inline const std::string PositionFactory<  float , 1 >::_PlyNames[] = { "x" };
template<> inline const std::string PositionFactory< double , 1 >::_PlyNames[] = { "x" };
template<> inline const std::string PositionFactory<  float , 2 >::_PlyNames[] = { "x" , "y" };
template<> inline const std::string PositionFactory< double , 2 >::_PlyNames[] = { "x" , "y" };
template<> inline const std::string PositionFactory<  float , 3 >::_PlyNames[] = { "x" , "y" , "z" };
template<> inline const std::string PositionFactory< double , 3 >::_PlyNames[] = { "x" , "y" , "z" };
template<> inline const std::string PositionFactory<  float , 4 >::_PlyNames[] = { "x" , "y" , "z" , "w" };
template<> inline const std::string PositionFactory< double , 4 >::_PlyNames[] = { "x" , "y" , "z" , "w" };

///////////////////
// NormalFactory //
///////////////////
template< typename Real , unsigned int Dim >
GregTurk::PlyProperty NormalFactory< Real , Dim >::plyReadProperty( unsigned int idx ) const
{
	if( idx>= plyReadNum() ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
}
template< typename Real , unsigned int Dim >
GregTurk::PlyProperty NormalFactory< Real , Dim >::plyWriteProperty( unsigned int idx ) const
{
	if( idx>=plyWriteNum() ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
}

template< typename Real , unsigned int Dim >
GregTurk::PlyProperty NormalFactory< Real , Dim >::plyStaticReadProperty( unsigned int idx ) const
{
	if( idx>=plyReadNum() ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( DataType , coords ) + sizeof(Real)*idx );
}
template< typename Real , unsigned int Dim >
GregTurk::PlyProperty NormalFactory< Real , Dim >::plyStaticWriteProperty( unsigned int idx ) const
{
	if( idx>=plyWriteNum() ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( DataType , coords ) + sizeof(Real)*idx );
}

template<> inline const std::string NormalFactory<  float , 1 >::_PlyNames[] = { "nx" };
template<> inline const std::string NormalFactory< double , 1 >::_PlyNames[] = { "nx" };
template<> inline const std::string NormalFactory<  float , 2 >::_PlyNames[] = { "nx" , "ny" };
template<> inline const std::string NormalFactory< double , 2 >::_PlyNames[] = { "nx" , "ny" };
template<> inline const std::string NormalFactory<  float , 3 >::_PlyNames[] = { "nx" , "ny" , "nz" };
template<> inline const std::string NormalFactory< double , 3 >::_PlyNames[] = { "nx" , "ny" , "nz" };
template<> inline const std::string NormalFactory<  float , 4 >::_PlyNames[] = { "nx" , "ny" , "nz" , "nw" };
template<> inline const std::string NormalFactory< double , 4 >::_PlyNames[] = { "nx" , "ny" , "nz" , "nw" };

////////////////////
// TextureFactory //
////////////////////
template< typename Real , unsigned int Dim >
GregTurk::PlyProperty TextureFactory< Real , Dim >::plyReadProperty( unsigned int idx ) const
{
	if( idx>= plyReadNum() ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
}
template< typename Real , unsigned int Dim >
GregTurk::PlyProperty TextureFactory< Real , Dim >::plyWriteProperty( unsigned int idx ) const
{
	if( idx>=plyWriteNum() ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
}

template< typename Real , unsigned int Dim >
GregTurk::PlyProperty TextureFactory< Real , Dim >::plyStaticReadProperty( unsigned int idx ) const
{
	if( idx>=plyReadNum() ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( DataType , coords ) + sizeof(Real)*idx );
}
template< typename Real , unsigned int Dim >
GregTurk::PlyProperty TextureFactory< Real , Dim >::plyStaticWriteProperty( unsigned int idx ) const
{
	if( idx>=plyWriteNum() ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( DataType , coords ) + sizeof(Real)*idx );
}

template<> inline const std::string TextureFactory<  float , 1 >::_PlyNames[] = { "u" };
template<> inline const std::string TextureFactory< double , 1 >::_PlyNames[] = { "u" };
template<> inline const std::string TextureFactory<  float , 2 >::_PlyNames[] = { "u" , "v" };
template<> inline const std::string TextureFactory< double , 2 >::_PlyNames[] = { "u" , "v" };
template<> inline const std::string TextureFactory<  float , 3 >::_PlyNames[] = { "u" , "v" , "w" };
template<> inline const std::string TextureFactory< double , 3 >::_PlyNames[] = { "u" , "v" , "w" };

/////////////////////
// RGBColorFactory //
/////////////////////
template< typename Real >
GregTurk::PlyProperty RGBColorFactory< Real >::plyReadProperty( unsigned int idx ) const
{
	if( idx>= plyReadNum() ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
}
template< typename Real >
GregTurk::PlyProperty RGBColorFactory< Real >::plyWriteProperty( unsigned int idx ) const
{
	if( idx>=plyWriteNum() ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
}

template< typename Real >
GregTurk::PlyProperty RGBColorFactory< Real >::plyStaticReadProperty( unsigned int idx ) const
{
	if( idx>=plyReadNum() ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( DataType , coords ) + sizeof(Real)*idx );
}
template< typename Real >
GregTurk::PlyProperty RGBColorFactory< Real >::plyStaticWriteProperty( unsigned int idx ) const
{
	if( idx>=plyWriteNum() ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( DataType , coords ) + sizeof(Real)*idx );
}

template< typename Real >
bool RGBColorFactory< Real >::plyValidReadProperties( const std::vector< GregTurk::PlyProperty > &plyProperties ) const
{
	for( int d=0 ; d<3 ; d++ )
	{
		std::string name1 = plyReadProperty(d).name , name2 = plyReadProperty(d+3).name;
		bool hasProperty = false;
		for( unsigned int j=0 ; j<plyProperties.size() ; j++ ) hasProperty |= plyProperties[j].name==name1 || plyProperties[j].name==name2;
		if( !hasProperty ) return false;
	}
	return true;
}


template< typename Real > const std::string RGBColorFactory< Real >::_PlyNames[] = { "red" , "green" , "blue" , "r" , "g" , "b" };

//////////////////////
// RGBAColorFactory //
//////////////////////
template< typename Real >
GregTurk::PlyProperty RGBAColorFactory< Real >::plyReadProperty( unsigned int idx ) const
{
	if( idx>= plyReadNum() ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
}
template< typename Real >
GregTurk::PlyProperty RGBAColorFactory< Real >::plyWriteProperty( unsigned int idx ) const
{
	if( idx>=plyWriteNum() ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
}

template< typename Real >
GregTurk::PlyProperty RGBAColorFactory< Real >::plyStaticReadProperty( unsigned int idx ) const
{
	if( idx>=plyReadNum() ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( DataType , coords ) + sizeof(Real)*idx );
}
template< typename Real >
GregTurk::PlyProperty RGBAColorFactory< Real >::plyStaticWriteProperty( unsigned int idx ) const
{
	if( idx>=plyWriteNum() ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( DataType , coords ) + sizeof(Real)*idx );
}

template< typename Real >
bool RGBAColorFactory< Real >::plyValidReadProperties( const std::vector< GregTurk::PlyProperty > &plyProperties ) const
{
	for( int d=0 ; d<4 ; d++ )
	{
		std::string name1 = plyReadProperty(d).name , name2 = plyReadProperty(d+4).name;
		bool hasProperty = false;
		for( unsigned int j=0 ; j<plyProperties.size() ; j++ ) hasProperty |= plyProperties[j].name==name1 || plyProperties[j].name==name2;
		if( !hasProperty ) return false;
	}
	return true;
}

template< typename Real > const std::string RGBAColorFactory< Real >::_PlyNames[] = { "red" , "green" , "blue" , "alpha" , "r" , "g" , "b" , "a" };

//////////////////
// ValueFactory //
//////////////////
template< typename Real >
GregTurk::PlyProperty ValueFactory< Real >::plyReadProperty( unsigned int idx ) const
{
	if( idx>= plyReadNum() ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
}
template< typename Real >
GregTurk::PlyProperty ValueFactory< Real >::plyWriteProperty( unsigned int idx ) const
{
	if( idx>=plyWriteNum() ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
}

template< typename Real >
GregTurk::PlyProperty ValueFactory< Real >::plyStaticReadProperty( unsigned int idx ) const
{
	if( idx>=plyReadNum() ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , 0 );
}
template< typename Real >
GregTurk::PlyProperty ValueFactory< Real >::plyStaticWriteProperty( unsigned int idx ) const
{
	if( idx>=plyWriteNum() ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _PlyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , 0 );
}

template< typename Real > const std::string ValueFactory< Real >::_PlyNames[] = { "value" };

///////////////////
// StaticFactory //
///////////////////
template< typename Real , unsigned int Dim >
StaticFactory< Real , Dim >::StaticFactory( const std::string plyNames[] , TypeOnDisk typeOnDisk ) : _typeOnDisk( typeOnDisk )
{
	for( unsigned int d=0 ; d<Dim ; d++ ) _plyNames[d] = plyNames[d];
}

template< typename Real , unsigned int Dim >
GregTurk::PlyProperty StaticFactory< Real , Dim >::plyReadProperty( unsigned int idx ) const
{
	if( idx>=Dim ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _plyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
}
template< typename Real , unsigned int Dim >
GregTurk::PlyProperty StaticFactory< Real , Dim >::plyWriteProperty( unsigned int idx ) const
{
	if( idx>=Dim ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _plyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , sizeof(Real)*idx );
}

template< typename Real , unsigned int Dim >
GregTurk::PlyProperty StaticFactory< Real , Dim >::plyStaticReadProperty( unsigned int idx ) const
{
	if( idx>=Dim ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _plyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( DataType , coords ) + sizeof(Real)*idx );
}
template< typename Real , unsigned int Dim >
GregTurk::PlyProperty StaticFactory< Real , Dim >::plyStaticWriteProperty( unsigned int idx ) const
{
	if( idx>=Dim ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _plyNames[idx] , ToPlyType( _typeOnDisk ) , PLY::Type< Real >() , (int)offsetof( DataType , coords ) + sizeof(Real)*idx );
}

template< typename Real , unsigned int Dim >
bool StaticFactory< Real , Dim >::operator == ( const StaticFactory< Real , Dim > &factory ) const
{
	if( _typeOnDisk!=factory._typeOnDisk ) return false;
	for( int d=0 ; d<Dim ; d++ ) if( _plyNames[d]!=factory._plyNames[d] ) return false;
	return true;
}

////////////////////
// DynamicFactory //
////////////////////
template< typename Real >
DynamicFactory< Real >::DynamicFactory( const std::vector< std::pair< std::string , TypeOnDisk > > &namesAndTypesOnDisk ) : _namesAndTypesOnDisk(namesAndTypesOnDisk)
{
}
template< typename Real >
DynamicFactory< Real >::DynamicFactory( const std::vector< GregTurk::PlyProperty > &plyProperties )
{
	_namesAndTypesOnDisk.resize( plyProperties.size() );
	for( int i=0 ; i<plyProperties.size() ; i++ ) _namesAndTypesOnDisk[i] = std::pair< std::string , TypeOnDisk >( plyProperties[i].name , FromPlyType( plyProperties[i].external_type ) );
}

template< typename Real >
bool DynamicFactory< Real >::readASCII( FILE *fp , DataType &dt ) const
{
	for( unsigned int i=0 ; i<dt.dim() ; i++ ) if( !DataIO< Real >::ReadASCII( fp , _namesAndTypesOnDisk[i].second , dt[i] ) ) return false;
	return true;
}
template< typename Real >
bool DynamicFactory< Real >::readBinary( FILE *fp , DataType &dt ) const
{
	for( unsigned int i=0 ; i<dt.dim() ; i++ ) if( !DataIO< Real >::ReadBinary( fp , _namesAndTypesOnDisk[i].second , dt[i] ) ) return false;
	return true;
}
template< typename Real >
void DynamicFactory< Real >::writeASCII( FILE *fp , const DataType &dt ) const
{
	for( unsigned int i=0 ; i<dt.dim() ; i++ ) DataIO< Real >::WriteASCII( fp , _namesAndTypesOnDisk[i].second , dt[i] );
}
template< typename Real >
void DynamicFactory< Real >::writeBinary( FILE *fp , const DataType &dt ) const
{
	for( unsigned int i=0 ; i<dt.dim() ; i++ ) DataIO< Real >::WriteBinary( fp , _namesAndTypesOnDisk[i].second , dt[i] );
}

template< typename Real >
GregTurk::PlyProperty DynamicFactory< Real >::plyReadProperty( unsigned int idx ) const
{
	if( idx>=plyReadNum() ) MK_ERROR_OUT( "read property out of bounds" );
	return GregTurk::PlyProperty( _namesAndTypesOnDisk[idx].first , ToPlyType( _namesAndTypesOnDisk[idx].second ) , PLY::Type< Real >() , sizeof(Real)*idx );
}
template< typename Real >
GregTurk::PlyProperty DynamicFactory< Real >::plyWriteProperty( unsigned int idx ) const
{
	if( idx>=plyWriteNum() ) MK_ERROR_OUT( "write property out of bounds" );
	return GregTurk::PlyProperty( _namesAndTypesOnDisk[idx].first , ToPlyType( _namesAndTypesOnDisk[idx].second ) , PLY::Type< Real >() , sizeof(Real)*idx );
}

template< typename Real >
bool DynamicFactory< Real >::operator == ( const DynamicFactory< Real > &factory ) const
{
	if( size()!=factory.size() ) return false;
	for( int i=0 ; i<size() ; i++ ) if( _namesAndTypesOnDisk[i].first!=factory._namesAndTypesOnDisk[i].first || _namesAndTypesOnDisk[i].second!=factory._namesAndTypesOnDisk[i].second ) return false;
	return true;
}

/////////////
// Factory //
/////////////
template< typename Real , typename ... Factories >
template< unsigned int I >
typename std::enable_if< I!=sizeof...(Factories) , GregTurk::PlyProperty >::type Factory< Real , Factories ... >::_plyReadProperty( unsigned int idx , size_t offset ) const
{
	if( idx<this->template get<I>().plyReadNum() )
	{
		GregTurk::PlyProperty prop = this->template get<I>().plyReadProperty(idx);
		prop.offset += (int)offset;
		return prop;
	}
	else return _plyReadProperty<I+1>( idx - this->template get<I>().plyReadNum() , offset + this->template get<I>().bufferSize() );
}

template< typename Real , typename ... Factories >
template< unsigned int I >
typename std::enable_if< I!=sizeof...(Factories) , GregTurk::PlyProperty >::type Factory< Real , Factories ... >::_plyWriteProperty( unsigned int idx , size_t offset ) const
{
	if( idx<this->template get<I>().plyWriteNum() )
	{
		GregTurk::PlyProperty prop = this->template get<I>().plyWriteProperty(idx);
		prop.offset += (int)offset;
		return prop;
	}
	else return _plyWriteProperty<I+1>( idx - this->template get<I>().plyWriteNum() , offset + this->template get<I>().bufferSize() );
}

template< typename Real , typename ... Factories >
template< unsigned int I >
typename std::enable_if< I!=sizeof...(Factories) , GregTurk::PlyProperty >::type Factory< Real , Factories ... >::_plyStaticReadProperty( unsigned int idx ) const
{
	if( idx<this->template get<I>().plyReadNum() )
	{
		DataType v;
		GregTurk::PlyProperty prop = this->template get<I>().plyStaticReadProperty( idx );
		prop.offset += (int)( (size_t)&v.template get<I>() - (size_t)&v );
		return prop;
	}
	else return _plyStaticReadProperty<I+1>( idx - this->template get<I>().plyReadNum() );
}

template< typename Real , typename ... Factories >
template< unsigned int I >
typename std::enable_if< I!=sizeof...(Factories) , GregTurk::PlyProperty >::type Factory< Real , Factories ... >::_plyStaticWriteProperty( unsigned int idx ) const
{
	if( idx<this->template get<I>().plyWriteNum() )
	{
		DataType v;
		GregTurk::PlyProperty prop = this->template get<I>().plyStaticWriteProperty( idx );
		prop.offset += (int)( (size_t)&v.template get<I>() - (size_t)&v );
		return prop;
	}
	else return _plyStaticWriteProperty<I+1>( idx - this->template get<I>().plyWriteNum() ) ;
}

