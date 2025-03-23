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

/////////////////////
// CmdLineReadable //
/////////////////////
inline CmdLineReadable::CmdLineReadable( const std::string &n ) : name(n) , set(false) {}

inline CmdLineReadable::~CmdLineReadable( void ){}

inline int CmdLineReadable::read( char ** , int ){ set = true ; return 0; }

template< typename Type >
Type CmdLineReadable::ToType( const std::string &str )
{
	Type type;
	std::stringstream( str ) >> type;
	return type;
}

template<>
std::string CmdLineReadable::ToType( const std::string &str ){ return str; }

//////////////////////
// CmdLineParameter //
//////////////////////
template< class Type > CmdLineParameter< Type >::CmdLineParameter( const std::string &name          ) : CmdLineReadable(name) { value = Type(); }

template< class Type > CmdLineParameter< Type >::CmdLineParameter( const std::string &name , Type v ) : CmdLineReadable(name) , value(v) {}

template< class Type >
int CmdLineParameter< Type >::read( char **argv , int argc )
{
	if( argc>0 )
	{
		value = ToType< Type >( argv[0] );
		set = true;
		return 1;
	}
	else return 0;
}

///////////////////////////
// CmdLineParameterArray //
///////////////////////////
template< class Type , int Dim >
CmdLineParameterArray< Type , Dim >::CmdLineParameterArray( const std::string &name , const Type* v ) : CmdLineReadable(name)
{
	if( v ) for( int i=0 ; i<Dim ; i++ ) values[i] = v[i];
	else    for( int i=0 ; i<Dim ; i++ ) values[i] = Type();
}

template< class Type , int Dim >
int CmdLineParameterArray< Type , Dim >::read( char **argv , int argc )
{
	if( argc>=Dim )
	{
		for( int i=0 ; i<Dim ; i++ ) values[i] = ToType< Type >( argv[i] );
		set = true;
		return Dim;
	}
	else return 0;
}

///////////////////////
// CmdLineParameters //
///////////////////////
template< class Type >
CmdLineParameters< Type >::CmdLineParameters( const std::string &name ) : CmdLineReadable(name) , values(NULL) , count(0) { }

template< class Type >
CmdLineParameters< Type >::~CmdLineParameters( void )
{
	if( values ) delete[] values;
	values = NULL;
	count = 0;
}

template< class Type >
int CmdLineParameters< Type >::read( char **argv , int argc )
{
	if( values ) delete[] values;
	values = NULL;

	if( argc>0 )
	{
		count = atoi(argv[0]);
		if( count<=0 || argc<=(int)count ) return 1;
		values = new Type[count];
		if( !values ) return 0;
		for( unsigned int i=0 ; i<count ; i++ ) values[i] = ToType< Type >( argv[i+1] );
		set = true;
		return count+1;
	}
	else return 0;
}

template< class Type >
void CmdLineParameters< Type >::resize( unsigned int sz )
{
	if( values ) delete[] values;
	values = NULL;
	count = sz;
	if( count ) values = new Type[count];
}

//////////////////////
// Helper functions //
//////////////////////

inline void CmdLineParse( int argc , char **argv , CmdLineReadable** params )
{
	while( argc>0 )
	{
		if( argv[0][0]=='-' && argv[0][1]=='-' )
		{
			CmdLineReadable* readable=NULL;
			for( int i=0 ; params[i]!=NULL && readable==NULL ; i++ ) if( params[i]->name==argv[0]+2 ) readable = params[i];
			if( readable )
			{
				int j = readable->read( argv+1 , argc-1 );
				argv += j , argc -= j;
			}
			else
			{
				MK_WARN( "Invalid option: " , argv[0] );
				for( int i=0 ; params[i]!=NULL ; i++ ) std::cerr << "\t--" << params[i]->name << std::endl;
			}
		}
		else MK_WARN( "Parameter name should be of the form --<name>: " , argv[0] );
		++argv , --argc;
	}
}

inline void CmdLineParse( int argc , char **argv , const std::vector< CmdLineReadable * > &params )
{
	while( argc>0 )
	{
		if( argv[0][0]=='-' && argv[0][1]=='-' )
		{
			CmdLineReadable* readable=NULL;
			for( int i=0 ; i<params.size() && readable==NULL ; i++ ) if( params[i]->name==argv[0]+2 ) readable = params[i];
			if( readable )
			{
				int j = readable->read( argv+1 , argc-1 );
				argv += j , argc -= j;
			}
			else
			{
				MK_WARN( "Invalid option: " , argv[0] );
				for( int i=0 ; i<params.size() ; i++ ) std::cerr << "\t--" << params[i]->name << std::endl;
			}
		}
		else MK_WARN( "Parameter name should be of the form --<name>: " , argv[0] );
		++argv , --argc;
	}
}

inline std::string ToUpper( const std::string &str )
{
	auto _ToUpper = []( char c ){ return c>='a' && c<='z' ? c+'A'-'a' : c; };
	std::string upper;
	upper.resize( str.size() );
	std::transform( str.begin() , str.end() , upper.begin() , _ToUpper );
	return upper;
}

inline std::string ToLower( const std::string &str )
{
	auto _ToLower = []( char c ){ return c>='A' && c<='Z' ? c+'a'-'A' : c; };
	std::string lower;
	lower.resize( str.size() );
	std::transform( str.begin() , str.end() , lower.begin() , _ToLower );
	return lower;
}

inline std::string GetFileExtension( const std::string &fileName )
{
	std::string ext;
	std::stringstream stream( fileName );
	while( std::getline( stream , ext , '.' ) ) ;
	return ext;
}

inline std::vector< std::string > ReadWords( const std::string &fileName )
{
	std::ifstream istream;
	istream.open( fileName );
	if( !istream ) MK_THROW( "Failed to open file for reading: " , fileName );
	std::vector< std::string > words;
	std::string word;
	while( istream >> word ) words.push_back( word );
	return words;
}

inline std::vector< std::string > GetWords( const std::string &str )
{
	std::stringstream sstream( str );
	std::vector< std::string > words;
	std::string word;
	while( sstream >> word ) words.push_back( word );
	return words;
}

inline std::vector< std::string > ReadLines( const std::string &fileName )
{
	std::ifstream istream;
	istream.open( fileName );
	if( !istream ) MK_THROW( "Failed to open file for reading: " , fileName );
	std::vector< std::string > lines;
	std::string line;
	while( std::getline( istream , line ) ) lines.push_back( line );
	return lines;
}
