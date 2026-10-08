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

//////////////////////
// Node::_StateInfo //
//////////////////////
struct Node::_StateInfo
{
	friend struct Node;

	enum struct NodeType
	{
		VARIABLE ,
		CONSTANT ,
		UNARY_OPERATOR ,
		BINARY_OPERATOR ,
		FUNCTION ,
		UNKNOWN ,
		L_PARENTH ,
		R_PARENTH ,
		COMMA
	};

	static inline const std::string NodeTypeNames[] =
	{
		"variable" ,
		"constant" ,
		"operator (unary)" ,
		"operator (binary)" ,
		"function" ,
		"unknown" ,
		"parentheses (left)" ,
		"parentheses (right)" ,
		"comma"
	};

	struct State
	{
		NodeType type;
		std::string name;

		State( NodeType type=NodeType::UNKNOWN , std::string name="" );
		Node operator()( const Node & node ) const;
		Node operator()( const Node & node1 , const Node & node2 ) const;
	};

	inline static const std::vector< std::vector< std::string > > UnaryOperators = { { "-" } };
	inline static const std::vector< std::vector< std::string > > BinaryOperators = { { "^" } , { "*" , "/" } , { "+" , "-" } };

	std::vector< State > state;

	_StateInfo( void );
	_StateInfo( std::string eqn , const std::vector< std::string > & vars );

	unsigned int openingParenth( unsigned int idx ) const;
	unsigned int closingParenth( unsigned int idx ) const;
	std::vector< _StateInfo > splitOnCommas( void ) const;


	_StateInfo sub( unsigned int begin , unsigned end ) const;

	void addParenths( void );
	void addFunctionParenths( void );
	void addUnaryOperatorParenths( const std::vector< std::string > & ops );
	void addBinaryOperatorParenths( const std::vector< std::string > & ops );

	static bool IsParenth( char c );
	static bool IsComma( char c );
	static bool IsOperator( char c );
	static bool IsOperator( std::string str );
	static bool IsOperator( char c , const std::vector< std::string > & ops );
	static bool IsOperator( std::string str , const std::vector< std::string > & ops );
	static bool IsUnaryOperator( char c );
	static bool IsUnaryOperator( std::string str );
	static bool IsBinaryOperator( char c );
	static bool IsBinaryOperator( std::string str );

	static void Print( const std::vector< State > & state , bool showType , bool indent=false );
	static unsigned int OpeningParenth( const std::vector< State > & state , unsigned int idx );
	static unsigned int ClosingParenth( const std::vector< State > & state , unsigned int idx );
};

/////////////////////////////
// Node::_StateInfo::State //
/////////////////////////////
Node::_StateInfo::State::State( NodeType type , std::string name ) : type(type) , name(name){}

Node Node::_StateInfo::State::operator()( const Node & node ) const
{
	switch( type )
	{
	case NodeType::UNARY_OPERATOR:
		if( name=="-" ) return - node;
		else MK_THROW( "Unrecognized unary operator: " , name );
		break;
	case NodeType::FUNCTION:
		if     ( name=="exp"  ) return Exp ( node );
		else if( name=="cos"  ) return Cos ( node );
		else if( name=="sin"  ) return Sin ( node );
		else if( name=="tan"  ) return Tan ( node );
		else if( name=="cosh" ) return Cosh( node );
		else if( name=="sinh" ) return Sinh( node );
		else if( name=="tanh" ) return Tanh( node );
		else if( name=="sqrt" ) return Sqrt( node );
		else if( name=="sgn" ) return Sgn( node );
		else if( name=="abs" ) return Abs( node );
		else if( name=="relu" ) return ReLU( node );
		else MK_THROW( "Unrecognized unary function: " , name );
		break;
	default:
		MK_THROW( "Node type is not unary operator or unary function: " , NodeTypeNames[ static_cast< int >( type ) ] );
	}
	return Node();
}

Node Node::_StateInfo::State::operator()( const Node & node1 , const Node & node2 ) const
{
	switch( type )
	{
	case NodeType::BINARY_OPERATOR:
		if     ( name=="+" ) return node1 + node2;
		else if( name=="-" ) return node1 + (-node2);
		else if( name=="*" ) return node1 * node2;
		else if( name=="/" ) return node1 * Pow( node2 , -1. );
		else if( name=="^" ) return Pow( node1 , node2 );
		else MK_THROW( "Unrecognized binary operator: " , name );
		break;
	case NodeType::FUNCTION:
		if     ( name=="min" ) return Min( node1 , node2 );
		else if( name=="max" ) return Max( node1 , node2 );
		else MK_THROW( "Unrecognized binary function: " , name );
		break;
	default:
		MK_THROW( "Node type is not binary operator: " , NodeTypeNames[ static_cast< int >( type ) ] );
	}
	return Node();
}

//////////////////////
// Node::_StateInfo //
//////////////////////

inline Node::_StateInfo::_StateInfo( void ){}

inline Node::_StateInfo::_StateInfo( std::string eqn , const std::vector< std::string > & vars )
{
	auto PrintNodeTypes = []( const std::vector< NodeType > & state , const std::vector< std::string > & tokens )
	{
		for( unsigned int i=0 ; i<state.size() ; i++ )
			std::cout << "\t" << NodeTypeNames[ static_cast< int >( state[i] ) ] << " : " << tokens[i] << std::endl;
	};

	auto PrintStates = []( const std::vector< State > & state , bool showType )
	{
		for( unsigned int i=0 ; i<state.size() ; i++ )
			if( showType ) std::cout << "\t" << NodeTypeNames[ static_cast< int >( state[i].type ) ] << " : " << state[i].name << std::endl;
			else           std::cout << " " << state[i].name;
		if( !showType ) std::cout << std::endl;
	};


	// Add spaces around operators and parentheses
	{
		std::string _eqn;
		for( unsigned int i=0 ; i<eqn.size() ; i++ )
		{
			if( IsOperator( eqn[i] ) || IsParenth( eqn[i] ) || IsComma( eqn[i] ) )
			{
				_eqn.push_back( ' ' );
				_eqn.push_back( eqn[i] );
				_eqn.push_back( ' ' );
			}
			else _eqn.push_back( eqn[i] );
		}
		eqn = _eqn;
	}

	// Tokenize on spaces
	std::vector< std::string > tokens;
	{
		auto RemoveLeadingWhiteSpace = []( const std::string str )
		{
			unsigned int idx = 0;
			while( idx<str.size() && std::isspace( str[idx] ) ) idx++;
			return str.substr( idx );
		};
		eqn = RemoveLeadingWhiteSpace( eqn );

		while( eqn.size() )
		{
			std::string token;
			for( unsigned int i=0 ; i<eqn.size() && !std::isspace( eqn[i] ) ; i++ ) token += eqn[i];
			if( !token.size() ) MK_THROW( "Expected non-empty token" );
			tokens.push_back( token );
			eqn = RemoveLeadingWhiteSpace( eqn.substr( token.size() ) );
		}
	}

	// Identify what we can
	std::vector< NodeType > _state( tokens.size() , NodeType::UNKNOWN );
	{
		int pCount = 0;
		for( unsigned int i=0 ; i<tokens.size() ; i++ )
		{
			if( IsOperator( tokens[i] ) ) _state[i] = NodeType::BINARY_OPERATOR;
			else if( tokens[i]=="(" || tokens[i]==")" )
			{
				if     ( tokens[i]=="(" ) pCount++ , _state[i] = NodeType::L_PARENTH;
				else if( tokens[i]==")" ) pCount-- , _state[i] = NodeType::R_PARENTH;
				if( pCount<0 ){ MK_THROW( "Negative parentheses count: " , eqn ); }
			}
			else if( tokens[i]=="," ) _state[i] = NodeType::COMMA;
		}
		if( pCount ){ MK_THROW( "Non-zero parentheses count" ); }
	}

	for( unsigned int i=0 ; i<tokens.size() ; i++ ) if( _state[i]==NodeType::UNKNOWN ) for( unsigned int j=0 ; j<vars.size() ; j++ ) if( tokens[i]==vars[j] )
		if( ( i==0 || _state[i-1]==NodeType::BINARY_OPERATOR || _state[i-1]==NodeType::L_PARENTH || _state[i-1]==NodeType::COMMA ) && ( i+1==tokens.size() || _state[i+1]==NodeType::BINARY_OPERATOR || _state[i+1]==NodeType::R_PARENTH || _state[i+1]==NodeType::COMMA ) )
		{
			_state[i] = NodeType::VARIABLE;
			break;
		}

	for( unsigned int i=0 ; i<_state.size() ; i++ )
	{
		if( _state[i]==NodeType::UNKNOWN )
		{
			try
			{
				if( tokens[i]=="Pi" ) state.emplace_back( NodeType::CONSTANT , tokens[i] );
				else
				{
					double foo = std::stod( tokens[i] );
					state.emplace_back( NodeType::CONSTANT , tokens[i] );
				}
			}
			catch( ... ){ state.emplace_back( NodeType::FUNCTION , tokens[i] ); }
		}
		else state.emplace_back( _state[i] , tokens[i] );
	}

	for( unsigned int i=0 ; i<state.size() ; i++ ) if( state[i].type==NodeType::FUNCTION && ( i+1>=state.size() || state[i+1].type!=NodeType::L_PARENTH ) )
		MK_THROW( "Expected left parenth after unary function: " , i , "] " , tokens[i] );


	addFunctionParenths();
	for( unsigned int i=0 ; i<UnaryOperators.size() ; i++ ) addUnaryOperatorParenths( UnaryOperators[i] );
	for( unsigned int i=0 ; i<BinaryOperators.size() ; i++ ) addBinaryOperatorParenths( BinaryOperators[i] );

	for( unsigned int i=0 ; i<state.size() ; i++ ) if( state[i].type==NodeType::UNKNOWN ) MK_THROW( "Failed to parse: " , i );
}

inline unsigned int Node::_StateInfo::openingParenth( unsigned int idx ) const
{
#ifdef NEW_EQUATION_PARSER
	return OpeningParenth( state , idx );
#else // !NEW_EQUATION_PARSER
	if( state[idx].type!=NodeType::R_PARENTH ) MK_THROW( "Expected right parenth: " , idx );
	int pCount = 1;
	for( unsigned int i=idx ; i>0 ; i-- )
	{
		if     ( state[i-1].type==NodeType::L_PARENTH ) pCount--;
		else if( state[i-1].type==NodeType::R_PARENTH ) pCount++;
		if( pCount==0 ) return i-1;
	}
	MK_THROW( "Could not find opening parenth: " , idx );
	return static_cast< unsigned int >(-1);
#endif // NEW_EQUATION_PARSER
}

inline unsigned int Node::_StateInfo::closingParenth( unsigned int idx ) const
{
#ifdef NEW_EQUATION_PARSER
	return ClosingParenth( state , idx );
#else // !NEW_EQUATION_PARSER
	if( state[idx].type!=NodeType::L_PARENTH ) MK_THROW( "Expected left parenth: " , idx );
	int pCount = 1;
	for( unsigned int i=idx+1 ; i<state.size() ; i++ )
	{
		if     ( state[i].type==NodeType::L_PARENTH ) pCount++;
		else if( state[i].type==NodeType::R_PARENTH ) pCount--;
		if( pCount==0 ) return i;
	}
	MK_THROW( "Could not find closing parenth: " , idx );
	return static_cast< unsigned int >(-1);
#endif // NEW_EQUATION_PARSER
}

inline Node::_StateInfo Node::_StateInfo::sub( unsigned int begin , unsigned end ) const
{
	_StateInfo stateInfo;
	for( unsigned int i=begin ; i<end ; i++ ) stateInfo.state.emplace_back( state[i] );
	return stateInfo;
}

inline void Node::_StateInfo::addParenths( void )
{
#ifdef NEW_EQUATION_PARSER
	if( state[0].type!=NodeType::L_PARENTH || closingParenth( 0 )!=state.size()-1 )
#else // !NEW_EQUATION_PARSER
	if( state[0].type!=NodeType::L_PARENTH )
#endif // NEW_EQUATION_PARSER
	{
		std::vector< State > _state;
		_state.emplace_back( NodeType::L_PARENTH , "(" );
		for( unsigned int j=0 ; j<state.size() ; j++ ) _state.push_back( state[j] );
		_state.emplace_back( NodeType::R_PARENTH , ")" );
		state = _state;
	}
}

inline void Node::_StateInfo::addFunctionParenths( void )
{
	for( unsigned int i=0 ; i<state.size() ; i++ )
		if( state[i].type==NodeType::FUNCTION && ( i==0 || ( i>0 && state[i-1].type!=NodeType::L_PARENTH ) ) )
		{
			unsigned int idx = closingParenth( i+1 );

			std::vector< State > _state;
			for( unsigned int j=0 ; j<i ; j++ ) _state.push_back( state[j] );
			_state.emplace_back( NodeType::L_PARENTH , "(" );
			for( unsigned int j=i ; j<idx ; j++ ) _state.push_back( state[j] );
			_state.emplace_back( NodeType::R_PARENTH , ")" );
			for( unsigned int j=idx ; j<state.size() ; j++ ) _state.push_back( state[j] );
			state = _state;
			addFunctionParenths();
		}
}

inline void Node::_StateInfo::addUnaryOperatorParenths( const std::vector< std::string > & ops )
{
	for( unsigned int i=0 ; i<state.size() ; i++ )
	{
		if( state[i].type==NodeType::BINARY_OPERATOR && IsOperator( state[i].name , ops ) )
		{
			if( i && ( state[i-1].type==NodeType::VARIABLE || state[i-1].type==NodeType::CONSTANT || state[i-1].type==NodeType::R_PARENTH ) ) continue;
			unsigned int end;

			if( state[i+1].type==NodeType::VARIABLE || state[i+1].type==NodeType::CONSTANT ) end = i+1;
			else if( state[i+1].type==NodeType::L_PARENTH ) end = closingParenth( i+1 );
			else MK_THROW( "Expected variable or left parenth" );

			if( i==0 || state[i-1].type!=NodeType::L_PARENTH || closingParenth(i-1)!=end+1 )
			{
				std::vector< State > _state;
				for( unsigned int j=0 ; j<i ; j++ ) _state.push_back( state[j] );
				_state.emplace_back( NodeType::L_PARENTH , "(" );
				for( unsigned int j=i ; j<=end ; j++ ) _state.push_back( state[j] );
				_state.emplace_back( NodeType::R_PARENTH , ")" );
				for( unsigned int j=end+1 ; j<state.size() ; j++ ) _state.push_back( state[j] );
				_state[i+1].type = NodeType::UNARY_OPERATOR;
				state = _state;
				addUnaryOperatorParenths( ops );
			}
			else state[i].type = NodeType::UNARY_OPERATOR;
		}
	}
}

inline void Node::_StateInfo::addBinaryOperatorParenths( const std::vector< std::string > & ops )
{
	for( unsigned int i=0 ; i<state.size() ; i++ ) if( state[i].type==NodeType::BINARY_OPERATOR && IsOperator( state[i].name , ops ) )
	{
		unsigned int begin , end;

		if( state[i-1].type==NodeType::VARIABLE || state[i-1].type==NodeType::CONSTANT ) begin = i-1;
		else if( state[i-1].type==NodeType::R_PARENTH ) begin = openingParenth( i-1 );
		else MK_THROW( "Expected variable or right parenth" );

		if( state[i+1].type==NodeType::VARIABLE || state[i+1].type==NodeType::CONSTANT ) end = i+1;
		else if( state[i+1].type==NodeType::L_PARENTH ) end = closingParenth( i+1 );
		else MK_THROW( "Expected variable or left parenth" );

		if( begin==0 || state[begin-1].type!=NodeType::L_PARENTH || closingParenth( begin-1 )!=end+1 )
		{
			std::vector< State > _state;
			for( unsigned int j=0 ; j<begin ; j++ ) _state.push_back( state[j] );
			_state.emplace_back( NodeType::L_PARENTH , "(" );
			for( unsigned int j=begin ; j<=end ; j++ ) _state.push_back( state[j] );
			_state.emplace_back( NodeType::R_PARENTH , ")" );
			for( unsigned int j=end+1 ; j<state.size() ; j++ ) _state.push_back( state[j] );
			state = _state;
			addBinaryOperatorParenths( ops );
		}
	}
}

inline bool Node::_StateInfo::IsComma( char c ){ return c==','; }

inline bool Node::_StateInfo::IsParenth( char c ){ return c=='(' || c==')'; }

inline bool Node::_StateInfo::IsOperator( char c ){ return IsOperator( std::string( 1 , c ) ); }
inline bool Node::_StateInfo::IsOperator( std::string str )
{
	return IsUnaryOperator( str ) || IsBinaryOperator( str );
}

inline bool Node::_StateInfo::IsOperator( char c , const std::vector< std::string > & ops ){ return IsOperator( std::string( 1 , c ) , ops ); }
inline bool Node::_StateInfo::IsOperator( std::string str , const std::vector< std::string > & ops )
{
	for( unsigned int i=0 ; i<ops.size() ; i++ ) if( str==ops[i] ) return true;
	return false;
}

inline bool Node::_StateInfo::IsUnaryOperator( char c ){ return IsUnaryOperator( std::string( 1 , c ) ); }
inline bool Node::_StateInfo::IsUnaryOperator( std::string str )
{
	for( unsigned int i=0 ; i<UnaryOperators.size() ; i++ ) if( IsOperator( str , UnaryOperators[i] ) ) return true;
	return false;
}

inline bool Node::_StateInfo::IsBinaryOperator( char c ){ return IsBinaryOperator( std::string( 1 , c ) ); }
inline bool Node::_StateInfo::IsBinaryOperator( std::string str )
{
	for( unsigned int i=0 ; i<BinaryOperators.size() ; i++ ) if( IsOperator( str , BinaryOperators[i] ) ) return true;
	return false;
}

inline unsigned int Node::_StateInfo::OpeningParenth( const std::vector< Node::_StateInfo::State > & state , unsigned int idx )
{
	if( state[idx].type!=NodeType::R_PARENTH ) MK_THROW( "Expected right parenth: " , idx );
	int pCount = 1;
	for( unsigned int i=idx ; i>0 ; i-- )
	{
		if     ( state[i-1].type==NodeType::L_PARENTH ) pCount--;
		else if( state[i-1].type==NodeType::R_PARENTH ) pCount++;
		if( pCount==0 ) return i-1;
	}
	MK_THROW( "Could not find opening parenth: " , idx );
	return static_cast< unsigned int >(-1);
}

inline unsigned int Node::_StateInfo::ClosingParenth( const std::vector< Node::_StateInfo::State > & state , unsigned int idx )
{
	if( state[idx].type!=NodeType::L_PARENTH ) MK_THROW( "Expected left parenth: " , idx );
	int pCount = 1;
	for( unsigned int i=idx+1 ; i<state.size() ; i++ )
	{
		if     ( state[i].type==NodeType::L_PARENTH ) pCount++;
		else if( state[i].type==NodeType::R_PARENTH ) pCount--;
		if( pCount==0 ) return i;
	}
	MK_THROW( "Could not find closing parenth: " , idx );
	return static_cast< unsigned int >(-1);
}
inline void Node::_StateInfo::Print( const std::vector< Node::_StateInfo::State > & state , bool showType , bool indent )
{
	unsigned int iCount = 0;
	for( unsigned int i=0 ; i<state.size() ; i++ )
	{
		if( state[i].type==Node::_StateInfo::NodeType::R_PARENTH ) iCount--;
		if( indent ) for( unsigned int j=0 ; j<iCount ; j++ ) std::cout << "..";

		if( showType ) std::cout << NodeTypeNames[ static_cast< int >( state[i].type ) ] << " : " << state[i].name << std::endl;
		else           std::cout << state[i].name;
		if( state[i].type==Node::_StateInfo::NodeType::L_PARENTH ) iCount++;
	}
	if( !showType ) std::cout << std::endl;
}

inline std::vector< Node::_StateInfo > Node::_StateInfo::splitOnCommas( void ) const
{
	unsigned int idx = 0;
	std::vector< unsigned int > commaIndices;
	while( idx<state.size() )
	{
		if( state[idx].type==Node::_StateInfo::NodeType::L_PARENTH ) idx = ClosingParenth( state , idx );
		else if( state[idx].type==_StateInfo::NodeType::COMMA ) commaIndices.push_back( idx );
		idx++;
	}

	std::vector< Node::_StateInfo > stateInfos( commaIndices.size()+1 );
	unsigned int start = 0;
	for( unsigned int i=0 ; i<=commaIndices.size() ; i++ )
	{
		unsigned int end = i<commaIndices.size() ? commaIndices[i] : static_cast< unsigned int >( state.size() );
		stateInfos[i].state.resize( end - start );
		for( unsigned int j=start ; j<end ; j++ ) stateInfos[i].state[j-start] = state[j];
		start = end+1;
	}

	return stateInfos;
}


//////////
// Node //
//////////

inline Node Node::Constant( double v )
{
	Node n;
	n._type = _NodeType::CONSTANT;
	n._value = v;
	return n;
}

inline Node Node::Variable( unsigned int i ){ return Node(i); }


inline Node::Node( void ) : _type( _NodeType::ZERO ) , _value(0) {}
inline Node::Node( unsigned int idx ) : _type( _NodeType::VARIABLE ) , _variableIndex( idx ) {}
inline Node::Node( std::string eqn , const std::vector< std::string > & vars , bool compress )
{
	_StateInfo stateInfo( eqn , vars );
	*this = _Parse( stateInfo , vars );
	if( compress ) this->compress();
}


inline Node Node::_Constant( double c ){ Node n ; n._type = _NodeType::CONSTANT ; n._value = c ; return n; }

inline Node Node::_DVariable( unsigned int vIndex , unsigned int dIndex )
{
	Node node;
	if( vIndex==dIndex )
	{
		node._type = _NodeType::CONSTANT;
		node._value = 1.;
	}
	else node._type = _NodeType::ZERO;
	return node;
}

inline double Node::_Evaluate( _NodeType type , const std::vector< double > & values )
{
	switch( type )
	{
	case _NodeType::ADDITION:
	{
		double val = 0.;
		for( unsigned int i=0 ; i<values.size() ; i++ ) val += values[i];
		return val;
	}
	case _NodeType::MULTIPLICATION:
	{
		double val = 1.;
		for( unsigned int i=0 ; i<values.size() ; i++ ) val *= values[i];
		return val;
	}
	case _NodeType::POWER:              return pow( values[0] , values[1] );
	case _NodeType::EXPONENTIAL:        return exp ( values[0] );
	case _NodeType::NATURAL_LOGARITHM:  return log ( values[0] );
	case _NodeType::COSINE:             return cos ( values[0] );
	case _NodeType::SINE:               return sin ( values[0] );
	case _NodeType::TANGENT:            return tan ( values[0] );
	case _NodeType::HYPERBOLIC_COSINE:  return cosh( values[0] );
	case _NodeType::HYPERBOLIC_SINE:    return sinh( values[0] );
	case _NodeType::HYPERBOLIC_TANGENT: return tanh( values[0] );
	case _NodeType::SIGN:               return values[0]<0 ? -1 : 1.;
	case _NodeType::ABSOLUTE_VALUE:     return fabs( values[0] );
	default: MK_THROW( "Node type is not a unary operator or function: " , _NodeTypeNames[ static_cast< unsigned int >(type) ] );
	}
	return 0;
}

inline Node Node::_Function( _NodeType type , const Node & n )
{
	Node node;
	node._children.push_back( n );
	node._type = type;
	return node;
}

inline Node Node::_Function( _NodeType type , const Node & node1 , const Node & node2 )
{
	Node node;
	node._children.push_back( node1 );
	node._children.push_back( node2 );
	node._type = type;
	return node;
}

inline Node Node::_Function( _NodeType type , const std::vector< Node > & nodes )
{
	Node node;
	node._children = nodes;
	node._type = type;
	return node;
}

inline Node Node::_DFunction( _NodeType type , const Node & n , const Node & d )
{
	switch( type )
	{
	case _NodeType::EXPONENTIAL:        return   d * Exp( n );
	case _NodeType::NATURAL_LOGARITHM:  return   d * Pow( n , -1. );
	case _NodeType::COSINE:             return - d * Sin( n );
	case _NodeType::SINE:               return   d * Cos( n );
	case _NodeType::TANGENT:            return   d * ( 1. + Pow( Tan(n) , 2 ) );
	case _NodeType::HYPERBOLIC_COSINE:  return   d * Sinh( n );
	case _NodeType::HYPERBOLIC_SINE:    return   d * Cosh( n );
	case _NodeType::HYPERBOLIC_TANGENT: return   d / Pow( Cosh(n) , 2 );
	case _NodeType::SIGN:               return   Node();
	case _NodeType::ABSOLUTE_VALUE:     return   Sgn( n ) * d;
	default: MK_THROW( "Node type is not a unary operator or function: " , _NodeTypeNames[ static_cast< unsigned int >(type) ] );
	}

	return Node();
}

inline Node Node::_DFunction( _NodeType type , const Node & node1 , const Node & dNode1 , const Node & node2 , const Node & dNode2 )
{
	switch( type )
	{
	case _NodeType::ADDITION:       return dNode1 + dNode2;
	case _NodeType::MULTIPLICATION: return ( dNode1 * node2 ) + ( node1 * dNode2 );
	case _NodeType::POWER:
		if( node2._type==_NodeType::CONSTANT ) return node2._value * dNode1 * Pow( node1 ,node2._value-1 );
		else                                   return Pow( node1 , node2 ) * ( Log( node1 ) * dNode2 + node2 / node1 * dNode1 );
	default: MK_THROW( "Node type is not a binary operator: " , _NodeTypeNames[ static_cast< unsigned int >(type) ] );
	}

	return Node();
}

inline Node Node::_DFunction( _NodeType type , const std::vector< Node > & nodes , const std::vector< Node > & dNodes )
{
	switch( type )
	{
	case _NodeType::ADDITION: return _Function( type , dNodes );
	case _NodeType::MULTIPLICATION:
	{
		std::vector< Node > summands( nodes.size( ) );
		{
			std::vector< Node > products( nodes.size() );
			for( unsigned int i=0 ; i<nodes.size() ; i++ )
			{
				for( unsigned int j=0 ; j<nodes.size() ; j++ )
					if( i==j ) products[j] = dNodes[j];
					else       products[j] =  nodes[j];
				summands[i] = _Function( _NodeType::MULTIPLICATION , products );
			}
		}
		return _Function( _NodeType::ADDITION , summands );
	}
	default: MK_THROW( "Node type is not an n-ary operator: " , _NodeTypeNames[ static_cast< unsigned int >(type) ] );
	}

	return Node();
}

inline Node Node::_Parse( _StateInfo stateInfo , const std::vector< std::string > & vars )
{
	auto GetConstantNode = []( std::string str )
	{
		double v;
		try
		{
			if( str=="Pi" ) v = M_PI;
			else v = std::stod( str );
		}
		catch( ... ){ MK_THROW( "Failed to convert constant to double: " , str ); }
		return _Constant( v );
	};
	auto GetVariableNode = [&]( std::string str )
	{
		unsigned int idx = static_cast< unsigned int >(-1);
		for( unsigned i=0 ; i<vars.size() ; i++ ) if( vars[i]==str ) idx = i;
		if( idx==static_cast< unsigned int >(-1) ) MK_THROW( "Could not find variable: " , str );
		return Node( idx );
	};

	if( stateInfo.state.size()==0 ) return Node();
	else if( stateInfo.state.size()==1 )
	{
		if     ( stateInfo.state[0].type==_StateInfo::NodeType::CONSTANT ) return GetConstantNode( stateInfo.state[0].name );
		else if( stateInfo.state[0].type==_StateInfo::NodeType::VARIABLE ) return GetVariableNode( stateInfo.state[0].name );
		else MK_THROW( "Expected constant or variable" );
	}
	else
	{
		// Check that there are enclosing parenths
		if( stateInfo.state.front().type!=_StateInfo::NodeType::L_PARENTH || stateInfo.state.back().type!=_StateInfo::NodeType::R_PARENTH ) MK_THROW( "Expected parenthesized state" );

		// Strip off the enclosing parenths
		stateInfo = stateInfo.sub( 1 , static_cast< unsigned int >( stateInfo.state.size()-1 ) );

		if( stateInfo.state.size()==1 ) return _Parse( stateInfo , vars );
		else
		{
			unsigned int idx = -1;

			// Find the function's symbol
			{
				unsigned int pCount = 0;
				for( unsigned int i=0 ; i<stateInfo.state.size() ; i++ )
				{
					if     ( stateInfo.state[i].type==_StateInfo::NodeType::L_PARENTH ) pCount++;
					else if( stateInfo.state[i].type==_StateInfo::NodeType::R_PARENTH ) pCount--;
					else if( pCount==0 )
					{
						if( stateInfo.state[i].type==_StateInfo::NodeType::FUNCTION || stateInfo.state[i].type==_StateInfo::NodeType::UNARY_OPERATOR || stateInfo.state[i].type==_StateInfo::NodeType::BINARY_OPERATOR )
						{
							if( idx!=-1 ) MK_THROW( "Expected a single operator or function" );
							else idx = i;
						}
					}
				}
			}
			if( idx==-1 )
			{
				std::cout << "State info:" << std::endl;
				_StateInfo::Print( stateInfo.state , true , true );
				MK_THROW( "Could not find operator or function" );
			}

			if( stateInfo.state[idx].type==_StateInfo::NodeType::UNARY_OPERATOR )
			{
				if( idx+1>stateInfo.state.size() ) MK_THROW( "Expected argument after operator" );
				unsigned int begin = idx+1;

				Node childNode;
				if     ( stateInfo.state[begin].type==_StateInfo::NodeType::CONSTANT  ) childNode = GetConstantNode( stateInfo.state[begin].name );
				else if( stateInfo.state[begin].type==_StateInfo::NodeType::VARIABLE  ) childNode = GetVariableNode( stateInfo.state[begin].name );
				else if( stateInfo.state[begin].type==_StateInfo::NodeType::L_PARENTH ) childNode = _Parse( stateInfo.sub( begin , stateInfo.closingParenth( begin )+1 ) , vars );
				else MK_THROW( "Expected constant, variable, or parentheses-enclosed equation" );
				return stateInfo.state[idx]( childNode );
			}
			else if( stateInfo.state[idx].type==_StateInfo::NodeType::BINARY_OPERATOR )
			{
				if( idx==0 ) MK_THROW( "Expected content before binary operator" );
				if( idx+1>stateInfo.state.size() ) MK_THROW( "Expected content after binary operator" );

				Node childNode1 , childNode2;
				{
					unsigned int end = idx-1;
					if     ( stateInfo.state[end].type==_StateInfo::NodeType::CONSTANT  ) childNode1 = GetConstantNode( stateInfo.state[end].name );
					else if( stateInfo.state[end].type==_StateInfo::NodeType::VARIABLE  ) childNode1 = GetVariableNode( stateInfo.state[end].name );
					else if( stateInfo.state[end].type==_StateInfo::NodeType::R_PARENTH ) childNode1 = _Parse( stateInfo.sub( stateInfo.openingParenth( end ) , end+1 ) , vars );
					else MK_THROW( "Expected constant, variable, or parenthesized expression" );
				}
				{
					unsigned int begin = idx+1;
					if     ( stateInfo.state[begin].type==_StateInfo::NodeType::CONSTANT  ) childNode2 = GetConstantNode( stateInfo.state[begin].name );
					else if( stateInfo.state[begin].type==_StateInfo::NodeType::VARIABLE  ) childNode2 = GetVariableNode( stateInfo.state[begin].name );
					else if( stateInfo.state[begin].type==_StateInfo::NodeType::L_PARENTH ) childNode2 = _Parse( stateInfo.sub( begin , stateInfo.closingParenth( begin )+1 ) , vars );
					else MK_THROW( "Expected constant, variable, or parenthesized expression" );
				}
				return stateInfo.state[idx]( childNode1 , childNode2 );
			}
			else if( stateInfo.state[idx].type==_StateInfo::NodeType::FUNCTION )
			{
#ifdef NEW_EQUATION_PARSER
				if( stateInfo.state[idx+1].type!=_StateInfo::NodeType::L_PARENTH ) MK_THROW( "Expected leading parenth" );

				// The function argument, with parentheses stripped off
				_StateInfo _stateInfo = stateInfo.sub( idx+2 , stateInfo.closingParenth( idx+1 ) );

				std::vector< _StateInfo > subStates = _stateInfo.splitOnCommas();
				if( subStates.size()==1 )
				{
					subStates[0].addParenths();
					return stateInfo.state[idx]( _Parse( subStates[0] , vars ) );
				}
				else if( subStates.size()==2 )
				{
					subStates[0].addParenths();
					subStates[1].addParenths();
					return stateInfo.state[idx]( _Parse( subStates[0] , vars ) , _Parse( subStates[1] , vars ) );
				}
				else MK_ERROR_OUT( "Only unary and binary functions supported: " , subStates.size() );

				Node childNode1 , childNode2;
				unsigned int _idx;
				{
					unsigned int start = idx+2;
					if     ( stateInfo.state[start].type==_StateInfo::NodeType::CONSTANT  ) childNode1 = GetConstantNode( stateInfo.state[start].name ) , _idx = start+1;
					else if( stateInfo.state[start].type==_StateInfo::NodeType::VARIABLE  ) childNode1 = GetVariableNode( stateInfo.state[start].name ) , _idx = start+1;
					else if( stateInfo.state[start].type==_StateInfo::NodeType::L_PARENTH )
					{
						unsigned int end = stateInfo.closingParenth( start );
						childNode1 = _Parse( stateInfo.sub( start , end+1 ) , vars );
						_idx = end+1;
					}
					else if( stateInfo.state[start].type==_StateInfo::NodeType::FUNCTION )
					{
						unsigned int end = static_cast< unsigned int >( stateInfo.state.size()-1 );
						childNode1 = _Parse( stateInfo.sub( start , end ) , vars );
						_idx = end+1;
					}
					else MK_THROW( "Expected constant, variable, or parenthesized expression" );
				}

				if( _idx>=stateInfo.state.size() || stateInfo.state[_idx].type!=_StateInfo::NodeType::COMMA ) return stateInfo.state[idx]( childNode1 );
				else
				{
					unsigned int start = _idx+1;
					if     ( stateInfo.state[start].type==_StateInfo::NodeType::CONSTANT  ) childNode2 = GetConstantNode( stateInfo.state[start].name ) , _idx = start+1;
					else if( stateInfo.state[start].type==_StateInfo::NodeType::VARIABLE  ) childNode2 = GetVariableNode( stateInfo.state[start].name ) , _idx = start+1;
					else if( stateInfo.state[start].type==_StateInfo::NodeType::L_PARENTH )
					{
						unsigned int end = stateInfo.closingParenth( start );
						childNode2 = _Parse( stateInfo.sub( start , end+1 ) , vars );
						_idx = end+1;
					}
					else MK_THROW( "Expected constant, variable, or parenthesized expression" );
					return stateInfo.state[idx]( childNode1 , childNode2 );
				}
#else // !NEW_EQUATION_PARSER
				if( stateInfo.state[idx+1].type!=_StateInfo::NodeType::L_PARENTH ) MK_THROW( "Expected leading parenth" );

				Node childNode1 , childNode2;
				unsigned int _idx;
				{
//					unsigned int start = idx+2;
					unsigned int start = idx+1;
					if     ( stateInfo.state[start].type==_StateInfo::NodeType::CONSTANT  ) childNode1 = GetConstantNode( stateInfo.state[start].name ) , _idx = start+1;
					else if( stateInfo.state[start].type==_StateInfo::NodeType::VARIABLE  ) childNode1 = GetVariableNode( stateInfo.state[start].name ) , _idx = start+1;
					else if( stateInfo.state[start].type==_StateInfo::NodeType::L_PARENTH )
					{
						unsigned int end = stateInfo.closingParenth( start );
						childNode1 = _Parse( stateInfo.sub( start , end+1 ) , vars );
						_idx = end+1;
					}
					else MK_THROW( "Expected constant, variable, or parenthesized expression" );
				}

				if( _idx>=stateInfo.state.size() || stateInfo.state[_idx].type!=_StateInfo::NodeType::COMMA ) return stateInfo.state[idx]( childNode1 );
				else
				{
					unsigned int start = _idx+1;
					if     ( stateInfo.state[start].type==_StateInfo::NodeType::CONSTANT  ) childNode2 = GetConstantNode( stateInfo.state[start].name ) , _idx = start+1;
					else if( stateInfo.state[start].type==_StateInfo::NodeType::VARIABLE  ) childNode2 = GetVariableNode( stateInfo.state[start].name ) , _idx = start+1;
					else if( stateInfo.state[start].type==_StateInfo::NodeType::L_PARENTH )
					{
						unsigned int end = stateInfo.closingParenth( start );
						childNode2 = _Parse( stateInfo.sub( start , end+1 ) , vars );
						_idx = end+1;
					}
					else MK_THROW( "Expected constant, variable, or parenthesized expression" );
					return stateInfo.state[idx]( childNode1 , childNode2 );
				}
#endif // NEW_EQUATION_PARSER
			}
			else MK_THROW( "Expected operator or function" );
		}
	}
	return Node();
}

inline std::ostream & operator << ( std::ostream & stream , const Node & node )
{
	Node::_Insert( stream , node , []( unsigned int idx ){ return std::string( "x" ) + std::to_string( idx ); } , true );
	return stream;
}

inline std::string Node::operator()( const std::vector< std::string > & varNames ) const
{
	std::stringstream sStream;
	_Insert( sStream , *this , [&]( unsigned int idx ){ return varNames[idx]; } ,true );
	return sStream.str();
}

inline bool Node::_isNegative( void ) const
{
	if( _type==_NodeType::CONSTANT ) return _value<0;
	else if( _type==_NodeType::MULTIPLICATION ) return _children[0]._isNegative();
	else return false;
}

inline bool Node::_divide( const Node & node )
{
	if( operator==( node ) )
	{
		*this = _Constant(1.);
		return true;
	}
	else if( _type==_NodeType::MULTIPLICATION )
	{
		for( unsigned int i=0 ; i<_children.size() ; i++ ) if( _children[i]._divide( node ) ) return true;
		return false;
	}
	else if( _type==_NodeType::ADDITION )
	{
		for( unsigned int i=0 ; i<_children.size() ; i++ ) if( !_children[i]._divide( node ) ) return false;
		return true;
	}
	else if( _type==_NodeType::POWER )
	{
		if( _children[1]._type==_NodeType::CONSTANT && _children[1]._value>=1 && _children[0]._isDivisible( node ) )
		{
			_children[1]._value -= 1.;
			return true;
		}
	}
	return false;
}

inline bool Node::_isDivisible( const Node & node ) const
{
	if( operator==( node ) ) return true;
	else if( _type==_NodeType::MULTIPLICATION )
	{
		for( unsigned int i=0 ; i<_children.size() ; i++ ) if( _children[i]._isDivisible( node ) ) return true;
		return false;
	}
	else if( _type==_NodeType::ADDITION )
	{
		for( unsigned int i=0 ; i<_children.size() ; i++ ) if( !_children[i]._isDivisible( node ) ) return false;
		return true;
	}
	else if( _type==_NodeType::POWER )
	{
		if( _children[1]._type==_NodeType::CONSTANT && _children[1]._value>=1 && _children[0]._isDivisible( node ) ) return true;
	}
	return false;
}
inline void Node::_Insert( std::ostream & stream , const Node & node , const std::function< std::string ( unsigned int ) > & varNameFunctor , bool processSign )
{
	switch( node._type )
	{
	case Node::_NodeType::ZERO:
		stream << 0;
		break;
	case Node::_NodeType::CONSTANT:
		if( processSign ) stream << node._value;
		else              stream << fabs( node._value );
		break;
	case Node::_NodeType::VARIABLE:
		stream << varNameFunctor( node._variableIndex );
		break;
	default:
		if( _IsFunction( node._type ) )
		{
			stream << _ToString( node._type ) << "(";
			_Insert( stream , node._children[0] , varNameFunctor , true );
			stream << ")";
		}
		else if( _IsOperator( node._type ) )
		{
			if( processSign && node._type==_NodeType::MULTIPLICATION && node._isNegative() ) stream << "-";
			__Insert( stream , node , varNameFunctor );
		}
		else MK_THROW( "Unrecognized node type: " , _NodeTypeNames[ static_cast< unsigned int >( node._type ) ] );
	}
}

inline void Node::__Insert( std::ostream & stream , const Node & node , const std::function< std::string ( unsigned int ) > & varNameFunctor )
{
	auto NeedsParenths = [&]( const Node & node )
	{
		return 
			node._type!=Node::_NodeType::CONSTANT &&
			node._type!=Node::_NodeType::VARIABLE &&
			!( node._type==Node::_NodeType::POWER && node._children[0]._type==Node::_NodeType::VARIABLE );
	};

	auto InsertNode = [&]( const Node & node )
	{
		if( NeedsParenths( node ) )
		{
			stream << "(";
			_Insert( stream , node , varNameFunctor , false );
			stream << ")";
		}
		else _Insert( stream , node , varNameFunctor , false );
	};

	bool first = true;
	for( size_t i=0 ; i<node._children.size() ; i++ )
	{
		if( node._type==_NodeType::MULTIPLICATION && node._children[i]._type==_NodeType::POWER && node._children[i]._children[1]._type==_NodeType::CONSTANT && node._children[i]._children[1]._isNegative() )
		{
			if( first ) stream << "1/";
			else        stream <<  "/";
			InsertNode( node._children[i]._children[0] );
			if( node._children[i]._children[1]._type==_NodeType::CONSTANT && node._children[i]._children[1]._value==-1 );
			else
			{
				stream << "^";
				InsertNode( node._children[i]._children[1] );
			}
			first = false;
		}
		else
		{
			if( node._type==_NodeType::MULTIPLICATION && node._children[i]._type==_NodeType::CONSTANT && node._children[i]._value==-1 ) continue;
			else if( node._type==_NodeType::ADDITION && node._children[i]._isNegative() ) stream << "-";
			else if( !first ) stream << _ToString( node._type );
			InsertNode( node._children[i] );
			first = false;
		}
	}
}

inline bool Node::_isConstant( void ) const
{
	bool isConstant = _type!=_NodeType::VARIABLE;
	for( unsigned int i=0 ; i<_children.size() ; i++ ) isConstant &= _children[i]._isConstant();
	return isConstant;
}

inline bool Node::__preCompress( void )
{
	auto IsZero             = []( const Node & n ){ return n._type==_NodeType::ZERO; };
	auto IsOne              = []( const Node & n ){ return n._type==_NodeType::CONSTANT && n._value==1; };
	auto IsMinusOne         = []( const Node & n ){ return n._type==_NodeType::CONSTANT && n._value==-1; };
	auto IsScalar           = []( const Node & n ){ return n._type==_NodeType::ZERO || n._type==_NodeType::CONSTANT; };
	auto IsVariableExponent = []( const Node & n ){ return n._type==_NodeType::POWER && n._children[0]._type==_NodeType::VARIABLE; };
	auto VariableIndex      = []( const Node & n )
	{
		if     ( n._type==_NodeType::VARIABLE )                                           return std::optional< unsigned int >( n._variableIndex );
		else if( n._type==_NodeType::POWER && n._children[0]._type==_NodeType::VARIABLE ) return std::optional< unsigned int >( n._children[0]._variableIndex  );
		else                                                                              return std::optional< unsigned int >( std::nullopt );
	};

	// Trim off constant sub-trees
	if( _isConstant() && _type!=_NodeType::ZERO && _type!=_NodeType::CONSTANT )
	{
		_value = operator()( nullptr );
		_type = _value==0 ? _NodeType::ZERO : _NodeType::CONSTANT;
		_children.resize(0);
		return true;
	}

	// Mark zero constants as zero
	if( _type==_NodeType::CONSTANT && !_value )
	{
		_type = _NodeType::ZERO;
		return true;
	}

	// Collapse trivial summation
	if( _type==_NodeType::ADDITION && _children.size()==0 )
	{
		*this = Node();
		return true;
	}

	// Collapse trivial multiplication
	if( _type==_NodeType::MULTIPLICATION && _children.size()==0 )
	{
		*this = _Constant(1.);
		return true;
	}

	// Collapse n-ary with a single argument to the argument
	if( _IsNAryOperator( _type ) && _children.size()==1 )
	{
		Node n = _children[0];
		*this = n;
		return true;
	}

	// Collapse n-ary children into parent
	if( _IsNAryOperator( _type ) )
	{
		for( unsigned int i=0 ; i<_children.size() ; i++ ) if( _children[i]._type==_type )
		{
			std::vector< Node > children;
			children.reserve( children.size() - 1 + _children[i].size() );
			for( unsigned int j=0 ; j<i ; j++ ) children.push_back( _children[j] );
			for( unsigned int j=0 ; j<_children[i]._children.size() ; j++ ) children.push_back( _children[i]._children[j] );
			for( unsigned int j=i+1 ; j<_children.size() ; j++ ) children.push_back( _children[j] );
			*this = _Function( _type , children );
			return true;
		}
	}

	if( _type==_NodeType::ADDITION )
	{
		// Collapse scalars and move to the front
		{
			unsigned int sCount = 0 , sIdx = -1;
			double s = 0;
			for( unsigned int i=0 ; i<_children.size() ; i++ ) if( IsScalar( _children[i] ) ) sCount++ , sIdx=i , s += _children[i]._value;
			if( sCount>1 )
			{
				unsigned int idx = 0;
				for( unsigned int i=0 ; i<_children.size() ; i++ ) if( !IsScalar( _children[i] ) ) _children[idx++] = _children[i];
				if( idx==0 )
				{
					*this = _Constant( s );
					return true;
				}

				Node n = *this;
				n._children.resize( idx );
				if( idx==1 ) n = _children[0];
				if( s==0 ) *this = n;
				else       *this = s+n;
				return true;
			}
			else if( sCount==1 && sIdx!=0 )
			{
				std::swap( _children[0] , _children[sIdx] );
				return true;
			}
			else if( sCount==1 && s==0 )
			{
				_children[0] = _children.back();
				_children.pop_back();
				return true;
			}
		}
	}

	if( _type==_NodeType::MULTIPLICATION )
	{
		// If any of the products are zero
		{
			bool hasZero = false;
			for( unsigned int i=0 ; i<_children.size() ; i++ ) hasZero |= IsZero( _children[i] );
			if( hasZero )
			{
				*this = Node();
				return true;
			}
		}

		// Collapse scalars and move to the front
		{
			unsigned int sCount = 0 , sIdx = -1;
			double s = 1.;
			for( unsigned int i=0 ; i<_children.size() ; i++ ) if( IsScalar( _children[i] ) ) sCount++ , sIdx=i , s*=_children[i]._value;;
			if( sCount>1 )
			{
				unsigned int idx = 0;
				for( unsigned int i=0 ; i<_children.size() ; i++ ) if( !IsScalar( _children[i] ) ) _children[idx++] = _children[i];

				Node n = *this;
				n._children.resize( idx );
				if     ( idx==0 ) n = Node();
				else if( idx==1 ) n = _children[0];
				if( s==1 ) *this = n;
				else       *this = s*n;
				return true;
			}
			else if( sCount==1 && sIdx!=0 )
			{
				std::swap( _children[0] , _children[sIdx] );
				return true;
			}
			else if( sCount==1 && s==1 )
			{
				_children[0] = _children.back();
				_children.pop_back();
				return true;
			}
		}


		// Collapse products involving the same variable
		for( unsigned int i=0 ; i<_children.size() ; i++ ) if( auto idx1=VariableIndex( _children[i] ) )
			for( unsigned int j=i+1 ; j<_children.size() ; j++ ) if( auto idx2=VariableIndex( _children[j] ) )
				if( *idx1==*idx2 )
				{
					Node n , v( *idx1 );
					if( IsVariableExponent( _children[i] ) && IsVariableExponent( _children[j] ) ) n = v^( _children[i]._children[1] + _children[j]._children[1] );
					else if( IsVariableExponent( _children[i] ) ) n = v^( _children[i]._children[1] + 1 );
					else if( IsVariableExponent( _children[j] ) ) n = v^( _children[j]._children[1] + 1 );
					else n = v^2;
					_children[i] = n;
					_children[j] = _children.back();
					_children.pop_back();
					return true;
				}
	}

	if( _type==_NodeType::POWER )
		if( _children[1]._type==_NodeType::ZERO )
		{
			*this = _Constant(1.);
			return true;
		}
		else if( _children[1]._type==_NodeType::CONSTANT && _children[1]._value==1. )
		{
			Node n = _children[0];
			*this = n;
			return true;
		}
		else if( _children[0]._type==_NodeType::POWER )
		{
			Node n = _children[0]._children[0] ^ ( _children[0]._children[1] * _children[1] );
			*this = n;
			return true;
		}
	return false;
}

inline bool Node::_isReciprocal( void ) const { return _type==_NodeType::POWER && _children[1]._type==_NodeType::CONSTANT && _children[1]._value==-1; }

inline bool Node::_hasNumerator( const Node & node ) const
{
	if( operator==(node) ) return true;
	else if( _type==_NodeType::POWER && _children[0]==node && _children[1]._type==_NodeType::CONSTANT && _children[1]._value>=1. ) return true;
	else if( _type==_NodeType::ADDITION )
	{
		for( const auto & child : _children ) if( !child._hasNumerator( node ) ) return false;
		return true;
	}
	else if( _type==_NodeType::MULTIPLICATION ) for( const auto & child : _children ) if( child._hasNumerator( node ) ) return true;
	return false;
}

inline bool Node::_hasDenominator( const Node & node ) const
{
	if( _type==_NodeType::POWER && _children[0]==node && _children[1]._type==_NodeType::CONSTANT && _children[1]._value<=-1. ) return true;
	else if( _type==_NodeType::MULTIPLICATION ) for( const auto & child : _children ) if( child._hasDenominator( node ) ) return true;
	else if( _type==_NodeType::ADDITION )
	{
		for( const auto & child : _children ) if( !child._hasDenominator( node ) ) return false;
		return true;
	}
	return false;
}

inline bool Node::_removeNumerator( const Node & node )
{
	if( operator==(node) ){ *this = _Constant(1.) ; return true; }
	else if( _type==_NodeType::POWER && _children[0]==node && _children[1]._type==_NodeType::CONSTANT && _children[1]._value>=1. ){ _children[1]._value -= 1. ; return true; }
	else if( _type==_NodeType::ADDITION )
	{
		auto copy = *this;
		for( auto & child : copy._children ) if( !child._removeNumerator( node ) ) return false;
		*this = copy;
		return true;
	}
	else if( _type==_NodeType::MULTIPLICATION )
	{
		for( auto & child : _children )	if( child._removeNumerator( node ) ) return true;
		return false;
	}
	return false;
}

inline bool Node::_removeDenominator( const Node & node )
{
	if( _type==_NodeType::POWER && _children[0]==node && _children[1]._type==_NodeType::CONSTANT && _children[1]._value<=-1. ){ _children[1]._value += 1. ; return true; }
	else if( _type==_NodeType::ADDITION )
	{
		auto copy = *this;
		for( auto & child : copy._children ) if( !child._removeDenominator( node ) ) return false;
		*this = copy;
		return true;
	}
	else if( _type==_NodeType::MULTIPLICATION )
	{
		for( auto & child : _children ) if( child._removeDenominator( node ) ) return true;
		return false;
	}
	return false;
}

inline bool Node::__postCompress( void )
{
	if( _type==_NodeType::MULTIPLICATION )
	{
		for( unsigned int i=0 ; i<_children.size() ; i++ )
		{
			// Check if a particular appears more than once in a product. If it does, compactify it into a power
			Node f , n = _children[i]._type==_NodeType::POWER ? _children[i]._children[0] : _children[i];
			unsigned int fCount = 0;
			for( unsigned int j=0 ; j<_children.size() ; j++ )
				if( _children[j]._type==_NodeType::POWER && _children[j]._children[0]==n ) f += _children[j]._children[1] , fCount++;
				else if( _children[j]==n ) f += 1 , fCount++;

			if( fCount>1 )
			{
				unsigned int idx=0;
				for( unsigned int j=0 ; j<_children.size() ; j++ )
					if( _children[j]._type==_NodeType::POWER && _children[j]._children[0]==n ) _children[j] = _Constant(1.) , idx = j;
					else if( _children[j]==n ) _children[j] = _Constant(1.) , idx = j;
				_children[idx] = _Function( _NodeType::POWER , n , f );
				_preCompress();
				_sort();
				return true;
			}
		}
	}
	if( _type==_NodeType::ADDITION )
	{
		for( unsigned int idx=0 ; auto child : _children )
		{
			// If all the summands have "child" as a numerator
			if( _hasNumerator( child ) )
			{
				// Remove "child" from the summands
				if( !_removeNumerator( child ) ) MK_ERROR_OUT( "Failed to remove child" );

				// Create a node which is the product of "child" with the current node with "child" removed
				Node n = (*this) * child;
				*this = n;

				_preCompress();
				_sort();
				return true;
			}
			else if( child._type==_NodeType::MULTIPLICATION )
			{
				for( auto & gChild : child._children )
				{
					// If all the summands have "gChild" as a numerator
					if( _hasNumerator( gChild ) )
					{
						if( !_removeNumerator( gChild ) ) MK_ERROR_OUT( "Failed to remove grand-child" );
						Node n = (*this) * gChild;
						*this = n;

						_preCompress();
						_sort();
						return true;
					}
				}
			}

			if( child._isReciprocal() && _hasDenominator( child._children[0] ) )
			{
				_removeDenominator( child._children[0] );
				Node n = *this;
				*this = n / child._children[0];
				_preCompress();
				_sort();
				return true;
			}
			else if( child._type==_NodeType::MULTIPLICATION )
			{
				for( auto & gChild : child._children )
				{
					if( gChild._isReciprocal() && _hasDenominator( gChild._children[0] ) )
					{
						_removeDenominator( gChild._children[0] );
						Node n = *this;
						*this = n / gChild._children[0];
						_preCompress();
						_sort();
						return true;
					}
				}
			}
			idx++;
		}
	}
	return false;
}

inline bool Node::_preCompress( void )
{
	bool compressed = false;
	for( unsigned int i=0 ; i<_children.size() ; i++ ) compressed |= _children[i]._preCompress();
	while( __preCompress() ) compressed = true;

	return compressed;
}

inline bool Node::_postCompress( void )
{
	bool compressed = false;
	for( unsigned int i=0 ; i<_children.size() ; i++ ) compressed |= _children[i]._postCompress();
	while( __postCompress() ) compressed = true;

	return compressed;
}

inline void Node::_sort( void )
{
	for( unsigned int i=0 ; i<_children.size() ; i++ ) _children[i]._sort();
	if( _IsSymmetricOperator( _type ) ) std::sort( _children.begin() , _children.end() );
}

inline void Node::_sanityCheck( void ) const
{
	for( const auto & child : _children ) child._sanityCheck();
	if     ( _IsFunction      ( _type ) && _children.size()!=1 ) MK_THROW( "Expected one child: "        , _children.size() );
	else if( _IsBinaryOperator( _type ) && _children.size()!=2 ) MK_THROW( "Expected two children "      , _children.size() );
	else if( _IsNAryOperator  ( _type ) && _children.size()==0 ) MK_THROW( "Expected non-zero children " , _children.size() );
}

inline void Node::compress( void )
{
	while( true )
	{
		while( _preCompress() );
		_sort();
		bool done = true;
		while( _postCompress() ) done = false;
		if( done ) break;
	}
}


inline void Node::replaceVariables( const Node * nodes )
{
	if( _type==Node::_NodeType::VARIABLE ) *this = nodes[ _variableIndex ];
	else for( unsigned int i=0 ; i<_children.size() ; i++ ) _children[i].replaceVariables( nodes );
}

template< unsigned int Dim >
double Node::operator()( Point< double , Dim > p ) const { return operator()( &p[0] ); }

inline double Node::operator()( const double * values ) const
{
	std::vector< double > _values;
	switch( _type )
	{
	case Node::_NodeType::ZERO:     return 0;
	case Node::_NodeType::CONSTANT: return _value;
	case Node::_NodeType::VARIABLE: return values[ _variableIndex ];
	default:
		if( _IsFunction( _type ) )
		{
			_values.resize(1);
			_values[0] = _children[0]( values );
			return _Evaluate( _type , _values );
		}
		else if( _IsBinaryOperator( _type ) )
		{
			_values.resize(2);
			_values[0] = _children[0]( values );
			_values[1] = _children[1]( values );
			return _Evaluate( _type , _values );
		}
		else if( _IsNAryOperator( _type ) )
		{
			_values.resize( _children.size() );
			for( unsigned int i=0 ; i<_children.size() ; i++ ) _values[i] = _children[i]( values );
			return _Evaluate( _type , _values );
		}
		else
		{
			MK_THROW( "Urecognized node type: " , _NodeTypeNames[ static_cast< unsigned int >( _type ) ] );
			return 0;
		}
	}
}

inline Node Node::d( unsigned int dIndex ) const
{
	Node node;
	if     ( _type==_NodeType::ZERO || _type==_NodeType::CONSTANT ) ;
	else if( _type==_NodeType::VARIABLE ) node = _DVariable( _variableIndex , dIndex );
	else if( _IsFunction( _type ) ) node = _DFunction( _type , _children[0] , _children[0].d( dIndex ) );
	else if( _IsBinaryOperator( _type ) ) node = _DFunction( _type , _children[0] , _children[0].d(dIndex) , _children[1] , _children[1].d(dIndex) );
	else if( _IsNAryOperator( _type ) )
	{
		std::vector< Node > dChildren( _children.size() );
		for( unsigned int i=0 ; i<_children.size() ; i++ ) dChildren[i] = _children[i].d( dIndex );
		node = _DFunction( _type , _children , dChildren );
	}
	else MK_THROW( "Unrecognized type: " , _NodeTypeNames[ static_cast< unsigned int >( _type ) ] );
	node.compress();
	return node;
}

inline unsigned int Node::size( void ) const
{
	unsigned int sz = 1;
	for( unsigned int i=0 ; i<_children.size() ; i++ ) sz += _children[i].size();
	return sz;
}

inline unsigned int Node::variables( void ) const
{
	std::set< unsigned int > variableIndexSet;
	_addVariableIndices( variableIndexSet );
	return static_cast< unsigned int >( variableIndexSet.size() );
}

inline void Node::_addVariableIndices( std::set< unsigned int > & variableIndexSet ) const
{
	if( _type==_NodeType::VARIABLE ) variableIndexSet.insert( _variableIndex );
	for( unsigned int i=0 ; i<_children.size() ; i++ ) _children[i]._addVariableIndices( variableIndexSet );
}

inline void Node::Trace( const EquationParser::Node & node )
{
	std::cout << EquationParser::Node::_NodeTypeNames[ static_cast< unsigned int >( node._type ) ] << " : " << node << std::endl;
	for( unsigned int i=0 ; i<node._children.size() ; i++ )
	{
		std::cout << "\t" << i << "] " << EquationParser::Node::_NodeTypeNames[ (int)node._children[i]._type ] << " : " << node._children[i] << std::endl;
		if( node._children[i]._type==EquationParser::Node::_NodeType::POWER ) std::cout << node._children[i]._children[1]._value << std::endl;
	}
	if( node._children.size()==0 );
	else if( node._children.size()==1 ) Trace( node._children[0] );
	else
	{
		unsigned int idx;
		std::cout << "Enter index: ";
		std::cin >> idx;
		Trace( node._children[idx] );
	}
}

inline unsigned int Node::_maxVarIndex( void ) const
{
	if( _type==_NodeType::VARIABLE ) return _variableIndex;
	else
	{
		unsigned int max = 0;
		for( unsigned int i=0 ; i<_children.size() ; i++ ) max = std::max< unsigned int >( max , _children[i]._maxVarIndex() );
		return max;
	}
}

inline double Node::SanityCheckCompression( const Node & node , unsigned int count , double radius )
{
	Node _node = node;
	_node.compress();

	double e2=0 , l2=0;
	unsigned int dim = node._maxVarIndex() + 1;
	std::vector< double > values( dim );
	for( unsigned int i=0 ; i<count ; i++ )
	{
		while( true )
		{
			for( unsigned int i=0 ; i<dim ; i++ ) values[i] = Random< double >() * 2. - 1.;
			double sum = 0;
			for( unsigned int i=0 ; i<dim ; i++ ) sum += values[i] * values[i];
			if( sum<1 ) break;
		}
		double v = node( &values[0] ) , _v = _node( &values[0] );
		e2 += ( v - _v ) * ( v - _v );
		l2 += v*v + _v*_v;
	}
	return sqrt( e2 / l2 );
}


inline void Node::SanityCheckParsing( const Node & node )
{
	node._sanityCheck();
	std::function< void ( const Node & , const Node & , unsigned int ) > Compare = [&Compare]( const Node & node1 , const Node & node2 , unsigned int offset )
	{
		if( node1!=node2 )
		{
			for( unsigned int i=0 ; i<offset ; i++ ) std::cout << "." ; std::cout << node1 << std::endl;
			for( unsigned int i=0 ; i<offset ; i++ ) std::cout << "." ; std::cout << node2 << std::endl;

			if( node1._type!=node2._type )
				MK_THROW( "Types differ: " , _NodeTypeNames[ static_cast< unsigned int >( node1._type ) ] , " != " , _NodeTypeNames[ static_cast< unsigned int >( node2._type ) ] );
			else if( node1._children.size()!=node2._children.size() )
				MK_THROW( "Number of children differ: " , node1._children.size() , " != " , node2._children.size() );
			else
			{
				if( !node1._children.size() )
				{
					std::cout << "[0] " << node1 << std::endl;
					std::cout << "[1] " << node2 << std::endl;
					MK_ERROR_OUT( "..." );
				}
				else for( unsigned int i=0 ; i<node1._children.size() ; i++ ) Compare( node1._children[i] , node2._children[i] , offset+1 );
			}
		}
	};
	std::vector< std::string > varNames( node._maxVarIndex()+1 );
	for( unsigned int i=0 ; i<varNames.size() ; i++ ) varNames[i] = std::string( "x" ) + std::to_string(i);
	Node _node = node;
	_node.compress();
	Node __node( node(varNames) , varNames );
	__node.compress();
	Compare( node , __node , 0 );
}

inline bool Node::_IsTerminal( _NodeType type )
{
	switch( type )
	{
	case _NodeType::ZERO:
	case _NodeType::CONSTANT:
	case _NodeType::VARIABLE:
		return true;
	default: return false;
	}
}

inline bool Node::_IsBinaryOperator( _NodeType type )
{
	switch( type )
	{
	case _NodeType::POWER:
		return true;
	default: return false;
	}
}

inline bool Node::_IsNAryOperator( _NodeType type )
{
	switch( type )
	{
	case _NodeType::ADDITION:
	case _NodeType::MULTIPLICATION:
		return true;
	default: return false;
	}
}

inline bool Node::_IsSymmetricOperator( _NodeType type )
{
	switch( type )
	{
	case _NodeType::ADDITION:
	case _NodeType::MULTIPLICATION:
		return true;
	default: return false;
	}
}

inline bool Node::_IsOperator( _NodeType type ){ return _IsBinaryOperator( type ) || _IsNAryOperator( type ); }

inline bool Node::_IsFunction( _NodeType type )
{
	switch( type )
	{
	case _NodeType::EXPONENTIAL:
	case _NodeType::COSINE:
	case _NodeType::SINE:
	case _NodeType::TANGENT:
	case _NodeType::HYPERBOLIC_COSINE:
	case _NodeType::HYPERBOLIC_SINE:
	case _NodeType::HYPERBOLIC_TANGENT:
	case _NodeType::SIGN:
	case _NodeType::ABSOLUTE_VALUE:
		return true;
	default: return false;
	}
}

inline std::string Node::_ToString( _NodeType type )
{
	switch( type )
	{
	case _NodeType::ADDITION:           return "+";
	case _NodeType::MULTIPLICATION:     return "*";
	case _NodeType::POWER:              return "^";
	case _NodeType::EXPONENTIAL:        return "exp";
	case _NodeType::COSINE:             return "cos";
	case _NodeType::SINE:               return "sin";
	case _NodeType::TANGENT:            return "tan";
	case _NodeType::HYPERBOLIC_COSINE:  return "cosh";
	case _NodeType::HYPERBOLIC_SINE:    return "sinh";
	case _NodeType::HYPERBOLIC_TANGENT: return "tanh";
	case _NodeType::SIGN:               return "sgn";
	case _NodeType::ABSOLUTE_VALUE:     return "abs";
	default:
		MK_THROW( "Unrecognized node type: " , _NodeTypeNames[ static_cast< unsigned int >( type ) ] );
		return "";
	}
}

inline bool Node::operator < ( const Node & n ) const
{
	if     ( _type!=n._type ) return static_cast< int >( _type ) < static_cast< int >( n._type );
	else if( _type==_NodeType::ZERO ) return false;
	else if( _type==_NodeType::CONSTANT ) return _value<n._value;
	else if( _type==_NodeType::VARIABLE ) return _variableIndex<n._variableIndex;
	else if( _type==_NodeType::POWER    )
	{
		if( n._children[1]<_children[1] || _children[1]<n._children[1] ) return n._children[1]<_children[1];
		else return _children[0]<n._children[0];
	}
	else if( _type==_NodeType::ADDITION || _type==_NodeType::MULTIPLICATION )
	{
		if( _children.size()!=n._children.size() ) return _children.size()<n._children.size();
		else for( unsigned int i=0 ; i<_children.size() ; i++ ) if( _children[i]<n._children[i] || n._children[i]<_children[i] ) return _children[i]<n._children[i];
	}
	else if( _children[0]<n._children[0] || n._children[0]<_children[0] ) return _children[0]<n._children[0];
	return false;
}

inline bool Node::operator == ( const Node & n ) const
{
	if( _type!=n._type ) return false;
	else
	{
		switch( _type )
		{
		case _NodeType::ZERO: return true;
		case _NodeType::CONSTANT: return _value==n._value;
		case _NodeType::VARIABLE: return _variableIndex==n._variableIndex;
		case _NodeType::ADDITION:
		case _NodeType::MULTIPLICATION:
			if( _children.size()!=n._children.size() ) return false;
			for( unsigned int i=0 ; i<_children.size() ; i++ ) if( !(_children[i]==n._children[i]) ) return false;
			return true;
		case _NodeType::POWER: return _children[0]==n._children[0] && _children[1]==n._children[1];
		case _NodeType::EXPONENTIAL:
		case _NodeType::NATURAL_LOGARITHM:
		case _NodeType::COSINE:
		case _NodeType::SINE:
		case _NodeType::TANGENT:
		case _NodeType::HYPERBOLIC_COSINE:
		case _NodeType::HYPERBOLIC_SINE:
		case _NodeType::HYPERBOLIC_TANGENT:
		case _NodeType::SIGN:
		case _NodeType::ABSOLUTE_VALUE:
			return _children[0]==n._children[0];
		default: MK_THROW( "Unrecognized type: " , _NodeTypeNames[ static_cast< unsigned int >( _type ) ] );
		}
	}
	return false;
}

inline bool Node::operator != ( const Node & n ) const { return ! ( operator==( n ) ); }


inline void Node::_SanityCheck( const Node & node1 , const Node & node2 , unsigned int tests , double eps )
{
	unsigned int vars = std::max< unsigned int >( node1.variables() , node2.variables() );
	std::vector< double > values( vars );
	for( unsigned int i=0 ; i<tests ; i++ )
	{
		for( unsigned int j=0 ; j<vars ; j++ ) values[j] = Random< double >();
		if( fabs( node1( &values[0] ) - node2( &values[0] ) )>eps ) MK_THROW( "Values don't match: " , node1( &values[0] ) , " != " , node2( &values[0] ) );
	}
}

/////////////////////////////
// Operators and Functions //
/////////////////////////////
inline Node & operator += ( Node & n , const Node & _n ){ return n = n + _n; }
inline Node & operator -= ( Node & n , const Node & _n ){ return n = n - _n; }
inline Node & operator += ( Node & n , double s ){ return n = n + s; }
inline Node & operator -= ( Node & n , double s ){ return n = n - s; }
inline Node & operator *= ( Node & n , const Node & _n ){ return n = n * _n; }
inline Node & operator /= ( Node & n , const Node & _n ){ return n = n / _n; }
inline Node & operator *= ( Node & n , double s ){ return n = n * s; }
inline Node & operator /= ( Node & n , double s ){ return n = n / s; }
inline Node operator - ( const Node & node ){ return node * -1.; }
inline Node operator + ( const Node & node1 , const Node & node2 ){ return Node::_Function( Node::_NodeType::ADDITION , node1 , node2 ); }
inline Node operator + ( const Node & node , double s ){ return Node::_Function( Node::_NodeType::ADDITION , node , Node::_Constant( s ) ); }
inline Node operator + ( double s , const Node & node ){ return Node::_Function( Node::_NodeType::ADDITION , Node::_Constant( s ) , node ); }
inline Node operator - ( const Node & node1 , const Node & node2 ){ return node1 + ( -node2 ); }
inline Node operator - ( const Node & node , double s ){ return node + (-s); }
inline Node operator - ( double s , const Node & node ){ return s + (-node); }
inline Node operator * ( const Node & node1 , const Node & node2 ){ return Node::_Function( Node::_NodeType::MULTIPLICATION , node1 , node2 ); }
inline Node operator * ( const Node & node , double s ){ return Node::_Function( Node::_NodeType::MULTIPLICATION , node , Node::_Constant( s ) ); }
inline Node operator * ( double s , const Node & node ){ return Node::_Function( Node::_NodeType::MULTIPLICATION , Node::_Constant( s ) , node ); }
inline Node operator / ( const Node & node1 , const Node & node2 ){ return node1 * Pow( node2 , -1. ); }
inline Node operator / ( const Node & node , double s ){ return node * (1./s); }
inline Node operator / ( double s , const Node & node ){ return s * Pow( node , -1. ); }
inline Node operator ^ ( const Node & node1 , const Node & node2 ){ return Node::_Function( Node::_NodeType::POWER , node1 , node2 ); }
inline Node operator ^ ( const Node & node , double s ){ return Node::_Function( Node::_NodeType::POWER , node , Node::_Constant( s ) ); }
inline Node operator ^ ( double s , const Node & node ){ return Node::_Function( Node::_NodeType::POWER , Node::_Constant( s ) , node ); }

inline Node Pow( const Node & n1 , const Node & n2 ){ return Node::_Function( Node::_NodeType::POWER , n1 , n2 ); }
inline Node Pow( const Node & n , double s ){ return Node::_Function( Node::_NodeType::POWER , n , Node::_Constant( s ) ); }
inline Node Pow( double s , const Node & n ){ return Node::_Function( Node::_NodeType::POWER , Node::_Constant( s ) , n ); }
inline Node Exp( const Node & n ){ return Node::_Function( Node::_NodeType::EXPONENTIAL , n ); }
inline Node Log( const Node & n ){ return Node::_Function( Node::_NodeType::NATURAL_LOGARITHM , n ); }
inline Node Cos( const Node & n ){ return Node::_Function( Node::_NodeType::COSINE , n ); }
inline Node Sin( const Node & n ){ return Node::_Function( Node::_NodeType::SINE , n ); }
inline Node Tan( const Node & n ){ return Node::_Function( Node::_NodeType::TANGENT , n ); }
inline Node Cosh( const Node & n ){ return Node::_Function( Node::_NodeType::HYPERBOLIC_COSINE , n ); }
inline Node Sinh( const Node & n ){ return Node::_Function( Node::_NodeType::HYPERBOLIC_SINE , n ); }
inline Node Tanh( const Node & n ){ return Node::_Function( Node::_NodeType::HYPERBOLIC_TANGENT , n ); }
inline Node Sqrt( const Node & n ){ return n ^ 0.5; }
inline Node Sgn( const Node & n ){ return Node::_Function( Node::_NodeType::SIGN , n ); }
inline Node Abs( const Node & n ){ return Node::_Function( Node::_NodeType::ABSOLUTE_VALUE , n ); }
inline Node ReLU( const Node & n ){ return ( n + Abs(n) ) / 2.; }
inline Node Min( const Node & n1 , const Node & n2 ){ return ( n1 + n2 ) / 2 - Abs( n1 - n2 ) / 2; }
inline Node Max( const Node & n1 , const Node & n2 ){ return ( n1 + n2 ) / 2 + Abs( n1 - n2 ) / 2; }

#if 1 // NEW_CODE
void Node::Analyze( const Node & node )
{
	std::function< void ( const Node & node ) > addChildren;
	std::vector< Node > nodes;

	addChildren = [&]( const Node & node ) -> void 
		{
			for( unsigned int i=0 ; i<node._children.size() ; i++ ) addChildren( node._children[i] );
			nodes.push_back( node );
		};

	addChildren( node );
	std::sort( nodes.begin() , nodes.end() );

	unsigned int count = 0;
	for( unsigned int i=1 ; i<nodes.size() ; i++ ) if( nodes[i]!=nodes[i-1] ) count++;
	std::cout << "Nodes: " << nodes.size() << " -> " << count << std::endl;
}
#endif // NEW_CODE

#ifdef EVALUATION_NODE
////////////////////
// EvaluationNode //
////////////////////
unsigned int EvaluationNode::size( void ) const { return static_cast< unsigned int >( _nodes.size() ); }
template< unsigned int Dim > double EvaluationNode::operator()( Point< double , Dim > p ) const { return operator()( &p[0] ); }
#if 0
struct CompressedNode
{
	CompressedNode( const Node & node );

	// Evaluates the the equation-tree at a prescribed set of values
	double operator()( const double * values ) const;
protected:
	struct _CompressedNode
	{
		Node * node;
		std::vector< _CompressedNode > children;
	};
	_CompressedNode _node;
	std::vector< Node > _nodes;
};
void Node::Analyze( const Node & node )
{
	std::function< void ( const Node & node ) > addChildren;
	std::vector< Node > nodes;

	addChildren = [&]( const Node & node ) -> void 
		{
			for( unsigned int i=0 ; i<node._children.size() ; i++ ) addChildren( node._children[i] );
			nodes.push_back( node );
		};

	addChildren( node );
	std::sort( nodes.begin() , nodes.end() );

	unsigned int count = 0;
	for( unsigned int i=1 ; i<nodes.size() ; i++ ) if( nodes[i]!=nodes[i-1] ) count++;
	std::cout << "Nodes: " << nodes.size() << " -> " << count << std::endl;

}

#endif
#endif // EVALUATION_NODE
