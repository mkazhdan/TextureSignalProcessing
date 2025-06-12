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

#ifndef CMD_LINE_PARSER_INCLUDED
#define CMD_LINE_PARSER_INCLUDED
#include <cstdarg>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <string.h>
#include "Exceptions.h"


namespace MishaK
{
	/** This class represents a named argument that can be read from the command line */
	class CmdLineReadable
	{
	public:
		/** Has the argument been set */
		bool set;

		/** The argument's name */
		std::string name;

		/** Constructor with the name of the argument */
		CmdLineReadable( const std::string &name );

		/** Destructor */
		virtual ~CmdLineReadable( void );

		/** Try to set the argument from the list of command line arguments.
		*** Returns thenumber of arguments ingested.*/ 
		virtual int read( char **argv , int argc );

		/** Transforms a string into the prescribed type */
		template< typename Type >
		static Type ToType( const std::string &str );
	};

	/** This templated class represents a named argument of the prescribed type */
	template< class Type >
	class CmdLineParameter : public CmdLineReadable
	{
	public:
		/** The value the parameter has been set to */
		Type value;

		/** Constructor with the name of the argument */
		CmdLineParameter( const std::string &name );

		/** Constructor with the name of the argument and the default value */
		CmdLineParameter( const std::string &name , Type v );

		/** Try to set the argument from the list of command line arguments.
		*** Returns thenumber of arguments ingested.*/ 
		int read( char **argv , int argc );
	};

	/** This templated class represents a named argument taking a fixed number of values of the prescribed type */
	template< class Type , int Dim >
	class CmdLineParameterArray : public CmdLineReadable
	{
	public:
		/** The values the parameter has been set to */
		Type values[Dim];

		/** Constructor with the name of the argument and the default values */
		CmdLineParameterArray( const std::string &name, const Type* v=NULL );

		/** Try to set the argument from the list of command line arguments.
		*** Returns thenumber of arguments ingested.*/ 
		int read( char **argv , int argc );
	};

	/** This templated class represents a named argument taking a variable number of values of of the prescribed type */
	template< class Type >
	class CmdLineParameters : public CmdLineReadable
	{
	public:
		/** The number of values the argument takes */
		unsigned int count;

		/** The values the parameter has been set to */
		Type *values;

		/** Constructor with the name of the argument */
		CmdLineParameters( const std::string &name );

		/** Destructor deallocating the array of values */
		~CmdLineParameters( void );

		/** Try to set the argument from the list of command line arguments.
		*** Returns thenumber of arguments ingested.*/ 
		int read( char **argv , int argc );

		/** Method for resizing the contents */
		void resize( unsigned int sz );
	};

	/** This function takes a list of arguments and tries to set the parameters.
	*** The last parameter must be a NULL pointer. */
	void CmdLineParse( int argc , char **argv, CmdLineReadable **params );

	/** This function takes a vectros of arguments and tries to set the parameters. */
	void CmdLineParse( int argc , char **argv, const std::vector< CmdLineReadable * > &params );

	/** Converts a string to upper case*/
	std::string ToUpper( const std::string &str );

	/** Converts a string to lower case*/
	std::string ToLower( const std::string &str );

	/** Returns the file extension */
	std::string GetFileExtension( const std::string &fileName );

	/** Returns and array of individual lines read from a file */
	std::vector< std::string > ReadLines( const std::string &fileName );

	/** Returns and array of individual words pulled from a string */
	std::vector< std::string > GetWords( const std::string &str );

	/** Returns and array of individual words read from a file */
	std::vector< std::string > ReadWords( const std::string &fileName );

#include "CmdLineParser.inl"
}
#endif // CMD_LINE_PARSER_INCLUDED
