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

#ifndef EXCEPTIONS_INCLUDED
#define EXCEPTIONS_INCLUDED

#define VERBOSE_MESSAGING
//#define TIMED_MESSAGING

#include <stdarg.h>
#include <exception>
#include <string>
#include <sstream>
#include <iostream>
#include <mutex>
#ifdef TIMED_MESSAGING
#include "Timer.h"
Timer waningTimer;
#endif // TIMED_MESSAGING


namespace MishaK
{
	template< typename ... Arguments > void _AddToMessageStream( std::stringstream &stream , Arguments ... arguments );
	inline void _AddToMessageStream( std::stringstream &stream ){ return; }
	template< typename Argument , typename ... Arguments > void _AddToMessageStream( std::stringstream &stream , Argument argument , Arguments ... arguments )
	{
		stream << argument;
		_AddToMessageStream( stream , arguments ... );
	}

#ifdef VERBOSE_MESSAGING
	template< typename ... Arguments >
	std::string MakeMessageString( std::string header , std::string fileName , int line , std::string functionName , Arguments ... arguments )
	{
		size_t headerSize = header.size();
		std::stringstream stream;

		// The first line is the header, the file name , and the line number
		stream << header << " " << fileName << " (Line " << line << ")" << std::endl;

		// Inset the second line by the size of the header and write the function name
		for( size_t i=0 ; i<=headerSize ; i++ ) stream << " ";
		stream << functionName << std::endl;

		// Inset the third line by the size of the header and write the rest
		for( size_t i=0 ; i<=headerSize ; i++ ) stream << " ";
		_AddToMessageStream( stream , arguments ... );

		return stream.str();
	}
	struct Exception : public std::exception
	{
		const char *what( void ) const noexcept { return _message.c_str(); }
		template< typename ... Args >
		Exception( const char *fileName , int line , const char *functionName , const char *format , Args ... args )
		{
			_message = MakeMessageString( "[EXCEPTION]" , fileName , line , functionName , format , args ... );
		}
	private:
		std::string _message;
	};

	template< typename ... Args > void Throw( const char *fileName , int line , const char *functionName , const char *format , Args ... args ){ throw Exception( fileName , line , functionName , format , args ... ); }
	template< typename ... Args >
	void Warn( const char *fileName , int line , const char *functionName , const char *format , Args ... args )
	{
		static std::mutex warnMutex;
		std::lock_guard< std::mutex > lock( warnMutex );
#ifdef TIMED_MESSAGING
		std::cerr << MakeMessageString( "[WARNING]" , fileName , line , functionName , format , args ... ,  " (Time=" , waningTimer.elapsed() , ")" ) << std::endl;
#else // !TIMED_MESSAGING
		std::cerr << MakeMessageString( "[WARNING]" , fileName , line , functionName , format , args ... ) << std::endl;
#endif // TIMED_MESSAGING
	}
	template< typename ... Args >
	void ErrorOut( const char *fileName , int line , const char *functionName , const char *format , Args ... args )
	{
		static std::mutex errorOutMutex;
		std::lock_guard< std::mutex > lock( errorOutMutex );
		std::cerr << MakeMessageString( "[ERROR]" , fileName , line , functionName , format , args ... ) << std::endl;
		exit( 0 );
	}
#else // !VERBOSE_MESSAGING
	template< typename ... Arguments >
	std::string MakeMessageString( std::string header , std::string functionName , Arguments ... arguments )
	{
		std::stringstream stream;

		// The first line is the header, the file name , and the line number
		stream << header << " " << functionName << ": ";

		_AddToMessageStream( stream , arguments ... );

		return stream.str();
	}

	struct Exception : public std::exception
	{
		const char *what( void ) const noexcept { return _message.c_str(); }
		template< typename ... Args >
		Exception( const char *functionName , const char *format , Args ... args )
		{
			_message = MakeMessageString( "[EXCEPTION]" , functionName , format , args ... );
		}
	private:
		std::string _message;
	};
	template< typename ... Args > void Throw( const char *functionName , const char *format , Args ... args ){ throw Exception( functionName , format , args ... ); }
	template< typename ... Args >
	void Warn( const char *functionName , const char *format , Args ... args )
	{
		static std::mutex warnMutex;
		std::lock_guard< std::mutex > lock( warnMutex );
		std::cerr << MakeMessageString( "[WARNING]" , functionName , format , args ... ) << std::endl;
	}
	template< typename ... Args >
	void ErrorOut( const char *functionName , const char *format , Args ... args )
	{
		static std::mutex errorOutMutex;
		std::lock_guard< std::mutex > lock( errorOutMutex );
		std::cerr << MakeMessageString( "[WARNING]" , functionName , format , args ... ) << std::endl;
		exit( 0 );
	}
#endif // VERBOSE_MESSAGING
}
#ifdef VERBOSE_MESSAGING
#ifndef MK_WARN
#define MK_WARN( ... ) MishaK::Warn( __FILE__ , __LINE__ , __FUNCTION__ , __VA_ARGS__ )
#endif // MK_WARN
#ifndef MK_WARN_ONCE
#define MK_WARN_ONCE( ... ) { static bool firstTime = true ; if( firstTime ) MishaK::Warn( __FILE__ , __LINE__ , __FUNCTION__ , __VA_ARGS__ ) ; firstTime = false; }
#endif // MK_WARN_ONCE
#ifndef MK_THROW
#define MK_THROW( ... ) MishaK::Throw( __FILE__ , __LINE__ , __FUNCTION__ , __VA_ARGS__ )
#endif // MK_THROW
#ifndef MK_ERROR_OUT
#define MK_ERROR_OUT( ... ) MishaK::ErrorOut( __FILE__ , __LINE__ , __FUNCTION__ , __VA_ARGS__ )
#endif // MK_ERROR_OUT
#else // !VERBOSE_MESSAGING
#ifndef MK_WARN
#define MK_WARN( ... ) MishaK::Warn( __FUNCTION__ , __VA_ARGS__ )
#endif // MK_WARN
#ifndef MK_WARN_ONCE
#define MK_WARN_ONCE( ... ) { static bool firstTime = true ; if( firstTime ) MishaK::Warn( __FUNCTION__ , __VA_ARGS__ ) ; firstTime = false; }
#endif // MK_WARN_ONCE
#ifndef MK_THROW
#define MK_THROW( ... ) MishaK::Throw( __FUNCTION__ , __VA_ARGS__ )
#endif // MK_THROW
#ifndef MK_ERROR_OUT
#define MK_ERROR_OUT( ... ) MishaK::ErrorOut( __FUNCTION__ , __VA_ARGS__ )
#endif // MK_ERROR_OUT
#endif // VERBOSE_MESSAGING

#endif // EXCEPTIONS_INCLUDED
