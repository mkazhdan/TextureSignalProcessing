/*
Copyright (c) 2018, Fabian Prada and Michael Kazhdan
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
#pragma once

#include <Eigen/Core>
#include <vectorclass/vectorclass.h>


template< class Real > struct ScalarPrecisionType{};
template< class Real > struct VectorPrecisionType{};

template<>
struct ScalarPrecisionType< float >
{
	typedef Eigen::VectorXf EigenVector;
	static float ZeroData( void ){ return 0.f; }
	static float  DotData( const float& v1 , const float& v2 ){ return v1*v2; }
};
template<>
struct ScalarPrecisionType< double >
{
	typedef Eigen::VectorXd EigenVector;
	static double ZeroData( void ){ return 0.; }
	static double  DotData( const double& v1 , const double& v2 ){ return v1*v2; }
};

template<>
struct VectorPrecisionType< float >
{
	typedef Vec4f VectorType;
	typedef Eigen::VectorXf EigenVector;
	static VectorType ZeroData( void ){ return VectorType( 0.f ); }
	static VectorType  SetData( float a , float b , float c ){ return VectorType( a , b , c , 0 ); }
	static float       DotData( const VectorType & v1 , const VectorType & v2 ){ return horizontal_add(v1*v2); }
};
template<>
struct VectorPrecisionType< double >
{
	typedef Vec4d VectorType;
	typedef Eigen::VectorXd EigenVector;
	static VectorType ZeroData( void ){ return VectorType( 0 ); }
	static VectorType  SetData( double a , double b , double c ){ return VectorType( a , b , c , 0 ); }
	static double      DotData( const VectorType & v1 , const VectorType & v2 ){ return horizontal_add(v1*v2); }
};


