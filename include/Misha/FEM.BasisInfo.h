/*
Copyright (c) 2015, Michael Kazhdan and Fabian Prada
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

template<> const unsigned int BasisInfo< BASIS_0_WHITNEY           >::CoefficientsPerElement = 1;
template<> const unsigned int BasisInfo< BASIS_1_CONFORMING        >::CoefficientsPerElement = 2;
template<> const unsigned int BasisInfo< BASIS_1_WHITNEY           >::CoefficientsPerElement = 1;
template<> const unsigned int BasisInfo< BASIS_1_TRIANGLE_CONSTANT >::CoefficientsPerElement = 2;
template<> const unsigned int BasisInfo< BASIS_2_WHITNEY           >::CoefficientsPerElement = 1;
template<> const unsigned int BasisInfo< BASIS_2_VERTEX_CONSTANT   >::CoefficientsPerElement = 1;

template<> const unsigned int BasisInfo< BASIS_0_WHITNEY           >::Coefficients = 3;
template<> const unsigned int BasisInfo< BASIS_1_CONFORMING        >::Coefficients = 6;
template<> const unsigned int BasisInfo< BASIS_1_WHITNEY           >::Coefficients = 3;
template<> const unsigned int BasisInfo< BASIS_1_TRIANGLE_CONSTANT >::Coefficients = 2;
template<> const unsigned int BasisInfo< BASIS_2_WHITNEY           >::Coefficients = 1;
template<> const unsigned int BasisInfo< BASIS_2_VERTEX_CONSTANT   >::Coefficients = 3;

template<> const unsigned int BasisInfo< BASIS_0_WHITNEY           >::Dimension = 0;
template<> const unsigned int BasisInfo< BASIS_1_CONFORMING        >::Dimension = 1;
template<> const unsigned int BasisInfo< BASIS_1_WHITNEY           >::Dimension = 1;
template<> const unsigned int BasisInfo< BASIS_1_TRIANGLE_CONSTANT >::Dimension = 1;
template<> const unsigned int BasisInfo< BASIS_2_WHITNEY           >::Dimension = 2;
template<> const unsigned int BasisInfo< BASIS_2_VERTEX_CONSTANT   >::Dimension = 2;

template<> const bool BasisInfo< BASIS_0_WHITNEY           >::Lumpable = true;
template<> const bool BasisInfo< BASIS_1_CONFORMING        >::Lumpable = false;
template<> const bool BasisInfo< BASIS_1_WHITNEY           >::Lumpable = true;
template<> const bool BasisInfo< BASIS_1_TRIANGLE_CONSTANT >::Lumpable = false;
template<> const bool BasisInfo< BASIS_2_WHITNEY           >::Lumpable = true;
template<> const bool BasisInfo< BASIS_2_VERTEX_CONSTANT   >::Lumpable = true;

template<> const unsigned int BasisInfo< BASIS_0_WHITNEY           >::ElementType = ELEMENT_VERTEX;
template<> const unsigned int BasisInfo< BASIS_1_CONFORMING        >::ElementType = ELEMENT_VERTEX;
template<> const unsigned int BasisInfo< BASIS_1_WHITNEY           >::ElementType = ELEMENT_EDGE;
template<> const unsigned int BasisInfo< BASIS_1_TRIANGLE_CONSTANT >::ElementType = ELEMENT_TRIANGLE;
template<> const unsigned int BasisInfo< BASIS_2_WHITNEY           >::ElementType = ELEMENT_TRIANGLE;
template<> const unsigned int BasisInfo< BASIS_2_VERTEX_CONSTANT   >::ElementType = ELEMENT_VERTEX;

template<> const bool BasisInfo< BASIS_0_WHITNEY           >::Singular = false;
template<> const bool BasisInfo< BASIS_1_CONFORMING        >::Singular = true;
template<> const bool BasisInfo< BASIS_1_WHITNEY           >::Singular = false;
template<> const bool BasisInfo< BASIS_1_TRIANGLE_CONSTANT >::Singular = false;
template<> const bool BasisInfo< BASIS_2_WHITNEY           >::Singular = false;
template<> const bool BasisInfo< BASIS_2_VERTEX_CONSTANT   >::Singular = false;
