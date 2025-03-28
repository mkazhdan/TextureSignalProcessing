/*
Copyright (c) 2024, Michael Kazhdan
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
#ifndef PRE_PROCESSING_INCLUDED
#define PRE_PROCESSING_INCLUDED

#define NEW_CODE					// General-purpose experimental code encapsulation

#define USE_EIGEN
#undef USE_CHOLMOD
#undef USE_EIGEN_PARDISO


//#define NO_OPEN_GL_VISUALIZATION		// Disable OpenGL visualization
//#define USE_TEXTURE_TRIANGLES		// Represent textures using a separate triangulation

#define INSERTION_EPSILON 1e-12		// Separation from interval end-points required for insertion


//#define DEBUG_BAD_LOOP
#ifdef DEBUG_BAD_LOOP
#include <functional>
std::function< bool ( unsigned int ) > DebugChartFunctor = [&]( unsigned int ){ return false; };
std::function< bool ( unsigned int ) > DebugCellFunctor = [&]( unsigned int ){ return false; };
bool DebugChart = false; 
bool DebugCell = false;
#endif // DEBUG_BAD_LOOP

#endif // PRE_PROCESSING_INCLUDED