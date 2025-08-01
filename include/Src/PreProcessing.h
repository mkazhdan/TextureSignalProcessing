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

// To do:
// -- Pull edges out of the divergence operator

// To do:
// -- Use analytic integration
// -- Remove GridChart::cellType and GridChart::texelType
// -- Change coefficient vector std::vector -> IndexVector (in IterativeSolvers.inl)
// -- Minimize static_cast< unsigned int >(...)
// 1. Modify code to distinguish between grid-based indexing and normalized coordinates
#ifndef PRE_PROCESSING_INCLUDED
#define PRE_PROCESSING_INCLUDED

#define NEW_CODE						// General-purpose experimental code encapsulation

//#define NO_OPEN_GL_VISUALIZATION		// Disable OpenGL visualization
//#define DEBUG_INDEXING				// Use separate classes to sanity check indexing
//#define SANITY_CHECK					// Enables sanity checks for debugging purposes

#define USE_EIGEN
//#define USE_CHOLMOD
//#define USE_EIGEN_PARDISO


#define INSERTION_EPSILON 1e-12			// Separation from interval end-points required for insertion
#define MIN_TEXEL_WEIGHT 1e-12
#define DEFAULT_JITTER 1e-6
#define ORTHOGONAL_PERTURBATION 1e-10
#define SANITY_PRECISION_EPSILON 1e-10
#define PRECISION_EPSILON 1e-10
#define PROLONGATION_EPSILON 1e-10
#define DIAGONAL_CLAMP 1e-10

#endif // PRE_PROCESSING_INCLUDED