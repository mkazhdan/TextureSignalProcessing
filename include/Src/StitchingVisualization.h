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

#ifndef STITCHING_VISUALIZATION
#define STITCHING_VISUALIZATION

#include <Src/TexturedMeshVisualization.h>
enum
{
	COLOR_TEXTURE,
	NORMAL_TEXTURE,
	TEXTURE_TYPE_COUNT
};

class StitchingVisualization : public TexturedMeshVisualization
{
public:
	StitchingVisualization();
	int visualizationMode;
	bool isBrushActive;
	int diskX, diskY;
	bool showDisk;
	
	//Single input
	bool showMask;
	GLuint compositeTextureBuffer;
	GLuint maskTextureBuffer;
	template< typename Real >
	void UpdateCompositeTextureBuffer( const Image< Point3D< Real > > &composite );
	template< typename Real >
	void UpdateMaskTextureBuffer( const Image< Point3D< Real > > &mask );


	//Multiple input
	int referenceIndex;
	std::vector<GLuint> referenceTextureBuffers;
	template< typename Real >
	void UpdateReferenceTextureBuffers( const std::vector< Image< Point3D< Real > > > &images );
	std::vector< GLuint > referenceConfidenceBuffers;
	template< typename Real >
	void UpdateReferenceConfidenceBuffers( const std::vector< Image< Real > > &confidences );

	unsigned char * colorTextureBuffer;
	int textureWidth;
	int textureHeight;
	template< typename Real >
	void UpdateTextureBuffer( const Image< Point3D< Real > > &image );
	void UpdateColorTextureBuffer();

	void display(void);
	void LoadGeometryData();
};

#include "StitchingVisualization.inl"
#endif // STITCHING_VISUALIZATION