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

void StitchingVisualization::LoadGeometryData() {
	glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, NULL);

	glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_FLOAT, 0, NULL);

	glBindBuffer(GL_ARRAY_BUFFER, coordinateBuffer);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	glTexCoordPointer(2, GL_FLOAT, 0, NULL);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, faceBuffer);
}

void StitchingVisualization::display( void )
{
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glEnable( GL_TEXTURE_2D );
	glEnable( GL_DEPTH_TEST );

	if( displayMode==ONE_REGION_DISPLAY )
	{
		if( showSource )
		{
			GLuint mBuffer = visualizationMode==MULTIPLE_INPUT_MODE ? referenceConfidenceBuffers[referenceIndex] : maskTextureBuffer;
			GLuint tBuffer = visualizationMode==MULTIPLE_INPUT_MODE ? referenceTextureBuffers[referenceIndex] : compositeTextureBuffer;
			setViewport( 0 );
			DrawRegion( showMesh , showMask ? mBuffer : tBuffer , false , false );
		}
		else
		{
#if 1
			// Show the mask / confidence
#ifdef SHOW_CUMULATIVE_CONFIDENCE
			GLuint mBuffer = visualizationMode==MULTIPLE_INPUT_MODE ? cumulativeConfidenceBuffer : maskTextureBuffer;
#else // !SHOW_CUMULATIVE_CONFIDENCE
			GLuint mBuffer = visualizationMode==MULTIPLE_INPUT_MODE ? referenceConfidenceBuffers[referenceIndex] : maskTextureBuffer;
#endif // SHOW_CUMULATIVE_CONFIDENCE
#else
			// Show the result
			GLuint mBuffer = visualizationMode==MULTIPLE_INPUT_MODE ? textureBuffer : compositeTextureBuffer;
#endif
			GLuint tBuffer = textureBuffer;
			setViewport( 0 );
			DrawRegion( showMesh , showMask ? mBuffer : tBuffer , false , false );
		}
	}
	else
	{
		// Show the output texture
		{
#if 1
			// Show the mask / confidence
#ifdef SHOW_CUMULATIVE_CONFIDENCE
			GLuint mBuffer = visualizationMode==MULTIPLE_INPUT_MODE ? cumulativeConfidenceBuffer : maskTextureBuffer;
#else // !SHOW_CUMULATIVE_CONFIDENCE
			GLuint mBuffer = visualizationMode==MULTIPLE_INPUT_MODE ? referenceConfidenceBuffers[referenceIndex] : maskTextureBuffer;
#endif // SHOW_CUMULATIVE_CONFIDENCE
#else
			// Show the result
			GLuint mBuffer = visualizationMode==MULTIPLE_INPUT_MODE ? textureBuffer : compositeTextureBuffer;
#endif
			GLuint tBuffer = textureBuffer;
			setViewport( 1 );
			DrawRegion( showMesh , showMask ? mBuffer : tBuffer , false , false );
		}

		// Show the input texture
		{
			GLuint mBuffer = visualizationMode==MULTIPLE_INPUT_MODE ? referenceConfidenceBuffers[referenceIndex] : maskTextureBuffer;
			GLuint tBuffer = visualizationMode==MULTIPLE_INPUT_MODE ? referenceTextureBuffers[referenceIndex] : compositeTextureBuffer;
			setViewport( 0 );
			DrawRegion( showMesh , showMask ? mBuffer : tBuffer , false , false );
		}
	}

	glDisable( GL_TEXTURE_2D );
	glBindTexture( GL_TEXTURE_2D , 0 );
	
	setViewport();

	if( showDisk && isBrushActive )
	{
		glDisable( GL_DEPTH_TEST );
		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		glOrtho(0, screenWidth, 0, screenHeight, -1, 1);

		glMatrixMode( GL_MODELVIEW );
		glLoadIdentity();

		glColor3f(0, 1.0, 0);
		GLUquadric* quad = gluNewQuadric();
		glTranslatef(diskX, screenHeight - diskY, 0);
		gluDisk(quad, 18, 22, 40, 3);
		gluDeleteQuadric(quad);
		glEnable( GL_DEPTH_TEST );
	}
}

void StitchingVisualization::UpdateColorTextureBuffer( void )
{
	if( !glIsBuffer(textureBuffer ) ) glGenTextures(1, &textureBuffer);
	glBindTexture( GL_TEXTURE_2D , textureBuffer );
	glTexImage2D( GL_TEXTURE_2D , 0 , GL_RGBA , textureWidth , textureHeight , 0 , GL_RGB , GL_UNSIGNED_BYTE , (GLvoid*)&colorTextureBuffer[0] );
	glBindTexture( GL_TEXTURE_2D , 0 );
}

template< typename Real >
void StitchingVisualization::UpdateCompositeTextureBuffer( const RegularGrid< 2 , Point3D< Real > > &composite )
{
	if( !glIsBuffer(compositeTextureBuffer) ) glGenTextures(1, &compositeTextureBuffer);
	glBindTexture(GL_TEXTURE_2D, compositeTextureBuffer);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	unsigned char * imValues = new unsigned char[composite.size() * 3];
	for( int j=0 ; j<composite.size() ; j++ ) for( int c=0 ; c<3 ; c++ ) imValues[3 * j + c] = (unsigned char)(composite[j][c] * 255.0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, composite.res(0), composite.res(1), 0, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)&imValues[0]);
	delete[] imValues;

	glBindTexture(GL_TEXTURE_2D, 0);
}

template< typename Real >
void StitchingVisualization::UpdateMaskTextureBuffer( const RegularGrid< 2 , Point3D< Real > > &mask )
{
	if( !glIsBuffer(maskTextureBuffer) ) glGenTextures(1, &maskTextureBuffer);
	glBindTexture(GL_TEXTURE_2D, maskTextureBuffer);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	unsigned char * imValues = new unsigned char[mask.size() * 3];
	for( int j=0 ; j<mask.size() ; j++ ) for( int c=0 ; c<3 ; c++ ) imValues[3 * j + c] = (unsigned char)(mask[j][c] * 255.0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, mask.res(0), mask.res(1), 0, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)&imValues[0]);
	delete[] imValues;

	glBindTexture(GL_TEXTURE_2D, 0);
}

template< typename Real >
void StitchingVisualization::UpdateReferenceTextureBuffers( const std::vector< RegularGrid< 2 , Point3D< Real > > > &images )
{
	referenceTextureBuffers.resize( images.size() );
	for( int i=0 ; i<referenceTextureBuffers.size() ; i++ )
	{
		if( !glIsBuffer( referenceTextureBuffers[i] ) ) glGenTextures( 1 , &referenceTextureBuffers[i] );

		glBindTexture(GL_TEXTURE_2D, referenceTextureBuffers[i]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		unsigned char * imValues = new unsigned char[images[i].size() * 3];
		for (int j = 0; j < images[i].size(); j++) for (int c = 0; c < 3; c++) imValues[3 * j + c] = (unsigned char)(images[i][j][c]*255.0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, images[i].res(0), images[i].res(1), 0, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)&imValues[0]);
		delete[] imValues;

		glBindTexture(GL_TEXTURE_2D, 0);
	}
}

template< typename Real >
void StitchingVisualization::UpdateReferenceConfidenceBuffers( const std::vector< RegularGrid< 2 , Real > > &confidences )
{
	referenceConfidenceBuffers.resize( confidences.size() );
	for( int i=0 ; i<referenceConfidenceBuffers.size() ; i++ )
	{
		if( !glIsBuffer( referenceConfidenceBuffers[i] ) ) glGenTextures( 1 , &referenceConfidenceBuffers[i] );

		glBindTexture( GL_TEXTURE_2D , referenceConfidenceBuffers[i] );
		glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_WRAP_S , GL_MIRRORED_REPEAT );
		glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_WRAP_T , GL_MIRRORED_REPEAT );

		glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_MIN_FILTER , GL_LINEAR );
		glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_MAG_FILTER , GL_LINEAR );

		unsigned char * imValues = new unsigned char[confidences[i].size() * 3];
		for( int j=0 ; j<confidences[i].size() ; j++ ) for( int c=0 ; c<3 ; c++ ) imValues[ 3*j + c ] = (unsigned char)(confidences[i][j]*255.0);
		glTexImage2D( GL_TEXTURE_2D , 0 , GL_RGBA , confidences[i].res(0) , confidences[i].res(1) , 0 , GL_RGB , GL_UNSIGNED_BYTE , (GLvoid*)&imValues[0] );
		delete[] imValues;

		glBindTexture( GL_TEXTURE_2D , 0 );
	}
#ifdef SHOW_CUMULATIVE_CONFIDENCE
	RegularGrid< 2 , Real > confidence( confidences[0].res() );
	for( size_t j=0 ; j<confidence.size() ; j++ ) confidence[j] = 0.;
	for( unsigned int i=0 ; i<confidences.size() ; i++ ) for( size_t j=0 ; j<confidence.size() ; j++ ) confidence[j] += confidences[i][j];
	Real max = 0;
	for( size_t j=0 ; j<confidence.size() ; j++ ) max = std::max< Real >( max , confidence[j] );
	for( size_t j=0 ; j<confidence.size() ; j++ ) confidence[j] /= max;
	if( !glIsBuffer( cumulativeConfidenceBuffer ) ) glGenTextures( 1 , &cumulativeConfidenceBuffer );

	glBindTexture( GL_TEXTURE_2D , cumulativeConfidenceBuffer );
	glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_WRAP_S , GL_MIRRORED_REPEAT );
	glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_WRAP_T , GL_MIRRORED_REPEAT );

	glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_MIN_FILTER , GL_LINEAR );
	glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_MAG_FILTER , GL_LINEAR );

	unsigned char * imValues = new unsigned char[confidence.size() * 3];
	for( int j=0 ; j<confidence.size() ; j++ ) for( int c=0 ; c<3 ; c++ ) imValues[ 3*j + c ] = (unsigned char)(confidence[j]*255.0);
	glTexImage2D( GL_TEXTURE_2D , 0 , GL_RGBA , confidence.res(0) , confidence.res(1) , 0 , GL_RGB , GL_UNSIGNED_BYTE , (GLvoid*)&imValues[0] );
	delete[] imValues;

	glBindTexture( GL_TEXTURE_2D , 0 );
#endif // SHOW_CUMULATIVE_CONFIDENCE
}

template< typename Real >
void StitchingVisualization::UpdateTextureBuffer( const RegularGrid< 2 , Point3D< Real > > &image )
{
	if (!glIsBuffer(textureBuffer)) {
		glGenTextures(1, &textureBuffer);
	}
	glBindTexture(GL_TEXTURE_2D, textureBuffer);

	// set basic parameters
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, image.res(0), image.res(1), 0, GL_RGB, GL_FLOAT, (GLvoid*)&image[0]);

	// Unbind the texture
	glBindTexture(GL_TEXTURE_2D, 0);

}

StitchingVisualization::StitchingVisualization( void )
{
	TexturedMeshVisualization();
	referenceIndex = 0;
	showDisk = true;
	isBrushActive = false;
	showMask = false;
	compositeTextureBuffer = 0;
	maskTextureBuffer = 0;
	showSource = false;

}