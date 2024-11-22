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

void TextureFilteringVisualization::LoadGeometryData() {
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
int TextureFilteringVisualization::slideBarWidth( void )
{
	GLint oldViewport[4] ,  viewport[4];
	glGetIntegerv( GL_VIEWPORT , oldViewport );
	setViewport( 0 );
	glGetIntegerv( GL_VIEWPORT , viewport );
	glViewport( oldViewport[0] , oldViewport[1] , oldViewport[2] , oldViewport[3] );
	return viewport[2];
}
void TextureFilteringVisualization::DrawSlideBar( void )
{
	glDisable( GL_LIGHTING );
	glDisable( GL_DEPTH_TEST );

	GLint viewport[4];
	glGetIntegerv( GL_VIEWPORT , viewport );
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	glOrtho( 0 , viewport[2] , 0 , viewport[3] , -1 ,  1 );

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glBegin( GL_QUADS );
	glColor3f( 0.0f , 0.0f , 1.0f ) ; glVertex2f( 20 , viewport[3] - 20 );
	glColor3f( 0.8f , 0.8f , 0.8f ) ; glVertex2f( viewport[2] / 2 , viewport[3] - 20 );
	glColor3f( 0.8f , 0.8f , 0.8f ) ; glVertex2f( viewport[2] / 2 , viewport[3] - 30 );
	glColor3f( 0.0f , 0.0f , 1.0f ) ; glVertex2f( 20 , viewport[3] - 30 );
	glEnd();

	glBegin(GL_QUADS);
	glColor3f( 0.8f , 0.8f , 0.8f ) ; glVertex2f( viewport[2] / 2  , viewport[3] - 20 );
	glColor3f( 1.0f , 0.0f , 0.0f ) ; glVertex2f( viewport[2] - 20 , viewport[3] - 20 );
	glColor3f( 1.0f , 0.0f , 0.0f ) ; glVertex2f( viewport[2] - 20 , viewport[3] - 30 );
	glColor3f( 0.8f , 0.8f , 0.8f ) ; glVertex2f( viewport[2] / 2  , viewport[3] - 30 );
	glEnd();

	float normalizedCursorPosition = 20 * (1.0 - slideBarCursorPosition) + ( viewport[2] - 20) * slideBarCursorPosition;

	glBegin(GL_QUADS);
	glColor3f(0.0, 0.0, 0.0);
	glVertex2f(normalizedCursorPosition - 12, viewport[3] - 18);
	glVertex2f(normalizedCursorPosition + 12, viewport[3] - 18);
	glVertex2f(normalizedCursorPosition + 12, viewport[3] - 32);
	glVertex2f(normalizedCursorPosition - 12, viewport[3] - 32);
	glEnd();

	glBegin(GL_QUADS);
	glColor3f(0.5, 0.5, 0.5);
	glVertex2f(normalizedCursorPosition - 10, viewport[3] - 20);
	glVertex2f(normalizedCursorPosition + 10, viewport[3] - 20);
	glVertex2f(normalizedCursorPosition + 10, viewport[3] - 30);
	glVertex2f(normalizedCursorPosition - 10, viewport[3] - 30);
	glEnd();
}

void TextureFilteringVisualization::display( void )
{
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	if( displayMode == ONE_REGION_DISPLAY )
	{
		glEnable(GL_TEXTURE_2D);
		glEnable(GL_DEPTH_TEST);

		setViewport(0);
		DrawRegion(showMesh, textureBuffer, textureType == NORMAL_TEXTURE);

		glDisable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 0);
	}
	else if (displayMode == TWO_REGION_DISPLAY) {
		SetGeometryCamera();
		LoadGeometryData();

		glEnable(GL_TEXTURE_2D);
		glEnable(GL_DEPTH_TEST);

		setViewport(0);
		glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glBindTexture(GL_TEXTURE_2D, textureBuffer);
		glDrawElements(GL_TRIANGLES, (GLsizei)(triangles.size() * 3), GL_UNSIGNED_INT, NULL);

		setViewport(1);
		SetLightingData();
		glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
		glBindTexture(GL_TEXTURE_2D, maskTextureBuffer);
		glDrawElements(GL_TRIANGLES, (GLsizei)(triangles.size() * 3), GL_UNSIGNED_INT, NULL);
		glDisable(GL_LIGHTING);
		glDisable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 0);
	}
	else if (displayMode == THREE_REGION_DISPLAY) {

		glEnable(GL_TEXTURE_2D);
		glEnable(GL_DEPTH_TEST);

		setViewport(0);
		DrawRegion(showMesh, textureBuffer, textureType == NORMAL_TEXTURE);

		setViewport(1);
		DrawRegion(false, textureBuffer, textureType == NORMAL_TEXTURE);
		//DrawRegion(showMesh, inputTextureBuffer, textureType == NORMAL_TEXTURE);

		setViewport(2);
		DrawRegion(showMesh, maskTextureBuffer, false, true);
		//DrawRegion(false, textureBuffer, textureType == NORMAL_TEXTURE);

		glDisable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 0);
	}

	else if( displayMode==FOUR_REGION_DISPLAY )
	{
		glEnable( GL_TEXTURE_2D );
		glEnable( GL_DEPTH_TEST );

		setViewport(0);
		DrawRegion( true , textureBuffer , textureType==NORMAL_TEXTURE );

		setViewport(1);
		DrawRegion( false , maskTextureBuffer );

		setViewport(2);
		DrawRegion( true , maskTextureBuffer , false , true );

		setViewport(3);
		DrawRegion( false , textureBuffer , textureType==NORMAL_TEXTURE );

		glDisable( GL_TEXTURE_2D );
		glBindTexture( GL_TEXTURE_2D, 0 );
	}

	if( showDisk && isBrushActive )
	{
		glDisable( GL_DEPTH_TEST );
		setViewport(0);
		GLint viewport[4];
		glGetIntegerv( GL_VIEWPORT , viewport );
		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		glOrtho( 0 , viewport[2] , 0 , viewport[3] , -1 ,  1 );

		glMatrixMode( GL_MODELVIEW );
		glLoadIdentity();

		if( positiveModulation) glColor3f( 1.0 , 0.0 , 0.0 );
		else                    glColor3f( 0.0 , 0.0 , 1.0 );
		GLUquadric* quad = gluNewQuadric();
		glTranslatef( diskX , screenHeight - diskY , 0);
		gluDisk( quad , 18 , 22 , 40 , 3 );
		gluDeleteQuadric( quad );
		glEnable( GL_DEPTH_TEST );
	}

	if( showSlideBar )
	{
		setViewport(0);
		DrawSlideBar();
	}
	setViewport();
}

void TextureFilteringVisualization::UpdateColorTextureBuffer() {
	if (!glIsBuffer(textureBuffer)) {
		glGenTextures(1, &textureBuffer);
	}
	glBindTexture(GL_TEXTURE_2D, textureBuffer);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, textureWidth, textureHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)&colorTextureBuffer[0]);
	glBindTexture(GL_TEXTURE_2D, 0);
}

void TextureFilteringVisualization::UpdateMaskTextureBuffer() {
	if (!glIsBuffer(maskTextureBuffer)) {
		glGenTextures(1, &maskTextureBuffer);
	}
	glBindTexture(GL_TEXTURE_2D, maskTextureBuffer);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, textureWidth, textureHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)&maskBufferValues[0]);
	glBindTexture(GL_TEXTURE_2D, 0);
}

template< typename Real >
void TextureFilteringVisualization::UpdateTextureBuffer( const Image< Point3D< Real > > &image )
{
	if( !glIsBuffer(textureBuffer) ) glGenTextures( 1 , &textureBuffer );
	glBindTexture(GL_TEXTURE_2D, textureBuffer);

	// set basic parameters
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, image.width(), image.height(), 0, GL_RGB, GL_FLOAT, (GLvoid*)&image[0]);

	// Unbind the texture
	glBindTexture(GL_TEXTURE_2D, 0);

}

TextureFilteringVisualization::TextureFilteringVisualization( void )
{
	TexturedMeshVisualization();

	showDisk = true;
	isBrushActive = false;
	showSlideBar = true;
	isSlideBarActive = false;
}