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

void StitchingVisualization::display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_TEXTURE_2D);
	glEnable(GL_DEPTH_TEST);

	setViewport(0);
	DrawRegion(showMesh, textureBuffer, false,false);

	setViewport(1);
	if(visualizationMode == MULTIPLE_INPUT_MODE) DrawRegion(showMesh, referenceTextureBuffers[referenceIndex], false, false);
	else DrawRegion(showMesh, compositeTextureBuffer, false, false); 

	glDisable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 0);
	
	setViewport();

	if (showDisk && isBrushActive)
	{
		glDisable(GL_DEPTH_TEST);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(0, screenWidth, 0, screenHeight, -1, 1);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		glColor3f(0, 1.0, 0);
		GLUquadric* quad = gluNewQuadric();
		glTranslatef(diskX, screenHeight - diskY, 0);
		gluDisk(quad, 18, 22, 40, 3);
		gluDeleteQuadric(quad);
		glEnable(GL_DEPTH_TEST);
	}
}

void StitchingVisualization::UpdateColorTextureBuffer() {
	if (!glIsBuffer(textureBuffer)) {
		glGenTextures(1, &textureBuffer);
	}
	glBindTexture(GL_TEXTURE_2D, textureBuffer);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, textureWidth, textureHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)&colorTextureBuffer[0]);
	glBindTexture(GL_TEXTURE_2D, 0);
}

void StitchingVisualization::UpdateCompositeTextureBuffer(const Image<Point3D<float>> & composite) {
	if (!glIsBuffer(compositeTextureBuffer)) {
		glGenTextures(1, &compositeTextureBuffer);
	}
	glBindTexture(GL_TEXTURE_2D, compositeTextureBuffer);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	unsigned char * imValues = new unsigned char[composite.size() * 3];
	for (int j = 0; j < composite.size(); j++) for (int c = 0; c < 3; c++) imValues[3 * j + c] = (unsigned char)(composite[j][c] * 255.0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, composite.width(), composite.height(), 0, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)&imValues[0]);
	delete imValues;

	glBindTexture(GL_TEXTURE_2D, 0);
}

void StitchingVisualization::UpdateReferenceTextureBuffers(const std::vector<Image<Point3D<float>>> & images) {
	referenceTextureBuffers.resize(images.size());
	for (int i = 0; i < referenceTextureBuffers.size(); i++) {
		if (!glIsBuffer(referenceTextureBuffers[i])) {
			glGenTextures(1, &referenceTextureBuffers[i]);
		}
		glBindTexture(GL_TEXTURE_2D, referenceTextureBuffers[i]);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		unsigned char * imValues = new unsigned char[images[i].size() * 3];
		for (int j = 0; j < images[i].size(); j++) for (int c = 0; c < 3; c++) imValues[3 * j + c] = (unsigned char)(images[i][j][c]*255.0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, images[i].width(), images[i].height(), 0, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)&imValues[0]);
		delete imValues;

		glBindTexture(GL_TEXTURE_2D, 0);
	}
}

void StitchingVisualization::UpdateTextureBuffer(const Image<Point3D<float>> & image) {
	if (!glIsBuffer(textureBuffer)) {
		glGenTextures(1, &textureBuffer);
	}
	glBindTexture(GL_TEXTURE_2D, textureBuffer);

	// set basic parameters
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, textureWidth, textureHeight, 0, GL_RGB, GL_FLOAT, (GLvoid*)&textureImage[0]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, image.width(), image.height(), 0, GL_RGB, GL_FLOAT, (GLvoid*)&image[0]);

	// Unbind the texture
	glBindTexture(GL_TEXTURE_2D, 0);

}

StitchingVisualization::StitchingVisualization(void)
{
	TexturedMeshVisualization();
	referenceIndex = 0;
	showDisk = true;
	isBrushActive = false;
}