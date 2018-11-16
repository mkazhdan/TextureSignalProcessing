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

float TexturedMeshVisualization::imageToScreenScale( void ) const
{
	return std::min< float >(float(screenWidth) / float(textureImage.width()), float(screenHeight) / float(textureImage.height())) * xForm.zoom;
}
Point< float , 2 > TexturedMeshVisualization::imageToScreen( float px , float py ) const
{
	float ip[] = { px, py }, ic[] = { float(textureImage.width()) / 2 + xForm.offset[0], float(textureImage.height()) / 2 - xForm.offset[1] }, sc[] = { float(screenWidth) / 2, float(screenHeight) / 2 };
	float scale = imageToScreenScale();
	Point< float, 2 > sp;
	sp[0] = sc[0] + (ip[0] - ic[0])*scale;
	sp[1] = sc[1] - (ip[1] - ic[1])*scale;
	return sp;
}

void TexturedMeshVisualization::SetupOffScreenBuffer( void )
{
	// The depth buffer texture
	glGenTextures( 1 , &offscreen_depth_texture );
	glBindTexture( GL_TEXTURE_2D , offscreen_depth_texture );
	glTexStorage2D( GL_TEXTURE_2D , 1 , GL_DEPTH_COMPONENT24 , screenWidth , screenHeight );
	glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_MAG_FILTER , GL_LINEAR );
	glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_MIN_FILTER , GL_LINEAR );
 	glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_WRAP_S , GL_CLAMP_TO_BORDER );
	glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_WRAP_T , GL_CLAMP_TO_BORDER );

	// The color buffer texture
	glGenTextures( 1 , &offscreen_color_texture );
	glBindTexture( GL_TEXTURE_2D , offscreen_color_texture );
	glTexStorage2D( GL_TEXTURE_2D , 1 , GL_RGBA8 , screenWidth , screenHeight );
	glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_MAG_FILTER , GL_LINEAR );
	glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_MIN_FILTER , GL_LINEAR );
	glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_WRAP_S , GL_CLAMP_TO_BORDER );
	glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_WRAP_T , GL_CLAMP_TO_BORDER );

	// Create and set up the FBO
	glGenFramebuffers( 1 , &offscreen_framebuffer_handle );
	glBindFramebuffer( GL_FRAMEBUFFER , offscreen_framebuffer_handle );
	glFramebufferTexture2D( GL_FRAMEBUFFER , GL_DEPTH_ATTACHMENT , GL_TEXTURE_2D , offscreen_depth_texture , 0 );
	glFramebufferTexture2D( GL_FRAMEBUFFER , GL_COLOR_ATTACHMENT0 , GL_TEXTURE_2D , offscreen_color_texture , 0 );
	GLenum drawBuffers[] = { GL_COLOR_ATTACHMENT0 };
	glDrawBuffers( 1 , drawBuffers );

	glBindFramebuffer( GL_FRAMEBUFFER , 0 );
}

void TexturedMeshVisualization::RenderOffScreenBuffer( Image< Point3D< float > > & image )
{
	if( !offscreen_framebuffer_handle ) SetupOffScreenBuffer();
	setViewport();
	glBindFramebuffer( GL_FRAMEBUFFER , offscreen_framebuffer_handle );
	display();
	glFlush();

	//Save color buffer to image
	Pointer(float) GLColorBuffer = AllocPointer< float >( sizeof(float) * 3 * screenWidth * screenHeight );
	glReadBuffer( GL_COLOR_ATTACHMENT0 );
	glReadPixels( 0, 0 , screenWidth , screenHeight , GL_RGB , GL_FLOAT , GLColorBuffer );
	glFinish();
	image.resize( screenWidth , screenHeight );
	for( int i=0 ; i<screenWidth ; i++ ) for( int j=0 ; j<screenHeight ; j++ )  for( int c=0 ; c<3 ; c++ )
		image(i,j)[c] = GLColorBuffer[ c + i*3 + (screenHeight-1-j) * screenWidth * 3];
	FreePointer( GLColorBuffer );
	glBindFramebuffer( GL_FRAMEBUFFER , 0 );
}


void TexturedMeshVisualization::WriteSceneConfigurationCallBack( Visualization* v , const char* prompt )
{
	const TexturedMeshVisualization* av = (TexturedMeshVisualization*)v;
	FILE * file;
	file = fopen( prompt , "wb" );
	fwrite( &av->screenWidth  , sizeof( int ) , 1 , file );
	fwrite( &av->screenHeight , sizeof( int ) , 1 , file );
	av->camera.write( file );
	fwrite( &av->zoom , sizeof(float) , 1 , file );
	fclose(file);
}

void TexturedMeshVisualization::ReadSceneConfigurationCallBack( Visualization *v , const char* prompt )
{
	TexturedMeshVisualization* av = (TexturedMeshVisualization*)v;
	FILE * file;
	file = fopen( prompt , "rb" );
	if( !file ) Miscellany::Throw( "Camera Configuration File Not Valid: %s" , prompt );
	else
	{
		fread( &av->screenWidth  , sizeof(int) , 1 , file );
		fread( &av->screenHeight , sizeof(int) , 1 , file );
		av->camera.read( file );
		fread( &av->zoom , sizeof(float) , 1 , file );
		fclose( file );
	}
}
void TexturedMeshVisualization::ToggleVectorFieldCallBack( Visualization* v , const char* )
{
	TexturedMeshVisualization* tmv = (TexturedMeshVisualization*)v;
	tmv->showVectorField = !tmv->showVectorField;
}
void TexturedMeshVisualization::ScreenshotCallBack(Visualization* v, const char* prompt) {
	Image< Point3D< float > > image;
	TexturedMeshVisualization* av = (TexturedMeshVisualization*)v;
	av->RenderOffScreenBuffer( image );
	image.write( prompt );
}

void TexturedMeshVisualization::UpdateVertexBuffer() {
	if (!glIsBuffer(vertexBuffer)) {
		glGenBuffers(1, &vertexBuffer);
	}
	if (!glIsBuffer(normalBuffer)) {
		glGenBuffers(1, &normalBuffer);
	}
	if (!glIsBuffer(colorBuffer)) {
		glGenBuffers(1, &colorBuffer);
	}
	if (!glIsBuffer(coordinateBuffer)) {
		glGenBuffers(1, &coordinateBuffer);
	}

	glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
	glBufferData(GL_ARRAY_BUFFER, 3 * triangles.size() * sizeof(Point3D<float>), &vertices[0], GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
	glBufferData(GL_ARRAY_BUFFER, 3 * triangles.size() * sizeof(Point3D<float>), &normals[0], GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, colorBuffer);
	glBufferData(GL_ARRAY_BUFFER, 3 * triangles.size() * sizeof(Point3D<float>), &colors[0], GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, coordinateBuffer);
	glBufferData(GL_ARRAY_BUFFER, 3 * triangles.size() * sizeof(Point2D<float>), &textureCoordinates[0], GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glGenVertexArrays(1, &vertexHandle);
	glBindVertexArray(vertexHandle);

	{//Texture
	 // Vertex position
		glEnableVertexAttribArray(0);
		glBindVertexBuffer(0, vertexBuffer, 0, sizeof(GLfloat) * 3);
		glVertexAttribFormat(0, 3, GL_FLOAT, GL_FALSE, 0);
		glVertexAttribBinding(0, 0);

		// Vertex texture
		glEnableVertexAttribArray(1);
		glBindVertexBuffer(1, coordinateBuffer, 0, sizeof(GLfloat) * 2);
		glVertexAttribFormat(1, 2, GL_FLOAT, GL_FALSE, 0);
		glVertexAttribBinding(1, 1);

		// Vertex Normal
		glEnableVertexAttribArray(2);
		glBindVertexBuffer(2, normalBuffer, 0, sizeof(GLfloat) * 3);
		glVertexAttribFormat(2, 3, GL_FLOAT, GL_FALSE, 0);
		glVertexAttribBinding(2, 2);
	}

	glBindVertexArray(0);
}

void TexturedMeshVisualization::UpdateFaceBuffer() {
	if (!glIsBuffer(faceBuffer)) {
		glGenBuffers(1, &faceBuffer);
	}
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, faceBuffer);
	TriangleIndex *_triangles = new TriangleIndex[triangles.size()];

	for (int i = 0; i<triangles.size(); i++) for (int j = 0; j<3; j++) _triangles[i][j] = 3 * i + j;
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, triangles.size() * sizeof(int) * 3, _triangles, GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	delete[] _triangles;
}

void TexturedMeshVisualization::UpdateTextureBuffer( void )
{
	if( !glIsBuffer( textureBuffer ) ) glGenTextures( 1 , &textureBuffer );

	int height = textureImage.height();
	int width = textureImage.width();
	unsigned char *clrs = new unsigned char[height*width * 3];
	for( int j=0 ; j<height ; j++ ) for( int i=0 ; i<width ; i++ ) for( int c=0 ; c<3 ; c++ )
		clrs[ 3*(width*j + i)+c ] = (unsigned char)( std::min< float >( (float)textureImage(i,j)[c]*255.f , 255.f ) );
	glBindTexture( GL_TEXTURE_2D , textureBuffer );
	glTexImage2D( GL_TEXTURE_2D , 0 , GL_RGBA , width , height , 0 , GL_RGB , GL_UNSIGNED_BYTE , (GLvoid*)&clrs[0] );
	glBindTexture( GL_TEXTURE_2D , 0 );

	delete[] clrs;
}


void TexturedMeshVisualization::SetLightingData( void )
{
	glLightModeli( GL_LIGHT_MODEL_LOCAL_VIEWER , GL_FALSE );
	glLightModeli( GL_LIGHT_MODEL_TWO_SIDE , GL_TRUE );
	glLightfv( GL_LIGHT0 , GL_AMBIENT  , &lightAmbient [0] );
	glLightfv( GL_LIGHT0 , GL_DIFFUSE  , &lightDiffuse [0] );
	glLightfv( GL_LIGHT0 , GL_SPECULAR , &lightSpecular[0] );

	lightPosition[0] = -camera.forward[0];
	lightPosition[1] = -camera.forward[1];
	lightPosition[2] = -camera.forward[2];
	lightPosition[3] = 0.f;

	glLightfv( GL_LIGHT0 , GL_POSITION , &lightPosition[0] );
	glEnable( GL_LIGHTING );
	glEnable( GL_LIGHT0 );

	glMaterialfv( GL_FRONT_AND_BACK , GL_AMBIENT   , &shapeAmbient [0]       );
	glMaterialfv( GL_FRONT_AND_BACK , GL_DIFFUSE   , &shapeDiffuse [0]       );
	glMaterialfv( GL_FRONT_AND_BACK , GL_SPECULAR  , &shapeSpecular[0]       );
	glMaterialf ( GL_FRONT_AND_BACK , GL_SHININESS ,  shapeSpecularShininess );
}

void TexturedMeshVisualization::SetGeometryCamera( void )
{
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	GLint viewport[4];
	glGetIntegerv( GL_VIEWPORT , viewport );
	gluPerspective( FOV , (float)viewport[2]/(float)viewport[3] , NEARZ , FARZ );
	//Draw Camera
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	gluLookAt( camera.position[0] , camera.position[1] , camera.position[2] , camera.position[0] + camera.forward[0] , camera.position[1] + camera.forward[1] , camera.position[2] + camera.forward[2] , camera.up[0] , camera.up[1] , camera.up[2] );
}

void TexturedMeshVisualization::SetTextureCamera( void )
{
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	{
		GLint viewport[4];
		glGetIntegerv( GL_VIEWPORT , viewport );
		float ar = (float)viewport[2]/(float)viewport[3] , ar_r = 1.f/ar;
		if( viewport[2]>viewport[3] ) glOrtho( -ar*0.5 , ar*0.5 , -     0.5 ,      0.5 , -1.f , 1.f );
		else                          glOrtho( -   0.5 ,    0.5 , -ar_r*0.5 , ar_r*0.5 , -1.f , 1.f );
	}

	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
}

void TexturedMeshVisualization::PhongShading( GLuint & textureBufferId )
{
	static bool firstTime = true;
	GLSLProgram * current_program = normalProgram;
	current_program->use();

	GLdouble projection[16] , modelview[16];
	glGetDoublev( GL_PROJECTION_MATRIX , projection );
	glGetDoublev( GL_MODELVIEW_MATRIX , modelview );

	current_program->setUniformMatrix< 4 >( "eye_projection" , projection );
	current_program->setUniformMatrix< 4 >( "world_to_eye" , modelview );
	current_program->setUniform< 3 >( "light_direction" , &camera.forward[0] , firstTime );
	current_program->setUniform< 3 >( "light_diffuse"   , &lightDiffuse  [0] , firstTime );
	current_program->setUniform< 3 >( "light_specular"  , &lightSpecular [0] , firstTime );
	current_program->setUniform< 3 >( "light_ambient"   , &lightAmbient  [0] , firstTime );
	current_program->setUniform< 3 >( "shape_diffuse"   , &shapeDiffuse  [0] , firstTime );
	current_program->setUniform< 3 >( "shape_specular"  , &shapeSpecular [0] , firstTime );
	current_program->setUniform< 3 >( "shape_ambient"   , &shapeAmbient  [0] , firstTime );
	current_program->setUniform( "shape_specular_shininess" , shapeSpecularShininess , firstTime );

	current_program->setUniform( "normal_texture" , 0 , firstTime );
	glActiveTexture( GL_TEXTURE0 );
	glBindTexture( GL_TEXTURE_2D , textureBufferId );

	glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE , GL_DECAL );
	glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_MAG_FILTER , useNearestSampling ? GL_NEAREST : GL_LINEAR );
	glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_MIN_FILTER , useNearestSampling ? GL_NEAREST : GL_LINEAR );

	//		glPolygonOffset( 30.f , 30.f );

	glBindVertexArray( vertexHandle );
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , faceBuffer );

	glEnable( GL_POLYGON_OFFSET_FILL );
	glPolygonOffset( polygonOffsetFactor , polygonOffsetUnits );
	glDrawElements( GL_TRIANGLES , 3 * (int)triangles.size() , GL_UNSIGNED_INT , 0 );
	glDisable( GL_POLYGON_OFFSET_FILL );

	glBindVertexArray(0);
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , 0 );
	glBindTexture( GL_TEXTURE_2D , 0 );
	glUseProgram(0);

	firstTime = false;
}

void TexturedMeshVisualization::DrawGeometry( GLuint& textureBufferId , bool phongShading , bool modulateLight )
{
	SetGeometryCamera();

	if( showVectorField ) glDisable( GL_TEXTURE_2D );
	else                  glEnable ( GL_TEXTURE_2D );

	if( phongShading ) PhongShading( textureBufferId );
	else
	{
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

		glBindTexture(GL_TEXTURE_2D, textureBufferId);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, useNearestSampling ? GL_NEAREST : GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, useNearestSampling ? GL_NEAREST : GL_LINEAR);
		if (modulateLight) SetLightingData();
		glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, modulateLight ? GL_MODULATE : GL_DECAL);

		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(polygonOffsetFactor, polygonOffsetUnits);
		glDrawElements(GL_TRIANGLES, (GLsizei)(triangles.size() * 3), GL_UNSIGNED_INT, NULL);
		glDisable(GL_POLYGON_OFFSET_FILL);
		if (modulateLight) glDisable(GL_LIGHTING);
		glBindTexture(GL_TEXTURE_2D, 0);

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}

	glDisable(GL_TEXTURE_2D);

	if( showVectorField )
	{
		static const float LengthScale = 1.f/40.f;
		static const float WidthScale = LengthScale / 5.f;
		glDisable( GL_LIGHTING );
		glBegin( GL_TRIANGLES );
		for( int i=0 ; i<vectorField.size() ; i++ )
		{
			Point3D< float > t;
			int tIdx = vectorField[i].tIdx;
			Point3D< float > v[] = { vertices[ triangles[tIdx][0] ] , vertices[ triangles[tIdx][1] ] , vertices[ triangles[tIdx][2] ] };
			Point3D< float > n = Point3D< float >::CrossProduct( v[1]-v[0] , v[2]-v[0] );
			Point3D< float > center = v[0] + (v[1]-v[0])*vectorField[i].p[0] + (v[2]-v[0])*vectorField[i].p[1];
			Point3D< float > v1 = (v[1]-v[0])*vectorField[i].v[0] + (v[2]-v[0])*vectorField[i].v[1];
			Point3D< float > v2 = Point3D< float >::CrossProduct( v1 , n );
			v2 /= (float)Length(v2);

			glColor3f( 0.35f , 0.35f , 0.35f );
			t = center - v2*WidthScale ; glVertex3f( t[0] , t[1] , t[2] );
			t = center + v2*WidthScale ; glVertex3f( t[0] , t[1] , t[2] );
			glColor3f( 1 , 1 , 1 );
			t = center + v1*LengthScale ; glVertex3f( t[0] , t[1] , t[2] );
		}
		glEnd();
	}

	if (showBoundaryEdges)
	{
		glDisable(GL_LIGHTING);
		glLineWidth(lineWidth);
		glColor3f(1.0, 1.0, 1.0);
		glDisable(GL_MULTISAMPLE);
		glEnable(GL_BLEND);
		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glDepthMask(GL_FALSE);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glBegin(GL_LINES);
		for (int i = 0; i < boundaryEdgeVertices.size(); i++) glVertex3f(boundaryEdgeVertices[i][0], boundaryEdgeVertices[i][1], boundaryEdgeVertices[i][2]);
		glEnd();
		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_BLEND);
		glDepthMask(GL_TRUE);
	}

	if (showEdges)
	{
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

		glDisable(GL_LIGHTING);
		glLineWidth(lineWidth);
		glColor3f(0.125, 0.125, 0.125);
		glDisable(GL_MULTISAMPLE);
		glEnable(GL_BLEND);
		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glDepthMask(GL_FALSE);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDrawElements(GL_TRIANGLES, (GLsizei)(triangles.size() * 3), GL_UNSIGNED_INT, NULL);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_BLEND);
		glDepthMask(GL_TRUE);

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_NORMAL_ARRAY);
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}
}

//void TexturedMeshVisualization::DrawGeometry(GLuint & textureBufferId, bool modulateLight) {
//	SetGeometryCamera();
//
//	glEnable(GL_TEXTURE_2D);
//
//	glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
//	glEnableClientState(GL_VERTEX_ARRAY);
//	glVertexPointer(3, GL_FLOAT, 0, NULL);
//
//	glBindBuffer(GL_ARRAY_BUFFER, normalBuffer);
//	glEnableClientState(GL_NORMAL_ARRAY);
//	glNormalPointer(GL_FLOAT, 0, NULL);
//
//	glBindBuffer(GL_ARRAY_BUFFER, coordinateBuffer);
//	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
//	glTexCoordPointer(2, GL_FLOAT, 0, NULL);
//
//	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, faceBuffer);
//
//	glBindTexture(GL_TEXTURE_2D, textureBufferId);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, useNearestSampling ? GL_NEAREST : GL_LINEAR);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, useNearestSampling ? GL_NEAREST : GL_LINEAR);
//	if (modulateLight) SetLightingData();
//	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, modulateLight ? GL_MODULATE : GL_DECAL);
//
//	glEnable(GL_POLYGON_OFFSET_FILL);
//	glPolygonOffset(polygonOffsetFactor, polygonOffsetUnits);
//	glDrawElements(GL_TRIANGLES, (GLsizei)(triangles.size() * 3), GL_UNSIGNED_INT, NULL);
//	glDisable(GL_POLYGON_OFFSET_FILL);
//	if (modulateLight) glDisable(GL_LIGHTING);
//	glBindTexture(GL_TEXTURE_2D, 0);
//
//
//	if (showEdges) {
//		glDisable(GL_LIGHTING);
//		glLineWidth(lineWidth);
//		glColor3f(0.125, 0.125, 0.125);
//		glDisable(GL_MULTISAMPLE);
//		glEnable(GL_BLEND);
//		glEnable(GL_LINE_SMOOTH);
//		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
//		glDepthMask(GL_FALSE);
//		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//
//		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
//		glDrawElements(GL_TRIANGLES, (GLsizei)(triangles.size() * 3), GL_UNSIGNED_INT, NULL);
//		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
//
//		glDisable(GL_LINE_SMOOTH);
//		glDisable(GL_BLEND);
//		glDepthMask(GL_TRUE);
//	}
//
//	glDisableClientState(GL_VERTEX_ARRAY);
//	glDisableClientState(GL_NORMAL_ARRAY);
//	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
//	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
//
//
//	glDisable(GL_TEXTURE_2D);
//
//
//	if (showBoundaryEdges) {
//		glDisable(GL_LIGHTING);
//		glLineWidth(lineWidth);
//		glColor3f(0.8, 0.8, 0.8);
//		glDisable(GL_MULTISAMPLE);
//		glEnable(GL_BLEND);
//		glEnable(GL_LINE_SMOOTH);
//		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
//		glDepthMask(GL_FALSE);
//		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//		glBegin(GL_LINES);
//		for (int i = 0; i < boundaryEdgeVertices.size(); i++) glVertex3f(boundaryEdgeVertices[i][0], boundaryEdgeVertices[i][1], boundaryEdgeVertices[i][2]);
//		glEnd();
//		glDisable(GL_LINE_SMOOTH);
//		glDisable(GL_BLEND);
//		glDepthMask(GL_TRUE);
//	}
//
//}
void TexturedMeshVisualization::DrawRegion( bool drawGeometry , GLuint & textureBufferId , bool phongShading , bool modulateLight )
{
	if( drawGeometry ) DrawGeometry( textureBufferId , phongShading , modulateLight );
	else               DrawTexture ( textureBufferId );
}

void TexturedMeshVisualization::DrawTexture( GLuint & textureBufferId )
{
	SetTextureCamera();

	glEnable(GL_TEXTURE_2D);

	glBindTexture(GL_TEXTURE_2D, textureBufferId);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);


	float screenScale = std::max<float>(screenWidth, screenHeight);

	Point2D<float> q;
	glBegin(GL_QUADS);
	{
		glTexCoord2d(0, 0);
		q = (Point2D<float>(-0.5, -0.5) + Point2D<float>(xForm.offset[0], xForm.offset[1]) / screenScale)*xForm.zoom;
		glVertex2f(q[0], q[1]);

		glTexCoord2d(1, 0);
		q = (Point2D<float>(0.5, -0.5) + Point2D<float>(xForm.offset[0], xForm.offset[1]) / screenScale)*xForm.zoom;
		glVertex2f(q[0], q[1]);

		glTexCoord2d(1, 1);
		q = (Point2D<float>(0.5, 0.5) + Point2D<float>(xForm.offset[0], xForm.offset[1]) / screenScale)*xForm.zoom;
		glVertex2f(q[0], q[1]);

		glTexCoord2d(0, 1);
		q = (Point2D<float>(-0.5, 0.5) + Point2D<float>(xForm.offset[0], xForm.offset[1]) / screenScale)*xForm.zoom;
		glVertex2f(q[0], q[1]);
	}
	glEnd();

	glDisable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 0);

	if (showEdges) {
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_LIGHTING);
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();

		glScalef(xForm.zoom, xForm.zoom, 1.0);
		glTranslatef(xForm.offset[0] / screenScale, xForm.offset[1] / screenScale, 0.f);
		glTranslatef(-0.5, -0.5, 0.f);

		glBindBuffer(GL_ARRAY_BUFFER, coordinateBuffer);
		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(2, GL_FLOAT, 0, NULL);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, faceBuffer);

		glLineWidth(lineWidth);
		glColor3f(0.125, 0.125, 0.125);
		glDisable(GL_MULTISAMPLE);
		glEnable(GL_BLEND);
		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glDepthMask(GL_FALSE);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDrawElements(GL_TRIANGLES, (GLsizei)(triangles.size() * 3), GL_UNSIGNED_INT, NULL);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_BLEND);
		glDepthMask(GL_TRUE);

		glDisableClientState(GL_VERTEX_ARRAY);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		glPopMatrix();
	}
}


void TexturedMeshVisualization::reshape( int w , int h ){ screenWidth = w , screenHeight = h; }
void TexturedMeshVisualization::setViewport( int whichRegion )
{
	switch( whichRegion )
	{
	case 0:
	{
		switch( displayMode )
		{
		case 1: glViewport( 0 , 0 , screenWidth , screenHeight ) ; break;
		case 2: glViewport( 0 , 0 , screenWidth /2 , screenHeight ) ; break;
		case 3: glViewport( 0 , 0 , screenWidth /3 , screenHeight ) ; break;
		case 4: glViewport( 0 , 0 , (2*screenWidth)/5 , screenHeight ) ; break;
		default: glViewport( 0 , 0 , screenWidth , screenHeight );
		};
	}
	break;
	case 1:
	{
		switch( displayMode )
		{
		case 2: glViewport( screenWidth/2 , 0 , screenWidth /2 , screenHeight ) ; break;
		case 3: glViewport( screenWidth/3 , 0 , screenWidth /3 , screenHeight ) ; break;
		case 4: glViewport( (2*screenWidth)/5 , 0 , screenWidth/5 , screenHeight/2 ) ; break;
		default: glViewport( 0 , 0 , screenWidth , screenHeight );
		};
	}
	break;
	case 2:
	{
		switch( displayMode )
		{
		case 3: glViewport( 2*screenWidth/3 , 0 , screenWidth /3 , screenHeight ) ; break;
		case 4: glViewport( (2*screenWidth)/5 , screenHeight/2 , screenWidth/5 , screenHeight/2 ) ; break;
		default: glViewport( 0 , 0 , screenWidth , screenHeight );
		};
	}
	break;
	case 3:
	{
		switch( displayMode )
		{
		case 4: glViewport( (3*screenWidth)/5 , 0 , (2*screenWidth)/5 , screenHeight ) ; break;
		default: glViewport( 0 , 0 , screenWidth , screenHeight );
		};
	}
	break;
	default: glViewport( 0 , 0 , screenWidth , screenHeight );
	};
}

void TexturedMeshVisualization::display( void )
{
	glClearColor(1, 1, 1, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if( displayMode == ONE_REGION_DISPLAY )
	{
		glEnable( GL_DEPTH_TEST );

		setViewport(0);
		DrawRegion(showMesh, textureBuffer, false, true);
	}
	else if( displayMode == TWO_REGION_DISPLAY )
	{
		glEnable( GL_DEPTH_TEST );

		setViewport(0);
		DrawRegion( true , textureBuffer , false , true );

		setViewport(1);
		DrawRegion( false , textureBuffer , false , true );

		setViewport();
	}
}


void TexturedMeshVisualization::keyboardFunc( unsigned char /*key*/ , int /*x*/ , int /*y*/ ){}

void TexturedMeshVisualization::mouseFunc( int button , int /*state*/ , int x , int y )
{
	if (!showMesh) {
		newX = x; newY = y;

		scaling = panning = false;
		if (button == GLUT_LEFT_BUTTON) panning = true;
		else if (button == GLUT_RIGHT_BUTTON) scaling = true;
	}
	else {
		newX = x; newY = y;

		rotating = scaling = panning = false;
		if (button == GLUT_LEFT_BUTTON)
			if (glutGetModifiers() & GLUT_ACTIVE_CTRL) panning = true;
			else                                        rotating = true;
		else if (button == GLUT_RIGHT_BUTTON) scaling = true;
	}
}

void TexturedMeshVisualization::motionFunc(int x, int y)
{
	if( !showMesh )
	{
		oldX = newX, oldY = newY, newX = x, newY = y;

		if (panning) xForm.offset[0] -= (newX - oldX) / imageToScreenScale(), xForm.offset[1] += (newY - oldY) / imageToScreenScale();
		else
		{
			float dz = (float)( pow( 1.1 , (float)(newY-oldY) / 8.f ) );
			xForm.zoom *= dz;
		}
		glutPostRedisplay();
	}
	else
	{
		oldX = newX, oldY = newY, newX = x, newY = y;
		int screenSize = std::min< int >(screenWidth, screenHeight);
		float rel_x = (newX - oldX) / (float)screenSize * 2;
		float rel_y = (newY - oldY) / (float)screenSize * 2;
		float pRight = -rel_x * zoom, pUp = rel_y * zoom;
		float pForward = rel_y * zoom;
		float rRight = rel_y, rUp = rel_x;

		//if (rotating) camera.rotateUp(rUp), camera.rotateRight(rRight);
		//else if (scaling) camera.translate(camera.forward *pForward);
		//else if (panning) camera.translate(camera.right * pRight + camera.up * pUp);
		
		if     ( rotating ) camera.rotateUp( -rUp ) , camera.rotateRight( -rRight );
		else if( scaling  ) camera.translate( camera.forward*pForward );
		else if( panning  ) camera.translate( -( camera.right*pRight + camera.up*pUp ) );

		glutPostRedisplay();
	}
}

TexturedMeshVisualization::TexturedMeshVisualization( bool hasVectorField )
{
	showVectorField = false;
	useNearestSampling = false;
	showEdges = false;
	showBoundaryEdges = false;
	showMesh = true;
	screenHeight = 800;
	screenWidth = 800;

	zoom = 1.f;
	rotating = scaling = panning = false;

	radius = 1.f;

	camera = Camera< float >( Point3D< float >( 0.f , 0.f , 2.f ) , Point3D< float >( 0.f , 0.f , -1.f ) , Point3D< float >( 0.f , 1.f , 0.f ) );

	callBacks.push_back(KeyboardCallBack(this, 'C', "read camera", "File Name", ReadSceneConfigurationCallBack));
	callBacks.push_back(KeyboardCallBack(this, 'c', "save camera", "File Name", WriteSceneConfigurationCallBack));
	callBacks.push_back(KeyboardCallBack(this, 'e', "show edges", ShowEdgesCallBack));
	callBacks.push_back(KeyboardCallBack(this, 'm', "toggle mesh-atlas", ToggleVisualizationMode));
	callBacks.push_back(KeyboardCallBack(this, 'K', "screenshot", "File Name", ScreenshotCallBack));
	callBacks.push_back(KeyboardCallBack(this, 'E', "show boundary", ShowBoundaryEdgesCallBack));
	callBacks.push_back(KeyboardCallBack(this, 'N', "use nearest", NearestSamplingCallBack));
	if( hasVectorField ) callBacks.push_back( KeyboardCallBack( this , 'V', "toggle vector field ", ToggleVectorFieldCallBack ) );


	lightDiffuse [0] = lightDiffuse [1] = lightDiffuse [2] = 1.f, lightDiffuse [3] = 1.f;
	lightAmbient [0] = lightAmbient [1] = lightAmbient [2] = 0.f, lightAmbient [3] = 1.f;
	lightSpecular[0] = lightSpecular[1] = lightSpecular[2] = 1.f, lightSpecular[3] = 1.f;
	lightPosition[3] = 0.f;
	shapeDiffuse [0] = shapeDiffuse [1] = shapeDiffuse [2] = 1.f , shapeDiffuse [3] = 1.f;
	shapeAmbient [0] = shapeAmbient [1] = shapeAmbient [2] = 0.f , shapeAmbient [3] = 1.f;
	shapeSpecular[0] = shapeSpecular[1] = shapeSpecular[2] = 1.f , shapeSpecular[3] = 1.f;
	shapeSpecularShininess = 128.f;
}


//(x,y) in [0,screenWidth-1] x [0,screenHeight-1]
//(i,j) in [0,textureImage.width() -1] x [0,textureImage.height() - 1]

Point< float, 2 > TexturedMeshVisualization::screenToImage(int x, int  y) {
	float ic[2] = { float(textureImage.width()) / 2 + xForm.offset[0], float(textureImage.height()) / 2 - xForm.offset[1] };
	float sc[2] = { float(screenWidth) / 2, float(screenHeight) / 2 };
	float scale = imageToScreenScale();
	Point< float, 2 > ip;
	ip[0] = ((float(x) - sc[0]) / scale) + ic[0];
	ip[1] = ((float(y) - sc[1]) / scale) + ic[1];
	return ip;
}

Point2D<float> TexturedMeshVisualization::selectImagePos(int x, int  y) {
	Point2D<float> imagePos;
	if (displayMode == ONE_REGION_DISPLAY) {
		imagePos = Point2D<float>(x / float(screenWidth), y / float(screenHeight));
	}
	else if (displayMode == TWO_REGION_DISPLAY) {
		imagePos = Point2D<float>((x - (screenWidth / 2)) / float(screenWidth / 2), y / float(screenHeight));
	}
	imagePos -= Point2D<float>(0.5, 0.5);

	//Reverse transformation
	float screenScale = std::max<float>(screenWidth, screenHeight);
	imagePos = imagePos / xForm.zoom - Point2D<float>(xForm.offset[0], xForm.offset[1]) / screenScale + Point2D<float>(0.5, 0.5);

	return imagePos;
}

bool TexturedMeshVisualization::select( int x , int  y , Point3D< float >& out ) //PERSPECTIVE
{
	setViewport( 0 );
	SetGeometryCamera();

	GLint viewport[4];
	glGetIntegerv( GL_VIEWPORT , viewport );

	GLdouble modelview[16];
	glGetDoublev( GL_MODELVIEW_MATRIX , modelview );

	GLdouble projection[16];
	glGetDoublev( GL_PROJECTION_MATRIX , projection );

	GLfloat winX = (GLfloat)x , winY = (GLfloat)y , winZ;
	winY = (float)viewport[3] - winY;

	glReadPixels( (int)winX , (int)winY , 1 , 1 , GL_DEPTH_COMPONENT , GL_FLOAT , &winZ );
	if( winZ==1.f ) return false;

	GLdouble posX , posY , posZ;
	gluUnProject( winX , winY , winZ , modelview , projection , viewport , &posX , &posY , &posZ );

	out = Point3D< float >( (float)posX , (float)posY , (float)posZ );

	setViewport();

	return true;
}