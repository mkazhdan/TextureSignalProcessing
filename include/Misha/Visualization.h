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

#ifndef VISUALIZATION_INCLUDED
#define VISUALIZATION_INCLUDED

#include <GL/glew.h>
#include <GL/glut.h>
#include <vector>
#include <algorithm>
#include "PNG.h"
#include "Array.h"
#include "Miscellany.h"

#define KEY_UPARROW		101
#define KEY_DOWNARROW	103
#define KEY_LEFTARROW	100
#define KEY_RIGHTARROW	102
#define KEY_PGUP		104
#define KEY_PGDN		105
#define KEY_CTRL_C        3
#define KEY_BACK_SPACE    8
#define KEY_ENTER        13
#define KEY_ESC          27

struct Visualization
{
protected:
	const double _MIN_FPS_TIME = 0.5;
	double _lastFPSTime;
	int _lastFPSCount;
	double _fps;
	int _currentFrame , _frameStride , _totalFrames;
	bool _exitAfterVideo;
public:
	int screenWidth , screenHeight;
	void *font , *promptFont;
	int fontHeight , promptFontHeight;
	bool showHelp , showInfo , showFPS;
	void (*promptCallBack)( Visualization* , const char* );
	char promptString[1024];
	int promptLength;
	char* snapshotName;
	char videoHeader[1024];
	bool flushImage;

	struct KeyboardCallBack
	{
		char key;
		char prompt[1024];
		char description[1024];
		void (*callBackFunction)( Visualization* , const char* );
		Visualization* visualization;
		KeyboardCallBack( Visualization* visualization , char key , const char* description , void (*callBackFunction)( Visualization* , const char* ) )
		{
			this->visualization = visualization;
			this->key = key;
			strcpy( this->description , description );
			prompt[0] = 0;
			this->callBackFunction = callBackFunction;
		}
		KeyboardCallBack( Visualization* visualization , char key , const char* description , const char* prompt , void ( *callBackFunction )( Visualization* , const char* ) )
		{
			this->visualization = visualization;
			this->key = key;
			strcpy( this->description , description );
			strcpy( this->prompt , prompt );
			this->callBackFunction = callBackFunction;
		}
	};

	std::vector< KeyboardCallBack > callBacks;
	std::vector< char* > info;
	Visualization( void )
	{
		callBacks.push_back( KeyboardCallBack( this , KEY_ESC    , "" , QuitCallBack ) );
		callBacks.push_back( KeyboardCallBack( this , KEY_CTRL_C , "" , QuitCallBack ) );
		callBacks.push_back( KeyboardCallBack( this , 'F' , "toggle fps"  , ToggleFPSCallBack ) );
		callBacks.push_back( KeyboardCallBack( this , 'H' , "toggle help" , ToggleHelpCallBack ) );
		callBacks.push_back( KeyboardCallBack( this , 'I' , "toggle info" , ToggleInfoCallBack ) );
		callBacks.push_back( KeyboardCallBack( this , 'i' , "save frame buffer" , "<Ouput image>" , SetFrameBufferCallBack ) );
		callBacks.push_back( KeyboardCallBack( this , 'v' , "save video" , "<Ouput header> [<stride>] <number of frames>" , SetVideoCallBack ) );
		snapshotName = NULL;
		flushImage = false;
		screenWidth = screenHeight = 512;
		font = GLUT_BITMAP_HELVETICA_12;
		fontHeight = 12;
		promptFont = GLUT_BITMAP_TIMES_ROMAN_24;
		promptFontHeight = 24;
		showHelp = showInfo = showFPS = true;
		promptCallBack = NULL;
		strcpy( promptString , "" );
		promptLength = 0;

		_lastFPSTime = Miscellany::Time();
		_lastFPSCount = 0;
		_fps = 0;
	}
	virtual void idle( void ){;}
	virtual void keyboardFunc( unsigned char /*key*/ , int /*x*/ , int /*y*/ ){;}
	virtual void specialFunc( int /*key*/ , int /*x*/ , int /*y*/ ){;}
	virtual void display( void ){;}
	virtual void reshape( int /*w*/ , int /*h*/ ){;}
	virtual void mouseFunc( int /*button*/ , int /*state*/ , int /*x*/ , int /*y*/ ){;}
	virtual void motionFunc( int /*x*/ , int /*y*/ ){;}
	virtual void passiveMotionFunc( int /*x*/ , int /*y*/ ){;}

	void Idle        ( void );
	void KeyboardFunc( unsigned char key , int x , int y );
	void SpecialFunc ( int key, int x, int y );
	void Display     ( void );
	void Reshape     ( int w , int h );
	void MouseFunc   ( int button , int state , int x , int y );
	void MotionFunc  ( int x , int y );
	void PassiveMotionFunc  ( int x , int y );

	static void           QuitCallBack( Visualization*   , const char* ){ exit( 0 ); }
	static void      ToggleFPSCallBack( Visualization* v , const char* ){ v->showFPS  = !v->showFPS ; }
	static void     ToggleHelpCallBack( Visualization* v , const char* ){ v->showHelp = !v->showHelp; }
	static void     ToggleInfoCallBack( Visualization* v , const char* ){ v->showInfo = !v->showInfo; }
	static void SetFrameBufferCallBack( Visualization* v , const char* prompt )
	{
		if( prompt )
		{
			v->snapshotName = new char[ strlen(prompt)+1 ];
			strcpy( v->snapshotName , prompt );
			v->flushImage = true;
		}
	}
	static void SetVideoCallBack( Visualization* v , const char* prompt )
	{
		if( prompt )
		{
			if     ( sscanf( prompt , " %s %d %d " , v->videoHeader , &v->_frameStride , &v->_totalFrames )==3 ) ;
			else if( sscanf( prompt , " %s %d "    , v->videoHeader ,                    &v->_totalFrames )==2 ) v->_frameStride=1;
			else  v->_totalFrames = 0;
			v->_currentFrame = 0;
			v->_exitAfterVideo = true;
		}
	}

	static void WriteLeftString( int x , int y , void* font , const char* format , ... );
	static int StringWidth( void* font , const char* format , ... );
	void writeLeftString( int x , int y , const char* format , ... ) const;
	void writeRightString( int x , int y , const char* format , ... ) const;
	int stringWidth( const char* format , ... ) const;

	void saveFrameBuffer( const char* fileName , int whichBuffer=GL_BACK );
};
struct VisualizationViewer
{
	static Visualization* visualization;
	static void Idle             ( void );
	static void KeyboardFunc     ( unsigned char key , int x , int y );
	static void SpecialFunc      ( int key, int x, int y );
	static void Display          ( void );
	static void Reshape          ( int w , int h );
	static void MouseFunc        ( int button , int state , int x , int y );
	static void MotionFunc       ( int x , int y );
	static void PassiveMotionFunc( int x , int y );
};
Visualization* VisualizationViewer::visualization = NULL;
void VisualizationViewer::Idle( void ){ visualization->Idle(); }
void VisualizationViewer::KeyboardFunc( unsigned char key , int x , int y ){ visualization->KeyboardFunc( key , x , y ); }
void VisualizationViewer::SpecialFunc( int key , int x , int y ){ visualization->SpecialFunc( key , x ,  y ); }
void VisualizationViewer::Display( void ){ visualization->Display(); }
void VisualizationViewer::Reshape( int w , int h ){ visualization->Reshape( w , h ); }
void VisualizationViewer::MouseFunc( int button , int state , int x , int y ){ visualization->MouseFunc( button , state , x , y ); }
void VisualizationViewer::MotionFunc( int x , int y ){ visualization->MotionFunc( x , y ); }
void Visualization::Reshape( int w , int h )
{
	reshape( w , h );
	screenWidth = w , screenHeight = h;
	glViewport( 0 , 0 , screenWidth , screenHeight );
}
void Visualization::MouseFunc( int button , int state , int x , int y ){ mouseFunc( button , state , x , y ); }
void Visualization::MotionFunc( int x , int y ){ motionFunc( x , y );}
void Visualization::PassiveMotionFunc( int x , int y ){ passiveMotionFunc( x , y );}
void Visualization::Idle( void )
{
	if( snapshotName )
	{
		if( flushImage )
		{
			flushImage = false;
			glutPostRedisplay();
			return;
		}
		else
		{
			saveFrameBuffer( snapshotName , GL_FRONT );
			delete[] snapshotName;
			snapshotName = NULL;
		}
	}
	else if( strlen(videoHeader) && _currentFrame<_totalFrames )
	{
		if( !(_currentFrame%_frameStride) )
		{
			char snapshotName[1024];
			sprintf( snapshotName , "%s.%04d.jpg" , videoHeader , _currentFrame/_frameStride );
			printf( "Frame: %d\n" , _currentFrame );
			saveFrameBuffer( snapshotName , GL_FRONT );
		}
		_currentFrame++;
		if( _currentFrame==_totalFrames && _exitAfterVideo ) exit( 0 );
	}
	idle();
}
void Visualization::KeyboardFunc( unsigned char key , int x , int y )
{
	if( promptCallBack )
	{
		size_t len = strlen( promptString );
		if( key==KEY_BACK_SPACE )
		{
			if( len>promptLength ) promptString[len-1] = 0;
		}
		else if( key==KEY_ENTER )
		{
			promptCallBack( this , promptString+promptLength );
			promptString[0] = 0;
			promptLength = 0;
			promptCallBack = NULL;
		}
		else if( key==KEY_CTRL_C )
		{
			promptString[0] = 0;
			promptLength = 0;
			promptCallBack = NULL;
		}
		else if( key>=32 && key<=126 ) // ' ' to '~'
		{
			promptString[ len ] = key;
			promptString[ len+1 ] = 0;
		}
		glutPostRedisplay();
		return;
	}
	switch( key )
	{
	case KEY_CTRL_C:
		exit( 0 );
		break;
	default:
		for( int i=0 ; i<callBacks.size() ; i++ ) if( callBacks[i].key==key )
		{
			if( strlen( callBacks[i].prompt ) )
			{
				sprintf( promptString , "%s: " , callBacks[i].prompt );
				promptLength = int( strlen( promptString ) );
				promptCallBack = callBacks[i].callBackFunction;
			}
			else (*callBacks[i].callBackFunction)( this , NULL );
			break;
		}
	}
	keyboardFunc( key , x , y );
	glutPostRedisplay();
}

void Visualization::SpecialFunc( int key , int x , int y ){ specialFunc( key , x , y );}
void Visualization::Display( void )
{
	glClearColor( 1 , 1 , 1 , 1 );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	display();

	_lastFPSCount++;
	double t = Miscellany::Time();
	if( t-_lastFPSTime > _MIN_FPS_TIME )
	{
		_fps = (double)_lastFPSCount / (t-_lastFPSTime);
		_lastFPSCount = 0;
		_lastFPSTime = t;
	}
	if( showFPS ) writeRightString( 5 , screenHeight - fontHeight - 5 , "%d x %d @ %.2f" , screenWidth , screenHeight , _fps );

	int offset = fontHeight/2;
	if( showHelp )
	{
		{
			GLint vp[4];
			glGetIntegerv( GL_VIEWPORT , vp );

			glMatrixMode( GL_PROJECTION );
			glPushMatrix();
			glLoadIdentity();
			glOrtho( vp[0] , vp[2] , vp[1] , vp[3] , 0 , 1 );

			glMatrixMode( GL_MODELVIEW );
			glPushMatrix();
			glLoadIdentity();

			int x=0 , y = offset;
			for( int i=0 ; i<callBacks.size() ; i++ ) if( strlen( callBacks[i].description ) ) x = std::max< int >( x , stringWidth( "\'%c\': %s" , callBacks[i].key , callBacks[i].description ) ) , y += fontHeight + offset;

			GLint srcAlpha , dstAlpha;
			glGetIntegerv( GL_BLEND_SRC_ALPHA , &srcAlpha );
			glGetIntegerv( GL_BLEND_DST_ALPHA , &dstAlpha );

			glEnable( GL_BLEND );
			glBlendFunc( GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA );
			glBegin( GL_QUADS );
			glColor4f( 1.f , 1.f , 1.f , 0.5f );
			glVertex2i( screenWidth-5 , 0 ) , glVertex2i( screenWidth-(x+15) , 0 ) , glVertex2i( screenWidth-(x+15) , y ) , glVertex2i( screenWidth-5 , y );
			glEnd();
			glBlendFunc( srcAlpha , dstAlpha );
			glDisable( GL_BLEND );
			glDisable( GL_DEPTH_TEST );
			glLineWidth( 2.f );
			glBegin( GL_LINE_LOOP );
			glColor4f( 0.f , 0.f , 0.f , 1.f );
			glVertex2i( screenWidth-5 , 0 ) , glVertex2i( screenWidth-(x+15) , 0 ) , glVertex2i( screenWidth-(x+15) , y ) , glVertex2i( screenWidth-5 , y );
			glEnd();

			glMatrixMode( GL_PROJECTION );
			glPopMatrix();

			glMatrixMode( GL_MODELVIEW );
			glPopMatrix();
		}

		{
			int y = offset , width = 0;

			for( int i=0 ; i<callBacks.size() ; i++ ) if( strlen( callBacks[i].description ) )
				width = std::max< int >( width , stringWidth( "\'%c\': %s" , callBacks[i].key , callBacks[i].description ) );
			for( int i=0 ; i<callBacks.size() ; i++ ) if( strlen( callBacks[i].description ) )
				writeLeftString( screenWidth - 10 - width , y , "\'%c\': %s" , callBacks[i].key , callBacks[i].description ) , y += fontHeight + offset;
		}
	}
#if 1
	if( showInfo && info.size() )
	{
		{
			GLint vp[4];
			glGetIntegerv( GL_VIEWPORT , vp );

			glMatrixMode( GL_MODELVIEW );
			glPushMatrix();
			glLoadIdentity();
			glMatrixMode( GL_PROJECTION );
			glPushMatrix();
			glLoadIdentity();
			glOrtho( vp[0] , vp[2] , vp[1] , vp[3] , 0 , 1 );
			int x=0 , y = offset;
			for( int i=0 ; i<info.size() ; i++ ) if( strlen( info[i] ) ) x = std::max< int >( x , glutBitmapLength( font , (unsigned char*) info[i] ) ) , y += fontHeight + offset;
			glEnable( GL_BLEND );
			GLint srcAlpha , dstAlpha;
			glGetIntegerv( GL_BLEND_SRC_ALPHA , &srcAlpha );
			glGetIntegerv( GL_BLEND_DST_ALPHA , &dstAlpha );
			glBlendFunc( GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA );
			glBegin( GL_QUADS );
			glColor4f( 1.f , 1.f , 1.f , 0.5f );
			glVertex2i( 5 , 0 ) , glVertex2i( x+15 , 0 ) , glVertex2i( x+15 , y ) , glVertex2i( 5 , y );
			glEnd();
			glBlendFunc( srcAlpha , dstAlpha );
			glDisable( GL_BLEND );
			glDisable( GL_DEPTH_TEST );
			glLineWidth( 2.f );
			glBegin( GL_LINE_LOOP );
			glColor4f( 0.f , 0.f , 0.f , 1.f );
			glVertex2i( 5 , 0 ) , glVertex2i( x+15 , 0 ) , glVertex2i( x+15 , y ) , glVertex2i( 5 , y );
			glEnd();

			glMatrixMode( GL_PROJECTION );
			glPopMatrix();

			glMatrixMode( GL_MODELVIEW );
			glPopMatrix();
		}
		{
			int y = offset;
			for( int i=0 ; i<info.size() ; i++ ) if( strlen( info[i] ) )
				writeLeftString( 10 , y , "%s" , info[i] ) , y += fontHeight + offset;
		}
	}
#else
	if( showInfo )
	{
		int y = offset;
		for( int i=0 ; i<info.size() ; i++ ) if( strlen( info[i] ) )
			writeLeftString( 10 , y , "%s" , info[i] ) , y += fontHeight + offset;
	}
#endif
	if( strlen( promptString ) )
	{
		void* _font = font;
		int _fontHeight = fontHeight;
		font = promptFont;
		fontHeight = promptFontHeight;

		int sw = StringWidth ( font , promptString );
		glColor4f( 1.f , 1.f , 1.f , 0.5 );
		glEnable( GL_BLEND );
		glBlendFunc( GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA );
		glBegin( GL_QUADS );
		{
			glVertex2f(     0 , screenHeight              );
			glVertex2f( sw+20 , screenHeight              );
			glVertex2f( sw+20 , screenHeight-fontHeight*2 );
			glVertex2f(     0 , screenHeight-fontHeight*2 );
		}
		glEnd();
		glDisable( GL_BLEND );
		glColor4f( 0 , 0 , 0 , 1 );
		glLineWidth( 2.f );
		glBegin( GL_LINE_LOOP );
		{
			glVertex2f(     0 , screenHeight              );
			glVertex2f( sw+20 , screenHeight              );
			glVertex2f( sw+20 , screenHeight-fontHeight*2 );
			glVertex2f(     0 , screenHeight-fontHeight*2 );
		}
		glEnd();
		double tm = Miscellany::Time();
		tm -= floor( tm );
		static char _promptString[1024];
		if( tm<0.5 ) sprintf( _promptString , "%s " , promptString );
		else         sprintf( _promptString , "%s_" , promptString );
		writeLeftString( 10 , screenHeight-fontHeight-fontHeight/2 , _promptString );
		font = _font;
		fontHeight = _fontHeight;
	}
	glutSwapBuffers();
}

void Visualization::WriteLeftString( int x , int y , void* font , const char* format , ... )
{
	static char str[1024];
	{
		va_list args;
		va_start( args , format );
		vsprintf( str , format , args );
		va_end( args );
	}

	GLint vp[4];

	glGetIntegerv( GL_VIEWPORT , vp );
	glPushMatrix();
	glLoadIdentity();
	glMatrixMode( GL_PROJECTION );
	glPushMatrix();
	glLoadIdentity();
	glOrtho( vp[0] , vp[2] , vp[1] , vp[3] , 0 , 1 );

	glDisable( GL_DEPTH_TEST );
	glDisable( GL_LIGHTING );
	glColor4f( 0 , 0 , 0 , 1 );
	glRasterPos2f( x , y  );
	int len = int( strlen( str ) );
	for( int i=0 ; i<len ; i++ ) glutBitmapCharacter( font , str[i] );
	glPopMatrix();
	glMatrixMode( GL_MODELVIEW );
	glPopMatrix();
}
int Visualization::StringWidth( void* font , const char* format , ... )
{
	static char str[1024];
	{
		va_list args;
		va_start( args , format );
		vsprintf( str , format , args );
		va_end( args );
	}
	return glutBitmapLength( font , (unsigned char*) str );
}
int Visualization::stringWidth( const char* format , ... ) const
{
	static char str[1024];
	{
		va_list args;
		va_start( args , format );
		vsprintf( str , format , args );
		va_end( args );
	}
	return glutBitmapLength( font , (unsigned char*) str );
}
void Visualization::writeLeftString( int x , int y , const char* format , ... ) const
{
	static char str[1024];
	{
		va_list args;
		va_start( args , format );
		vsprintf( str , format , args );
		va_end( args );
	}
	WriteLeftString( x , y , font , str );
}
void Visualization::writeRightString( int x , int y , const char* format , ... ) const
{
	static char str[1024];
	{
		va_list args;
		va_start( args , format );
		vsprintf( str , format , args );
		va_end( args );
	}
	WriteLeftString( screenWidth-x-glutBitmapLength( font , (unsigned char*) str ) , y , font  ,str );
}
void Visualization::saveFrameBuffer( const char* fileName , int whichBuffer )
{
	Pointer( float ) pixels = AllocPointer< float >( sizeof(float) * 3 * screenWidth * screenHeight );
	Pointer( unsigned char ) _pixels = AllocPointer< unsigned char >( sizeof(unsigned char) * 3 * screenWidth * screenHeight );
	glReadBuffer( whichBuffer );
	glReadPixels( 0 , 0 , screenWidth , screenHeight , GL_RGB , GL_FLOAT , pixels );
	for( int j=0 ; j<screenHeight ; j++ ) for( int i=0 ; i<screenWidth ; i++ ) for( int c=0 ; c<3 ; c++ )
	{
		int ii = int( pixels[ c + i * 3 + ( screenHeight - 1 - j ) * screenWidth * 3 ]*256 );
		if( ii<  0 ) ii =   0;
		if( ii>255 ) ii = 255;
		_pixels[ c + i * 3 + j * screenWidth * 3 ] = (unsigned char)ii;
	}
	FreePointer( pixels );
	if( !ImageWriter::Write( fileName , _pixels , screenWidth , screenHeight , 3 ) ) Miscellany::Warn( "Failed to write image: %s" , fileName );
	FreePointer( _pixels );
}

#endif // VISUALIZATION_INCLUDED