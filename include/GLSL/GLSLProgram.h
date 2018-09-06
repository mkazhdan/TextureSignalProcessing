//Code from : OpenGL 4 Shading Language Cookbook Second Edition, by David Wolff
// Modified to only support simple vertex/fragment shaders
#ifndef GLSLPROGRAM_H
#define GLSLPROGRAM_H

#include <GL/glew.h>
#include <string>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <iostream>

#pragma warning( disable : 4290 )

class GLSLProgramException : public std::runtime_error
{
public:
	GLSLProgramException( const std::string & msg ) : std::runtime_error( msg ) { }
};

namespace GLSLShader
{
	enum GLSLShaderType
	{
		VERTEX          = GL_VERTEX_SHADER ,
		FRAGMENT        = GL_FRAGMENT_SHADER ,
		GEOMETRY        = GL_GEOMETRY_SHADER ,
		TESS_CONTROL    = GL_TESS_CONTROL_SHADER ,
		TESS_EVALUATION = GL_TESS_EVALUATION_SHADER ,
		COMPUTE         = GL_COMPUTE_SHADER
	};
};

class GLSLProgram
{
	template< unsigned int Dim > static void glUniformiv      ( GLint location , GLsizei count ,                       const GLint*   value );
	template< unsigned int Dim > static void glUniformfv      ( GLint location , GLsizei count ,                       const GLfloat* value );
	template< unsigned int Dim > static void glUniformMatrixfv( GLint location , GLsizei count , GLboolean transpose , const GLfloat* value );
public:
	int  handle;
	bool linked;
	std::string vertex_shader_src;
	std::string fragment_shader_src;

	static bool fileExists( const std::string & fileName );
	static std::string getExtension( const char * fileName );
	static GLSLShader::GLSLShaderType getShaderType( const char* fileName ) throw( GLSLProgramException );

	// Make these private in order to make the object non-copyable
	GLSLProgram( const GLSLProgram & other ) { }
	GLSLProgram & operator=( const GLSLProgram &other ){ return *this; }

	GLSLProgram( void );
	GLSLProgram( const char* vs_filename , const char* fs_filename ) throw( GLSLProgramException );
	GLSLProgram( const std::string& vs_src , const std::string& fs_src );

	~GLSLProgram( void );

	void   compileShader( const char *fileName )                                                                     throw( GLSLProgramException );
	void   compileShader( const char * fileName      , GLSLShader::GLSLShaderType type )                             throw( GLSLProgramException );
	void   compileShader( const std::string & source , GLSLShader::GLSLShaderType type , const char *fileName=NULL ) throw( GLSLProgramException );

	void   link( void )     throw( GLSLProgramException );
	void   validate( void ) throw( GLSLProgramException );
	void   use( void )      throw( GLSLProgramException );

	int    getHandle( void );
	bool   isLinked( void );

	void   bindAttribLocation( GLuint location , const char * name );
	void   bindFragDataLocation( GLuint location , const char * name );

	void setUniform( const char *name , int    val , bool showWarning=true );
	void setUniform( const char *name , float  val , bool showWarning=true );
	void setUniform( const char *name , double val , bool showWarning=true );
	template< unsigned int > void setUniform( const char* name , const int*    v , bool showWarning=true );
	template< unsigned int > void setUniform( const char* name , const float*  v , bool showWarning=true );
	template< unsigned int > void setUniform( const char* name , const double* v , bool showWarning=true );
	template< unsigned int > void setUniformMatrix( const char* name , const float*  m , bool showWarning=true );
	template< unsigned int > void setUniformMatrix( const char* name , const double* m , bool showWarning=true );

	void   printActiveUniforms( void );
	void   printActiveUniformBlocks( void );
	void   printActiveAttribs( void );

	const char * getTypeString( GLenum type );
	void setup( void );
};
template< unsigned int Dim > void GLSLProgram::glUniformiv( GLint location  , GLsizei count , const GLint* value ){ fprintf( stderr , "[ERROR] glUniform[%d]iv undefined\n" , Dim ) , exit( 0 ); }
template<> void GLSLProgram::glUniformiv< 1 >( GLint location  , GLsizei count , const GLint* value ){ glUniform1iv( location , count , value ); }
template<> void GLSLProgram::glUniformiv< 2 >( GLint location  , GLsizei count , const GLint* value ){ glUniform2iv( location , count , value ); }
template<> void GLSLProgram::glUniformiv< 3 >( GLint location  , GLsizei count , const GLint* value ){ glUniform3iv( location , count , value ); }
template<> void GLSLProgram::glUniformiv< 4 >( GLint location  , GLsizei count , const GLint* value ){ glUniform4iv( location , count , value ); }
template< unsigned int Dim > void GLSLProgram::glUniformfv( GLint location  , GLsizei count , const GLfloat* value ){ fprintf( stderr , "[ERROR] glUniform[%d]fv undefined\n" , Dim ) , exit( 0 ); }
template<> void GLSLProgram::glUniformfv< 1 >( GLint location  , GLsizei count , const GLfloat* value ){ glUniform1fv( location , count , value ); }
template<> void GLSLProgram::glUniformfv< 2 >( GLint location  , GLsizei count , const GLfloat* value ){ glUniform2fv( location , count , value ); }
template<> void GLSLProgram::glUniformfv< 3 >( GLint location  , GLsizei count , const GLfloat* value ){ glUniform3fv( location , count , value ); }
template<> void GLSLProgram::glUniformfv< 4 >( GLint location  , GLsizei count , const GLfloat* value ){ glUniform4fv( location , count , value ); }
template< unsigned int Dim > void GLSLProgram::glUniformMatrixfv( GLint location  , GLsizei count , GLboolean transpose , const GLfloat* value ){ fprintf( stderr , "[ERROR] glUniformMatrix[%d]fv undefined\n" , Dim ) , exit( 0 ); }
template<> void GLSLProgram::glUniformMatrixfv< 2 >( GLint location  , GLsizei count , GLboolean transpose , const GLfloat* value ){ glUniformMatrix2fv( location , count , transpose , value ); }
template<> void GLSLProgram::glUniformMatrixfv< 3 >( GLint location  , GLsizei count , GLboolean transpose , const GLfloat* value ){ glUniformMatrix3fv( location , count , transpose , value ); }
template<> void GLSLProgram::glUniformMatrixfv< 4 >( GLint location  , GLsizei count , GLboolean transpose , const GLfloat* value ){ glUniformMatrix4fv( location , count , transpose , value ); }

namespace GLSLShaderInfo
{
	struct shader_file_extension
	{
		const char *ext;
		GLSLShader::GLSLShaderType type;
	};

	struct shader_file_extension extensions[] =
	{
		{ ".vs"   , GLSLShader::VERTEX },
		{ ".vert" , GLSLShader::VERTEX },
		{ ".gs"   , GLSLShader::GEOMETRY },
		{ ".geom" , GLSLShader::GEOMETRY },
		{ ".tcs"  , GLSLShader::TESS_CONTROL },
		{ ".tes"  , GLSLShader::TESS_EVALUATION },
		{ ".fs"   , GLSLShader::FRAGMENT },
		{ ".frag" , GLSLShader::FRAGMENT },
		{ ".cs"   , GLSLShader::COMPUTE }
	};
}

GLSLProgram::GLSLProgram( void ) : handle( 0 ) , linked( false ) {}
GLSLProgram::GLSLProgram( const char* vs_filename , const char* fs_filename ) throw( GLSLProgramException ) : handle( 0 ) , linked( false )
{
	// Check that the first shader is a vertex shader and that the file exists. If it does, read in the code.
	if( getShaderType( vs_filename )!=GLSLShader::VERTEX ) throw GLSLProgramException( std::string( "Expected vertex shader: " ) + std::string( vs_filename ) );
	if( !fileExists( vs_filename ) ) throw GLSLProgramException( std::string( "Vertex shader: " ) + std::string( vs_filename ) + " not found." );
	{
		std::ifstream inFile( vs_filename , std::ios::in );
		if( !inFile ) throw GLSLProgramException( std::string( "Unable to open vertex shader: " ) + std::string( vs_filename ) );
		std::stringstream code;
		code << inFile.rdbuf();
		inFile.close();
		vertex_shader_src = code.str();
	}

	// Check that the second shader is a fragment shader and that the file exists. If it does, read in the code.
	if( getShaderType( fs_filename )!=GLSLShader::FRAGMENT ) throw GLSLProgramException( std::string( "Expected fragment shader: " ) + std::string( fs_filename ) );
	if( !fileExists( fs_filename ) ) throw GLSLProgramException( std::string( "Fragment shader: " ) + std::string( fs_filename ) + " not found." );
	{
		std::ifstream inFile( fs_filename , std::ios::in );
		if( !inFile ) throw GLSLProgramException( std::string( "Unable to open fragment shader: " ) + std::string( fs_filename ) );
		std::stringstream code;
		code << inFile.rdbuf();
		inFile.close();
		fragment_shader_src = code.str();
	}
}
GLSLProgram::GLSLProgram( const std::string& vs_src , const std::string& fs_src ) : vertex_shader_src( vs_src ) , fragment_shader_src( fs_src ) , handle( 0 ) , linked( false ) {}

GLSLProgram::~GLSLProgram( void )
{
	if( handle==0 ) return;

	// Query the number of attached shaders
	GLint numShaders = 0;
	glGetProgramiv( handle , GL_ATTACHED_SHADERS , &numShaders );

	// Get the shader names
	GLuint * shaderNames = new GLuint[numShaders];
	glGetAttachedShaders( handle, numShaders , NULL , shaderNames );

	// Delete the shaders
	for( int i=0 ; i<numShaders ; i++ ) glDeleteShader( shaderNames[i] );

	// Delete the program
	glDeleteProgram( handle );

	delete[] shaderNames;
}

GLSLShader::GLSLShaderType GLSLProgram::getShaderType( const char* fileName ) throw( GLSLProgramException )
{
	int numExts = sizeof( GLSLShaderInfo::extensions ) / sizeof( GLSLShaderInfo::shader_file_extension );

	std::string ext = getExtension( fileName );
	bool matchFound = false;
	for( int i=0 ; i<numExts ; i++ ) if( ext==GLSLShaderInfo::extensions[i].ext ) return GLSLShaderInfo::extensions[i].type;
	throw GLSLProgramException( std::string( "Unrecognized extension: " ) + ext );
}

std::string GLSLProgram::getExtension( const char * name )
{
	std::string nameStr( name );
	
	size_t loc = nameStr.find_last_of('.');
	if( loc!=std::string::npos ) return nameStr.substr( loc , std::string::npos );
	return "";
}

void GLSLProgram::compileShader( const std::string & source , GLSLShader::GLSLShaderType type , const char * fileName ) throw( GLSLProgramException )
{
	if( handle<=0 )
	{
		handle = glCreateProgram();
		if( handle==0 ) throw GLSLProgramException( "Unable to create shader program." );
	}

	GLuint shaderHandle = glCreateShader( type );

	const char * c_code = source.c_str();
	glShaderSource( shaderHandle , 1 , &c_code , NULL );

	// Compile the shader
	glCompileShader( shaderHandle );

	// Check for errors
	int result;
	glGetShaderiv( shaderHandle , GL_COMPILE_STATUS , &result );
	if( GL_FALSE==result )
	{
		// Compile failed, get log
		int length = 0;
		std::string logString;
		glGetShaderiv( shaderHandle , GL_INFO_LOG_LENGTH , &length );
		if( length>0 )
		{
			char * c_log = new char[length];
			int written = 0;
			glGetShaderInfoLog( shaderHandle , length , &written , c_log );
			logString = c_log;
			delete[] c_log;
		}
		std::string msg;
		if( fileName ) msg = std::string( fileName ) + ": shader compliation failed\n";
		else           msg = "Shader compilation failed.\n";
		msg += logString;

		throw GLSLProgramException( msg );
	}
	else glAttachShader( handle , shaderHandle );
}

void GLSLProgram::link( void ) throw( GLSLProgramException )
{
	if( linked ) return;
	if( handle<=0 ) throw GLSLProgramException( "Program has not been compiled." );

	glLinkProgram( handle );

	int status = 0;
	glGetProgramiv( handle , GL_LINK_STATUS , &status );
	if( GL_FALSE==status )
	{
		// Store log and return false
		int length = 0;
		std::string logString;

		glGetProgramiv( handle , GL_INFO_LOG_LENGTH , &length );

		if( length>0 )
		{
			char * c_log = new char[length];
			int written = 0;
			glGetProgramInfoLog( handle , length , &written , c_log );
			logString = c_log;
			delete[] c_log;
		}
		throw GLSLProgramException( std::string( "Program link failed:\n" ) + logString );
	}
	else linked = true;
}

void GLSLProgram::use( void ) throw( GLSLProgramException )
{
	if( handle<=0 || (!linked) ) throw GLSLProgramException( "Shader has not been linked" );
	glUseProgram( handle );
}

int GLSLProgram::getHandle( void ){ return handle; }
bool GLSLProgram::isLinked( void ){ return linked; }

void GLSLProgram::bindAttribLocation( GLuint location , const char * name ){ glBindAttribLocation(handle, location, name); }
void GLSLProgram::bindFragDataLocation( GLuint location , const char * name ){ glBindFragDataLocation( handle , location , name ); }
/////////////////////////////
// GLSLProgram::setUniform //
/////////////////////////////
void GLSLProgram::setUniform( const char *name , int val , bool showWarning )
{
	GLint loc = glGetUniformLocation( handle , name );
	if( loc>=0 ) glUniform1i( loc , val );
	else if( showWarning ) fprintf( stderr , "[WARNING] Non-existant uniform: %s\n" , name );
}
void GLSLProgram::setUniform( const char *name , float val , bool showWarning )
{
	GLint loc = glGetUniformLocation( handle , name );
	if( loc>=0 ) glUniform1f( loc , val );
	else if( showWarning ) fprintf( stderr , "[WARNING] Non-existant uniform: %s\n" , name );
}
void GLSLProgram::setUniform( const char *name , double val , bool showWarning )
{
	GLint loc = glGetUniformLocation( handle , name );
	if( loc>=0 ) glUniform1f( loc , val );
	else if( showWarning ) fprintf( stderr , "[WARNING] Non-existant uniform: %s\n" , name );
}
template< unsigned int Dim >
void GLSLProgram::setUniform( const char* name , const int* v , bool showWarning )
{
	GLint loc = glGetUniformLocation( handle , name );
	if( loc>=0 ) glUniformiv< Dim >( loc , 1 , v );
	else if( showWarning ) fprintf( stderr , "[WARNING] Non-existant uniform: %s\n" , name );
}
template< unsigned int Dim >
void GLSLProgram::setUniform( const char* name , const float* v , bool showWarning )
{
	GLint loc = glGetUniformLocation( handle , name );
	if( loc>=0 ) glUniformfv< Dim >( loc , 1 , v );
	else if( showWarning ) fprintf( stderr , "[WARNING] Non-existant uniform: %s\n" , name );
}
template< unsigned int Dim >
void GLSLProgram::setUniform( const char* name , const double* v , bool showWarning )
{
	float _v[Dim];
	for( int i=0 ; i<Dim ; i++ ) _v[i] = (float)v[i];
	setUniform< Dim >( name , _v , showWarning );
}
template< unsigned int Dim >
void GLSLProgram::setUniformMatrix( const char* name , const float* m , bool showWarning )
{
	GLint loc = glGetUniformLocation( handle , name );
	if( loc>=0 ) glUniformMatrixfv< Dim >( loc , 1 , GL_FALSE , m );
	else if( showWarning ) fprintf( stderr , "[WARNING] Non-existant uniform: %s\n" , name );
};
template< unsigned int Dim >
void GLSLProgram::setUniformMatrix( const char* name , const double* m , bool showWarning )
{
	float _m[Dim*Dim];
	for( int i=0 ; i<Dim*Dim ; i++ ) _m[i] = (float)m[i];
	setUniformMatrix< Dim >( name , _m , showWarning );
}



void GLSLProgram::validate( void ) throw( GLSLProgramException )
{
	if( !isLinked() ) throw GLSLProgramException("Program is not linked");

	GLint status;
	glValidateProgram( handle );
	glGetProgramiv( handle , GL_VALIDATE_STATUS , &status );

	if( GL_FALSE==status )
	{
		// Store log and return false
		int length = 0;
		std::string logString;

		glGetProgramiv( handle , GL_INFO_LOG_LENGTH , &length );

		if( length>0 )
		{
			char * c_log = new char[length];
			int written = 0;
			glGetProgramInfoLog( handle , length , &written , c_log );
			logString = c_log;
			delete[] c_log;
		}
		throw GLSLProgramException( std::string( "Program failed to validate\n" ) + logString );
	}
}

bool GLSLProgram::fileExists( const std::string & fileName )
{
	struct stat info;
	int ret = -1;

	ret = stat(fileName.c_str(), &info);
	return 0 == ret;
}

void GLSLProgram::setup( void )
{
	try
	{
		compileShader(   vertex_shader_src , GLSLShader::VERTEX   );
		compileShader( fragment_shader_src , GLSLShader::FRAGMENT );
		link();
		validate();
	}
	catch( GLSLProgramException &e )
	{
		std::cerr << e.what() << std::endl;
		exit( EXIT_FAILURE );
	}
}

void GLSLProgram::printActiveUniforms( void )
{
	GLint numUniforms = 0;
	glGetProgramInterfaceiv(handle, GL_UNIFORM, GL_ACTIVE_RESOURCES, &numUniforms);

	GLenum properties[] = { GL_NAME_LENGTH, GL_TYPE, GL_LOCATION, GL_BLOCK_INDEX };

	printf("Active uniforms:\n");
	for (int i = 0; i < numUniforms; ++i) {
		GLint results[4];
		glGetProgramResourceiv(handle, GL_UNIFORM, i, 4, properties, 4, NULL, results);

		if (results[3] != -1) continue;  // Skip uniforms in blocks 
		GLint nameBufSize = results[0] + 1;
		char * name = new char[nameBufSize];
		glGetProgramResourceName(handle, GL_UNIFORM, i, nameBufSize, NULL, name);
		printf("%-5d %s (%s)\n", results[2], name, getTypeString(results[1]));
		delete[] name;
	}
}

void GLSLProgram::printActiveUniformBlocks( void )
{
	GLint numBlocks = 0;

	glGetProgramInterfaceiv(handle, GL_UNIFORM_BLOCK, GL_ACTIVE_RESOURCES, &numBlocks);
	GLenum blockProps[] = { GL_NUM_ACTIVE_VARIABLES, GL_NAME_LENGTH };
	GLenum blockIndex[] = { GL_ACTIVE_VARIABLES };
	GLenum props[] = { GL_NAME_LENGTH, GL_TYPE, GL_BLOCK_INDEX };

	for (int block = 0; block < numBlocks; ++block) {
		GLint blockInfo[2];
		glGetProgramResourceiv(handle, GL_UNIFORM_BLOCK, block, 2, blockProps, 2, NULL, blockInfo);
		GLint numUnis = blockInfo[0];

		char * blockName = new char[blockInfo[1] + 1];
		glGetProgramResourceName(handle, GL_UNIFORM_BLOCK, block, blockInfo[1] + 1, NULL, blockName);
		printf("Uniform block \"%s\":\n", blockName);
		delete[] blockName;

		GLint * unifIndexes = new GLint[numUnis];
		glGetProgramResourceiv(handle, GL_UNIFORM_BLOCK, block, 1, blockIndex, numUnis, NULL, unifIndexes);

		for (int unif = 0; unif < numUnis; ++unif) {
			GLint uniIndex = unifIndexes[unif];
			GLint results[3];
			glGetProgramResourceiv(handle, GL_UNIFORM, uniIndex, 3, props, 3, NULL, results);

			GLint nameBufSize = results[0] + 1;
			char * name = new char[nameBufSize];
			glGetProgramResourceName(handle, GL_UNIFORM, uniIndex, nameBufSize, NULL, name);
			printf("    %s (%s)\n", name, getTypeString(results[1]));
			delete[] name;
		}

		delete[] unifIndexes;
	}
}

void GLSLProgram::printActiveAttribs( void )
{
	GLint numAttribs;
	glGetProgramInterfaceiv(handle, GL_PROGRAM_INPUT, GL_ACTIVE_RESOURCES, &numAttribs);

	GLenum properties[] = { GL_NAME_LENGTH, GL_TYPE, GL_LOCATION };

	printf("Active attributes:\n");
	for (int i = 0; i < numAttribs; ++i) {
		GLint results[3];
		glGetProgramResourceiv(handle, GL_PROGRAM_INPUT, i, 3, properties, 3, NULL, results);

		GLint nameBufSize = results[0] + 1;
		char * name = new char[nameBufSize];
		glGetProgramResourceName(handle, GL_PROGRAM_INPUT, i, nameBufSize, NULL, name);
		printf("%-5d %s (%s)\n", results[2], name, getTypeString(results[1]));
		delete[] name;
	}
}

const char * GLSLProgram::getTypeString(GLenum type) {
	// There are many more types than are covered here, but
	// these are the most common in these examples.
	switch( type )
	{
	case GL_FLOAT:        return "float";
	case GL_FLOAT_VEC2:   return "vec2";
	case GL_FLOAT_VEC3:   return "vec3";
	case GL_FLOAT_VEC4:   return "vec4";
	case GL_DOUBLE:       return "double";
	case GL_INT:          return "int";
	case GL_INT_VEC2:     return "ivec2";
	case GL_INT_VEC3:     return "ivec3";
	case GL_UNSIGNED_INT: return "unsigned int";
	case GL_BOOL:         return "bool";
	case GL_FLOAT_MAT2:   return "mat2";
	case GL_FLOAT_MAT3:   return "mat3";
	case GL_FLOAT_MAT4:   return "mat4";
	default:              return "?";
	}
}


#endif // GLSLPROGRAM_H