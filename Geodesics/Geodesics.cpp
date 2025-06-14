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

#include <Src/PreProcessing.h>

#include <Misha/CmdLineParser.h> 
#include <Misha/Miscellany.h>
#include <Misha/Exceptions.h>
#include <Misha/FEM.h>
#include <Misha/MultiThreading.h>
#include <Src/Hierarchy.h>
#include <Src/SimpleTriangleMesh.h>
#include <Src/Basis.h>
#include <Src/Solver.h>
#include <Src/QuadratureIntegration.h>
#include <Src/Operators.h>
#include <Src/Padding.h>
#include <Src/TexturedMeshVisualization.h>

using namespace MishaK;
using namespace MishaK::TSP;

CmdLineParameter< std::string >
	Input( "in" ) ,
	CameraConfig( "camera" );

CmdLineParameter< unsigned int >
	Width( "width" , 1024 ) ,
	Height( "height" , 1024 ) ,
	Levels( "levels" , 4 ) ,
	DisplayMode( "display" , TWO_REGION_DISPLAY ) ,
	MatrixQuadrature( "mQuadrature" , 6 ) ,
	VectorFieldQuadrature( "vfQuadrature" , 6 ) ,
	MultigridBlockHeight ( "mBlockH" ,  16 ) ,
	MultigridBlockWidth  ( "mBlockW" , 128 ) ,
	MultigridPaddedHeight( "mPadH"   ,   0 ) , 
	MultigridPaddedWidth ( "mPadW"   ,   2 ) ,
	RandomJitter( "jitter" , 0 );

CmdLineParameter< double >
	DiffusionInterpolationWeight( "interpolation" , 1e3 ) ,
	CollapseEpsilon( "collapse" , 0 );

CmdLineReadable
	Serial( "serial" ) , 
	Verbose( "verbose" ) ,
	NoHelp( "noHelp" ) ,
	DetailVerbose( "detail" ) ,
	UseDirectSolver( "useDirectSolver" ) ,
	Double( "double" ) ,
	Nearest( "nearest" ) ,
	PreciseIntegration( "preciseIntegration" );


CmdLineReadable* params[] =
{
	&Input , &Width , &Height , &DiffusionInterpolationWeight , &CameraConfig , &Levels , &UseDirectSolver , &Serial , &DisplayMode , &MultigridBlockHeight , &MultigridBlockWidth , &MultigridPaddedHeight , &MultigridPaddedWidth ,
	&Verbose , &DetailVerbose ,
	&RandomJitter ,
	&Double ,
	&MatrixQuadrature , &VectorFieldQuadrature ,
	&PreciseIntegration ,
	&NoHelp ,
	&Nearest ,
	&CollapseEpsilon ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n", ex );
	printf( "\t --%s <input mesh>\n" , Input.name.c_str() );
	printf( "\t[--%s <texture width>=%d]\n" , Width.name.c_str() , Width.value );
	printf( "\t[--%s <texture height>=%d]\n" , Height.name.c_str() , Height.value );
	printf( "\t[--%s <diffusion interpolation weight>=%f]\n" , DiffusionInterpolationWeight.name.c_str() , DiffusionInterpolationWeight.value );
	printf( "\t[--%s <system matrix quadrature points per triangle>=%d]\n" , MatrixQuadrature.name.c_str() , MatrixQuadrature.value );
	printf( "\t[--%s <normalized vector field quadrature points per triangle>=%d]\n" , VectorFieldQuadrature.name.c_str() , VectorFieldQuadrature.value );
	printf( "\t[--%s]\n" , PreciseIntegration.name.c_str() );
	printf( "\t[--%s]\n" , UseDirectSolver.name.c_str() );
	printf( "\t[--%s <jittering seed>]\n" , RandomJitter.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );

	printf( "\t[--%s <camera configuration file>]\n" , CameraConfig.name.c_str() );
	printf( "\t[--%s <hierarchy levels>=%d]\n" , Levels.name.c_str() , static_cast< unsigned int >(Levels.value) );
	printf( "\t[--%s]\n" , DetailVerbose.name.c_str() );
	printf( "\t[--%s <display mode>=%d]\n" , DisplayMode.name.c_str() , DisplayMode.value );
	printf( "\t\t%d] One Region \n" , ONE_REGION_DISPLAY );
	printf( "\t\t%d] Two Region \n" , TWO_REGION_DISPLAY );
	printf( "\t[--%s <multigrid block width>=%d]\n"   , MultigridBlockWidth.name.c_str()   , MultigridBlockWidth.value   );
	printf( "\t[--%s <multigrid block height>=%d]\n"  , MultigridBlockHeight.name.c_str()  , MultigridBlockHeight.value  );
	printf( "\t[--%s <multigrid padded width>=%d]\n"  , MultigridPaddedWidth.name.c_str()  , MultigridPaddedWidth.value  );
	printf( "\t[--%s <multigrid padded height>=%d]\n" , MultigridPaddedHeight.name.c_str() , MultigridPaddedHeight.value );
	printf( "\t[--%s <collapse epsilon>=%g]\n" , CollapseEpsilon.name.c_str() , CollapseEpsilon.value );
	printf( "\t[--%s]\n" , Nearest.name.c_str() );
	printf( "\t[--%s]\n" , Serial.name.c_str() );
	printf( "\t[--%s]\n" , NoHelp.name.c_str() );
}

template< typename PreReal , typename Real >
class Geodesics
{
public:
	static TexturedTriangleMesh< PreReal > mesh;
	static int textureWidth;
	static int textureHeight;
	static Real diffusionInterpolationWeight;
	static Real geodesicInterpolationWeight;
	static unsigned int levels;
	static int updateCount;
	
	static int steps;
	static char stepsString[];

	static Padding padding;
	
	static std::vector< Point3D< float > > textureNodePositions;

	static HierarchicalSystem< PreReal , Real > hierarchy;

	static std::vector< TextureNodeInfo< PreReal > > textureNodes;
	static Image< int > nodeIndex;

	static SparseMatrix< Real , int > mass;
	static SparseMatrix< Real , int > stiffness;

	static SparseMatrix< Real , int > smoothImpulseMatrix;
	static SparseMatrix< Real , int > geodesicDistanceMatrix;

	static unsigned int impulseTexel;

	static ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< PreReal , 2 > > > parameterMetric;

	static Real smoothImpulseRange;
	static Real geodesicDistanceRange;

	//Impulse Smoothing
	static std::vector< SystemCoefficients< Real > > multigridSmoothImpulseCoefficients;
	static std::vector< MultigridLevelVariables< Real > > multigridSmoothImpulseVariables;
	
	//Geodesic Distance
	static std::vector< SystemCoefficients< Real > > multigridGeodesicDistanceCoefficients;
	static std::vector< MultigridLevelVariables< Real > > multigridGeodesicDistanceVariables;

#if defined( USE_CHOLMOD )
	typedef CholmodCholeskySolver< Real , 1 > DirectSolver;
#elif defined( USE_EIGEN )
	typedef EigenCholeskySolver< Real , 1 > DirectSolver;
#elif defined( USE_EIGEN_PARDISO )
	typedef EigenPardisoSolver< Real , 1 > DirectSolver;
#else
#error "[ERROR] No solver defined!"
#endif

	static VCycleSolvers< DirectSolver > smoothImpulseSolvers;
	static VCycleSolvers< DirectSolver > geodesicDistanceSolvers;
	static DirectSolver fineSmoothImpulseSolver;
	static DirectSolver fineGeodesicDistanceSolver;

	static std::vector<MultigridLevelIndices<Real>> multigridIndices;

	static GradientIntegrator< Real > gradientIntegrator;
	static typename GradientIntegrator< Real >::template Scratch< Real , Real > gradientIntegratorScratch;

	static unsigned char * outputBuffer;

	//Visulization
	static TexturedMeshVisualization visualization;
	static int mouseX, mouseY;
	static bool mouseSelectionActive;

	static void UpdateOutputBuffer( const std::vector< Real >& solution );

	static void ToggleUpdateCallBack( Visualization* v , const char* prompt );
	static void IncrementUpdateCallBack( Visualization* v , const char* prompt );
	static void ExportTextureCallBack( Visualization* v , const char* prompt );
	static void Init( void );
	static void InitializeVisualization( int width , int height );
	static void ComputeExactSolution( void);
	static void UpdateSolution( void );
	static void InitializeSystem( int width , int height );

	static void Display( void ){ visualization.Display(); }
	static void MouseFunc( int button , int state , int x , int y );
	static void MotionFunc( int x , int y );
	static void Reshape( int w , int h ){ visualization.Reshape(w,h); }
	static void KeyboardFunc( unsigned char key , int x , int y ){ visualization.KeyboardFunc( key , x , y ); }
	static void Idle( void );
};

template< typename PreReal , typename Real > TexturedTriangleMesh< PreReal >								Geodesics< PreReal , Real >::mesh;
template< typename PreReal , typename Real > int															Geodesics< PreReal , Real >::textureWidth;
template< typename PreReal , typename Real > int															Geodesics< PreReal , Real >::textureHeight;
template< typename PreReal , typename Real > TexturedMeshVisualization										Geodesics< PreReal , Real >::visualization;
template< typename PreReal , typename Real > int															Geodesics< PreReal , Real >::mouseX = -1;
template< typename PreReal , typename Real > int															Geodesics< PreReal , Real >::mouseY = -1;
template< typename PreReal , typename Real > bool															Geodesics< PreReal , Real >::mouseSelectionActive = false;
template< typename PreReal , typename Real > Padding														Geodesics< PreReal , Real >::padding;

template< typename PreReal , typename Real > ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< PreReal , 2 > > >	Geodesics< PreReal , Real >::parameterMetric;

template< typename PreReal , typename Real > SparseMatrix< Real , int >										Geodesics< PreReal , Real >::mass;
template< typename PreReal , typename Real > SparseMatrix< Real , int >										Geodesics< PreReal , Real >::stiffness;

template< typename PreReal , typename Real > SparseMatrix< Real , int >										Geodesics< PreReal , Real >::smoothImpulseMatrix;
template< typename PreReal , typename Real > SparseMatrix< Real , int >										Geodesics< PreReal , Real >::geodesicDistanceMatrix;

template< typename PreReal , typename Real > Real															Geodesics< PreReal , Real >::diffusionInterpolationWeight;
template< typename PreReal , typename Real > Real															Geodesics< PreReal , Real >::geodesicInterpolationWeight = 1e-6;

template< typename PreReal , typename Real > std::vector< TextureNodeInfo< PreReal > >						Geodesics< PreReal , Real >::textureNodes;
template< typename PreReal , typename Real > Image<int>														Geodesics< PreReal , Real >::nodeIndex;

template< typename PreReal , typename Real > int															Geodesics< PreReal , Real >::steps;
template< typename PreReal , typename Real > char															Geodesics< PreReal , Real >::stepsString[1024];

template< typename PreReal , typename Real > unsigned int													Geodesics< PreReal , Real >::levels;
template< typename PreReal , typename Real > HierarchicalSystem< PreReal , Real >							Geodesics< PreReal , Real >::hierarchy;

template< typename PreReal , typename Real > unsigned char *												Geodesics< PreReal , Real >::outputBuffer;
template< typename PreReal , typename Real > std::vector< MultigridLevelIndices< Real > >					Geodesics< PreReal , Real >::multigridIndices;

//Impulse Smoothing
template< typename PreReal , typename Real > std::vector< SystemCoefficients< Real > >						Geodesics< PreReal , Real >::multigridSmoothImpulseCoefficients;
template< typename PreReal , typename Real > std::vector<MultigridLevelVariables<Real>>						Geodesics< PreReal , Real >::multigridSmoothImpulseVariables;
template< typename PreReal , typename Real > VCycleSolvers< typename Geodesics< PreReal , Real >::DirectSolver > Geodesics< PreReal , Real >::smoothImpulseSolvers;

//Geodesic Distance
template< typename PreReal , typename Real > std::vector< SystemCoefficients< Real > >						Geodesics< PreReal , Real >::multigridGeodesicDistanceCoefficients;
template< typename PreReal , typename Real > std::vector<MultigridLevelVariables<Real>>						Geodesics< PreReal , Real >::multigridGeodesicDistanceVariables;
template< typename PreReal , typename Real > VCycleSolvers< typename Geodesics< PreReal , Real >::DirectSolver > Geodesics< PreReal , Real >::geodesicDistanceSolvers;

template< typename PreReal , typename Real > typename Geodesics< PreReal , Real >::DirectSolver				Geodesics< PreReal , Real >::fineSmoothImpulseSolver;
template< typename PreReal , typename Real > typename Geodesics< PreReal , Real >::DirectSolver				Geodesics< PreReal , Real >::fineGeodesicDistanceSolver;

//Samples
template< typename PreReal , typename Real > unsigned int													Geodesics< PreReal , Real >::impulseTexel = static_cast< unsigned int >(-1);
template< typename PreReal , typename Real > std::vector<Point3D< float > >									Geodesics< PreReal , Real >::textureNodePositions;

template< typename PreReal , typename Real > Real															Geodesics< PreReal , Real >::smoothImpulseRange;
template< typename PreReal , typename Real > Real															Geodesics< PreReal , Real >::geodesicDistanceRange;

template< typename PreReal , typename Real > GradientIntegrator< Real >										Geodesics< PreReal , Real >::gradientIntegrator;
template< typename PreReal , typename Real > typename GradientIntegrator< Real >::template Scratch< Real , Real > Geodesics< PreReal , Real >::gradientIntegratorScratch;

template< typename PreReal , typename Real > int															Geodesics< PreReal , Real >::updateCount = -1;


template< typename PreReal , typename Real >
void Geodesics< PreReal , Real >::ComputeExactSolution( void )
{
	//(1) Smoothing impulse	
	solve( fineSmoothImpulseSolver , multigridSmoothImpulseVariables[0].x , multigridSmoothImpulseVariables[0].rhs );

	//(1) Integrating vector field	
	std::vector< Real >& fineGeodesicDistanceRHS = multigridGeodesicDistanceVariables[0].rhs;

	auto VectorFunction = []( Point2D< Real > v , SquareMatrix< Real , 2 > tensor )
		{
			Point2D< Real > _v = tensor * v;
			Real len2 = Point2D< Real >::Dot( v , _v );
			if( len2>0 ) return -v / (Real)sqrt( len2 );
			else         return -v;
		};

	gradientIntegrator( multigridSmoothImpulseVariables[0].x , VectorFunction , gradientIntegratorScratch , fineGeodesicDistanceRHS );

	//(3) Update geodesic distance solution	
	solve( fineGeodesicDistanceSolver , multigridGeodesicDistanceVariables[0].x , fineGeodesicDistanceRHS );

	Real expectedMinDistance = multigridGeodesicDistanceVariables[0].x[impulseTexel];

	ThreadPool::ParallelFor
		(
			0 , multigridGeodesicDistanceVariables[0].x.size() ,
			[&]( unsigned int , size_t i ){ multigridGeodesicDistanceVariables[0].x[i] -= expectedMinDistance; }
		);
}

template< typename PreReal , typename Real >
void Geodesics< PreReal , Real >::UpdateOutputBuffer( const std::vector< Real > &solution )
{
	//Real expectedMaxDistance = solution[impulseTexel];

	Real attenuationRadius = (Real)0.2; //value within 0 and 0.5
	Real attenuationStart  = (Real)0.8 - attenuationRadius;
	ThreadPool::ParallelFor
		(
			0 , textureNodes.size() ,
			[&]( unsigned int , size_t i )
			{
				int ci = textureNodes[i].ci;
				int cj = textureNodes[i].cj;
				int offset = 3 * (textureWidth*cj + ci);
				Real value = solution[i];
				Point3D< float > red( 255.f , 64.f , 64.f ) , blue( 64.f , 64.f , 255.f ) , green( 64.f , 255.f , 64.f ); ;
				Point3D< float > color = value < 1.f ? red * (1.f - value) + blue * value : blue * (2.f - value) + green * ( value-1.f );
				float scaledDistance = value * 60.f;
				int closestInteger = floor(scaledDistance);

				float residual = scaledDistance - float(closestInteger);
				if( closestInteger%2!=0 ) residual = 1.f - residual;
				float attenuationWeight = ( residual-attenuationStart ) / ( 2.f * attenuationRadius );
				attenuationWeight = std::min< float >( std::max< float >( 0 , attenuationWeight ), 1.f );
				attenuationWeight = attenuationWeight*attenuationWeight*(3.f - 2.f * attenuationWeight );
				color = color * ( 1.f-attenuationWeight );

				outputBuffer[offset + 0] = (unsigned char)(color[0]);
				outputBuffer[offset + 1] = (unsigned char)(color[1]);
				outputBuffer[offset + 2] = (unsigned char)(color[2]);
			}
		);

	glBindTexture(GL_TEXTURE_2D, visualization.textureBuffer);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, textureWidth, textureHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)&outputBuffer[0]);
	glBindTexture(GL_TEXTURE_2D, 0);

	glutPostRedisplay();
}

template< typename PreReal , typename Real >
void Geodesics< PreReal , Real >::Idle( void )
{
	visualization.Idle();
	int selectedTexel = -1;
	if( mouseSelectionActive )
	{
		Point3D< float > selectedPoint;
		if( visualization.showMesh )
		{
			visualization.select(mouseX, mouseY, selectedPoint);
			float minDistance = FLT_MAX;
			for( int i=0 ; i<textureNodePositions.size() ; i++ )
			{
				float squaredDistance = Point3D< float >::SquareNorm( textureNodePositions[i]-selectedPoint );
				if( squaredDistance<minDistance ) minDistance = squaredDistance , selectedTexel = i;
			}
		}
		else
		{
			Point2D< float > ip = visualization.selectImagePos(mouseX, mouseY);
			//printf("Texture Coord %f %f \n", ip[0], ip[1]);
			int i = floor(ip[0] * float(nodeIndex.res(0)) - 0.5f);
			int j = floor((1.0 - ip[1]) * float(nodeIndex.res(1)) - 0.5f);
			//printf("Image pos %d %d \n", i, j);
			if( i>=0 && i<(int)nodeIndex.res(0) && j>=0 && j<(int)nodeIndex.res(1) ) selectedTexel = nodeIndex(i, j);
		}

		if( impulseTexel!=selectedTexel && selectedTexel!=-1 )
		{
			impulseTexel = selectedTexel;
			steps = 0;
			memset( &multigridSmoothImpulseVariables[0].rhs[0] , 0 , multigridSmoothImpulseVariables[0].rhs.size() * sizeof(Real) );
			multigridSmoothImpulseVariables[0].rhs[impulseTexel] = 1.0;
			memset( &multigridSmoothImpulseVariables[0].x[0] , 0 , multigridSmoothImpulseVariables[0].x.size() * sizeof(Real) );
			memset( &multigridGeodesicDistanceVariables[0].x[0] , 0 , multigridGeodesicDistanceVariables[0].x.size() * sizeof(Real) );
		}
	}

	if( impulseTexel!=-1 )
	{
		if( updateCount && !visualization.promptCallBack )
		{
			UpdateSolution();
			if( updateCount>0 ) updateCount--;
			steps++;
			sprintf( stepsString , "Steps: %d" , steps );
			UpdateOutputBuffer( multigridGeodesicDistanceVariables[0].x );
		}
	}
	if( !visualization.promptCallBack ) glutPostRedisplay();	
}

template< typename PreReal , typename Real>
void Geodesics< PreReal , Real >::MouseFunc( int button , int state , int x , int y )
{
	visualization.newX = x; visualization.newY = y;
	visualization.rotating = visualization.scaling = visualization.panning = false;

	if (state == GLUT_DOWN && glutGetModifiers() & GLUT_ACTIVE_SHIFT) {

		int selectedTexel = -1;
		if (visualization.showMesh)
		{
			Point3D< float > selectedPoint;
			if (visualization.select(x, y, selectedPoint))
			{
				mouseSelectionActive = true;
				mouseX = x;
				mouseY = y;
				float minDistance = FLT_MAX;
				for( int i=0 ; i<textureNodePositions.size() ; i++ )
				{
					float squaredDistance = Point3D< float >::SquareNorm(textureNodePositions[i] - selectedPoint);
					if( squaredDistance<minDistance ) minDistance = squaredDistance , selectedTexel = i;
				}
			}
		}
		else {
			Point2D<float> ip = visualization.selectImagePos(x, y);
			int i = floor(ip[0] * float(nodeIndex.res(0)) - 0.5f);
			int j = floor((1.0 - ip[1])*float(nodeIndex.res(1)) - 0.5f);
			if( i>=0 && i<(int)nodeIndex.res(0) && j >= 0 && j<(int)nodeIndex.res(1) )
			{
				mouseSelectionActive = true;
				selectedTexel = nodeIndex(i, j);
			}
		}
		if( selectedTexel!=-1 && selectedTexel!=impulseTexel )
		{
			impulseTexel = selectedTexel;
			steps = 0;
			memset(&multigridSmoothImpulseVariables[0].rhs[0], 0, multigridSmoothImpulseVariables[0].rhs.size() * sizeof(Real));
			multigridSmoothImpulseVariables[0].rhs[impulseTexel] = 1.0;
			memset(&multigridSmoothImpulseVariables[0].x[0], 0, multigridSmoothImpulseVariables[0].x.size() * sizeof(Real));
			memset(&multigridGeodesicDistanceVariables[0].x[0], 0, multigridGeodesicDistanceVariables[0].x.size() * sizeof(Real));
		}
		if( impulseTexel!=-1 )
		{
			if( UseDirectSolver.set )
			{
				ComputeExactSolution();
				UpdateOutputBuffer( multigridGeodesicDistanceVariables[0].x );
			}
			glutPostRedisplay();
		}
	}
	else
	{
		mouseSelectionActive = false;
		if( ( button==GLUT_LEFT_BUTTON || button==GLUT_RIGHT_BUTTON ) && glutGetModifiers() & GLUT_ACTIVE_CTRL ) visualization.panning = true;
		else if( button==GLUT_LEFT_BUTTON  ) visualization.rotating = true;
		else if( button==GLUT_RIGHT_BUTTON ) visualization.scaling  = true;
	}
}

template< typename PreReal , typename Real >
void Geodesics< PreReal , Real >::MotionFunc( int x , int y )
{
	if (mouseSelectionActive) {
		mouseX = x;
		mouseY = y;
		glutPostRedisplay();
	}
	else
	{
		if( !visualization.showMesh )
		{
			visualization.oldX = visualization.newX, visualization.oldY = visualization.newY, visualization.newX = x, visualization.newY = y;

			if (visualization.panning) visualization.xForm.offset[0] -= (visualization.newX - visualization.oldX) / visualization.imageToScreenScale(), visualization.xForm.offset[1] += (visualization.newY - visualization.oldY) / visualization.imageToScreenScale();
			else
			{
				float dz = (float)pow( 1.1 , (float)( visualization.newY-visualization.oldY ) / 8.f );
				visualization.xForm.zoom *= dz;
			}

		}
		else
		{
			visualization.oldX = visualization.newX, visualization.oldY = visualization.newY, visualization.newX = x, visualization.newY = y;
			int screenSize = std::min< int >(visualization.screenWidth, visualization.screenHeight);
			float rel_x = (visualization.newX - visualization.oldX) / (float)screenSize * 2;
			float rel_y = (visualization.newY - visualization.oldY) / (float)screenSize * 2;

			float pRight = rel_x * visualization.zoom, pUp = -rel_y * visualization.zoom;
			float pForward = rel_y * visualization.zoom;
			float rRight = -rel_y, rUp = -rel_x;

			if     ( visualization.rotating ) visualization.camera.rotateUp( -rUp ) , visualization.camera.rotateRight( -rRight );
			else if( visualization.scaling  ) visualization.camera.translate( visualization.camera.forward*pForward );
			else if( visualization.panning  ) visualization.camera.translate( -( visualization.camera.right*pRight + visualization.camera.up*pUp ) );
		}
		glutPostRedisplay();
	}
}

template< typename PreReal , typename Real>
void Geodesics< PreReal , Real >::ExportTextureCallBack( Visualization * /*v*/ , const char* prompt )
{

	Image< Point3D< float > > outputImage;
	outputImage.resize( textureWidth , textureHeight );
	for( int i=0 ; i<outputImage.size() ; i++ ) outputImage[i] = Point3D< float >( outputBuffer[3*i] , outputBuffer[3*i+1] , outputBuffer[3*i+2]) / float(255.0);
	padding.unpad( outputImage );
	WriteImage< 8 >( outputImage , prompt );
}

template< typename PreReal , typename Real >
void Geodesics< PreReal , Real >::UpdateSolution( void )
{
	// (1) Update smoothed input solution
	VCycle( multigridSmoothImpulseVariables , multigridSmoothImpulseCoefficients , multigridIndices , smoothImpulseSolvers , false , false );

	// (2) Integrate normalized vector field
	auto VectorFunction = []( Point2D< Real > v , SquareMatrix< Real , 2 > tensor )
		{
			Point2D< Real > _v = tensor * v;
			Real len2 = Point2D<Real>::Dot( v , _v );
			if( len2>0 ) return -v / (Real)sqrt( len2 );
			else         return -v;
		};

	gradientIntegrator( multigridSmoothImpulseVariables[0].x , VectorFunction , gradientIntegratorScratch , multigridGeodesicDistanceVariables[0].rhs );

	// (3) Update geodesic distance solution	
	VCycle( multigridGeodesicDistanceVariables , multigridGeodesicDistanceCoefficients , multigridIndices, geodesicDistanceSolvers , false , false );

	Real expectedMinDistance = multigridGeodesicDistanceVariables[0].x[impulseTexel];

	ThreadPool::ParallelFor
	(
		0 , multigridGeodesicDistanceVariables[0].x.size() ,
		[&]( unsigned int , size_t i ){  multigridGeodesicDistanceVariables[0].x[i] -= expectedMinDistance; }
	);
}

template< typename PreReal , typename Real >
void Geodesics< PreReal , Real >::InitializeSystem( int width , int height )
{
	Miscellany::PerformanceMeter pMeter( '.' );
	ExplicitIndexVector< ChartIndex , AtlasChart< PreReal > > atlasCharts;
	MultigridBlockInfo multigridBlockInfo( MultigridBlockWidth.value , MultigridBlockHeight.value , MultigridPaddedWidth.value , MultigridPaddedHeight.value );
	InitializeHierarchy( mesh , width , height , levels , textureNodes , hierarchy , atlasCharts , multigridBlockInfo );
	if( Verbose.set ) std::cout << pMeter( "Hierarchy" ) << std::endl;

	//Initialize node index
	nodeIndex.resize( width , height );
	for( int i=0 ; i<nodeIndex.size() ; i++ ) nodeIndex[i] = -1;
	for( int i=0 ; i<textureNodes.size() ; i++ ) nodeIndex( textureNodes[i].ci , textureNodes[i].cj ) = i;

	BoundaryProlongationData< Real > boundaryProlongation;
	InitializeBoundaryProlongationData( hierarchy.gridAtlases[0] , boundaryProlongation );

	MassAndStiffnessOperators< Real > massAndStiffnessOperators;

	InitializeMetric( mesh , EMBEDDING_METRIC , atlasCharts , parameterMetric );

	pMeter.reset();
	OperatorInitializer::Initialize( MatrixQuadrature.value , massAndStiffnessOperators , hierarchy.gridAtlases[0] , parameterMetric , atlasCharts , gradientIntegrator , VectorFieldQuadrature.value , !PreciseIntegration.set );
	gradientIntegratorScratch = gradientIntegrator.template getScratch< Real , Real >();
	if( Verbose.set ) std::cout << pMeter( "System" ) << std::endl;

	if( UseDirectSolver.set )
	{
		FullMatrixConstruction( hierarchy.gridAtlases[0] , massAndStiffnessOperators.massCoefficients , mass );
		FullMatrixConstruction( hierarchy.gridAtlases[0] , massAndStiffnessOperators.stiffnessCoefficients , stiffness );
		smoothImpulseMatrix = mass * diffusionInterpolationWeight + stiffness;
		geodesicDistanceMatrix = mass * geodesicInterpolationWeight + stiffness;
		if( Verbose.set ) std::cout << pMeter( "Assembled matrice" ) << std::endl;
	}

//////////////////////////////////// Initialize multigrid indices

	multigridIndices.resize(levels);
	for( unsigned int i=0 ; i<levels ; i++ )
	{
		const GridAtlas< PreReal , Real > &gridAtlas = hierarchy.gridAtlases[i];
		multigridIndices[i].threadTasks = gridAtlas.threadTasks;
		multigridIndices[i].boundaryToCombined = gridAtlas.indexConverter.boundaryToCombined();
		multigridIndices[i].segmentedLines = gridAtlas.segmentedLines;
		multigridIndices[i].rasterLines = gridAtlas.rasterLines;
		multigridIndices[i].restrictionLines = gridAtlas.restrictionLines;
		multigridIndices[i].prolongationLines = gridAtlas.prolongationLines;
		if( i<levels-1 ) multigridIndices[i].boundaryRestriction = hierarchy.boundaryRestriction[i];
	}
	if( Verbose.set ) std::cout << pMeter( "MG indices" ) << std::endl;

//////////////////////////////////// Initialize multigrid coefficients

	UpdateLinearSystem( diffusionInterpolationWeight , (Real)1. , hierarchy , multigridSmoothImpulseCoefficients , massAndStiffnessOperators , smoothImpulseSolvers , fineSmoothImpulseSolver , smoothImpulseMatrix , DetailVerbose.set , true , UseDirectSolver.set );
	UpdateLinearSystem(  geodesicInterpolationWeight , (Real)1. , hierarchy , multigridGeodesicDistanceCoefficients , massAndStiffnessOperators , geodesicDistanceSolvers , fineGeodesicDistanceSolver , geodesicDistanceMatrix , DetailVerbose.set , true , UseDirectSolver.set );
	if( Verbose.set ) std::cout << pMeter( "MG coefficients" ) << std::endl;

//////////////////////////////////// Initialize multigrid variables

	multigridSmoothImpulseVariables.resize(levels);
	for( unsigned int i=0 ; i<levels ; i++ )
	{
		MultigridLevelVariables<Real> & variables = multigridSmoothImpulseVariables[i];
		variables.x.resize( static_cast< unsigned int >( hierarchy.gridAtlases[i].endCombinedTexelIndex ) );
		variables.rhs.resize( static_cast< unsigned int >( hierarchy.gridAtlases[i].endCombinedTexelIndex ) );
		variables.residual.resize( static_cast< unsigned int >( hierarchy.gridAtlases[i].endCombinedTexelIndex ) );
		variables.boundary_rhs.resize( hierarchy.gridAtlases[i].indexConverter.numBoundary() );
		variables.boundary_value.resize( hierarchy.gridAtlases[i].indexConverter.numBoundary() );
		variables.variable_boundary_value.resize( hierarchy.gridAtlases[i].indexConverter.numBoundary() );
	}

	multigridGeodesicDistanceVariables.resize( levels );
	for( unsigned int i=0 ; i<levels ; i++ )
	{
		MultigridLevelVariables<Real> & variables = multigridGeodesicDistanceVariables[i];
		variables.x.resize( static_cast< unsigned int >( hierarchy.gridAtlases[i].endCombinedTexelIndex ) );
		variables.rhs.resize( static_cast< unsigned int >( hierarchy.gridAtlases[i].endCombinedTexelIndex ) );
		variables.residual.resize( static_cast< unsigned int >( hierarchy.gridAtlases[i].endCombinedTexelIndex ) );
		variables.boundary_rhs.resize( hierarchy.gridAtlases[i].indexConverter.numBoundary() );
		variables.boundary_value.resize( hierarchy.gridAtlases[i].indexConverter.numBoundary() );
		variables.variable_boundary_value.resize( hierarchy.gridAtlases[i].indexConverter.numBoundary() );
	}
	if( Verbose.set ) std::cout << pMeter( "MG variables" ) << std::endl;
}

 template< typename PreReal , typename Real>
void Geodesics< PreReal , Real >::InitializeVisualization( int width , int height )
{
	outputBuffer = new unsigned char[ height*width* 3];
	memset( outputBuffer , 204 , height * width * 3 * sizeof(unsigned char) );

	unsigned int tCount = (unsigned int)mesh.surface.triangles.size();

	visualization.triangles.resize( tCount );
	visualization.vertices.resize( 3*tCount );
	visualization.textureCoordinates.resize( 3*tCount );
	visualization.normals.resize( 3*tCount );
	visualization.colors.resize( 3*tCount , Point3D< float >( 0.75f , 0.75f , 0.75f ) );

	for( unsigned int t=0 , idx=0 ; t<tCount ; t++ )
	{
		Point3D< float > n = mesh.surfaceTriangle( t ).normal();
		n /= Point3D< float >::Length( n );

		for( int k=0 ; k<3 ; k++ , idx++ )
		{
			visualization.triangles[t][k] = idx;
			visualization.vertices[idx] = mesh.surface.vertices[ mesh.surface.triangles[t][k] ];
			visualization.normals[idx] = n;
			visualization.textureCoordinates[idx] = mesh.texture.vertices[ mesh.texture.triangles[t][k] ];
		}
	}
	
	std::vector< unsigned int > boundaryHalfEdges = mesh.texture.boundaryHalfEdges();

	for( unsigned int e=0 ; e<boundaryHalfEdges.size() ; e++ )
	{
		SimplexIndex< 1 > eIndex = mesh.surface.edgeIndex( boundaryHalfEdges[e] );
		for( unsigned int i=0 ; i<2 ; i++ ) visualization.chartBoundaryVertices.push_back( Point3D< float >( mesh.surface.vertices[ eIndex[i] ] ) );
	}

	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization,  's' , "export texture" , "Output Texture" , ExportTextureCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , ' ' , "toggle update" , ToggleUpdateCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , '+' , "increment update" , IncrementUpdateCallBack ) );
	visualization.info.push_back( stepsString );

	visualization.UpdateVertexBuffer();
	visualization.UpdateFaceBuffer();

	visualization.textureImage.resize( textureWidth , textureHeight );
	for( int i=0 ; i<textureWidth*textureHeight ; i++ ) visualization.textureImage[i] = Point3D< float >( 0.5f ,  0.5f ,  0.5f );
	for( int i=0 ; i<textureNodes.size() ; i++ )
	{
		int ci = textureNodes[i].ci;
		int cj = textureNodes[i].cj;
		visualization.textureImage(ci,cj) = Point3D< float >( 0.8f , 0.8f , 0.8f );
	}
	visualization.UpdateTextureBuffer();
}

template< typename PreReal , typename Real >
void Geodesics< PreReal , Real >::ToggleUpdateCallBack( Visualization * /*v*/ , const char * /*prompt*/ )
{
	if( updateCount ) updateCount = 0;
	else              updateCount = -1;
}

template< typename PreReal , typename Real >
void Geodesics< PreReal , Real >::IncrementUpdateCallBack( Visualization * /*v*/ , const char * /*prompt*/ )
{
	if( updateCount<0 ) updateCount = 1;
	else updateCount++;
}

template< typename PreReal , typename Real >
void Geodesics< PreReal , Real >::Init( void )
{
	sprintf( stepsString , "Steps: 0" );
//	levels = std::max< unsigned int >( Levels.value , 1 );
	levels = Levels.value;
	diffusionInterpolationWeight = DiffusionInterpolationWeight.value;
	textureWidth = Width.value;
	textureHeight = Height.value;

	mesh.read( Input.value , DetailVerbose.set , CollapseEpsilon.value );

	if( RandomJitter.set )
	{
		if( RandomJitter.value ) srand( RandomJitter.value );
		else                     srand( time(NULL) );
		PreReal jitterScale = (PreReal)1e-3 / std::max< int >( textureWidth , textureHeight );
		for( int i=0 ; i<mesh.texture.vertices.size() ; i++ ) mesh.texture.vertices[i] += Point2D< PreReal >( (PreReal)1. - Random< PreReal >()*2 , (PreReal)1. - Random<PreReal>()*2 ) * jitterScale;
	}

	{
		padding = Padding::Init( textureWidth , textureHeight , mesh.texture.vertices , DetailVerbose.set );
		padding.pad( textureWidth , textureHeight , mesh.texture.vertices );
		textureWidth  += padding.width();
		textureHeight += padding.height();
	}

	//Define centroid and scale for visualization
	Point3D< PreReal > centroid = mesh.surface.centroid();
	PreReal radius = mesh.surface.boundingRadius( centroid );
	for( unsigned int i=0 ; i<mesh.surface.vertices.size() ; i++ ) mesh.surface.vertices[i] = ( mesh.surface.vertices[i]-centroid ) / radius;

	Miscellany::PerformanceMeter pMeter( '.' );
	InitializeSystem( textureWidth , textureHeight );

	if( Verbose.set )
	{
		std::cout << pMeter( "Initialized" ) << std::endl;
		std::cout << "Resolution: " << textureNodes.size() << " / " << textureWidth << " x " << textureHeight << std::endl;
	}

	//Assign position to exterior nodes using barycentric-exponential map
	{
		FEM::RiemannianMesh< PreReal , unsigned int > rMesh( GetPointer( mesh.surface.triangles ) , mesh.surface.triangles.size() );
		rMesh.setMetricFromEmbedding( GetPointer( mesh.surface.vertices ) );
		rMesh.makeUnitArea();
		Pointer( FEM::CoordinateXForm< PreReal > ) xForms = rMesh.getCoordinateXForms();

		for( unsigned int i=0 ; i<textureNodes.size() ; i++ ) if( textureNodes[i].tID!=AtlasMeshTriangleIndex(-1) && !textureNodes[i].isInterior )
		{
			FEM::HermiteSamplePoint< PreReal > _p;
			_p.tIdx = static_cast< unsigned int >( textureNodes[i].tID );
			_p.p = Point2D< PreReal >( (PreReal)1./3 , (PreReal)1./3 );
			_p.v = textureNodes[i].barycentricCoords - _p.p;

#ifdef SANITY_CHECK
			rMesh.exp( xForms , _p , 0 , false );
#else // !SANITY_CHECK
			rMesh.exp( xForms , _p );
#endif // SANITY_CHECK

			textureNodes[i].tID = AtlasMeshTriangleIndex( _p.tIdx );
			textureNodes[i].barycentricCoords = _p.p;
		}
	}

	textureNodePositions.resize( textureNodes.size() );
	for( int i=0 ; i<textureNodePositions.size() ; i++ ) textureNodePositions[i] = mesh.surface( textureNodes[i] );
}

template< typename PreReal , typename Real >
void _main( int argc, char *argv[] )
{
	Geodesics< PreReal , Real >::Init();

	glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
	Geodesics< PreReal , Real >::visualization.displayMode = DisplayMode.value;
	if     ( DisplayMode.value==ONE_REGION_DISPLAY ) Geodesics< PreReal , Real >::visualization.screenWidth =  800 , Geodesics< PreReal , Real >::visualization.screenHeight = 800;
	else if( DisplayMode.value==TWO_REGION_DISPLAY ) Geodesics< PreReal , Real >::visualization.screenWidth = 1440 , Geodesics< PreReal , Real >::visualization.screenHeight = 720;
	Geodesics< PreReal , Real >::visualization.useNearestSampling = Nearest.set;

	glutInitWindowSize( Geodesics< PreReal , Real >::visualization.screenWidth , Geodesics< PreReal , Real >::visualization.screenHeight );

	glutInit( &argc , argv );
	char windowName[1024];
	sprintf( windowName , "Goedsics" );
	glutCreateWindow( windowName );
	if( glewInit()!=GLEW_OK ) MK_THROW( "glewInit failed" );
	glutDisplayFunc ( Geodesics< PreReal , Real >::Display );
	glutReshapeFunc ( Geodesics< PreReal , Real >::Reshape );
	glutMouseFunc   ( Geodesics< PreReal , Real >::MouseFunc );
	glutMotionFunc  ( Geodesics< PreReal , Real >::MotionFunc );
	glutKeyboardFunc( Geodesics< PreReal , Real >::KeyboardFunc );
	if( !UseDirectSolver.set ) glutIdleFunc( Geodesics< PreReal , Real >::Idle );
	if( CameraConfig.set ) Geodesics< PreReal , Real >::visualization.ReadSceneConfigurationCallBack( &Geodesics< PreReal , Real >::visualization , CameraConfig.value.c_str() );
	Geodesics< PreReal , Real >::InitializeVisualization( Geodesics< PreReal , Real >::textureWidth , Geodesics< PreReal , Real >::textureHeight );
	glutMainLoop();
}

int main( int argc , char* argv[] )
{
	CmdLineParse( argc-1 , argv+1 , params );
	if( !Input.set ){ ShowUsage( argv[0] ) ; return EXIT_FAILURE; }
	if( Serial.set ) ThreadPool::ParallelizationType = ThreadPool::ParallelType::NONE;
	if( !NoHelp.set )
	{
		printf( "+---------------------------------------------+\n" );
		printf( "| Interface Controls:                         |\n" );
		printf( "|    [Left Mouse]:                 rotate     |\n" );
		printf( "|    [Right Mouse]:                zoom       |\n" );
		printf( "|    [Left/Right Mouse] + [CTRL]:  pan        |\n" );
		printf( "|    [Left/Right Mouse] + [SHIFT]: set source |\n" );
		printf( "+---------------------------------------------+\n" );
	}
	try
	{
		if( Double.set ) _main< double , double >( argc , argv );
		else             _main< double , float  >( argc , argv );
	}
	catch( Exception &e )
	{
		printf( "%s\n" , e.what() );
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
