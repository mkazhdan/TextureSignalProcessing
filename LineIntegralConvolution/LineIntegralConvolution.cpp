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
#include <Misha/FEM.h>
#include <Src/Hierarchy.h>
#include <Src/Basis.h>
#include <Src/ChartDecomposition.h>
#include <Src/HSV.h>
#include <Src/Solver.h>
#include <Src/MassAndStiffness.h>
#include <Src/Padding.h>
#include <Src/TexturedMeshVisualization.h>

using namespace MishaK;

CmdLineParameter< std::string > Input( "in" );
CmdLineParameter< std::string > Output( "out" );
CmdLineParameter< int   > OutputVCycles( "outVCycles" , 10 );
CmdLineReadable MinimalCurvature( "minimal" );
CmdLineReadable Double( "double" );
CmdLineParameter< std::string > InVectorField( "inVF" );
CmdLineParameter< std::string > OutVectorField( "outVF" );
CmdLineParameter< int   > Width( "width" , 2048 );
CmdLineParameter< int   > Height( "height" , 2048 );
CmdLineParameter< float > LICInterpolationWeight( "licInterpolation" , 1e4 );
CmdLineParameter< float > SharpeningInterpolationWeight( "sharpInterpolation" , 1e4 );
CmdLineParameter< float > SharpeningGradientModulation( "sharpModulation" , 100 );
CmdLineParameter< float > AnisotropyExponent( "aExp" , 0.f );
CmdLineParameter< int   > NormalSmoothingIterations( "nIters" , 2 );
CmdLineParameter< float > NormalSmoothingInterpolation( "nInterpolation" , 1e3f );
CmdLineParameter< int   > Levels( "levels" , 4 );
CmdLineParameter< int   > MatrixQuadrature( "mQuadrature" , 6 );


CmdLineParameter< std::string > CameraConfig("camera");
CmdLineParameter< int   > DisplayMode("display", TWO_REGION_DISPLAY);


CmdLineParameter< int   > MultigridBlockHeight("mBlockH", 16);
CmdLineParameter< int   > MultigridBlockWidth("mBlockW", 128);
CmdLineParameter< int   > MultigridPaddedHeight("mPadH", 0);
CmdLineParameter< int   > MultigridPaddedWidth("mPadW", 2);

CmdLineParameter< int   > RandomJitter( "jitter" , 0 );
CmdLineReadable Verbose("verbose");
CmdLineReadable DetailVerbose("detail");
CmdLineReadable UseDirectSolver("useDirectSolver");
CmdLineReadable IntrinsicVectorField( "intrinsicVF" );
CmdLineReadable Serial( "serial" );
CmdLineReadable NoHelp( "noHelp" );

CmdLineParameter< double > CollapseEpsilon( "collapse" , 0 );

CmdLineReadable* params[] =
{
	&Input , &Output , &MinimalCurvature , &InVectorField , &OutVectorField , &IntrinsicVectorField , &Width,&Height , &LICInterpolationWeight , &SharpeningInterpolationWeight , &SharpeningGradientModulation , &CameraConfig, &Levels,&UseDirectSolver,&Serial,&DisplayMode,&MultigridBlockHeight,&MultigridBlockWidth,&MultigridPaddedHeight,&MultigridPaddedWidth,&Verbose,
	&DetailVerbose , &RandomJitter ,
	&Double ,
	&MatrixQuadrature ,
	&OutputVCycles ,
	&NoHelp , &AnisotropyExponent ,
	&NormalSmoothingIterations , &NormalSmoothingInterpolation ,
	&CollapseEpsilon ,
	NULL
};

void ShowUsage(const char* ex)
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , Input.name.c_str() );
	printf( "\t[--%s <output texture>\n" , Output.name.c_str() );
	printf( "\t[--%s <output v-cycles>=%d]\n" , OutputVCycles.name.c_str() , OutputVCycles.value );
	printf( "\t[--%s <input vector field file>\n" , InVectorField.name.c_str() );
	printf( "\t[--%s <output vector field file>\n" , OutVectorField.name.c_str() );
	printf( "\t[--%s <LIC interpolation weight>=%f]\n" , LICInterpolationWeight.name.c_str() , LICInterpolationWeight.value );
	printf( "\t[--%s <sharpening interpolation weight>=%f]\n" , SharpeningInterpolationWeight.name.c_str() , SharpeningInterpolationWeight.value   );
	printf( "\t[--%s <sharpening gradient modulation>=%f]\n" , SharpeningGradientModulation.name.c_str() , SharpeningGradientModulation.value );
	printf( "\t[--%s <texture width>=%d]\n" , Width.name.c_str()  , Width.value  );
	printf( "\t[--%s <texture height>=%d]\n", Height.name.c_str() , Height.value );
	printf( "\t[--%s <system matrix quadrature points per triangle>=%d]\n" , MatrixQuadrature.name.c_str() , MatrixQuadrature.value );
	printf( "\t[--%s]\n" , IntrinsicVectorField.name.c_str() );
	printf( "\t[--%s]\n" , MinimalCurvature.name.c_str() );
	printf( "\t[--%s]\n" , UseDirectSolver.name.c_str() );
	printf( "\t[--%s <jittering seed>]\n" , RandomJitter.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );

	printf( "\t[--%s <camera configuration file>\n" , CameraConfig.name.c_str());
	printf( "\t[--%s <hierarchy levels>=%d]\n" , Levels.name.c_str() , Levels.value );
	printf( "\t[--%s]\n" , DetailVerbose.name.c_str() );
	printf( "\t[--%s <display mode>=%d]\n" , DisplayMode.name.c_str() , DisplayMode.value );
	printf( "\t\t%d] One Region \n", ONE_REGION_DISPLAY);
	printf( "\t\t%d] Two Region \n", TWO_REGION_DISPLAY);

	printf( "\t[--%s <multigrid block width>=%d]\n"   , MultigridBlockWidth.name.c_str()   , MultigridBlockWidth.value   );
	printf( "\t[--%s <multigrid block height>=%d]\n"  , MultigridBlockHeight.name.c_str()  , MultigridBlockHeight.value  );
	printf( "\t[--%s <multigrid padded width>=%d]\n"  , MultigridPaddedWidth.name.c_str()  , MultigridPaddedWidth.value  );
	printf( "\t[--%s <multigrid padded height>=%d]\n" , MultigridPaddedHeight.name.c_str() , MultigridPaddedHeight.value );
	printf( "\t[--%s <normal smoothing iterations>=%d]\n" , NormalSmoothingIterations.name.c_str() , NormalSmoothingIterations.value );
	printf( "\t[--%s <normal smoothing interpolation>=%f]\n" , NormalSmoothingInterpolation.name.c_str() , NormalSmoothingInterpolation.value );
	printf( "\t[--%s <anisotropy exponent>=%f]\n" , AnisotropyExponent.name.c_str() , AnisotropyExponent.value );
	printf( "\t[--%s <collapse epsilon>=%g]\n" , CollapseEpsilon.name.c_str() , CollapseEpsilon.value );
	printf( "\t[--%s]\n" , Serial.name.c_str() );
	printf( "\t[--%s]\n" , NoHelp.name.c_str() );
	printf( "\t[--%s]\n" , Double.name.c_str() );
}

template< typename PreReal , typename Real >
class LineConvolution
{
public:
	static Real sharpeningGradientModulation;
	static Real sharpeningInterpolationWeight;
	static Real licInterpolationWeight;

	static OrientedTexturedTriangleMesh< PreReal > mesh;
	static int textureWidth;
	static int textureHeight;
	static int levels;

	static int steps;
	static char stepsString[];

	static Padding padding;

	static HierarchicalSystem< PreReal , Real > hierarchy;
	static std::vector< BilinearElementIndex > bilinearElementIndices;

	static std::vector< TextureNodeInfo< PreReal > > textureNodes;
	static Image< int > nodeIndex;

	static SparseMatrix< Real , int > anisotropicMass;
	static SparseMatrix< Real , int > anisotropicStiffness;
	static SparseMatrix< Real , int > mass;
	static SparseMatrix< Real , int > stiffness;
	static SparseMatrix< Real , int > lineConvolutionMatrix;
	static SparseMatrix< Real , int > modulationMatrix;

	static int impulseTexel;

	static std::vector< AtlasChart< PreReal > > atlasCharts;
	static std::vector<std::vector< SquareMatrix< PreReal , 2 > > > parameterMetric;

	static Real lineConvolutionRange;
	static Real modulationRange;

	static std::vector< Point3D< Real > > randSignal;

	static std::vector< Point3D< Real > > mass_x0;
	static std::vector< Point3D< Real > > stiffness_x0;

	//Impulse Smoothing
	static std::vector< SystemCoefficients< Real > > multigridLineConvolutionCoefficients;
	static std::vector< MultigridLevelVariables< Point3D< Real > > > multigridLineConvolutionVariables;

	//Geodesic Distance
	static std::vector< SystemCoefficients< Real > > multigridModulationCoefficients;
	static std::vector< MultigridLevelVariables< Point3D< Real > > > multigridModulationVariables;

#if defined( USE_CHOLMOD )
	typedef CholmodCholeskySolver< Real , 3 > DirectSolver;
#elif defined( USE_EIGEN )
	typedef EigenCholeskySolver< Real , 3 > DirectSolver;
#elif defined( USE_EIGEN_PARDISO )
	typedef EigenPardisoSolver< Real , 3 > DirectSolver;
#else
#error "[ERROR] No solver defined!"
#endif

	static VCycleSolvers< DirectSolver > lineConvolutionSolvers;
	static VCycleSolvers< DirectSolver > modulationSolvers;

	static DirectSolver fineLineConvolutionSolver;
	static DirectSolver fineModulationSolver;

	static std::vector<MultigridLevelIndices<Real>> multigridIndices;

	static SparseMatrix< Real , int > coarseBoundaryFineBoundaryProlongation;
	static SparseMatrix< Real , int > fineBoundaryCoarseBoundaryRestriction;
	static std::vector< Real > coarseBoundaryValues;
	static std::vector< Real > coarseBoundaryRHS;
	static std::vector< Real > fineBoundaryValues;
	static std::vector< Real > fineBoundaryRHS;

	// Anisotropic Linear Operators
	static SystemCoefficients< Real > anisoMassCoefficients;
	static SystemCoefficients< Real > anisoStiffnessCoefficients;

	// Isotropic Linear Operators
	static SystemCoefficients< Real > massCoefficients;
	static SystemCoefficients< Real > stiffnessCoefficients;

	static unsigned char * outputBuffer;

	//Visulization
	static TexturedMeshVisualization visualization;

	static void SetOutputBuffer( const std::vector< Point3D< Real > > & solution );
	static void UpdateOutputBuffer( const std::vector< Point3D< Real > > & solution );

	static void SharpeningInterpolationWeightCallBack(Visualization* v, const char* prompt);
	static void SharpeningGradientModulationCallBack(Visualization* v, const char* prompt);
	static void LICInterpolationWeightCallBack(Visualization* v, const char* prompt);

	static int updateCount;

	static void ToggleUpdateCallBack(Visualization* v, const char* prompt);
	static void IncrementUpdateCallBack( Visualization* v , const char* prompt );
	static void ExportTextureCallBack(Visualization* v, const char* prompt);

	static void Init( void );
	static void InitializeVisualization( void );
	static void ComputeExactSolution( bool verbose= false );
	static void UpdateSolution( bool verbose=false , bool detailVerbose=false );
	static void InitializeSystem( const FEM::RiemannianMesh< PreReal , unsigned int >& rMesh , int width , int height );
	static void Reset( void );

	static void Display(void) { visualization.Display(); }
	static void MouseFunc(int button, int state, int x, int y);
	static void MotionFunc(int x, int y);
	static void Reshape(int w, int h) { visualization.Reshape(w, h); }
	static void KeyboardFunc(unsigned char key, int x, int y) { visualization.KeyboardFunc( key , x , y ); }
	static void Idle();
};

template< typename PreReal , typename Real > Real														LineConvolution< PreReal , Real >::sharpeningGradientModulation;
template< typename PreReal , typename Real > Real														LineConvolution< PreReal , Real >::sharpeningInterpolationWeight;
template< typename PreReal , typename Real > Real														LineConvolution< PreReal , Real >::licInterpolationWeight;

template< typename PreReal , typename Real > OrientedTexturedTriangleMesh< PreReal >					LineConvolution< PreReal , Real >::mesh;
template< typename PreReal , typename Real > int														LineConvolution< PreReal , Real >::textureWidth;
template< typename PreReal , typename Real > int														LineConvolution< PreReal , Real >::textureHeight;

template< typename PreReal , typename Real > TexturedMeshVisualization									LineConvolution< PreReal , Real >::visualization( true );

template< typename PreReal , typename Real > std::vector< AtlasChart< PreReal > >						LineConvolution< PreReal , Real >::atlasCharts;
template< typename PreReal , typename Real > std::vector< std::vector< SquareMatrix< PreReal , 2 > > >	LineConvolution< PreReal , Real >::parameterMetric;

template< typename PreReal , typename Real > Padding													LineConvolution< PreReal , Real >::padding;
template< typename PreReal , typename Real > SparseMatrix< Real , int >									LineConvolution< PreReal , Real >::anisotropicMass;
template< typename PreReal , typename Real > SparseMatrix< Real , int >									LineConvolution< PreReal , Real >::anisotropicStiffness;
template< typename PreReal , typename Real > SparseMatrix< Real , int >									LineConvolution< PreReal , Real >::lineConvolutionMatrix;
template< typename PreReal , typename Real > SparseMatrix< Real , int >									LineConvolution< PreReal , Real >::modulationMatrix;
template< typename PreReal , typename Real > SparseMatrix< Real , int >									LineConvolution< PreReal , Real >::mass;
template< typename PreReal , typename Real > SparseMatrix< Real , int >									LineConvolution< PreReal , Real >::stiffness;
template< typename PreReal , typename Real > std::vector< TextureNodeInfo< PreReal > >					LineConvolution< PreReal , Real >::textureNodes;
template< typename PreReal , typename Real > Image<int>													LineConvolution< PreReal , Real >::nodeIndex;
template< typename PreReal , typename Real > std::vector< BilinearElementIndex >						LineConvolution< PreReal , Real >::bilinearElementIndices;

template< typename PreReal , typename Real > int														LineConvolution< PreReal , Real >::steps;
template< typename PreReal , typename Real > char														LineConvolution< PreReal , Real >::stepsString[1024];
template< typename PreReal , typename Real > int														LineConvolution< PreReal , Real >::levels;
template< typename PreReal , typename Real > HierarchicalSystem< PreReal , Real >						LineConvolution< PreReal , Real >::hierarchy;

template< typename PreReal , typename Real > unsigned char *											LineConvolution< PreReal , Real >::outputBuffer;
template< typename PreReal , typename Real > std::vector< MultigridLevelIndices< Real > >				LineConvolution< PreReal , Real >::multigridIndices;

//Impulse Smoothing
template< typename PreReal , typename Real > std::vector< SystemCoefficients< Real > >					LineConvolution< PreReal , Real >::multigridLineConvolutionCoefficients;
template< typename PreReal , typename Real > std::vector< MultigridLevelVariables< Point3D< Real > > >	LineConvolution< PreReal , Real >::multigridLineConvolutionVariables;
template< typename PreReal , typename Real > VCycleSolvers< typename LineConvolution< PreReal , Real >::DirectSolver >		LineConvolution< PreReal , Real >::lineConvolutionSolvers;

//Geodesic Distance
template< typename PreReal , typename Real > std::vector< SystemCoefficients< Real > >					LineConvolution< PreReal , Real >::multigridModulationCoefficients;
template< typename PreReal , typename Real > std::vector< MultigridLevelVariables< Point3D< Real > > >	LineConvolution< PreReal , Real >::multigridModulationVariables;
template< typename PreReal , typename Real > VCycleSolvers< typename LineConvolution< PreReal , Real >::DirectSolver >		LineConvolution< PreReal , Real >::modulationSolvers;

template< typename PreReal , typename Real > typename LineConvolution< PreReal , Real >::DirectSolver	LineConvolution< PreReal , Real >::fineLineConvolutionSolver;
template< typename PreReal , typename Real > typename LineConvolution< PreReal , Real >::DirectSolver	LineConvolution< PreReal , Real >::fineModulationSolver;

template< typename PreReal , typename Real > std::vector< Point3D< Real > >								LineConvolution< PreReal , Real >::randSignal;
template< typename PreReal , typename Real > std::vector< Point3D< Real > >								LineConvolution< PreReal , Real >::mass_x0;
template< typename PreReal , typename Real > std::vector< Point3D< Real > >								LineConvolution< PreReal , Real >::stiffness_x0;

template< typename PreReal , typename Real > int														LineConvolution< PreReal , Real >::impulseTexel = -1;

template< typename PreReal , typename Real > Real														LineConvolution< PreReal , Real >::lineConvolutionRange;
template< typename PreReal , typename Real > Real														LineConvolution< PreReal , Real >::modulationRange;

template< typename PreReal , typename Real > SparseMatrix< Real , int >									LineConvolution< PreReal , Real >::coarseBoundaryFineBoundaryProlongation;
template< typename PreReal , typename Real > SparseMatrix< Real , int >									LineConvolution< PreReal , Real >::fineBoundaryCoarseBoundaryRestriction;

template< typename PreReal , typename Real > std::vector< Real >										LineConvolution< PreReal , Real >::coarseBoundaryValues;
template< typename PreReal , typename Real > std::vector< Real >										LineConvolution< PreReal , Real >::coarseBoundaryRHS;
template< typename PreReal , typename Real > std::vector< Real >										LineConvolution< PreReal , Real >::fineBoundaryValues;
template< typename PreReal , typename Real > std::vector< Real >										LineConvolution< PreReal , Real >::fineBoundaryRHS;

template< typename PreReal , typename Real > SystemCoefficients< Real >									LineConvolution< PreReal , Real >::anisoMassCoefficients;
template< typename PreReal , typename Real > SystemCoefficients< Real >									LineConvolution< PreReal , Real >::anisoStiffnessCoefficients;
template< typename PreReal , typename Real > SystemCoefficients< Real >									LineConvolution< PreReal , Real >::massCoefficients;
template< typename PreReal , typename Real > SystemCoefficients< Real >									LineConvolution< PreReal , Real >::stiffnessCoefficients;

template< typename PreReal , typename Real > int														LineConvolution< PreReal , Real >::updateCount = 0;

template< typename PreReal , typename Real >
void LineConvolution< PreReal , Real >::ComputeExactSolution( bool verbose )
{
	Miscellany::Timer timer;

	// (1) Line Convolution	
	// RHS = Mass * randSignal * licInterpolationWeight
	MultiplyBySystemMatrix_NoReciprocals( anisoMassCoefficients , hierarchy.gridAtlases[0].boundaryGlobalIndex , hierarchy.gridAtlases[0].rasterLines , randSignal , multigridLineConvolutionVariables[0].rhs );
	ThreadPool::ParallelFor( 0 , textureNodes.size() , [&]( unsigned int , size_t i ){ multigridLineConvolutionVariables[0].rhs[i] *= licInterpolationWeight; } );

	timer.reset();
	solve( fineLineConvolutionSolver , multigridLineConvolutionVariables[0].x , multigridLineConvolutionVariables[0].rhs );
	if( verbose ) printf( "Line convolution %.4f\n" , timer.elapsed() );

	//(2) Compute modulation RHS
	mass_x0.resize( textureNodes.size() );
	MultiplyBySystemMatrix_NoReciprocals( massCoefficients , hierarchy.gridAtlases[0].boundaryGlobalIndex , hierarchy.gridAtlases[0].rasterLines , multigridLineConvolutionVariables[0].x , mass_x0 );

	stiffness_x0.resize( textureNodes.size() );
	MultiplyBySystemMatrix_NoReciprocals( stiffnessCoefficients , hierarchy.gridAtlases[0].boundaryGlobalIndex , hierarchy.gridAtlases[0].rasterLines , multigridLineConvolutionVariables[0].x , stiffness_x0 );

	ThreadPool::ParallelFor( 0 , textureNodes.size() , [&]( unsigned int , size_t i ){ multigridModulationVariables[0].rhs[i] = mass_x0[i] * sharpeningInterpolationWeight + stiffness_x0[i] * sharpeningGradientModulation; } );

	//(3) Modulation
	if( verbose ) timer.reset();
	solve( fineModulationSolver , multigridModulationVariables[0].x , multigridModulationVariables[0].rhs );
	if( verbose ) printf( "Modulation %.4f \n" , timer.elapsed() );
}

template< typename PreReal , typename  Real >
void LineConvolution< PreReal , Real >::SetOutputBuffer( const std::vector< Point3D< Real > >& solution )
{
	ThreadPool::ParallelFor
		(
			0 , textureNodes.size() ,
			[&]( unsigned int , size_t i )
			{
				int ci = textureNodes[i].ci;
				int cj = textureNodes[i].cj;
				int offset = 3 * (textureWidth*cj + ci);
				outputBuffer[offset+0] = (unsigned char)( std::min< float >( std::max< float >( 0 , solution[i][0] ) , 1.f )*255.f );
				outputBuffer[offset+1] = (unsigned char)( std::min< float >( std::max< float >( 0 , solution[i][1] ) , 1.f )*255.f );
				outputBuffer[offset+2] = (unsigned char)( std::min< float >( std::max< float >( 0 , solution[i][2] ) , 1.f )*255.f );
			}
		);
}

template< typename PreReal , typename  Real >
void LineConvolution< PreReal , Real >::UpdateOutputBuffer( const std::vector< Point3D< Real > >& solution )
{
	SetOutputBuffer( solution );

	glBindTexture( GL_TEXTURE_2D , visualization.textureBuffer );
	glTexImage2D( GL_TEXTURE_2D , 0 , GL_RGBA , textureWidth , textureHeight , 0 , GL_RGB , GL_UNSIGNED_BYTE , (GLvoid*)&outputBuffer[0] );
	glBindTexture( GL_TEXTURE_2D , 0 );
	glutPostRedisplay();
}

template< typename PreReal , typename  Real >
void LineConvolution< PreReal , Real >::Idle( void )
{
	if( updateCount && !visualization.promptCallBack )
	{
		UpdateSolution();
		if( updateCount>0 ) updateCount--;
		steps++;
		sprintf( stepsString , "Steps: %d" , steps );
	}
	UpdateOutputBuffer( multigridModulationVariables[0].x );
}

template< typename PreReal , typename  Real >
void LineConvolution< PreReal , Real >::MouseFunc( int button , int /*state*/ , int x , int y )
{

	visualization.newX = x; visualization.newY = y;
	visualization.rotating = visualization.scaling = visualization.panning = false;

	if( ( button==GLUT_LEFT_BUTTON || button==GLUT_RIGHT_BUTTON ) && glutGetModifiers() & GLUT_ACTIVE_CTRL) visualization.panning = true;
	else if( button==GLUT_LEFT_BUTTON  ) visualization.rotating = true;
	else if( button==GLUT_RIGHT_BUTTON ) visualization.scaling  = true;
}
template< typename PreReal , typename Real >
void LineConvolution< PreReal , Real >::MotionFunc(int x, int y) {

	if( !visualization.showMesh )
	{
		visualization.oldX = visualization.newX, visualization.oldY = visualization.newY, visualization.newX = x, visualization.newY = y;
		if( visualization.panning ) visualization.xForm.offset[0] -= ( visualization.newX-visualization.oldX ) / visualization.imageToScreenScale() , visualization.xForm.offset[1] += ( visualization.newY-visualization.oldY ) / visualization.imageToScreenScale();
		else
		{
			float dz = (float)pow( 1.1 , (float)( visualization.newY-visualization.oldY ) / 8.f );
			visualization.xForm.zoom *= dz;
		}

	}
	else
	{
		visualization.oldX = visualization.newX , visualization.oldY = visualization.newY , visualization.newX = x , visualization.newY = y;
		int screenSize = std::min< int >( visualization.screenWidth , visualization.screenHeight );
		float rel_x = (float)( visualization.newX - visualization.oldX ) / screenSize * 2;
		float rel_y = (float)( visualization.newY - visualization.oldY ) / screenSize * 2;

		float pRight = rel_x * visualization.zoom, pUp = -rel_y * visualization.zoom;
		float pForward = rel_y * visualization.zoom;
		float rRight = -rel_y, rUp = -rel_x;

		if     ( visualization.rotating ) visualization.camera.rotateUp( -rUp ) , visualization.camera.rotateRight( -rRight );
		else if( visualization.scaling  ) visualization.camera.translate( visualization.camera.forward*pForward);
		else if( visualization.panning  ) visualization.camera.translate( -( visualization.camera.right*pRight + visualization.camera.up*pUp ) );
	}
	glutPostRedisplay();
}
template< typename PreReal , typename Real >
void LineConvolution< PreReal , Real >::SharpeningInterpolationWeightCallBack( Visualization * /*v*/ , const char* prompt )
{
	for( int i=0 ; i<multigridLineConvolutionVariables[0].x.size() ; i++) multigridLineConvolutionVariables[0].x[i] *= 0;
	for( int i=0 ; i<multigridModulationVariables[0].x.size() ; i++) multigridModulationVariables[0].x[i] *= 0;

	sharpeningInterpolationWeight = atof(prompt);

	if( UseDirectSolver.set ) modulationMatrix = mass * sharpeningInterpolationWeight + stiffness;

	UpdateLinearSystem( sharpeningInterpolationWeight , (Real)1. , hierarchy , multigridModulationCoefficients , massCoefficients , stiffnessCoefficients , modulationSolvers , fineModulationSolver , modulationMatrix , DetailVerbose.set , false , UseDirectSolver.set );
	Reset();
	if( UseDirectSolver.set ) UpdateOutputBuffer( multigridModulationVariables[0].x );
}
template< typename PreReal , typename Real >
void LineConvolution< PreReal , Real >::LICInterpolationWeightCallBack( Visualization * /*v*/ , const char* prompt )
{
	for( int i=0 ; i<multigridLineConvolutionVariables[0].x.size() ; i++ ) multigridLineConvolutionVariables[0].x[i] *= 0;
	for( int i=0 ; i<multigridModulationVariables[0].x.size() ; i++ ) multigridModulationVariables[0].x[i] *= 0;

	licInterpolationWeight = atof(prompt);

	if( UseDirectSolver.set ) lineConvolutionMatrix = anisotropicMass * licInterpolationWeight + anisotropicStiffness;

	UpdateLinearSystem( licInterpolationWeight , (Real)1. , hierarchy , multigridLineConvolutionCoefficients , anisoMassCoefficients , anisoStiffnessCoefficients , lineConvolutionSolvers , fineLineConvolutionSolver , lineConvolutionMatrix , DetailVerbose.set , false , UseDirectSolver.set );

	Reset();
	if( UseDirectSolver.set ) UpdateOutputBuffer( multigridModulationVariables[0].x );
}
template< typename PreReal , typename Real >
void LineConvolution< PreReal , Real >::SharpeningGradientModulationCallBack( Visualization * /*v*/ , const char* prompt )
{
	for( int i=0 ; i<multigridLineConvolutionVariables[0].x.size() ; i++ ) multigridLineConvolutionVariables[0].x[i] *= 0;
	for( int i=0 ; i<multigridModulationVariables[0].x.size() ; i++ ) multigridModulationVariables[0].x[i] *= 0;

	sharpeningGradientModulation = atof(prompt);
	Reset();
	if( UseDirectSolver.set ) UpdateOutputBuffer( multigridModulationVariables[0].x );
}

template< typename PreReal , typename Real >
void LineConvolution< PreReal , Real >::ToggleUpdateCallBack( Visualization * /*v*/ , const char * /*prompt*/ )
{
	if( updateCount ) updateCount =  0;
	else              updateCount = -1;
}

template< typename PreReal , typename Real >
void LineConvolution< PreReal , Real >::IncrementUpdateCallBack( Visualization * /*v*/ , const char * /*prompt*/ )
{
	if( updateCount<0 ) updateCount = 1;
	else updateCount++;
}

template< typename PreReal , typename Real >
void LineConvolution< PreReal , Real >::ExportTextureCallBack( Visualization * /*v*/ , const char *prompt )
{
	Image< Point3D< float > > outputImage;
	outputImage.resize( textureWidth , textureHeight );
	for( int i=0 ; i<outputImage.size() ; i++ ) outputImage[i] = Point3D< float >( outputBuffer[3*i] , outputBuffer[3*i+1] , outputBuffer[3*i+2] ) / 255.f;
	padding.unpad( outputImage );
	WriteImage< 8 >( outputImage , prompt );
}

template< typename PreReal , typename Real >
void LineConvolution< PreReal , Real >::UpdateSolution( bool verbose , bool detailVerbose )
{
	Miscellany::Timer timer;

	// (1) Update smoothed input solution
	MultiplyBySystemMatrix_NoReciprocals( anisoMassCoefficients , hierarchy.gridAtlases[0].boundaryGlobalIndex , hierarchy.gridAtlases[0].rasterLines , randSignal , multigridLineConvolutionVariables[0].rhs );
	ThreadPool::ParallelFor( 0 , textureNodes.size() , [&]( unsigned int , size_t i ){ multigridLineConvolutionVariables[0].rhs[i] *= licInterpolationWeight; } );
	if( verbose ) timer.reset();
	VCycle( multigridLineConvolutionVariables , multigridLineConvolutionCoefficients , multigridIndices , lineConvolutionSolvers , detailVerbose , detailVerbose );
	if( verbose ) printf( "Smoothing impulse %.4f\n" , timer.elapsed() );

	// (2) Compute modulation RHS
	mass_x0.resize(textureNodes.size());
	MultiplyBySystemMatrix_NoReciprocals( massCoefficients , hierarchy.gridAtlases[0].boundaryGlobalIndex , hierarchy.gridAtlases[0].rasterLines , multigridLineConvolutionVariables[0].x , mass_x0 );

	stiffness_x0.resize(textureNodes.size());
	MultiplyBySystemMatrix_NoReciprocals( stiffnessCoefficients , hierarchy.gridAtlases[0].boundaryGlobalIndex , hierarchy.gridAtlases[0].rasterLines , multigridLineConvolutionVariables[0].x , stiffness_x0 );

	ThreadPool::ParallelFor( 0 , textureNodes.size() , [&]( unsigned int , size_t i ){ multigridModulationVariables[0].rhs[i] = mass_x0[i] * sharpeningInterpolationWeight + stiffness_x0[i] * sharpeningGradientModulation; } );

	// (3) Update geodesic distance solution	
	if( verbose ) timer.reset();
	VCycle( multigridModulationVariables , multigridModulationCoefficients , multigridIndices , modulationSolvers , detailVerbose , detailVerbose );
	if( verbose ) printf( "Solving geodesic distance %.4f\n" , timer.elapsed() );
}

template< typename PreReal , typename Real >
void LineConvolution< PreReal , Real >::InitializeSystem( const FEM::RiemannianMesh< PreReal , unsigned int >& rMesh , int width , int height )
{
	Miscellany::Timer timer;
	MultigridBlockInfo multigridBlockInfo(MultigridBlockWidth.value, MultigridBlockHeight.value, MultigridPaddedWidth.value, MultigridPaddedHeight.value, 0);
	InitializeHierarchy( mesh , width , height , levels , textureNodes , bilinearElementIndices , hierarchy , atlasCharts , multigridBlockInfo , true , DetailVerbose.set );
	if( Verbose.set ) printf( "\tInitialized hierarchy: %.2f(s)\n" , timer.elapsed() );

	//Initialize node index
	nodeIndex.resize(width, height);
	for( int i=0 ; i<nodeIndex.size() ; i++ ) nodeIndex[i] = -1;
	for( int i=0 ; i<textureNodes.size() ; i++ )
	{
		if( textureNodes[i].ci<0 || textureNodes[i].ci>textureWidth-1 || textureNodes[i].cj<0 || textureNodes[i].cj>textureHeight-1 )
			MK_THROW( "Invalid node! " , textureNodes[i].ci , " " , textureNodes[i].cj );
		nodeIndex(textureNodes[i].ci, textureNodes[i].cj) = i;
	}

	BoundaryProlongationData< Real > boundaryProlongation;
	InitializeBoundaryProlongationData( hierarchy.gridAtlases[0] , boundaryProlongation );

	//////////////////////////////////// Initialize multigrid indices
	multigridIndices.resize( levels );
	for( int i=0 ; i<levels ; i++ )
	{
		const GridAtlas< PreReal , Real > &gridAtlas = hierarchy.gridAtlases[i];
		multigridIndices[i].threadTasks = gridAtlas.threadTasks;
		multigridIndices[i].boundaryGlobalIndex = gridAtlas.boundaryGlobalIndex;
		multigridIndices[i].segmentedLines = gridAtlas.segmentedLines;
		multigridIndices[i].rasterLines = gridAtlas.rasterLines;
		multigridIndices[i].restrictionLines = gridAtlas.restrictionLines;
		multigridIndices[i].prolongationLines = gridAtlas.prolongationLines;
		if( i<levels-1 ) multigridIndices[i].boundaryRestriction = hierarchy.boundaryRestriction[i];
	}

	//////////////////////////////////// Initialize multigrid coefficients

	//////////////////////////////////// 	Line Convolution coefficients
	{
		std::vector< Point2D< PreReal > > vectorField;
		if( InVectorField.set )
		{
			if( IntrinsicVectorField.set )
			{
				ReadVector( vectorField , InVectorField.value );
				if( vectorField.size()!=mesh.triangles.size() ) MK_THROW( "Triangle and vector counts don't match: " , mesh.triangles.size() , " != " , vectorField.size() );
			}
			else
			{
				std::vector< Point3D< PreReal > > _vectorField;
				ReadVector( _vectorField , InVectorField.value );
				if( _vectorField.size()!=mesh.triangles.size() ) MK_THROW( "Triangle and vector counts don't match: " , mesh.triangles.size() , " != " , _vectorField.size() );
				vectorField.resize( _vectorField.size() );
				ThreadPool::ParallelFor
					(
						0 , mesh.triangles.size() ,
						[&]( unsigned int , size_t i )
						{
							Point3D< PreReal > v[] = { mesh.vertices[ mesh.triangles[i][0] ] , mesh.vertices[ mesh.triangles[i][1] ] , mesh.vertices[ mesh.triangles[i][2] ] };
							Point3D< PreReal > d[] = { v[1]-v[0] , v[2]-v[0] };
							SquareMatrix< PreReal , 2 > Dot;
							for( int j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ ) Dot(j,k) = Point3D< PreReal >::Dot( d[j] , d[k] );
							Point2D< PreReal > dot( Point3D< PreReal >::Dot( d[0] , _vectorField[i] ) , Point3D< PreReal >::Dot( d[1] , _vectorField[i] ) );
							vectorField[i] = Dot.inverse() * dot;
						}
					);
			}
		}
		else
		{
			// Compute the principal curvatures
			std::vector< PrincipalCurvature< PreReal > > principalCurvatures;
			mesh.updateNormals();
			// Smooth the normals
			{
				Miscellany::Timer tmr;

				SparseMatrix< PreReal , int > M , _M = rMesh.template massMatrix< FEM::BASIS_0_WHITNEY >() , _S = rMesh.template stiffnessMatrix< FEM::BASIS_0_WHITNEY >();
				M.resize( 2*mesh.vertices.size() );
				ThreadPool::ParallelFor
					(
						0 , mesh.vertices.size() ,
						[&]( unsigned int , size_t i )
						{
							for( int ii=0 ; ii<2 ; ii++ )
							{
								M.SetRowSize( 2*i+ii , 2*_M.rowSizes[i] );
								for( int j=0 ; j<_M.rowSizes[i] ; j++ ) for( int jj=0 ; jj<2 ; jj++ ) M[2*i+ii][2*j+jj].N = _M[i][j].N*2+jj;
							}
						}
					);
				std::vector< Point3D< PreReal > > tangents( mesh.vertices.size()*2 );
				std::vector< PreReal > b( mesh.vertices.size()*2 ) , o( mesh.vertices.size()*2 );

				typedef EigenSolverCholeskyLDLt< PreReal , typename SparseMatrix< PreReal , int >::RowIterator > Solver;
				Solver solver( M , true );

				for( int iter=0 ; iter<NormalSmoothingIterations.value ; iter++ )
				{

					// Set the tangent directions
					ThreadPool::ParallelFor
					(
						0 , mesh.vertices.size() ,
						[&]( unsigned int , size_t i )
						{
							Point3D< PreReal > v( 1 , 0 , 0 );
							if( fabs( Point3D< PreReal >::Dot( v , mesh.normals[i] ) )>0.99 ) v = Point3D< PreReal >( 0 , 1 , 0 );
							tangents[2*i+0] = Point3D< PreReal >::CrossProduct( mesh.normals[i] , v               ) ; tangents[2*i+0] /= Point3D< PreReal >::Length( tangents[2*i+0] );
							tangents[2*i+1] = Point3D< PreReal >::CrossProduct( mesh.normals[i] , tangents[2*i+0] ) ; tangents[2*i+1] /= Point3D< PreReal >::Length( tangents[2*i+1] );
						}
					);

					// Solve for the tangent offsets minimizing the dirichlet energy:
					// E( o1 , o2 ) = || \sum o[i] * T[i] ||^2 + e * || \nabla( \sum n[i] + o[i] * T[i] ) ||^2
					//              = o^t * T^t * M * T * o + e * [ o^t * T^t * S * T * o + 2 * o^t * T^t * S * n + n^t * S * n ]
					// \nabla E = 0:
					// 0 = T^t * ( M + e * S ) * T * o + e * T^t * S * n
					{
						ThreadPool::ParallelFor
						(
							0 , mesh.vertices.size() ,
							[&]( unsigned int , size_t i )
							{
								for( int ii=0 ; ii<2 ; ii++ ) 
								{
									b[2*i+ii] = 0;
									for( int j=0 ; j<_M.rowSizes[i] ; j++ )
									{
										for( int jj=0 ; jj<2 ; jj++ ) M[2*i+ii][2*j+jj].Value = ( _M[i][j].Value*NormalSmoothingInterpolation.value + _S[i][j].Value ) * Point3D< Real >::Dot( tangents[2*i+ii] , tangents[ 2*_M[i][j].N+jj ] );
										b[2*i+ii] -= _S[i][j].Value * Point3D< Real >::Dot( mesh.normals[ _S[i][j].N ] , tangents[2*i+ii] );
									}
								}
							}
						);
					}
					{
						solver.update( M );
						solver.solve( GetPointer( b ) , GetPointer( o ) );

						ThreadPool::ParallelFor
						(
							0 , mesh.vertices.size() ,
							[&]( unsigned int , size_t i )
							{
								mesh.normals[i] += tangents[2*i+0] * o[2*i+0] + tangents[2*i+1] * o[2*i+1] , mesh.normals[i] /= Point3D< PreReal >::Length( mesh.normals[i] );
							}
						);
					}
				}		
				if( Verbose.set ) printf( "\tSmoothed normals: %.2f(s)\n" , tmr.elapsed() );
			}
			InitializePrincipalCurvatureDirection( mesh , mesh.normals , principalCurvatures );
			mesh.updateNormals();

			// Set the vector-field to the principal curvature direction times the umbilicity
			vectorField.resize( principalCurvatures.size() );
			ThreadPool::ParallelFor
				(
					0 , principalCurvatures.size() ,
					[&]( unsigned int , size_t t ){  vectorField[t] = principalCurvatures[t].dirs[ MinimalCurvature.set ? 0 : 1 ] * ( principalCurvatures[t].values[1] - principalCurvatures[t].values[0] ); }
				);
		}
		// Normalize the vector-field to have unit-norm
		{
		std::vector< SquareMatrix< PreReal , 2 > > embeddingMetric;
			InitializeEmbeddingMetric( mesh , true , embeddingMetric );
			{
				PreReal norm = 0 , area = 0;
				for( int t=0 ; t<embeddingMetric.size() ; t++ )
				{
					PreReal a = (PreReal)sqrt( embeddingMetric[t].determinant() ) / 2.;
					norm += Point2D< PreReal >::Dot( vectorField[t] , embeddingMetric[t]*vectorField[t] ) * a;
					area += a;
				}
				norm = sqrt( norm / area );
				for( int t=0 ; t<embeddingMetric.size() ; t++ ) vectorField[t] /= (Real)norm;
			}
		}

		if( OutVectorField.set )
		{
#if 1
			std::cerr << "[WARNING] Forcing extrinsic output" << std::endl;
			std::vector< Point3D< PreReal > > _vectorField( vectorField.size() );
			ThreadPool::ParallelFor
				(
					0 , mesh.triangles.size() ,
					[&]( unsigned int , size_t i )
					{
						Point3D< PreReal > v[] = { mesh.vertices[ mesh.triangles[i][0] ] , mesh.vertices[ mesh.triangles[i][1] ] , mesh.vertices[ mesh.triangles[i][2] ] };
						_vectorField[i] = (v[1]-v[0]) * vectorField[i][0] + (v[2]-v[0]) * vectorField[i][1];
					}
				);
			WriteVector( _vectorField , OutVectorField.value );
#else
			if( IntrinsicVectorField.set ) WriteVector( vectorField , OutVectorField.value );
			else
			{
				std::vector< Point3D< PreReal > > _vectorField( vectorField.size() );
				ThreadPool::ParallelFor
				(
					0 , mesh.triangles.size() ,
					[&]( unsigned int , size_t i )
					{
						Point3D< PreReal > v[] = { mesh.vertices[ mesh.triangles[i][0] ] , mesh.vertices[ mesh.triangles[i][1] ] , mesh.vertices[ mesh.triangles[i][2] ] };
						_vectorField[i] = (v[1]-v[0]) * vectorField[i][0] + (v[2]-v[0]) * vectorField[i][1];
					}
				);
				WriteVector( _vectorField , OutVectorField.value );
			}
#endif
		}
		{
			std::vector< FEM::SamplePoint< PreReal > > randomSamples = rMesh.randomSamples( 5e5 );
			visualization.vectorField.resize( randomSamples.size() );
			for( int i=0 ; i<randomSamples.size() ; i++ )
			{
				visualization.vectorField[i].tIdx = randomSamples[i].tIdx;
				visualization.vectorField[i].p = Point2D< float >( randomSamples[i].p );
				visualization.vectorField[i].v = Point2D< float >( vectorField[ randomSamples[i].tIdx ] );
			}
		}

		auto LengthToAnisotropy = [&]( PreReal len )
		{
			// g <- g + gOrtho * anisotropy 
			// 0 -> 0
			// 1 -> 1e5
			// infty -> infty
			return (PreReal)( pow( len , AnisotropyExponent.value ) * 1e5 );
		};
		InitializeAnisotropicMetric( mesh , atlasCharts , vectorField , LengthToAnisotropy , parameterMetric );

		std::vector< Point3D< Real > > __inputSignal;
		std::vector< Real > __texelToCellCoeffs;
		SparseMatrix< Real , int> __boundaryCellBasedStiffnessRHSMatrix[3];

		timer.reset();
		{
			switch( MatrixQuadrature.value )
			{
			case  1: InitializeMassAndStiffness< 1>( anisoMassCoefficients , anisoStiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case  3: InitializeMassAndStiffness< 3>( anisoMassCoefficients , anisoStiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case  6: InitializeMassAndStiffness< 6>( anisoMassCoefficients , anisoStiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case 12: InitializeMassAndStiffness<12>( anisoMassCoefficients , anisoStiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case 24: InitializeMassAndStiffness<24>( anisoMassCoefficients , anisoStiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case 32: InitializeMassAndStiffness<32>( anisoMassCoefficients , anisoStiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			default: MK_THROW( "Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles" );
			}
		}
		if( Verbose.set ) printf( "\tInitialized mass and stiffness: %.2f(s)\n" , timer.elapsed() );

		if( UseDirectSolver.set )
		{
			Miscellany::Timer tmr;
			FullMatrixConstruction( hierarchy.gridAtlases[0] , anisoMassCoefficients , anisotropicMass);
			FullMatrixConstruction( hierarchy.gridAtlases[0] , anisoStiffnessCoefficients , anisotropicStiffness);
			lineConvolutionMatrix = anisotropicMass * licInterpolationWeight + anisotropicStiffness;
			printf( "Assembling matrices =  %.4f\n" , tmr.elapsed() );
		}

		timer.reset();
		UpdateLinearSystem( licInterpolationWeight , (Real)1. , hierarchy , multigridLineConvolutionCoefficients , anisoMassCoefficients , anisoStiffnessCoefficients , lineConvolutionSolvers , fineLineConvolutionSolver , lineConvolutionMatrix , DetailVerbose.set , true , UseDirectSolver.set );
		if( Verbose.set ) printf( "\tInitialized multigrid coefficients: %.2f(s)\n" , timer.elapsed() );
	}

	//////////////////////////////////// 	Modulation coefficients
	SparseMatrix< Real, int > modMatrix;
	{
		InitializeMetric( mesh , EMBEDDING_METRIC , atlasCharts , parameterMetric );

		std::vector< Point3D< Real > > __inputSignal;
		std::vector< Real > __texelToCellCoeffs;
		SparseMatrix< Real , int > __boundaryCellBasedStiffnessRHSMatrix[3];

		timer.reset();
		{
			switch( MatrixQuadrature.value )
			{
			case  1: InitializeMassAndStiffness< 1>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case  3: InitializeMassAndStiffness< 3>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case  6: InitializeMassAndStiffness< 6>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case 12: InitializeMassAndStiffness<12>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case 24: InitializeMassAndStiffness<24>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case 32: InitializeMassAndStiffness<32>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			default: MK_THROW( "Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles" );
			}
		}
		if( Verbose.set ) printf( "\tInitialized mass and stiffness: %.2f(s)\n" , timer.elapsed() );

		if( UseDirectSolver.set )
		{
			Miscellany::Timer tmr;
			FullMatrixConstruction( hierarchy.gridAtlases[0] , massCoefficients , mass);
			FullMatrixConstruction( hierarchy.gridAtlases[0] , stiffnessCoefficients , stiffness);
			modMatrix = mass * sharpeningInterpolationWeight + stiffness;
			printf( "Assembling matrices =  %.4f\n" , tmr.elapsed() );
		}

		timer.reset();
		UpdateLinearSystem( sharpeningInterpolationWeight , (Real)1. , hierarchy , multigridModulationCoefficients , massCoefficients , stiffnessCoefficients , modulationSolvers , fineModulationSolver , modMatrix , DetailVerbose.set , true , UseDirectSolver.set );
		if( Verbose.set ) printf( "\tInitialized multigrid coefficients: %.2f(s)\n" , timer.elapsed() );
	}

	//////////////////////////////////// Initialize multigrid variables
	multigridLineConvolutionVariables.resize(levels);
	for (int i = 0; i < levels; i++) {
		MultigridLevelVariables< Point3D< Real > >& variables = multigridLineConvolutionVariables[i];
		variables.x.resize(hierarchy.gridAtlases[i].numTexels);
		variables.rhs.resize(hierarchy.gridAtlases[i].numTexels);
		variables.residual.resize(hierarchy.gridAtlases[i].numTexels);
		variables.boundary_rhs.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
		variables.boundary_value.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
		variables.variable_boundary_value.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
	}

	multigridModulationVariables.resize(levels);
	for (int i = 0; i < levels; i++) {
		MultigridLevelVariables< Point3D< Real > >& variables = multigridModulationVariables[i];
		variables.x.resize(hierarchy.gridAtlases[i].numTexels);
		variables.rhs.resize(hierarchy.gridAtlases[i].numTexels);
		variables.residual.resize(hierarchy.gridAtlases[i].numTexels);
		variables.boundary_rhs.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
		variables.boundary_value.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
		variables.variable_boundary_value.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
	}

	randSignal.resize( textureNodes.size() );

	for( int i=0 ; i<randSignal.size() ; i++ )
	{
		Point3D< float > randomColor = HSV2RGB( Random< float >() , 1.f , 1.f );
		randSignal[i] = Point3D< Real >( randomColor[0] , randomColor[1] , randomColor[2] );
	}
	Reset();
}

template< typename PreReal , typename  Real >
void LineConvolution< PreReal , Real >::Reset( void )
{
	for( int i=0 ; i<multigridLineConvolutionVariables[0].x.size() ; i++) multigridLineConvolutionVariables[0].x[i] *= 0;
	for( int i=0 ; i<multigridModulationVariables[0].x.size() ; i++) multigridModulationVariables[0].x[i] *= 0;

	if( UseDirectSolver.set ) ComputeExactSolution();
	else for( int i=0 ; i<multigridModulationVariables[0].x.size() ; i++) multigridModulationVariables[0].x[i] = randSignal[i];

	steps = 0;
}

template< typename PreReal , typename  Real >
void LineConvolution< PreReal , Real >::InitializeVisualization( void )
{
	int tCount = (int)mesh.triangles.size();

	visualization.triangles.resize( tCount );
	visualization.vertices.resize( 3*tCount );
	visualization.colors.resize( 3*tCount , Point3D< float >( 0.75f , 0.75f , 0.75f ) );
	visualization.textureCoordinates.resize( 3*tCount );
	visualization.normals.resize( 3*tCount );


	for( int i=0 ; i<tCount ; i++ ) for( int k=0 ; k<3 ; k++ ) visualization.triangles[i][k] = 3*i+k;

	for (int i = 0; i<tCount; i++) {
		for (int j = 0; j < 3; j++) {
			visualization.vertices[3 * i + j] = mesh.vertices[mesh.triangles[i][j]];
			visualization.normals[3 * i + j] = mesh.normals[mesh.triangles[i][j]];
			visualization.textureCoordinates[3 * i + j] = mesh.textureCoordinates[3 * i + j];
		}
	}



	std::vector<int> boundaryEdges;
	mesh.initializeBoundaryEdges( boundaryEdges );

	for (int e = 0; e < boundaryEdges.size(); e++) {
		int tIndex = boundaryEdges[e] / 3;
		int kIndex = boundaryEdges[e] % 3;
		for (int c = 0; c < 2; c++)
		{
			Point3D< float > v = Point3D< float >( mesh.vertices[mesh.triangles[tIndex][(kIndex + c) % 3]] );
			visualization.boundaryEdgeVertices.push_back(v);
		}
	}

	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'f' , "lic interpolation weight" , "LIC Interpolation Weight" , LICInterpolationWeightCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'g' , "sharpening gradient modulation" , "Sharpening Gradient Modulation" , SharpeningGradientModulationCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'y' , "sharpening interpolation weight" , "Sharpening Interpolation Weight" , SharpeningInterpolationWeightCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 's' , "export texture" , "Output Texture" , ExportTextureCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , ' ' , "toggle update" , ToggleUpdateCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , '+' , "increment update" , IncrementUpdateCallBack ) );
	visualization.info.push_back( stepsString );

	visualization.UpdateVertexBuffer();
	visualization.UpdateFaceBuffer();
	visualization.UpdateTextureBuffer();

	UpdateOutputBuffer( multigridModulationVariables[0].x );
}
template< typename PreReal , typename Real >
void LineConvolution< PreReal , Real >::Init( void )
{
	sprintf( stepsString , "Steps: 0" );
	levels = Levels.value;
	textureWidth = Width.value;
	textureHeight = Height.value;
	sharpeningGradientModulation = SharpeningGradientModulation.value;
	sharpeningInterpolationWeight = SharpeningInterpolationWeight.value;
	licInterpolationWeight = LICInterpolationWeight.value;

	mesh.read( Input.value , DetailVerbose.set , CollapseEpsilon.value );

	if( RandomJitter.set )
	{
		if( RandomJitter.value ) srand( RandomJitter.value );
		else                     srand( time(NULL) );
		std::vector< Point2D< PreReal > >randomOffset( mesh.vertices.size() );
		PreReal jitterScale = (PreReal)1e-3 / std::max< int >( textureWidth , textureHeight );
		for( int i=0 ; i<randomOffset.size() ; i++ ) randomOffset[i] = Point2D< PreReal >( (PreReal)1. - Random< PreReal >()*2 , (PreReal)1. - Random< PreReal >()*2 ) * jitterScale;
		for( int i=0 ; i<mesh.triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ ) mesh.textureCoordinates[ 3*i+k ] += randomOffset[ mesh.triangles[i][k] ];
	}

	{
		padding = Padding::Init( textureWidth , textureHeight , mesh.textureCoordinates , DetailVerbose.set );
		padding.pad( textureWidth , textureHeight , mesh.textureCoordinates );
		textureWidth  += padding.width();
		textureHeight += padding.height();
	}

	// Define centroid and scale for visualization
	Point3D< PreReal > centroid;
	for( int i=0 ; i<mesh.vertices.size() ; i++ ) centroid += mesh.vertices[i];
	centroid /= (int)mesh.vertices.size();
	PreReal radius = 0;
	for( int i=0 ; i<mesh.vertices.size() ; i++ ) radius = std::max< PreReal >( radius , Point3D< PreReal >::Length( mesh.vertices[i]-centroid ) );
	for( int i=0 ; i<mesh.vertices.size() ; i++ ) mesh.vertices[i] = ( mesh.vertices[i]-centroid ) / radius;

	Miscellany::Timer timer;
	FEM::RiemannianMesh< PreReal , unsigned int > rMesh( GetPointer( mesh.triangles ) , mesh.triangles.size() );
	rMesh.setMetricFromEmbedding( GetPointer( mesh.vertices ) );
	rMesh.makeUnitArea();

	InitializeSystem( rMesh , textureWidth , textureHeight );
	if( Verbose.set )
	{
		printf( "Resolution: %d / %d x %d\n" , (int)textureNodes.size() , textureWidth , textureHeight );
		printf( "Initialized system %.2f(s)\n" , timer.elapsed() );
		printf( "Peak Memory (MB): %d\n" , Miscellany::MemoryInfo::PeakMemoryUsageMB() );
	}

	//Assign position to exterior nodes using barycentric-exponential map
	{
		Pointer( FEM::CoordinateXForm< PreReal > ) xForms = rMesh.getCoordinateXForms();

		for( int i=0 ; i<textureNodes.size() ; i++ ) if( textureNodes[i].tID!=-1 && !textureNodes[i].isInterior )
		{
			FEM::HermiteSamplePoint< PreReal > _p;
			_p.tIdx = textureNodes[i].tID;
			_p.p = Point2D< PreReal >( (PreReal)1./3 , (PreReal)1./3 );
			_p.v = textureNodes[i].barycentricCoords - _p.p;

			rMesh.exp(xForms, _p);

			textureNodes[i].tID = _p.tIdx;
			textureNodes[i].barycentricCoords = _p.p;
		}
		DeletePointer( xForms );		
	}

	outputBuffer = new unsigned char[ textureHeight*textureWidth*3 ];
	memset( outputBuffer , 204 , textureHeight*textureWidth*3*sizeof( unsigned char ) );
}

template< typename PreReal , typename Real>
void _main( int argc , char* argv[] )
{
	LineConvolution< PreReal , Real >::Init();
	if( !Output.set )
	{
		glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
		LineConvolution< PreReal , Real >::visualization.displayMode = DisplayMode.value;
		if     ( DisplayMode.value==ONE_REGION_DISPLAY ) LineConvolution< PreReal , Real >::visualization.screenWidth =  800 , LineConvolution< PreReal , Real >::visualization.screenHeight = 800;
		else if( DisplayMode.value==TWO_REGION_DISPLAY ) LineConvolution< PreReal , Real >::visualization.screenWidth = 1440 , LineConvolution< PreReal , Real >::visualization.screenHeight = 720;
		glutInitWindowSize( LineConvolution< PreReal , Real >::visualization.screenWidth , LineConvolution< PreReal , Real >::visualization.screenHeight );
		glutInit( &argc , argv );
		char windowName[1024];
		sprintf( windowName , "Line Integral Convolution" );
		glutCreateWindow( windowName );
		if( glewInit()!=GLEW_OK ) MK_THROW( "glewInit failed" );
		glutDisplayFunc ( LineConvolution< PreReal , Real >::Display );
		glutReshapeFunc ( LineConvolution< PreReal , Real >::Reshape );
		glutMouseFunc   ( LineConvolution< PreReal , Real >::MouseFunc );
		glutMotionFunc  ( LineConvolution< PreReal , Real >::MotionFunc );
		glutKeyboardFunc( LineConvolution< PreReal , Real >::KeyboardFunc );
		if( !UseDirectSolver.set ) glutIdleFunc( LineConvolution< PreReal , Real >::Idle );
		if( CameraConfig.set ) LineConvolution< PreReal , Real >::visualization.ReadSceneConfigurationCallBack( &LineConvolution< PreReal , Real >::visualization , CameraConfig.value.c_str() );
		LineConvolution< PreReal , Real >::InitializeVisualization();
		glutMainLoop();
	}
	else
	{
		if( UseDirectSolver.set ) LineConvolution< PreReal , Real >::ComputeExactSolution();
		else for( int i=0 ; i<OutputVCycles.value ; i++ ) LineConvolution< PreReal , Real >::UpdateSolution();
		LineConvolution< PreReal , Real >::SetOutputBuffer( LineConvolution< PreReal , Real >::multigridModulationVariables[0].x );
		LineConvolution< PreReal , Real >::ExportTextureCallBack( &LineConvolution< PreReal , Real >::visualization , Output.value.c_str() );
	}
}

int main( int argc , char *argv[] )
{
	CmdLineParse( argc-1 , argv+1 , params );
	if( !Input.set ) { ShowUsage(argv[0]); return EXIT_FAILURE; }
	if( Serial.set ) ThreadPool::ParallelizationType = ThreadPool::ParallelType::NONE;
	if( !NoHelp.set && !Output.set )
	{
		printf( "+-----------------------------------------------+\n" );
		printf( "| Interface Controls:                           |\n" );
		printf( "|    [Left Mouse]:                 rotate       |\n" );
		printf( "|    [Right Mouse]:                zoom         |\n" );
		printf( "|    [Left/Right Mouse] + [CTRL]:  pan          |\n" );
		printf( "|    [SPACE]:                      start solver |\n" );
		printf( "+-----------------------------------------------+\n" );
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
