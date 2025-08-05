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

#ifdef USE_EIGEN_PARDISO
#include <Eigen/PardisoSupport>
#endif // USE_EIGEN_PARDISO

#include <Misha/CmdLineParser.h> 
#include <Misha/Miscellany.h>
#include <Misha/Exceptions.h>
#include <Misha/FEM.h>
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

const float StripeRates[] = { 0.062f , 0.062f };
const float    DotRates[] = { 0.0367f , 0.0649f };
const float DefaultStripesSamplesFraction = 0.01f;
const float DefaultDotsSamplesFraction = 0.001f;

CmdLineParameter< std::string >
	Input( "in" ) ,
	Output( "out" ) ,
	CameraConfig( "camera" ) ,
	VectorField( "inVF" );


CmdLineParameter< unsigned int >
	OutputSteps( "outSteps" , 1000 ) ,
	Width( "width" , 512 ) ,
	Height( "height" , 512 ) ,
	DisplayMode( "display" , TWO_REGION_DISPLAY ) ,
	MatrixQuadrature( "mQuadrature" , 6 ) ,
	RHSQuadrature( "rhsQuadrature" , 3 ) ,
	MultigridBlockHeight ( "mBlockH" ,  16 ) ,
	MultigridBlockWidth  ( "mBlockW" , 128 ) ,
	MultigridPaddedHeight( "mPadH"   ,   0 ) ,
	MultigridPaddedWidth ( "mPadW"   ,   2 ) ,
	RandomJitter( "jitter" , 0 ) ,
	Levels( "levels" , 4 );

CmdLineParameterArray< float , 2 >
	FeedKillRates( "fk" , StripeRates );

CmdLineParameter< float >
	Speed( "speed" , 10.f ) ,
	DiffusionScale( "diff" , 1.f ) ,
	SamplesFraction( "samples" ) ,
	AnisotropyScale( "aScl" , 1.f ) ,
	AnisotropyExponent( "aExp" , 0.f );

CmdLineReadable
	Verbose( "verbose" ) , 
	NearestSampling( "nearest" ) , 
	NoHelp( "noHelp" ) ,
	DetailVerbose( "detail" ) ,
	UseDirectSolver( "useDirectSolver" ) ,
	Double( "double" ) ,
	ApproximateIntegration( "approximateIntegration" ) ,
	Dots( "dots" ) ,
	Serial( "serial" ) ,
	Run( "run" ) ,
	IntrinsicVectorField( "intrinsicVF" );

CmdLineParameter< double >
	CollapseEpsilon( "collapse" , 0 );

CmdLineReadable* params[] =
{
	&Input , &Output , &OutputSteps , &Width , &Height , &Speed , &FeedKillRates , &DiffusionScale , &SamplesFraction , &CameraConfig , &Levels , &UseDirectSolver , &Serial , &DisplayMode , &MultigridBlockHeight , &MultigridBlockWidth , &MultigridPaddedHeight , &MultigridPaddedWidth ,
	&Verbose , &DetailVerbose ,
	&RandomJitter ,
	&Double ,
	&MatrixQuadrature , &RHSQuadrature ,
	&ApproximateIntegration , &Dots ,
	&NoHelp ,
	&VectorField , &IntrinsicVectorField , &AnisotropyScale , &AnisotropyExponent , 
	&CollapseEpsilon ,
	&Run ,
	&NearestSampling ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n", ex );
	printf( "\t --%s <input mesh>\n" , Input.name.c_str() );
	printf( "\t[--%s <output texture>\n" , Output.name.c_str() );
	printf( "\t[--%s <output steps>=%d]\n" , OutputSteps.name.c_str() , OutputSteps.value );
	printf( "\t[--%s <texture width>=%d]\n" , Width.name.c_str() , Width.value );
	printf( "\t[--%s <texture height>=%d]\n" , Height.name.c_str() , Height.value );
	printf( "\t[--%s <time-step>=%f]\n" , Speed.name.c_str() , Speed.value );
	printf( "\t[--%s <feed/kill rates>=%f %f]\n" , FeedKillRates.name.c_str() , FeedKillRates.values[0] , FeedKillRates.values[1] );
	printf( "\t[--%s <diffusion scale>=%f]\n" , DiffusionScale.name.c_str() , DiffusionScale.value );
	printf( "\t[--%s <samples fraction>=%f / %f]\n" , SamplesFraction.name.c_str() , DefaultStripesSamplesFraction , DefaultDotsSamplesFraction );
	printf( "\t[--%s <system matrix quadrature points per triangle>=%d]\n" , MatrixQuadrature.name.c_str() , MatrixQuadrature.value );
	printf( "\t[--%s <right-hand-side quadrature points per triangle>=%d]\n" , RHSQuadrature.name.c_str() , RHSQuadrature.value );
	printf( "\t[--%s]\n" , ApproximateIntegration.name.c_str() );
	printf( "\t[--%s]\n" , UseDirectSolver.name.c_str() );
	printf( "\t[--%s <jittering seed>]\n" , RandomJitter.name.c_str() );
	printf( "\t[--%s]\n" , Dots.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );

	printf( "\t[--%s <camera configuration file>]\n" , CameraConfig.name.c_str() );
	printf( "\t[--%s <hierarchy levels>=%d]\n" , Levels.name.c_str() , Levels.value );
	printf( "\t[--%s]\n" , DetailVerbose.name.c_str() );
	printf( "\t[--%s <display mode>=%d]\n" , DisplayMode.name.c_str() , DisplayMode.value );
	printf( "\t\t%d] One Region \n" , ONE_REGION_DISPLAY );
	printf( "\t\t%d] Two Region \n" , TWO_REGION_DISPLAY );
	printf( "\t[--%s <multigrid block width>=%d]\n"   , MultigridBlockWidth.name.c_str()   , MultigridBlockWidth.value   );
	printf( "\t[--%s <multigrid block height>=%d]\n"  , MultigridBlockHeight.name.c_str()  , MultigridBlockHeight.value  );
	printf( "\t[--%s <multigrid padded width>=%d]\n"  , MultigridPaddedWidth.name.c_str()  , MultigridPaddedWidth.value  );
	printf( "\t[--%s <multigrid padded height>=%d]\n" , MultigridPaddedHeight.name.c_str() , MultigridPaddedHeight.value );

	printf( "\t[--%s <input vector field>]\n" , VectorField.name.c_str() );
	printf( "\t[--%s <anisotropy scale>=%f]\n" , AnisotropyScale.name.c_str() , AnisotropyScale.value );
	printf( "\t[--%s <anisotropy exponent>=%f]\n" , AnisotropyExponent.name.c_str() , AnisotropyExponent.value );
	printf( "\t[--%s <collapse epsilon>=%g]\n" , CollapseEpsilon.name.c_str() , CollapseEpsilon.value );
	printf( "\t[--%s]\n" , IntrinsicVectorField.name.c_str() );
	printf( "\t[--%s]\n" , Run.name.c_str() );
	printf( "\t[--%s]\n" , Serial.name.c_str() );
	printf( "\t[--%s]\n" , NearestSampling.name.c_str() );

	printf( "\t[--%s]\n" , NoHelp.name.c_str() );
}

template< typename PreReal , typename Real >
class GrayScottReactionDiffusion
{
public:
	static std::vector< FEM::SamplePoint< PreReal > > randomSamples;
	static TexturedTriangleMesh< PreReal > mesh;
	static int textureWidth;
	static int textureHeight;
	static Real diffusionRates[2];
	static Real kill;
	static Real feed;
	static Real speed;
	static unsigned int levels;
	static int steps;
	static char stepsString[];
	static int whichConcentration;

	static Padding padding;

	static std::vector< Point3D< float > > textureNodePositions;

	static HierarchicalSystem< PreReal , Real > hierarchy;

	static std::vector< TextureNodeInfo< PreReal > > textureNodes;
	static RegularGrid< 2 , int > nodeIndex;

	static SparseMatrix< Real , int > mass;
	static SparseMatrix< Real , int > stiffness;
	static SparseMatrix< Real , int > systemMatrices[2];

	static unsigned int seedTexel;

	static std::vector< SystemCoefficients< Real > > multigridCoefficients[2];
	static std::vector< MultigridLevelVariables< Real > > multigridVariables[2];

#ifdef USE_EIGEN_PARDISO
	using EigenSolver = Eigen::PardisoLDLT< Eigen::SparseMatrix< double > >;
#else // !USE_EIGEN_PARDISO
	using EigenSolver = Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > >;
#endif // USE_EIGEN_PARDISO
	static VCycleSolvers< EigenSolver > vCycleSolvers[2];
	static EigenSolverWrapper< EigenSolver > fineSolvers[2];

	static std::vector< MultigridLevelIndices< Real > > multigridIndices;

	static ScalarIntegrator< Real > scalarIntegrator;
	static typename ScalarIntegrator< Real >::template Scratch< Point< Real , 2 > , Point< Real , 2 > > scalarIntegratorScratch;

	// Linear Operators
	static MassAndStiffnessOperators< Real > massAndStiffnessOperators;

	static unsigned char * outputBuffer;

	//Visulization
	static TexturedMeshVisualization visualization;
	static int mouseX , mouseY;
	static bool mouseSelectionActive;

	static void SetOutputBuffer( const std::vector< Real >& solution );
	static void UpdateOutputBuffer( const std::vector< Real >& solution );
	static void UpdateOutputBuffer( void );

	static int updateCount;

	static void SetConcentration1CallBack( Visualization* v , const char* prompt );
	static void SetConcentration2CallBack( Visualization* v , const char* prompt );
	static void HideConcentrationsCallBack( Visualization *v , const char *prompt );
	static void ToggleUpdateCallBack( Visualization* v , const char* prompt );
	static void IncrementUpdateCallBack( Visualization* v , const char* prompt );
	static void ExportTextureCallBack( Visualization* v , const char* prompt );
	static void Init( void );
	static void InitializeVisualization( void );
	static void SetRightHandSide( void );
	static void UpdateExactSolution( bool verbose=false );
	static void UpdateApproximateSolution( bool verbose=false , bool detailVerbose=false );
	static void InitializeSystem( int width , int height );
	static void InitializeConcentrations( void );

	static void Display( void ){ visualization.Display(); }
	static void MouseFunc( int button , int state , int x , int y );
	static void MotionFunc( int x , int y );
	static void Reshape( int w , int h ){ visualization.Reshape(w,h); }
	static void KeyboardFunc( unsigned char key , int x , int y ){ visualization.KeyboardFunc( key , x , y ); }
	static void Idle( void );
};

template< typename PreReal , typename Real > std::vector< FEM::SamplePoint< PreReal > >								GrayScottReactionDiffusion< PreReal , Real >::randomSamples;
template< typename PreReal , typename Real > TexturedTriangleMesh< PreReal >										GrayScottReactionDiffusion< PreReal , Real >::mesh;
template< typename PreReal , typename Real > int																	GrayScottReactionDiffusion< PreReal , Real >::textureWidth;
template< typename PreReal , typename Real > int																	GrayScottReactionDiffusion< PreReal , Real >::textureHeight;
template< typename PreReal , typename Real > TexturedMeshVisualization												GrayScottReactionDiffusion< PreReal , Real >::visualization;
template< typename PreReal , typename Real > int																	GrayScottReactionDiffusion< PreReal , Real >::mouseX = -1;
template< typename PreReal , typename Real > int																	GrayScottReactionDiffusion< PreReal , Real >::mouseY = -1;
template< typename PreReal , typename Real > bool																	GrayScottReactionDiffusion< PreReal , Real >::mouseSelectionActive = false;
template< typename PreReal , typename Real > Padding																GrayScottReactionDiffusion< PreReal , Real >::padding;

template< typename PreReal , typename Real > SparseMatrix< Real , int >												GrayScottReactionDiffusion< PreReal , Real >::mass;
template< typename PreReal , typename Real > SparseMatrix< Real , int >												GrayScottReactionDiffusion< PreReal , Real >::stiffness;
template< typename PreReal , typename Real > SparseMatrix< Real , int >												GrayScottReactionDiffusion< PreReal , Real >::systemMatrices[2];

template< typename PreReal , typename Real > Real																	GrayScottReactionDiffusion< PreReal , Real >::diffusionRates[2];
template< typename PreReal , typename Real > Real																	GrayScottReactionDiffusion< PreReal , Real >::speed;
template< typename PreReal , typename Real > Real																	GrayScottReactionDiffusion< PreReal , Real >::kill;
template< typename PreReal , typename Real > Real																	GrayScottReactionDiffusion< PreReal , Real >::feed;

template< typename PreReal , typename Real > std::vector< TextureNodeInfo< PreReal > >								GrayScottReactionDiffusion< PreReal , Real >::textureNodes;
template< typename PreReal , typename Real > RegularGrid< 2 , int >													GrayScottReactionDiffusion< PreReal , Real >::nodeIndex;

template< typename PreReal , typename Real > int																	GrayScottReactionDiffusion< PreReal , Real >::steps;
template< typename PreReal , typename Real > char																	GrayScottReactionDiffusion< PreReal , Real >::stepsString[1024];
template< typename PreReal , typename Real > unsigned int															GrayScottReactionDiffusion< PreReal , Real >::levels;
template< typename PreReal , typename Real > HierarchicalSystem< PreReal , Real >									GrayScottReactionDiffusion< PreReal , Real >::hierarchy;

template< typename PreReal , typename Real > unsigned char *														GrayScottReactionDiffusion< PreReal , Real >::outputBuffer;
template< typename PreReal , typename Real > std::vector< MultigridLevelIndices< Real > >							GrayScottReactionDiffusion< PreReal , Real >::multigridIndices;

template< typename PreReal , typename Real > std::vector< SystemCoefficients< Real > >								GrayScottReactionDiffusion< PreReal , Real >::multigridCoefficients[2];
template< typename PreReal , typename Real > std::vector< MultigridLevelVariables< Real > >							GrayScottReactionDiffusion< PreReal , Real >::multigridVariables[2];
template< typename PreReal , typename Real > VCycleSolvers< typename GrayScottReactionDiffusion< PreReal , Real >::EigenSolver >		GrayScottReactionDiffusion< PreReal , Real >::vCycleSolvers[2];
template< typename PreReal , typename Real > EigenSolverWrapper< typename GrayScottReactionDiffusion< PreReal , Real >::EigenSolver >	GrayScottReactionDiffusion< PreReal , Real >::fineSolvers[2];



//Samples

template< typename PreReal , typename Real > unsigned int															GrayScottReactionDiffusion< PreReal , Real >::seedTexel = -1;
template< typename PreReal , typename Real > std::vector< Point3D< float > >										GrayScottReactionDiffusion< PreReal , Real >::textureNodePositions;

template< typename PreReal , typename Real > ScalarIntegrator< Real >												GrayScottReactionDiffusion< PreReal , Real >::scalarIntegrator;
template< typename PreReal , typename Real > typename ScalarIntegrator< Real >::template Scratch< Point< Real , 2 > , Point< Real , 2 > > GrayScottReactionDiffusion< PreReal , Real >::scalarIntegratorScratch;

template< typename PreReal , typename Real > MassAndStiffnessOperators< Real >										GrayScottReactionDiffusion< PreReal , Real >::massAndStiffnessOperators;

template< typename PreReal , typename Real > int																	GrayScottReactionDiffusion< PreReal , Real >::whichConcentration = 1;
template< typename PreReal , typename Real > int																	GrayScottReactionDiffusion< PreReal , Real >::updateCount = 0;


template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::SetRightHandSide( void )
{
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// ( D + d * diff1 * S ) a[t+d] = D( a[t] + d * ( - a[t] * b[t] * b[t] + feed * ( 1 - a[t] ) ) )    //
	// ( D + d * diff2 * S ) b[t+d] = D( b[t] + d * (   a[t] * b[t] * b[t] - ( kill + feed ) * b[t] ) ) //
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	// [MK] Should pull this out and initialize it once
	static std::vector< Point2D< Real > > ab_x( multigridVariables[0][0].x.size() );
	static std::vector< Point2D< Real > > ab_rhs( multigridVariables[0][0].rhs.size() );
	for( unsigned int ab=0 ; ab<2 ; ab++ ) ThreadPool::ParallelFor( 0 , ab_x.size() , [&]( unsigned int , size_t i ){ ab_x[i][ab] = multigridVariables[ab][0].x[i]; } );
	for( unsigned int ab=0 ; ab<2 ; ab++ ) massAndStiffnessOperators.mass( multigridVariables[ab][0].x , multigridVariables[ab][0].rhs );

	auto ABFunction = [&]( Point2D< Real > ab , SquareMatrix< Real , 2 > )
		{
			return Point2D< Real >
				(
					(Real)( speed * ( - ab[0] * ab[1] * ab[1] + feed * ( 1 - ab[0] ) ) ) ,
					(Real)( speed * (   ab[0] * ab[1] * ab[1] - ( kill + feed ) * ab[1] ) )
				);
		};

	scalarIntegrator( ab_x , ABFunction , scalarIntegratorScratch , ab_rhs );
	for( unsigned int ab=0 ; ab<2 ; ab++ ) ThreadPool::ParallelFor( 0 , ab_rhs.size() , [&]( size_t i ){ multigridVariables[ab][0].rhs[i] += ab_rhs[i][ab]; } );
}

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::UpdateExactSolution( bool verbose )
{
	Miscellany::PerformanceMeter pMeter( '.' );

	// (1) Compute the right-hand-sides
	{
		SetRightHandSide();
		if( verbose ) std::cout << pMeter( "Integrating" ) << std::endl;
	}

	// (2) Solve the linear systems
	{
		for( int ab=0 ; ab<2 ; ab++ )
		{
			fineSolvers[ab].solve( multigridVariables[ab][0].x , multigridVariables[ab][0].rhs );
			for( int i=0 ; i<multigridVariables[ab][0].x.size() ; i++ ) multigridVariables[ab][0].x[i] = std::max< Real >( multigridVariables[ab][0].x[i] , 0 );
		}
		if( verbose ) std::cout << pMeter( "Direct solve" ) << std::endl;
	}
}

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::UpdateApproximateSolution( bool verbose , bool detailVerbose )
{
	Miscellany::PerformanceMeter pMeter( '.' );

	// Compute the right-hand-sides
	{
		SetRightHandSide();
		if( verbose ) std::cout << pMeter( "Integrated" ) << std::endl;
	}

	// Solve the linear systems
	{
		for( int ab=0 ; ab<2 ; ab++ )
		{
			VCycle( multigridVariables[ab] , multigridCoefficients[ab] , multigridIndices , vCycleSolvers[ab] , 2 , detailVerbose , detailVerbose );
			ThreadPool::ParallelFor( 0 , multigridVariables[ab][0].x.size() , [&]( unsigned int , size_t i ){ multigridVariables[ab][0].x[i] = std::max< Real >( multigridVariables[ab][0].x[i] , 0 ); } );
		}
		if( verbose ) std::cout << pMeter( "V-Cycle solve" ) << std::endl;
	}
}

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::SetOutputBuffer( const std::vector< Real > & solution )
{
	ThreadPool::ParallelFor
		(
			0 , textureNodes.size() ,
			[&]( unsigned int , size_t i )
			{
				int ci = textureNodes[i].ci;
				int cj = textureNodes[i].cj;
				int offset = textureWidth*cj + ci;
				float value = (float)solution[i];
				value = std::max< float >( 0.f , std::min< Real >( 1.f , value*2.f ) );
				value = 1.f - value;
				value *= 0.75f;

				unsigned char color = (unsigned char)( value*255.f );
				outputBuffer[offset] = color;
			}
		);
}

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::UpdateOutputBuffer( const std::vector< Real > & solution )
{
	SetOutputBuffer( solution );
	glBindTexture( GL_TEXTURE_2D , visualization.textureBuffer );
	glTexImage2D( GL_TEXTURE_2D , 0 , GL_RGBA , textureWidth , textureHeight , 0 , GL_LUMINANCE , GL_UNSIGNED_BYTE , (GLvoid*)&outputBuffer[0] );
	glBindTexture( GL_TEXTURE_2D , 0 );

	glutPostRedisplay();
}

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::UpdateOutputBuffer( void )
{
	std::vector< Real > solution( textureNodes.size() , 0 );
	UpdateOutputBuffer( solution );
}

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::Idle( void )
{
	visualization.Idle();
	if( updateCount && !visualization.promptCallBack )
	{
		if( UseDirectSolver.set ) UpdateExactSolution();
		else UpdateApproximateSolution();

		if( updateCount>0 ) updateCount--;
		steps++;
		sprintf( stepsString , "Steps: %d" , steps );
	}
	if( whichConcentration<0 ) UpdateOutputBuffer();
	else UpdateOutputBuffer( multigridVariables[whichConcentration][0].x );
}

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::MouseFunc( int button , int state , int x , int y )
{
	visualization.newX = x , visualization.newY = y;
	visualization.rotating = visualization.scaling = visualization.panning = false;

	if( state==GLUT_DOWN && glutGetModifiers() & GLUT_ACTIVE_SHIFT )
	{
		seedTexel = -1;
		if( visualization.showMesh )
		{
			Point3D< float > selectedPoint;
			if( visualization.select( x , y , selectedPoint ) )
			{
				mouseSelectionActive = true;
				mouseX = x , mouseY = y;
				float minDistance = FLT_MAX;
				for( int i=0 ; i<textureNodePositions.size() ; i++ )
				{
					float squaredDistance = Point3D< float >::SquareNorm( textureNodePositions[i]-selectedPoint );
					if( squaredDistance<minDistance ) minDistance = squaredDistance , seedTexel = i;
				}
			}
		}
		else
		{
			Point2D< float > ip = visualization.selectImagePos( x , y );
			int i = floor( ip[0] * float( nodeIndex.res(0)) - 0.5f );
			int j = floor( (1.0-ip[1])*float( nodeIndex.res(1) ) - 0.5f );
			if( i>=0 && i<(int)nodeIndex.res(0) && j>=0 && j<(int)nodeIndex.res(1) ) mouseSelectionActive = true , seedTexel = nodeIndex(i,j);
		}
		InitializeConcentrations();
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
void GrayScottReactionDiffusion< PreReal , Real >::MotionFunc( int x , int y )
{
	if( mouseSelectionActive )
	{
		mouseX = x , mouseY = y;
		glutPostRedisplay();
	}
	else
	{
		if( !visualization.showMesh )
		{
			visualization.oldX = visualization.newX , visualization.oldY = visualization.newY , visualization.newX = x , visualization.newY = y;
			if( visualization.panning ) visualization.xForm.offset[0] -= ( visualization.newX - visualization.oldX ) / visualization.imageToScreenScale() , visualization.xForm.offset[1] += ( visualization.newY - visualization.oldY ) / visualization.imageToScreenScale();
			else
			{
				float dz = (float)pow( 1.1 , (float)( visualization.newY - visualization.oldY ) / 8.f );
				visualization.xForm.zoom *= dz;
			}
		}
		else
		{
			visualization.oldX = visualization.newX , visualization.oldY = visualization.newY , visualization.newX = x , visualization.newY = y;
			int screenSize = std::min< int >( visualization.screenWidth , visualization.screenHeight );
			float rel_x = ( visualization.newX - visualization.oldX ) / (float)screenSize * 2;
			float rel_y = ( visualization.newY - visualization.oldY ) / (float)screenSize * 2;

			float pRight = rel_x * visualization.zoom , pUp = -rel_y * visualization.zoom;
			float pForward = rel_y * visualization.zoom;
			float rRight = -rel_y , rUp = -rel_x;

			if     ( visualization.rotating ) visualization.camera.rotateUp( -rUp ) , visualization.camera.rotateRight( -rRight );
			else if( visualization.scaling  ) visualization.camera.translate( visualization.camera.forward*pForward);
			else if( visualization.panning  ) visualization.camera.translate( -( visualization.camera.right*pRight + visualization.camera.up*pUp ) );
		}
		glutPostRedisplay();
	}
}

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::ExportTextureCallBack( Visualization * /*v*/ , const char* prompt )
{
	RegularGrid< 2 , Point3D< Real > > outputImage;
	outputImage.resize( textureWidth , textureHeight );
	for( int i=0 ; i<outputImage.size() ; i++ ) outputImage[i] = Point3D< Real >( outputBuffer[i] , outputBuffer[i] , outputBuffer[i] ) / (Real)255.;
	padding.unpad( outputImage );
	WriteImage< 8 >( outputImage , prompt );
}

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::InitializeConcentrations( void )
{
	steps = 0;
	for( int i=0 ; i<multigridVariables[0][0].x.size() ; i++ ) multigridVariables[0][0].x[i] = 1;
	for( int i=0 ; i<multigridVariables[1][0].x.size() ; i++ ) multigridVariables[1][0].x[i] = 0;
	if( seedTexel!=-1 ) multigridVariables[1][0].x[seedTexel] = 1;
	else
	{
		for( int i=0 ; i<randomSamples.size() ; i++ )
		{
			int tIdx = randomSamples[i].tIdx;
			Point2D< Real > p = randomSamples[i].p;
			Point2D< Real > t = mesh.textureTriangle( tIdx )( Point2D< PreReal >( p ) );
			t[0] *= nodeIndex.res(0) , t[1] *= nodeIndex.res(1);
			int idx = nodeIndex( (int)floor( t[0] +0.5 ) , (int)floor( t[1] +0.5 ) );
			if( idx>=0 && idx<multigridVariables[1][0].x.size() ) multigridVariables[1][0].x[idx] = 1;
			else MK_WARN( "Bad random texel: " , t[0] , " " , t[1] , ": " , idx );
		}
	}
}

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::InitializeSystem( int width , int height )
{
	Miscellany::PerformanceMeter pMeter( '.' );
	ExplicitIndexVector< ChartIndex , AtlasChart< PreReal > > atlasCharts;
	ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< PreReal , 2 > > > parameterMetric;
	MultigridBlockInfo multigridBlockInfo( MultigridBlockWidth.value , MultigridBlockHeight.value , MultigridPaddedWidth.value , MultigridPaddedHeight.value );
	InitializeHierarchy( mesh , width , height , levels , textureNodes , hierarchy , atlasCharts , multigridBlockInfo );
	if( Verbose.set ) std::cout << pMeter( "Hierarchy" ) << std::endl;

	//Initialize node index
	nodeIndex.resize( width , height );
	for( int i=0 ; i<nodeIndex.size() ; i++ ) nodeIndex[i] = -1;
	for( int i=0 ; i<textureNodes.size() ; i++ )
	{
		if( nodeIndex( textureNodes[i].ci , textureNodes[i].cj )!=-1 ) if( false ) MK_WARN( "Multiple nodes mapped to pixel " , textureNodes[i].ci , " " , textureNodes[i].cj );
		nodeIndex( textureNodes[i].ci , textureNodes[i].cj ) = i;
	}

	BoundaryProlongationData< Real > boundaryProlongation;
	InitializeBoundaryProlongationData( hierarchy.gridAtlases[0] , boundaryProlongation );

	if( VectorField.set )
	{
		ExplicitIndexVector< AtlasMeshTriangleIndex , Point2D< PreReal > > vectorField;
		// Read in the vector field
		if( IntrinsicVectorField.set )
		{
			ReadVector( ( std::vector< Point2D< PreReal > > & )vectorField , VectorField.value );
			if( vectorField.size()!=mesh.numTriangles() ) MK_THROW( "Triangle and vector counts don't match: " , mesh.numTriangles() , " != " , vectorField.size() );
		}
		else
		{
			ExplicitIndexVector< AtlasMeshTriangleIndex , Point3D< PreReal > > _vectorField;
			ReadVector( ( std::vector< Point2D< PreReal > > & )_vectorField , VectorField.value );
			if( _vectorField.size()!=mesh.numTriangles() ) MK_THROW( "Triangle and vector counts don't match: " , mesh.numTriangles() , " != " , _vectorField.size() );
			vectorField.resize( _vectorField.size() );
			ThreadPool::ParallelFor
				(
					0 , mesh.numTriangles() ,
					[&]( unsigned int , size_t i )
					{
						Simplex< PreReal , 3 , 2 > s = mesh.surfaceTriangle( static_cast< unsigned int >(i) );
						Point3D< PreReal > d[] = { s[1]-s[0] , s[2]-s[0] };
						SquareMatrix< PreReal , 2 > Dot;
						for( int j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; k++ ) Dot(j,k) = Point3D< PreReal >::Dot( d[j] , d[k] );
						Point2D< PreReal > dot( Point3D< PreReal >::Dot( d[0] , _vectorField[ AtlasMeshTriangleIndex(i) ] ) , Point3D< PreReal >::Dot( d[1] , _vectorField[ AtlasMeshTriangleIndex(i) ] ) );
						vectorField[ AtlasMeshTriangleIndex(i) ] = Dot.inverse() * dot;
					}
				);
		}
		// Normalize the vector-field to have unit-norm
		{
			ExplicitIndexVector< AtlasMeshTriangleIndex , SquareMatrix< PreReal , 2 > > embeddingMetric;
			InitializeEmbeddingMetric( mesh , true , embeddingMetric );
			{
				PreReal norm = 0 , area = 0;
				for( int t=0 ; t<embeddingMetric.size() ; t++ )
				{
					PreReal a = (PreReal)sqrt( embeddingMetric[ AtlasMeshTriangleIndex(t) ].determinant() ) / 2.;
					norm += Point2D< PreReal >::Dot( vectorField[ AtlasMeshTriangleIndex(t) ] , embeddingMetric[ AtlasMeshTriangleIndex(t) ] * vectorField[ AtlasMeshTriangleIndex(t) ] ) * a;
					area += a;
				}
				norm = sqrt( norm / area );
				for( int t=0 ; t<embeddingMetric.size() ; t++ ) vectorField[ AtlasMeshTriangleIndex(t) ] /= norm;
			}
		}
		// Initialize the metric from the vector field
		auto LengthToAnisotropy = [&]( PreReal len )
		{
			// g <- g + gOrtho * anisotropy 
			// 0 -> 0
			// 1 -> 1e5
			// infty -> infty
			return (PreReal)( pow( len , AnisotropyExponent.value ) * AnisotropyScale.value );
		};
		InitializeAnisotropicMetric( mesh , atlasCharts , vectorField , LengthToAnisotropy , parameterMetric );
	}
	else
	{
		// Initialize the metric from the embedding
		InitializeMetric( mesh , EMBEDDING_METRIC , atlasCharts , parameterMetric );
	}

	// Scale the metric so that the area is equal to the resolution
	for( unsigned int i=0 ; i<parameterMetric.size() ; i++ ) for( unsigned int j=0 ; j<parameterMetric[ ChartIndex(i) ].size() ; j++ ) parameterMetric[ ChartIndex(i) ][ ChartMeshTriangleIndex(j) ] *= textureNodes.size() / 2;

	pMeter.reset();
	OperatorInitializer::Initialize( MatrixQuadrature.value , massAndStiffnessOperators , hierarchy.gridAtlases[0] , parameterMetric , atlasCharts , scalarIntegrator , RHSQuadrature.value , ApproximateIntegration.set );
	scalarIntegratorScratch = scalarIntegrator.template getScratch< Point< Real , 2 > , Point< Real , 2 > >();
	if( Verbose.set ) std::cout << pMeter( "Mass and stiffness" ) << std::endl;

	if( UseDirectSolver.set )
	{
		FullMatrixConstruction( hierarchy.gridAtlases[0] , massAndStiffnessOperators.massCoefficients , mass );
		FullMatrixConstruction( hierarchy.gridAtlases[0] , massAndStiffnessOperators.stiffnessCoefficients , stiffness );
		systemMatrices[0] = mass + stiffness * diffusionRates[0] * speed;
		systemMatrices[1] = mass + stiffness * diffusionRates[1] * speed;
		if( Verbose.set ) std::cout << pMeter( "Assembled" ) << std::endl;
	}

	//////////////////////////////////// Initialize multigrid indices

	multigridIndices.resize( levels );
	for( unsigned int i=0 ; i<levels ; i++ )
	{
		const typename GridAtlas<>::IndexConverter & indexConverter = hierarchy.gridAtlases[i].indexConverter;
		const GridAtlas< PreReal , Real > &gridAtlas = hierarchy.gridAtlases[i];
		multigridIndices[i].threadTasks = gridAtlas.threadTasks;
		multigridIndices[i].boundaryToCombined = indexConverter.boundaryToCombined();
		multigridIndices[i].segmentedLines = gridAtlas.segmentedLines;
		multigridIndices[i].rasterLines = gridAtlas.rasterLines;
		multigridIndices[i].restrictionLines = gridAtlas.restrictionLines;
		multigridIndices[i].prolongationLines = gridAtlas.prolongationLines;
		if( i<levels-1 ) multigridIndices[i].boundaryRestriction = hierarchy.boundaryRestriction[i];
	}

	//////////////////////////////////// Initialize multigrid coefficients

	pMeter.reset();
	for( int ab=0 ; ab<2 ; ab++ ) UpdateLinearSystem( (Real)1. , diffusionRates[ab] * speed , hierarchy , multigridCoefficients[ab] , massAndStiffnessOperators , vCycleSolvers[ab] , fineSolvers[ab] , systemMatrices[ab] , DetailVerbose.set , true , UseDirectSolver.set );
	if( Verbose.set ) std::cout << pMeter( "Initialized MG" ) << std::endl;

	//////////////////////////////////// Initialize multigrid variables

	for( unsigned int ab=0 ; ab<2 ; ab++ )
	{
		multigridVariables[ab].resize( levels );
		for( unsigned int i=0 ; i<levels ; i++ )
		{
			const typename GridAtlas<>::IndexConverter & indexConverter = hierarchy.gridAtlases[i].indexConverter;
			MultigridLevelVariables< Real >& variables = multigridVariables[ab][i];
			variables.x.resize( indexConverter.numCombined() );
			variables.rhs.resize( indexConverter.numCombined() );
			variables.residual.resize( indexConverter.numCombined() );
			variables.boundary_rhs.resize( indexConverter.numBoundary() );
			variables.boundary_value.resize( indexConverter.numBoundary() );
			variables.variable_boundary_value.resize( indexConverter.numBoundary() );
		}
	}
}

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::InitializeVisualization( void )
{
	unsigned int tCount = (unsigned int)mesh.numTriangles();

	visualization.triangles.resize( tCount );
	visualization.vertices.resize( 3*tCount );
	visualization.colors.resize( 3*tCount , Point3D< float >( 0.75f , 0.75f , 0.75f ) );
	visualization.textureCoordinates.resize( 3*tCount );
	visualization.normals.resize( 3*tCount );

	for( unsigned int t=0 , idx=0 ; t<tCount ; t++ )
	{
		Simplex< PreReal , 3 , 2 > sSimplex = mesh.surfaceTriangle(t);
		Simplex< PreReal , 2 , 2 > tSimplex = mesh.textureTriangle(t);
		Point3D< float > n = sSimplex.normal();
		n /= Point3D< float >::Length( n );

		for( int k=0 ; k<3 ; k++ , idx++ )
		{
			visualization.triangles[t][k] = idx;
			visualization.vertices[idx] = sSimplex[k];
			visualization.normals[idx] = n;
			visualization.textureCoordinates[idx] = tSimplex[k];
		}
	}

	std::vector< unsigned int > boundaryHalfEdges = mesh.texture.boundaryHalfEdges();

	for( int e=0 ; e<boundaryHalfEdges.size() ; e++ )
	{
		SimplexIndex< 1 > eIndex = mesh.surface.edgeIndex( boundaryHalfEdges[e] );
		for( int i=0 ; i<2 ; i++ ) visualization.chartBoundaryVertices.push_back( mesh.surface.vertices[ eIndex[i] ] );
	}

	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization,  's' , "export texture" , "Output Texture" , ExportTextureCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , ' ' , "toggle update" , ToggleUpdateCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , '+' , "increment update" , IncrementUpdateCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , '2' , "show concentration 2" , SetConcentration2CallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , '1' , "show concentration 1" , SetConcentration1CallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , '0' , "hide concetrations" , HideConcentrationsCallBack ) );

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

	visualization.info.push_back( stepsString );
}

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::SetConcentration1CallBack( Visualization* , const char* ){ whichConcentration = 0; }

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::SetConcentration2CallBack( Visualization* , const char* ){ whichConcentration = 1; }

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::HideConcentrationsCallBack( Visualization* , const char* ){ whichConcentration = -1; }

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::ToggleUpdateCallBack( Visualization* , const char* )
{
	if( updateCount ) updateCount = 0;
	else              updateCount = -1;
}

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::IncrementUpdateCallBack( Visualization* , const char* )
{
	if( updateCount<0 ) updateCount = 1;
	else updateCount++;
}

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::Init( void )
{
	sprintf( stepsString , "Steps: 0" );
	levels = std::max< unsigned int >( Levels.value , 1 );
	feed = FeedKillRates.values[0];
	kill = FeedKillRates.values[1];
	speed = Speed.value;
	diffusionRates[0] = 1.0;
	diffusionRates[1] = 0.5;
	for( int ab=0 ; ab<2 ; ab++ ) diffusionRates[ab] *= DiffusionScale.value / 10.;
	textureWidth = Width.value;
	textureHeight = Height.value;
	mesh.read( Input.value , DetailVerbose.set , CollapseEpsilon.value );
	
	if( RandomJitter.set )
	{
		if( RandomJitter.value ) srand( RandomJitter.value );
		else                     srand( time(NULL) );
		PreReal jitterScale = (PreReal)1e-3 / std::max< int >( textureWidth , textureHeight );
		for( int i=0 ; i<mesh.texture.vertices.size() ; i++ ) mesh.texture.vertices[i] += Point2D< PreReal >( (PreReal)1. - Random< PreReal >()*2 , (PreReal)1. - Random< PreReal >()*2 ) * jitterScale;
	}

	{
		padding = Padding::Init( textureWidth , textureHeight , mesh.texture.vertices , DetailVerbose.set );
		padding.pad( textureWidth , textureHeight , mesh.texture.vertices );
		textureWidth  += padding.width();
		textureHeight += padding.height();
	}

	// Define centroid and scale for visualization
	Point3D< PreReal > centroid = mesh.surface.centroid();
	PreReal radius = mesh.surface.boundingRadius( centroid );
	for( unsigned int i=0 ; i<mesh.surface.vertices.size() ; i++ ) mesh.surface.vertices[i] = ( mesh.surface.vertices[i]-centroid ) / radius;

	Miscellany::PerformanceMeter pMeter( '.' );
	InitializeSystem( textureWidth , textureHeight );
	if( Verbose.set ) std::cout << pMeter( "Initialized" ) << std::endl;
	if( Verbose.set ) printf( "Resolution: %d / %d x %d\n" , (int)textureNodes.size() , textureWidth , textureHeight );

	// Assign position to exterior nodes using barycentric-exponential map
	{
		FEM::RiemannianMesh< PreReal , unsigned int > rMesh( GetPointer( mesh.surface.triangles ) , mesh.surface.triangles.size() );
		rMesh.setMetricFromEmbedding( GetPointer( mesh.surface.vertices ) );
		rMesh.makeUnitArea();
		Pointer( FEM::CoordinateXForm< PreReal > ) xForms = rMesh.getCoordinateXForms();

		for( unsigned int i=0 ; i<textureNodes.size() ; i++ ) if( textureNodes[i].tID!=AtlasMeshTriangleIndex(-1) && !textureNodes[i].isInterior )
		{
			FEM::HermiteSamplePoint< PreReal > _p;
			_p.tIdx = static_cast< unsigned int >(textureNodes[i].tID);
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
		const int SAMPLES = (int)( multigridVariables[1][0].x.size() * SamplesFraction.value );
		randomSamples = rMesh.randomSamples( SAMPLES );
	}

	textureNodePositions.resize( textureNodes.size() );
	for( int i=0 ; i<textureNodePositions.size() ; i++ ) textureNodePositions[i] = mesh.surface( textureNodes[i] );

	InitializeConcentrations();

	outputBuffer = new unsigned char[ textureHeight*textureWidth];
	memset( outputBuffer , 204 , textureHeight * textureWidth * sizeof(unsigned char) );
}

template< typename PreReal , typename Real >
void _main( int argc , char* argv[] )
{
	GrayScottReactionDiffusion< PreReal , Real >::Init();
	if( Run.set ) GrayScottReactionDiffusion< PreReal , Real >::updateCount = -1;
	if( !Output.set )
	{
		glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
		GrayScottReactionDiffusion< PreReal , Real >::visualization.displayMode = DisplayMode.value;
		if     ( DisplayMode.value==ONE_REGION_DISPLAY ) GrayScottReactionDiffusion< PreReal , Real >::visualization.screenWidth =  800 , GrayScottReactionDiffusion< PreReal , Real >::visualization.screenHeight = 800;
		else if( DisplayMode.value==TWO_REGION_DISPLAY ) GrayScottReactionDiffusion< PreReal , Real >::visualization.screenWidth = 1440 , GrayScottReactionDiffusion< PreReal , Real >::visualization.screenHeight = 720;
		GrayScottReactionDiffusion< PreReal , Real >::visualization.useNearestSampling = NearestSampling.set;
		glutInitWindowSize( GrayScottReactionDiffusion< PreReal , Real >::visualization.screenWidth , GrayScottReactionDiffusion< PreReal , Real >::visualization.screenHeight );

		glutInit( &argc , argv );
		char windowName[1024];
		sprintf( windowName , "Gray-Scott Reaction Diffusion");
		glutCreateWindow( windowName );
		if( glewInit()!=GLEW_OK ) MK_THROW( "glewInit failed" );
		glutDisplayFunc ( GrayScottReactionDiffusion< PreReal , Real >::Display );
		glutReshapeFunc ( GrayScottReactionDiffusion< PreReal , Real >::Reshape );
		glutMouseFunc   ( GrayScottReactionDiffusion< PreReal , Real >::MouseFunc );
		glutMotionFunc  ( GrayScottReactionDiffusion< PreReal , Real >::MotionFunc );
		glutKeyboardFunc( GrayScottReactionDiffusion< PreReal , Real >::KeyboardFunc );
		glutIdleFunc    ( GrayScottReactionDiffusion< PreReal , Real >::Idle );
		if( CameraConfig.set ) GrayScottReactionDiffusion< PreReal , Real >::visualization.ReadSceneConfigurationCallBack( &GrayScottReactionDiffusion< PreReal , Real >::visualization , CameraConfig.value.c_str() );
		GrayScottReactionDiffusion< PreReal , Real >::InitializeVisualization();
		glutMainLoop();
	}
	else
	{
		Miscellany::Timer timer;
		for( unsigned int i=0 ; i<OutputSteps.value ; i++ )
		{
			if( Verbose.set ) printf( "%d / %d \r" , i+1 , OutputSteps.value );
			if( UseDirectSolver.set ) GrayScottReactionDiffusion< PreReal , Real >::UpdateExactSolution();
			else                      GrayScottReactionDiffusion< PreReal , Real >::UpdateApproximateSolution();
		}
		if( Verbose.set )
		{
			printf( "\n" );
			double total_time = timer.elapsed();
			printf( "Reaction-diffusion total time / time per iteration: %.2f(s) / %.4f(s)\n" , total_time , total_time / OutputSteps.value );
		}
		GrayScottReactionDiffusion< PreReal , Real >::SetOutputBuffer( GrayScottReactionDiffusion< PreReal , Real >::multigridVariables[1][0].x );
		GrayScottReactionDiffusion< PreReal , Real >::ExportTextureCallBack( &GrayScottReactionDiffusion< PreReal , Real >::visualization , Output.value.c_str() );
	}
}

int main( int argc , char* argv[] )
{
	CmdLineParse( argc-1 , argv+1 , params );
	if( !Input.set ){ ShowUsage( argv[0] ) ; return EXIT_FAILURE; }
	if( !SamplesFraction.set )
	{
		if( Dots.set ) SamplesFraction.value = DefaultDotsSamplesFraction;
		else           SamplesFraction.value = DefaultStripesSamplesFraction;
	}
	if( Dots.set && !FeedKillRates.set ) FeedKillRates.values[0] = DotRates[0] , FeedKillRates.values[1] = DotRates[1];
	if( Serial.set ) ThreadPool::ParallelizationType = ThreadPool::ParallelType::NONE;
	if( !NoHelp.set && !Output.set )
	{
		printf( "+------------------------------------------------+\n" );
		printf( "| Interface Controls:                            |\n" );
		printf( "|    [Left Mouse]:                 rotate        |\n" );
		printf( "|    [Right Mouse]:                zoom          |\n" );
		printf( "|    [Left/Right Mouse] + [CTRL]:  pan           |\n" );
		printf( "|    [Left/Right Mouse] + [SHIFT]: set seed      |\n" );
		printf( "|    [SPACE]:                      start process |\n" );
		printf( "+------------------------------------------------+\n" );
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