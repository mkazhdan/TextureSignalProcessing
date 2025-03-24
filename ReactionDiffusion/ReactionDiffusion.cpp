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
#include <Src/Hierarchy.h>
#include <Src/SimpleTriangleMesh.h>
#include <Src/Basis.h>
#include <Src/Solver.h>
#include <Src/QuadratureIntergration.h>
#include <Src/MassAndStiffness.h>
#include <Src/Padding.h>
#include <Src/TexturedMeshVisualization.h>

using namespace MishaK;

const float StripeRates[] = { 0.062f , 0.062f };
const float    DotRates[] = { 0.0367f , 0.0649f };
const float DefaultStripesSamplesFraction = 0.01f;
const float DefaultDotsSamplesFraction = 0.001f;

CmdLineParameter< std::string > Input( "in" );
CmdLineParameter< std::string > Output( "out" );
CmdLineParameter< int   > OutputSteps( "outSteps" , 1000 );
CmdLineParameter< int   > Width( "width" , 512 );
CmdLineParameter< int   > Height( "height" , 512 );
CmdLineParameter< float > Speed( "speed" , 10.f );
CmdLineParameterArray< float , 2 > FeedKillRates( "fk" , StripeRates );
CmdLineParameter< float > DiffusionScale( "diff" , 1.f );
CmdLineParameter< float > SamplesFraction( "samples" );
CmdLineParameter< int   > Levels( "levels" , 4 );
CmdLineParameter< std::string > CameraConfig( "camera" );
CmdLineParameter< int   > DisplayMode( "display" , TWO_REGION_DISPLAY );
CmdLineParameter< int   > MatrixQuadrature( "mQuadrature" , 6 );
CmdLineParameter< int   > RHSQuadrature( "rhsQuadrature" , 3 );

CmdLineParameter< int   > MultigridBlockHeight ( "mBlockH" ,  16 );
CmdLineParameter< int   > MultigridBlockWidth  ( "mBlockW" , 128 );
CmdLineParameter< int   > MultigridPaddedHeight( "mPadH"   ,   0 );
CmdLineParameter< int   > MultigridPaddedWidth ( "mPadW"   ,   2 );

CmdLineParameter< int   > RandomJitter( "jitter" , 0 );
CmdLineReadable Verbose( "verbose" );
CmdLineReadable NoHelp( "noHelp" );
CmdLineReadable DetailVerbose( "detail" );
CmdLineReadable UseDirectSolver( "useDirectSolver" );
CmdLineReadable Double( "double" );
CmdLineReadable ApproximateIntegration( "approximateIntegration" );
CmdLineReadable Dots( "dots" );
CmdLineReadable Serial( "serial" );

CmdLineParameter< std::string > VectorField( "inVF" );
CmdLineParameter< float > AnisotropyScale( "aScl" , 1.f );
CmdLineParameter< float > AnisotropyExponent( "aExp" , 0.f );
CmdLineReadable IntrinsicVectorField( "intrinsicVF" );

#ifdef NEW_CODE
CmdLineParameter< double > CollapseEpsilon( "collapse" , 0 );
#endif // NEW_CODE

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
#ifdef NEW_CODE
	&CollapseEpsilon ,
#endif // NEW_CODE
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
#ifdef NEW_CODE
	printf( "\t[--%s <collapse epsilon>=%g]\n" , CollapseEpsilon.name.c_str() , CollapseEpsilon.value );
#endif // NEW_CODE
	printf( "\t[--%s]\n" , IntrinsicVectorField.name.c_str() );
	printf( "\t[--%s]\n" , Serial.name.c_str() );

	printf( "\t[--%s]\n" , NoHelp.name.c_str() );
}

template< typename PreReal , typename Real >
class GrayScottReactionDiffusion
{
public:
	static std::vector< FEM::SamplePoint< PreReal > > randomSamples;
	static OrientedTexturedTriangleMesh< PreReal > mesh;
	static int textureWidth;
	static int textureHeight;
	static Real diffusionRates[2];
	static Real kill;
	static Real feed;
	static Real speed;
	static int levels;
	static int steps;
	static char stepsString[];
	static int whichConcentration;

	static Padding padding;

	static std::vector< Point3D< float > > textureNodePositions;

	static HierarchicalSystem< PreReal , Real > hierarchy;

	static std::vector< BilinearElementIndex > bilinearElementIndices;

	static std::vector< TextureNodeInfo< PreReal > > textureNodes;
	static Image< int > nodeIndex;

	static SparseMatrix< Real , int > mass;
	static SparseMatrix< Real , int > stiffness;
	static SparseMatrix< Real , int > systemMatrices[2];

	static int seedTexel;

	static std::vector< AtlasChart< PreReal > > atlasCharts;
	static std::vector< std::vector< SquareMatrix< PreReal , 2 > > > parameterMetric;

	static std::vector< SystemCoefficients< Real > > multigridCoefficients[2];
	static std::vector< MultigridLevelVariables< Real > > multigridVariables[2];

#if defined( USE_CHOLMOD )
	typedef CholmodCholeskySolver< Real , 1 > DirectSolver;
#elif defined( USE_EIGEN )
	typedef EigenCholeskySolver< Real , 1 > DirectSolver;
#elif defined( USE_EIGEN_PARDISO )
	typedef EigenPardisoSolver< Real , 1 > DirectSolver;
#else
#error "[ERROR] No solver defined!"
#endif

	static VCycleSolvers< DirectSolver > vCycleSolvers[2];
	static DirectSolver fineSolvers[2];

	static std::vector< MultigridLevelIndices< Real > > multigridIndices;

	static SparseMatrix< Real , int > coarseBoundaryFineBoundaryProlongation;
	static SparseMatrix< Real , int > fineBoundaryCoarseBoundaryRestriction;
	static std::vector< Point2D< Real > > coarseBoundaryValues;
	static std::vector< Point2D< Real > > coarseBoundaryRHS;
	static std::vector< Point2D< Real > > fineBoundaryValues;
	static std::vector< Point2D< Real > > fineBoundaryRHS;

	//Samples
	static ScalarElementSamples< Real > scalarSamples;
	static std::vector< InteriorCellLine > interiorCellLines;
	static std::vector< std::pair< int , int > > interiorCellLineIndex;

	//Linear Operators
	static SystemCoefficients< Real > massCoefficients;
	static SystemCoefficients< Real > stiffnessCoefficients;

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
template< typename PreReal , typename Real > OrientedTexturedTriangleMesh< PreReal >								GrayScottReactionDiffusion< PreReal , Real >::mesh;
template< typename PreReal , typename Real > int																	GrayScottReactionDiffusion< PreReal , Real >::textureWidth;
template< typename PreReal , typename Real > int																	GrayScottReactionDiffusion< PreReal , Real >::textureHeight;
template< typename PreReal , typename Real > TexturedMeshVisualization												GrayScottReactionDiffusion< PreReal , Real >::visualization;
template< typename PreReal , typename Real > int																	GrayScottReactionDiffusion< PreReal , Real >::mouseX = -1;
template< typename PreReal , typename Real > int																	GrayScottReactionDiffusion< PreReal , Real >::mouseY = -1;
template< typename PreReal , typename Real > bool																	GrayScottReactionDiffusion< PreReal , Real >::mouseSelectionActive = false;
template< typename PreReal , typename Real > Padding																GrayScottReactionDiffusion< PreReal , Real >::padding;

template< typename PreReal , typename Real > std::vector< AtlasChart< PreReal > >									GrayScottReactionDiffusion< PreReal , Real> ::atlasCharts;
template< typename PreReal , typename Real > std::vector< std::vector< SquareMatrix< PreReal , 2 > > >				GrayScottReactionDiffusion< PreReal , Real >::parameterMetric;

template< typename PreReal , typename Real > SparseMatrix< Real , int >												GrayScottReactionDiffusion< PreReal , Real >::mass;
template< typename PreReal , typename Real > SparseMatrix< Real , int >												GrayScottReactionDiffusion< PreReal , Real >::stiffness;
template< typename PreReal , typename Real > SparseMatrix< Real , int >												GrayScottReactionDiffusion< PreReal , Real >::systemMatrices[2];

template< typename PreReal , typename Real > Real																	GrayScottReactionDiffusion< PreReal , Real >::diffusionRates[2];
template< typename PreReal , typename Real > Real																	GrayScottReactionDiffusion< PreReal , Real >::speed;
template< typename PreReal , typename Real > Real																	GrayScottReactionDiffusion< PreReal , Real >::kill;
template< typename PreReal , typename Real > Real																	GrayScottReactionDiffusion< PreReal , Real >::feed;

template< typename PreReal , typename Real > std::vector< TextureNodeInfo< PreReal > >								GrayScottReactionDiffusion< PreReal , Real >::textureNodes;
template< typename PreReal , typename Real > Image< int >															GrayScottReactionDiffusion< PreReal , Real >::nodeIndex;
template< typename PreReal , typename Real > std::vector< BilinearElementIndex >									GrayScottReactionDiffusion< PreReal , Real >::bilinearElementIndices;

template< typename PreReal , typename Real > int																	GrayScottReactionDiffusion< PreReal , Real >::steps;
template< typename PreReal , typename Real > char																	GrayScottReactionDiffusion< PreReal , Real >::stepsString[1024];
template< typename PreReal , typename Real > int																	GrayScottReactionDiffusion< PreReal , Real >::levels;
template< typename PreReal , typename Real > HierarchicalSystem< PreReal , Real >									GrayScottReactionDiffusion< PreReal , Real >::hierarchy;

template< typename PreReal , typename Real > unsigned char *														GrayScottReactionDiffusion< PreReal , Real >::outputBuffer;
template< typename PreReal , typename Real > std::vector< MultigridLevelIndices< Real > >							GrayScottReactionDiffusion< PreReal , Real >::multigridIndices;

template< typename PreReal , typename Real > std::vector< SystemCoefficients< Real > >								GrayScottReactionDiffusion< PreReal , Real >::multigridCoefficients[2];
template< typename PreReal , typename Real > std::vector< MultigridLevelVariables< Real > >							GrayScottReactionDiffusion< PreReal , Real >::multigridVariables[2];
template< typename PreReal , typename Real > VCycleSolvers< typename GrayScottReactionDiffusion< PreReal , Real >::DirectSolver >	GrayScottReactionDiffusion< PreReal , Real >::vCycleSolvers[2];
template< typename PreReal , typename Real > typename GrayScottReactionDiffusion< PreReal , Real >::DirectSolver	GrayScottReactionDiffusion< PreReal , Real >::fineSolvers[2];



//Samples
template< typename PreReal , typename Real > ScalarElementSamples< Real >											GrayScottReactionDiffusion< PreReal , Real >::scalarSamples;
template< typename PreReal , typename Real > std::vector< InteriorCellLine >										GrayScottReactionDiffusion< PreReal , Real >::interiorCellLines;
template< typename PreReal , typename Real > std::vector< std::pair< int , int > >									GrayScottReactionDiffusion< PreReal , Real >::interiorCellLineIndex;

template< typename PreReal , typename Real > int																	GrayScottReactionDiffusion< PreReal , Real >::seedTexel = -1;
template< typename PreReal , typename Real > std::vector< Point3D< float > >										GrayScottReactionDiffusion< PreReal , Real >::textureNodePositions;

template< typename PreReal , typename Real > SparseMatrix< Real , int >												GrayScottReactionDiffusion< PreReal , Real >::coarseBoundaryFineBoundaryProlongation;
template< typename PreReal , typename Real > SparseMatrix< Real , int >												GrayScottReactionDiffusion< PreReal , Real >::fineBoundaryCoarseBoundaryRestriction;

template< typename PreReal , typename Real > std::vector< Point2D< Real > >											GrayScottReactionDiffusion< PreReal , Real >::coarseBoundaryValues;
template< typename PreReal , typename Real > std::vector< Point2D< Real > >											GrayScottReactionDiffusion< PreReal , Real >::coarseBoundaryRHS;
template< typename PreReal , typename Real > std::vector< Point2D< Real > >											GrayScottReactionDiffusion< PreReal , Real >::fineBoundaryValues;
template< typename PreReal , typename Real > std::vector< Point2D< Real > >											GrayScottReactionDiffusion< PreReal , Real >::fineBoundaryRHS;

template< typename PreReal , typename Real > SystemCoefficients< Real >												GrayScottReactionDiffusion< PreReal , Real >::massCoefficients;
template< typename PreReal , typename Real > SystemCoefficients< Real >												GrayScottReactionDiffusion< PreReal , Real >::stiffnessCoefficients;

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
	for( int ab=0 ; ab<2 ; ab++ ) ThreadPool::ParallelFor( 0 , ab_x.size() , [&]( unsigned int , size_t i ){ ab_x[i][ab] = multigridVariables[ab][0].x[i]; } );
	const std::vector< int >& boundaryGlobalIndex = hierarchy.gridAtlases[0].boundaryGlobalIndex;

	ThreadPool::ParallelFor( 0 , boundaryGlobalIndex.size() , [&]( unsigned int , size_t i ){coarseBoundaryValues[i] = ab_x[ boundaryGlobalIndex[i] ]; } );
	coarseBoundaryFineBoundaryProlongation.Multiply( &coarseBoundaryValues[0] , &fineBoundaryValues[0] );

	for( int ab=0 ; ab<2 ; ab++ ) MultiplyBySystemMatrix_NoReciprocals( massCoefficients , hierarchy.gridAtlases[0].boundaryGlobalIndex , hierarchy.gridAtlases[0].rasterLines , multigridVariables[ab][0].x , multigridVariables[ab][0].rhs );
	auto ABFunction = [&]( Point2D< Real > ab , SquareMatrix< Real , 2 > )
	{
		return Point2D< Real >
			(
				(Real)( speed * ( - ab[0] * ab[1] * ab[1] + feed * ( 1 - ab[0] ) ) ) ,
				(Real)( speed * (   ab[0] * ab[1] * ab[1] - ( kill + feed ) * ab[1] ) )
			);
	};
	memset( &ab_rhs[0] , 0 , ab_rhs.size() * sizeof( Point2D< Real > ) );
	memset( &fineBoundaryRHS[0] , 0 , fineBoundaryRHS.size() * sizeof( Point2D< Real > ) );
	Integrate< Real >( interiorCellLines , scalarSamples , ab_x , fineBoundaryValues , ABFunction , ab_rhs , fineBoundaryRHS );
	fineBoundaryCoarseBoundaryRestriction.Multiply( &fineBoundaryRHS[0] , &coarseBoundaryRHS[0] );
	for( int ab=0 ; ab<2 ; ab++ )
	{
		ThreadPool::ParallelFor( 0 , ab_rhs.size() , [&]( unsigned int , size_t i ){ multigridVariables[ab][0].rhs[i] += ab_rhs[i][ab]; } );
		ThreadPool::ParallelFor( 0 , boundaryGlobalIndex.size() , [&]( unsigned int , size_t i ){ multigridVariables[ab][0].rhs[ boundaryGlobalIndex[i] ] += coarseBoundaryRHS[i][ab]; } );
	}
}

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::UpdateExactSolution( bool verbose )
{
	// (1) Compute the right-hand-sides
	{
		Miscellany::Timer timer;
		SetRightHandSide();
		if( verbose ) printf( "Integrating normalized vector field %.4f\n" , timer.elapsed() );
	}

	// (2) Solve the linear systems
	{
		Miscellany::Timer timer;
		for( int ab=0 ; ab<2 ; ab++ )
		{
			solve( fineSolvers[ab] , multigridVariables[ab][0].x , multigridVariables[ab][0].rhs );
			for( int i=0 ; i<multigridVariables[ab][0].x.size() ; i++ ) multigridVariables[ab][0].x[i] = std::max< Real >( multigridVariables[ab][0].x[i] , 0 );
		}
		if( verbose ) printf( "Performed direct solve %.4f\n" , timer.elapsed() );
	}
}

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::UpdateApproximateSolution( bool verbose , bool detailVerbose )
{
	double rhsTime , vCycleTime;
	// Compute the right-hand-sides
	{
		Miscellany::Timer timer;
		SetRightHandSide();
		rhsTime = timer.elapsed();
	}

	// Solve the linear systems
	{
		Miscellany::Timer timer;
		for( int ab=0 ; ab<2 ; ab++ )
		{
			VCycle( multigridVariables[ab] , multigridCoefficients[ab] , multigridIndices , vCycleSolvers[ab] , detailVerbose , detailVerbose );
			ThreadPool::ParallelFor( 0 , multigridVariables[ab][0].x.size() , [&]( unsigned int , size_t i ){ multigridVariables[ab][0].x[i] = std::max< Real >( multigridVariables[ab][0].x[i] , 0 ); } );
		}
		vCycleTime = timer.elapsed();
	}
	if( verbose ) printf( "Integrated constraints / performed v-cycle: %.3f / %.3f(s)\n" , rhsTime , vCycleTime );
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
		else                      UpdateApproximateSolution();

		if( updateCount>0 ) updateCount--;
		steps++;
		sprintf( stepsString , "Steps: %d" , steps );
	}
	if( whichConcentration<0 ) UpdateOutputBuffer();
	else                       UpdateOutputBuffer( multigridVariables[whichConcentration][0].x );
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
	Image< Point3D< Real > > outputImage;
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
			Point2D< Real > t =
				mesh.textureCoordinates[ tIdx*3 + 0 ] * (Real)( 1. - p[0] - p[1] ) +
				mesh.textureCoordinates[ tIdx*3 + 1 ] * (Real)(      p[0]        ) +
				mesh.textureCoordinates[ tIdx*3 + 2 ] * (Real)(             p[1] );
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
	Miscellany::Timer timer;
	MultigridBlockInfo multigridBlockInfo( MultigridBlockWidth.value , MultigridBlockHeight.value , MultigridPaddedWidth.value , MultigridPaddedHeight.value , 0 );

	InitializeHierarchy( mesh , width , height , levels , textureNodes , bilinearElementIndices , hierarchy , atlasCharts , multigridBlockInfo , true , DetailVerbose.set );
	if( Verbose.set ) printf( "\tInitialized hierarchy: %.2f(s)\n" , timer.elapsed() );

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

	std::vector< Point3D< Real > > __inputSignal;
	std::vector< Real > __texelToCellCoeffs;
	SparseMatrix< Real , int > __boundaryCellBasedStiffnessRHSMatrix[3];

	if( VectorField.set )
	{
		std::vector< Point2D< PreReal > > vectorField;
		// Read in the vector field
		if( IntrinsicVectorField.set )
		{
			ReadVector( vectorField , VectorField.value );
			if( vectorField.size()!=mesh.triangles.size() ) MK_THROW( "Triangle and vector counts don't match: " , mesh.triangles.size() , " != " , vectorField.size() );
		}
		else
		{
			std::vector< Point3D< PreReal > > _vectorField;
			ReadVector( _vectorField , VectorField.value );
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
				for( int t=0 ; t<embeddingMetric.size() ; t++ ) vectorField[t] /= norm;
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
	for( int i=0 ; i<parameterMetric.size() ; i++ ) for( int j=0 ; j<parameterMetric[i].size() ; j++ ) parameterMetric[i][j] *= textureNodes.size() / 2;

	timer.reset();
	{
		switch( MatrixQuadrature.value )
		{
		case 1:  InitializeMassAndStiffness< 1>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
		case 3:  InitializeMassAndStiffness< 3>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
		case 6:  InitializeMassAndStiffness< 6>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
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
		FullMatrixConstruction( hierarchy.gridAtlases[0] , massCoefficients , mass );
		FullMatrixConstruction( hierarchy.gridAtlases[0] , stiffnessCoefficients , stiffness );
		systemMatrices[0] = mass + stiffness * diffusionRates[0] * speed;
		systemMatrices[1] = mass + stiffness * diffusionRates[1] * speed;
		if( Verbose.set ) printf( "\tAssembled matrices: %.2f(s)\n" , tmr.elapsed() );
	}

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

	timer.reset();
	for( int ab=0 ; ab<2 ; ab++ ) UpdateLinearSystem( (Real)1. , diffusionRates[ab] * speed , hierarchy , multigridCoefficients[ab] , massCoefficients , stiffnessCoefficients , vCycleSolvers[ab] , fineSolvers[ab] , systemMatrices[ab] , DetailVerbose.set , true , UseDirectSolver.set );
	if( Verbose.set ) printf( "\tInitialized multigrid coefficients: %.2f(s)\n" , timer.elapsed() );

	//////////////////////////////////// Initialize multigrid variables

	for( int ab=0 ; ab<2 ; ab++ )
	{
		multigridVariables[ab].resize( levels );
		for( int i=0 ; i<levels ; i++ )
		{
			MultigridLevelVariables< Real >& variables = multigridVariables[ab][i];
		 	variables.x.resize( hierarchy.gridAtlases[i].numTexels );
			variables.rhs.resize( hierarchy.gridAtlases[i].numTexels );
		 	variables.residual.resize( hierarchy.gridAtlases[i].numTexels );
			variables.boundary_rhs.resize( hierarchy.gridAtlases[i].boundaryGlobalIndex.size() );
			variables.boundary_value.resize( hierarchy.gridAtlases[i].boundaryGlobalIndex.size() );
			variables.variable_boundary_value.resize( hierarchy.gridAtlases[i].boundaryGlobalIndex.size() );
		}
	}


	//////////////////////////////////// Initialize cell samples

	InitializeGridAtlasInteriorCellLines( hierarchy.gridAtlases[0].gridCharts , interiorCellLines , interiorCellLineIndex );
	if( interiorCellLineIndex.size()!=hierarchy.gridAtlases[0].numInteriorCells ) MK_THROW( "Inconsistent number of interior cells! Expected " , hierarchy.gridAtlases[0].numInteriorCells , " . Result " , interiorCellLineIndex.size() , "." );

	coarseBoundaryFineBoundaryProlongation = boundaryProlongation.coarseBoundaryFineBoundaryProlongation;
	fineBoundaryCoarseBoundaryRestriction = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction;
	std::vector< int > fineBoundaryIndex = boundaryProlongation.fineBoundaryIndex;
	int numFineBoundarNodes = boundaryProlongation.numFineBoundarNodes;

	scalarSamples.resize( interiorCellLines.size() );

	timer.reset();
	{
		switch( RHSQuadrature.value )
		{
		case  1: InitializeIntegration<  1 >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , scalarSamples , ApproximateIntegration.set ) ; break;
		case  3: InitializeIntegration<  3 >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , scalarSamples , ApproximateIntegration.set ) ; break;
		case  6: InitializeIntegration<  6 >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , scalarSamples , ApproximateIntegration.set ) ; break;
		case 12: InitializeIntegration< 12 >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , scalarSamples , ApproximateIntegration.set ) ; break;
		case 24: InitializeIntegration< 24 >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , scalarSamples , ApproximateIntegration.set ) ; break;
		case 32: InitializeIntegration< 32 >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , scalarSamples , ApproximateIntegration.set ) ; break;
		default: MK_THROW( "Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles" );
		}
	}
	if( Verbose.set ) printf( "\tInitialized vector field integration: %.2f(s)\n" , timer.elapsed() );
	coarseBoundaryValues.resize( hierarchy.gridAtlases[0].numTexels - hierarchy.gridAtlases[0].numDeepTexels );
	coarseBoundaryRHS.resize   ( hierarchy.gridAtlases[0].numTexels - hierarchy.gridAtlases[0].numDeepTexels );
	fineBoundaryValues.resize( numFineBoundarNodes );
	fineBoundaryRHS.resize   ( numFineBoundarNodes );

	scalarSamples.sort();
}

template< typename PreReal , typename Real >
void GrayScottReactionDiffusion< PreReal , Real >::InitializeVisualization( void )
{
	int tCount = (int)mesh.triangles.size();

	visualization.triangles.resize( tCount );
	visualization.vertices.resize( 3*tCount );
	visualization.colors.resize( 3*tCount , Point3D< float >( 0.75f , 0.75f , 0.75f ) );
	visualization.textureCoordinates.resize( 3*tCount );
	visualization.normals.resize( 3*tCount );


	for( int i=0 ; i<tCount ; i++ ) for( int k=0 ; k<3 ; k++ ) visualization.triangles[i][k] = 3*i+k;

	for( int i=0 ; i<tCount ; i++ ) for ( int j=0 ; j<3 ; j++ )
	{
		visualization.vertices          [3*i+j] = mesh.vertices[ mesh.triangles[i][j] ];
		visualization.normals           [3*i+j] = mesh.normals [ mesh.triangles[i][j] ];
		visualization.textureCoordinates[3*i+j] = mesh.textureCoordinates[3*i+j];
	}

	std::vector< int > boundaryEdges;
	mesh.initializeBoundaryEdges( boundaryEdges );

	for( int e=0 ; e<boundaryEdges.size() ; e++ )
	{
		int tIndex = boundaryEdges[e] / 3;
		int kIndex = boundaryEdges[e] % 3;
		for( int c=0 ; c<2 ; c++ )
		{
			Point3D< float > v = Point3D< float >( mesh.vertices[ mesh.triangles[tIndex][ (kIndex+c)%3 ] ] );
			visualization.boundaryEdgeVertices.push_back(v);
		}
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
	levels = Levels.value;
	feed = FeedKillRates.values[0];
	kill = FeedKillRates.values[1];
	speed = Speed.value;
	diffusionRates[0] = 1.0;
	diffusionRates[1] = 0.5;
	for( int ab=0 ; ab<2 ; ab++ ) diffusionRates[ab] *= DiffusionScale.value / 10.;
	textureWidth = Width.value;
	textureHeight = Height.value;
#ifdef NEW_CODE
	mesh.read( Input.value , DetailVerbose.set , CollapseEpsilon.value );
#else // !NEW_CODE
	mesh.read( Input.value , DetailVerbose.set );
#endif // NEW_CODE
	
	if( RandomJitter.set )
	{
		if( RandomJitter.value ) srand( RandomJitter.value );
		else                     srand( time(NULL) );
		std::vector< Point2D< PreReal > > randomOffset( mesh.vertices.size() );
		PreReal jitterScale = (PreReal)1e-3 / std::max< int >( textureWidth , textureHeight );
		for( int i=0 ; i<randomOffset.size() ; i++ ) randomOffset[i] = Point2D< PreReal >( (PreReal)1. - Random< PreReal >()*2 , (PreReal)1 - Random< PreReal >()*2 )*jitterScale;
		for( int i=0 ; i<mesh.triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ ) mesh.textureCoordinates[3*i+k] += randomOffset[ mesh.triangles[i][k] ];
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
	InitializeSystem( textureWidth , textureHeight );

	if( Verbose.set )
	{
		printf( "Resolution: %d / %d x %d\n" , (int)textureNodes.size() , textureWidth , textureHeight );
		printf( "Initialized system %.2f(s)\n" , timer.elapsed() );
		printf( "Peak Memory (MB): %d\n" , Miscellany::MemoryInfo::PeakMemoryUsageMB() );
	}

	//Assign position to exterior nodes using barycentric-exponential map
	{
		FEM::RiemannianMesh< PreReal , unsigned int > rMesh( GetPointer( mesh.triangles ) , mesh.triangles.size() );

		rMesh.setMetricFromEmbedding( GetPointer( mesh.vertices ) );
		rMesh.makeUnitArea();
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
		const int SAMPLES = (int)( multigridVariables[1][0].x.size() * SamplesFraction.value );
		randomSamples = rMesh.randomSamples( SAMPLES );
	}

	textureNodePositions.resize( textureNodes.size() );
	for( int i=0 ; i<textureNodePositions.size() ; i++ )
	{
		Point2D< PreReal > barycentricCoords = textureNodes[i].barycentricCoords;
		int tID = textureNodes[i].tID;
		Point3D< PreReal > p =
			mesh.vertices[ mesh.triangles[tID][0] ] * ( (PreReal)1.-barycentricCoords[0]-barycentricCoords[1] ) +
			mesh.vertices[ mesh.triangles[tID][1] ] *               barycentricCoords[0]                        +
			mesh.vertices[ mesh.triangles[tID][2] ] *                                    barycentricCoords[1]   ;

		textureNodePositions[i] = Point3D< float >( p );
	}

	InitializeConcentrations();

	outputBuffer = new unsigned char[ textureHeight*textureWidth];
	memset( outputBuffer , 204 , textureHeight * textureWidth * sizeof(unsigned char) );
}

template< typename PreReal , typename Real >
void _main( int argc , char* argv[] )
{
	GrayScottReactionDiffusion< PreReal , Real >::Init();
	if( !Output.set )
	{
		glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
		GrayScottReactionDiffusion< PreReal , Real >::visualization.displayMode = DisplayMode.value;
		if     ( DisplayMode.value==ONE_REGION_DISPLAY ) GrayScottReactionDiffusion< PreReal , Real >::visualization.screenWidth =  800 , GrayScottReactionDiffusion< PreReal , Real >::visualization.screenHeight = 800;
		else if( DisplayMode.value==TWO_REGION_DISPLAY ) GrayScottReactionDiffusion< PreReal , Real >::visualization.screenWidth = 1440 , GrayScottReactionDiffusion< PreReal , Real >::visualization.screenHeight = 720;
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
		for( int i=0 ; i<OutputSteps.value ; i++ )
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