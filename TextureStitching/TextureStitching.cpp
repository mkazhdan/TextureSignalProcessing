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

enum
{
	SINGLE_INPUT_MODE,
	MULTIPLE_INPUT_MODE,
	INPUT_MODE_COUNT
};

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
#include <Src/MassAndStiffness.h>
#include <Src/Padding.h>
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
#include <Src/StitchingVisualization.h>
#endif // NO_OPEN_GL_VISUALIZATION

using namespace MishaK;

CmdLineParameterArray< std::string , 2 > In( "in" );
CmdLineParameter< std::string > InMask( "mask" );
CmdLineParameter< std::string > InputLowFrequency( "inLow" );
CmdLineParameter< std::string > Output( "out" );
CmdLineParameter< int   > OutputVCycles( "outVCycles" , 6 );
CmdLineParameter< float > InterpolationWeight( "interpolation" , 1e2 );
CmdLineParameter< int   > Levels( "levels" , 4 );
CmdLineParameter< int   > MatrixQuadrature( "mQuadrature" , 6 );

CmdLineParameter< int   > MultigridBlockHeight ( "mBlockH" ,  16 );
CmdLineParameter< int   > MultigridBlockWidth  ( "mBlockW" , 128 );
CmdLineParameter< int   > MultigridPaddedHeight( "mPadH"   ,   0 );
CmdLineParameter< int   > MultigridPaddedWidth ( "mPadW"   ,   2 );

CmdLineParameter< int    > RandomJitter( "jitter" , 0 );
CmdLineParameter< int    > ChartMaskErode( "erode" , 0 );
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
CmdLineParameter< std::string > CameraConfig( "camera" );
#endif // NO_OPEN_GL_VISUALIZATION
CmdLineReadable UseDirectSolver( "useDirectSolver" );
CmdLineReadable Verbose( "verbose" );
CmdLineReadable NoHelp( "noHelp" );
CmdLineReadable DetailVerbose( "detail" );
CmdLineReadable Double( "double" );
CmdLineReadable MultiInput( "multi" );
CmdLineReadable Serial( "serial" );
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
CmdLineReadable Nearest( "nearest" );
#endif // NO_OPEN_GL_VISUALIZATION
CmdLineParameter< double > CollapseEpsilon( "collapse" , 0 );

CmdLineReadable* params[] =
{
	&In , &InMask , &Output , &InterpolationWeight , &Levels , &UseDirectSolver , &Serial, &Verbose ,
	&InputLowFrequency ,
	&DetailVerbose , &MultigridBlockHeight , &MultigridBlockWidth , &MultigridPaddedHeight , &MultigridPaddedWidth , &RandomJitter ,
	&Double , &MatrixQuadrature , &OutputVCycles , &NoHelp , &MultiInput ,
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	&CameraConfig ,
#endif // NO_OPEN_GL_VISUALIZATION
	&ChartMaskErode ,
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	&Nearest ,
#endif // NO_OPEN_GL_VISUALIZATION
	&CollapseEpsilon ,
	NULL
};

void ShowUsage( const char *ex )
{
	printf( "Usage %s:\n" , ex );

	printf( "\t --%s <input mesh and texels>\n" , In.name.c_str() );
	printf( "\t[--%s <input mask>]\n" , InMask.name.c_str() );
	printf( "\t[--%s <input low-frequency texture>\n" , InputLowFrequency.name.c_str() );
#ifdef NO_OPEN_GL_VISUALIZATION
	printf( "\t --%s <output texture>\n" , Output.name.c_str() );
#else // !NO_OPEN_GL_VISUALIZATION
	printf( "\t[--%s <output texture>]\n" , Output.name.c_str() );
#endif // NO_OPEN_GL_VISUALIZATION
	printf( "\t[--%s <chart mask erosion radius>=%d]\n" , ChartMaskErode.name.c_str() , ChartMaskErode.value );
	printf( "\t[--%s <output v-cycles>=%d]\n" , OutputVCycles.name.c_str() , OutputVCycles.value );
	printf( "\t[--%s <interpolation weight>=%f]\n" , InterpolationWeight.name.c_str() , InterpolationWeight.value );
	printf( "\t[--%s <system matrix quadrature points per triangle>=%d]\n" , MatrixQuadrature.name.c_str(), MatrixQuadrature.value );
	printf( "\t[--%s]\n" , UseDirectSolver.name.c_str() );
	printf( "\t[--%s]\n" , MultiInput.name.c_str() );
	printf( "\t[--%s <jittering seed>]\n" , RandomJitter.name.c_str() );
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	printf( "\t[--%s]\n" , Nearest.name.c_str() );
#endif // NO_OPEN_GL_VISUALIZATION
	printf( "\t[--%s]\n" , Verbose.name.c_str() );

#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	printf( "\t[--%s <camera configuration file>]\n" , CameraConfig.name.c_str() );
#endif // NO_OPEN_GL_VISUALIZATION
	printf( "\t[--%s <hierarchy levels>=%d]\n" , Levels.name.c_str(), Levels.value );
	printf( "\t[--%s <multigrid block width>=%d]\n" , MultigridBlockWidth.name.c_str() , MultigridBlockWidth.value );
	printf( "\t[--%s <multigrid block height>=%d]\n" , MultigridBlockHeight.name.c_str() , MultigridBlockHeight.value );
	printf( "\t[--%s <multigrid padded width>=%d]\n" , MultigridPaddedWidth.name.c_str() , MultigridPaddedWidth.value );
	printf( "\t[--%s <multigrid padded height>=%d]\n" , MultigridPaddedHeight.name.c_str() , MultigridPaddedHeight.value );
	printf( "\t[--%s <collapse epsilon>=%g]\n" , CollapseEpsilon.name.c_str() , CollapseEpsilon.value );
	printf( "\t[--%s]\n", Serial.name.c_str() );
	printf( "\t[--%s]\n", DetailVerbose.name.c_str() );
	printf( "\t[--%s]\n" , NoHelp.name.c_str() );
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
class Stitching
{
public:
	static int inputMode;
	static TexturedTriangleMesh< PreReal > mesh;
	static int textureWidth;
	static int textureHeight;
	static Real interpolationWeight;
	static int levels;

	static HierarchicalSystem< PreReal , Real > hierarchy;
	static bool rhsUpdated;
	static bool positiveModulation;

	// Single input mode
	static Image< int > inputMask;
	static Image< Point3D< Real > > lowFrequencyTexture;
	static Image< Point3D< Real > > inputComposition;
	static Image< Point3D< Real > > inputColorMask;

	// Multiple input mode
	static int numTextures;
	static std::vector< Image< Real > > inputConfidence;
	static std::vector< Image< Point3D< Real > > > inputTextures;
	static std::vector< std::vector< Point3D< Real > > > partialTexelValues;
	static std::vector< std::vector< Point3D< Real > > > partialEdgeValues;

	static Image< Point3D< Real > > filteredTexture;

#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	// UI
	static char interpolationStr[1024];
	static char referenceTextureStr[1024];
#endif // NO_OPEN_GL_VISUALIZATION

	static std::vector< Point3D< float > >textureNodePositions;
	static std::vector< Point3D< float > >textureEdgePositions;

	static std::vector< AtlasChart< PreReal > > atlasCharts;

	static std::vector< BilinearElementIndex > bilinearElementIndices;

	static std::vector< TextureNodeInfo< PreReal > > textureNodes;

	static SparseMatrix< Real , int > mass;
	static SparseMatrix< Real , int > stiffness;
	static SparseMatrix< Real , int > stitchingMatrix;

	static std::vector< Point3D< Real > > texelMass;
	static std::vector< Point3D< Real > > texelDivergence;

	static std::vector< SystemCoefficients< Real > > multigridStitchingCoefficients;
	static std::vector< MultigridLevelVariables< Point3D< Real > > > multigridStitchingVariables;
	static std::vector< MultigridLevelIndices< Real > > multigridIndices;

#if defined( USE_CHOLMOD )
	typedef CholmodCholeskySolver< Real , 3 > DirectSolver;
#elif defined( USE_EIGEN )
	typedef EigenCholeskySolver< Real , 3 > DirectSolver;
#elif defined( USE_EIGEN_PARDISO )
	typedef EigenPardisoSolver< Real , 3 > DirectSolver;
#else
#error "[ERROR] No solver defined!"
#endif

	static VCycleSolvers< DirectSolver > vCycleSolvers;
	static DirectSolver directSolver;

	static std::map< SimplexIndex< 1 > , unsigned int > edgeIndex;
	static std::vector< SimplexIndex< 1 > > edgePairs;

	static SparseMatrix< Real , int > boundaryDivergenceMatrix;

	static std::vector< Real > deepDivergenceCoefficients;
	static std::vector< DivegenceRasterLine > divergenceRasterLines;

	static std::vector< bool > unobservedTexel;
	static std::vector< Point3D< Real > > texelValues;
	static std::vector< Point3D< Real > > edgeValues;

	// Linear Operators
	static SystemCoefficients< Real > massCoefficients;
	static SystemCoefficients< Real > stiffnessCoefficients;

	// Stitching UI
	static int textureIndex;

	static int steps;
	static char stepsString[];

	static Padding padding;

#ifdef NO_OPEN_GL_VISUALIZATION
	static unsigned int updateCount;
	static void WriteTexture( const char *fileName );
#else // !NO_OPEN_GL_VISUALIZATION
	// Visulization
	static StitchingVisualization visualization;
	static unsigned int updateCount;

	static void ToggleForwardReferenceTextureCallBack ( Visualization *v , const char *prompt );
	static void ToggleBackwardReferenceTextureCallBack( Visualization *v , const char *prompt );
	static void ToggleMaskCallBack                    ( Visualization *v , const char *prompt );
	static void ToggleUpdateCallBack                  ( Visualization *v , const char *prompt );
	static void IncrementUpdateCallBack               ( Visualization *v , const char *prompt );
	static void ExportTextureCallBack                 ( Visualization *v , const char *prompt );
	static void InterpolationWeightCallBack           ( Visualization *v , const char *prompt );
#endif // NO_OPEN_GL_VISUALIZATION

	static Image< Point3D< unsigned char > > GetChartMask( void );
	static void LoadTextures( void );
	static void LoadMasks( void );
	static void ParseImages( void );
	static void SetUpSystem( void );
	static void SolveSystem( void );
	static void Init( void );
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	static void InitializeVisualization( void );
#endif // NO_OPEN_GL_VISUALIZATION
	static void UpdateSolution( bool verbose=false , bool detailVerbose=false );
	static void ComputeExactSolution( bool verbose=false );
	static void InitializeSystem( int width , int height );
	static void _InitializeSystem( std::vector< std::vector< SquareMatrix< PreReal , 2 > > > &parameterMetric , BoundaryProlongationData< Real > &boundaryProlongation );

#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	static void UpdateFilteredColorTexture( const std::vector< Point3D< Real > > &solution );
#endif // NO_OPEN_GL_VISUALIZATION
	static void UpdateFilteredTexture( const std::vector< Point3D< Real > > &solution );

#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	static void Display( void ){ visualization.Display(); }
	static void MouseFunc( int button , int state , int x , int y );
	static void MotionFunc( int x , int y );
	static void Reshape( int w , int h ) { visualization.Reshape(w,h); }
	static void KeyboardFunc( unsigned char key , int x , int y ) { visualization.KeyboardFunc(key,x,y); }
	static void Idle( void );
#endif // NO_OPEN_GL_VISUALIZATION
};


#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
template< typename PreReal , typename Real , unsigned int TextureBitDepth > char															Stitching< PreReal , Real , TextureBitDepth >::referenceTextureStr[1024];
template< typename PreReal , typename Real , unsigned int TextureBitDepth > char															Stitching< PreReal , Real , TextureBitDepth >::interpolationStr[1024];
#endif // NO_OPEN_GL_VISUALIZATION

template< typename PreReal , typename Real , unsigned int TextureBitDepth > int																Stitching< PreReal , Real , TextureBitDepth >::inputMode;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > TexturedTriangleMesh< PreReal >									Stitching< PreReal , Real , TextureBitDepth >::mesh;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > int																Stitching< PreReal , Real , TextureBitDepth >::textureWidth;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > int																Stitching< PreReal , Real , TextureBitDepth >::textureHeight;
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
template< typename PreReal , typename Real , unsigned int TextureBitDepth > StitchingVisualization											Stitching< PreReal , Real , TextureBitDepth >::visualization;
#endif // NO_OPEN_GL_VISUALIZATION
template< typename PreReal , typename Real , unsigned int TextureBitDepth > SparseMatrix< Real , int >										Stitching< PreReal , Real , TextureBitDepth >::mass;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > SparseMatrix< Real , int >										Stitching< PreReal , Real , TextureBitDepth >::stiffness;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > SparseMatrix< Real , int >										Stitching< PreReal , Real , TextureBitDepth >::stitchingMatrix;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > Real															Stitching< PreReal , Real , TextureBitDepth >::interpolationWeight;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< TextureNodeInfo< PreReal > >						Stitching< PreReal , Real , TextureBitDepth >::textureNodes;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< BilinearElementIndex >								Stitching< PreReal , Real , TextureBitDepth >::bilinearElementIndices;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > int																Stitching< PreReal , Real , TextureBitDepth >::steps;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > char															Stitching< PreReal , Real , TextureBitDepth >::stepsString[1024];
template< typename PreReal , typename Real , unsigned int TextureBitDepth > int																Stitching< PreReal , Real , TextureBitDepth >::levels;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > HierarchicalSystem< PreReal , Real >							Stitching< PreReal , Real , TextureBitDepth >::hierarchy;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > bool															Stitching< PreReal , Real , TextureBitDepth >::rhsUpdated = true;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > bool															Stitching< PreReal , Real , TextureBitDepth >::positiveModulation = true;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > Image< Point3D< Real > >										Stitching< PreReal , Real , TextureBitDepth >::filteredTexture;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > Image< int >												    Stitching< PreReal , Real , TextureBitDepth >::inputMask;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > Image< Point3D< Real > >										Stitching< PreReal , Real , TextureBitDepth >::lowFrequencyTexture;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > Image< Point3D< Real > >										Stitching< PreReal , Real , TextureBitDepth >::inputComposition;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > Image< Point3D< Real > >										Stitching< PreReal , Real , TextureBitDepth >::inputColorMask;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > int																Stitching< PreReal , Real , TextureBitDepth >::numTextures;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Image< Real > >									Stitching< PreReal , Real , TextureBitDepth >::inputConfidence;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Image< Point3D< Real > > >							Stitching< PreReal , Real , TextureBitDepth >::inputTextures;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Point3D< Real > >									Stitching< PreReal , Real , TextureBitDepth >::texelMass;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Point3D< Real > >									Stitching< PreReal , Real , TextureBitDepth >::texelDivergence;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< SystemCoefficients< Real > >						Stitching< PreReal , Real , TextureBitDepth >::multigridStitchingCoefficients;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< MultigridLevelVariables< Point3D< Real > > >		Stitching< PreReal , Real , TextureBitDepth >::multigridStitchingVariables;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< MultigridLevelIndices< Real > >					Stitching< PreReal , Real , TextureBitDepth >::multigridIndices;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > VCycleSolvers< typename Stitching< PreReal , Real , TextureBitDepth >::DirectSolver >		Stitching< PreReal , Real , TextureBitDepth >::vCycleSolvers;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > typename Stitching< PreReal , Real , TextureBitDepth >::DirectSolver				Stitching< PreReal , Real , TextureBitDepth >::directSolver;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< AtlasChart< PreReal > >							Stitching< PreReal , Real , TextureBitDepth >::atlasCharts;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Point3D< float > >									Stitching< PreReal , Real , TextureBitDepth >::textureNodePositions;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Point3D< float > >									Stitching< PreReal , Real , TextureBitDepth >::textureEdgePositions;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > Padding															Stitching< PreReal , Real , TextureBitDepth >::padding;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > SystemCoefficients< Real >										Stitching< PreReal , Real , TextureBitDepth >::massCoefficients;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > SystemCoefficients< Real >										Stitching< PreReal , Real , TextureBitDepth >::stiffnessCoefficients;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > unsigned int													Stitching< PreReal , Real , TextureBitDepth >::updateCount = static_cast< unsigned int >(-1);

template< typename PreReal , typename Real , unsigned int TextureBitDepth >  std::map< SimplexIndex< 1 > , unsigned int >					Stitching< PreReal , Real , TextureBitDepth >::edgeIndex;
template< typename PreReal , typename Real , unsigned int TextureBitDepth >  std::vector< SimplexIndex< 1 > >								Stitching< PreReal , Real , TextureBitDepth >::edgePairs;

template< typename PreReal , typename Real , unsigned int TextureBitDepth >  SparseMatrix< Real , int >										Stitching< PreReal , Real , TextureBitDepth >::boundaryDivergenceMatrix;
template< typename PreReal , typename Real , unsigned int TextureBitDepth >  std::vector< Real >											Stitching< PreReal , Real , TextureBitDepth >::deepDivergenceCoefficients;
template< typename PreReal , typename Real , unsigned int TextureBitDepth >  std::vector< DivegenceRasterLine >  							Stitching< PreReal , Real , TextureBitDepth >::divergenceRasterLines;

template< typename PreReal , typename Real , unsigned int TextureBitDepth >  std::vector< bool >											Stitching< PreReal , Real , TextureBitDepth >::unobservedTexel;
template< typename PreReal , typename Real , unsigned int TextureBitDepth >  std::vector< Point3D< Real > >									Stitching< PreReal , Real , TextureBitDepth >::texelValues;
template< typename PreReal , typename Real , unsigned int TextureBitDepth >  std::vector< Point3D< Real > >									Stitching< PreReal , Real , TextureBitDepth >::edgeValues;
template< typename PreReal , typename Real , unsigned int TextureBitDepth >  std::vector< std::vector< Point3D< Real > > >					Stitching< PreReal , Real , TextureBitDepth >::partialTexelValues;
template< typename PreReal , typename Real , unsigned int TextureBitDepth >  std::vector< std::vector< Point3D< Real > > >					Stitching< PreReal , Real , TextureBitDepth >::partialEdgeValues;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > int																Stitching< PreReal , Real , TextureBitDepth >::textureIndex = 0;

#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::UpdateFilteredColorTexture( const std::vector< Point3D< Real > > &solution )
{
	ThreadPool::ParallelFor
		(
			0 , textureNodes.size() ,
			[&]( unsigned int , size_t i )
			{
				int ci = textureNodes[i].ci;
				int cj = textureNodes[i].cj;
				int offset = 3 * ( textureWidth*cj + ci );
				for( int c=0 ; c<3 ; c++ )
				{
					Real value = std::min< Real >( (Real)1 , std::max< Real >( (Real)0. , solution[i][c] ) );
					visualization.colorTextureBuffer[offset + c] = (unsigned char)(value*255.0);
				}
			}
		);
}
#endif // NO_OPEN_GL_VISUALIZATION

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::UpdateFilteredTexture( const std::vector< Point3D< Real > > &solution )
{
	ThreadPool::ParallelFor
		(
			0 , textureNodes.size() ,
			[&]( unsigned int , size_t i )
			{
				int ci = textureNodes[i].ci , cj = textureNodes[i].cj;
				filteredTexture(ci,cj) = solution[i];
			}
		);
}

#ifdef NO_OPEN_GL_VISUALIZATION
template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::WriteTexture( const char *fileName )
{
	UpdateFilteredTexture( multigridStitchingVariables[0].x );
	Image< Point3D< Real > > outputTexture = filteredTexture;
	padding.unpad( outputTexture );
	WriteImage< TextureBitDepth >( outputTexture , fileName );
}

#else // !NO_OPEN_GL_VISUALIZATION
template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::Idle( void )
{
	visualization.Idle();

	float radius = 0.03;
	float radiusSquared = radius * radius;
	if( inputMode==MULTIPLE_INPUT_MODE )
	{
		if( visualization.isBrushActive )
		{
			Point3D< float > selectedPoint;
			bool validSelection = false;
			if( visualization.showMesh ) validSelection = visualization.select( visualization.diskX , visualization.diskY , selectedPoint );
			if( validSelection )
			{
				ThreadPool::ParallelFor( 0 , textureNodePositions.size() , [&]( unsigned int , size_t i ){ if( Point3D< float >::SquareNorm( textureNodePositions[i]-selectedPoint )<radiusSquared ) texelValues[i] = partialTexelValues[ textureIndex ][i]; } );
				ThreadPool::ParallelFor( 0 , textureNodePositions.size() , [&]( unsigned int , size_t i ){ if( Point3D< float >::SquareNorm( textureEdgePositions[i]-selectedPoint )<radiusSquared )  edgeValues[i] =  partialEdgeValues[ textureIndex ][i]; } );
				rhsUpdated = false;
			}
		}
	}

	if( updateCount && !UseDirectSolver.set && !visualization.promptCallBack )
	{
		UpdateSolution();
		UpdateFilteredColorTexture( multigridStitchingVariables[0].x );
		visualization.UpdateColorTextureBuffer();
		if( updateCount>0 ) updateCount--;
		steps++;
		sprintf( stepsString , "Steps: %d" , steps );
	}

	glutPostRedisplay();
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::MouseFunc( int button , int state , int x , int y )
{
	visualization.newX = x; visualization.newY = y;
	visualization.rotating = visualization.scaling = visualization.panning = false;
	visualization.isBrushActive = false;

	if( state==GLUT_DOWN && glutGetModifiers() & GLUT_ACTIVE_SHIFT )
	{
		visualization.isBrushActive = true;
		visualization.diskX = x;
		visualization.diskY = y;

		if     ( button==GLUT_RIGHT_BUTTON ) positiveModulation = true;
		else if( button==GLUT_LEFT_BUTTON  ) positiveModulation = false;
	}
	else
	{
		if( visualization.showMesh )
		{
			visualization.newX = x , visualization.newY = y;

			visualization.rotating = visualization.scaling = visualization.panning = false;
			if( ( button==GLUT_LEFT_BUTTON || button==GLUT_RIGHT_BUTTON ) && glutGetModifiers() & GLUT_ACTIVE_CTRL ) visualization.panning = true;
			else if( button==GLUT_LEFT_BUTTON  ) visualization.rotating = true;
			else if( button==GLUT_RIGHT_BUTTON ) visualization.scaling  = true;
		}
		else
		{
			if( button==GLUT_LEFT_BUTTON  ) visualization.panning = true;
			else if( button==GLUT_RIGHT_BUTTON ) visualization.scaling  = true;
		}

	}

	glutPostRedisplay();
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::MotionFunc( int x , int y )
{
	if( visualization.isBrushActive )
	{
		visualization.diskX = x;
		visualization.diskY = y;
	}
	else
	{
		if( visualization.showMesh )
		{
			visualization.oldX = visualization.newX , visualization.oldY = visualization.newY , visualization.newX = x , visualization.newY = y;
			int screenSize = std::min< int >( visualization.screenWidth , visualization.screenHeight );
			float rel_x = ( visualization.newX-visualization.oldX ) / (float)screenSize * 2;
			float rel_y = ( visualization.newY-visualization.oldY ) / (float)screenSize * 2;

			float pRight   =  rel_x * visualization.zoom , pUp = -rel_y * visualization.zoom;
			float pForward =  rel_y * visualization.zoom;
			float rRight   = -rel_y , rUp = -rel_x;

			if     ( visualization.rotating ) visualization.camera.rotateUp( -rUp ) , visualization.camera.rotateRight( -rRight );
			else if( visualization.scaling  ) visualization.camera.translate( visualization.camera.forward*pForward );
			else if( visualization.panning  ) visualization.camera.translate( -( visualization.camera.right*pRight + visualization.camera.up*pUp ) );
		}
		else
		{
			visualization.oldX = visualization.newX , visualization.oldY = visualization.newY , visualization.newX = x , visualization.newY = y;
			if( visualization.panning ) visualization.xForm.offset[0] += ( visualization.newX - visualization.oldX ) , visualization.xForm.offset[1] -= ( visualization.newY - visualization.oldY );
			else if( visualization.scaling )
			{
				float dz = (float) pow( 1.1 , (double)( visualization.newY - visualization.oldY)/8 );
				visualization.xForm.zoom *= dz;
			}
		}
	}
	glutPostRedisplay();
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth > void Stitching< PreReal , Real , TextureBitDepth >::ToggleMaskCallBack( Visualization * /*v*/ , const char * )
{
	visualization.showMask = !visualization.showMask;
	glutPostRedisplay();
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth > void Stitching< PreReal , Real , TextureBitDepth >::ToggleForwardReferenceTextureCallBack( Visualization * /*v*/ , const char* )
{
	textureIndex = ( textureIndex+1 ) % numTextures;
	visualization.referenceIndex = textureIndex;
	sprintf( referenceTextureStr , "Reference Texture: %02d of %02d \n" , textureIndex , numTextures );
	glutPostRedisplay();
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth > void Stitching< PreReal , Real , TextureBitDepth >::ToggleBackwardReferenceTextureCallBack( Visualization * /*v*/ , const char* )
{
	textureIndex = ( textureIndex+numTextures-1 ) % numTextures;
	visualization.referenceIndex = textureIndex;
	sprintf( referenceTextureStr, "Reference Texture: %02d of %02d \n" , textureIndex , numTextures );
	glutPostRedisplay();
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth > void Stitching< PreReal , Real , TextureBitDepth >::ToggleUpdateCallBack( Visualization * /*v*/ , const char * )
{
	if( updateCount ) updateCount = 0;
	else              updateCount = -1;
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::IncrementUpdateCallBack( Visualization * /*v*/ , const char * )
{
	if( updateCount<0 ) updateCount = 1;
	else                updateCount++;
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::ExportTextureCallBack( Visualization * /*v*/ , const char *prompt )
{
	UpdateFilteredTexture( multigridStitchingVariables[0].x );
	Image< Point3D< Real > > outputTexture = filteredTexture;
	padding.unpad( outputTexture );
	WriteImage< TextureBitDepth >( outputTexture , prompt );
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void  Stitching< PreReal , Real , TextureBitDepth >::InterpolationWeightCallBack( Visualization * /*v*/ , const char *prompt )
{
	interpolationWeight = atof(prompt);
	if( UseDirectSolver.set ) stitchingMatrix = mass*interpolationWeight + stiffness;
	Miscellany::Timer timer;
	UpdateLinearSystem( interpolationWeight , (Real)1. , hierarchy , multigridStitchingCoefficients , massCoefficients , stiffnessCoefficients , vCycleSolvers , directSolver , stitchingMatrix , DetailVerbose.set , false , UseDirectSolver.set );
	if( Verbose.set ) printf( "\tInitialized multigrid coefficients: %.2f(s)\n" , timer.elapsed() );

	ThreadPool::ParallelFor( 0 , multigridStitchingVariables[0].rhs.size() , [&]( unsigned int , size_t i ){ multigridStitchingVariables[0].rhs[i] = texelMass[i] * interpolationWeight + texelDivergence[i]; } );

	if( UseDirectSolver.set ) ComputeExactSolution(Verbose.set);
	else for( int i=0 ; i<OutputVCycles.value ; i++ ) VCycle( multigridStitchingVariables , multigridStitchingCoefficients , multigridIndices , vCycleSolvers , false , false );

	UpdateFilteredColorTexture( multigridStitchingVariables[0].x );
	visualization.UpdateColorTextureBuffer();
	sprintf( interpolationStr , "Interpolation weight: %e\n" , interpolationWeight );
}
#endif // NO_OPEN_GL_VISUALIZATION


template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::ComputeExactSolution( bool verbose )
{
	Miscellany::Timer timer;
	solve( directSolver , multigridStitchingVariables[0].x , multigridStitchingVariables[0].rhs );
	if( verbose ) printf( "Solving time =  %.4f\n" , timer.elapsed() );
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::UpdateSolution( bool verbose , bool detailVerbose )
{
	if( !rhsUpdated )
	{
		Miscellany::Timer timer;

		MultiplyBySystemMatrix_NoReciprocals( massCoefficients , hierarchy.gridAtlases[0].indexConverter , hierarchy.gridAtlases[0].rasterLines , texelValues , texelMass );
		ComputeDivergence( edgeValues , texelDivergence , deepDivergenceCoefficients , boundaryDivergenceMatrix , divergenceRasterLines );

		ThreadPool::ParallelFor( 0 , textureNodes.size() , [&]( unsigned int , size_t i ){ multigridStitchingVariables[0].rhs[i] = texelMass[i] * interpolationWeight + texelDivergence[i]; } );

		if( verbose ) printf( "RHS update time %.4f\n" , timer.elapsed() );
		rhsUpdated = true;
	}

	VCycle( multigridStitchingVariables , multigridStitchingCoefficients , multigridIndices , vCycleSolvers , verbose , detailVerbose );
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::_InitializeSystem( std::vector< std::vector< SquareMatrix< PreReal , 2 > > > &parameterMetric , BoundaryProlongationData< Real > &boundaryProlongation )
{
	// Unused parameters
	std::vector< Point3D< Real > > inputSignal;
	std::vector< Real > texelToCellCoeffs;
	SparseMatrix< Real , int > boundaryCellBasedStiffnessRHSMatrix[3];

	{
		switch( MatrixQuadrature.value )
		{
		case  1: InitializeMassAndStiffness< 1>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix , true , edgeIndex , boundaryDivergenceMatrix , deepDivergenceCoefficients ) ; break;
		case  3: InitializeMassAndStiffness< 3>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix , true , edgeIndex , boundaryDivergenceMatrix , deepDivergenceCoefficients ) ; break;
		case  6: InitializeMassAndStiffness< 6>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix , true , edgeIndex , boundaryDivergenceMatrix , deepDivergenceCoefficients ) ; break;
		case 12: InitializeMassAndStiffness<12>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix , true , edgeIndex , boundaryDivergenceMatrix , deepDivergenceCoefficients ) ; break;
		case 24: InitializeMassAndStiffness<24>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix , true , edgeIndex , boundaryDivergenceMatrix , deepDivergenceCoefficients ) ; break;
		case 32: InitializeMassAndStiffness<32>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix , true , edgeIndex , boundaryDivergenceMatrix , deepDivergenceCoefficients ) ; break;
		default: MK_THROW( "Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles" );
		}
	}
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::InitializeSystem( int width , int height )
{
	Miscellany::Timer timer;
	MultigridBlockInfo multigridBlockInfo( MultigridBlockWidth.value , MultigridBlockHeight.value , MultigridPaddedWidth.value , MultigridPaddedHeight.value );
	InitializeHierarchy( mesh , width , height , levels , textureNodes , bilinearElementIndices , hierarchy , atlasCharts , multigridBlockInfo , false );
	if( Verbose.set ) printf( "\tInitialized hierarchy: %.2f(s)\n" , timer.elapsed() );

	BoundaryProlongationData< Real > boundaryProlongation;
	InitializeBoundaryProlongationData( hierarchy.gridAtlases[0] , boundaryProlongation );

	std::vector< Point3D< Real > > inputSignal;
	std::vector< Real > texelToCellCoeffs;

	std::vector< std::vector< SquareMatrix< PreReal , 2 > > > parameterMetric;
	InitializeMetric( mesh , EMBEDDING_METRIC , atlasCharts , parameterMetric );
	InitializeIntraChartEdgeIndexing( hierarchy.gridAtlases[0].gridCharts , edgeIndex );

	if( Verbose.set ) timer.reset();
	_InitializeSystem( parameterMetric , boundaryProlongation );
	if( Verbose.set ) printf( "\tInitialized system: %.2f(s)\n" , timer.elapsed() );
	
	InitializeDivergenceRasteLines( edgeIndex , hierarchy.gridAtlases[0].rasterLines , divergenceRasterLines );
	
	edgePairs.resize( edgeIndex.size() );
	for( auto edgeIter=edgeIndex.begin() ; edgeIter!=edgeIndex.end() ; edgeIter++ ) edgePairs[ (*edgeIter).second ] = (*edgeIter).first;

	texelMass.resize( textureNodes.size() );
	texelDivergence.resize( textureNodes.size() );

	if( UseDirectSolver.set )
	{
		FullMatrixConstruction( hierarchy.gridAtlases[0] , massCoefficients , mass );
		FullMatrixConstruction( hierarchy.gridAtlases[0] , stiffnessCoefficients , stiffness );
		stitchingMatrix = mass*interpolationWeight + stiffness;
	}

	multigridIndices.resize(levels);
	for( unsigned int i=0 ; i<levels ; i++ )
	{
		const typename GridAtlas<>::IndexConverter & indexConverter = hierarchy.gridAtlases[i].indexConverter;
		const GridAtlas< PreReal , Real > &gridAtlas = hierarchy.gridAtlases[i];
		multigridIndices[i].threadTasks = gridAtlas.threadTasks;
		multigridIndices[i].boundaryToSupported = indexConverter.boundaryToSupported();
		multigridIndices[i].segmentedLines = gridAtlas.segmentedLines;
		multigridIndices[i].rasterLines = gridAtlas.rasterLines;
		multigridIndices[i].restrictionLines = gridAtlas.restrictionLines;
		multigridIndices[i].prolongationLines = gridAtlas.prolongationLines;
		if( i<levels-1 ) multigridIndices[i].boundaryRestriction = hierarchy.boundaryRestriction[i];
	}

	if( Verbose.set ) timer.reset();
	UpdateLinearSystem( interpolationWeight , (Real)1. , hierarchy , multigridStitchingCoefficients , massCoefficients , stiffnessCoefficients , vCycleSolvers , directSolver , stitchingMatrix , DetailVerbose.set, true, UseDirectSolver.set );
	if( Verbose.set ) printf( "\tInitialized multigrid coefficients: %.2f(s)\n" , timer.elapsed() );

	multigridStitchingVariables.resize(levels);
	for( unsigned int i=0 ; i<levels ; i++ )
	{
		const typename GridAtlas<>::IndexConverter & indexConverter = hierarchy.gridAtlases[i].indexConverter;
		MultigridLevelVariables< Point3D< Real > >& variables = multigridStitchingVariables[i];
		variables.x.resize( hierarchy.gridAtlases[i].numTexels );
		variables.rhs.resize( hierarchy.gridAtlases[i].numTexels );
		variables.residual.resize( hierarchy.gridAtlases[i].numTexels );
		variables.boundary_rhs.resize( indexConverter.numBoundary() );
		variables.boundary_value.resize( indexConverter.numBoundary() );
		variables.variable_boundary_value.resize( indexConverter.numBoundary() );
	}
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::SetUpSystem( void )
{
	texelMass.resize( textureNodes.size() );
	MultiplyBySystemMatrix_NoReciprocals( massCoefficients , hierarchy.gridAtlases[0].indexConverter , hierarchy.gridAtlases[0].rasterLines , texelValues , texelMass );

	ComputeDivergence( edgeValues , texelDivergence , deepDivergenceCoefficients , boundaryDivergenceMatrix , divergenceRasterLines );

	ThreadPool::ParallelFor
		(
			0 , texelValues.size() ,
			[&]( unsigned int , size_t i )
			{
				multigridStitchingVariables[0].x[i] = texelValues[i];
				multigridStitchingVariables[0].rhs[i] = texelMass[i] * interpolationWeight + texelDivergence[i];
			}
		);

	filteredTexture.resize( textureWidth , textureHeight );
	for( int i=0 ; i<filteredTexture.size() ; i++ ) filteredTexture[i] = Point3D< Real >( (Real)0.5 , (Real)0.5 , (Real)0.5 );

	UpdateFilteredTexture( multigridStitchingVariables[0].x );
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::SolveSystem( void )
{
	if( UseDirectSolver.set ) ComputeExactSolution( Verbose.set );
	else for( int i=0 ; i<OutputVCycles.value ; i++ ) VCycle( multigridStitchingVariables , multigridStitchingCoefficients , multigridIndices , vCycleSolvers , false , false );	
	UpdateFilteredTexture( multigridStitchingVariables[0].x );
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
Image< Point3D< unsigned char > > Stitching< PreReal , Real , TextureBitDepth >::GetChartMask( void )
{
	auto IsBlack = []( Point3D< unsigned char > c ){ return !c[0] && !c[1] && !c[2]; };

	Image< Point3D< unsigned char > > chartMask;
	chartMask.resize( textureWidth , textureHeight );
	for( unsigned int i=0 ; i<(unsigned int)textureWidth ; i++ ) for( unsigned int j=0 ; j<(unsigned int)textureHeight ; j++ ) chartMask(i,j) = Point3D< unsigned char >(0,0,0);

	std::map< int , Point3D< unsigned char > > chartColors;
	auto ColorAlreadyUsed = [&]( Point3D< unsigned char > c )
		{
			for( auto iter=chartColors.begin() ; iter!=chartColors.end() ; iter++ )
				if( c[0]==iter->second[0] && c[1]==iter->second[1] && c[2]==iter->second[2] ) return true;
			return false;
		};

	for( unsigned int i=0 ; i<textureNodes.size() ; i++ ) if( chartColors.find( textureNodes[i].chartID )==chartColors.end() )
	{
		Point3D< unsigned char > c;
		while( IsBlack(c) || ColorAlreadyUsed(c) ) c[0] = rand()%256 , c[1] = rand()%256 , c[2] = rand()%256;
		chartColors[ textureNodes[i].chartID ] = c;
	}

	for( int i=0 ; i<textureNodes.size() ; i++ ) chartMask( textureNodes[i].ci , textureNodes[i].cj ) = chartColors[ textureNodes[i].chartID ];
	for( int e=0 ; e<ChartMaskErode.value ; e++ )
	{
		Image< Point3D< unsigned char > > _chartMask = chartMask;
		for( int i=0 ; i<textureWidth ; i++ ) for( int j=0 ; j<textureHeight ; j++ ) if( chartMask(i,j)[0] || chartMask(i,j)[1] || chartMask(i,j)[2] )
			for( int di=-1 ; di<=1 ; di++ ) for( int dj=-1 ; dj<=1 ; dj++ )
				if( i+di>=0 && i+di<textureWidth && j+dj>0 && j+dj<textureHeight ) if( IsBlack( _chartMask(i+di,j+dj) ) ) chartMask(i,j) = Point3D< unsigned char >();
	}
	return chartMask;
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::LoadTextures( void )
{
	if( inputMode==MULTIPLE_INPUT_MODE )
	{
		bool countingTextures = true;
		numTextures = 0;
		while( countingTextures )
		{
			char textureName[256];
			sprintf( textureName , In.values[1].c_str() , numTextures );
			FILE * file = fopen( textureName , "r" );
			if( file ) numTextures++;
			else countingTextures = false;
		}

		inputTextures.resize( numTextures );
		ThreadPool::ParallelFor
			(
				0 , numTextures ,
				[&]( unsigned int , size_t i )
				{
					char textureName[256];
					sprintf( textureName , In.values[1].c_str() , i );
					ReadImage< TextureBitDepth >( inputTextures[i] , textureName );
				}
			);

		textureWidth = inputTextures[0].res(0);
		textureHeight = inputTextures[0].res(1);

		if( Verbose.set ) printf( "Texture count: %d\n" , numTextures );
	}
	else
	{
		ReadImage< TextureBitDepth >( inputComposition , In.values[1] );
		textureWidth = inputComposition.res(0);
		textureHeight = inputComposition.res(1);
	}
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::LoadMasks( void )
{
	if( inputMode==MULTIPLE_INPUT_MODE )
	{
		inputConfidence.resize( numTextures );
		ThreadPool::ParallelFor
			(
				0 , numTextures ,
				[&]( unsigned int , size_t i )
				{
					char confidenceName[256];
					sprintf( confidenceName , InMask.value.c_str() , i );
					Image< Point3D< Real > > textureConfidence;
					ReadImage< 8 >( textureConfidence , confidenceName );
					inputConfidence[i].resize( textureWidth , textureHeight );
					for( int p=0 ; p<textureConfidence.size() ; p++ ) inputConfidence[i][p] = Point3D< Real >::Dot( textureConfidence[p] , Point3D< Real >( (Real)1./3 , (Real)1./3 , (Real)1./3 ) );
				}
			);
	}
	else
	{
		Image< Point3D< unsigned char > > textureConfidence;
		if( InMask.set ) ReadImage< 8 >( textureConfidence , InMask.value );
		else textureConfidence = GetChartMask();

		inputColorMask.resize( textureConfidence.res() );
		for( unsigned int i=0 ; i<textureConfidence.res(0) ; i++ ) for( unsigned int j=0 ; j<textureConfidence.res(1) ; j++ ) for( unsigned int c=0 ; c<3 ; c++ )
			inputColorMask(i,j)[c] = ( (Real)textureConfidence(i,j)[c] )/255;

		inputMask.resize( textureWidth , textureHeight );
		for( int p=0 ; p<textureConfidence.size() ; p++ )
		{
			int index = ( (int)textureConfidence[p][0] ) * 256 * 256 + ( (int)textureConfidence[p][1] ) * 256 + ( (int)textureConfidence[p][2] );
			inputMask[p] = index==0 ? -1 : index;
		}
	}
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::ParseImages( void )
{
	int numNodes = (int)textureNodes.size();
	int numEdges = (int)edgeIndex.size();
	unobservedTexel.resize( numNodes );
	texelValues.resize( numNodes );
 	edgeValues.resize( numEdges );

	if( inputMode==MULTIPLE_INPUT_MODE )
	{
		std::vector< Real > texelWeight;
		std::vector< Real > edgeWeight;
		texelWeight.resize( numNodes , 0 );
		edgeWeight.resize( numEdges , 0 );

		partialTexelValues.resize( numTextures );
		for( int i=0 ; i<numTextures ; i++ ) partialTexelValues[i].resize( numNodes );
		partialEdgeValues.resize( numTextures );
		for( int i=0 ; i<numTextures ; i++ ) partialEdgeValues[i].resize( edgeIndex.size() );

		for( int textureIter=0 ; textureIter<numTextures ; textureIter++ )
		{
			Image< Point3D< Real > > textureValues = InputLowFrequency.set ? lowFrequencyTexture : inputTextures[textureIter];
			Image< Real > textureConfidence = inputConfidence[textureIter];

			for( int i=0 ; i<textureNodes.size() ; i++ )
			{
				Real weight = textureConfidence( textureNodes[i].ci , textureNodes[i].cj );
				Point3D< Real > value = textureValues( textureNodes[i].ci , textureNodes[i].cj );
				partialTexelValues[textureIter][i] = value;
				texelValues[i] += value * weight;
				texelWeight[i] += weight;
			}

			for( int e=0 ; e<edgePairs.size() ; e++ )
			{

				Real weight[2];
				Point3D< Real > value[2];
				for( int k=0 ; k<2 ; k++ )
				{
					unsigned int i = edgePairs[e][k];
					weight[k] = textureConfidence( textureNodes[i].ci , textureNodes[i].cj );
					value[k] = textureValues( textureNodes[i].ci , textureNodes[i].cj );
				}
				Real eWeight = weight[0] * weight[1];
				Point3D<Real> eValue = value[1] - value[0];
				partialEdgeValues[textureIter][e] = eWeight > 0 ? eValue : Point3D<Real>();
				edgeValues[e] += eValue*eWeight;
				edgeWeight[e] += eWeight;
			}

		}
		for( int i=0 ; i<numNodes ; i++ )
		{
			unobservedTexel[i] = texelWeight[i]==0;
			if( texelWeight[i]>0 ) texelValues[i] /= texelWeight[i];
		}

		for( int i=0 ; i<numEdges ; i++ ) if( edgeWeight[i]>0 ) edgeValues[i] /= edgeWeight[i];
	}
	else
	{
		for( int i=0 ; i<textureNodes.size() ; i++ )
		{
			unobservedTexel[i] = inputMask( textureNodes[i].ci , textureNodes[i].cj )==-1;
			texelValues[i] = InputLowFrequency.set ? lowFrequencyTexture( textureNodes[i].ci , textureNodes[i].cj ) : inputComposition( textureNodes[i].ci , textureNodes[i].cj );
		}
		for( int e=0 ; e<edgePairs.size() ; e++ )
		{
			const SimplexIndex< 1 > &edgeCorners = edgePairs[e];
			unsigned int ci[] = { textureNodes[ edgeCorners[0] ].ci , textureNodes[ edgeCorners[1] ].ci };
			unsigned int cj[] = { textureNodes[ edgeCorners[0] ].cj , textureNodes[ edgeCorners[1] ].cj };
			if( inputMask( ci[0] , cj[0] )!=-1 && inputMask( ci[0] , cj[0] )==inputMask( ci[1] , cj[1] ) ) edgeValues[e] = inputComposition( ci[1] , cj[1] ) - inputComposition( ci[0] , cj[0] );
			else edgeValues[e] = Point3D< Real >(0, 0, 0);
		}
	}

	// Set unobserved texel values to be the average of observed ones
	Point3D< Real > avgObservedTexel;
	int observedTexelCount = 0;
	for( int i=0 ; i<unobservedTexel.size() ; i++ )
	{
		if( !unobservedTexel[i] )
		{
			avgObservedTexel += texelValues[i];
			observedTexelCount++;
		}
	}
	avgObservedTexel /= (Real)observedTexelCount;
	for( int i=0 ; i<unobservedTexel.size() ; i++ ) if( unobservedTexel[i] ) texelValues[i] = avgObservedTexel;

	ThreadPool::ParallelFor( 0 , textureNodes.size() , [&]( unsigned int , size_t i ){ multigridStitchingVariables[0].x[i] = texelValues[i]; } );
}

#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::InitializeVisualization( void )
{
	sprintf( referenceTextureStr , "Reference Texture: %02d of %02d\n" , textureIndex,numTextures );
	sprintf( interpolationStr , "Interpolation: %.2e\n" , interpolationWeight );

	visualization.textureWidth = textureWidth;
	visualization.textureHeight = textureHeight;

	visualization.colorTextureBuffer = new unsigned char[textureHeight*textureWidth * 3];
	memset( visualization.colorTextureBuffer , 204 , textureHeight * textureWidth * 3 * sizeof(unsigned char) );

	unsigned int tCount = (unsigned int)mesh.numTriangles();

	visualization.triangles.resize( tCount );
	visualization.vertices.resize( 3*tCount );
	visualization.colors.resize( 3*tCount , Point3D< double >( 0.75 , 0.75 , 0.75 ) );
	visualization.textureCoordinates.resize( 3*tCount );
	visualization.normals.resize( 3*tCount );


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

	for( int e=0 ; e<boundaryHalfEdges.size() ; e++ )
	{
		SimplexIndex< 1 > eIndex = mesh.surface.edgeIndex( boundaryHalfEdges[e] );
		for( int i=0 ; i<2 ; i++ ) visualization.chartBoundaryVertices.push_back( mesh.surface.vertices[ eIndex[i] ] );
	}

	if( inputMode==SINGLE_INPUT_MODE ) visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'M' , "toggle mask" , ToggleMaskCallBack ) );
	else                               visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'M' , "toggle weights" , ToggleMaskCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 's' , "export texture" , "Output Texture" , ExportTextureCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'y' , "interpolation weight" , "Interpolation Weight" , InterpolationWeightCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , ' ' , "toggle update" , ToggleUpdateCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , '+' , "increment update" , IncrementUpdateCallBack ) );
	
	if( inputMode==MULTIPLE_INPUT_MODE )
	{
		visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 't' , "toggle reference" , ToggleForwardReferenceTextureCallBack ) );
		visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'T' , "toggle reference" , ToggleBackwardReferenceTextureCallBack ) );
	}
	
	visualization.info.push_back( stepsString );
	
	visualization.info.push_back( interpolationStr );
	if( inputMode==MULTIPLE_INPUT_MODE ) visualization.info.push_back( referenceTextureStr );

	visualization.UpdateVertexBuffer();
	visualization.UpdateFaceBuffer();
	visualization.UpdateTextureBuffer( filteredTexture );

	if( inputMode==MULTIPLE_INPUT_MODE ) visualization.UpdateReferenceTextureBuffers( inputTextures );
	if( inputMode==MULTIPLE_INPUT_MODE ) visualization.UpdateReferenceConfidenceBuffers( inputConfidence );
	if( inputMode==SINGLE_INPUT_MODE )   visualization.UpdateCompositeTextureBuffer( inputComposition );
	if( inputMode==SINGLE_INPUT_MODE )   visualization.UpdateMaskTextureBuffer( inputColorMask );

}
#endif // NO_OPEN_GL_VISUALIZATION

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::Init( void )
{
	sprintf( stepsString , "Steps: 0" );
	levels = std::max< int >( Levels.value , 1 );
	interpolationWeight = InterpolationWeight.value;

	mesh.read( In.values[0] , DetailVerbose.set , CollapseEpsilon.value );
	if( InputLowFrequency.set ) ReadImage< TextureBitDepth >( lowFrequencyTexture , InputLowFrequency.value );

	// Define centroid and scale for visualization
	Point3D< PreReal > centroid = mesh.surface.centroid();
	PreReal radius = mesh.surface.boundingRadius( centroid );
	for( int i=0 ; i<mesh.surface.vertices.size() ; i++ ) mesh.surface.vertices[i] = ( mesh.surface.vertices[i] - centroid ) / radius;

	if( RandomJitter.set )
	{
		if( RandomJitter.value ) srand( RandomJitter.value );
		else                     srand( (unsigned int)time(NULL) );
		PreReal jitterScale = (PreReal)1e-3 / std::max< int >( textureWidth , textureHeight );

		for( int i=0 ; i<mesh.texture.vertices.size() ; i++ ) mesh.texture.vertices[i] += Point2D< PreReal >( (PreReal)1. - Random< PreReal >()*2 , (PreReal)1. - Random< PreReal >()*2 ) * jitterScale;
	}

	{
		padding = Padding::Init( textureWidth , textureHeight , mesh.texture.vertices , DetailVerbose.set );
		padding.pad( textureWidth , textureHeight , mesh.texture.vertices );
		if( inputMode==MULTIPLE_INPUT_MODE ) for( int i=0 ; i<numTextures ; i++ ) padding.pad( inputTextures[i] ) , padding.pad( inputConfidence[i] );
		else padding.pad( inputComposition ) , padding.pad( inputMask );

		textureWidth  += padding.width();
		textureHeight += padding.height();
	}

	Miscellany::Timer timer;
	InitializeSystem( textureWidth , textureHeight );
	if( Verbose.set )
	{
		printf( "Resolution: %d / %d x %d\n" , (int)textureNodes.size() , textureWidth , textureHeight );
		printf( "Initialized system %.2f(s)\n" , timer.elapsed() );
		printf( "Peak Memory (MB): %d\n" , Miscellany::MemoryInfo::PeakMemoryUsageMB() );
	}

	// Assign position to exterior nodes using barycentric-exponential map
	{
		FEM::RiemannianMesh< PreReal , unsigned int > rMesh( GetPointer( mesh.surface.triangles) , mesh.surface.triangles.size() );
		rMesh.setMetricFromEmbedding( GetPointer( mesh.surface.vertices ) );
		rMesh.makeUnitArea();
		Pointer( FEM::CoordinateXForm< PreReal > ) xForms = rMesh.getCoordinateXForms();

#ifdef NEW_INDEXING
		for( int i=0 ; i<textureNodes.size() ; i++ ) if( textureNodes[i].tID!=AtlasMeshTriangleIndex(-1) && !textureNodes[i].isInterior )
#else // !NEW_INDEXING
		for( int i=0 ; i<textureNodes.size() ; i++ ) if( textureNodes[i].tID!=-1 && !textureNodes[i].isInterior )
#endif // NEW_INDEXING
		{
			FEM::HermiteSamplePoint< PreReal > _p;
#ifdef NEW_INDEXING
			_p.tIdx = static_cast< unsigned int >( textureNodes[i].tID );
#else // !NEW_INDEXING
			_p.tIdx = textureNodes[i].tID;
#endif // NEW_INDEXING
			_p.p = Point2D< PreReal >( (PreReal)1./3 , (PreReal)1./3 );
			_p.v = textureNodes[i].barycentricCoords - _p.p;

			rMesh.exp( xForms , _p );

#ifdef NEW_INDEXING
			textureNodes[i].tID = AtlasMeshTriangleIndex( _p.tIdx );
#else // !NEW_INDEXING
			textureNodes[i].tID = _p.tIdx;
#endif // NEW_INDEXING
			textureNodes[i].barycentricCoords = _p.p;
		}
	}

	textureNodePositions.resize( textureNodes.size() );
	for( int i=0 ; i<textureNodePositions.size() ; i++ ) textureNodePositions[i] = mesh.surface( textureNodes[i] );

	textureEdgePositions.resize( edgePairs.size() );
	for( int i=0 ; i<edgePairs.size() ; i++ ) textureEdgePositions[i] = ( textureNodePositions[ edgePairs[i][0] ] + textureNodePositions[ edgePairs[i][1] ] ) / 2;

	{
		unsigned int multiChartTexelCount = 0;
		Image< int > texelId;
		texelId.resize( textureWidth , textureHeight );
		for( int i=0 ; i<texelId.size() ; i++ ) texelId[i] = -1;
		for( int i=0 ; i<textureNodes.size() ; i++ )
		{
			int ci = textureNodes[i].ci , cj = textureNodes[i].cj;
			if( texelId(ci,cj)!=-1 ) multiChartTexelCount++;
			texelId(ci,cj) = i;
		}
		if( multiChartTexelCount ) MK_WARN( "Non-zero multi-chart texels: " , multiChartTexelCount );
	}
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void _main( int argc , char *argv[] )
{
	Stitching< PreReal , Real , TextureBitDepth >::inputMode = MultiInput.set ? MULTIPLE_INPUT_MODE : SINGLE_INPUT_MODE;
	Stitching< PreReal , Real , TextureBitDepth >::updateCount = 0;

	Stitching< PreReal , Real , TextureBitDepth >::LoadTextures();
	Stitching< PreReal , Real , TextureBitDepth >::Init();
	Stitching< PreReal , Real , TextureBitDepth >::LoadMasks();
	Stitching< PreReal , Real , TextureBitDepth >::ParseImages();
	Stitching< PreReal , Real , TextureBitDepth >::SetUpSystem();

#ifdef NO_OPEN_GL_VISUALIZATION
	Stitching< PreReal , Real , TextureBitDepth >::SolveSystem();
	Stitching< PreReal , Real , TextureBitDepth >::WriteTexture( Output.value.c_str() );
#else // !NO_OPEN_GL_VISUALIZATION
	if( !Output.set )
	{
		glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
		Stitching< PreReal , Real , TextureBitDepth >::visualization.visualizationMode = Stitching< PreReal , Real , TextureBitDepth >::inputMode;
		Stitching< PreReal , Real , TextureBitDepth >::visualization.displayMode = TWO_REGION_DISPLAY;
		Stitching< PreReal , Real , TextureBitDepth >::visualization.screenWidth = 1600;
		Stitching< PreReal , Real , TextureBitDepth >::visualization.screenHeight = 800;
		Stitching< PreReal , Real , TextureBitDepth >::visualization.useNearestSampling = Nearest.set;

		glutInitWindowSize( Stitching< PreReal , Real , TextureBitDepth >::visualization.screenWidth , Stitching< PreReal , Real , TextureBitDepth >::visualization.screenHeight );
		glutInit( &argc , argv );
		char windowName[1024];
		sprintf( windowName , "Stitching" );
		glutCreateWindow( windowName );
		if( glewInit()!=GLEW_OK ) MK_THROW( "glewInit failed" );
		glutDisplayFunc ( Stitching< PreReal , Real , TextureBitDepth >::Display );
		glutReshapeFunc ( Stitching< PreReal , Real , TextureBitDepth >::Reshape );
		glutMouseFunc   ( Stitching< PreReal , Real , TextureBitDepth >::MouseFunc );
		glutMotionFunc  ( Stitching< PreReal , Real , TextureBitDepth >::MotionFunc );
		glutKeyboardFunc( Stitching< PreReal , Real , TextureBitDepth >::KeyboardFunc );
		glutIdleFunc    ( Stitching< PreReal , Real , TextureBitDepth >::Idle );
		if( CameraConfig.set ) Stitching< PreReal , Real , TextureBitDepth >::visualization.ReadSceneConfigurationCallBack( &Stitching< PreReal , Real , TextureBitDepth >::visualization , CameraConfig.value.c_str() );
		Stitching< PreReal , Real , TextureBitDepth >::InitializeVisualization();
		glutMainLoop();
	}
	else
	{
		Stitching< PreReal , Real , TextureBitDepth >::SolveSystem();
		Stitching< PreReal , Real , TextureBitDepth >::ExportTextureCallBack( &Stitching< PreReal , Real , TextureBitDepth >::visualization , Output.value.c_str() );
	}
#endif // NO_OPEN_GL_VISUALIZATION
}

template< typename PreReal , typename Real >
void _main( int argc , char *argv[] , unsigned int bitDepth )
{
	switch( bitDepth )
	{
	case  8: return _main< PreReal , Real ,  8 >( argc , argv );
	case 16: return _main< PreReal , Real , 16 >( argc , argv );
	case 32: return _main< PreReal , Real , 32 >( argc , argv );
	case 64: return _main< PreReal , Real , 64 >( argc , argv );
	default: MK_THROW( "Only bit depths of 8, 16, 32, and 64 supported: " , bitDepth );
	}
}


int main( int argc , char* argv[] )
{
	CmdLineParse( argc-1 , argv+1 , params );
#ifdef NO_OPEN_GL_VISUALIZATION
	if( !In.set || !Output.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
#else // !NO_OPEN_GL_VISUALIZATION
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
#endif // NO_OPEN_GL_VISUALIZATION
	if( MultiInput.set && !InMask.set ) MK_THROW( "Input mask required for multi-input" );

	unsigned int bitDepth;
	{
		unsigned int width , height , channels;
		if( MultiInput.set )
		{
			unsigned int numTextures = 0;
			while( true )
			{
				char textureName[256];
				sprintf( textureName , In.values[1].c_str() , numTextures );
				FILE * file = fopen( textureName , "r" );
				if( file )
				{
					fclose( file );
					unsigned int _width , _height , _channels , _bitDepth;
					ImageReader< 8 >::GetInfo( textureName , _width , _height , _channels , _bitDepth );
					if( !numTextures ) width = _width , height = _height , channels = _channels , bitDepth = _bitDepth;
					else if( width!=_width || height!=_height || channels!=_channels || bitDepth!=_bitDepth )
						MK_THROW( "Image properties don't match: (" , width , " " , height , " " , channels , " " , bitDepth , ") != (" , _width , " " , _height , " " , _channels , " " , _bitDepth , ")" );
					numTextures++;
				}
				else break;
			}
		}
		else ImageReader< 8 >::GetInfo( In.values[1].c_str() , width , height , channels , bitDepth );
	}

	if( Serial.set ) ThreadPool::ParallelizationType = ThreadPool::ParallelType::NONE;
	if( !NoHelp.set && !Output.set )
	{
		printf( "+----------------------------------------------------------------------------+\n" );
		printf( "| Interface Controls:                                                        |\n" );
		printf( "|    [Left Mouse]:                rotate                                     |\n" );
		printf( "|    [Right Mouse]:               zoom                                       |\n" );
		printf( "|    [Left/Right Mouse] + [CTRL]: pan                                        |\n" );
		if( MultiInput.set )
		{
			printf( "|    [Left Mouse] + [SHIFT]:      mark region to in-paint                    |\n" );
			printf( "|    't':                         toggle textures forward                    |\n" );
			printf( "|    'T':                         toggle textures backward                   |\n" );
		}
		else printf( "|    [SPACE]:                     start solver                               |\n" );
		printf( "|    'y':                         prescribe interpolation weight             |\n" );
		printf( "+----------------------------------------------------------------------------+\n" );
	}
	try
	{
		if( Double.set ) _main< double , double >( argc , argv , bitDepth );
		else             _main< double , float  >( argc , argv , bitDepth );
	}
	catch( Exception &e )
	{
		printf( "%s\n" , e.what() );
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
