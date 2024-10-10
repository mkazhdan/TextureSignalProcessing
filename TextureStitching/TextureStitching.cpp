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


#include <Misha/CmdLineParser.h> 
#include <Misha/Miscellany.h>
#include <Misha/FEM.h>
#include <Src/Hierarchy.h>
#include <Src/SimpleMesh.h>
#include <Src/Basis.h>
#include <Src/Solver.h>
#include <Src/MassAndStiffness.h>
#include <Src/Padding.h>
#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
#include <Src/StitchingVisualization.h>
#endif // NO_VISUALIZATION

cmdLineParameterArray< char * , 3 > In( "in" );
cmdLineParameter< char* > Output( "out" );
cmdLineParameter< int   > OutputVCycles( "outVCycles" , 6 );
cmdLineParameter< float > InterpolationWeight( "interpolation" , 1e2 );
cmdLineParameter< int   > Threads( "threads", omp_get_num_procs() );
cmdLineParameter< int   > Levels( "levels" , 4 );
cmdLineParameter< int   > MatrixQuadrature( "mQuadrature" , 6 );

cmdLineParameter< int   > MultigridBlockHeight ( "mBlockH" ,  16 );
cmdLineParameter< int   > MultigridBlockWidth  ( "mBlockW" , 128 );
cmdLineParameter< int   > MultigridPaddedHeight( "mPadH"   ,   0 );
cmdLineParameter< int   > MultigridPaddedWidth ( "mPadW"   ,   2 );

cmdLineReadable RandomJitter( "jitter" );
cmdLineParameter< int    > OutputChartMaskErode( "outChartMaskErode" , 0 );
cmdLineParameter< char * > OutputChartMask( "outChartMask" );
#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
cmdLineParameter< char* > CameraConfig( "camera" );
#endif // NO_VISUALIZATION
cmdLineReadable UseDirectSolver( "useDirectSolver" );
cmdLineReadable Verbose( "verbose" );
cmdLineReadable NoHelp( "noHelp" );
cmdLineReadable DetailVerbose( "detail" );
cmdLineReadable Double( "double" );
cmdLineReadable MultiInput( "multi" );
cmdLineReadable* params[] =
{
	&In , &Output , &InterpolationWeight , &Levels , &UseDirectSolver , &Threads, &Verbose ,
	&DetailVerbose , &MultigridBlockHeight , &MultigridBlockWidth , &MultigridPaddedHeight , &MultigridPaddedWidth , &RandomJitter ,
	&Double , &MatrixQuadrature , &OutputVCycles , &NoHelp , &MultiInput ,
#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
	&CameraConfig ,
#endif // NO_VISUALIZATION
	&OutputChartMaskErode ,
	&OutputChartMask ,
	NULL
};

void ShowUsage(const char* ex)
{
	printf("Usage %s:\n", ex);

	printf( "\t --%s <input mesh, texels, and mask>\n" , In.name );
#ifdef NO_VISUALIZATION
	printf( "\t --%s <output texture>\n" , Output.name );
#else // !NO_VISUALIZATION
	printf( "\t[--%s <output texture>]\n" , Output.name );
#endif // NO_VISUALIZATION
	printf( "\t[--%s <output chart mask image>]\n" , OutputChartMask.name );
	printf( "\t[--%s <output chart mask erosion radius>=%d]\n" , OutputChartMaskErode.name );
	printf( "\t[--%s <output v-cycles>=%d]\n" , OutputVCycles.name , OutputVCycles.value );
	printf( "\t[--%s <interpolation weight>=%f]\n" , InterpolationWeight.name , InterpolationWeight.value );
	printf( "\t[--%s <system matrix quadrature points per triangle>=%d]\n" , MatrixQuadrature.name, MatrixQuadrature.value );
	printf( "\t[--%s]\n" , UseDirectSolver.name );
	printf( "\t[--%s]\n" , MultiInput.name );
	printf( "\t[--%s]\n" , RandomJitter.name );
	printf( "\t[--%s]\n" , Verbose.name );

#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
	printf( "\t[--%s <camera configuration file>]\n" , CameraConfig.name );
#endif // NO_VISUALIZATION
	printf( "\t[--%s <hierarchy levels>=%d]\n" , Levels.name, Levels.value );
	printf( "\t[--%s <threads>=%d]\n", Threads.name , Threads.value );
	printf( "\t[--%s <multigrid block width>=%d]\n" , MultigridBlockWidth.name , MultigridBlockWidth.value );
	printf( "\t[--%s <multigrid block height>=%d]\n" , MultigridBlockHeight.name , MultigridBlockHeight.value );
	printf( "\t[--%s <multigrid padded width>=%d]\n" , MultigridPaddedWidth.name , MultigridPaddedWidth.value );
	printf( "\t[--%s <multigrid padded height>=%d]\n" , MultigridPaddedHeight.name , MultigridPaddedHeight.value );
	printf( "\t[--%s]\n", DetailVerbose.name );
	printf( "\t[--%s]\n" , NoHelp.name );
}

template< typename PreReal , typename Real >
class Stitching
{
public:
	static int inputMode;
	static TexturedMesh< PreReal > mesh;
	static int textureWidth;
	static int textureHeight;
	static Real interpolationWeight;
	static int levels;

	static HierarchicalSystem< PreReal , Real > hierarchy;
	static bool rhsUpdated;
	static bool positiveModulation;

	// Single input mode
	static Image< int > inputMask;
	static Image< Point3D< Real > > inputComposition;
	static Image< Point3D< Real > > inputColorMask;

	// Multiple input mode
	static int numTextures;
	static std::vector< Image< Real > > inputConfidence;
	static std::vector< Image< Point3D< Real > > > inputTextures;
	static std::vector< std::vector< Point3D< Real > > > partialTexelValues;
	static std::vector< std::vector< Point3D< Real > > > partialEdgeValues;

	static Image< Point3D< Real > > filteredTexture;

#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
	// UI
	static char interpolationStr[1024];
	static char referenceTextureStr[1024];
#endif // NO_VISUALIZATION

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
#elif defined( USE_EIGEN_SIMPLICIAL )
	typedef EigenCholeskySolver< Real , 3 > DirectSolver;
#elif defined( USE_EIGEN_PARDISO )
	typedef EigenPardisoSolver< Real , 3 > DirectSolver;
#else
#error "[ERROR] No solver defined!"
#endif

	static VCycleSolvers< DirectSolver > vCycleSolvers;
	static DirectSolver directSolver;

	static std::unordered_map< unsigned long long , int > edgeIndex;
	static std::vector< std::pair< int , int > > edgePairs;

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

#ifdef NO_VISUALIZATION
	static int updateCount;
	static void WriteTexture( const char *fileName );
#else // !NO_VISUALIZATION
	// Visulization
	static StitchingVisualization visualization;
	static int updateCount;

	static void ToggleForwardReferenceTextureCallBack ( Visualization *v , const char *prompt );
	static void ToggleBackwardReferenceTextureCallBack( Visualization *v , const char *prompt );
	static void ToggleMaskCallBack                    ( Visualization *v , const char *prompt );
	static void ToggleUpdateCallBack                  ( Visualization *v , const char *prompt );
	static void IncrementUpdateCallBack               ( Visualization *v , const char *prompt );
	static void ExportTextureCallBack                 ( Visualization *v , const char *prompt );
	static void InterpolationWeightCallBack           ( Visualization *v , const char *prompt );
#endif // NO_VISUALIZATION

	static void LoadImages( void );
	static void ParseImages( void );
	static void SetUpSystem( void );
	static void SolveSystem( void );
	static void Init( void );
#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
	static void InitializeVisualization( void );
#endif // NO_VISUALIZATION
	static void UpdateSolution( bool verbose=false , bool detailVerbose=false );
	static void ComputeExactSolution( bool verbose=false );
	static void InitializeSystem( int width , int height );
	static void _InitializeSystem( std::vector< std::vector< SquareMatrix< PreReal , 2 > > > &parameterMetric , BoundaryProlongationData< Real > &boundaryProlongation );

#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
	static void UpdateFilteredColorTexture( const std::vector< Point3D< Real > > &solution );
#endif // NO_VISUALIZATION
	static void UpdateFilteredTexture( const std::vector< Point3D< Real > > &solution );

#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
	static void Display( void ){ visualization.Display(); }
	static void MouseFunc( int button , int state , int x , int y );
	static void MotionFunc( int x , int y );
	static void Reshape( int w , int h ) { visualization.Reshape(w,h); }
	static void KeyboardFunc( unsigned char key , int x , int y ) { visualization.KeyboardFunc(key,x,y); }
	static void Idle( void );
#endif // NO_VISUALIZATION
};

#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
template< typename PreReal , typename Real > char															Stitching< PreReal , Real >::referenceTextureStr[1024];
template< typename PreReal , typename Real > char															Stitching< PreReal , Real >::interpolationStr[1024];
#endif // NO_VISUALIZATION

template< typename PreReal , typename Real > int															Stitching< PreReal , Real >::inputMode;
template< typename PreReal , typename Real > TexturedMesh< PreReal >										Stitching< PreReal , Real >::mesh;
template< typename PreReal , typename Real > int															Stitching< PreReal , Real >::textureWidth;
template< typename PreReal , typename Real > int															Stitching< PreReal , Real >::textureHeight;
#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
template< typename PreReal , typename Real > StitchingVisualization											Stitching< PreReal , Real >::visualization;
#endif // NO_VISUALIZATION
template< typename PreReal , typename Real > SparseMatrix< Real , int >										Stitching< PreReal , Real >::mass;
template< typename PreReal , typename Real > SparseMatrix< Real , int >										Stitching< PreReal , Real >::stiffness;
template< typename PreReal , typename Real > SparseMatrix< Real , int >										Stitching< PreReal , Real >::stitchingMatrix;

template< typename PreReal , typename Real > Real															Stitching< PreReal , Real >::interpolationWeight;

template< typename PreReal , typename Real > std::vector< TextureNodeInfo< PreReal > >						Stitching< PreReal , Real >::textureNodes;
template< typename PreReal , typename Real > std::vector< BilinearElementIndex >							Stitching< PreReal , Real >::bilinearElementIndices;

template< typename PreReal , typename Real > int															Stitching< PreReal , Real >::steps;
template< typename PreReal , typename Real > char															Stitching< PreReal , Real >::stepsString[1024];
template< typename PreReal , typename Real > int															Stitching< PreReal , Real >::levels;
template< typename PreReal , typename Real > HierarchicalSystem< PreReal , Real >							Stitching< PreReal , Real >::hierarchy;

template< typename PreReal , typename Real > bool															Stitching< PreReal , Real >::rhsUpdated = true;
template< typename PreReal , typename Real > bool															Stitching< PreReal , Real >::positiveModulation = true;

template< typename PreReal , typename Real> Image< Point3D< Real > >										Stitching< PreReal , Real >::filteredTexture;

template< typename PreReal , typename Real > Image< int >												    Stitching< PreReal , Real >::inputMask;
template< typename PreReal , typename Real > Image< Point3D< Real > >										Stitching< PreReal , Real >::inputComposition;
template< typename PreReal , typename Real > Image< Point3D< Real > >										Stitching< PreReal , Real >::inputColorMask;

template< typename PreReal , typename Real > int															Stitching< PreReal , Real >::numTextures;
template< typename PreReal , typename Real > std::vector< Image< Real > >									Stitching< PreReal , Real >::inputConfidence;
template< typename PreReal , typename Real > std::vector< Image< Point3D< Real > > >						Stitching< PreReal , Real >::inputTextures;
template< typename PreReal , typename Real > std::vector< Point3D< Real > >									Stitching< PreReal , Real >::texelMass;
template< typename PreReal , typename Real > std::vector< Point3D< Real > >									Stitching< PreReal , Real >::texelDivergence;

template< typename PreReal , typename Real > std::vector< SystemCoefficients< Real > >						Stitching< PreReal , Real >::multigridStitchingCoefficients;
template< typename PreReal , typename Real > std::vector< MultigridLevelVariables< Point3D< Real > > >		Stitching< PreReal , Real >::multigridStitchingVariables;
template< typename PreReal , typename Real > std::vector< MultigridLevelIndices< Real > >					Stitching< PreReal , Real >::multigridIndices;

template< typename PreReal , typename Real > VCycleSolvers< typename Stitching< PreReal , Real >::DirectSolver >		Stitching< PreReal , Real >::vCycleSolvers;
template< typename PreReal , typename Real > typename Stitching< PreReal , Real >::DirectSolver				Stitching< PreReal , Real >::directSolver;

template< typename PreReal , typename Real > std::vector< AtlasChart< PreReal > >							Stitching< PreReal , Real >::atlasCharts;

template< typename PreReal , typename Real > std::vector< Point3D< float > >								Stitching< PreReal , Real >::textureNodePositions;
template< typename PreReal , typename Real > std::vector< Point3D< float > >								Stitching< PreReal , Real >::textureEdgePositions;

template< typename PreReal , typename Real > Padding														Stitching< PreReal , Real >::padding;

template< typename PreReal , typename Real > SystemCoefficients< Real >										Stitching< PreReal , Real >::massCoefficients;
template< typename PreReal , typename Real > SystemCoefficients< Real >										Stitching< PreReal , Real >::stiffnessCoefficients;

template< typename PreReal , typename Real > int															Stitching< PreReal , Real >::updateCount = -1;

template< typename PreReal , typename Real >  std::unordered_map< unsigned long long , int >				Stitching< PreReal , Real >::edgeIndex;
template< typename PreReal , typename Real >  std::vector< std::pair< int , int > >							Stitching< PreReal , Real >::edgePairs;

template< typename PreReal , typename Real >  SparseMatrix< Real , int >									Stitching< PreReal , Real >::boundaryDivergenceMatrix;
template< typename PreReal , typename Real >  std::vector< Real >											Stitching< PreReal , Real >::deepDivergenceCoefficients;
template< typename PreReal , typename Real >  std::vector< DivegenceRasterLine >  							Stitching< PreReal , Real >::divergenceRasterLines;

template< typename PreReal , typename Real >  std::vector< bool >											Stitching< PreReal , Real >::unobservedTexel;
template< typename PreReal , typename Real >  std::vector< Point3D< Real > >								Stitching< PreReal , Real >::texelValues;
template< typename PreReal , typename Real >  std::vector< Point3D< Real > >								Stitching< PreReal , Real >::edgeValues;
template< typename PreReal , typename Real >  std::vector< std::vector< Point3D< Real > > >					Stitching< PreReal , Real >::partialTexelValues;
template< typename PreReal , typename Real >  std::vector< std::vector< Point3D< Real > > >					Stitching< PreReal , Real >::partialEdgeValues;
template< typename PreReal , typename Real > int															Stitching< PreReal , Real >::textureIndex = 0;


#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
template< typename PreReal , typename Real >
void Stitching< PreReal , Real >::UpdateFilteredColorTexture( const std::vector< Point3D< Real > > &solution )
{
#pragma omp parallel for
	for( int i=0 ; i<textureNodes.size() ; i++ )
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
}
#endif // NO_VISUALIZATION

template< typename PreReal , typename Real >
void Stitching< PreReal , Real >::UpdateFilteredTexture( const std::vector< Point3D< Real > > &solution )
{
#pragma omp parallel for
	for( int i=0 ; i<textureNodes.size() ; i++ )
	{
		int ci = textureNodes[i].ci , cj = textureNodes[i].cj;
		filteredTexture(ci,cj) = solution[i];
	}
}

#ifdef NO_VISUALIZATION
template< typename PreReal , typename Real >
void Stitching< PreReal , Real >::WriteTexture( const char *fileName )
{
	UpdateFilteredTexture( multigridStitchingVariables[0].x );
	Image< Point3D< Real > > outputTexture = filteredTexture;
	if( padding.nonTrivial ) UnpadImage( padding , outputTexture );
	outputTexture.write( fileName );
}

#else // !NO_VISUALIZATION
template< typename PreReal , typename Real >
void Stitching< PreReal , Real >::Idle( void )
{
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
#pragma omp parallel for
				for( int i=0 ; i<textureNodePositions.size() ; i++ ) if( Point3D< float >::SquareNorm( textureNodePositions[i]-selectedPoint )<radiusSquared ) texelValues[i] = partialTexelValues[textureIndex][i];
#pragma omp parallel for
				for( int i=0 ; i<textureEdgePositions.size() ; i++ ) if( Point3D< float >::SquareNorm( textureEdgePositions[i]-selectedPoint )<radiusSquared ) edgeValues[i] = partialEdgeValues[textureIndex][i];
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

template< typename PreReal , typename Real >
void Stitching< PreReal , Real >::MouseFunc( int button , int state , int x , int y )
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
	}

	glutPostRedisplay();
}

template< typename PreReal , typename Real >
void Stitching< PreReal , Real >::MotionFunc( int x , int y )
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
			if( visualization.panning ) visualization.xForm.offset[0] -= ( visualization.newX - visualization.oldX ) / visualization.imageToScreenScale() , visualization.xForm.offset[1] += ( visualization.newY - visualization.oldY ) / visualization.imageToScreenScale();
			else
			{
				float dz = (float) pow( 1.1 , (double)( visualization.newY - visualization.oldY)/8 );
				visualization.xForm.zoom *= dz;
			}
		}
	}
	glutPostRedisplay();
}

template< typename PreReal , typename Real > void Stitching< PreReal , Real >::ToggleMaskCallBack( Visualization * /*v*/ , const char * )
{
	visualization.showMask = !visualization.showMask;
	glutPostRedisplay();
}

template< typename PreReal , typename Real > void Stitching< PreReal , Real >::ToggleForwardReferenceTextureCallBack( Visualization * /*v*/ , const char* )
{
	textureIndex = ( textureIndex+1 ) % numTextures;
	visualization.referenceIndex = textureIndex;
	sprintf( referenceTextureStr , "Reference Texture: %02d of %02d \n" , textureIndex , numTextures );
	glutPostRedisplay();
}

template< typename PreReal , typename Real > void Stitching< PreReal , Real >::ToggleBackwardReferenceTextureCallBack( Visualization * /*v*/ , const char* )
{
	textureIndex = ( textureIndex+numTextures-1 ) % numTextures;
	visualization.referenceIndex = textureIndex;
	sprintf( referenceTextureStr, "Reference Texture: %02d of %02d \n" , textureIndex , numTextures );
	glutPostRedisplay();
}

template< typename PreReal , typename Real > void Stitching< PreReal , Real >::ToggleUpdateCallBack( Visualization * /*v*/ , const char * )
{
	if( updateCount ) updateCount = 0;
	else              updateCount = -1;
}

template< typename PreReal , typename Real >
void Stitching< PreReal , Real >::IncrementUpdateCallBack( Visualization * /*v*/ , const char * )
{
	if( updateCount<0 ) updateCount = 1;
	else                updateCount++;
}

template< typename PreReal , typename Real >
void Stitching< PreReal , Real >::ExportTextureCallBack( Visualization * /*v*/ , const char *prompt )
{
	UpdateFilteredTexture( multigridStitchingVariables[0].x );
	Image< Point3D< Real > > outputTexture = filteredTexture;
	if( padding.nonTrivial ) UnpadImage( padding , outputTexture );
	outputTexture.write( prompt );
}

template< typename PreReal , typename Real >
void  Stitching< PreReal , Real >::InterpolationWeightCallBack( Visualization * /*v*/ , const char *prompt )
{
	interpolationWeight = atof(prompt);
	if( UseDirectSolver.set ) stitchingMatrix = mass*interpolationWeight + stiffness;
	Miscellany::Timer timer;
	UpdateLinearSystem( interpolationWeight , (Real)1. , hierarchy , multigridStitchingCoefficients , massCoefficients , stiffnessCoefficients , vCycleSolvers , directSolver , stitchingMatrix , DetailVerbose.set , false , UseDirectSolver.set );
	if( Verbose.set ) printf( "\tInitialized multigrid coefficients: %.2f(s)\n" , timer.elapsed() );

#pragma omp parallel for
	for( int i=0 ; i<multigridStitchingVariables[0].rhs.size() ; i++ ) multigridStitchingVariables[0].rhs[i] = texelMass[i] * interpolationWeight + texelDivergence[i];

	if( UseDirectSolver.set ) ComputeExactSolution(Verbose.set);
	else for( int i=0 ; i<OutputVCycles.value ; i++ ) VCycle( multigridStitchingVariables , multigridStitchingCoefficients , multigridIndices , vCycleSolvers , false , false );

	UpdateFilteredColorTexture( multigridStitchingVariables[0].x );
	visualization.UpdateColorTextureBuffer();
	sprintf( interpolationStr , "Interpolation weight: %e\n" , interpolationWeight );
}
#endif // NO_VISUALIZATION


template< typename PreReal , typename Real >
void Stitching< PreReal , Real >::ComputeExactSolution( bool verbose )
{
	Miscellany::Timer timer;
	solve( directSolver , multigridStitchingVariables[0].x , multigridStitchingVariables[0].rhs );
	if( verbose ) printf( "Solving time =  %.4f\n" , timer.elapsed() );
}

template< typename PreReal , typename Real >
void Stitching< PreReal , Real >::UpdateSolution( bool verbose , bool detailVerbose )
{
	if( !rhsUpdated )
	{
		Miscellany::Timer timer;

		MultiplyBySystemMatrix_NoReciprocals( massCoefficients , hierarchy.gridAtlases[0].boundaryGlobalIndex , hierarchy.gridAtlases[0].rasterLines , texelValues , texelMass );
		ComputeDivergence( edgeValues , texelDivergence , deepDivergenceCoefficients , boundaryDivergenceMatrix , divergenceRasterLines );

#pragma omp parallel for
		for( int i=0 ; i<textureNodes.size() ; i++ ) multigridStitchingVariables[0].rhs[i] = texelMass[i] * interpolationWeight + texelDivergence[i];

		if( verbose ) printf( "RHS update time %.4f\n" , timer.elapsed() );
		rhsUpdated = true;
	}

	VCycle( multigridStitchingVariables , multigridStitchingCoefficients , multigridIndices , vCycleSolvers , verbose , detailVerbose );
}

template< typename PreReal , typename Real >
void Stitching< PreReal , Real >::_InitializeSystem( std::vector< std::vector< SquareMatrix< PreReal , 2 > > > &parameterMetric , BoundaryProlongationData< Real > &boundaryProlongation )
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
		default: Miscellany::Throw( "Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles" );
		}
	}
}

template< typename PreReal , typename Real >
void Stitching< PreReal , Real >::InitializeSystem( int width , int height )
{
	Miscellany::Timer timer;
	MultigridBlockInfo multigridBlockInfo( MultigridBlockWidth.value , MultigridBlockHeight.value , MultigridPaddedWidth.value , MultigridPaddedHeight.value , 0 );
	InitializeHierarchy( mesh , width , height , levels , textureNodes , bilinearElementIndices , hierarchy , atlasCharts , multigridBlockInfo , true , DetailVerbose.set );
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
	for( auto edgeIter=edgeIndex.begin() ; edgeIter!=edgeIndex.end() ; edgeIter++ )
	{
		unsigned long long edgeKey = (*edgeIter).first;
		int edgeId = (*edgeIter).second;

		unsigned long edgeCorners[2];
		GetMeshEdgeIndices( edgeKey , edgeCorners[0] , edgeCorners[1] );
		edgePairs[edgeId] = std::pair< int , int >( edgeCorners[0] , edgeCorners[1] );
	}

	texelMass.resize( textureNodes.size() );
	texelDivergence.resize( textureNodes.size() );

	if( UseDirectSolver.set )
	{
		FullMatrixConstruction( hierarchy.gridAtlases[0] , massCoefficients , mass );
		FullMatrixConstruction( hierarchy.gridAtlases[0] , stiffnessCoefficients , stiffness );
		stitchingMatrix = mass*interpolationWeight + stiffness;
	}

	multigridIndices.resize(levels);
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

	if( Verbose.set ) timer.reset();
	UpdateLinearSystem( interpolationWeight , (Real)1. , hierarchy , multigridStitchingCoefficients , massCoefficients , stiffnessCoefficients , vCycleSolvers , directSolver , stitchingMatrix , DetailVerbose.set, true, UseDirectSolver.set );
	if( Verbose.set ) printf( "\tInitialized multigrid coefficients: %.2f(s)\n" , timer.elapsed() );

	multigridStitchingVariables.resize(levels);
	for( int i=0 ; i<levels ; i++ )
	{
		MultigridLevelVariables< Point3D< Real > >& variables = multigridStitchingVariables[i];
		variables.x.resize( hierarchy.gridAtlases[i].numTexels );
		variables.rhs.resize( hierarchy.gridAtlases[i].numTexels );
		variables.residual.resize( hierarchy.gridAtlases[i].numTexels );
		variables.boundary_rhs.resize( hierarchy.gridAtlases[i].boundaryGlobalIndex.size() );
		variables.boundary_value.resize( hierarchy.gridAtlases[i].boundaryGlobalIndex.size() );
		variables.variable_boundary_value.resize( hierarchy.gridAtlases[i].boundaryGlobalIndex.size() );
	}
}
template< typename PreReal , typename Real >
void Stitching< PreReal , Real >::SetUpSystem( void )
{
	texelMass.resize( textureNodes.size() );
	MultiplyBySystemMatrix_NoReciprocals( massCoefficients , hierarchy.gridAtlases[0].boundaryGlobalIndex , hierarchy.gridAtlases[0].rasterLines , texelValues , texelMass );

	ComputeDivergence( edgeValues , texelDivergence , deepDivergenceCoefficients , boundaryDivergenceMatrix , divergenceRasterLines );

#pragma omp parallel for
	for( int i=0 ; i<texelValues.size() ; i++ )
	{
		multigridStitchingVariables[0].x[i] = texelValues[i];
		multigridStitchingVariables[0].rhs[i] = texelMass[i] * interpolationWeight + texelDivergence[i];
	}

	filteredTexture.resize( textureWidth , textureHeight );
	for( int i=0 ; i<filteredTexture.size() ; i++ ) filteredTexture[i] = Point3D< Real >( (Real)0.5 , (Real)0.5 , (Real)0.5 );

	UpdateFilteredTexture( multigridStitchingVariables[0].x );
}

template< typename PreReal , typename Real >
void Stitching< PreReal , Real >::SolveSystem( void )
{
	if( UseDirectSolver.set ) ComputeExactSolution( Verbose.set );
	else for( int i=0 ; i<OutputVCycles.value ; i++ ) VCycle( multigridStitchingVariables , multigridStitchingCoefficients , multigridIndices , vCycleSolvers , false , false );	
	UpdateFilteredTexture( multigridStitchingVariables[0].x );
}

template< typename PreReal , typename Real >
void Stitching< PreReal , Real >::LoadImages( void )
{
	if( inputMode==MULTIPLE_INPUT_MODE )
	{
		bool countingTextures = true;
		numTextures = 0;
		while( countingTextures )
		{
			char textureName[256];
			sprintf( textureName , In.values[1] , numTextures );
			FILE * file = fopen( textureName , "r" );
			if( file ) numTextures++;
			else countingTextures = false;
		}

		inputTextures.resize( numTextures );
#pragma omp parallel for
		for( int i=0 ; i<numTextures ; i++ )
		{
			char textureName[256];
			sprintf( textureName , In.values[1] , i );
			inputTextures[i].read( textureName );
		}

		textureWidth = inputTextures[0].width();
		textureHeight = inputTextures[0].height();

		if( Verbose.set ) printf( "Texture count: %d\n" , numTextures );

		inputConfidence.resize( numTextures );
#pragma omp parallel for
		for( int i=0 ; i<numTextures ; i++ )
		{
			char confidenceName[256];
			sprintf( confidenceName , In.values[2] , i );
			Image< Point3D< Real > > textureConfidence;
			textureConfidence.read( confidenceName );
			inputConfidence[i].resize( textureWidth , textureHeight );
			for( int p=0 ; p<textureConfidence.size() ; p++ ) inputConfidence[i][p] = Point3D< Real >::Dot( textureConfidence[p] , Point3D< Real >( (Real)1./3 , (Real)1./3 , (Real)1./3 ) );
		}
	}
	else
	{
		inputComposition.read( In.values[1] );
		textureWidth = inputComposition.width();
		textureHeight = inputComposition.height();
		inputColorMask.read( In.values[2] );

		Image< Point3D< unsigned char > > textureConfidence;
		textureConfidence.read( In.values[2] );

		inputMask.resize( textureWidth , textureHeight );
		for( int p=0 ; p<textureConfidence.size() ; p++ )
		{
			int index = ( (int)textureConfidence[p][0] ) * 256 * 256 + ( (int)textureConfidence[p][1] ) * 256 + ( (int)textureConfidence[p][2] );
			inputMask[p] = index==0 ? -1 : index;
		}
	}
}

template< typename PreReal , typename Real >
void Stitching< PreReal , Real >::ParseImages( void )
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
			Image< Point3D< Real > > textureValues = inputTextures[textureIter];
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
				int edgeCorners[2] = { edgePairs[e].first , edgePairs[e].second };

				Real weight[2];
				Point3D< Real > value[2];
				for( int k=0 ; k<2 ; k++ )
				{
					int i = edgeCorners[k];
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
			texelValues[i] = inputComposition( textureNodes[i].ci , textureNodes[i].cj );
		}
		for( int e=0 ; e<edgePairs.size() ; e++ )
		{
			int edgeCorners[2] = { edgePairs[e].first , edgePairs[e].second };
			int ci[] = { textureNodes[ edgeCorners[0] ].ci , textureNodes[ edgeCorners[1] ].ci };
			int cj[] = { textureNodes[ edgeCorners[0] ].cj , textureNodes[ edgeCorners[1] ].cj };
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

#pragma omp parallel for
	for( int i=0 ; i<textureNodes.size() ; i++ ) multigridStitchingVariables[0].x[i] = texelValues[i];
}

#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
template< typename PreReal , typename Real >
void Stitching< PreReal , Real >::InitializeVisualization( void )
{
	sprintf( referenceTextureStr , "Reference Texture: %02d of %02d\n" , textureIndex,numTextures );
	sprintf( interpolationStr , "Interpolation: %.2e\n" , interpolationWeight );

	visualization.textureWidth = textureWidth;
	visualization.textureHeight = textureHeight;

	visualization.colorTextureBuffer = new unsigned char[textureHeight*textureWidth * 3];
	memset( visualization.colorTextureBuffer , 204 , textureHeight * textureWidth * 3 * sizeof(unsigned char) );



	int tCount = (int)mesh.triangles.size();

	visualization.triangles.resize( tCount );
	visualization.vertices.resize( 3*tCount );
	visualization.colors.resize( 3*tCount , Point3D< double >( 0.75 , 0.75 , 0.75 ) );
	visualization.textureCoordinates.resize( 3*tCount );
	visualization.normals.resize( 3*tCount );


	for( int i=0 ; i<tCount ; i++ ) for( int k=0 ; k<3 ; k++ ) visualization.triangles[i][k] = 3*i + k;

	for( int i=0 ; i<tCount ; i++ ) for( int j=0 ; j<3 ; j++ )
	{
		visualization.vertices[ 3*i+j ] = mesh.vertices[ mesh.triangles[i][j] ];
		visualization.normals[ 3*i+j ] = mesh.normals[ mesh.triangles[i][j] ];
		visualization.textureCoordinates[ 3*i+j ] = mesh.textureCoordinates[ 3*i+j ];
	}

	std::vector< int > boundaryEdges;
	mesh.initializeBoundaryEdges( boundaryEdges );

	for( int e=0 ; e<boundaryEdges.size() ; e++ )
	{
		int tIndex = boundaryEdges[e] / 3;
		int kIndex = boundaryEdges[e] % 3;
		for( int c=0 ; c<2 ; c++ )
		{
			Point3D< double > v = mesh.vertices[ mesh.triangles[tIndex][ (kIndex+c)%3 ] ];
			visualization.boundaryEdgeVertices.push_back( v );
		}
	}

	if( inputMode==SINGLE_INPUT_MODE ) visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'M' , "toggle mask" , ToggleMaskCallBack ) );
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
	if( inputMode==SINGLE_INPUT_MODE )   visualization.UpdateCompositeTextureBuffer( inputComposition );
	if( inputMode==SINGLE_INPUT_MODE )   visualization.UpdateMaskTextureBuffer( inputColorMask );

}
#endif // NO_VISUALIZATION

template< typename PreReal , typename Real >
void Stitching< PreReal , Real >::Init( void )
{
	sprintf( stepsString , "Steps: 0" );
	levels = std::max< int >( Levels.value , 1 );
	interpolationWeight = InterpolationWeight.value;

	mesh.read( In.values[0] , NULL , DetailVerbose.set );

	// Define centroid and scale for visualization
	Point3D< PreReal > centroid;
	for( int i=0 ; i<mesh.vertices.size() ; i++ ) centroid += mesh.vertices[i];
	centroid /= (int)mesh.vertices.size();
	PreReal radius = 0;
	for( int i=0 ; i<mesh.vertices.size() ; i++ ) radius = std::max< PreReal >( radius , Point3D< PreReal >::Length( mesh.vertices[i]-centroid ) );
	for( int i=0 ; i<mesh.vertices.size() ; i++ ) mesh.vertices[i] = ( mesh.vertices[i]-centroid ) / radius;

	for( int i=0 ; i<mesh.textureCoordinates.size() ; i++ ) mesh.textureCoordinates[i][1] = 1.0-mesh.textureCoordinates[i][1];

	if( RandomJitter.set )
	{
		srand( time(NULL) );
		std::vector< Point2D< PreReal > > randomOffset( mesh.vertices.size() );
		PreReal jitterScale = (PreReal)1e-3 / std::max< int >( textureWidth , textureHeight );
		for( int i=0 ; i<randomOffset.size() ; i++ ) randomOffset[i] = Point2D< PreReal >( (PreReal)1. - Random< PreReal >()*2 , (PreReal)1. - Random< PreReal >()*2 ) * jitterScale;
		for( int i=0 ; i<mesh.triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ ) mesh.textureCoordinates[ 3*i+k ] += randomOffset[ mesh.triangles[i][k] ];
	}

	ComputePadding( padding , textureWidth , textureHeight , mesh.textureCoordinates , DetailVerbose.set );
	if( padding.nonTrivial )
	{
		PadTextureCoordinates( padding , textureWidth , textureHeight , mesh.textureCoordinates );
		if( inputMode==MULTIPLE_INPUT_MODE ) for( int i=0 ; i<numTextures ; i++ ) PadImage( padding , inputTextures[i] ) , PadImage( padding , inputConfidence[i] );
		else PadImage( padding , inputComposition ) , PadImage( padding , inputMask );
		
		textureWidth  += ( padding.left   + padding.right );
		textureHeight += ( padding.bottom + padding.top );		
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
		FEM::RiemannianMesh< PreReal > rMesh( GetPointer( mesh.triangles) , mesh.triangles.size() );
		rMesh.setMetricFromEmbedding( GetPointer(mesh.vertices) );
		rMesh.makeUnitArea();
		Pointer( FEM::CoordinateXForm< PreReal > ) xForms = rMesh.getCoordinateXForms();

		for( int i=0 ; i<textureNodes.size() ; i++ ) if( textureNodes[i].tID!=-1 && !textureNodes[i].isInterior )
		{
			FEM::HermiteSamplePoint< PreReal > _p;
			_p.tIdx = textureNodes[i].tID;
			_p.p = Point2D< PreReal >( (PreReal)1./3 , (PreReal)1./3 );
			_p.v = textureNodes[i].barycentricCoords - _p.p;

			rMesh.exp( xForms , _p );

			textureNodes[i].tID = _p.tIdx;
			textureNodes[i].barycentricCoords = _p.p;
		}
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

	textureEdgePositions.resize( edgePairs.size() );
	for( int i=0 ; i<edgePairs.size() ; i++ ) textureEdgePositions[i] = ( textureNodePositions[ edgePairs[i].first ] + textureNodePositions[ edgePairs[i].second ] ) / 2;

	if( OutputChartMask.set )
	{
		auto IsBlack = []( Point3D< unsigned char > c ){ return !c[0] && !c[1] && !c[2]; };

		Image< Point3D< unsigned char > > mask;
		mask.resize( textureWidth , textureHeight );
		for( unsigned int i=0 ; i<textureWidth ; i++ ) for( unsigned int j=0 ; j<textureHeight ; j++ ) mask(i,j) = Point3D< unsigned char >(0,0,0);

		std::map< int , Point3D< unsigned char > > chartColors;
		for( unsigned int i=0 ; i<textureNodes.size() ; i++ ) if( chartColors.find( textureNodes[i].chartID )==chartColors.end() )
		{
			Point3D< unsigned char > c;
			while( IsBlack(c) ) c[0] = rand()%256 , c[1] = rand()%256 , c[2] = rand()%256;
			chartColors[ textureNodes[i].chartID ] = c;
		}

		for( int i=0 ; i<textureNodes.size() ; i++ ) mask( textureNodes[i].ci , textureNodes[i].cj ) = chartColors[ textureNodes[i].chartID ];
		for( int e=0 ; e<OutputChartMaskErode.value ; e++ )
		{
			Image< Point3D< unsigned char > > _mask = mask;
			for( int i=0 ; i<textureWidth ; i++ ) for( int j=0 ; j<textureHeight ; j++ ) if( mask(i,j)[0] || mask(i,j)[1] || mask(i,j)[2] )
				for( int di=-1 ; di<=1 ; di++ ) for( int dj=-1 ; dj<=1 ; dj++ )
					if( i+di>=0 && i+di<textureWidth && j+dj>0 && j+dj<textureHeight ) if( IsBlack( _mask(i+di,j+dj) ) ) mask(i,j) = Point3D< unsigned char >();
		}
		mask.write( OutputChartMask.value );
	}
	{
		int multiChartTexelCount = 0;
		Image< int > texelId;
		texelId.resize( textureWidth , textureHeight );
		for( int i=0 ; i<texelId.size() ; i++ ) texelId[i] = -1;
		for( int i=0 ; i<textureNodes.size() ; i++ )
		{
			int ci = textureNodes[i].ci , cj = textureNodes[i].cj;
			if( texelId(ci,cj)!=-1 ) multiChartTexelCount++;
			texelId(ci,cj) = i;
		}
		if( multiChartTexelCount ) Miscellany::Warn( "%d texels belong to multiple charts" , multiChartTexelCount );
	}
}

template< typename PreReal , typename Real >
void _main( int argc , char *argv[] )
{
	Stitching< PreReal , Real >::inputMode = MultiInput.set ? MULTIPLE_INPUT_MODE : SINGLE_INPUT_MODE;
	Stitching< PreReal , Real >::updateCount = MultiInput.set ? -1 : 0;

	Stitching< PreReal , Real >::LoadImages();
	Stitching< PreReal , Real >::Init();
	Stitching< PreReal , Real >::ParseImages();
	Stitching< PreReal , Real >::SetUpSystem();

#ifdef NO_VISUALIZATION
	Stitching< PreReal , Real >::SolveSystem();
	Stitching< PreReal , Real >::WriteTexture( Output.value );
#else // !NO_VISUALIZATION
	if( !Output.set )
	{
		glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
		Stitching< PreReal , Real >::visualization.visualizationMode = Stitching< PreReal , Real >::inputMode;
		Stitching< PreReal , Real >::visualization.displayMode = TWO_REGION_DISPLAY;
		Stitching< PreReal , Real >::visualization.screenWidth = 1600;
		Stitching< PreReal , Real >::visualization.screenHeight = 800;

		glutInitWindowSize( Stitching< PreReal , Real >::visualization.screenWidth , Stitching< PreReal , Real >::visualization.screenHeight );
		glutInit( &argc , argv );
		char windowName[1024];
		sprintf( windowName , "Stitching" );
		glutCreateWindow( windowName );
		if( glewInit()!=GLEW_OK ) Miscellany::Throw( "glewInit failed" );
		glutDisplayFunc ( Stitching< PreReal , Real >::Display );
		glutReshapeFunc ( Stitching< PreReal , Real >::Reshape );
		glutMouseFunc   ( Stitching< PreReal , Real >::MouseFunc );
		glutMotionFunc  ( Stitching< PreReal , Real >::MotionFunc );
		glutKeyboardFunc( Stitching< PreReal , Real >::KeyboardFunc );
		glutIdleFunc    ( Stitching< PreReal , Real >::Idle );
		if( CameraConfig.set ) Stitching< PreReal , Real >::visualization.ReadSceneConfigurationCallBack( &Stitching< PreReal , Real >::visualization , CameraConfig.value );
		Stitching< PreReal , Real >::InitializeVisualization();
		glutMainLoop();
	}
	else
	{
		Stitching< PreReal , Real >::SolveSystem();
		Stitching< PreReal , Real >::ExportTextureCallBack( &Stitching< PreReal , Real >::visualization , Output.value );
	}
#endif // NO_VISUALIZATION
}

int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , argv+1 , params );
#ifdef NO_VISUALIZATION
	if( !In.set || !Output.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
#else // !NO_VISUALIZATION
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
#endif // NO_VISUALIZATION
	omp_set_num_threads( Threads.value );
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
		if( Double.set ) _main< double , double >( argc , argv );
		else             _main< double , float  >( argc , argv );
	}
	catch( Miscellany::Exception &e )
	{
		printf( "%s\n" , e.what() );
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
