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

#include <Misha/CmdLineParser.h> 
#include <Misha/Miscellany.h>
#include <Misha/FEM.h>
#include <Src/Hierarchy.h>
#include <Src/SimpleMesh.h>
#include <Src/Basis.h>
#include <Src/Solver.h>
#include <Src/MassAndStiffness.h>
#include <Src/InteriorTexelToCellLines.inl>
#include <Src/Padding.h>
#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
#include <Src/TextureFilteringVisualization.h>

const std::string vertex_shader_src =
#include <Shaders/normal_texture_vertex.vs>
;
const std::string fragment_shader_src =
#include <Shaders/normal_texture_fragment.fs>
;
#endif //NO_VISUALIZATION

cmdLineParameterArray< char* , 2 > Input( "in" );
cmdLineParameter< char* > Output( "out" );
cmdLineParameter< int   > OutputVCycles( "outVCycles" , 6 );
cmdLineParameter< float > InterpolationWeight( "interpolation" , 1e3 );
cmdLineParameter< float > GradientModulation( "modulation" , 1.0 );
#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
cmdLineParameter< int   > DisplayMode("display",FOUR_REGION_DISPLAY);
#endif // NO_VISUALIZATION
cmdLineParameter< int   > Threads("threads", omp_get_num_procs());
cmdLineParameter< int   > Levels("levels", 4);
cmdLineParameter< int   > MatrixQuadrature( "mQuadrature" , 6 );

cmdLineParameter< int   > MultigridBlockHeight("mBlockH", 16);
cmdLineParameter< int   > MultigridBlockWidth ("mBlockW", 128);
cmdLineParameter< int   > MultigridPaddedHeight("mPadH", 0);
cmdLineParameter< int   > MultigridPaddedWidth("mPadW", 2);

cmdLineReadable RandomJitter("jitter");
cmdLineParameter< char* > CameraConfig("camera");
cmdLineReadable UseDirectSolver("useDirectSolver");
cmdLineReadable Verbose( "verbose" );
cmdLineReadable NoHelp( "noHelp" );
cmdLineReadable DetailVerbose( "detail" );
cmdLineReadable Double( "double" );
cmdLineReadable* params[] =
{
	&Input , &Output , &InterpolationWeight , &GradientModulation , &CameraConfig , &Levels , &UseDirectSolver , &Threads  , &Verbose ,
	&DetailVerbose , &MultigridBlockHeight , &MultigridBlockWidth , &MultigridPaddedHeight , &MultigridPaddedWidth , &RandomJitter ,
#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
	&DisplayMode ,
#endif // NO_VISUALIZATION
	&Double ,
	&MatrixQuadrature ,
	&OutputVCycles ,
	&NoHelp ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );

	printf( "\t --%s <input mesh and texture>\n" , Input.name );
#ifdef NO_VISUALIZATION
	printf( "\t --%s <output texture>\n" , Output.name );
#else // !NO_VISUALIZATION
	printf( "\t[--%s <output texture>]\n" , Output.name );
#endif // NO_VISUALIZATION
	printf( "\t[--%s <output v-cycles>=%d]\n" , OutputVCycles.name , OutputVCycles.value );
	printf( "\t[--%s <interpolation weight>=%f]\n" , InterpolationWeight.name , InterpolationWeight.value );
	printf( "\t[--%s <gradient modulation>=%f]\n" , GradientModulation.name , GradientModulation.value );
	printf( "\t[--%s <system matrix quadrature points per triangle>=%d]\n" , MatrixQuadrature.name , MatrixQuadrature.value );
	printf( "\t[--%s]\n" , UseDirectSolver.name );
	printf( "\t[--%s]\n" , RandomJitter.name );
	printf( "\t[--%s]\n" , Verbose.name );

	printf( "\t[--%s <camera configuration file>]\n" , CameraConfig.name );
	printf( "\t[--%s <hierarchy levels>=%d]\n" , Levels.name , Levels.value );
	printf( "\t[--%s <threads>=%d]\n" , Threads.name , Threads.value );
	printf( "\t[--%s]\n" , DetailVerbose.name );
#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
	printf( "\t[--%s <display mode>=%d]\n" , DisplayMode.name , DisplayMode.value );
	printf( "\t\t%d] One Region \n"   , ONE_REGION_DISPLAY   );
	printf( "\t\t%d] Two Region \n"   , TWO_REGION_DISPLAY   );
	printf( "\t\t%d] Three Region \n" , THREE_REGION_DISPLAY );
	printf( "\t\t%d] Four Region \n"  , FOUR_REGION_DISPLAY  );
#endif // NO_VISUALIZATION
	printf( "\t[--%s <multigrid block width>=%d]\n"    , MultigridBlockWidth.name    , MultigridBlockWidth.value    );
	printf( "\t[--%s <multigrid block height>=%d]\n"   , MultigridBlockHeight.name   , MultigridBlockHeight.value   );
	printf( "\t[--%s <multigrid padded width>=%d]\n"   , MultigridPaddedWidth.name   , MultigridPaddedWidth.value   );
	printf( "\t[--%s <multigrid padded height>=%d]\n"  , MultigridPaddedHeight.name  , MultigridPaddedHeight.value  );
	printf( "\t[--%s]\n" , NoHelp.name );
}

enum
{
	INPUT_TEXTURE,
	OUTPUT_TEXTURE,
	TEXTURE_COUNT
};

template< typename PreReal , typename Real>
class TextureFilter
{
public:
	static TexturedMesh< PreReal > mesh;
	static int textureWidth;
	static int textureHeight;
	static Real interpolationWeight;
	static Real gradientModulation;
	static int levels;

	static HierarchicalSystem< PreReal , Real > hierarchy;
	static bool gradientModulationUpdated;
	static bool positiveModulation;

	static Image< Point3D< Real > > filteredTexture;

#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
	//UI
	static char gradientModulationStr[1024];
	static char interpolationStr[1024];
#endif // NO_VISUALIZATION

	static std::vector<Real> uniformTexelModulationMask;
	static std::vector<Real> cellModulationMask;
	static std::vector<Real> uniformCellModulationMask;
	static std::vector<Real> texelStiffness[3];

	static std::vector< Point3D< float > > cellCenterPositions;
	static std::vector< Point3D< float> >textureNodePositions;

	static std::vector< AtlasChart< PreReal > > atlasCharts;

	static std::vector< BilinearElementIndex > bilinearElementIndices;
	
	static std::vector< TextureNodeInfo< PreReal > > textureNodes;

	static SparseMatrix< Real , int > mass;
	static SparseMatrix< Real , int > stiffness;
	static SparseMatrix< Real , int > filteringMatrix;

	//RHS computation
	static std::vector<InteriorTexelToCellLine> interiorTexelToCellLines;
	static std::vector< Point3D< Real > > interiorTexelToCellCoeffs;
	static SparseMatrix<Real, int> boundaryCellBasedStiffnessRHSMatrix[3];
	static std::vector<Real> boundaryTexelStiffness[3];
	static std::vector< Point3D< Real > > texelModulatedStiffness;

	static std::vector< Point3D< Real > > mass_x0;
	static std::vector< Point3D< Real > > stiffness_x0;

	static std::vector< SystemCoefficients< Real > > multigridFilteringCoefficients;
	static std::vector< MultigridLevelVariables< Point3D< Real > > > multigridFilteringVariables;
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


	//Linear Operators
	static SystemCoefficients< Real > massCoefficients;
	static SystemCoefficients< Real > stiffnessCoefficients;

	static int steps;
	static char stepsString[];

	static Padding padding;

#ifdef NO_VISUALIZATION
	static int updateCount;
	static void WriteTexture( const char* prompt );
#else // !NO_VISUALIZATION
	//Visulization
	static TextureFilteringVisualization visualization;
	static int updateCount;

	static void ToggleUpdateCallBack( Visualization* v , const char* prompt );
	static void IncrementUpdateCallBack( Visualization* v , const char* prompt );
	static void ExportTextureCallBack(Visualization* v, const char* prompt);

	static void GradientModulationCallBack( Visualization* v , const char* prompt );
	static void InterpolationWeightCallBack( Visualization* v , const char* prompt );
#endif // NO_VISUALIZATION

	static void Init( void );
#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
	static void InitializeVisualization();
#endif // NO_VISUALIZATION
	static void UpdateSolution(bool verbose = false, bool detailVerbose = false);
	static void ComputeExactSolution( bool verbose=false );
	static void InitializeSystem( int width , int height );
	static void _InitializeSystem( std::vector<std::vector< SquareMatrix< PreReal , 2 > > > &parameterMetric , BoundaryProlongationData< Real > &boundaryProlongation , std::vector< Point3D< Real > > &inputSignal , std::vector< Real >& texelToCellCoeffs );

#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
	static void UpdateFilteredColorTexture( const std::vector< Point3D< Real > >& solution );
#endif // NO_VISUALIZATION
	static void UpdateFilteredTexture( const std::vector< Point3D< Real > >& solution );

#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
	static void UpdateMaskTexture();

	static void Display( void ){ visualization.Display(); }
	static void MouseFunc( int button , int state , int x , int y );
	static void MotionFunc( int x , int y );
	static void Reshape( int w , int h ){ visualization.Reshape(w,h); }
	static void KeyboardFunc(unsigned char key, int x, int y){ visualization.KeyboardFunc( key , x , y ); }
	static void Idle( void );
#endif // NO_VISUALIZATION
};

#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
template< typename PreReal , typename Real > char																TextureFilter< PreReal , Real >::gradientModulationStr[1024];
template< typename PreReal , typename Real > char																TextureFilter< PreReal , Real >::interpolationStr[1024];
#endif // NO_VISUALIZATION

template< typename PreReal , typename Real > TexturedMesh< PreReal >											TextureFilter< PreReal , Real >::mesh;
template< typename PreReal , typename Real > int																TextureFilter< PreReal , Real >::textureWidth;
template< typename PreReal , typename Real > int																TextureFilter< PreReal , Real >::textureHeight;
#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
template< typename PreReal , typename Real > TextureFilteringVisualization										TextureFilter< PreReal , Real >::visualization;
#endif // NO_VISUALIZATION
template< typename PreReal , typename Real > SparseMatrix< Real , int >											TextureFilter< PreReal , Real >::mass;
template< typename PreReal , typename Real > SparseMatrix< Real , int >											TextureFilter< PreReal , Real >::stiffness;
template< typename PreReal , typename Real > SparseMatrix< Real , int >											TextureFilter< PreReal , Real >::filteringMatrix;

template< typename PreReal , typename Real > Real																TextureFilter< PreReal , Real >::interpolationWeight;
template< typename PreReal , typename Real > Real																TextureFilter< PreReal , Real >::gradientModulation;

template< typename PreReal , typename Real > std::vector< TextureNodeInfo< PreReal > >							TextureFilter< PreReal , Real >::textureNodes;
template< typename PreReal , typename Real > std::vector< BilinearElementIndex >								TextureFilter< PreReal , Real >::bilinearElementIndices;

template< typename PreReal , typename Real > int																TextureFilter< PreReal , Real >::steps;
template< typename PreReal , typename Real > char																TextureFilter< PreReal , Real >::stepsString[1024];
template< typename PreReal , typename Real > int																TextureFilter< PreReal , Real >::levels;
template< typename PreReal , typename Real > HierarchicalSystem< PreReal , Real >								TextureFilter< PreReal , Real >::hierarchy;


template< typename PreReal , typename Real > bool																TextureFilter< PreReal , Real >::gradientModulationUpdated = true;
template< typename PreReal , typename Real > bool																TextureFilter< PreReal , Real >::positiveModulation = true;

//UI
template< typename PreReal , typename Real > std::vector< Real >												TextureFilter< PreReal , Real >::texelStiffness[3];
template< typename PreReal , typename Real > std::vector< Real >												TextureFilter< PreReal , Real >::cellModulationMask;
template< typename PreReal , typename Real > std::vector< Real >												TextureFilter< PreReal , Real >::uniformCellModulationMask;

template< typename PreReal , typename Real > Image< Point3D< Real > >											TextureFilter< PreReal , Real >::filteredTexture;
template< typename PreReal , typename Real > std::vector< Point3D< Real > >										TextureFilter< PreReal , Real >::stiffness_x0;
template< typename PreReal , typename Real > std::vector< Point3D< Real > >										TextureFilter< PreReal , Real >::mass_x0;

template< typename PreReal , typename Real > std::vector< SystemCoefficients< Real > >							TextureFilter< PreReal , Real >::multigridFilteringCoefficients;
template< typename PreReal , typename Real > std::vector< MultigridLevelVariables< Point3D< Real > > >			TextureFilter< PreReal , Real >::multigridFilteringVariables;
template< typename PreReal , typename Real > std::vector<MultigridLevelIndices<Real>>							TextureFilter< PreReal , Real >::multigridIndices;

template< typename PreReal , typename Real > VCycleSolvers< typename TextureFilter< PreReal , Real >::DirectSolver >		TextureFilter< PreReal , Real >::vCycleSolvers;
template< typename PreReal , typename Real > typename TextureFilter< PreReal , Real >::DirectSolver				TextureFilter< PreReal , Real >::directSolver;

template< typename PreReal , typename Real > std::vector< AtlasChart< PreReal > >								TextureFilter< PreReal , Real >::atlasCharts;

template< typename PreReal , typename Real > std::vector<InteriorTexelToCellLine>								TextureFilter< PreReal , Real >::interiorTexelToCellLines;
template< typename PreReal , typename Real > std::vector< Point3D< Real > >										TextureFilter< PreReal , Real >::interiorTexelToCellCoeffs;
template< typename PreReal , typename Real > SparseMatrix<Real, int>											TextureFilter< PreReal , Real >::boundaryCellBasedStiffnessRHSMatrix[3];
template< typename PreReal , typename Real > std::vector<Real>													TextureFilter< PreReal , Real >::boundaryTexelStiffness[3];
template< typename PreReal , typename Real > std::vector< Point3D< Real > >										TextureFilter< PreReal , Real >::texelModulatedStiffness;

template< typename PreReal , typename Real > std::vector< Point3D< float > >									TextureFilter< PreReal , Real >::cellCenterPositions;
template< typename PreReal , typename Real > std::vector< Point3D< float > >									TextureFilter< PreReal , Real >::textureNodePositions;
template< typename PreReal , typename Real > std::vector< Real >												TextureFilter< PreReal , Real >::uniformTexelModulationMask;

template< typename PreReal , typename Real > Padding															TextureFilter< PreReal , Real >::padding;

template< typename PreReal , typename Real > SystemCoefficients< Real >											TextureFilter< PreReal , Real >::massCoefficients;
template< typename PreReal , typename Real > SystemCoefficients< Real >											TextureFilter< PreReal , Real >::stiffnessCoefficients;
template< typename PreReal , typename Real > int																TextureFilter< PreReal , Real >::updateCount = -1;


#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
template< typename PreReal , typename Real >
void TextureFilter< PreReal , Real >::UpdateMaskTexture( void )
{
#pragma omp parallel for
	for (int i = 0; i < textureNodes.size(); i++){
		Real texelModulationValue = uniformTexelModulationMask[i];
		if (texelModulationValue != 0.5){
			Point3D< Real > newColor;
			if( texelModulationValue>0.5 )
			{
				texelModulationValue = 2.0 * texelModulationValue - 1.0;
				newColor = Point3D< Real >( (Real)1. , (Real)0. , (Real)0. ) * texelModulationValue + Point3D< Real >( (Real)0.8 , (Real)0.8 , (Real)0.8 ) * ( (Real)1.0 - texelModulationValue);
			}
			else
			{
				texelModulationValue = 2.0 *texelModulationValue;
				newColor = Point3D< Real >( (Real)0. , (Real)0. , (Real)1. ) * ( (Real)1. - texelModulationValue ) + Point3D< Real >( (Real)0.8 , (Real)0.8 , (Real)0.8 ) * texelModulationValue;
			}
			int ci = textureNodes[i].ci;
			int cj = textureNodes[i].cj;
			int offset = 3 * (textureWidth*cj + ci);
			visualization.maskBufferValues[offset + 0] = (unsigned char)(newColor[0] * 255.0);
			visualization.maskBufferValues[offset + 1] = (unsigned char)(newColor[1] * 255.0);
			visualization.maskBufferValues[offset + 2] = (unsigned char)(newColor[2] * 255.0);
		}
	}

	visualization.UpdateMaskTextureBuffer();
}

template< typename PreReal , typename Real >
void TextureFilter< PreReal , Real >::UpdateFilteredColorTexture( const std::vector< Point3D< Real > > & solution )
{
#pragma omp parallel for
	for (int i = 0; i < textureNodes.size(); i++) {
		int ci = textureNodes[i].ci;
		int cj = textureNodes[i].cj;
		int offset = 3 * (textureWidth*cj + ci);
		for (int c = 0; c < 3; c++)
		{
			Real value = std::min< Real >( (Real)1. , std::max< Real >( 0 , solution[i][c] ) );
			visualization.colorTextureBuffer[offset + c] = (unsigned char)(value*255.0);
		}
	}
}
#endif // NO_VISUALIZATION

template< typename PreReal , typename Real >
void TextureFilter< PreReal , Real >::UpdateFilteredTexture( const std::vector< Point3D< Real > >& solution )
{
#pragma omp parallel for
	for( int i=0 ; i<textureNodes.size() ; i++ )
	{
		int ci = textureNodes[i].ci , cj = textureNodes[i].cj;
		filteredTexture(ci,cj) = Point3D< Real >( solution[i][0] , solution[i][1] , solution[i][2] );
	}
}

#ifdef NO_VISUALIZATION
template< typename PreReal , typename Real >
void TextureFilter< PreReal , Real >::WriteTexture( const char* fileName )
{
	UpdateFilteredTexture( multigridFilteringVariables[0].x );
	Image< Point3D< Real > > outputTexture = filteredTexture;
	if( padding.nonTrivial ) UnpadImage( padding , outputTexture );

	char* ext = GetFileExtension( fileName );
	if( !strcasecmp( ext , "normap" ) ) WriteBinaryImage( outputTexture , fileName );
	else outputTexture.write( fileName );
	delete[] ext;
}
#else // !NO_VISUALIZATION
template< typename PreReal , typename Real >
void TextureFilter< PreReal , Real >::Idle( void )
{
	if( visualization.promptCallBack ) visualization.showSlideBar = false;
	else                               visualization.showSlideBar = true;

	auto RescaleFunction = []( float x ){ return x*2.f; };

	float radius = 0.1f;
	float modulationVartiation = 0.2f;
	if( visualization.isBrushActive )
	{
		Point3D< float > selectedPoint;
		bool validSelection = false;
		if( visualization.showMesh ) validSelection = visualization.select( visualization.diskX , visualization.diskY , selectedPoint );

		if( validSelection )
		{			
			float modulationSign = positiveModulation ? 1.f : -1.f;
			modulationVartiation *= modulationSign;
			for( int i=0 ; i<bilinearElementIndices.size() ; i++ )
			{
				float distanceRatio = Point3D< float >::Length( cellCenterPositions[i]-selectedPoint) / radius;
				float factor = 1.f - distanceRatio;
				factor = factor < 0 ? 0 : factor*factor*(-2.0*factor + 3.0);
				float uniformModulationMaskValue = std::max< float >( 0.f , std::min< float >( 1.f , uniformCellModulationMask[i] + modulationVartiation * factor ) );
				uniformCellModulationMask[i] = uniformModulationMaskValue;
				cellModulationMask[i] = RescaleFunction( uniformModulationMaskValue );
			}
			if( true )
			{
				for( int i=0 ; i<textureNodePositions.size() ; i++ )
				{
					float distanceRatio = Point3D< float >::Length( textureNodePositions[i] - selectedPoint ) / radius;
					float factor = 1.f - distanceRatio;
					factor = factor < 0 ? 0 : factor*factor*( -2.f * factor + 3.f );
					float modulationMaskValue = std::max< float >( 0 , std::min< float >( 1.f , uniformTexelModulationMask[i] + modulationVartiation * factor ) );
					uniformTexelModulationMask[i] = modulationMaskValue;
				}
				UpdateMaskTexture();
			}
			steps = 0;
		}
	}
	else if( visualization.isSlideBarActive )
	{
		if( visualization.slideBarCursorOldPosition!=visualization.slideBarCursorPosition )
		{
			float diff = (float)( visualization.slideBarCursorPosition - visualization.slideBarCursorOldPosition );
			visualization.slideBarCursorOldPosition = visualization.slideBarCursorPosition;

			for( int i=0 ; i<bilinearElementIndices.size() ; i++ )
			{
				float uniformModulationMaskValue = std::max< float >( 0.f , std::min< float >( 1.f , uniformCellModulationMask[i] + diff ) );
				uniformCellModulationMask[i] = uniformModulationMaskValue;
				cellModulationMask[i] = RescaleFunction( uniformModulationMaskValue );
			}

			if( true )
			{
				for( int i=0 ; i<textureNodePositions.size() ; i++ )
				{
					float modulationMaskValue = std::max< float >( 0.f , std::min< float >( 1.f , uniformTexelModulationMask[i] + diff ) );
					uniformTexelModulationMask[i] = modulationMaskValue;
				}
				UpdateMaskTexture();
			}
			steps = 0;
		}
	}

	if( updateCount && !UseDirectSolver.set && !visualization.promptCallBack )
	{
		UpdateSolution();

		if( visualization.textureType==COLOR_TEXTURE )
		{
			UpdateFilteredColorTexture( multigridFilteringVariables[0].x );
			visualization.UpdateColorTextureBuffer();
		}
		else
		{
			UpdateFilteredTexture( multigridFilteringVariables[0].x );
			visualization.UpdateTextureBuffer( filteredTexture );
		}
		if( updateCount>0 ) updateCount--;
		steps++;
		sprintf( stepsString , "Steps: %d" , steps );
	}

	glutPostRedisplay();
}

template< typename PreReal , typename Real >
void TextureFilter< PreReal , Real >::MouseFunc( int button , int state , int x , int y )
{
	if( state==GLUT_UP && UseDirectSolver.set && ( visualization.isBrushActive || visualization.isSlideBarActive ) )
	{
		CellStiffnessToTexelStiffness< Real , 3 >( cellModulationMask , interiorTexelToCellLines , interiorTexelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix , boundaryTexelStiffness , hierarchy.gridAtlases[0].boundaryGlobalIndex , texelModulatedStiffness );
		int numTexels = (int)multigridFilteringVariables[0].rhs.size();
#pragma omp parallel for
		for( int i=0 ; i<numTexels ; i++ ) multigridFilteringVariables[0].rhs[i] = mass_x0[i]*interpolationWeight + texelModulatedStiffness[i];
		ComputeExactSolution( DetailVerbose.set );
		if( visualization.textureType==COLOR_TEXTURE )
		{
			UpdateFilteredColorTexture( multigridFilteringVariables[0].x );
			visualization.UpdateColorTextureBuffer();
		}
		else
		{
			UpdateFilteredTexture( multigridFilteringVariables[0].x );
			visualization.UpdateTextureBuffer( filteredTexture );
		}
	}
	visualization.newX = x; visualization.newY = y;
	visualization.rotating = visualization.scaling = visualization.panning = false;
	visualization.isSlideBarActive = false;
	visualization.isBrushActive = false;

	if( state==GLUT_DOWN && glutGetModifiers() & GLUT_ACTIVE_SHIFT )
	{
		visualization.isBrushActive = true;
		gradientModulationUpdated = false;
		visualization.diskX = x;
		visualization.diskY = y;

		if     ( button==GLUT_RIGHT_BUTTON ) positiveModulation = true;
		else if( button==GLUT_LEFT_BUTTON  ) positiveModulation = false;
	}
	else if( visualization.showSlideBar && x>10 && x<visualization.slideBarWidth()-10 && y>18 && y<32 ) // Slide bar update
	{
		visualization.isSlideBarActive = true;
		gradientModulationUpdated = false;
		float slideBarCursorPosition = (float)( x-20.f ) / ( visualization.slideBarWidth()-40 );
		slideBarCursorPosition = std::min< float >(std::max< float >( slideBarCursorPosition , 0.f ) , 1.f );
		visualization.slideBarCursorPosition = slideBarCursorPosition;
	}
	else
	{
		if( visualization.showMesh )
		{
			visualization.newX = x; visualization.newY = y;

			visualization.rotating = visualization.scaling = visualization.panning = false;
			if( ( button==GLUT_LEFT_BUTTON || button==GLUT_RIGHT_BUTTON ) && glutGetModifiers() & GLUT_ACTIVE_CTRL ) visualization.panning = true;
			else if( button==GLUT_LEFT_BUTTON  ) visualization.rotating = true;
			else if( button==GLUT_RIGHT_BUTTON ) visualization.scaling  = true;
		}
	}

	glutPostRedisplay();
}

template< typename PreReal , typename Real >
void TextureFilter< PreReal , Real >::MotionFunc( int x , int y )
{

	if( visualization.isBrushActive )
	{
		gradientModulationUpdated = false;
		visualization.diskX = x;
		visualization.diskY = y;
	}
	else if( visualization.showSlideBar && visualization.isSlideBarActive )
	{
		gradientModulationUpdated = false;
		float slideBarCursorPosition = (float)( x-20.f ) / ( visualization.slideBarWidth()-40 );
		slideBarCursorPosition = std::min< float >(std::max< float >( slideBarCursorPosition , 0.f ) , 1.f );
		visualization.slideBarCursorPosition = slideBarCursorPosition;
	}
	else
	{
		if( visualization.showMesh )
		{
			visualization.oldX = visualization.newX, visualization.oldY = visualization.newY, visualization.newX = x, visualization.newY = y;
			int screenSize = std::min< int >(visualization.screenWidth, visualization.screenHeight);
			float rel_x = (visualization.newX - visualization.oldX) / (float)screenSize * 2;
			float rel_y = (visualization.newY - visualization.oldY) / (float)screenSize * 2;

			float pRight = rel_x * visualization.zoom, pUp = -rel_y * visualization.zoom;
			float pForward = rel_y * visualization.zoom;
			float rRight = -rel_y, rUp = -rel_x;

			if     ( visualization.rotating ) visualization.camera.rotateUp( -rUp ) , visualization.camera.rotateRight( -rRight );
			else if( visualization.scaling  ) visualization.camera.translate( visualization.camera.forward*pForward);
			else if( visualization.panning  ) visualization.camera.translate( -( visualization.camera.right*pRight + visualization.camera.up*pUp ) );
		}
		else {
			visualization.oldX = visualization.newX, visualization.oldY = visualization.newY, visualization.newX = x, visualization.newY = y;

			if (visualization.panning) visualization.xForm.offset[0] -= (visualization.newX - visualization.oldX) / visualization.imageToScreenScale(), visualization.xForm.offset[1] += (visualization.newY - visualization.oldY) / visualization.imageToScreenScale();
			else
			{
				float dz = (float)pow( 1.1f , (float)( visualization.newY-visualization.oldY ) / 8 );
				visualization.xForm.zoom *= dz;
			}
		}
	}
	glutPostRedisplay();
}

template< typename PreReal , typename Real >
void TextureFilter< PreReal , Real >::ToggleUpdateCallBack( Visualization * /*v*/ , const char * /*prompt*/ )
{
	if( updateCount ) updateCount = 0;
	else              updateCount = -1;
}

template< typename PreReal , typename Real >
void TextureFilter< PreReal , Real >::IncrementUpdateCallBack( Visualization * /*v*/ , const char * /*prompt*/ )
{
	if( updateCount<0 ) updateCount = 1;
	else updateCount++;
}

template< typename PreReal , typename Real >
void TextureFilter< PreReal , Real >::ExportTextureCallBack( Visualization * /*v*/ , const char* prompt )
{
	UpdateFilteredTexture( multigridFilteringVariables[0].x );
	Image< Point3D< Real > > outputTexture = filteredTexture;
	if( padding.nonTrivial ) UnpadImage( padding , outputTexture );

	char* ext = GetFileExtension( prompt );
	if( !strcasecmp( ext , "normap" ) ) WriteBinaryImage( outputTexture , prompt );
	else
	{
		if( visualization.textureType==NORMAL_TEXTURE )
		{
			for( int i=0 ; i<outputTexture.size() ; i++ )
			{
				outputTexture[i] /= Point3D< Real >::Length( outputTexture[i] );
				outputTexture[i] = outputTexture[i] * 0.5f + Point3D< Real >( (Real)0.5 , (Real)0.5 , (Real)0.5 );
			}
		}
		outputTexture.write( prompt );
	}
	delete[] ext;
}

template< typename PreReal , class Real >
void  TextureFilter< PreReal , Real >::GradientModulationCallBack( Visualization * /*v*/ , const char* prompt )
{
	gradientModulation = atof( prompt );
#pragma omp parallel for
	for( int i=0 ; i<multigridFilteringVariables[0].rhs.size() ; i++ ) multigridFilteringVariables[0].rhs[i] = mass_x0[i] * interpolationWeight + stiffness_x0[i] * gradientModulation;

	if( UseDirectSolver.set )
	{
		ComputeExactSolution();
		if( visualization.textureType==COLOR_TEXTURE )
		{
			UpdateFilteredColorTexture( multigridFilteringVariables[0].x );
			visualization.UpdateColorTextureBuffer();
		}
		else
		{
			UpdateFilteredTexture( multigridFilteringVariables[0].x );
			visualization.UpdateTextureBuffer( filteredTexture );
		}
	}
	sprintf( gradientModulationStr , "Gradient modulation: %e\n" , gradientModulation );
}

template< typename PreReal , class Real >
void  TextureFilter< PreReal , Real >::InterpolationWeightCallBack( Visualization * /*v*/ , const char* prompt )
{
	interpolationWeight = atof(prompt);
	if( UseDirectSolver.set ) filteringMatrix = mass*interpolationWeight + stiffness;
	Miscellany::Timer timer;
	UpdateLinearSystem( interpolationWeight , (Real)1. , hierarchy , multigridFilteringCoefficients , massCoefficients , stiffnessCoefficients , vCycleSolvers , directSolver , filteringMatrix , DetailVerbose.set , false , UseDirectSolver.set );
	if( Verbose.set ) printf( "\tInitialized multigrid coefficients: %.2f(s)\n" , timer.elapsed() );

#pragma omp parallel for
	for( int i=0 ; i<multigridFilteringVariables[0].rhs.size() ; i++ ) multigridFilteringVariables[0].rhs[i] = mass_x0[i]*interpolationWeight + stiffness_x0[i] * gradientModulation;

	if( UseDirectSolver.set )
	{
		ComputeExactSolution();
		if( visualization.textureType==COLOR_TEXTURE )
		{
			UpdateFilteredColorTexture( multigridFilteringVariables[0].x );
			visualization.UpdateColorTextureBuffer();
		}
		else
		{
			UpdateFilteredTexture( multigridFilteringVariables[0].x );
			visualization.UpdateTextureBuffer( filteredTexture );
		}
	}
	sprintf( interpolationStr , "Interpolation weight: %e\n" , interpolationWeight );
}
#endif // NO_VISUALIZATION

template< typename PreReal , typename Real >
void TextureFilter< PreReal , Real >::ComputeExactSolution( bool verbose )
{
	Miscellany::Timer timer;
	solve( directSolver , multigridFilteringVariables[0].x , multigridFilteringVariables[0].rhs );
	if( verbose ) printf( "Solving time =  %.4f\n" , timer.elapsed() );
}

template< typename PreReal , typename Real >
void TextureFilter< PreReal , Real >::UpdateSolution( bool verbose , bool detailVerbose )
{
	if( !gradientModulationUpdated )
	{
		int numTexels = (int)multigridFilteringVariables[0].rhs.size();

		Miscellany::Timer timer;
		CellStiffnessToTexelStiffness< Real , 3 >(cellModulationMask, interiorTexelToCellLines, interiorTexelToCellCoeffs, boundaryCellBasedStiffnessRHSMatrix, boundaryTexelStiffness, hierarchy.gridAtlases[0].boundaryGlobalIndex, texelModulatedStiffness);
#pragma omp parallel for
		for( int i=0 ; i<numTexels ; i++ ) multigridFilteringVariables[0].rhs[i] = mass_x0[i]*interpolationWeight + texelModulatedStiffness[i];

		if( verbose ) printf( "RHS update time %.4f\n" , timer.elapsed() );	
		gradientModulationUpdated = true;
	}

	VCycle( multigridFilteringVariables , multigridFilteringCoefficients , multigridIndices , vCycleSolvers , verbose , detailVerbose );
}

template< typename PreReal , typename Real >
void TextureFilter< PreReal , Real >::_InitializeSystem( std::vector< std::vector< SquareMatrix< PreReal , 2 > > > &parameterMetric , BoundaryProlongationData< Real > &boundaryProlongation , std::vector< Point3D< Real > > &inputSignal , std::vector< Real > &texelToCellCoeffs )
{
	Miscellany::Timer timer;
	{
		switch( MatrixQuadrature.value )
		{
		case  1: InitializeMassAndStiffness< 1>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , true , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix ) ; break;
		case  3: InitializeMassAndStiffness< 3>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , true , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix ) ; break;
		case  6: InitializeMassAndStiffness< 6>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , true , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix ) ; break;
		case 12: InitializeMassAndStiffness<12>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , true , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix ) ; break;
		case 24: InitializeMassAndStiffness<24>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , true , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix ) ; break;
		case 32: InitializeMassAndStiffness<32>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , true , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix ) ; break;
		default: Miscellany::Throw( "Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles" );
		}
	}
	if( Verbose.set ) printf( "\tInitialized mass and stiffness: %.2f(s)\n" , timer.elapsed() );
}

template< typename PreReal , typename Real >
void TextureFilter< PreReal , Real >::InitializeSystem( int width , int height )
{
	Miscellany::Timer timer;

	MultigridBlockInfo multigridBlockInfo(MultigridBlockWidth.value, MultigridBlockHeight.value,MultigridPaddedWidth.value,MultigridPaddedHeight.value, 0);
	InitializeHierarchy( mesh , width , height , levels , textureNodes , bilinearElementIndices , hierarchy , atlasCharts , multigridBlockInfo , true , DetailVerbose.set );
	if( Verbose.set ) printf( "\tInitialized hierarchy: %.2f(s)\n" , timer.elapsed() );

	BoundaryProlongationData< Real > boundaryProlongation;
	InitializeBoundaryProlongationData( hierarchy.gridAtlases[0] , boundaryProlongation );

	std::vector< Point3D< Real > > _x0;
	_x0.resize(textureNodes.size());
	
	for( int i=0 ; i<textureNodes.size() ; i++ )
	{
		Point3D< Real > texelValue = mesh.texture(textureNodes[i].ci, textureNodes[i].cj);
		_x0[i] = Point3D< Real >( texelValue[0] , texelValue[1] , texelValue[2] );
	}

	std::vector< Point3D< Real > > inputSignal(textureNodes.size());
	for( int i=0 ; i<textureNodes.size() ; i++ ) inputSignal[i] = mesh.texture( textureNodes[i].ci , textureNodes[i].cj );
	std::vector< Real > texelToCellCoeffs;

	timer.reset();
	std::vector< std::vector< SquareMatrix< PreReal , 2 > > > parameterMetric;
	InitializeMetric( mesh , EMBEDDING_METRIC , atlasCharts , parameterMetric );
	_InitializeSystem( parameterMetric , boundaryProlongation , inputSignal , texelToCellCoeffs );

	interiorTexelToCellCoeffs.resize(4 * hierarchy.gridAtlases[0].numDeepTexels);
	for( int i=0 ; i<4*hierarchy.gridAtlases[0].numDeepTexels ; i++ ) interiorTexelToCellCoeffs[i] = Point3D< Real >( Real(texelToCellCoeffs[3*i+0]) , Real(texelToCellCoeffs[3*i+1]) , Real(texelToCellCoeffs[3*i+2]) );
	
	InitializeInteriorTexelToCellLines( interiorTexelToCellLines , hierarchy.gridAtlases[0] );

	for (int c = 0; c < 3; c++) boundaryTexelStiffness[c].resize(hierarchy.gridAtlases[0].boundaryGlobalIndex.size());
	texelModulatedStiffness.resize(hierarchy.gridAtlases[0].numTexels);

	if( UseDirectSolver.set )
	{
		FullMatrixConstruction( hierarchy.gridAtlases[0] , massCoefficients , mass );
		FullMatrixConstruction( hierarchy.gridAtlases[0] , stiffnessCoefficients , stiffness );
		filteringMatrix  = mass*interpolationWeight + stiffness;
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

	timer.reset();
	UpdateLinearSystem( interpolationWeight , (Real)1. , hierarchy , multigridFilteringCoefficients , massCoefficients , stiffnessCoefficients , vCycleSolvers , directSolver , filteringMatrix , DetailVerbose.set , true , UseDirectSolver.set );
	if( Verbose.set ) printf( "\tInitialized multigrid coefficients: %.2f(s)\n" , timer.elapsed() );

	multigridFilteringVariables.resize(levels);
	for( int i=0 ; i<levels ; i++ )
	{
		MultigridLevelVariables< Point3D< Real > >& variables = multigridFilteringVariables[i];
		variables.x.resize(hierarchy.gridAtlases[i].numTexels);
		variables.rhs.resize(hierarchy.gridAtlases[i].numTexels);
		variables.residual.resize(hierarchy.gridAtlases[i].numTexels);
		variables.boundary_rhs.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
		variables.boundary_value.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
		variables.variable_boundary_value.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
	}

	mass_x0.resize(textureNodes.size());
	MultiplyBySystemMatrix_NoReciprocals( massCoefficients , hierarchy.gridAtlases[0].boundaryGlobalIndex , hierarchy.gridAtlases[0].rasterLines , _x0 , mass_x0 );

	stiffness_x0.resize(textureNodes.size());
	MultiplyBySystemMatrix_NoReciprocals( stiffnessCoefficients , hierarchy.gridAtlases[0].boundaryGlobalIndex , hierarchy.gridAtlases[0].rasterLines , _x0 , stiffness_x0 );

#pragma omp parallel for
	for (int i = 0; i <_x0.size(); i++){
		multigridFilteringVariables[0].x[i] = _x0[i];
		multigridFilteringVariables[0].rhs[i] = mass_x0[i] * interpolationWeight + stiffness_x0[i] * gradientModulation;
	}

	filteredTexture.resize(width, height);
	for (int i = 0; i < filteredTexture.size(); i++) filteredTexture[i] = Point3D< Real >( (Real)0.5 , (Real)0.5 , (Real)0.5 );


	if( UseDirectSolver.set ) ComputeExactSolution( Verbose.set );
	UpdateFilteredTexture( multigridFilteringVariables[0].x );
}

#ifdef NO_VISUALIZATION
#else // !NO_VISUALIZATION
template< typename PreReal , typename Real >
void TextureFilter< PreReal , Real >::InitializeVisualization( void )
{
	sprintf( gradientModulationStr , "Gradient modulation: %.2e\n" , gradientModulation );
	sprintf( interpolationStr , "Interpolation: %.2e\n" , interpolationWeight );

	visualization.textureWidth = textureWidth;
	visualization.textureHeight = textureHeight;

	visualization.colorTextureBuffer = new unsigned char[textureHeight*textureWidth * 3];
	memset(visualization.colorTextureBuffer, 128, textureHeight * textureWidth * 3 * sizeof(unsigned char));



	int tCount = (int)mesh.triangles.size();

	visualization.triangles.resize(tCount);
	visualization.vertices.resize(3 * tCount);
	visualization.colors.resize( 3*tCount , Point3D< float >( 0.75f , 0.75f , 0.75f ) );
	visualization.textureCoordinates.resize( 3*tCount );
	visualization.normals.resize( 3*tCount );


	for (int i = 0; i < tCount; i++) for (int k = 0; k < 3; k++) visualization.triangles[i][k] = 3 * i + k;

	for (int i = 0; i<tCount; i++){
		for (int j = 0; j < 3; j++){
			visualization.vertices[3 * i + j] = mesh.vertices[mesh.triangles[i][j]];
			visualization.normals[3 * i + j] = mesh.normals[mesh.triangles[i][j]];
			visualization.textureCoordinates[3 * i + j] = mesh.textureCoordinates[3 * i + j];
		}
	}

	std::vector<int> boundaryEdges;
	mesh.initializeBoundaryEdges( boundaryEdges);

	for( int e=0 ; e<boundaryEdges.size() ; e++ )
	{
		int tIndex = boundaryEdges[e] / 3;
		int kIndex = boundaryEdges[e] % 3;
		for (int c = 0; c < 2; c++)
		{
			Point3D< PreReal > v = mesh.vertices[ mesh.triangles[tIndex][ (kIndex+c)%3 ] ];
			visualization.boundaryEdgeVertices.push_back(v);
		}
	}

	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 's' , "export texture" , "Output Texture", ExportTextureCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'y' , "interpolation weight" , "Interpolation Weight" , InterpolationWeightCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'g' , "gradient modulation" , "Gradient Modulation" , GradientModulationCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , ' ' , "toggle update" , ToggleUpdateCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , '+' , "increment update" , IncrementUpdateCallBack ) );
	visualization.info.push_back( stepsString );

	visualization.info.push_back( gradientModulationStr );
	visualization.info.push_back( interpolationStr );

	char* ext = GetFileExtension( Input.values[1] );
	if( !strcasecmp( ext , "normap" ) )
	{
		visualization.textureType = NORMAL_TEXTURE;
		visualization.normalProgram = new GLSLProgram( vertex_shader_src , fragment_shader_src );
		visualization.normalProgram->bindAttribLocation( 0 , "vertex_position" );
		visualization.normalProgram->bindAttribLocation( 1 , "vertex_texture" );
		visualization.normalProgram->bindAttribLocation( 2 , "vertex_normal" );
//		visualization.normalProgram->bindFragDataLocation( 0 , "frag_color" );
		visualization.normalProgram->setup();
	}
	else visualization.textureType = COLOR_TEXTURE;
	delete[] ext;

	visualization.UpdateVertexBuffer();
	visualization.UpdateFaceBuffer();
	visualization.UpdateTextureBuffer( filteredTexture );

	visualization.maskBufferValues = new unsigned char[textureHeight*textureWidth * 3];
	memset(visualization.maskBufferValues, 128, textureHeight * textureWidth * 3 * sizeof(unsigned char));
	for (int i = 0; i < textureNodes.size(); i++) {
			int ci = textureNodes[i].ci;
			int cj = textureNodes[i].cj;
			int offset = 3 * (textureWidth*cj + ci);
			visualization.maskBufferValues[offset + 0] = (unsigned char)(0.8 * 255.0);
			visualization.maskBufferValues[offset + 1] = (unsigned char)(0.8 * 255.0);
			visualization.maskBufferValues[offset + 2] = (unsigned char)(0.8 * 255.0);
	}
	visualization.UpdateMaskTextureBuffer();
}
#endif // NO_VISUALIZATION

template< typename PreReal , typename Real >
void TextureFilter< PreReal , Real >::Init( void )
{
	sprintf( stepsString , "Steps: 0" );
	levels = std::max<int>(Levels.value,1);
	interpolationWeight = InterpolationWeight.value;
	gradientModulation = GradientModulation.value;

	mesh.read( Input.values[0] , Input.values[1] , DetailVerbose.set );

	textureWidth = mesh.texture.width();
	textureHeight = mesh.texture.height();

	//Define centroid and scale for visualization
	Point3D< PreReal > centroid;
	for( int i=0 ; i<mesh.vertices.size() ; i++ ) centroid += mesh.vertices[i];
	centroid /= (int)mesh.vertices.size();
	PreReal radius = 0;
	for( int i=0 ; i<mesh.vertices.size() ; i++ ) radius = std::max< PreReal >( radius , Point3D< PreReal >::Length( mesh.vertices[i] - centroid ) );
	for (int i = 0; i < mesh.vertices.size(); i++) mesh.vertices[i] = (mesh.vertices[i] - centroid) / radius;

	if( true ) for( int i=0 ; i<mesh.textureCoordinates.size() ; i++ ) mesh.textureCoordinates[i][1] = (PreReal)1. - mesh.textureCoordinates[i][1];

	if( RandomJitter.set )
	{
		srand( time(NULL) );
		std::vector< Point2D< PreReal > >randomOffset( mesh.vertices.size() );
		PreReal jitterScale = (PreReal)1e-3 / std::max< int >( textureWidth , textureHeight );
		for( int i=0 ; i<randomOffset.size() ; i++ ) randomOffset[i] = Point2D< PreReal >( (PreReal)1. - Random< PreReal >()*2 , (PreReal)1. - Random< PreReal >()*2 )*jitterScale;
		for( int i=0 ; i<mesh.triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ ) mesh.textureCoordinates[ 3*i+k ] += randomOffset[ mesh.triangles[i][k] ];
	}

	ComputePadding( padding , textureWidth , textureHeight , mesh.textureCoordinates , DetailVerbose.set );
	if( padding.nonTrivial )
	{
		PadTextureCoordinates(padding, textureWidth, textureHeight, mesh.textureCoordinates);
		PadImage( padding , mesh.texture );
		textureWidth += (padding.left + padding.right);
		textureHeight += (padding.bottom + padding.top);
	}

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
		FEM::RiemannianMesh< PreReal > rMesh( GetPointer( mesh.triangles ) , mesh.triangles.size() );
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
	}

	textureNodePositions.resize(textureNodes.size());
	for( int i=0 ; i<textureNodePositions.size() ; i++ )
	{
		Point2D< PreReal > barycentricCoords = textureNodes[i].barycentricCoords;
		int tID = textureNodes[i].tID;
		Point3D< PreReal > p =
			mesh.vertices[mesh.triangles[tID][0]] * ( (PreReal)1. - barycentricCoords[0] - barycentricCoords[1] ) +
			mesh.vertices[mesh.triangles[tID][1]] *                 barycentricCoords[0]                          +
			mesh.vertices[mesh.triangles[tID][2]] *                                        barycentricCoords[1]   ;
		textureNodePositions[i] = Point3D< float >( p );
	}

	uniformTexelModulationMask.resize( textureNodes.size() , 0.5 );

	for (int c = 0; c < 3; c++) texelStiffness[c].resize(textureNodes.size());

	cellModulationMask.resize( bilinearElementIndices.size() , 1 );
	uniformCellModulationMask.resize( bilinearElementIndices.size() , 0.5 );

	cellCenterPositions.resize( bilinearElementIndices.size() );
	for( int i=0 ; i<bilinearElementIndices.size() ; i++ )
		cellCenterPositions[i] = Point3D< float >
		(
			textureNodePositions[ bilinearElementIndices[i][0] ] +
			textureNodePositions[ bilinearElementIndices[i][1] ] +
			textureNodePositions[ bilinearElementIndices[i][2] ] +
			textureNodePositions[ bilinearElementIndices[i][3] ]
		) / 4.f;

	if( true )
	{
		int multiChartTexelCount = 0;
		Image< int > texelId;
		texelId.resize(textureWidth, textureHeight);
		for( int i=0 ; i<texelId.size() ; i++ ) texelId[i] = -1;
		for( int i=0 ; i<textureNodes.size() ; i++)
		{
			int ci = textureNodes[i].ci , cj = textureNodes[i].cj;
			if( texelId(ci,cj)!=-1 ) multiChartTexelCount++;
			texelId(ci,cj) = i;
		}
		if( multiChartTexelCount ) Miscellany::Warn( "%d texels belong to multiple charts!\n" , multiChartTexelCount );
	}
}

template< typename PreReal , typename Real >
void _main( int argc , char *argv[] )
{
	TextureFilter< PreReal , Real >::Init();

#ifdef NO_VISUALIZATION
	if( UseDirectSolver.set ) TextureFilter< PreReal , Real >::ComputeExactSolution( Verbose.set );
	else for ( int i=0 ; i<OutputVCycles.value ; i++ ) TextureFilter< PreReal , Real >::UpdateSolution();
	TextureFilter< PreReal , Real >::WriteTexture( Output.value );
#else // !NO_VISUALIZATION
	if( !Output.set )
	{
		glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
		TextureFilter< PreReal , Real >::visualization.displayMode = DisplayMode.value;
		if     ( DisplayMode.value==ONE_REGION_DISPLAY   ) TextureFilter< PreReal , Real >::visualization.screenWidth =  800 , TextureFilter< PreReal , Real >::visualization.screenHeight = 800;
		else if( DisplayMode.value==TWO_REGION_DISPLAY   ) TextureFilter< PreReal , Real >::visualization.screenWidth = 1600 , TextureFilter< PreReal , Real >::visualization.screenHeight = 800;
		else if( DisplayMode.value==THREE_REGION_DISPLAY ) TextureFilter< PreReal , Real >::visualization.screenWidth = 1200 , TextureFilter< PreReal , Real >::visualization.screenHeight = 800;
		else if( DisplayMode.value==FOUR_REGION_DISPLAY  ) TextureFilter< PreReal , Real >::visualization.screenWidth = 1500 , TextureFilter< PreReal , Real >::visualization.screenHeight = 600;
		glutInitWindowSize( TextureFilter< PreReal , Real >::visualization.screenWidth , TextureFilter< PreReal , Real >::visualization.screenHeight );
		glutInit( &argc , argv );
		char windowName[1024];
		sprintf( windowName , "Texture Filtering" );
		glutCreateWindow( windowName );
		if( glewInit()!=GLEW_OK ) Miscellany::Throw( "glewInit failed" );
		glutDisplayFunc ( TextureFilter< PreReal , Real >::Display );
		glutReshapeFunc ( TextureFilter< PreReal , Real >::Reshape );
		glutMouseFunc   ( TextureFilter< PreReal , Real >::MouseFunc );
		glutMotionFunc  ( TextureFilter< PreReal , Real >::MotionFunc );
		glutKeyboardFunc( TextureFilter< PreReal , Real >::KeyboardFunc) ;
		glutIdleFunc    ( TextureFilter< PreReal , Real >::Idle );
		if( CameraConfig.set ) TextureFilter< PreReal , Real >::visualization.ReadSceneConfigurationCallBack( &TextureFilter< PreReal , Real >::visualization , CameraConfig.value );
		TextureFilter< PreReal , Real >::InitializeVisualization();
		TextureFilter< PreReal , Real >::visualization.showSlideBar = true;
		glutMainLoop(); 
	}
	else
	{
		if( UseDirectSolver.set ) TextureFilter< PreReal , Real >::ComputeExactSolution( Verbose.set );
		else for ( int i=0 ; i<OutputVCycles.value ; i++ ) TextureFilter< PreReal , Real >::UpdateSolution();
		TextureFilter< PreReal , Real >::ExportTextureCallBack( &TextureFilter< PreReal , Real >::visualization , Output.value );
	}
#endif // NO_VISUALIZATION
}

int main(int argc, char* argv[])
{
	cmdLineParse( argc-1 , argv+1 , params );
#ifdef NO_VISUALIZATION
	if( !Input.set || !Output.set ) { ShowUsage( argv[0] ) ; return EXIT_FAILURE; }
#else // !NO_VISUALIZATION
	if( !Input.set ) { ShowUsage( argv[0] ) ; return EXIT_FAILURE; }
#endif // NO_VISUALIZATION
	omp_set_num_threads( Threads.value );
	if( !NoHelp.set && !Output.set )
	{
		printf( "+----------------------------------------------------------------------+\n" );
		printf( "| Interface Controls:                                                  |\n" );
		printf( "|    [Left Mouse]:                rotate                               |\n" );
		printf( "|    [Right Mouse]:               zoom                                 |\n" );
		printf( "|    [Left/Right Mouse] + [CTRL]: pan                                  |\n" );
		printf( "|    [Left Mouse] + [SHIFT]:      dampen gradients                     |\n" );
		printf( "|    [Right Mouse] + [SHIFT]:     amplify gradients                    |\n" );
		printf( "|    'g':                         prescribe global gradient modulation |\n" );
		printf( "|    'y':                         prescribe interpolation weight       |\n" );
		printf( "+----------------------------------------------------------------------+\n" );
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
