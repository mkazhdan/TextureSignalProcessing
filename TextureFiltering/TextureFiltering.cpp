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
#include <Src/MassAndStiffness.h>
#include <Src/InteriorTexelToCellLines.h>
#include <Src/Padding.h>
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
#include <Src/TextureFilteringVisualization.h>

const std::string vertex_shader_src =
#include <Shaders/normal_texture_vertex.vs>
;
const std::string fragment_shader_src =
#include <Shaders/normal_texture_fragment.fs>
;
#endif // NO_OPEN_GL_VISUALIZATION

using namespace MishaK;

CmdLineParameterArray< std::string , 2 > Input( "in" );
CmdLineParameter< std::string > InputLowFrequency( "inLow" );
CmdLineParameter< std::string > Output( "out" );
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
CmdLineParameter< std::string > Snapshot( "snapshot" );
#endif // NO_OPEN_GL_VISUALIZATION
CmdLineParameter< int   > OutputVCycles( "outVCycles" , 6 );
CmdLineParameter< float > InterpolationWeight( "interpolation" , 1e3 );
CmdLineParameter< float > GradientModulation( "modulation" , 1.0 );
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
CmdLineParameter< int   > DisplayMode( "display" , FOUR_REGION_DISPLAY );
#endif // 
CmdLineParameter< int   > Levels("levels", 4);
CmdLineParameter< int   > MatrixQuadrature( "mQuadrature" , 6 );

CmdLineParameter< int   > MultigridBlockHeight ( "mBlockH" , 16 );
CmdLineParameter< int   > MultigridBlockWidth  ( "mBlockW" , 128 );
CmdLineParameter< int   > MultigridPaddedHeight( "mPadH" , 0 );
CmdLineParameter< int   > MultigridPaddedWidth ( "mPadW" , 2 );

CmdLineParameter< int   > RandomJitter( "jitter" , 0 );
CmdLineReadable Paused( "paused" );
CmdLineParameter< std::string > CameraConfig( "camera" );
CmdLineReadable UseDirectSolver( "useDirectSolver" );
CmdLineReadable Verbose( "verbose" );
CmdLineReadable NoHelp( "noHelp" );
CmdLineReadable DetailVerbose( "detail" );
CmdLineReadable Double( "double" );
CmdLineReadable Seamless( "seamless" );
CmdLineReadable Serial( "serial" );
CmdLineReadable ColorAsNormal( "colorAsNormal" );
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
CmdLineReadable Nearest( "nearest" );
#endif // NO_OPEN_GL_VISUALIZATION

CmdLineParameter< double > CollapseEpsilon( "collapse" , 0 );

CmdLineReadable* params[] =
{
	&Input , &Output , &InterpolationWeight , &GradientModulation , &CameraConfig , &Levels , &UseDirectSolver , &Serial  , &Verbose ,
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	&Snapshot ,
#endif // NO_OPEN_GL_VISUALIZATION
	&InputLowFrequency ,
	&DetailVerbose , &MultigridBlockHeight , &MultigridBlockWidth , &MultigridPaddedHeight , &MultigridPaddedWidth , &RandomJitter ,
	&Paused ,
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	&DisplayMode ,
#endif // NO_OPEN_GL_VISUALIZATION
	&Double ,
	&MatrixQuadrature ,
	&OutputVCycles ,
	&Seamless ,
	&ColorAsNormal ,
	&CollapseEpsilon ,
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	&Nearest ,
#endif // NO_OPEN_GL_VISUALIZATION
	&NoHelp ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );

	printf( "\t --%s <input mesh and texture>\n" , Input.name.c_str() );
	printf( "\t[--%s <input low-frequency texture>\n" , InputLowFrequency.name.c_str() );
#ifdef NO_OPEN_GL_VISUALIZATION
	printf( "\t --%s <output texture>\n" , Output.name.c_str() );
#else // !NO_OPEN_GL_VISUALIZATION
	printf( "\t[--%s <output texture>]\n" , Output.name.c_str() );
	printf( "\t[--%s <snapshot name>]\n" , Snapshot.name.c_str() );
#endif // NO_OPEN_GL_VISUALIZATION
	printf( "\t[--%s <output v-cycles>=%d]\n" , OutputVCycles.name.c_str() , OutputVCycles.value );
	printf( "\t[--%s <interpolation weight>=%f]\n" , InterpolationWeight.name.c_str() , InterpolationWeight.value );
	printf( "\t[--%s <gradient modulation>=%f]\n" , GradientModulation.name.c_str() , GradientModulation.value );
	printf( "\t[--%s <system matrix quadrature points per triangle>=%d]\n" , MatrixQuadrature.name.c_str() , MatrixQuadrature.value );
	printf( "\t[--%s]\n" , Seamless.name.c_str() );
	printf( "\t[--%s]\n" , UseDirectSolver.name.c_str() );
	printf( "\t[--%s <jittering seed>]\n" , RandomJitter.name.c_str() );
	printf( "\t[--%s]\n" , ColorAsNormal.name.c_str() );
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	printf( "\t[--%s]\n" , Nearest.name.c_str() );
#endif // NO_OPEN_GL_VISUALIZATION
	printf( "\t[--%s]\n" , Verbose.name.c_str() );

	printf( "\t[--%s <camera configuration file>]\n" , CameraConfig.name.c_str() );
	printf( "\t[--%s <hierarchy levels>=%d]\n" , Levels.name.c_str() , Levels.value );
	printf( "\t[--%s]\n" , DetailVerbose.name.c_str() );
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	printf( "\t[--%s <display mode>=%d]\n" , DisplayMode.name.c_str() , DisplayMode.value );
	printf( "\t\t%d] One Region \n"   , ONE_REGION_DISPLAY   );
	printf( "\t\t%d] Two Region \n"   , TWO_REGION_DISPLAY   );
	printf( "\t\t%d] Three Region \n" , THREE_REGION_DISPLAY );
	printf( "\t\t%d] Four Region \n"  , FOUR_REGION_DISPLAY  );
#endif // NO_OPEN_GL_VISUALIZATION
	printf( "\t[--%s <multigrid block width>=%d]\n"    , MultigridBlockWidth.name.c_str()    , MultigridBlockWidth.value    );
	printf( "\t[--%s <multigrid block height>=%d]\n"   , MultigridBlockHeight.name.c_str()   , MultigridBlockHeight.value   );
	printf( "\t[--%s <multigrid padded width>=%d]\n"   , MultigridPaddedWidth.name.c_str()   , MultigridPaddedWidth.value   );
	printf( "\t[--%s <multigrid padded height>=%d]\n"  , MultigridPaddedHeight.name.c_str()  , MultigridPaddedHeight.value  );
	printf( "\t[--%s <collapse epsilon>=%g]\n" , CollapseEpsilon.name.c_str() , CollapseEpsilon.value );
	printf( "\t[--%s]\n" , Serial.name.c_str() );
	printf( "\t[--%s]\n" , NoHelp.name.c_str() );
	printf( "\t[--%s]\n" , Paused.name.c_str() );
}

enum
{
	INPUT_TEXTURE,
	OUTPUT_TEXTURE,
	TEXTURE_COUNT
};

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
class TextureFilter
{
public:
	static TexturedTriangleMesh< PreReal > mesh;
	static int textureWidth;
	static int textureHeight;
	static Real interpolationWeight;
	static Real gradientModulation;
	static int levels;

	static HierarchicalSystem< PreReal , Real > hierarchy;
	static bool gradientModulationUpdated;
	static bool positiveModulation;

	static Image< Point3D< Real > > filteredTexture;
	static Image< Point3D< float > > highFrequencyTexture;
	static Image< Point3D< float > > lowFrequencyTexture;

#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	//UI
	static char gradientModulationStr[1024];
	static char interpolationStr[1024];
#endif // NO_OPEN_GL_VISUALIZATION

	static std::vector< Real > uniformTexelModulationMask;
	static std::vector< Real > cellModulationMask;
	static std::vector< Real > uniformCellModulationMask;
	static std::vector< Real > texelStiffness[3];

	static std::vector< Point3D< float > > cellCenterPositions;
	static std::vector< Point3D< float > > textureNodePositions;

	static IndexVector< ChartIndex , AtlasChart< PreReal > > atlasCharts;

#ifdef NEW_CODE
	static IndexVector< AtlasCombinedCellIndex , BilinearElementIndex< AtlasCombinedTexelIndex > > bilinearElementIndices;
#else // !NEW_CODE
	static IndexVector< AtlasCombinedCellIndex , BilinearElementIndex< unsigned int > > bilinearElementIndices;
#endif // NEW_CODE
	
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
#elif defined( USE_EIGEN )
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

#ifdef NO_OPEN_GL_VISUALIZATION
	static int updateCount;
	static void WriteTexture( const char* prompt );
#else // !NO_OPEN_GL_VISUALIZATION
	// Visualization
	static TextureFilteringVisualization visualization;
	static int updateCount;

	static void ToggleUpdateCallBack( Visualization* v , const char* prompt );
	static void IncrementUpdateCallBack( Visualization* v , const char* prompt );
	static void ExportTextureCallBack(Visualization* v, const char* prompt);

	static void GradientModulationCallBack( Visualization* v , const char* prompt );
	static void InterpolationWeightCallBack( Visualization* v , const char* prompt );
#endif // NO_OPEN_GL_VISUALIZATION

	static void Init( void );
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	static void InitializeVisualization();
#endif // NO_OPEN_GL_VISUALIZATION
	static void UpdateSolution(bool verbose = false, bool detailVerbose = false);
	static void ComputeExactSolution( bool verbose=false );
	static void InitializeSystem( int width , int height );
	static void _InitializeSystem( IndexVector< ChartIndex , IndexVector< ChartMeshTriangleIndex , SquareMatrix< PreReal , 2 > > > &parameterMetric , BoundaryProlongationData< Real > &boundaryProlongation , std::vector< Point3D< Real > > &inputSignal , std::vector< Real >& texelToCellCoeffs );

#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	static void UpdateFilteredColorTexture( const std::vector< Point3D< Real > >& solution );
#endif // NO_OPEN_GL_VISUALIZATION
	static void UpdateFilteredTexture( const std::vector< Point3D< Real > >& solution );

#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	static void UpdateMaskTexture();

	static void Display( void ){ visualization.Display(); }
	static void MouseFunc( int button , int state , int x , int y );
	static void MotionFunc( int x , int y );
	static void Reshape( int w , int h ){ visualization.Reshape(w,h); }
	static void KeyboardFunc(unsigned char key, int x, int y){ visualization.KeyboardFunc( key , x , y ); }
	static void Idle( void );
#endif // NO_OPEN_GL_VISUALIZATION
};

#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
template< typename PreReal , typename Real , unsigned int TextureBitDepth > char																TextureFilter< PreReal , Real , TextureBitDepth >::gradientModulationStr[1024];
template< typename PreReal , typename Real , unsigned int TextureBitDepth > char																TextureFilter< PreReal , Real , TextureBitDepth >::interpolationStr[1024];
#endif // NO_OPEN_GL_VISUALIZATION

template< typename PreReal , typename Real , unsigned int TextureBitDepth > TexturedTriangleMesh< PreReal >										TextureFilter< PreReal , Real , TextureBitDepth >::mesh;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > int																	TextureFilter< PreReal , Real , TextureBitDepth >::textureWidth;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > int																	TextureFilter< PreReal , Real , TextureBitDepth >::textureHeight;
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
template< typename PreReal , typename Real , unsigned int TextureBitDepth > TextureFilteringVisualization										TextureFilter< PreReal , Real , TextureBitDepth >::visualization;
#endif // NO_OPEN_GL_VISUALIZATION
template< typename PreReal , typename Real , unsigned int TextureBitDepth > SparseMatrix< Real , int >											TextureFilter< PreReal , Real , TextureBitDepth >::mass;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > SparseMatrix< Real , int >											TextureFilter< PreReal , Real , TextureBitDepth >::stiffness;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > SparseMatrix< Real , int >											TextureFilter< PreReal , Real , TextureBitDepth >::filteringMatrix;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > Real																TextureFilter< PreReal , Real , TextureBitDepth >::interpolationWeight;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > Real																TextureFilter< PreReal , Real , TextureBitDepth >::gradientModulation;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< TextureNodeInfo< PreReal > >							TextureFilter< PreReal , Real , TextureBitDepth >::textureNodes;
#ifdef NEW_CODE
template< typename PreReal , typename Real , unsigned int TextureBitDepth > IndexVector< AtlasCombinedCellIndex , BilinearElementIndex< AtlasCombinedTexelIndex > >	TextureFilter< PreReal , Real , TextureBitDepth >::bilinearElementIndices;
#else // !NEW_CODE
template< typename PreReal , typename Real , unsigned int TextureBitDepth > IndexVector< AtlasCombinedCellIndex , BilinearElementIndex< unsigned int > >	TextureFilter< PreReal , Real , TextureBitDepth >::bilinearElementIndices;
#endif // NEW_CODE

template< typename PreReal , typename Real , unsigned int TextureBitDepth > int																	TextureFilter< PreReal , Real , TextureBitDepth >::steps;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > char																TextureFilter< PreReal , Real , TextureBitDepth >::stepsString[1024];
template< typename PreReal , typename Real , unsigned int TextureBitDepth > int																	TextureFilter< PreReal , Real , TextureBitDepth >::levels;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > HierarchicalSystem< PreReal , Real >								TextureFilter< PreReal , Real , TextureBitDepth >::hierarchy;


template< typename PreReal , typename Real , unsigned int TextureBitDepth > bool																TextureFilter< PreReal , Real , TextureBitDepth >::gradientModulationUpdated = true;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > bool																TextureFilter< PreReal , Real , TextureBitDepth >::positiveModulation = true;

//UI
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Real >													TextureFilter< PreReal , Real , TextureBitDepth >::texelStiffness[3];
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Real >													TextureFilter< PreReal , Real , TextureBitDepth >::cellModulationMask;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Real >													TextureFilter< PreReal , Real , TextureBitDepth >::uniformCellModulationMask;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > Image< Point3D< Real > >											TextureFilter< PreReal , Real , TextureBitDepth >::filteredTexture;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > Image< Point3D< float > >											TextureFilter< PreReal , Real , TextureBitDepth >::highFrequencyTexture;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > Image< Point3D< float > >											TextureFilter< PreReal , Real , TextureBitDepth >::lowFrequencyTexture;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Point3D< Real > >										TextureFilter< PreReal , Real , TextureBitDepth >::stiffness_x0;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Point3D< Real > >										TextureFilter< PreReal , Real , TextureBitDepth >::mass_x0;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< SystemCoefficients< Real > >							TextureFilter< PreReal , Real , TextureBitDepth >::multigridFilteringCoefficients;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< MultigridLevelVariables< Point3D< Real > > >			TextureFilter< PreReal , Real , TextureBitDepth >::multigridFilteringVariables;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector<MultigridLevelIndices<Real>>							TextureFilter< PreReal , Real , TextureBitDepth >::multigridIndices;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > VCycleSolvers< typename TextureFilter< PreReal , Real , TextureBitDepth >::DirectSolver >		TextureFilter< PreReal , Real , TextureBitDepth >::vCycleSolvers;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > typename TextureFilter< PreReal , Real , TextureBitDepth >::DirectSolver						TextureFilter< PreReal , Real , TextureBitDepth >::directSolver;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > IndexVector< ChartIndex , AtlasChart< PreReal > >					TextureFilter< PreReal , Real , TextureBitDepth >::atlasCharts;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector<InteriorTexelToCellLine>								TextureFilter< PreReal , Real , TextureBitDepth >::interiorTexelToCellLines;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Point3D< Real > >										TextureFilter< PreReal , Real , TextureBitDepth >::interiorTexelToCellCoeffs;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > SparseMatrix<Real, int>												TextureFilter< PreReal , Real , TextureBitDepth >::boundaryCellBasedStiffnessRHSMatrix[3];
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector<Real>													TextureFilter< PreReal , Real , TextureBitDepth >::boundaryTexelStiffness[3];
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Point3D< Real > >										TextureFilter< PreReal , Real , TextureBitDepth >::texelModulatedStiffness;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Point3D< float > >										TextureFilter< PreReal , Real , TextureBitDepth >::cellCenterPositions;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Point3D< float > >										TextureFilter< PreReal , Real , TextureBitDepth >::textureNodePositions;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Real >													TextureFilter< PreReal , Real , TextureBitDepth >::uniformTexelModulationMask;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > Padding																TextureFilter< PreReal , Real , TextureBitDepth >::padding;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > SystemCoefficients< Real >											TextureFilter< PreReal , Real , TextureBitDepth >::massCoefficients;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > SystemCoefficients< Real >											TextureFilter< PreReal , Real , TextureBitDepth >::stiffnessCoefficients;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > int																	TextureFilter< PreReal , Real , TextureBitDepth >::updateCount = -1;

template< typename Real >
Image< Point3D< Real > > ColorToNormal( const Image< Point3D< Real > > &color )
{
	Image< Point3D< Real > > normal( color.res() );
	for( unsigned int y=0 ; y<(unsigned int)color.res(1) ; y++ ) for( unsigned int x=0 ; x<(unsigned int)color.res(0) ; x++ )
		normal(x,y) = color(x,y) * static_cast< Real >( 2 ) - Point3D< Real >(1.f,1.f,1.f);
	return normal;
}

template< typename Real >
Image< Point3D< Real > > NormalToColor( const Image< Point3D< Real > > &normal )
{
	Image< Point3D< Real > > color( normal.res() );
	for( unsigned int y=0 ; y<(unsigned int)normal.res(1) ; y++ ) for( unsigned int x=0 ; x<(unsigned int)normal.res(0) ; x++ )
		color(x,y) = ( normal(x,y) + Point3D< Real >( static_cast< Real >( 1. ) , static_cast< Real >( 1. ) , static_cast< Real >( 1. ) ) ) / static_cast< Real >( 2 );
	return color;
}

#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void TextureFilter< PreReal , Real , TextureBitDepth >::UpdateMaskTexture( void )
{
	ThreadPool::ParallelFor
		(
			0 , textureNodes.size() ,
			[&]( unsigned int , size_t i )
			{
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
		);

	visualization.UpdateMaskTextureBuffer();
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void TextureFilter< PreReal , Real , TextureBitDepth >::UpdateFilteredColorTexture( const std::vector< Point3D< Real > > & solution )
{
	ThreadPool::ParallelFor
		(
			0 , textureNodes.size() ,
			[&]( unsigned int , size_t i )
			{
				int ci = textureNodes[i].ci;
				int cj = textureNodes[i].cj;
				int offset = 3 * (textureWidth*cj + ci);
				for (int c = 0; c < 3; c++)
				{
					Real value = std::min< Real >( (Real)1. , std::max< Real >( 0 , solution[i][c] ) );
					visualization.colorTextureBuffer[offset + c] = (unsigned char)(value*255.0);
				}
			}
		);
}
#endif // NO_OPEN_GL_VISUALIZATION

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void TextureFilter< PreReal , Real , TextureBitDepth >::UpdateFilteredTexture( const std::vector< Point3D< Real > >& solution )
{
	ThreadPool::ParallelFor
		(
			0 , textureNodes.size() ,
			[&]( unsigned int , size_t i )
			{
				int ci = textureNodes[i].ci , cj = textureNodes[i].cj;
				filteredTexture(ci,cj) = Point3D< Real >( solution[i][0] , solution[i][1] , solution[i][2] );
			}
		);
}

#ifdef NO_OPEN_GL_VISUALIZATION
template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void TextureFilter< PreReal , Real , TextureBitDepth >::WriteTexture( const char* fileName )
{
	UpdateFilteredTexture( multigridFilteringVariables[0].x );
	Image< Point3D< Real > > outputTexture = filteredTexture;
	padding.unpad( outputTexture );

	std::string ext = ToLower( GetFileExtension( fileName ) );
	if( ext==std::string( "normap" ) ) WriteBinaryImage( outputTexture , fileName );
	else
	{
		if( ColorAsNormal.set ) outputTexture = NormalToColor( outputTexture );
		WriteImage< TextureBitDepth >( outputTexture , fileName );
	}
}
#else // !NO_OPEN_GL_VISUALIZATION
template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void TextureFilter< PreReal , Real , TextureBitDepth >::Idle( void )
{
	if( visualization.promptCallBack || visualization.snapshotName ) visualization.showHelp = visualization.showInfo = visualization.showSlideBar = false;
	else                                                             visualization.showSlideBar = true;

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

	visualization.Idle();
	glutPostRedisplay();
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void TextureFilter< PreReal , Real , TextureBitDepth >::MouseFunc( int button , int state , int x , int y )
{
	if( state==GLUT_UP && UseDirectSolver.set && ( visualization.isBrushActive || visualization.isSlideBarActive ) )
	{
		CellStiffnessToTexelStiffness< Real , 3 >( cellModulationMask , interiorTexelToCellLines , interiorTexelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix , boundaryTexelStiffness , hierarchy.gridAtlases[0].indexConverter , texelModulatedStiffness );
		int numTexels = (int)multigridFilteringVariables[0].rhs.size();
		ThreadPool::ParallelFor( 0 , numTexels , [&]( unsigned int , size_t i ){ multigridFilteringVariables[0].rhs[i] = mass_x0[i]*interpolationWeight + texelModulatedStiffness[i]; } );
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

		visualization.positiveModulation = positiveModulation;
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

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void TextureFilter< PreReal , Real , TextureBitDepth >::MotionFunc( int x , int y )
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

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void TextureFilter< PreReal , Real , TextureBitDepth >::ToggleUpdateCallBack( Visualization * /*v*/ , const char * /*prompt*/ )
{
	if( updateCount ) updateCount = 0;
	else              updateCount = -1;
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void TextureFilter< PreReal , Real , TextureBitDepth >::IncrementUpdateCallBack( Visualization * /*v*/ , const char * /*prompt*/ )
{
	if( updateCount<0 ) updateCount = 1;
	else updateCount++;
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void TextureFilter< PreReal , Real , TextureBitDepth >::ExportTextureCallBack( Visualization * /*v*/ , const char* prompt )
{
	UpdateFilteredTexture( multigridFilteringVariables[0].x );
	Image< Point3D< Real > > outputTexture = filteredTexture;
	padding.unpad( outputTexture );

	std::string ext = ToLower( GetFileExtension( prompt ) );
	if( ext==std::string( "normap" ) ) WriteBinaryImage( outputTexture , prompt );
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
		WriteImage< 8 >( outputTexture , prompt );
	}
}

template< typename PreReal , class Real , unsigned int TextureBitDepth >
void  TextureFilter< PreReal , Real , TextureBitDepth >::GradientModulationCallBack( Visualization * /*v*/ , const char* prompt )
{
	gradientModulation = atof( prompt );
	ThreadPool::ParallelFor( 0 , multigridFilteringVariables[0].rhs.size() , [&]( unsigned int , size_t i ){ multigridFilteringVariables[0].rhs[i] = mass_x0[i] * interpolationWeight + stiffness_x0[i] * gradientModulation; } );

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

template< typename PreReal , class Real , unsigned int TextureBitDepth >
void  TextureFilter< PreReal , Real , TextureBitDepth >::InterpolationWeightCallBack( Visualization * /*v*/ , const char* prompt )
{
	interpolationWeight = atof(prompt);
	if( UseDirectSolver.set ) filteringMatrix = mass*interpolationWeight + stiffness;
	Miscellany::Timer timer;
	UpdateLinearSystem( interpolationWeight , (Real)1. , hierarchy , multigridFilteringCoefficients , massCoefficients , stiffnessCoefficients , vCycleSolvers , directSolver , filteringMatrix , DetailVerbose.set , false , UseDirectSolver.set );
	if( Verbose.set ) printf( "\tInitialized multigrid coefficients: %.2f(s)\n" , timer.elapsed() );

	ThreadPool::ParallelFor( 0 ,multigridFilteringVariables[0].rhs.size() , [&]( unsigned int , size_t i ){ multigridFilteringVariables[0].rhs[i] = mass_x0[i]*interpolationWeight + stiffness_x0[i] * gradientModulation; } );

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
#endif // NO_OPEN_GL_VISUALIZATION

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void TextureFilter< PreReal , Real , TextureBitDepth >::ComputeExactSolution( bool verbose )
{
	Miscellany::Timer timer;
	solve( directSolver , multigridFilteringVariables[0].x , multigridFilteringVariables[0].rhs );
	if( verbose ) printf( "Solving time =  %.4f\n" , timer.elapsed() );
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void TextureFilter< PreReal , Real , TextureBitDepth >::UpdateSolution( bool verbose , bool detailVerbose )
{
	if( !gradientModulationUpdated )
	{
		int numTexels = (int)multigridFilteringVariables[0].rhs.size();

		Miscellany::Timer timer;
		CellStiffnessToTexelStiffness< Real , 3 >( cellModulationMask , interiorTexelToCellLines , interiorTexelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix , boundaryTexelStiffness , hierarchy.gridAtlases[0].indexConverter , texelModulatedStiffness );
		ThreadPool::ParallelFor( 0 , numTexels , [&]( unsigned int , size_t i ){ multigridFilteringVariables[0].rhs[i] = mass_x0[i]*interpolationWeight + texelModulatedStiffness[i]; } );

		if( verbose ) printf( "RHS update time %.4f\n" , timer.elapsed() );	
		gradientModulationUpdated = true;
	}

	VCycle( multigridFilteringVariables , multigridFilteringCoefficients , multigridIndices , vCycleSolvers , verbose , detailVerbose );
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void TextureFilter< PreReal , Real , TextureBitDepth >::_InitializeSystem( IndexVector< ChartIndex , IndexVector< ChartMeshTriangleIndex , SquareMatrix< PreReal , 2 > > > &parameterMetric , BoundaryProlongationData< Real > &boundaryProlongation , std::vector< Point3D< Real > > &inputSignal , std::vector< Real > &texelToCellCoeffs )
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
		default: MK_THROW( "Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles" );
		}
	}
	if( Verbose.set ) printf( "\tInitialized mass and stiffness: %.2f(s)\n" , timer.elapsed() );
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void TextureFilter< PreReal , Real , TextureBitDepth >::InitializeSystem( int width , int height )
{
	Miscellany::Timer timer;
	MultigridBlockInfo multigridBlockInfo( MultigridBlockWidth.value , MultigridBlockHeight.value , MultigridPaddedWidth.value , MultigridPaddedHeight.value );
	InitializeHierarchy( mesh , width , height , levels , textureNodes , bilinearElementIndices , hierarchy , atlasCharts , multigridBlockInfo , false );
	if( Verbose.set ) printf( "\tInitialized hierarchy: %.2f(s)\n" , timer.elapsed() );

	BoundaryProlongationData< Real > boundaryProlongation;
	InitializeBoundaryProlongationData( hierarchy.gridAtlases[0] , boundaryProlongation );

	std::vector< Point3D< Real > > low_x0( textureNodes.size() ) , high_x0( textureNodes.size() );
	for( int i=0 ; i<textureNodes.size() ; i++ ) high_x0[i] = highFrequencyTexture( textureNodes[i].ci , textureNodes[i].cj ) , low_x0[i] = lowFrequencyTexture( textureNodes[i].ci , textureNodes[i].cj );

	std::vector< Point3D< Real > > inputSignal( textureNodes.size() );
	for( int i=0 ; i<textureNodes.size() ; i++ ) inputSignal[i] = highFrequencyTexture( textureNodes[i].ci , textureNodes[i].cj );

	std::vector< Real > texelToCellCoeffs;

	timer.reset();
	IndexVector< ChartIndex , IndexVector< ChartMeshTriangleIndex , SquareMatrix< PreReal , 2 > > > parameterMetric;
	InitializeMetric( mesh , EMBEDDING_METRIC , atlasCharts , parameterMetric );
	_InitializeSystem( parameterMetric , boundaryProlongation , inputSignal , texelToCellCoeffs );

	interiorTexelToCellCoeffs.resize( 4 * static_cast< unsigned int >(hierarchy.gridAtlases[0].endInteriorTexelIndex) );
	for( unsigned int i=0 ; i<interiorTexelToCellCoeffs.size() ; i++ ) interiorTexelToCellCoeffs[i] = Point3D< Real >( Real(texelToCellCoeffs[3*i+0]) , Real(texelToCellCoeffs[3*i+1]) , Real(texelToCellCoeffs[3*i+2]) );
	
	InitializeInteriorTexelToCellLines( interiorTexelToCellLines , hierarchy.gridAtlases[0] );

	for( unsigned int c=0 ; c<3 ; c++ ) boundaryTexelStiffness[c].resize( hierarchy.gridAtlases[0].indexConverter.numBoundary() );
#ifdef NEW_CODE
	texelModulatedStiffness.resize( static_cast< unsigned int >(hierarchy.gridAtlases[0].endCombinedTexelIndex));
#else // !NEW_CODE
	texelModulatedStiffness.resize(hierarchy.gridAtlases[0].numTexels);
#endif // NEW_CODE

	if( UseDirectSolver.set )
	{
		FullMatrixConstruction( hierarchy.gridAtlases[0] , massCoefficients , mass );
		FullMatrixConstruction( hierarchy.gridAtlases[0] , stiffnessCoefficients , stiffness );
		filteringMatrix  = mass*interpolationWeight + stiffness;
	}

	multigridIndices.resize(levels);
	for( unsigned int i=0 ; i<levels ; i++ )
	{
		const typename GridAtlas<>::IndexConverter & indexConverter = hierarchy.gridAtlases[i].indexConverter;
		const GridAtlas< PreReal , Real > &gridAtlas = hierarchy.gridAtlases[i];
		multigridIndices[i].threadTasks = gridAtlas.threadTasks;
#ifdef NEW_CODE
		multigridIndices[i].boundaryToCombined = indexConverter.boundaryToCombined();
#else // !NEW_CODE
		multigridIndices[i].boundaryToSupported = indexConverter.boundaryToSupported();
#endif // NEW_CODE
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
	for( unsigned int i=0 ; i<levels ; i++ )
	{
		const typename GridAtlas<>::IndexConverter & indexConverter = hierarchy.gridAtlases[i].indexConverter;
		MultigridLevelVariables< Point3D< Real > >& variables = multigridFilteringVariables[i];
#ifdef NEW_CODE
		variables.x.resize( indexConverter.numCombined() );
		variables.rhs.resize( indexConverter.numCombined() );
		variables.residual.resize( indexConverter.numCombined() );
#else // !NEW_CODE
		variables.x.resize(hierarchy.gridAtlases[i].numTexels);
		variables.rhs.resize(hierarchy.gridAtlases[i].numTexels);
		variables.residual.resize(hierarchy.gridAtlases[i].numTexels);
#endif // NEW_CODE
		variables.boundary_rhs.resize( indexConverter.numBoundary() );
		variables.boundary_value.resize( indexConverter.numBoundary() );
		variables.variable_boundary_value.resize( indexConverter.numBoundary() );
	}

	mass_x0.resize( textureNodes.size() );
	MultiplyBySystemMatrix_NoReciprocals( massCoefficients , hierarchy.gridAtlases[0].indexConverter , hierarchy.gridAtlases[0].rasterLines , low_x0 , mass_x0 );

	stiffness_x0.resize(textureNodes.size());
	MultiplyBySystemMatrix_NoReciprocals( stiffnessCoefficients , hierarchy.gridAtlases[0].indexConverter , hierarchy.gridAtlases[0].rasterLines , high_x0 , stiffness_x0 );

	ThreadPool::ParallelFor
		(
			0 , low_x0.size() ,
			[&]( unsigned int , size_t i )
			{
				multigridFilteringVariables[0].x[i] = low_x0[i];
				multigridFilteringVariables[0].rhs[i] = mass_x0[i] * interpolationWeight + stiffness_x0[i] * gradientModulation;
			}
		);

	filteredTexture.resize(width, height);
	for (int i = 0; i < filteredTexture.size(); i++) filteredTexture[i] = Point3D< Real >( (Real)0.5 , (Real)0.5 , (Real)0.5 );


	if( UseDirectSolver.set ) ComputeExactSolution( Verbose.set );
	UpdateFilteredTexture( multigridFilteringVariables[0].x );
}

#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void TextureFilter< PreReal , Real , TextureBitDepth >::InitializeVisualization( void )
{
	sprintf( gradientModulationStr , "Gradient modulation: %.2e\n" , gradientModulation );
	sprintf( interpolationStr , "Interpolation: %.2e\n" , interpolationWeight );

	visualization.textureWidth = textureWidth;
	visualization.textureHeight = textureHeight;

	visualization.colorTextureBuffer = new unsigned char[textureHeight*textureWidth * 3];
	memset(visualization.colorTextureBuffer, 128, textureHeight * textureWidth * 3 * sizeof(unsigned char));


	unsigned int tCount = (unsigned int)mesh.numTriangles();

	visualization.triangles.resize(tCount);
	visualization.vertices.resize(3 * tCount);
	visualization.colors.resize( 3*tCount , Point3D< float >( 0.75f , 0.75f , 0.75f ) );
	visualization.textureCoordinates.resize( 3*tCount );
	visualization.normals.resize( 3*tCount );


	for( unsigned int t=0 , idx=0 ; t<tCount ; t++ )
	{
		Simplex< PreReal , 3 , 2 > sTriangle = mesh.surfaceTriangle(t);
		Simplex< PreReal , 2 , 2 > tTriangle = mesh.textureTriangle(t);
		Point3D< float > n = sTriangle.normal();
		n /= Point3D< float >::Length( n );

		for( int k=0 ; k<3 ; k++ , idx++ )
		{
			visualization.triangles[t][k] = idx;
			visualization.vertices[idx] = sTriangle[k];
			visualization.normals[idx] = n;
			visualization.textureCoordinates[idx] = tTriangle[k];
		}
	}

	std::vector< unsigned int > boundaryHalfEdges = mesh.texture.boundaryHalfEdges();

	for( int e=0 ; e<boundaryHalfEdges.size() ; e++ )
	{
		SimplexIndex< 1 > eIndex = mesh.surface.edgeIndex( boundaryHalfEdges[e] );
		for( int i=0 ; i<2 ; i++ ) visualization.chartBoundaryVertices.push_back( mesh.surface.vertices[ eIndex[i] ] );
	}

	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 's' , "export texture" , "Output Texture", ExportTextureCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'y' , "interpolation weight" , "Interpolation Weight" , InterpolationWeightCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'g' , "gradient modulation" , "Gradient Modulation" , GradientModulationCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , ' ' , "toggle update" , ToggleUpdateCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , '+' , "increment update" , IncrementUpdateCallBack ) );
	visualization.info.push_back( stepsString );

	visualization.info.push_back( gradientModulationStr );
	visualization.info.push_back( interpolationStr );

	std::string ext = ToLower( GetFileExtension( Input.values[1] ) );
	if( ext==std::string( "normap" ) || ColorAsNormal.set )
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
#endif // NO_OPEN_GL_VISUALIZATION

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void TextureFilter< PreReal , Real , TextureBitDepth >::Init( void )
{
	sprintf( stepsString , "Steps: 0" );
	levels = std::max<int>(Levels.value,1);
	interpolationWeight = InterpolationWeight.value;
	gradientModulation = GradientModulation.value;


	mesh.read( Input.values[0].c_str() , DetailVerbose.set , CollapseEpsilon.value );

	{
		std::string ext = ToLower( GetFileExtension( Input.values[1] ) );
		if( ext==std::string( "normap" ) ) ReadBinaryImage( highFrequencyTexture , Input.values[1] );
		else
		{
			ReadImage< TextureBitDepth >( highFrequencyTexture , Input.values[1] );
			if( ColorAsNormal.set ) highFrequencyTexture = ColorToNormal( highFrequencyTexture );
		}
	}
	if( InputLowFrequency.set )
	{
		std::string ext = ToLower( GetFileExtension( InputLowFrequency.value.c_str() ) );
		if( ext==std::string( "normap" ) ) ReadBinaryImage( lowFrequencyTexture , InputLowFrequency.value );
		else
		{
			ReadImage< TextureBitDepth >( lowFrequencyTexture , InputLowFrequency.value );
			if( ColorAsNormal.set ) lowFrequencyTexture = ColorToNormal( lowFrequencyTexture );
		}
	}
	else lowFrequencyTexture = highFrequencyTexture;
	if( lowFrequencyTexture.res(0)!=highFrequencyTexture.res(0) || lowFrequencyTexture.res(1)!=highFrequencyTexture.res(1) )
		MK_THROW( "Low/high texture resolutions don't match: " , lowFrequencyTexture.res(0) , " x " , lowFrequencyTexture.res(1) , " != " , highFrequencyTexture.res(0) , " x " , highFrequencyTexture.res(1) );

	textureWidth = highFrequencyTexture.res(0);
	textureHeight = highFrequencyTexture.res(1);

	// Define centroid and scale for visualization
	Point3D< PreReal > centroid = mesh.surface.centroid();
	PreReal radius = mesh.surface.boundingRadius( centroid );
	for( int i=0 ; i<mesh.surface.vertices.size() ; i++ ) mesh.surface.vertices[i] = ( mesh.surface.vertices[i] - centroid ) / radius;

	// Apply a random jitter to the texture coordinates
	if( RandomJitter.set )
	{
		if( RandomJitter.value ) srand( RandomJitter.value );
		else                     srand( (unsigned int)time(NULL) );
		PreReal jitterScale = (PreReal)1e-3 / std::max< int >( textureWidth , textureHeight );
		for( int i=0 ; i<mesh.texture.vertices.size() ; i++ ) mesh.texture.vertices[i] += Point2D< PreReal >( (PreReal)1. - Random< PreReal >()*2 , (PreReal)1. - Random< PreReal >()*2 ) * jitterScale;
	}

	// Pad the texture and texture coordinates
	{
		padding = Padding::Init( textureWidth , textureHeight , mesh.texture.vertices , DetailVerbose.set );
		padding.pad( textureWidth , textureHeight , mesh.texture.vertices );
		padding.pad( highFrequencyTexture );
		padding.pad( lowFrequencyTexture );
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

			rMesh.exp(xForms, _p);
			
			textureNodes[i].tID = AtlasMeshTriangleIndex( _p.tIdx );
			textureNodes[i].barycentricCoords = _p.p;
		}
	}

	textureNodePositions.resize( textureNodes.size() );
	for( int i=0 ; i<textureNodePositions.size() ; i++ ) textureNodePositions[i] = mesh.surface( textureNodes[i] );

	uniformTexelModulationMask.resize( textureNodes.size() , 0.5 );

	for (int c = 0; c < 3; c++) texelStiffness[c].resize(textureNodes.size());

	cellModulationMask.resize( bilinearElementIndices.size() , 1 );
	uniformCellModulationMask.resize( bilinearElementIndices.size() , 0.5 );

	cellCenterPositions.resize( bilinearElementIndices.size() );
	for( int i=0 ; i<bilinearElementIndices.size() ; i++ )
		cellCenterPositions[i] = Point3D< float >
		(
#ifdef NEW_CODE
			textureNodePositions[ static_cast< unsigned int >( bilinearElementIndices[ AtlasCombinedCellIndex(i) ][0] ) ] +
			textureNodePositions[ static_cast< unsigned int >( bilinearElementIndices[ AtlasCombinedCellIndex(i) ][1] ) ] +
			textureNodePositions[ static_cast< unsigned int >( bilinearElementIndices[ AtlasCombinedCellIndex(i) ][2] ) ] +
			textureNodePositions[ static_cast< unsigned int >( bilinearElementIndices[ AtlasCombinedCellIndex(i) ][3] ) ]
#else // !NEW_CODE
			textureNodePositions[ bilinearElementIndices[ AtlasCombinedCellIndex(i) ][0] ] +
			textureNodePositions[ bilinearElementIndices[ AtlasCombinedCellIndex(i) ][1] ] +
			textureNodePositions[ bilinearElementIndices[ AtlasCombinedCellIndex(i) ][2] ] +
			textureNodePositions[ bilinearElementIndices[ AtlasCombinedCellIndex(i) ][3] ]
#endif // NEW_CODE
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
		if( multiChartTexelCount ) MK_WARN( "Non-zero multi-chart texels: " , multiChartTexelCount );
	}
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void _main( int argc , char *argv[] )
{
	TextureFilter< PreReal , Real , TextureBitDepth >::Init();

	TextureFilter< PreReal , Real , TextureBitDepth >::updateCount = Paused.set ?  0 : -1;

#ifdef NO_OPEN_GL_VISUALIZATION
	if( UseDirectSolver.set ) TextureFilter< PreReal , Real , TextureBitDepth >::ComputeExactSolution( Verbose.set );
	else for ( int i=0 ; i<OutputVCycles.value ; i++ ) TextureFilter< PreReal , Real , TextureBitDepth >::UpdateSolution();
	TextureFilter< PreReal , Real , TextureBitDepth >::WriteTexture( Output.value.c_str() );
#else // !NO_OPEN_GL_VISUALIZATION
	if( !Output.set )
	{
		glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
		TextureFilter< PreReal , Real , TextureBitDepth >::visualization.displayMode = DisplayMode.value;
		if     ( DisplayMode.value==ONE_REGION_DISPLAY   ) TextureFilter< PreReal , Real , TextureBitDepth >::visualization.screenWidth =  800 , TextureFilter< PreReal , Real , TextureBitDepth >::visualization.screenHeight = 800;
		else if( DisplayMode.value==TWO_REGION_DISPLAY   ) TextureFilter< PreReal , Real , TextureBitDepth >::visualization.screenWidth = 1600 , TextureFilter< PreReal , Real , TextureBitDepth >::visualization.screenHeight = 800;
		else if( DisplayMode.value==THREE_REGION_DISPLAY ) TextureFilter< PreReal , Real , TextureBitDepth >::visualization.screenWidth = 1200 , TextureFilter< PreReal , Real , TextureBitDepth >::visualization.screenHeight = 800;
		else if( DisplayMode.value==FOUR_REGION_DISPLAY  ) TextureFilter< PreReal , Real , TextureBitDepth >::visualization.screenWidth = 1500 , TextureFilter< PreReal , Real , TextureBitDepth >::visualization.screenHeight = 600;
		TextureFilter< PreReal , Real , TextureBitDepth >::visualization.useNearestSampling = Nearest.set;
		if( Snapshot.set ) TextureFilter< PreReal , Real , TextureBitDepth >::visualization.setSnapshot( Snapshot.value.c_str() , true );
		glutInitWindowSize( TextureFilter< PreReal , Real , TextureBitDepth >::visualization.screenWidth , TextureFilter< PreReal , Real , TextureBitDepth >::visualization.screenHeight );
		glutInit( &argc , argv );
		char windowName[1024];
		sprintf( windowName , "Texture Filtering" );
		glutCreateWindow( windowName );
		if( glewInit()!=GLEW_OK ) MK_THROW( "glewInit failed" );
		glutDisplayFunc ( TextureFilter< PreReal , Real , TextureBitDepth >::Display );
		glutReshapeFunc ( TextureFilter< PreReal , Real , TextureBitDepth >::Reshape );
		glutMouseFunc   ( TextureFilter< PreReal , Real , TextureBitDepth >::MouseFunc );
		glutMotionFunc  ( TextureFilter< PreReal , Real , TextureBitDepth >::MotionFunc );
		glutKeyboardFunc( TextureFilter< PreReal , Real , TextureBitDepth >::KeyboardFunc) ;
		glutIdleFunc    ( TextureFilter< PreReal , Real , TextureBitDepth >::Idle );
		if( CameraConfig.set ) TextureFilter< PreReal , Real , TextureBitDepth >::visualization.ReadSceneConfigurationCallBack( &TextureFilter< PreReal , Real , TextureBitDepth >::visualization , CameraConfig.value.c_str() );
		TextureFilter< PreReal , Real , TextureBitDepth >::InitializeVisualization();
		TextureFilter< PreReal , Real , TextureBitDepth >::visualization.showSlideBar = true;
		glutMainLoop(); 
	}
	else
	{
		if( UseDirectSolver.set ) TextureFilter< PreReal , Real , TextureBitDepth >::ComputeExactSolution( Verbose.set );
		else for ( int i=0 ; i<OutputVCycles.value ; i++ ) TextureFilter< PreReal , Real , TextureBitDepth >::UpdateSolution();
		TextureFilter< PreReal , Real , TextureBitDepth >::ExportTextureCallBack( &TextureFilter< PreReal , Real , TextureBitDepth >::visualization , Output.value.c_str() );
	}
#endif // NO_OPEN_GL_VISUALIZATION

}

template< typename PreReal , typename Real >
void _main( int argc , char *argv[] , unsigned int bitDepth )
{
	switch( bitDepth )
	{
	case (unsigned int)-1:
	case  8: return _main< PreReal , Real ,  8 >( argc , argv );
	case 16: return _main< PreReal , Real , 16 >( argc , argv );
	case 32: return _main< PreReal , Real , 32 >( argc , argv );
	case 64: return _main< PreReal , Real , 64 >( argc , argv );
	default: MK_THROW( "Only bit depths of 8, 16, 32, and 64 supported: " , bitDepth );
	}
}

int main(int argc, char* argv[])
{
	CmdLineParse( argc-1 , argv+1 , params );
#ifdef NO_OPEN_GL_VISUALIZATION
	if( !Input.set || !Output.set ) { ShowUsage( argv[0] ) ; return EXIT_FAILURE; }
#else // !NO_OPEN_GL_VISUALIZATION
	if( !Input.set ) { ShowUsage( argv[0] ) ; return EXIT_FAILURE; }
#endif // NO_OPEN_GL_VISUALIZATION

	unsigned int bitDepth = -1;
	std::string ext = ToLower( GetFileExtension( Input.values[1] ) );
	if( ext!=std::string( "normap" ) )
	{
		unsigned int width , height , channels;
		ImageReader< 8 >::GetInfo( Input.values[1].c_str() , width , height , channels , bitDepth );
	}

	if( Serial.set ) ThreadPool::ParallelizationType = ThreadPool::ParallelType::NONE;
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
