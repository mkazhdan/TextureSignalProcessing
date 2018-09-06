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
#include <Src/SimpleMesh.h>
#include <Src/Basis.h>
#include <Misha/FEM.h>
#include <Src/Solver.h>
#include <Src/Hierarchy.h>
#include <Src/MassAndStiffness.h>
#include <Src/InteriorTexelToCellLines.inl>
#include <Src/Padding.h>
#include <Src/TextureFilteringVisualization.h>

const std::string vertex_shader_src =
#include <Shaders/normal_texture_vertex.vs>
;
const std::string fragment_shader_src =
#include <Shaders/normal_texture_fragment.fs>
;

cmdLineParameterArray< char* , 2 > Input( "in" );
cmdLineParameter< char* > Output( "out" );
cmdLineParameter< int   > OutputVCycles( "outVCycles" , 6 );
cmdLineParameter< float > InterpolationWeight( "interpolation" , 1e3 );
cmdLineParameter< float > GradientModulation( "modulation" , 1.0 );
cmdLineParameter< int   > DisplayMode("display",FOUR_REGION_DISPLAY);
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
	&Input , &Output , &InterpolationWeight , &GradientModulation , &CameraConfig , &Levels , &UseDirectSolver , &Threads , &DisplayMode , &Verbose ,
	&DetailVerbose , &MultigridBlockHeight , &MultigridBlockWidth , &MultigridPaddedHeight , &MultigridPaddedWidth , &RandomJitter ,
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
	printf( "\t[--%s <output texture>]\n" , Output.name );
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
	printf( "\t[--%s <display mode>=%d]\n" , DisplayMode.name , DisplayMode.value );
	printf( "\t\t%d] One Region \n"   , ONE_REGION_DISPLAY   );
	printf( "\t\t%d] Two Region \n"   , TWO_REGION_DISPLAY   );
	printf( "\t\t%d] Three Region \n" , THREE_REGION_DISPLAY );
	printf( "\t\t%d] Four Region \n"  , FOUR_REGION_DISPLAY  );
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

template<class Real>
class TextureFilter
{
public:
	static TexturedMesh mesh;
	static int textureWidth;
	static int textureHeight;
	static double interpolationWeight;
	static double gradientModulation;
	static int levels;

	static HierarchicalSystem hierarchy;
	static bool gradientModulationUpdated;
	static bool positiveModulation;

	static Image<Point3D<float>> filteredTexture;
	//UI
	static char gradientModulationStr[1024];
	static char interpolationStr[1024];

	static std::vector<Real> uniformTexelModulationMask;
	static std::vector<Real> cellModulationMask;
	static std::vector<Real> uniformCellModulationMask;
	static std::vector<Real> texelStiffness[3];
	static std::vector<Point3D<float>> cellCenterPositions;

	static std::vector<Point3D<float>>textureNodePositions;

	static std::vector<AtlasChart> atlasCharts;

	static std::vector< BilinearElementIndex > bilinearElementIndices;
	
	static std::vector<TextureNodeInfo> textureNodes;

	static SparseMatrix<double,int> mass;
	static SparseMatrix<double,int> stiffness;
	static SparseMatrix< double, int > filteringMatrix;
	//RHS computation
	static std::vector<InteriorTexelToCellLine> interiorTexelToCellLines;
	static std::vector< Point3D< Real > > interiorTexelToCellCoeffs;
	static SparseMatrix<Real, int> boundaryCellBasedStiffnessRHSMatrix[3];
	static std::vector<Real> boundaryTexelStiffness[3];
	static std::vector< Point3D< Real > > texelModulatedStiffness;

	static std::vector< Point3D< Real > > mass_x0;
	static std::vector< Point3D< Real > > stiffness_x0;

	static std::vector< MultigridLevelCoefficients< Real > > multigridFilteringCoefficients;
	static std::vector< MultigridLevelVariables< Point3D< Real > > > multigridFilteringVariables;
	static std::vector< MultigridLevelIndices< Real > > multigridIndices;

#if defined( USE_CHOLMOD )
	typedef  std::vector< CholmodCholeskySolver< Real , 3 > > BoundarySolverType;
	typedef  CholmodCholeskySolver< Real , 3 > CoarseSolverType;
	typedef  CholmodCholeskySolver< Real , 3 > DirectSolverType;
#elif defined( USE_EIGEN_SIMPLICIAL )
	typedef  std::vector< EigenCholeskySolver< Real , 3 > > BoundarySolverType;
	typedef  EigenCholeskySolver< Real , 3 > CoarseSolverType;
	typedef  EigenCholeskySolver< Real , 3 > DirectSolverType;
#elif defined( USE_EIGEN_PARDISO )
	typedef  std::vector< EigenPardisoSolver< Real , 3 > > BoundarySolverType;
	typedef  EigenPardisoSolver< Real , 3 > CoarseSolverType;
	typedef  EigenPardisoSolver< Real , 3 > DirectSolverType;
#else
#error "[ERROR] No solver defined!"
#endif

	static BoundarySolverType boundarySolver;
	static CoarseSolverType coarseSolver;
	static DirectSolverType directSolver;


	//Linear Operators
	static std::vector<double> deepMassCoefficients;
	static std::vector<double> deepStiffnessCoefficients;
	static SparseMatrix<double, int> boundaryBoundaryMassMatrix;
	static SparseMatrix<double, int> boundaryBoundaryStiffnessMatrix;
	static SparseMatrix<double, int> boundaryDeepMassMatrix;
	static SparseMatrix<double, int> boundaryDeepStiffnessMatrix;

	static int steps;
	static char stepsString[];

	static Padding padding;

	//Visulization
	static TextureFilteringVisualization visualization;
	static int updateCount;

	static void ToggleUpdateCallBack( Visualization* v , const char* prompt );
	static void IncrementUpdateCallBack( Visualization* v , const char* prompt );
	static void ExportTextureCallBack(Visualization* v, const char* prompt);

	static void GradientModulationCallBack( Visualization* v , const char* prompt );
	static void InterpolationWeightCallBack( Visualization* v , const char* prompt );

	static int Init();
	static void InitializeVisualization();
	static int UpdateSolution(bool verbose = false, bool detailVerbose = false);
	static void ComputeExactSolution( bool verbose=false );
	static int InitializeSystem(const int width, const int height);
	static int _InitializeSystem( std::vector<std::vector<SquareMatrix<double, 2>>>& parameterMetric , BoundaryProlongationData& boundaryProlongation , std::vector<Point3D<double>>& inputSignal , std::vector<double>& texelToCellCoeffs );

	static void UpdateFilteredColorTexture( const std::vector< Point3D< Real > >& solution );
	static void UpdateFilteredTexture( const std::vector< Point3D< Real > >& solution );

	static void UpdateMaskTexture();


	static void Display( void ){ visualization.Display(); }
	static void MouseFunc( int button , int state , int x , int y );
	static void MotionFunc( int x , int y );
	static void Reshape( int w , int h ){ visualization.Reshape(w,h); }
	static void KeyboardFunc(unsigned char key, int x, int y){ visualization.KeyboardFunc( key , x , y ); }
	static void Idle( void );
};

template< class Real > char TextureFilter< Real >::gradientModulationStr[1024];
template< class Real > char TextureFilter< Real >::interpolationStr[1024];

template<class Real> TexturedMesh												TextureFilter<Real>::mesh;
template<class Real> int														TextureFilter<Real>::textureWidth;
template<class Real> int														TextureFilter<Real>::textureHeight;
template<class Real> TextureFilteringVisualization								TextureFilter<Real>::visualization;
template<class Real> SparseMatrix<double,int>									TextureFilter<Real>::mass;
template<class Real> SparseMatrix<double,int>									TextureFilter<Real>::stiffness;
template<class Real> SparseMatrix<double, int>									TextureFilter<Real>::filteringMatrix;

template<class Real> double														TextureFilter<Real>::interpolationWeight;
template<class Real> double														TextureFilter<Real>::gradientModulation;

template<class Real> std::vector<TextureNodeInfo>								TextureFilter<Real>::textureNodes;
template<class Real> std::vector< BilinearElementIndex >						TextureFilter<Real>::bilinearElementIndices;

template< class Real > int														TextureFilter< Real >::steps;
template< class Real > char														TextureFilter< Real >::stepsString[1024];
template<class Real> int														TextureFilter<Real>::levels;
template<class Real> HierarchicalSystem											TextureFilter<Real>::hierarchy;


template<class Real> bool TextureFilter<Real>::gradientModulationUpdated = true;
template<class Real> bool TextureFilter<Real>::positiveModulation = true;

//UI
template<class Real> std::vector<Point3D<float>>										TextureFilter<Real>::cellCenterPositions;
template<class Real> std::vector<Real>													TextureFilter<Real>::texelStiffness[3];
template<class Real> std::vector<Real>													TextureFilter<Real>::cellModulationMask;
template<class Real> std::vector<Real>													TextureFilter<Real>::uniformCellModulationMask;

template<class Real> Image<Point3D<float>>												TextureFilter<Real>::filteredTexture;
template<class Real> std::vector< Point3D< Real > >										TextureFilter<Real>::stiffness_x0;
template<class Real> std::vector< Point3D< Real > >										TextureFilter<Real>::mass_x0;

template<class Real> std::vector<MultigridLevelCoefficients<Real>>						TextureFilter<Real>::multigridFilteringCoefficients;
template<class Real> std::vector< MultigridLevelVariables< Point3D< Real > > >			TextureFilter<Real>::multigridFilteringVariables;
template<class Real> std::vector<MultigridLevelIndices<Real>>							TextureFilter<Real>::multigridIndices;


template<class Real> typename TextureFilter<Real>::BoundarySolverType					TextureFilter<Real>::boundarySolver;
template<class Real> typename TextureFilter<Real>::CoarseSolverType						TextureFilter<Real>::coarseSolver;
template<class Real> typename TextureFilter<Real>::DirectSolverType						TextureFilter<Real>::directSolver;

template<class Real> std::vector<AtlasChart>											TextureFilter<Real>::atlasCharts;

template<class Real> std::vector<InteriorTexelToCellLine>								TextureFilter<Real>::interiorTexelToCellLines;
template<class Real> std::vector< Point3D< Real > >										TextureFilter<Real>::interiorTexelToCellCoeffs;
template<class Real> SparseMatrix<Real, int>											TextureFilter<Real>::boundaryCellBasedStiffnessRHSMatrix[3];
template<class Real> std::vector<Real>													TextureFilter<Real>::boundaryTexelStiffness[3];
template<class Real> std::vector< Point3D< Real > >										TextureFilter<Real>::texelModulatedStiffness;

template<class Real> std::vector<Point3D<float>>										TextureFilter<Real>::textureNodePositions;
template<class Real> std::vector<Real>													TextureFilter<Real>::uniformTexelModulationMask;

template<class Real> Padding															TextureFilter<Real>::padding;


template<class Real>  std::vector<double>												TextureFilter<Real>::deepMassCoefficients;
template<class Real>  std::vector<double>												TextureFilter<Real>::deepStiffnessCoefficients;
template<class Real>  SparseMatrix<double, int>											TextureFilter<Real>::boundaryBoundaryMassMatrix;
template<class Real>  SparseMatrix<double, int>											TextureFilter<Real>::boundaryBoundaryStiffnessMatrix;
template<class Real>  SparseMatrix<double, int>											TextureFilter<Real>::boundaryDeepMassMatrix;
template<class Real>  SparseMatrix<double, int>											TextureFilter<Real>::boundaryDeepStiffnessMatrix;
template< class Real > int																TextureFilter< Real >::updateCount = -1;

template<class Real>
void TextureFilter<Real>::UpdateMaskTexture(){

#pragma omp parallel for
	for (int i = 0; i < textureNodes.size(); i++){
		float texelModulationValue = uniformTexelModulationMask[i];
		if (texelModulationValue != 0.5){
			Point3D<float> newColor;
			if (texelModulationValue > 0.5) {
				texelModulationValue = 2.0 * texelModulationValue - 1.0;
				newColor = Point3D<float>(1.f, 0.f, 0.f) * texelModulationValue + Point3D<float>(0.8f, 0.8f, 0.8f) * (1.0 - texelModulationValue);
			}
			else {
				texelModulationValue = 2.0 *texelModulationValue;
				newColor = Point3D<float>(0.f, 0.f, 1.f) *(1.0 - texelModulationValue) + Point3D<float>(0.8f, 0.8f, 0.8f) * texelModulationValue;
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

template<class Real>
void TextureFilter< Real >::UpdateFilteredColorTexture( const std::vector< Point3D< Real > > & solution )
{
#pragma omp parallel for
	for (int i = 0; i < textureNodes.size(); i++) {
		int ci = textureNodes[i].ci;
		int cj = textureNodes[i].cj;
		int offset = 3 * (textureWidth*cj + ci);
		for (int c = 0; c < 3; c++) {
			double value = std::min<double>(1.0, std::max<double>(0, solution[i][c]));
			visualization.colorTextureBuffer[offset + c] = (unsigned char)(value*255.0);
		}
	}
}

template< class Real >
void TextureFilter< Real >::UpdateFilteredTexture( const std::vector< Point3D< Real > >& solution )
{
#pragma omp parallel for
	for( int i=0 ; i<textureNodes.size() ; i++ )
	{
		int ci = textureNodes[i].ci , cj = textureNodes[i].cj;
		filteredTexture(ci,cj) = Point3D< float >( solution[i][0] , solution[i][1] , solution[i][2] );
	}
}

template< class Real >
void TextureFilter< Real >::Idle( void )
{
	if( visualization.promptCallBack ) visualization.showSlideBar = false;
	else                               visualization.showSlideBar = true;

	auto RescaleFunction = []( double x )
	{
		//return x*(2.0 *x + 1);
		return 2.0 *x;
	};

	float radius = 0.1;
	float modulationVartiation = 0.2;
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
				float distanceRatio = Point3D< float >::Length(cellCenterPositions[i] - selectedPoint) / radius;
				float factor = 1.0 - distanceRatio;
				factor = factor < 0 ? 0 : factor*factor*(-2.0*factor + 3.0);
				Real uniformModulationMaskValue = std::max< Real >( 0 , std::min< Real >( 1.0 , uniformCellModulationMask[i] + modulationVartiation * factor ) );
				uniformCellModulationMask[i] = uniformModulationMaskValue;
				cellModulationMask[i] = RescaleFunction( uniformModulationMaskValue );
			}
			if( true )
			{
				for( int i=0 ; i<textureNodePositions.size() ; i++ )
				{
					float distanceRatio = Point3D< float >::Length( textureNodePositions[i] - selectedPoint ) / radius;
					float factor = 1.0 - distanceRatio;
					factor = factor < 0 ? 0 : factor*factor*(-2.0*factor + 3.0);
					Real modulationMaskValue = std::max< Real >( 0 , std::min< Real >( 1.0 , uniformTexelModulationMask[i] + modulationVartiation * factor ) );
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
			Real diff = (Real)(visualization.slideBarCursorPosition - visualization.slideBarCursorOldPosition);
			visualization.slideBarCursorOldPosition = visualization.slideBarCursorPosition;

			for( int i=0 ; i<bilinearElementIndices.size() ; i++ )
			{
				Real uniformModulationMaskValue = std::max< Real >( 0 , std::min< Real >( 1.0 , uniformCellModulationMask[i] + diff ) );
				uniformCellModulationMask[i] = uniformModulationMaskValue;
				cellModulationMask[i] = RescaleFunction( uniformModulationMaskValue );
			}

			if( true )
			{
				for( int i=0 ; i<textureNodePositions.size() ; i++ )
				{
					Real modulationMaskValue = std::max< Real >( 0 , std::min< Real >( 1.0 , uniformTexelModulationMask[i] + diff ) );
					uniformTexelModulationMask[i] = modulationMaskValue;
				}
				UpdateMaskTexture();
			}
			steps = 0;
		}
	}

	if( updateCount && !UseDirectSolver.set && !visualization.promptCallBack )
	{
		if( !UpdateSolution() ) fprintf( stderr , "[ERROR] Updated solution failed!\n" );

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

template< class Real >
void TextureFilter< Real >::MouseFunc( int button , int state , int x , int y )
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
		double slideBarCursorPosition = double(x - 20) / double(visualization.slideBarWidth() - 40);
		slideBarCursorPosition = std::min<double>(std::max<double>(slideBarCursorPosition, 0), 1.0);
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

template<class Real>
void TextureFilter<Real>::MotionFunc( int x , int y )
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
		double slideBarCursorPosition = double(x - 20) / double(visualization.slideBarWidth() - 40);
		slideBarCursorPosition = std::min<double>(std::max<double>(slideBarCursorPosition, 0), 1.0);
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

			int imageSize = std::min< int >(visualization.screenWidth, visualization.screenHeight);
			if (visualization.panning) visualization.xForm.offset[0] -= (visualization.newX - visualization.oldX) / visualization.imageToScreenScale(), visualization.xForm.offset[1] += (visualization.newY - visualization.oldY) / visualization.imageToScreenScale();
			else
			{
				float dz = float(pow(1.1, double(visualization.newY - visualization.oldY) / 8));
				visualization.xForm.zoom *= dz;
			}
		}
	}
	glutPostRedisplay();
}

template< class Real > void TextureFilter< Real >::ToggleUpdateCallBack( Visualization* v , const char* prompt )
{
	if( updateCount ) updateCount = 0;
	else              updateCount = -1;
}
template< class Real > void TextureFilter< Real >::IncrementUpdateCallBack( Visualization* v , const char* prompt )
{
	if( updateCount<0 ) updateCount = 1;
	else updateCount++;
}

template<class Real>
void TextureFilter<Real>::ExportTextureCallBack( Visualization* v , const char* prompt )
{
	UpdateFilteredTexture( multigridFilteringVariables[0].x );
	Image< Point3D< float > > outputTexture = filteredTexture;
	if( padding.nonTrivial ) UnpadImage( padding , outputTexture );

	char* ext = GetFileExtension( prompt );
	if( !strcasecmp( ext , "normap" ) ) WriteBinaryImage( outputTexture , prompt );
	else
	{
		if( visualization.textureType==NORMAL_TEXTURE )
		{
			for( int i=0 ; i<outputTexture.size() ; i++ )
			{
				outputTexture[i] /= Point3D<float>::Length( outputTexture[i] );
				outputTexture[i] = outputTexture[i] * 0.5f + Point3D< float >( 0.5f , 0.5f , 0.5f );
			}
		}
		outputTexture.write( prompt );
	}
	delete[] ext;
}

template< class Real >
void  TextureFilter< Real >::GradientModulationCallBack( Visualization* v , const char* prompt )
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

template< class Real >
void  TextureFilter< Real >::InterpolationWeightCallBack( Visualization* v , const char* prompt )
{
	interpolationWeight = atof(prompt);
	if( UseDirectSolver.set ) filteringMatrix = mass*interpolationWeight + stiffness;
	clock_t t_begin = clock();
	if( !UpdateLinearSystem( interpolationWeight , 1.0 , hierarchy , multigridFilteringCoefficients ,
		deepMassCoefficients , deepStiffnessCoefficients ,
		boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix ,
		boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix ,
		coarseSolver , boundarySolver , directSolver ,
		filteringMatrix , DetailVerbose.set , false , UseDirectSolver.set ) ) {
		fprintf( stderr , "[ERROR] Failed system update!\n" );
	}
	if( Verbose.set ) printf( "\tInitialized multigrid coefficients: %.2f(s)\n" , double(clock() - t_begin) / CLOCKS_PER_SEC);

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


template<class Real>
void TextureFilter<Real>::ComputeExactSolution( bool verbose )
{
	clock_t t_begin = clock();
	solve( directSolver , multigridFilteringVariables[0].x , multigridFilteringVariables[0].rhs );
	if( verbose ) printf( "Solving time =  %.4f\n" , double(clock() - t_begin) / CLOCKS_PER_SEC );
}

template<class Real>
int TextureFilter<Real>::UpdateSolution( bool verbose , bool detailVerbose )
{
	if( !gradientModulationUpdated )
	{
		int numTexels = (int)multigridFilteringVariables[0].rhs.size();
		clock_t p_begin = clock();

		CellStiffnessToTexelStiffness< Real , 3 >(cellModulationMask, interiorTexelToCellLines, interiorTexelToCellCoeffs, boundaryCellBasedStiffnessRHSMatrix, boundaryTexelStiffness, hierarchy.gridAtlases[0].boundaryGlobalIndex, texelModulatedStiffness);
#pragma omp parallel for
		for( int i=0 ; i<numTexels ; i++ ) multigridFilteringVariables[0].rhs[i] = mass_x0[i]*interpolationWeight + texelModulatedStiffness[i];

		if( verbose ) printf("RHS update time %.4f  \n", double(clock() - p_begin) / CLOCKS_PER_SEC);	
		gradientModulationUpdated = true;
	}

	VCycle( multigridFilteringVariables , multigridFilteringCoefficients , multigridIndices , boundarySolver , coarseSolver , verbose , detailVerbose );

	return 1;
}

template<>
int TextureFilter< float >::_InitializeSystem( std::vector<std::vector<SquareMatrix<double, 2>>>& parameterMetric , BoundaryProlongationData& boundaryProlongation , std::vector<Point3D<double>>& inputSignal , std::vector<double>& texelToCellCoeffs )
{
	SparseMatrix<double, int> _boundaryCellBasedStiffnessRHSMatrix[3];
	clock_t t_begin = clock();
	{
		int ret = 0;
		switch( MatrixQuadrature.value )
		{
		case 1:
			ret = InitializeMassAndStiffness< 1>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts,boundaryProlongation , true , inputSignal , texelToCellCoeffs , _boundaryCellBasedStiffnessRHSMatrix );
			break;
		case 3:
			ret = InitializeMassAndStiffness< 3>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts,boundaryProlongation , true , inputSignal , texelToCellCoeffs , _boundaryCellBasedStiffnessRHSMatrix );
			break;
		case 6:
			ret = InitializeMassAndStiffness< 6>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts,boundaryProlongation , true , inputSignal , texelToCellCoeffs , _boundaryCellBasedStiffnessRHSMatrix );
			break;
		case 12:
			ret = InitializeMassAndStiffness<12>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts,boundaryProlongation , true , inputSignal , texelToCellCoeffs , _boundaryCellBasedStiffnessRHSMatrix );
			break;
		case 24:
			ret = InitializeMassAndStiffness<24>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts,boundaryProlongation , true , inputSignal , texelToCellCoeffs , _boundaryCellBasedStiffnessRHSMatrix );
			break;
		case 32:
			ret = InitializeMassAndStiffness<32>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts,boundaryProlongation , true , inputSignal , texelToCellCoeffs , _boundaryCellBasedStiffnessRHSMatrix );
			break;
		default:
			fprintf( stderr , "[ERROR] Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles\n" );
		}
		if( !ret )
		{
			fprintf( stderr , "[ERROR] Failed intialization!\n" );
			return 0;
		}
	}
	if( Verbose.set ) printf( "\tInitialized mass and stiffness: %.2f(s)\n" , double(clock() - t_begin) / CLOCKS_PER_SEC);
	for( int c=0 ; c<3 ; c++ ) boundaryCellBasedStiffnessRHSMatrix[c] = _boundaryCellBasedStiffnessRHSMatrix[c];
	return 1;
}
template<>
int TextureFilter< double >::_InitializeSystem( std::vector<std::vector<SquareMatrix<double, 2>>>& parameterMetric , BoundaryProlongationData& boundaryProlongation , std::vector<Point3D<double>>& inputSignal , std::vector<double>& texelToCellCoeffs )
{
	clock_t t_begin = clock();
	{
		int ret = 0;
		switch( MatrixQuadrature.value )
		{
		case 1:
			InitializeMassAndStiffness< 1>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , true , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix );
			break;
		case 3:
			InitializeMassAndStiffness< 3>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , true , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix );
			break;
		case 6:
			InitializeMassAndStiffness< 6>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , true , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix );
			break;
		case 12:
			InitializeMassAndStiffness<12>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , true , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix );
			break;
		case 24:
			InitializeMassAndStiffness<24>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , true , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix );
			break;
		case 32:
			InitializeMassAndStiffness<32>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , true , inputSignal , texelToCellCoeffs , boundaryCellBasedStiffnessRHSMatrix );
			break;
		default:
			fprintf( stderr , "[ERROR] Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles\n" );
		}
		if( !ret )
		{
			fprintf( stderr , "[ERROR] Failed intialization!\n" );
			return 0;
		}
	}
	if( Verbose.set ) printf( "\tInitialized mass and stiffness: %.2f(s)\n" , double(clock() - t_begin) / CLOCKS_PER_SEC);
	return 1;
}
template<class Real>
int TextureFilter<Real>::InitializeSystem( const int width , const int height )
{
	clock_t t_begin;

	t_begin = clock();
	MultigridBlockInfo multigridBlockInfo(MultigridBlockWidth.value, MultigridBlockHeight.value,MultigridPaddedWidth.value,MultigridPaddedHeight.value, 0);
	if( !InitializeHierarchy( mesh , width , height , levels , textureNodes , bilinearElementIndices , hierarchy , atlasCharts , multigridBlockInfo , true , DetailVerbose.set ) )
	{
		printf("ERROR : Failed intialization! \n");
		return 0;
	}
	if( Verbose.set ) printf( "\tInitialized hierarchy: %.2f(s)\n" , double(clock() - t_begin) / CLOCKS_PER_SEC);

	BoundaryProlongationData boundaryProlongation;
	if (!InitializeBoundaryProlongationData(hierarchy.gridAtlases[0], boundaryProlongation)){
		printf("ERROR : Failed boundary prolongation! \n");
		return 0;
	}

	std::vector< Point3D< Real > > _x0;
	_x0.resize(textureNodes.size());
	
	for (int i = 0; i < textureNodes.size(); i++) {
		Point3D<float> texelValue = mesh.texture(textureNodes[i].ci, textureNodes[i].cj);
		_x0[i] = Point3D< Real >( texelValue[0] , texelValue[1] , texelValue[2] );
	}

	std::vector<Point3D<double>> inputSignal(textureNodes.size());
	for (int i = 0; i < textureNodes.size(); i++) inputSignal[i] = mesh.texture(textureNodes[i].ci, textureNodes[i].cj);
	std::vector<double> texelToCellCoeffs;


	t_begin = clock();

	std::vector<std::vector<SquareMatrix<double, 2>>> parameterMetric;
	if (!InitializeMetric(mesh, EMBEDDING_METRIC, atlasCharts, parameterMetric)) {
		printf("ERROR: Unable to initialize metric \n");
		return 0;
	}

	if( !_InitializeSystem( parameterMetric , boundaryProlongation , inputSignal , texelToCellCoeffs ) ) return 0;

	interiorTexelToCellCoeffs.resize(4 * hierarchy.gridAtlases[0].numDeepTexels);
	for( int i=0 ; i<4*hierarchy.gridAtlases[0].numDeepTexels ; i++ ) interiorTexelToCellCoeffs[i] = Point3D< Real >( Real(texelToCellCoeffs[3*i+0]) , Real(texelToCellCoeffs[3*i+1]) , Real(texelToCellCoeffs[3*i+2]) );
	
	if (!InitializeInteriorTexelToCellLines(interiorTexelToCellLines, hierarchy.gridAtlases[0])) {
		printf("ERROR: Interior texel to cell not initialized! \n");
		return 0;
	}

	for (int c = 0; c < 3; c++) boundaryTexelStiffness[c].resize(hierarchy.gridAtlases[0].boundaryGlobalIndex.size());
	texelModulatedStiffness.resize(hierarchy.gridAtlases[0].numTexels);

	if( UseDirectSolver.set )
	{
		FullMatrixConstruction(hierarchy.gridAtlases[0], deepMassCoefficients, boundaryBoundaryMassMatrix, boundaryDeepMassMatrix, mass);
		FullMatrixConstruction(hierarchy.gridAtlases[0], deepStiffnessCoefficients, boundaryBoundaryStiffnessMatrix, boundaryDeepStiffnessMatrix, stiffness);
		filteringMatrix  = mass*interpolationWeight + stiffness;
	}

	multigridIndices.resize(levels);
	for( int i=0 ; i<levels ; i++ )
	{
		const GridAtlas & gridAtlas = hierarchy.gridAtlases[i];
		multigridIndices[i].threadTasks = gridAtlas.threadTasks;
		multigridIndices[i].boundaryGlobalIndex = gridAtlas.boundaryGlobalIndex;
		multigridIndices[i].segmentedLines = gridAtlas.segmentedLines;
		multigridIndices[i].rasterLines = gridAtlas.rasterLines;
		multigridIndices[i].restrictionLines = gridAtlas.restrictionLines;
		multigridIndices[i].prolongationLines = gridAtlas.prolongationLines;
		if( i<levels-1 ) multigridIndices[i].boundaryRestriction = hierarchy.boundaryRestriction[i];
	}

	t_begin = clock();
	if (!UpdateLinearSystem( interpolationWeight , 1.0 , hierarchy , multigridFilteringCoefficients,
		deepMassCoefficients, deepStiffnessCoefficients,
		boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix,
		boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix,
		coarseSolver, boundarySolver, directSolver,
		filteringMatrix, DetailVerbose.set, true, UseDirectSolver.set)) {
		printf("ERROR : Failed system update! \n");
		return 0;
	}
	if( Verbose.set ) printf( "\tInitialized multigrid coefficients: %.2f(s)\n" , double(clock() - t_begin) / CLOCKS_PER_SEC );

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
	MultiplyBySystemMatrix_NoReciprocals(deepMassCoefficients, boundaryDeepMassMatrix, boundaryBoundaryMassMatrix, hierarchy.gridAtlases[0].boundaryGlobalIndex, hierarchy.gridAtlases[0].rasterLines, _x0, mass_x0);

	stiffness_x0.resize(textureNodes.size());
	MultiplyBySystemMatrix_NoReciprocals(deepStiffnessCoefficients, boundaryDeepStiffnessMatrix, boundaryBoundaryStiffnessMatrix, hierarchy.gridAtlases[0].boundaryGlobalIndex, hierarchy.gridAtlases[0].rasterLines, _x0, stiffness_x0);

#pragma omp parallel for
	for (int i = 0; i <_x0.size(); i++){
		multigridFilteringVariables[0].x[i] = _x0[i];
		multigridFilteringVariables[0].rhs[i] = mass_x0[i] * interpolationWeight + stiffness_x0[i] * gradientModulation;
	}

	filteredTexture.resize(width, height);
	for (int i = 0; i < filteredTexture.size(); i++) filteredTexture[i] = Point3D<float>(0.5f, 0.5f, 0.5f);


	if( UseDirectSolver.set ) ComputeExactSolution( Verbose.set );
	UpdateFilteredTexture( multigridFilteringVariables[0].x );
	return 1;
}

template<class Real>
void TextureFilter<Real>::InitializeVisualization( void )
{
	visualization.textureWidth = textureWidth;
	visualization.textureHeight = textureHeight;

	visualization.colorTextureBuffer = new unsigned char[textureHeight*textureWidth * 3];
	memset(visualization.colorTextureBuffer, 128, textureHeight * textureWidth * 3 * sizeof(unsigned char));



	int tCount = (int)mesh.triangles.size();

	visualization.triangles.resize(tCount);
	visualization.vertices.resize(3 * tCount);
	visualization.colors.resize(3 * tCount, Point3D<double>(0.75, 0.75, 0.75));
	visualization.textureCoordinates.resize(3 * tCount);
	visualization.normals.resize(3 * tCount);


	for (int i = 0; i < tCount; i++) for (int k = 0; k < 3; k++) visualization.triangles[i][k] = 3 * i + k;

	for (int i = 0; i<tCount; i++){
		for (int j = 0; j < 3; j++){
			visualization.vertices[3 * i + j] = mesh.vertices[mesh.triangles[i][j]];
			visualization.normals[3 * i + j] = mesh.normals[mesh.triangles[i][j]];
			visualization.textureCoordinates[3 * i + j] = mesh.textureCoordinates[3 * i + j];
		}
	}

	std::vector<int> boundaryEdges;
	if (!InitializeBoundaryEdges(mesh, boundaryEdges)) {
		printf("Unable to initialize boundary edges! \n");
	}

	for (int e = 0; e < boundaryEdges.size(); e++) {
		int tIndex = boundaryEdges[e] / 3;
		int kIndex = boundaryEdges[e] % 3;
		for (int c = 0; c < 2; c++) {
			Point3D<double> v = mesh.vertices[mesh.triangles[tIndex][(kIndex + c) % 3]];
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
		visualization.normalProgram->setup();
	}
	else visualization.textureType = COLOR_TEXTURE;
	delete[] ext;

	visualization.UpdateVertexBuffer();
	visualization.UpdateFaceBuffer();
	visualization.UpdateTextureBuffer(filteredTexture);

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

template< class Real >
int TextureFilter< Real >::Init( void )
{
	sprintf( stepsString , "Steps: 0" );
	levels = std::max<int>(Levels.value,1);
	interpolationWeight = InterpolationWeight.value;
	gradientModulation = GradientModulation.value;

	sprintf( gradientModulationStr , "Gradient modulation: %.2e\n" , gradientModulation );
	sprintf( interpolationStr , "Interpolation: %.2e\n" , interpolationWeight );

	if( !ReadTexturedMesh( mesh , Input.values[0] , Input.values[1] , DetailVerbose.set ) )
	{
		printf( "Unable to read mesh data: %s %s\n" , Input.values[0] , Input.values[1] );
		return 0;
	}

	textureWidth = mesh.texture.width();
	textureHeight = mesh.texture.height();

	//Define centroid and scale for visualization
	Point3D<double> centroid(0.f, 0.f, 0.f);
	for (int i = 0; i < mesh.vertices.size(); i++) centroid += mesh.vertices[i];
	centroid /= double(mesh.vertices.size());
	double radius = 0.f;
	for (int i = 0; i < mesh.vertices.size(); i++) radius = std::max<double>(radius, Point3D<double>::Length(mesh.vertices[i] - centroid));
	for (int i = 0; i < mesh.vertices.size(); i++) mesh.vertices[i] = (mesh.vertices[i] - centroid) / radius;

	if (1) for (int i = 0; i < mesh.textureCoordinates.size(); i++)mesh.textureCoordinates[i][1] = 1.0 - mesh.textureCoordinates[i][1];

	if( RandomJitter.set )
	{
		srand( time(NULL) );
		std::vector<Point2D < double >>randomOffset(mesh.vertices.size());
		double jitterScale = 1e-3 / double(std::max<int>(textureWidth, textureHeight));
		for (int i = 0; i < randomOffset.size(); i++) randomOffset[i] = Point2D < double >(1.0 - 2.0 * double(rand()) / double(RAND_MAX), 1.0 - 2.0 *  double(rand()) / double(RAND_MAX))*jitterScale;
		for (int i = 0; i < mesh.triangles.size(); i++) for (int k = 0; k < 3; k++)mesh.textureCoordinates[3 * i + k] += randomOffset[mesh.triangles[i][k]];
	}

	ComputePadding( padding , textureWidth , textureHeight , mesh.textureCoordinates , DetailVerbose.set );
	if( padding.nonTrivial )
	{
		PaddTextureCoordinates(padding, textureWidth, textureHeight, mesh.textureCoordinates);
		PaddImage( padding , mesh.texture );
		textureWidth += (padding.left + padding.right);
		textureHeight += (padding.bottom + padding.top);
	}

	clock_t t = clock();
	if( !InitializeSystem( textureWidth , textureHeight ) )
	{
		printf( "Unable to initialize system\n" );
		return 0;
	}
	if( Verbose.set )
	{
		printf( "Resolution: %d / %d x %d\n" , (int)textureNodes.size() , textureWidth , textureHeight );
		printf( "Initialized system %.2f(s)\n" , double(clock()-t) / CLOCKS_PER_SEC );
		printf( "Peak Memory (MB): %d\n" , Miscellany::MemoryInfo::PeakMemoryUsageMB() );
	}

	//Assign position to exterior nodes using barycentric-exponential map
	{
		FEM::RiemannianMesh< double > rMesh( GetPointer( mesh.triangles ) , mesh.triangles.size() );
		rMesh.setMetricFromEmbedding( GetPointer( mesh.vertices ) );
		rMesh.makeUnitArea();
		Pointer( FEM::CoordinateXForm< double > ) xForms = rMesh.getCoordinateXForms();

		for( int i=0 ; i<textureNodes.size() ; i++ ) if( textureNodes[i].tId!=-1 && !textureNodes[i].isInterior )
		{
			FEM::HermiteSamplePoint< double > _p;
			_p.tIdx = textureNodes[i].tId;
			_p.p = Point2D< double >( 1./3 , 1./3 );
			_p.v = textureNodes[i].barycentricCoords - _p.p;

			rMesh.exp(xForms, _p);
			
			textureNodes[i].tId = _p.tIdx;
			textureNodes[i].barycentricCoords = _p.p;
		}
	}

	textureNodePositions.resize(textureNodes.size());
	for (int i = 0; i < textureNodePositions.size(); i++){
		Point2D<double> barincetricCoords = textureNodes[i].barycentricCoords;
		int tId = textureNodes[i].tId;
		Point3D<float> surfacePosition = mesh.vertices[mesh.triangles[tId][0]] * (1.0 - barincetricCoords[0] - barincetricCoords[1]) +
										 mesh.vertices[mesh.triangles[tId][1]] * barincetricCoords[0] +
										 mesh.vertices[mesh.triangles[tId][2]] * barincetricCoords[1];
		textureNodePositions[i] = surfacePosition;
	}

	uniformTexelModulationMask.resize(textureNodes.size(), 0.5);

	for (int c = 0; c < 3; c++) texelStiffness[c].resize(textureNodes.size());

	cellModulationMask.resize( bilinearElementIndices.size() , 1 );
	uniformCellModulationMask.resize( bilinearElementIndices.size() , 0.5 );

	cellCenterPositions.resize( bilinearElementIndices.size() );
	for (int i = 0; i < bilinearElementIndices.size(); i++){
		cellCenterPositions[i] = (textureNodePositions[ bilinearElementIndices[i][0] ] +
								  textureNodePositions[ bilinearElementIndices[i][1] ] +
								  textureNodePositions[ bilinearElementIndices[i][2] ] +
								  textureNodePositions[ bilinearElementIndices[i][3] ]) / 4.0;
	}

	if( true )
	{
		int multiChartTexelCount = 0;
		Image< int > texelId;
		texelId.resize(textureWidth, textureHeight);
		for( int i=0 ; i<texelId.size() ; i++ ) texelId[i] = -1;
		for( int i=0 ; i<textureNodes.size() ; i++)
		{
			int ci = textureNodes[i].ci , cj = textureNodes[i].cj;
			if( texelId(ci,cj)!=-1 )
			{
				if( false ) fprintf( stderr , "[WARNING] Texel (%d %d) belong to multiple charts!\n" , ci , cj );
				multiChartTexelCount++;
			}
			texelId(ci,cj) = i;
		}
		if( multiChartTexelCount ) fprintf( stderr , "[WARNING] %d texels belong to multiple charts!\n" , multiChartTexelCount );
	}

	return 1;
}

template<class Real>
int _main(int argc, char* argv[])
{
	if( !TextureFilter< Real >::Init() ) return 0;

	if( !Output.set )
	{
		glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
		TextureFilter<Real>::visualization.displayMode = DisplayMode.value;
		if     ( DisplayMode.value==ONE_REGION_DISPLAY   ) TextureFilter< Real >::visualization.screenWidth =  800 , TextureFilter< Real >::visualization.screenHeight = 800;
		else if( DisplayMode.value==TWO_REGION_DISPLAY   ) TextureFilter< Real >::visualization.screenWidth = 1600 , TextureFilter< Real >::visualization.screenHeight = 800;
		else if( DisplayMode.value==THREE_REGION_DISPLAY ) TextureFilter< Real >::visualization.screenWidth = 1200 , TextureFilter< Real >::visualization.screenHeight = 800;
		else if( DisplayMode.value==FOUR_REGION_DISPLAY  ) TextureFilter< Real >::visualization.screenWidth = 1500 , TextureFilter< Real >::visualization.screenHeight = 600;
		glutInitWindowSize(TextureFilter<Real>::visualization.screenWidth, TextureFilter<Real>::visualization.screenHeight);
		glutInit(&argc, argv);
		char windowName[1024];
		sprintf(windowName, "Texture Filtering");
		glutCreateWindow(windowName);
		if (glewInit() != GLEW_OK) fprintf(stderr, "[ERROR] glewInit failed\n"), exit(0);
		glutDisplayFunc(TextureFilter<Real>::Display);
		glutReshapeFunc(TextureFilter<Real>::Reshape);
		glutMouseFunc(TextureFilter<Real>::MouseFunc);
		glutMotionFunc(TextureFilter<Real>::MotionFunc);
		glutKeyboardFunc(TextureFilter<Real>::KeyboardFunc);
		glutIdleFunc( TextureFilter< Real >::Idle );
		if (CameraConfig.set) TextureFilter<Real>::visualization.ReadSceneConfigurationCallBack(&TextureFilter<Real>::visualization, CameraConfig.value);
		TextureFilter<Real>::InitializeVisualization();
		TextureFilter<Real>::visualization.showSlideBar = true;
		glutMainLoop(); 
	}
	else
	{
		if( UseDirectSolver.set ) TextureFilter< Real >::ComputeExactSolution( Verbose.set );
		else for ( int i=0 ; i<OutputVCycles.value ; i++ ) TextureFilter< Real >::UpdateSolution();
		TextureFilter< Real >::ExportTextureCallBack( &TextureFilter< Real >::visualization , Output.value );
	}

	return 0;
}

int main(int argc, char* argv[])
{
	cmdLineParse( argc-1 , argv+1 , params );
	if( !Input.set ) { ShowUsage( argv[0] ) ; return EXIT_FAILURE; }
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
	if( Double.set ) _main< double >( argc , argv );
	else             _main< float  >( argc , argv );
	return 0;
}
