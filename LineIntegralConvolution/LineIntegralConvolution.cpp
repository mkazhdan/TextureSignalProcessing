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
#include <Src/Basis.h>
#include <Misha/FEM.h>
#include <Src/ChartDecomposition.h>
#include <Src/HSV.h>
#include <Src/Solver.h>
#include <Src/Hierarchy.h>
#include <Src/MassAndStiffness.h>
#include <Src/Padding.h>
#include <Src/TexturedMeshVisualization.h>
#include <Misha/Miscellany.h>

cmdLineParameter< char* > Input( "in" );
cmdLineParameter< char* > Output( "out" );
cmdLineParameter< int   > OutputVCycles( "outVCycles" , 10 );
cmdLineReadable MinimalCurvature( "minimal" );
cmdLineReadable Double( "double" );
cmdLineParameter< char* > VectorField( "vf" );
cmdLineParameter< int   > Width( "width" , 2048 );
cmdLineParameter< int   > Height( "height" , 2048 );
cmdLineParameter< float > LICInterpolationWeight( "licInterpolation" , 1e4 );
cmdLineParameter< float > SharpeningInterpolationWeight( "sharpInterpolation" , 1e4 );
cmdLineParameter< float > SharpeningGradientModulation( "sharpModulation" , 100 );
cmdLineParameter< float > AnisoExponent( "aExp" , 0.f );
cmdLineParameter< int   > Levels( "levels" , 4 );
cmdLineParameter< int   > MatrixQuadrature( "mQuadrature" , 6 );


cmdLineParameter< char* > CameraConfig("camera");
cmdLineParameter< int   > Threads("threads", omp_get_num_procs());
cmdLineParameter< int   > DisplayMode("display", TWO_REGION_DISPLAY);


cmdLineParameter< int   > MultigridBlockHeight("mBlockH", 16);
cmdLineParameter< int   > MultigridBlockWidth("mBlockW", 128);
cmdLineParameter< int   > MultigridPaddedHeight("mPadH", 0);
cmdLineParameter< int   > MultigridPaddedWidth("mPadW", 2);

cmdLineReadable RandomJitter( "jitter" );
cmdLineReadable Verbose("verbose");
cmdLineReadable DetailVerbose("detail");
cmdLineReadable UseDirectSolver("useDirectSolver");
cmdLineReadable IntrinsicVectorField( "intrinsicVF" );
cmdLineReadable NoHelp( "noHelp" );

cmdLineReadable* params[] =
{
	&Input , &Output , &MinimalCurvature , &VectorField , &IntrinsicVectorField , &Width,&Height , &LICInterpolationWeight , &SharpeningInterpolationWeight , &SharpeningGradientModulation , &CameraConfig, &Levels,&UseDirectSolver,&Threads,&DisplayMode,&MultigridBlockHeight,&MultigridBlockWidth,&MultigridPaddedHeight,&MultigridPaddedWidth,&Verbose,
	&DetailVerbose , &RandomJitter ,
	&Double ,
	&MatrixQuadrature ,
	&OutputVCycles ,
	&NoHelp , &AnisoExponent ,
	NULL
};

void ShowUsage(const char* ex)
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input mesh>\n" , Input.name );
	printf( "\t[--%s <output texture>\n" , Output.name );
	printf( "\t[--%s <output v-cycles>=%d]\n" , OutputVCycles.name , OutputVCycles.value );
	printf( "\t[--%s <vector field file>\n" , VectorField.name );
	printf( "\t[--%s <LIC interpolation weight>=%f]\n" , LICInterpolationWeight.name , LICInterpolationWeight.value );
	printf( "\t[--%s <sharpening interpolation weight>=%f]\n" , SharpeningInterpolationWeight.name , SharpeningInterpolationWeight.value   );
	printf( "\t[--%s <sharpening gradient modulation>=%f]\n" , SharpeningGradientModulation.name , SharpeningGradientModulation.value );
	printf( "\t[--%s <texture width>=%d]\n" , Width.name  , Width.value  );
	printf( "\t[--%s <texture height>=%d]\n", Height.name , Height.value );
	printf( "\t[--%s <system matrix quadrature points per triangle>=%d]\n" , MatrixQuadrature.name , MatrixQuadrature.value );
	printf( "\t[--%s]\n" , IntrinsicVectorField.name );
	printf( "\t[--%s]\n" , MinimalCurvature.name );
	printf( "\t[--%s]\n" , UseDirectSolver.name );
	printf( "\t[--%s]\n" , RandomJitter.name );
	printf( "\t[--%s]\n" , Verbose.name );

	printf( "\t[--%s <camera configuration file>\n" , CameraConfig.name);
	printf( "\t[--%s <hierarchy levels>=%d]\n" , Levels.name , Levels.value );
	printf( "\t[--%s <threads>=%d]\n" , Threads.name , Threads.value );
	printf( "\t[--%s]\n" , DetailVerbose.name );
	printf( "\t[--%s <display mode>=%d]\n" , DisplayMode.name , DisplayMode.value );
	printf( "\t\t%d] One Region \n", ONE_REGION_DISPLAY);
	printf( "\t\t%d] Two Region \n", TWO_REGION_DISPLAY);

	printf( "\t[--%s <multigrid block width>=%d]\n"   , MultigridBlockWidth.name   , MultigridBlockWidth.value   );
	printf( "\t[--%s <multigrid block height>=%d]\n"  , MultigridBlockHeight.name  , MultigridBlockHeight.value  );
	printf( "\t[--%s <multigrid padded width>=%d]\n"  , MultigridPaddedWidth.name  , MultigridPaddedWidth.value  );
	printf( "\t[--%s <multigrid padded height>=%d]\n" , MultigridPaddedHeight.name , MultigridPaddedHeight.value );
	printf( "\t[--%s <anisotropy exponent>=%f]\n" , AnisoExponent.name , AnisoExponent.value );
	printf( "\t[--%s]\n" , NoHelp.name );
}

template<class Real>
class LineConvolution
{
public:
	static Real sharpeningGradientModulation;
	static Real sharpeningInterpolationWeight;
	static Real licInterpolationWeight;

	static TexturedMesh mesh;
	static int textureWidth;
	static int textureHeight;
	static int levels;
	static std::vector<Point3D<float>> textureNodePositions;

	static Padding padding;

	static HierarchicalSystem hierarchy;
	static std::vector< BilinearElementIndex > bilinearElementIndices;

	static std::vector<TextureNodeInfo> textureNodes;
	static Image<int> nodeIndex;

	static SparseMatrix<double, int> anisotropicMass;
	static SparseMatrix<double, int> anisotropicStiffness;

	static SparseMatrix<double, int> mass;
	static SparseMatrix<double, int> stiffness;


	static SparseMatrix< double, int > lineConvolutionMatrix;
	static SparseMatrix< double, int > modulationMatrix;

	static int impulseTexel;

	static std::vector<AtlasChart> atlasCharts;
	static std::vector<std::vector<SquareMatrix<double, 2>>> parameterMetric;

	static double lineConvolutionRange;
	static double modulationRange;

	static std::vector< Point3D< Real > > randSignal;

	static std::vector< Point3D< Real > > mass_x0;
	static std::vector< Point3D< Real > > stiffness_x0;

	//Impulse Smoothing
	static std::vector< MultigridLevelCoefficients< Real > > multigridLineConvolutionCoefficients;
	static std::vector< MultigridLevelVariables< Point3D< Real > > > multigridLineConvolutionVariables;

	//Geodesic Distance
	static std::vector< MultigridLevelCoefficients< Real > > multigridModulationCoefficients;
	static std::vector< MultigridLevelVariables< Point3D< Real > > > multigridModulationVariables;

#if defined( USE_CHOLMOD )
	typedef  std::vector< CholmodCholeskySolver< Real , 3 > > BoundarySolverType;
	typedef  CholmodCholeskySolver< Real , 3 >  CoarseSolverType;
	typedef  CholmodCholeskySolver< Real , 3 > DirectSolverType;
#elif defined( USE_EIGEN_SIMPLICIAL )
	typedef  std::vector< EigenCholeskySolver< Real , 3 > > BoundarySolverType;
	typedef  EigenCholeskySolver< Real , 3 >  CoarseSolverType;
	typedef  EigenCholeskySolver< Real , 3 > DirectSolverType;
#elif defined( USE_EIGEN_PARDISO )
	typedef  std::vector< EigenPardisoSolver< Real , 3 > > BoundarySolverType;
	typedef  EigenPardisoSolver< Real , 3 > CoarseSolverType;
	typedef  EigenPardisoSolver< Real , 3 > DirectSolverType;
#else
#error "[ERROR] No solver defined!"
#endif

	static BoundarySolverType	boundaryLineConvolutionSolver;
	static BoundarySolverType	boundaryModulationSolver;

	static CoarseSolverType coarseLineConvolutionSolver;
	static CoarseSolverType coarseModulationSolver;

	static DirectSolverType fineLineConvolutionSolver;
	static DirectSolverType fineModulationSolver;

	static std::vector<MultigridLevelIndices<Real>> multigridIndices;

	static SparseMatrix< Real , int > coarseBoundaryFineBoundaryProlongation;
	static SparseMatrix< Real , int > fineBoundaryCoarseBoundaryRestriction;
	static std::vector< Real > coarseBoundaryValues;
	static std::vector< Real > coarseBoundaryRHS;
	static std::vector< Real > fineBoundaryValues;
	static std::vector< Real > fineBoundaryRHS;

	//Anisotropic Linear Operators
	static std::vector< double > anisoDeepMassCoefficients;
	static std::vector< double > anisoDeepStiffnessCoefficients;
	static SparseMatrix< double , int > anisoBoundaryBoundaryMassMatrix;
	static SparseMatrix< double , int > anisoBoundaryBoundaryStiffnessMatrix;
	static SparseMatrix< double , int > anisoBoundaryDeepMassMatrix;
	static SparseMatrix< double , int > anisoBoundaryDeepStiffnessMatrix;

	//Isotropic Linear Operators
	static std::vector< double > deepMassCoefficients;
	static std::vector< double > deepStiffnessCoefficients;
	static SparseMatrix< double , int> boundaryBoundaryMassMatrix;
	static SparseMatrix< double , int> boundaryBoundaryStiffnessMatrix;
	static SparseMatrix< double , int> boundaryDeepMassMatrix;
	static SparseMatrix< double , int> boundaryDeepStiffnessMatrix;

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

	static int Init();
	static void InitializeVisualization( const int width , const int height );
	static void ComputeExactSolution( bool verbose= false );
	static int UpdateSolution( bool verbose=false , bool detailVerbose=false );
	static int InitializeSystem( const int width , const int height , double scale );

	static void Display(void) { visualization.Display(); }
	static void MouseFunc(int button, int state, int x, int y);
	static void MotionFunc(int x, int y);
	static void Reshape(int w, int h) { visualization.Reshape(w, h); }
	static void KeyboardFunc(unsigned char key, int x, int y) { visualization.KeyboardFunc(key, x, y); }
	static void Idle();
};

template<class Real> Real														LineConvolution<Real>::sharpeningGradientModulation;
template<class Real> Real														LineConvolution<Real>::sharpeningInterpolationWeight;
template<class Real> Real														LineConvolution<Real>::licInterpolationWeight;

template<class Real> TexturedMesh												LineConvolution<Real>::mesh;
template<class Real> int														LineConvolution<Real>::textureWidth;
template<class Real> int														LineConvolution<Real>::textureHeight;

template<class Real> TexturedMeshVisualization									LineConvolution<Real>::visualization;

template<class Real> std::vector<AtlasChart>									LineConvolution<Real>::atlasCharts;
template<class Real> std::vector<std::vector<SquareMatrix<double, 2>>>			LineConvolution<Real>::parameterMetric;

template<class Real> Padding													LineConvolution<Real>::padding;

template<class Real> SparseMatrix<double, int>									LineConvolution<Real>::anisotropicMass;
template<class Real> SparseMatrix<double, int>									LineConvolution<Real>::anisotropicStiffness;


template<class Real> SparseMatrix<double, int>									LineConvolution<Real>::lineConvolutionMatrix;
template<class Real> SparseMatrix<double, int>									LineConvolution<Real>::modulationMatrix;

template<class Real> SparseMatrix<double, int>									LineConvolution<Real>::mass;
template<class Real> SparseMatrix<double, int>									LineConvolution<Real>::stiffness;


template<class Real> std::vector<TextureNodeInfo>								LineConvolution<Real>::textureNodes;
template<class Real> Image<int>													LineConvolution<Real>::nodeIndex;
template<class Real> std::vector< BilinearElementIndex >						LineConvolution<Real>::bilinearElementIndices;

template<class Real> int														LineConvolution<Real>::levels;
template<class Real> HierarchicalSystem											LineConvolution<Real>::hierarchy;

template<class Real> unsigned char *											LineConvolution<Real>::outputBuffer;
template<class Real> std::vector<MultigridLevelIndices<Real>>					LineConvolution<Real>::multigridIndices;

//Impulse Smoothing
template<class Real> std::vector<MultigridLevelCoefficients<Real>>				LineConvolution<Real>::multigridLineConvolutionCoefficients;
template<class Real> std::vector< MultigridLevelVariables< Point3D< Real > > >	LineConvolution<Real>::multigridLineConvolutionVariables;
template<class Real> typename LineConvolution<Real>::CoarseSolverType			LineConvolution<Real>::coarseLineConvolutionSolver;

//Geodesic Distance
template<class Real> std::vector<MultigridLevelCoefficients<Real>>				LineConvolution<Real>::multigridModulationCoefficients;
template<class Real> std::vector< MultigridLevelVariables< Point3D< Real > > >	LineConvolution<Real>::multigridModulationVariables;
template<class Real> typename LineConvolution<Real>::CoarseSolverType						LineConvolution<Real>::coarseModulationSolver;

template<class Real> typename LineConvolution<Real>::DirectSolverType						LineConvolution<Real>::fineLineConvolutionSolver;
template<class Real> typename LineConvolution<Real>::DirectSolverType						LineConvolution<Real>::fineModulationSolver;

template<class Real>  typename LineConvolution<Real>::BoundarySolverType					LineConvolution<Real>::boundaryLineConvolutionSolver;
template<class Real>  typename LineConvolution<Real>::BoundarySolverType					LineConvolution<Real>::boundaryModulationSolver;

template<class Real> std::vector< Point3D< Real > >		LineConvolution<Real>::randSignal;

template<class Real> std::vector< Point3D< Real > >		LineConvolution<Real>::mass_x0;
template<class Real> std::vector< Point3D< Real > >		LineConvolution<Real>::stiffness_x0;

template<class Real> static std::vector< Point3D< Real > > mass_x0;
template<class Real> static std::vector< Point3D< Real > > stiffness_x0;

template<class Real> int															LineConvolution<Real>::impulseTexel = -1;
template<class Real> std::vector<Point3D<float>>									LineConvolution<Real>::textureNodePositions;

template<class Real> double															LineConvolution<Real>::lineConvolutionRange;
template<class Real> double															LineConvolution<Real>::modulationRange;

template<class Real> SparseMatrix<Real, int>										LineConvolution<Real>::coarseBoundaryFineBoundaryProlongation;
template<class Real> SparseMatrix<Real, int>										LineConvolution<Real>::fineBoundaryCoarseBoundaryRestriction;

template<class Real> std::vector<Real>												LineConvolution<Real>::coarseBoundaryValues;
template<class Real> std::vector<Real>												LineConvolution<Real>::coarseBoundaryRHS;
template<class Real> std::vector<Real>												LineConvolution<Real>::fineBoundaryValues;
template<class Real> std::vector<Real>												LineConvolution<Real>::fineBoundaryRHS;

template<class Real> std::vector<double>											LineConvolution<Real>::anisoDeepMassCoefficients;
template<class Real> std::vector<double>											LineConvolution<Real>::anisoDeepStiffnessCoefficients;
template<class Real> SparseMatrix<double, int>										LineConvolution<Real>::anisoBoundaryBoundaryMassMatrix;
template<class Real> SparseMatrix<double, int>										LineConvolution<Real>::anisoBoundaryBoundaryStiffnessMatrix;
template<class Real> SparseMatrix<double, int>										LineConvolution<Real>::anisoBoundaryDeepMassMatrix;
template<class Real> SparseMatrix<double, int>										LineConvolution<Real>::anisoBoundaryDeepStiffnessMatrix;

template<class Real> std::vector<double>											LineConvolution<Real>::deepMassCoefficients;
template<class Real> std::vector<double>											LineConvolution<Real>::deepStiffnessCoefficients;
template<class Real> SparseMatrix<double, int>										LineConvolution<Real>::boundaryBoundaryMassMatrix;
template<class Real> SparseMatrix<double, int>										LineConvolution<Real>::boundaryBoundaryStiffnessMatrix;
template<class Real> SparseMatrix<double, int>										LineConvolution<Real>::boundaryDeepMassMatrix;
template<class Real> SparseMatrix<double, int>										LineConvolution<Real>::boundaryDeepStiffnessMatrix;

template< class Real > int															LineConvolution<Real>::updateCount = 0;
template< class Real >
void LineConvolution<Real>::ComputeExactSolution( bool verbose )
{
	clock_t begin;

	// (1) Line Convolution	
	// RHS = Mass * randSignal * licInterpolationWeight
	MultiplyBySystemMatrix_NoReciprocals( anisoDeepMassCoefficients , anisoBoundaryDeepMassMatrix , anisoBoundaryBoundaryMassMatrix , hierarchy.gridAtlases[0].boundaryGlobalIndex , hierarchy.gridAtlases[0].rasterLines , randSignal , multigridLineConvolutionVariables[0].rhs );
#pragma omp parallel for
	for( int i=0 ; i<textureNodes.size() ; i++ ) multigridLineConvolutionVariables[0].rhs[i] *= licInterpolationWeight;

	begin = clock();
	solve( fineLineConvolutionSolver , multigridLineConvolutionVariables[0].x , multigridLineConvolutionVariables[0].rhs );
	if( verbose ) printf( "Line convolution %.4f\n" , double( clock()-begin ) / CLOCKS_PER_SEC );

	//(2) Compute modulation RHS
	mass_x0.resize( textureNodes.size() );
	MultiplyBySystemMatrix_NoReciprocals( deepMassCoefficients , boundaryDeepMassMatrix , boundaryBoundaryMassMatrix , hierarchy.gridAtlases[0].boundaryGlobalIndex , hierarchy.gridAtlases[0].rasterLines , multigridLineConvolutionVariables[0].x , mass_x0 );

	stiffness_x0.resize(textureNodes.size());
	MultiplyBySystemMatrix_NoReciprocals( deepStiffnessCoefficients , boundaryDeepStiffnessMatrix , boundaryBoundaryStiffnessMatrix , hierarchy.gridAtlases[0].boundaryGlobalIndex , hierarchy.gridAtlases[0].rasterLines , multigridLineConvolutionVariables[0].x , stiffness_x0 );

#pragma omp parallel for
	for( int i=0 ; i<textureNodes.size() ; i++ ) multigridModulationVariables[0].rhs[i] = mass_x0[i] * sharpeningInterpolationWeight + stiffness_x0[i] * sharpeningGradientModulation;

	//(3) Modulation
	if( verbose ) begin = clock();
	solve( fineModulationSolver , multigridModulationVariables[0].x , multigridModulationVariables[0].rhs );
	if (verbose ) printf( "Modulation %.4f \n" , double( clock()-begin ) / CLOCKS_PER_SEC);
}

template<class Real>
void LineConvolution<Real>::SetOutputBuffer( const std::vector< Point3D< Real > >& solution )
{
#pragma omp parallel for
	for( int i=0 ; i<textureNodes.size() ; i++ )
	{
		int ci = textureNodes[i].ci;
		int cj = textureNodes[i].cj;
		int offset = 3 * (textureWidth*cj + ci);
		outputBuffer[offset+0] = (unsigned char)( std::min< double >( std::max< double >( 0 , solution[i][0] ) , 1.0 )*255.0 );
		outputBuffer[offset+1] = (unsigned char)( std::min< double >( std::max< double >( 0 , solution[i][1] ) , 1.0 )*255.0 );
		outputBuffer[offset+2] = (unsigned char)( std::min< double >( std::max< double >( 0 , solution[i][2] ) , 1.0 )*255.0 );
	}
}

template<class Real>
void LineConvolution<Real>::UpdateOutputBuffer( const std::vector< Point3D< Real > >& solution )
{
	SetOutputBuffer( solution );

	glBindTexture(GL_TEXTURE_2D, visualization.textureBuffer);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, textureWidth, textureHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)&outputBuffer[0]);
	glBindTexture(GL_TEXTURE_2D, 0);
	glutPostRedisplay();
}

template< class Real >
void LineConvolution<Real>::Idle( void )
{
	if( updateCount )
	{
		if( !UpdateSolution() ) fprintf( stderr , "[ERROR] Updated solution failed!\n" );
		else if( updateCount>0 ) updateCount--;
	}
	UpdateOutputBuffer( multigridModulationVariables[0].x );
}

template<class Real>
void LineConvolution<Real>::MouseFunc(int button, int state, int x, int y) {

	visualization.newX = x; visualization.newY = y;
	visualization.rotating = visualization.scaling = visualization.panning = false;

	if( ( button==GLUT_LEFT_BUTTON || button==GLUT_RIGHT_BUTTON ) && glutGetModifiers() & GLUT_ACTIVE_CTRL) visualization.panning = true;
	else if( button==GLUT_LEFT_BUTTON  ) visualization.rotating = true;
	else if( button==GLUT_RIGHT_BUTTON ) visualization.scaling  = true;
}

template<class Real>
void LineConvolution<Real>::MotionFunc(int x, int y) {

	if (!visualization.showMesh) {
		visualization.oldX = visualization.newX, visualization.oldY = visualization.newY, visualization.newX = x, visualization.newY = y;

		int imageSize = std::min< int >(visualization.screenWidth, visualization.screenHeight);
		if (visualization.panning) visualization.xForm.offset[0] -= (visualization.newX - visualization.oldX) / visualization.imageToScreenScale(), visualization.xForm.offset[1] += (visualization.newY - visualization.oldY) / visualization.imageToScreenScale();
		else
		{
			float dz = float(pow(1.1, double(visualization.newY - visualization.oldY) / 8));
			visualization.xForm.zoom *= dz;
		}

	}
	else {
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
	glutPostRedisplay();
}


template<class Real>
void LineConvolution<Real>::SharpeningInterpolationWeightCallBack( Visualization* v , const char* prompt )
{
	for( int i=0 ; i<multigridLineConvolutionVariables[0].x.size() ; i++) multigridLineConvolutionVariables[0].x[i] *= 0;
	for( int i=0 ; i<multigridModulationVariables[0].x.size() ; i++) multigridModulationVariables[0].x[i] *= 0;

	sharpeningInterpolationWeight = atof(prompt);

	if( UseDirectSolver.set ) modulationMatrix = mass * sharpeningInterpolationWeight + stiffness;

	if (!UpdateLinearSystem( sharpeningInterpolationWeight , 1.0 , hierarchy , multigridModulationCoefficients ,
		deepMassCoefficients, deepStiffnessCoefficients,
		boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix,
		boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix,
		coarseModulationSolver, boundaryModulationSolver, fineModulationSolver,
		modulationMatrix, DetailVerbose.set, false, UseDirectSolver.set)) {
		printf("ERROR : Failed system update! \n");
	}

	if( UseDirectSolver.set )
	{
		ComputeExactSolution();
		UpdateOutputBuffer( multigridModulationVariables[0].x );
	}
}

template< class Real >
void LineConvolution< Real >::LICInterpolationWeightCallBack( Visualization* v , const char* prompt )
{
	for( int i=0 ; i<multigridLineConvolutionVariables[0].x.size() ; i++) multigridLineConvolutionVariables[0].x[i] *= 0;
	for( int i=0 ; i<multigridModulationVariables[0].x.size() ; i++) multigridModulationVariables[0].x[i] *= 0;

	licInterpolationWeight = atof(prompt);

	if( UseDirectSolver.set ) lineConvolutionMatrix = anisotropicMass * licInterpolationWeight + anisotropicStiffness;

	if (!UpdateLinearSystem( licInterpolationWeight , 1.0 , hierarchy , multigridLineConvolutionCoefficients ,
		anisoDeepMassCoefficients, anisoDeepStiffnessCoefficients,
		anisoBoundaryBoundaryMassMatrix, anisoBoundaryBoundaryStiffnessMatrix,
		anisoBoundaryDeepMassMatrix, anisoBoundaryDeepStiffnessMatrix,
		coarseLineConvolutionSolver, boundaryLineConvolutionSolver, fineLineConvolutionSolver,
		lineConvolutionMatrix, DetailVerbose.set, false, UseDirectSolver.set)) {
		printf("ERROR : Failed system update! \n");
	}

	if( UseDirectSolver.set )
	{
		ComputeExactSolution();
		UpdateOutputBuffer( multigridModulationVariables[0].x );
	}


}

template< class Real >
void LineConvolution< Real >::SharpeningGradientModulationCallBack( Visualization* v , const char* prompt )
{
	for( int i=0 ; i<multigridLineConvolutionVariables[0].x.size() ; i++ ) multigridLineConvolutionVariables[0].x[i] *= 0;
	for( int i=0 ; i<multigridModulationVariables[0].x.size() ; i++ ) multigridModulationVariables[0].x[i] *= 0;

	sharpeningGradientModulation = atof(prompt);
	if( UseDirectSolver.set )
	{
		ComputeExactSolution();
		UpdateOutputBuffer( multigridModulationVariables[0].x );
	}
}

template< class Real > void LineConvolution< Real >::ToggleUpdateCallBack( Visualization* v , const char* prompt )
{
	if( updateCount ) updateCount = 0;
	else              updateCount = -1;
}
template< class Real > void LineConvolution< Real >::IncrementUpdateCallBack( Visualization* v , const char* prompt )
{
	if( updateCount<0 ) updateCount = 1;
	else updateCount++;
}

template<class Real>
void LineConvolution<Real>::ExportTextureCallBack( Visualization* v , const char* prompt )
{
	Image< Point3D< float > > outputImage;
	outputImage.resize( textureWidth , textureHeight );
	for( int i=0 ; i<outputImage.size() ; i++ ) outputImage[i] = Point3D< float >( outputBuffer[3*i] , outputBuffer[3*i+1] , outputBuffer[3*i+2] ) / 255.f;
	if( padding.nonTrivial ) UnpadImage( padding , outputImage );
	outputImage.write( prompt );
}

template<class Real>
int LineConvolution< Real >::UpdateSolution( bool verbose , bool detailVerbose )
{
	clock_t begin;

	//(1)Update smoothed input solution
	MultiplyBySystemMatrix_NoReciprocals( anisoDeepMassCoefficients , anisoBoundaryDeepMassMatrix , anisoBoundaryBoundaryMassMatrix , hierarchy.gridAtlases[0].boundaryGlobalIndex , hierarchy.gridAtlases[0].rasterLines , randSignal , multigridLineConvolutionVariables[0].rhs );
#pragma omp parallel for
	for( int i=0 ; i<textureNodes.size() ; i++ ) multigridLineConvolutionVariables[0].rhs[i] *= licInterpolationWeight;
	if (verbose) begin = clock();
	VCycle(multigridLineConvolutionVariables, multigridLineConvolutionCoefficients, multigridIndices, boundaryLineConvolutionSolver, coarseLineConvolutionSolver, detailVerbose, detailVerbose);
	if (verbose) printf("Smoothing impulse %.4f \n", double(clock() - begin) / CLOCKS_PER_SEC);

	//(2) Compute modulation RHS
	mass_x0.resize(textureNodes.size());
	MultiplyBySystemMatrix_NoReciprocals(deepMassCoefficients, boundaryDeepMassMatrix, boundaryBoundaryMassMatrix, hierarchy.gridAtlases[0].boundaryGlobalIndex, hierarchy.gridAtlases[0].rasterLines, multigridLineConvolutionVariables[0].x, mass_x0);

	stiffness_x0.resize(textureNodes.size());
	MultiplyBySystemMatrix_NoReciprocals(deepStiffnessCoefficients, boundaryDeepStiffnessMatrix, boundaryBoundaryStiffnessMatrix, hierarchy.gridAtlases[0].boundaryGlobalIndex, hierarchy.gridAtlases[0].rasterLines, multigridLineConvolutionVariables[0].x, stiffness_x0);

#pragma omp parallel for
	for( int i=0 ; i<textureNodes.size() ; i++ ) multigridModulationVariables[0].rhs[i] = mass_x0[i] * sharpeningInterpolationWeight + stiffness_x0[i] * sharpeningGradientModulation;

	//(3) Update geodesic distance solution	
	if (verbose) begin = clock();
	VCycle(multigridModulationVariables, multigridModulationCoefficients, multigridIndices, boundaryModulationSolver, coarseModulationSolver, detailVerbose, detailVerbose);
	if (verbose) printf("Solving geodesic distance %.4f \n", double(clock() - begin) / CLOCKS_PER_SEC);

	return 1;
}


template<class Real>
int LineConvolution<Real>::InitializeSystem( const int width , const int height , double scale )
{
	clock_t t_begin;
		
	t_begin = clock();
	MultigridBlockInfo multigridBlockInfo(MultigridBlockWidth.value, MultigridBlockHeight.value, MultigridPaddedWidth.value, MultigridPaddedHeight.value, 0);
	if( !InitializeHierarchy( mesh , width , height , levels , textureNodes , bilinearElementIndices , hierarchy , atlasCharts , multigridBlockInfo , true , DetailVerbose.set ) )
	{
		fprintf( stderr , "[ERROR] Failed intialization!\n" );
		return 0;
	}
	if( Verbose.set ) printf( "\tInitialized hierarchy: %.2f(s)\n" , double(clock() - t_begin) / CLOCKS_PER_SEC);

	//Initialize node index
	nodeIndex.resize(width, height);
	for (int i = 0; i < nodeIndex.size(); i++)nodeIndex[i] = -1;
	for (int i = 0; i < textureNodes.size(); i++) {
		if (textureNodes[i].ci < 0 || textureNodes[i].ci > textureWidth - 1 || textureNodes[i].cj < 0 || textureNodes[i].cj > textureHeight - 1) {
			printf("Error: Invalid node! %d %d \n", textureNodes[i].ci, textureNodes[i].cj);
			return 0;
		}
		if (nodeIndex(textureNodes[i].ci, textureNodes[i].cj) != -1) {
			if (0)printf("[WARNING] Multiple nodes mapped to pixel %d %d!\n", textureNodes[i].ci, textureNodes[i].cj);
		}
		nodeIndex(textureNodes[i].ci, textureNodes[i].cj) = i;
	}

	BoundaryProlongationData boundaryProlongation;
	if (!InitializeBoundaryProlongationData(hierarchy.gridAtlases[0], boundaryProlongation)) {
		printf("ERROR : Failed boundary prolongation! \n");
		return 0;
	}

	//////////////////////////////////// Initialize multigrid indices
	multigridIndices.resize(levels);
	for (int i = 0; i < levels; i++){
		const GridAtlas & gridAtlas = hierarchy.gridAtlases[i];
		multigridIndices[i].threadTasks = gridAtlas.threadTasks;
		multigridIndices[i].boundaryGlobalIndex = gridAtlas.boundaryGlobalIndex;
		multigridIndices[i].segmentedLines = gridAtlas.segmentedLines;
		multigridIndices[i].rasterLines = gridAtlas.rasterLines;
		multigridIndices[i].restrictionLines = gridAtlas.restrictionLines;
		multigridIndices[i].prolongationLines = gridAtlas.prolongationLines;
		if (i < levels - 1) {
			multigridIndices[i].boundaryRestriction = hierarchy.boundaryRestriction[i];
		}
	}

	//////////////////////////////////// Initialize multigrid coefficients

	//////////////////////////////////// 	Line Convolution coefficients
	{
		std::vector< Point2D< double > > vectorField;
		if( VectorField.set )
		{
			if( IntrinsicVectorField.set )
			{
				if( !ReadVector( vectorField , VectorField.value ) ){ fprintf( stderr , "[ERROR] Unable to read vector field: %s\n" , VectorField.value ) ; return 0; }
				if( vectorField.size()!=mesh.triangles.size() ){ fprintf( stderr , "[ERROR] Triangle and vector counts don't match: %d != %d\n" , (int)mesh.triangles.size() , (int)vectorField.size() ) ; return 0; }
			}
			else
			{
				std::vector< Point3D< double > > _vectorField;
				if( !ReadVector( _vectorField , VectorField.value ) ){ fprintf( stderr , "[ERROR] Unable to read vector field: %s\n" , VectorField.value ) ; return 0; }
				if( _vectorField.size()!=mesh.triangles.size() ){ fprintf( stderr , "[ERROR] Triangle and vector counts don't match: %d != %d\n" , (int)mesh.triangles.size() , (int)_vectorField.size() ) ; return 0; }
				vectorField.resize( _vectorField.size() );
#pragma omp parallel for
				for( int i=0 ; i<mesh.triangles.size() ; i++ )
				{
					Point3D< double > v[] = { mesh.vertices[ mesh.triangles[i][0] ] , mesh.vertices[ mesh.triangles[i][1] ] , mesh.vertices[ mesh.triangles[i][2] ] };
					Point3D< double > d[] = { (v[1]-v[0])*scale , (v[2]-v[0])*scale };
					SquareMatrix< double , 2 > Dot;
					for( int j=0 ; j<2 ; j++ ) for( int k=0 ; k<2 ; j++ ) Dot( j , k ) = Point3D< double >::Dot( d[j] , d[k] );
					Point2D< double > dot( Point3D< double >::Dot( d[0] , _vectorField[i] ) , Point3D< double >::Dot( d[1] , _vectorField[i] ) );
					vectorField[i] = Dot.inverse() * dot;
				}
			}
		}
		else
		{
			// Compute the principal curvatures
			std::vector< PrincipalCurvature< double > > principalCurvatures;
			UpdateNormals( mesh );
			Eigen::SparseMatrix< double > meshMassMatrix;
			Eigen::SparseMatrix< double > meshStiffnessMatrix;
			InitializeMeshMatrices( mesh , meshMassMatrix , meshStiffnessMatrix );
			Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > > meshSolver( meshMassMatrix + meshStiffnessMatrix*1e-4 );
			SmoothSignal( meshMassMatrix , meshSolver , mesh.normals , true );

			InitializePrincipalCurvatureDirection( mesh , mesh.normals , principalCurvatures );

			// Set the vector-field to the principal curvature direction times the umbilicity
			vectorField.resize( principalCurvatures.size() );
#pragma omp parallel for
			for( int t=0 ; t<principalCurvatures.size() ; t++ ) vectorField[t] = principalCurvatures[t].dirs[ MinimalCurvature.set ? 0 : 1 ] * ( principalCurvatures[t].values[1] - principalCurvatures[t].values[0] );
		}
		{
			std::vector< SquareMatrix< double , 2 > > embeddingMetric;

			InitializeEmbeddingMetric( mesh , true , embeddingMetric );
			// Make the vector-field a unit vector-field
			if( false )
			{
#pragma omp parallel for
				for( int t=0 ; t<vectorField.size() ; t++ )
				{
					double len2 = Point2D< double >::Dot( vectorField[t] , embeddingMetric[t]*vectorField[t] );
					if( len2>0 ) vectorField[t] /= (Real)sqrt( len2 );
					else         vectorField[t] *= 0;
				}
			}
			// Normalize the scale of the vector-field
			{
				double norm = 0 , area = 0;
				for( int t=0 ; t<embeddingMetric.size() ; t++ )
				{
					double a = sqrt( embeddingMetric[t].determinant() ) / 2.;
					norm += Point2D< double >::Dot( vectorField[t] , embeddingMetric[t]*vectorField[t] ) * a;
					area += a;
				}
				norm = sqrt( norm / area );
				for( int t=0 ; t<embeddingMetric.size() ; t++ ) vectorField[t] /= norm;
			}
		}
		auto LengthToAnisotropy = [&]( double len )
		{
			// g <- g + gOrtho * anisotropy 
			// 0 -> 0
			// 1 -> 1e5
			// infty -> infty
			return pow( len , AnisoExponent.value ) * 1e5;
		};
		if( !InitializeAnisotropicMetric( mesh , atlasCharts , vectorField , LengthToAnisotropy , parameterMetric ) )
		{
			printf( "[ERROR] Unable to initialize metric \n");
			return 0;
		} 

		std::vector<Point3D<double>> __inputSignal;
		std::vector<double> __texelToCellCoeffs;
		SparseMatrix<double, int> __boundaryCellBasedStiffnessRHSMatrix[3];

		t_begin = clock();
		{
			int ret = 0;
			switch( MatrixQuadrature.value )
			{
			case  1: ret = InitializeMassAndStiffness< 1>( anisoDeepMassCoefficients , anisoDeepStiffnessCoefficients , anisoBoundaryBoundaryMassMatrix , anisoBoundaryBoundaryStiffnessMatrix , anisoBoundaryDeepMassMatrix , anisoBoundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case  3: ret = InitializeMassAndStiffness< 3>( anisoDeepMassCoefficients , anisoDeepStiffnessCoefficients , anisoBoundaryBoundaryMassMatrix , anisoBoundaryBoundaryStiffnessMatrix , anisoBoundaryDeepMassMatrix , anisoBoundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case  6: ret = InitializeMassAndStiffness< 6>( anisoDeepMassCoefficients , anisoDeepStiffnessCoefficients , anisoBoundaryBoundaryMassMatrix , anisoBoundaryBoundaryStiffnessMatrix , anisoBoundaryDeepMassMatrix , anisoBoundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case 12: ret = InitializeMassAndStiffness<12>( anisoDeepMassCoefficients , anisoDeepStiffnessCoefficients , anisoBoundaryBoundaryMassMatrix , anisoBoundaryBoundaryStiffnessMatrix , anisoBoundaryDeepMassMatrix , anisoBoundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case 24: ret = InitializeMassAndStiffness<24>( anisoDeepMassCoefficients , anisoDeepStiffnessCoefficients , anisoBoundaryBoundaryMassMatrix , anisoBoundaryBoundaryStiffnessMatrix , anisoBoundaryDeepMassMatrix , anisoBoundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case 32: ret = InitializeMassAndStiffness<32>( anisoDeepMassCoefficients , anisoDeepStiffnessCoefficients , anisoBoundaryBoundaryMassMatrix , anisoBoundaryBoundaryStiffnessMatrix , anisoBoundaryDeepMassMatrix , anisoBoundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			default: fprintf( stderr , "[ERROR] Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles\n" );
			}
			if( !ret )
			{
				fprintf( stderr , "[ERROR] Failed intialization!\n" );
				return 0;
			}
		}
		if( Verbose.set ) printf( "\tInitialized mass and stiffness: %.2f(s)\n" , double(clock() - t_begin) / CLOCKS_PER_SEC);

		if( UseDirectSolver.set )
		{
			clock_t t_begin;
			t_begin = clock();
			FullMatrixConstruction(hierarchy.gridAtlases[0], anisoDeepMassCoefficients, anisoBoundaryBoundaryMassMatrix, anisoBoundaryDeepMassMatrix, anisotropicMass);
			FullMatrixConstruction(hierarchy.gridAtlases[0], anisoDeepStiffnessCoefficients, anisoBoundaryBoundaryStiffnessMatrix, anisoBoundaryDeepStiffnessMatrix, anisotropicStiffness);
			lineConvolutionMatrix = anisotropicMass * licInterpolationWeight + anisotropicStiffness;
			printf("Assembling matrices =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
		}

		t_begin = clock();
		if(!UpdateLinearSystem( licInterpolationWeight , 1.0 , hierarchy , multigridLineConvolutionCoefficients ,
			anisoDeepMassCoefficients, anisoDeepStiffnessCoefficients,
			anisoBoundaryBoundaryMassMatrix, anisoBoundaryBoundaryStiffnessMatrix,
			anisoBoundaryDeepMassMatrix, anisoBoundaryDeepStiffnessMatrix,
			coarseLineConvolutionSolver, boundaryLineConvolutionSolver,fineLineConvolutionSolver,
			lineConvolutionMatrix, DetailVerbose.set, true, UseDirectSolver.set)){
			printf("ERROR : Failed system update! \n");
			return 0;
		}
		if( Verbose.set ) printf( "\tInitialized multigrid coefficients: %.2f(s)\n" , double(clock() - t_begin) / CLOCKS_PER_SEC);
	}

	//////////////////////////////////// 	Modulation coefficients
	SparseMatrix< double, int > modulationMatrix;
	{
		if (!InitializeMetric(mesh, EMBEDDING_METRIC, atlasCharts, parameterMetric)) {
			printf("ERROR: Unable to initialize metric \n");
			return 0;
		}

		std::vector<Point3D<double>> __inputSignal;
		std::vector<double> __texelToCellCoeffs;
		SparseMatrix<double, int> __boundaryCellBasedStiffnessRHSMatrix[3];

		t_begin = clock();
		{
			int ret = 0;
			switch( MatrixQuadrature.value )
			{
			case 1:
				ret = InitializeMassAndStiffness< 1>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix );
				break;
			case 3:
				ret = InitializeMassAndStiffness< 3>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix );
				break;
			case 6:
				ret = InitializeMassAndStiffness< 6>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix );
				break;
			case 12:
				ret = InitializeMassAndStiffness<12>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix );
				break;
			case 24:
				ret = InitializeMassAndStiffness<24>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix );
				break;
			case 32:
				ret = InitializeMassAndStiffness<32>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix );
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

		if (UseDirectSolver.set) {
			clock_t t_begin;
			t_begin = clock();
			FullMatrixConstruction(hierarchy.gridAtlases[0], deepMassCoefficients, boundaryBoundaryMassMatrix, boundaryDeepMassMatrix, mass);
			FullMatrixConstruction(hierarchy.gridAtlases[0], deepStiffnessCoefficients, boundaryBoundaryStiffnessMatrix, boundaryDeepStiffnessMatrix, stiffness);
			modulationMatrix = mass * sharpeningInterpolationWeight + stiffness;
			printf("Assembling matrices =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
		}

		t_begin = clock();
		if (!UpdateLinearSystem( sharpeningInterpolationWeight , 1.0 , hierarchy , multigridModulationCoefficients ,
			deepMassCoefficients, deepStiffnessCoefficients,
			boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix,
			boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix,
			coarseModulationSolver, boundaryModulationSolver, fineModulationSolver,
			modulationMatrix, DetailVerbose.set, true, UseDirectSolver.set)){
			printf("ERROR : Failed system update! \n");
			return 0;
		}
		if( Verbose.set ) printf( "\tInitialized multigrid coefficients: %.2f(s)\n" , double(clock() - t_begin) / CLOCKS_PER_SEC);
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

	int numTexels = hierarchy.gridAtlases[0].numTexels;
	int numFineNodes = hierarchy.gridAtlases[0].numFineNodes;

	randSignal.resize( textureNodes.size() );

	for ( int i=0 ; i<randSignal.size() ; i++ )
	{
		Point3D< float > randomColor = HSV2RGB( double( rand() ) / double(RAND_MAX), 1 , 1 );
		randSignal[i] = Point3D< Real >( randomColor[0] , randomColor[1] , randomColor[2] );
	}

	for( int i=0 ; i<multigridLineConvolutionVariables[0].x.size() ; i++) multigridLineConvolutionVariables[0].x[i] *= 0;
	for( int i=0 ; i<multigridModulationVariables[0].x.size() ; i++) multigridModulationVariables[0].x[i] *= 0;

	if( UseDirectSolver.set ) ComputeExactSolution();
	else for( int i=0 ; i<multigridModulationVariables[0].x.size() ; i++) multigridModulationVariables[0].x[i] = randSignal[i];

	return 1;
}

template< class Real >
void LineConvolution< Real >::InitializeVisualization( const int width , const int height )
{
	int tCount = (int)mesh.triangles.size();

	visualization.triangles.resize( tCount );
	visualization.vertices.resize( 3*tCount );
	visualization.colors.resize( 3*tCount , Point3D< double >( 0.75 , 0.75 , 0.75 ) );
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

	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'f' , "lic interpolation weight" , "LIC Interpolation Weight" , LICInterpolationWeightCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'g' , "sharpening gradient modulation" , "Sharpening Gradient Modulation" , SharpeningGradientModulationCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'y' , "sharpening interpolation weight" , "Sharpening Interpolation Weight" , SharpeningInterpolationWeightCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 's' , "export texture" , "Output Texture" , ExportTextureCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , ' ' , "toggle update" , ToggleUpdateCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , '+' , "increment update" , IncrementUpdateCallBack ) );

	visualization.UpdateVertexBuffer();
	visualization.UpdateFaceBuffer();
	visualization.UpdateTextureBuffer();

	UpdateOutputBuffer( multigridModulationVariables[0].x );
}

template<class Real>
int LineConvolution<Real>::Init( void )
{
	levels = Levels.value;
	textureWidth = Width.value;
	textureHeight = Height.value;
	sharpeningGradientModulation = SharpeningGradientModulation.value;
	sharpeningInterpolationWeight = SharpeningInterpolationWeight.value;
	licInterpolationWeight = LICInterpolationWeight.value;

	if( !ReadTexturedMesh( mesh , Input.value , NULL , DetailVerbose.set ) )
	{
		printf("Unable to read mesh data\n");
		return 0;
	}
	if (1) for (int i = 0; i < mesh.textureCoordinates.size(); i++)mesh.textureCoordinates[i][1] = 1.0 - mesh.textureCoordinates[i][1];

	if( RandomJitter.set )
	{
		srand( time( NULL ) );
		std::vector<Point2D < double >>randomOffset( mesh.vertices.size() );
		double jitterScale = 1e-3 / double(std::max<int>(textureWidth, textureHeight));
		for( int i=0 ; i<randomOffset.size() ; i++ ) randomOffset[i] = Point2D < double >(1.0 - 2.0 * double(rand()) / double(RAND_MAX), 1.0 - 2.0 *  double(rand()) / double(RAND_MAX))*jitterScale;
		for( int i=0 ; i<mesh.triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ ) mesh.textureCoordinates[ 3*i+k ] += randomOffset[ mesh.triangles[i][k] ];
	}

	ComputePadding( padding , textureWidth , textureHeight , mesh.textureCoordinates , DetailVerbose.set );
	if( padding.nonTrivial )
	{
		PaddTextureCoordinates(padding, textureWidth, textureHeight, mesh.textureCoordinates);
		textureWidth += (padding.left + padding.right);
		textureHeight += (padding.bottom + padding.top);
	}

	//Define centroid and scale for visualization
	Point3D< double > centroid(0.f, 0.f, 0.f);
	for( int i=0 ; i<mesh.vertices.size() ;  i++ ) centroid += mesh.vertices[i];
	centroid /= double(mesh.vertices.size());
	double radius = 0.f;
	for( int i=0 ; i< mesh.vertices.size() ; i++ ) radius = std::max<double>(radius, Point3D<double>::Length(mesh.vertices[i] - centroid));
	for( int i=0 ; i< mesh.vertices.size() ; i++ ) mesh.vertices[i] = (mesh.vertices[i] - centroid) / radius;

	clock_t t = clock();
	if( !InitializeSystem( textureWidth , textureHeight , radius ) ){ printf("Unable to initialize system\n") ; return 0; }
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

	outputBuffer = new unsigned char[ textureHeight*textureWidth*3 ];
	memset( outputBuffer , 204 , textureHeight*textureWidth*3*sizeof( unsigned char ) );

	return 1;
}

template<class Real>
int _main(int argc, char* argv[])
{
	if( !LineConvolution< Real >::Init() ) return 0;
	if( !Output.set )
	{
		glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
		LineConvolution< Real >::visualization.displayMode = DisplayMode.value;
		if     ( DisplayMode.value==ONE_REGION_DISPLAY ) LineConvolution<Real>::visualization.screenWidth =  800 , LineConvolution<Real>::visualization.screenHeight = 800;
		else if( DisplayMode.value==TWO_REGION_DISPLAY ) LineConvolution<Real>::visualization.screenWidth = 1440 , LineConvolution<Real>::visualization.screenHeight = 720;
		glutInitWindowSize( LineConvolution< Real >::visualization.screenWidth , LineConvolution< Real >::visualization.screenHeight );
		glutInit( &argc , argv );
		char windowName[1024];
		sprintf( windowName , "Line Integral Convolution" );
		glutCreateWindow( windowName );
		if( glewInit()!=GLEW_OK ) fprintf( stderr , "[ERROR] glewInit failed\n" ) , exit(0);
		glutDisplayFunc( LineConvolution< Real >::Display );
		glutReshapeFunc( LineConvolution< Real >::Reshape );
		glutMouseFunc( LineConvolution< Real >::MouseFunc );
		glutMotionFunc( LineConvolution< Real >::MotionFunc );
		glutKeyboardFunc( LineConvolution< Real >::KeyboardFunc );
		if( !UseDirectSolver.set )glutIdleFunc( LineConvolution< Real >::Idle );
		if( CameraConfig.set ) LineConvolution< Real >::visualization.ReadSceneConfigurationCallBack( &LineConvolution< Real >::visualization , CameraConfig.value );
		LineConvolution< Real >::InitializeVisualization( LineConvolution< Real >::textureWidth , LineConvolution< Real >::textureHeight );
		glutMainLoop();
	}
	else
	{
		if( UseDirectSolver.set ) LineConvolution< Real >::ComputeExactSolution();
		else for( int i=0 ; i<OutputVCycles.value ; i++ ) LineConvolution< Real >::UpdateSolution();
		LineConvolution< Real >::SetOutputBuffer( LineConvolution< Real >::multigridModulationVariables[0].x );
		LineConvolution< Real >::ExportTextureCallBack( &LineConvolution<Real>::visualization , Output.value );
	}
	return 1;
}

int main(int argc, char* argv[])
{
	cmdLineParse(argc - 1, argv + 1, params);
	if( !Input.set ) { ShowUsage(argv[0]); return EXIT_FAILURE; }
	omp_set_num_threads( Threads.value );
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
	if( Double.set ) _main< double >( argc , argv );
	else             _main< float  >( argc , argv );
	return 0;
}
