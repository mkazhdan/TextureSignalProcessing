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
#include <Src/SimpleMesh.h>
#include <Src/Basis.h>
#include <Misha/FEM.h>
#include <Src/Solver.h>
#include <Src/Hierarchy.h>
#include <Src/MassAndStiffness.h>
#include <Src/Padding.h>
#include <Src/StitchingVisualization.h>

cmdLineParameterArray< char * , 3 > In( "in" );
cmdLineParameter< char* > Output("out");
cmdLineParameter< int   > OutputVCycles("outVCycles", 6);
cmdLineParameter< float > InterpolationWeight("interpolation", 1e2);
cmdLineParameter< int   > Threads("threads", omp_get_num_procs());
cmdLineParameter< int   > Levels("levels", 4);
cmdLineParameter< int   > MatrixQuadrature("mQuadrature", 6);

cmdLineParameter< int   > MultigridBlockHeight("mBlockH", 16);
cmdLineParameter< int   > MultigridBlockWidth("mBlockW", 128);
cmdLineParameter< int   > MultigridPaddedHeight("mPadH", 0);
cmdLineParameter< int   > MultigridPaddedWidth("mPadW", 2);

cmdLineReadable RandomJitter("jitter");
cmdLineParameter< char* > CameraConfig("camera");
cmdLineReadable UseDirectSolver("useDirectSolver");
cmdLineReadable Verbose("verbose");
cmdLineReadable NoHelp("noHelp");
cmdLineReadable DetailVerbose("detail");
cmdLineReadable Double("double");
cmdLineReadable MultiInput( "multi" );
cmdLineReadable* params[] =
{
	&In , &Output , &InterpolationWeight , &CameraConfig , &Levels , &UseDirectSolver , &Threads, &Verbose ,
	&DetailVerbose , &MultigridBlockHeight , &MultigridBlockWidth , &MultigridPaddedHeight , &MultigridPaddedWidth , &RandomJitter ,
	&Double , &MatrixQuadrature , &OutputVCycles , &NoHelp , &MultiInput ,
	NULL
};

void ShowUsage(const char* ex)
{
	printf("Usage %s:\n", ex);

	printf( "\t --%s <input mesh, texels , and mask>\n" , In.name );
	printf( "\t[--%s <output texture>]\n" , Output.name );
	printf( "\t[--%s <output v-cycles>=%d]\n" , OutputVCycles.name , OutputVCycles.value );
	printf( "\t[--%s <interpolation weight>=%f]\n" , InterpolationWeight.name , InterpolationWeight.value );
	printf( "\t[--%s <system matrix quadrature points per triangle>=%d]\n" , MatrixQuadrature.name, MatrixQuadrature.value );
	printf( "\t[--%s]\n" , UseDirectSolver.name );
	printf( "\t[--%s]\n" , MultiInput.name );
	printf( "\t[--%s]\n" , RandomJitter.name );
	printf( "\t[--%s]\n" , Verbose.name );

	printf( "\t[--%s <camera configuration file>]\n" , CameraConfig.name );
	printf( "\t[--%s <hierarchy levels>=%d]\n" , Levels.name, Levels.value );
	printf( "\t[--%s <threads>=%d]\n", Threads.name , Threads.value );
	printf( "\t[--%s <multigrid block width>=%d]\n" , MultigridBlockWidth.name , MultigridBlockWidth.value );
	printf( "\t[--%s <multigrid block height>=%d]\n" , MultigridBlockHeight.name , MultigridBlockHeight.value );
	printf( "\t[--%s <multigrid padded width>=%d]\n" , MultigridPaddedWidth.name , MultigridPaddedWidth.value );
	printf( "\t[--%s <multigrid padded height>=%d]\n" , MultigridPaddedHeight.name , MultigridPaddedHeight.value );
	printf( "\t[--%s]\n", DetailVerbose.name );
	printf( "\t[--%s]\n" , NoHelp.name );
}

template<class Real>
class Stitching
{
public:
	static int inputMode;
	static TexturedMesh mesh;
	static int textureWidth;
	static int textureHeight;
	static double interpolationWeight;
	static int levels;

	static HierarchicalSystem hierarchy;
	static bool rhsUpdated;
	static bool positiveModulation;

	// Single input mode
	static Image<int> inputMask;
	static Image<Point3D<float>> inputComposition;

	// Multiple input mode
	static int numTextures;
	static std::vector<Image<float>> inputConfidence;
	static std::vector<Image<Point3D<float>>> inputTextures;
	static std::vector<std::vector<Point3D<Real>>> partialTexelValues;
	static std::vector<std::vector<Point3D<Real>>> partialEdgeValues;

	static Image<Point3D<float>> filteredTexture;
	// UI
	static char interpolationStr[1024];
	static char referenceTextureStr[1024];

	static std::vector<Point3D<float>>textureNodePositions;
	static std::vector<Point3D<float>>textureEdgePositions;

	static std::vector<AtlasChart> atlasCharts;

	static std::vector< BilinearElementIndex > bilinearElementIndices;

	static std::vector<TextureNodeInfo> textureNodes;

	static SparseMatrix<double, int> mass;
	static SparseMatrix<double, int> stiffness;
	static SparseMatrix< double, int > stitchingMatrix;

	static std::vector< Point3D< Real > > texelMass;
	static std::vector< Point3D< Real > > texelDivergence;

	static std::vector< MultigridLevelCoefficients< Real > > multigridStitchingCoefficients;
	static std::vector< MultigridLevelVariables< Point3D< Real > > > multigridStitchingVariables;
	static std::vector< MultigridLevelIndices< Real > > multigridIndices;

#if defined( USE_CHOLMOD )
	typedef  std::vector< CholmodCholeskySolver< Real, 3 > > BoundarySolverType;
	typedef  CholmodCholeskySolver< Real, 3 > CoarseSolverType;
	typedef  CholmodCholeskySolver< Real, 3 > DirectSolverType;
#elif defined( USE_EIGEN_SIMPLICIAL )
	typedef  std::vector< EigenCholeskySolver< Real, 3 > > BoundarySolverType;
	typedef  EigenCholeskySolver< Real, 3 > CoarseSolverType;
	typedef  EigenCholeskySolver< Real, 3 > DirectSolverType;
#elif defined( USE_EIGEN_PARDISO )
	typedef  std::vector< EigenPardisoSolver< Real, 3 > > BoundarySolverType;
	typedef  EigenPardisoSolver< Real, 3 > CoarseSolverType;
	typedef  EigenPardisoSolver< Real, 3 > DirectSolverType;
#else
#error "[ERROR] No solver defined!"
#endif

	static BoundarySolverType boundarySolver;
	static CoarseSolverType coarseSolver;
	static DirectSolverType directSolver;

	static std::unordered_map<unsigned long long, int> edgeIndex;
	static std::vector<std::pair<int,int>> edgePairs;

	static SparseMatrix< Real, int > boundaryDivergenceMatrix;

	static std::vector< Real > deepDivergenceCoefficients;
	static std::vector<DivegenceRasterLine> divergenceRasterLines;

	static std::vector<bool> unobservedTexel;
	static std::vector<Point3D<Real>> texelValues;
	static std::vector<Point3D<Real>> edgeValues;

	// Linear Operators
	static std::vector<double> deepMassCoefficients;
	static std::vector<double> deepStiffnessCoefficients;
	static SparseMatrix<double, int> boundaryBoundaryMassMatrix;
	static SparseMatrix<double, int> boundaryBoundaryStiffnessMatrix;
	static SparseMatrix<double, int> boundaryDeepMassMatrix;
	static SparseMatrix<double, int> boundaryDeepStiffnessMatrix;

	// Stitching UI
	static int textureIndex;

	static int steps;
	static char stepsString[];

	static Padding padding;

	// Visulization
	static StitchingVisualization visualization;
	static int updateCount;

	static void ToggleForwardReferenceTextureCallBack(Visualization* v, const char* prompt);
	static void ToggleBackwardReferenceTextureCallBack(Visualization* v, const char* prompt);
	static void ToggleUpdateCallBack(Visualization* v, const char* prompt);
	static void IncrementUpdateCallBack(Visualization* v, const char* prompt);
	static void ExportTextureCallBack(Visualization* v, const char* prompt);

	static void InterpolationWeightCallBack(Visualization* v, const char* prompt);

	static void LoadImages();
	static void ParseImages();
	static void SolveSytem();
	static int Init();
	static void InitializeVisualization();
	static int UpdateSolution(bool verbose = false, bool detailVerbose = false);
	static void ComputeExactSolution(bool verbose = false);
	static int InitializeSystem(const int width, const int height);
	static int _InitializeSystem(std::vector<std::vector<SquareMatrix<double, 2>>>& parameterMetric, BoundaryProlongationData& boundaryProlongation);

	static void UpdateFilteredColorTexture(const std::vector< Point3D< Real > >& solution);
	static void UpdateFilteredTexture(const std::vector< Point3D< Real > >& solution);

	static void Display(void) { visualization.Display(); }
	static void MouseFunc(int button, int state, int x, int y);
	static void MotionFunc(int x, int y);
	static void Reshape(int w, int h) { visualization.Reshape(w, h); }
	static void KeyboardFunc(unsigned char key, int x, int y) { visualization.KeyboardFunc(key, x, y); }
	static void Idle(void);
};
template< class Real > char Stitching< Real >::referenceTextureStr[1024];
template< class Real > char Stitching< Real >::interpolationStr[1024];


template<class Real> int														Stitching<Real>::inputMode;
template<class Real> TexturedMesh												Stitching<Real>::mesh;
template<class Real> int														Stitching<Real>::textureWidth;
template<class Real> int														Stitching<Real>::textureHeight;
template<class Real> StitchingVisualization										Stitching<Real>::visualization;
template<class Real> SparseMatrix<double, int>									Stitching<Real>::mass;
template<class Real> SparseMatrix<double, int>									Stitching<Real>::stiffness;
template<class Real> SparseMatrix<double, int>									Stitching<Real>::stitchingMatrix;

template<class Real> double														Stitching<Real>::interpolationWeight;

template<class Real> std::vector<TextureNodeInfo>								Stitching<Real>::textureNodes;
template<class Real> std::vector< BilinearElementIndex >						Stitching<Real>::bilinearElementIndices;

template< class Real > int														Stitching< Real >::steps;
template< class Real > char														Stitching< Real >::stepsString[1024];
template<class Real> int														Stitching<Real>::levels;
template<class Real> HierarchicalSystem											Stitching<Real>::hierarchy;


template<class Real> bool Stitching<Real>::rhsUpdated = true;
template<class Real> bool Stitching<Real>::positiveModulation = true;

template<class Real> Image<Point3D<float>>												Stitching<Real>::filteredTexture;


template<class Real> Image<int>														    Stitching<Real>::inputMask;
template<class Real> Image<Point3D<float>>												Stitching<Real>::inputComposition;

template<class Real> int																Stitching<Real>::numTextures;
template<class Real> std::vector<Image<float>>											Stitching<Real>::inputConfidence;
template<class Real> std::vector<Image<Point3D<float>>>									Stitching<Real>::inputTextures;
template<class Real> std::vector< Point3D< Real > >										Stitching<Real>::texelMass;
template<class Real> std::vector< Point3D< Real > >										Stitching<Real>::texelDivergence;

template<class Real> std::vector<MultigridLevelCoefficients<Real>>						Stitching<Real>::multigridStitchingCoefficients;
template<class Real> std::vector< MultigridLevelVariables< Point3D< Real > > >			Stitching<Real>::multigridStitchingVariables;
template<class Real> std::vector<MultigridLevelIndices<Real>>							Stitching<Real>::multigridIndices;


template<class Real> typename Stitching<Real>::BoundarySolverType					Stitching<Real>::boundarySolver;
template<class Real> typename Stitching<Real>::CoarseSolverType						Stitching<Real>::coarseSolver;
template<class Real> typename Stitching<Real>::DirectSolverType						Stitching<Real>::directSolver;

template<class Real> std::vector<AtlasChart>											Stitching<Real>::atlasCharts;

template<class Real> std::vector<Point3D<float>>										Stitching<Real>::textureNodePositions;
template<class Real> std::vector<Point3D<float>>										Stitching<Real>::textureEdgePositions;

template<class Real> Padding															Stitching<Real>::padding;


template<class Real>  std::vector<double>												Stitching<Real>::deepMassCoefficients;
template<class Real>  std::vector<double>												Stitching<Real>::deepStiffnessCoefficients;
template<class Real>  SparseMatrix<double, int>											Stitching<Real>::boundaryBoundaryMassMatrix;
template<class Real>  SparseMatrix<double, int>											Stitching<Real>::boundaryBoundaryStiffnessMatrix;
template<class Real>  SparseMatrix<double, int>											Stitching<Real>::boundaryDeepMassMatrix;
template<class Real>  SparseMatrix<double, int>											Stitching<Real>::boundaryDeepStiffnessMatrix;
template< class Real > int																Stitching< Real >::updateCount = -1;

template<class Real>  std::unordered_map<unsigned long long, int>						Stitching< Real >::edgeIndex;
template<class Real>  std::vector<std::pair<int, int>>									Stitching< Real >::edgePairs;

template<class Real>  SparseMatrix< Real, int >											Stitching< Real >::boundaryDivergenceMatrix;
template<class Real>  std::vector< Real >												Stitching< Real >::deepDivergenceCoefficients;
template<class Real>  std::vector<DivegenceRasterLine>  								Stitching< Real >::divergenceRasterLines;

template<class Real>  std::vector<bool>													Stitching< Real >::unobservedTexel;
template<class Real>  std::vector<Point3D<Real>>										Stitching< Real >::texelValues;
template<class Real>  std::vector<Point3D<Real>>										Stitching< Real >::edgeValues;
template<class Real>  std::vector<std::vector<Point3D<Real>>>							Stitching< Real >::partialTexelValues;
template<class Real>  std::vector<std::vector<Point3D<Real>>>							Stitching< Real >::partialEdgeValues;
template<class Real> int																Stitching< Real >::textureIndex = 0;

template<class Real>
void Stitching< Real >::UpdateFilteredColorTexture(const std::vector< Point3D< Real > > & solution)
{
#pragma omp parallel for
	for (int i = 0; i < textureNodes.size(); i++) {
		int ci = textureNodes[i].ci;
		int cj = textureNodes[i].cj;
		int offset = 3 * (textureWidth*cj + ci);
		for( int c=0 ; c<3 ; c++ )
		{
			double value = std::min<double>(1.0, std::max<double>(0, solution[i][c]));
			visualization.colorTextureBuffer[offset + c] = (unsigned char)(value*255.0);
		}
	}
}

template< class Real >
void Stitching< Real >::UpdateFilteredTexture(const std::vector< Point3D< Real > >& solution)
{
#pragma omp parallel for
	for (int i = 0; i<textureNodes.size(); i++)
	{
		int ci = textureNodes[i].ci, cj = textureNodes[i].cj;
		filteredTexture(ci, cj) = Point3D< float >(solution[i][0], solution[i][1], solution[i][2]);
	}
}

template< class Real >
void Stitching< Real >::Idle(void)
{
	float radius = 0.03;
	float radiusSquared = radius * radius;
	if( inputMode==MULTIPLE_INPUT_MODE )
	{
		if (visualization.isBrushActive) {
			Point3D< float > selectedPoint;
			bool validSelection = false;
			if (visualization.showMesh) validSelection = visualization.select(visualization.diskX, visualization.diskY, selectedPoint);
			if (validSelection){
#pragma omp parallel for
				for (int i = 0; i < textureNodePositions.size(); i++)if (Point3D< float >::SquareNorm(textureNodePositions[i] - selectedPoint) < radiusSquared) texelValues[i] = partialTexelValues[textureIndex][i];
#pragma omp parallel for
				for (int i = 0; i < textureEdgePositions.size(); i++)if (Point3D< float >::SquareNorm(textureEdgePositions[i] - selectedPoint) < radiusSquared) edgeValues[i] = partialEdgeValues[textureIndex][i];
				rhsUpdated = false;
			}
		}
	}

	if (updateCount && !UseDirectSolver.set && !visualization.promptCallBack){
		if (!UpdateSolution()) fprintf(stderr, "[ERROR] Updated solution failed!\n");
		UpdateFilteredColorTexture(multigridStitchingVariables[0].x);
		visualization.UpdateColorTextureBuffer();
		if (updateCount>0) updateCount--;
		steps++;
		sprintf(stepsString, "Steps: %d", steps);
	}

	glutPostRedisplay();
}

template< class Real >
void Stitching< Real >::MouseFunc(int button, int state, int x, int y)
{
	visualization.newX = x; visualization.newY = y;
	visualization.rotating = visualization.scaling = visualization.panning = false;
	visualization.isBrushActive = false;

	if (state == GLUT_DOWN && glutGetModifiers() & GLUT_ACTIVE_SHIFT)
	{
		visualization.isBrushActive = true;
		visualization.diskX = x;
		visualization.diskY = y;

		if (button == GLUT_RIGHT_BUTTON) positiveModulation = true;
		else if (button == GLUT_LEFT_BUTTON) positiveModulation = false;
	}
	else
	{
		if (visualization.showMesh)
		{
			visualization.newX = x; visualization.newY = y;

			visualization.rotating = visualization.scaling = visualization.panning = false;
			if ((button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) && glutGetModifiers() & GLUT_ACTIVE_CTRL) visualization.panning = true;
			else if (button == GLUT_LEFT_BUTTON) visualization.rotating = true;
			else if (button == GLUT_RIGHT_BUTTON) visualization.scaling = true;
		}
	}

	glutPostRedisplay();
}

template<class Real>
void Stitching<Real>::MotionFunc(int x, int y)
{
	if (visualization.isBrushActive)
	{
		visualization.diskX = x;
		visualization.diskY = y;
	}
	else
	{
		if (visualization.showMesh)
		{
			visualization.oldX = visualization.newX, visualization.oldY = visualization.newY, visualization.newX = x, visualization.newY = y;
			int screenSize = std::min< int >(visualization.screenWidth, visualization.screenHeight);
			float rel_x = (visualization.newX - visualization.oldX) / (float)screenSize * 2;
			float rel_y = (visualization.newY - visualization.oldY) / (float)screenSize * 2;

			float pRight = rel_x * visualization.zoom, pUp = -rel_y * visualization.zoom;
			float pForward = rel_y * visualization.zoom;
			float rRight = -rel_y, rUp = -rel_x;

			if (visualization.rotating) visualization.camera.rotateUp(-rUp), visualization.camera.rotateRight(-rRight);
			else if (visualization.scaling) visualization.camera.translate(visualization.camera.forward*pForward);
			else if (visualization.panning) visualization.camera.translate(-(visualization.camera.right*pRight + visualization.camera.up*pUp));
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

template< class Real > void Stitching< Real >::ToggleForwardReferenceTextureCallBack(Visualization* v, const char* ) {
	textureIndex = (textureIndex + 1) % numTextures;
	visualization.referenceIndex = textureIndex;
	sprintf(referenceTextureStr, "Reference Texture: %02d of %02d \n", textureIndex, numTextures);
	glutPostRedisplay();
}
template< class Real > void Stitching< Real >::ToggleBackwardReferenceTextureCallBack(Visualization* v, const char* ) {
	textureIndex = (textureIndex + numTextures - 1) % numTextures;
	visualization.referenceIndex = textureIndex;
	sprintf(referenceTextureStr, "Reference Texture: %02d of %02d \n", textureIndex, numTextures);
	glutPostRedisplay();
}

template< class Real > void Stitching< Real >::ToggleUpdateCallBack(Visualization* v, const char* prompt)
{
	if (updateCount) updateCount = 0;
	else              updateCount = -1;
}
template< class Real > void Stitching< Real >::IncrementUpdateCallBack(Visualization* v, const char* prompt)
{
	if (updateCount<0) updateCount = 1;
	else updateCount++;
}

template<class Real>
void Stitching<Real>::ExportTextureCallBack(Visualization* v, const char* prompt)
{
	UpdateFilteredTexture(multigridStitchingVariables[0].x);
	Image< Point3D< float > > outputTexture = filteredTexture;
	if (padding.nonTrivial) UnpadImage(padding, outputTexture);
	outputTexture.write(prompt);
}

template< class Real >
void  Stitching< Real >::InterpolationWeightCallBack(Visualization* v, const char* prompt)
{
	interpolationWeight = atof(prompt);
	if (UseDirectSolver.set) stitchingMatrix = mass*interpolationWeight + stiffness;
	clock_t t_begin = clock();
	if (!UpdateLinearSystem(interpolationWeight, 1.0, hierarchy, multigridStitchingCoefficients,
		deepMassCoefficients, deepStiffnessCoefficients,
		boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix,
		boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix,
		coarseSolver, boundarySolver, directSolver,
		stitchingMatrix, DetailVerbose.set, false, UseDirectSolver.set)) {
		fprintf(stderr, "[ERROR] Failed system update!\n");
	}
	if (Verbose.set) printf("\tInitialized multigrid coefficients: %.2f(s)\n", double(clock() - t_begin) / CLOCKS_PER_SEC);

#pragma omp parallel for
	for (int i = 0; i<multigridStitchingVariables[0].rhs.size(); i++) multigridStitchingVariables[0].rhs[i] = texelMass[i] * interpolationWeight + texelDivergence[i];

	if (UseDirectSolver.set) ComputeExactSolution(Verbose.set);
	else for (int i = 0; i < OutputVCycles.value; i++)VCycle(multigridStitchingVariables, multigridStitchingCoefficients, multigridIndices, boundarySolver, coarseSolver, false, false);

	UpdateFilteredColorTexture(multigridStitchingVariables[0].x);
	visualization.UpdateColorTextureBuffer();
	sprintf(interpolationStr, "Interpolation weight: %e\n", interpolationWeight);
}


template<class Real>
void Stitching<Real>::ComputeExactSolution(bool verbose)
{
	clock_t t_begin = clock();
	solve(directSolver, multigridStitchingVariables[0].x, multigridStitchingVariables[0].rhs);
	if (verbose) printf("Solving time =  %.4f\n", double(clock() - t_begin) / CLOCKS_PER_SEC);
}

template<class Real>
int Stitching<Real>::UpdateSolution(bool verbose, bool detailVerbose)
{
	if (!rhsUpdated)
	{
		int numTexels = (int)multigridStitchingVariables[0].rhs.size();
		clock_t p_begin = clock();

		MultiplyBySystemMatrix_NoReciprocals(deepMassCoefficients, boundaryDeepMassMatrix, boundaryBoundaryMassMatrix, hierarchy.gridAtlases[0].boundaryGlobalIndex, hierarchy.gridAtlases[0].rasterLines, texelValues, texelMass);
		ComputeDivergence(edgeValues, texelDivergence, deepDivergenceCoefficients, boundaryDivergenceMatrix, divergenceRasterLines);

#pragma omp parallel for
		for (int i = 0; i <textureNodes.size(); i++) multigridStitchingVariables[0].rhs[i] = texelMass[i] * interpolationWeight + texelDivergence[i];

		if (verbose) printf("RHS update time %.4f  \n", double(clock() - p_begin) / CLOCKS_PER_SEC);
		rhsUpdated = true;
	}

	VCycle(multigridStitchingVariables, multigridStitchingCoefficients, multigridIndices, boundarySolver, coarseSolver, verbose, detailVerbose);

	return 1;
}

template<>
int Stitching< float >::_InitializeSystem(std::vector<std::vector<SquareMatrix<double, 2>>>& parameterMetric, BoundaryProlongationData& boundaryProlongation)
{
	// Unused parameters
	std::vector<Point3D<double>> inputSignal;
	std::vector<double> texelToCellCoeffs;
	SparseMatrix<double, int> boundaryCellBasedStiffnessRHSMatrix[3];

	SparseMatrix< double, int > _boundaryDivergenceMatrix;
	std::vector< double >  _deepDivergenceCoefficients;
	clock_t t_begin = clock();
	{
		int ret = 0;
		switch (MatrixQuadrature.value)
		{
		case 1:
			ret = InitializeMassAndStiffness< 1>(deepMassCoefficients, deepStiffnessCoefficients, boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix, boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix, hierarchy, parameterMetric, atlasCharts, boundaryProlongation, false, inputSignal, texelToCellCoeffs, boundaryCellBasedStiffnessRHSMatrix, true, edgeIndex, _boundaryDivergenceMatrix, _deepDivergenceCoefficients);
			break;
		case 3:
			ret = InitializeMassAndStiffness< 3>(deepMassCoefficients, deepStiffnessCoefficients, boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix, boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix, hierarchy, parameterMetric, atlasCharts, boundaryProlongation, false, inputSignal, texelToCellCoeffs, boundaryCellBasedStiffnessRHSMatrix, true, edgeIndex, _boundaryDivergenceMatrix, _deepDivergenceCoefficients);
			break;
		case 6:
			ret = InitializeMassAndStiffness< 6>(deepMassCoefficients, deepStiffnessCoefficients, boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix, boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix, hierarchy, parameterMetric, atlasCharts, boundaryProlongation, false, inputSignal, texelToCellCoeffs, boundaryCellBasedStiffnessRHSMatrix, true, edgeIndex, _boundaryDivergenceMatrix, _deepDivergenceCoefficients);
			break;
		case 12:
			ret = InitializeMassAndStiffness<12>(deepMassCoefficients, deepStiffnessCoefficients, boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix, boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix, hierarchy, parameterMetric, atlasCharts, boundaryProlongation, false, inputSignal, texelToCellCoeffs, boundaryCellBasedStiffnessRHSMatrix, true, edgeIndex, _boundaryDivergenceMatrix, _deepDivergenceCoefficients);
			break;
		case 24:
			ret = InitializeMassAndStiffness<24>(deepMassCoefficients, deepStiffnessCoefficients, boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix, boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix, hierarchy, parameterMetric, atlasCharts, boundaryProlongation, false, inputSignal, texelToCellCoeffs, boundaryCellBasedStiffnessRHSMatrix, true, edgeIndex, _boundaryDivergenceMatrix, _deepDivergenceCoefficients);
			break;
		case 32:
			ret = InitializeMassAndStiffness<32>(deepMassCoefficients, deepStiffnessCoefficients, boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix, boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix, hierarchy, parameterMetric, atlasCharts, boundaryProlongation, false, inputSignal, texelToCellCoeffs, boundaryCellBasedStiffnessRHSMatrix, true, edgeIndex, _boundaryDivergenceMatrix, _deepDivergenceCoefficients);
			break;
		default:
			fprintf(stderr, "[ERROR] Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles\n");
		}
		if (!ret)
		{
			fprintf(stderr, "[ERROR] Failed intialization!\n");
			return 0;
		}
	}

	if (Verbose.set) printf("\tInitialized mass, stiffness, and divergence: %.2f(s)\n", double(clock() - t_begin) / CLOCKS_PER_SEC);
	boundaryDivergenceMatrix = _boundaryDivergenceMatrix;
	deepDivergenceCoefficients.resize(_deepDivergenceCoefficients.size());
	for (int i = 0; i < _deepDivergenceCoefficients.size(); i++) deepDivergenceCoefficients[i] = static_cast<float>(_deepDivergenceCoefficients[i]);

	return 1;
}
template<>
int Stitching< double >::_InitializeSystem(std::vector<std::vector<SquareMatrix<double, 2>>>& parameterMetric, BoundaryProlongationData& boundaryProlongation)
{
	// Unused parameters
	std::vector<Point3D<double>> inputSignal;
	std::vector<double> texelToCellCoeffs;
	SparseMatrix<double, int> boundaryCellBasedStiffnessRHSMatrix[3];

	clock_t t_begin = clock();
	{
		int ret = 0;
		switch (MatrixQuadrature.value)
		{
		case 1:
			ret = InitializeMassAndStiffness< 1>(deepMassCoefficients, deepStiffnessCoefficients, boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix, boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix, hierarchy, parameterMetric, atlasCharts, boundaryProlongation, false, inputSignal, texelToCellCoeffs, boundaryCellBasedStiffnessRHSMatrix, true, edgeIndex, boundaryDivergenceMatrix, deepDivergenceCoefficients);
			break;
		case 3:
			ret = InitializeMassAndStiffness< 3>(deepMassCoefficients, deepStiffnessCoefficients, boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix, boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix, hierarchy, parameterMetric, atlasCharts, boundaryProlongation, false, inputSignal, texelToCellCoeffs, boundaryCellBasedStiffnessRHSMatrix, true, edgeIndex, boundaryDivergenceMatrix, deepDivergenceCoefficients);
			break;
		case 6:
			ret = InitializeMassAndStiffness< 6>(deepMassCoefficients, deepStiffnessCoefficients, boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix, boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix, hierarchy, parameterMetric, atlasCharts, boundaryProlongation, false, inputSignal, texelToCellCoeffs, boundaryCellBasedStiffnessRHSMatrix, true, edgeIndex, boundaryDivergenceMatrix, deepDivergenceCoefficients);
			break;
		case 12:
			ret = InitializeMassAndStiffness<12>(deepMassCoefficients, deepStiffnessCoefficients, boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix, boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix, hierarchy, parameterMetric, atlasCharts, boundaryProlongation, false, inputSignal, texelToCellCoeffs, boundaryCellBasedStiffnessRHSMatrix, true, edgeIndex, boundaryDivergenceMatrix, deepDivergenceCoefficients);
			break;
		case 24:
			ret = InitializeMassAndStiffness<24>(deepMassCoefficients, deepStiffnessCoefficients, boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix, boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix, hierarchy, parameterMetric, atlasCharts, boundaryProlongation, false, inputSignal, texelToCellCoeffs, boundaryCellBasedStiffnessRHSMatrix, true, edgeIndex, boundaryDivergenceMatrix, deepDivergenceCoefficients);
			break;
		case 32:
			ret = InitializeMassAndStiffness<32>(deepMassCoefficients, deepStiffnessCoefficients, boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix, boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix, hierarchy, parameterMetric, atlasCharts, boundaryProlongation, false, inputSignal, texelToCellCoeffs, boundaryCellBasedStiffnessRHSMatrix, true, edgeIndex, boundaryDivergenceMatrix, deepDivergenceCoefficients);
			break;
		default:
			fprintf(stderr, "[ERROR] Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles\n");
		}
		if (!ret)
		{
			fprintf(stderr, "[ERROR] Failed intialization!\n");
			return 0;
		}
	}

	if (Verbose.set) printf("\tInitialized mass and stiffness: %.2f(s)\n", double(clock() - t_begin) / CLOCKS_PER_SEC);
	return 1;
}
template<class Real>
int Stitching<Real>::InitializeSystem(const int width, const int height)
{
	clock_t t_begin;

	if (Verbose.set) t_begin = clock();
	MultigridBlockInfo multigridBlockInfo(MultigridBlockWidth.value, MultigridBlockHeight.value, MultigridPaddedWidth.value, MultigridPaddedHeight.value, 0);
	if (!InitializeHierarchy(mesh, width, height, levels, textureNodes, bilinearElementIndices, hierarchy, atlasCharts, multigridBlockInfo, true, DetailVerbose.set))
	{
		printf("ERROR : Failed intialization! \n");
		return 0;
	}
	if (Verbose.set) printf("\tInitialized hierarchy: %.2f(s)\n", double(clock() - t_begin) / CLOCKS_PER_SEC);

	BoundaryProlongationData boundaryProlongation;
	if (!InitializeBoundaryProlongationData(hierarchy.gridAtlases[0], boundaryProlongation)) {
		printf("ERROR : Failed boundary prolongation! \n");
		return 0;
	}

	std::vector<Point3D<double>> inputSignal;
	std::vector<double> texelToCellCoeffs;

	std::vector<std::vector<SquareMatrix<double, 2>>> parameterMetric;
	if (!InitializeMetric(mesh, EMBEDDING_METRIC, atlasCharts, parameterMetric)) {
		printf("ERROR: Unable to initialize metric \n");
		return 0;
	}

	if (!InitializeIntraChartEdgeIndexing(hierarchy.gridAtlases[0].gridCharts, edgeIndex)) {
		printf("ERROR: Unable to initialize intra chart edge indices \n");
		return 0;
	}

	if (Verbose.set) t_begin = clock();
	if (!_InitializeSystem(parameterMetric, boundaryProlongation)) return 0;
	if (Verbose.set) printf("\tInitialized system: %.2f(s)\n", double(clock() - t_begin) / CLOCKS_PER_SEC);
	
	if (!InitializeDivergenceRasteLines(edgeIndex, hierarchy.gridAtlases[0].rasterLines, divergenceRasterLines)) {
		printf("ERROR: Unable to initialize divergence raster lines! \n");
		return 0;
	}
	
	edgePairs.resize(edgeIndex.size());
	for (auto edgeIter = edgeIndex.begin(); edgeIter != edgeIndex.end(); edgeIter++) {
		unsigned long long edgeKey = (*edgeIter).first;
		int edgeId = (*edgeIter).second;

		unsigned long edgeCorners[2];
		GetMeshEdgeIndices(edgeKey, edgeCorners[0], edgeCorners[1]);
		edgePairs[edgeId] = std::pair<int, int>(edgeCorners[0], edgeCorners[1]);
	}

	texelMass.resize(textureNodes.size());
	texelDivergence.resize(textureNodes.size());

	if (UseDirectSolver.set)
	{
		FullMatrixConstruction(hierarchy.gridAtlases[0], deepMassCoefficients, boundaryBoundaryMassMatrix, boundaryDeepMassMatrix, mass);
		FullMatrixConstruction(hierarchy.gridAtlases[0], deepStiffnessCoefficients, boundaryBoundaryStiffnessMatrix, boundaryDeepStiffnessMatrix, stiffness);
		stitchingMatrix = mass*interpolationWeight + stiffness;
	}

	multigridIndices.resize(levels);
	for (int i = 0; i<levels; i++)
	{
		const GridAtlas & gridAtlas = hierarchy.gridAtlases[i];
		multigridIndices[i].threadTasks = gridAtlas.threadTasks;
		multigridIndices[i].boundaryGlobalIndex = gridAtlas.boundaryGlobalIndex;
		multigridIndices[i].segmentedLines = gridAtlas.segmentedLines;
		multigridIndices[i].rasterLines = gridAtlas.rasterLines;
		multigridIndices[i].restrictionLines = gridAtlas.restrictionLines;
		multigridIndices[i].prolongationLines = gridAtlas.prolongationLines;
		if (i<levels - 1) multigridIndices[i].boundaryRestriction = hierarchy.boundaryRestriction[i];
	}

	if (Verbose.set) t_begin = clock();
	if (!UpdateLinearSystem(interpolationWeight, 1.0, hierarchy, multigridStitchingCoefficients,
		deepMassCoefficients, deepStiffnessCoefficients,
		boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix,
		boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix,
		coarseSolver, boundarySolver, directSolver,
		stitchingMatrix, DetailVerbose.set, true, UseDirectSolver.set)) {
		printf("ERROR : Failed system update! \n");
		return 0;
	}
	if (Verbose.set) printf("\tInitialized multigrid coefficients: %.2f(s)\n", double(clock() - t_begin) / CLOCKS_PER_SEC);

	multigridStitchingVariables.resize(levels);
	for (int i = 0; i<levels; i++)
	{
		MultigridLevelVariables< Point3D< Real > >& variables = multigridStitchingVariables[i];
		variables.x.resize(hierarchy.gridAtlases[i].numTexels);
		variables.rhs.resize(hierarchy.gridAtlases[i].numTexels);
		variables.residual.resize(hierarchy.gridAtlases[i].numTexels);
		variables.boundary_rhs.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
		variables.boundary_value.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
		variables.variable_boundary_value.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
	}

	return 1;
}

template<class Real>
void Stitching<Real>::SolveSytem() {

	texelMass.resize(textureNodes.size());
	MultiplyBySystemMatrix_NoReciprocals(deepMassCoefficients, boundaryDeepMassMatrix, boundaryBoundaryMassMatrix, hierarchy.gridAtlases[0].boundaryGlobalIndex, hierarchy.gridAtlases[0].rasterLines, texelValues, texelMass);

	ComputeDivergence(edgeValues, texelDivergence, deepDivergenceCoefficients, boundaryDivergenceMatrix, divergenceRasterLines);

#pragma omp parallel for
	for (int i = 0; i <texelValues.size(); i++) {
		multigridStitchingVariables[0].x[i] = texelValues[i];
		multigridStitchingVariables[0].rhs[i] = texelMass[i] * interpolationWeight + texelDivergence[i];
	}

	filteredTexture.resize(textureWidth, textureHeight);
	for (int i = 0; i < filteredTexture.size(); i++) filteredTexture[i] = Point3D<float>(0.5f, 0.5f, 0.5f);
	
	if (UseDirectSolver.set) ComputeExactSolution(Verbose.set);
	else for (int i = 0; i < OutputVCycles.value; i++)VCycle(multigridStitchingVariables, multigridStitchingCoefficients, multigridIndices, boundarySolver, coarseSolver, false, false);
	
	UpdateFilteredTexture(multigridStitchingVariables[0].x);

}

template<class Real>
void Stitching<Real>::LoadImages() {

	if ( inputMode==MULTIPLE_INPUT_MODE )
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
			Image< Point3D< float > > textureConfidence;
			textureConfidence.read( confidenceName );
			inputConfidence[i].resize( textureWidth , textureHeight );
			for( int p=0 ; p<textureConfidence.size() ; p++ ) inputConfidence[i][p] = Point3D< float >::Dot( textureConfidence[p] , Point3D< float >( 1.f/3 , 1.f/3 , 1.f/3 ) );
		}
	}
	else
	{
		inputComposition.read( In.values[1] );
		textureWidth = inputComposition.width();
		textureHeight = inputComposition.height();

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

template<class Real>
void Stitching<Real>::ParseImages() {

	int width = textureWidth;
	int height = textureHeight;

	int numNodes = (int)textureNodes.size();
	int numEdges = (int)edgeIndex.size();
	unobservedTexel.resize(numNodes);
	texelValues.resize(numNodes);
	edgeValues.resize(numEdges);

	if (inputMode == MULTIPLE_INPUT_MODE){
		std::vector<Real> texelWeight;
		std::vector<Real> edgeWeight;
		texelWeight.resize(numNodes, 0);
		edgeWeight.resize(numEdges, 0);

		partialTexelValues.resize(numTextures);
		for (int i = 0; i < numTextures; i++) partialTexelValues[i].resize(numNodes);
		partialEdgeValues.resize(numTextures);
		for (int i = 0; i < numTextures; i++) partialEdgeValues[i].resize(edgeIndex.size());

		for (int textureIter = 0; textureIter < numTextures; textureIter++) {
			Image<Point3D<float>> textureValues = inputTextures[textureIter];
			Image<float> textureConfidence = inputConfidence[textureIter];

			for (int i = 0; i < textureNodes.size(); i++) {
				Real weight = textureConfidence(textureNodes[i].ci, textureNodes[i].cj);
				Point3D<Real> value = textureValues(textureNodes[i].ci, textureNodes[i].cj);
				partialTexelValues[textureIter][i] = value;
				texelValues[i] += value * weight;
				texelWeight[i] += weight;
			}

			for (int e = 0; e < edgePairs.size(); e++) {
				int edgeCorners[2] = { edgePairs[e].first,edgePairs[e].second };

				Real weight[2];
				Point3D<Real> value[2];
				for (int k = 0;k < 2; k++) {
					int i = edgeCorners[k];
					weight[k] = textureConfidence(textureNodes[i].ci, textureNodes[i].cj);
					value[k] = textureValues(textureNodes[i].ci, textureNodes[i].cj);
				}
				Real eWeight = weight[0] * weight[1];
				Point3D<Real> eValue = value[1] - value[0];
				partialEdgeValues[textureIter][e] = eWeight > 0 ? eValue : Point3D<Real>();
				edgeValues[e] += eValue*eWeight;
				edgeWeight[e] += eWeight;
			}

		}
		for (int i = 0; i < numNodes; i++) {
			unobservedTexel[i] = (texelWeight[i] == 0);
			if (texelWeight[i] > 0) texelValues[i] /= texelWeight[i];
		}

		for (int i = 0; i < numEdges; i++) {
			if (edgeWeight[i] > 0) edgeValues[i] /= edgeWeight[i];
		}
	}
	else {
		for (int i = 0; i < textureNodes.size(); i++) {
			unobservedTexel[i] = (inputMask(textureNodes[i].ci, textureNodes[i].cj) == -1);
			texelValues[i] = inputComposition(textureNodes[i].ci, textureNodes[i].cj);
		}
		for (int e = 0; e < edgePairs.size(); e++) {
			int edgeCorners[2] = { edgePairs[e].first,edgePairs[e].second };
			int ci[2] = { textureNodes[edgeCorners[0]].ci ,textureNodes[edgeCorners[1]].ci };
			int cj[2] = { textureNodes[edgeCorners[0]].cj ,textureNodes[edgeCorners[1]].cj };
			if (inputMask(ci[0], cj[0]) != -1 && inputMask(ci[0], cj[0]) == inputMask(ci[1], cj[1])) {
				edgeValues[e] = inputComposition(ci[1], cj[1]) - inputComposition(ci[0], cj[0]);
			}
			else {
				edgeValues[e] = Point3D<Real>(0, 0, 0);
			}
		}
	}

	// Set unobserved texel values to be the average of observed ones
	Point3D<Real> avgObservedTexel;
	int observedTexelCount = 0;
	for (int i = 0; i < unobservedTexel.size(); i++) {
		if (!unobservedTexel[i]) {
			avgObservedTexel += texelValues[i];
			observedTexelCount++;
		}
	}
	avgObservedTexel /= Real(observedTexelCount);
	for (int i = 0; i < unobservedTexel.size(); i++)if (unobservedTexel[i])texelValues[i] = avgObservedTexel;

#pragma omp parallel for
	for (int i = 0; i <textureNodes.size(); i++) multigridStitchingVariables[0].x[i] = texelValues[i];

}

template<class Real>
void Stitching<Real>::InitializeVisualization(void)
{
	visualization.textureWidth = textureWidth;
	visualization.textureHeight = textureHeight;

	visualization.colorTextureBuffer = new unsigned char[textureHeight*textureWidth * 3];
	memset(visualization.colorTextureBuffer, 204, textureHeight * textureWidth * 3 * sizeof(unsigned char));



	int tCount = (int)mesh.triangles.size();

	visualization.triangles.resize(tCount);
	visualization.vertices.resize(3 * tCount);
	visualization.colors.resize(3 * tCount, Point3D<double>(0.75, 0.75, 0.75));
	visualization.textureCoordinates.resize(3 * tCount);
	visualization.normals.resize(3 * tCount);


	for (int i = 0; i < tCount; i++) for (int k = 0; k < 3; k++) visualization.triangles[i][k] = 3 * i + k;

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

	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 's', "export texture", "Output Texture", ExportTextureCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'y', "interpolation weight", "Interpolation Weight", InterpolationWeightCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, ' ', "toggle update", ToggleUpdateCallBack));
	visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, '+', "increment update", IncrementUpdateCallBack));
	
	if(inputMode == MULTIPLE_INPUT_MODE)
	{
		visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 't', "toggle reference", ToggleForwardReferenceTextureCallBack));
		visualization.callBacks.push_back(Visualization::KeyboardCallBack(&visualization, 'T', "toggle reference", ToggleBackwardReferenceTextureCallBack));
	}
	
	visualization.info.push_back(stepsString);
	
	visualization.info.push_back(interpolationStr);
	if (inputMode == MULTIPLE_INPUT_MODE) visualization.info.push_back(referenceTextureStr);

	visualization.UpdateVertexBuffer();
	visualization.UpdateFaceBuffer();
	visualization.UpdateTextureBuffer(filteredTexture);

	if (inputMode == MULTIPLE_INPUT_MODE) visualization.UpdateReferenceTextureBuffers(inputTextures);
	if (inputMode == SINGLE_INPUT_MODE) visualization.UpdateCompositeTextureBuffer(inputComposition);
}

template< class Real >
int Stitching< Real >::Init(void)
{
	sprintf(stepsString, "Steps: 0");
	levels = std::max<int>(Levels.value, 1);
	interpolationWeight = InterpolationWeight.value;
	sprintf(referenceTextureStr, "Reference Texture: %02d of %02d \n", textureIndex,numTextures);
	sprintf(interpolationStr, "Interpolation: %.2e\n", interpolationWeight);

	if( !ReadTexturedMesh( mesh, In.values[0] , NULL , DetailVerbose.set ) )
	{
		printf("Unable to read mesh data\n");
		return 0;
	}

	// Define centroid and scale for visualization
	Point3D<double> centroid(0.f, 0.f, 0.f);
	for (int i = 0; i < mesh.vertices.size(); i++) centroid += mesh.vertices[i];
	centroid /= double(mesh.vertices.size());
	double radius = 0.f;
	for (int i = 0; i < mesh.vertices.size(); i++) radius = std::max<double>(radius, Point3D<double>::Length(mesh.vertices[i] - centroid));
	for (int i = 0; i < mesh.vertices.size(); i++) mesh.vertices[i] = (mesh.vertices[i] - centroid) / radius;

	if (1) for (int i = 0; i < mesh.textureCoordinates.size(); i++)mesh.textureCoordinates[i][1] = 1.0 - mesh.textureCoordinates[i][1];

	if (RandomJitter.set)
	{
		srand(time(NULL));
		std::vector<Point2D < double >>randomOffset(mesh.vertices.size());
		double jitterScale = 1e-3 / double(std::max<int>(textureWidth, textureHeight));
		for (int i = 0; i < randomOffset.size(); i++) randomOffset[i] = Point2D < double >(1.0 - 2.0 * double(rand()) / double(RAND_MAX), 1.0 - 2.0 *  double(rand()) / double(RAND_MAX))*jitterScale;
		for (int i = 0; i < mesh.triangles.size(); i++) for (int k = 0; k < 3; k++)mesh.textureCoordinates[3 * i + k] += randomOffset[mesh.triangles[i][k]];
	}



	ComputePadding(padding, textureWidth, textureHeight, mesh.textureCoordinates, DetailVerbose.set);
	if (padding.nonTrivial)
	{
		PadTextureCoordinates(padding, textureWidth, textureHeight, mesh.textureCoordinates);
		if (inputMode == MULTIPLE_INPUT_MODE) {
			for (int i = 0; i < numTextures; i++) {
				PadImage(padding, inputTextures[i]);
				PadImage(padding, inputConfidence[i]);
			}
		}
		else {
			PadImage(padding, inputComposition);
			PadImage(padding, inputMask);
		}
		
		textureWidth += (padding.left + padding.right);
		textureHeight += (padding.bottom + padding.top);		
	}

	clock_t t = clock();
	if (!InitializeSystem(textureWidth, textureHeight))
	{
		printf("Unable to initialize system\n");
		return 0;
	}
	if (Verbose.set)
	{
		printf("Resolution: %d / %d x %d\n", (int)textureNodes.size(), textureWidth, textureHeight);
		printf("Initialized system %.2f(s)\n", double(clock() - t) / CLOCKS_PER_SEC);
		printf("Peak Memory (MB): %d\n", Miscellany::MemoryInfo::PeakMemoryUsageMB());
	}

	// Assign position to exterior nodes using barycentric-exponential map
	{
		FEM::RiemannianMesh< double > rMesh(GetPointer(mesh.triangles), mesh.triangles.size());
		rMesh.setMetricFromEmbedding(GetPointer(mesh.vertices));
		rMesh.makeUnitArea();
		Pointer(FEM::CoordinateXForm< double >) xForms = rMesh.getCoordinateXForms();

		for (int i = 0; i<textureNodes.size(); i++) if (textureNodes[i].tId != -1 && !textureNodes[i].isInterior)
		{
			FEM::HermiteSamplePoint< double > _p;
			_p.tIdx = textureNodes[i].tId;
			_p.p = Point2D< double >(1. / 3, 1. / 3);
			_p.v = textureNodes[i].barycentricCoords - _p.p;

			rMesh.exp(xForms, _p);

			textureNodes[i].tId = _p.tIdx;
			textureNodes[i].barycentricCoords = _p.p;
		}
	}

	textureNodePositions.resize(textureNodes.size());
	for (int i = 0; i < textureNodePositions.size(); i++) {
		Point2D<double> barincetricCoords = textureNodes[i].barycentricCoords;
		int tId = textureNodes[i].tId;
		Point3D<float> surfacePosition = mesh.vertices[mesh.triangles[tId][0]] * (1.0 - barincetricCoords[0] - barincetricCoords[1]) +
			mesh.vertices[mesh.triangles[tId][1]] * barincetricCoords[0] +
			mesh.vertices[mesh.triangles[tId][2]] * barincetricCoords[1];
		textureNodePositions[i] = surfacePosition;
	}

	textureEdgePositions.resize(edgePairs.size());
	for (int i = 0; i < edgePairs.size(); i++) {
		textureEdgePositions[i] = (textureNodePositions[edgePairs[i].first] + textureNodePositions[edgePairs[i].second]) / 2.0;
	}

	if (true)
	{
		int multiChartTexelCount = 0;
		Image< int > texelId;
		texelId.resize(textureWidth, textureHeight);
		for (int i = 0; i<texelId.size(); i++) texelId[i] = -1;
		for (int i = 0; i<textureNodes.size(); i++)
		{
			int ci = textureNodes[i].ci, cj = textureNodes[i].cj;
			if (texelId(ci, cj) != -1)
			{
				if (false) fprintf(stderr, "[WARNING] Texel (%d %d) belong to multiple charts!\n", ci, cj);
				multiChartTexelCount++;
			}
			texelId(ci, cj) = i;
		}
		if (multiChartTexelCount) fprintf(stderr, "[WARNING] %d texels belong to multiple charts!\n", multiChartTexelCount);
	}

	return 1;
}

template<class Real>
int _main(int argc, char* argv[])
{
	Stitching< Real >::inputMode = MultiInput.set ? MULTIPLE_INPUT_MODE : SINGLE_INPUT_MODE;

	Stitching< Real >::LoadImages();
	if( !Stitching< Real >::Init() ) return 0;
	Stitching< Real >::ParseImages();
	Stitching< Real >::SolveSytem();

	if( !Output.set )
	{
		glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
		Stitching<Real>::visualization.visualizationMode = Stitching< Real >::inputMode;
		Stitching<Real>::visualization.displayMode = TWO_REGION_DISPLAY;
		Stitching< Real >::visualization.screenWidth = 1600, Stitching< Real >::visualization.screenHeight = 800;

		glutInitWindowSize(Stitching<Real>::visualization.screenWidth, Stitching<Real>::visualization.screenHeight);
		glutInit(&argc, argv);
		char windowName[1024];
		sprintf(windowName, "Stitching");
		glutCreateWindow(windowName);
		if (glewInit() != GLEW_OK) fprintf(stderr, "[ERROR] glewInit failed\n"), exit(0);
		glutDisplayFunc(Stitching<Real>::Display);
		glutReshapeFunc(Stitching<Real>::Reshape);
		glutMouseFunc(Stitching<Real>::MouseFunc);
		glutMotionFunc(Stitching<Real>::MotionFunc);
		glutKeyboardFunc(Stitching<Real>::KeyboardFunc);
		glutIdleFunc(Stitching< Real >::Idle);
		if (CameraConfig.set) Stitching<Real>::visualization.ReadSceneConfigurationCallBack(&Stitching<Real>::visualization, CameraConfig.value);
		Stitching<Real>::InitializeVisualization();
		glutMainLoop();
	}
	else
	{
		Stitching< Real >::ExportTextureCallBack(&Stitching< Real >::visualization, Output.value);
	}

	return 0;
}

int main(int argc, char* argv[])
{
	cmdLineParse( argc-1 , argv+1 , params );
	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}
	omp_set_num_threads( Threads.value );
	if( !NoHelp.set && !Output.set )
	{
		printf( "+----------------------------------------------------------------------------+\n" );
		printf( "| Interface Controls:                                                        |\n" );
		printf( "|    [Left Mouse]:                rotate                                     |\n" );
		printf( "|    [Right Mouse]:               zoom                                       |\n" );
		printf( "|    [Left/Right Mouse] + [CTRL]: pan                                        |\n" );
		printf( "|    [Left Mouse] + [SHIFT]:      mark region to in-paint (multi-mode only)  |\n" );
		printf( "|    't':                         toggle textures forward (multi-mode only)  |\n" );
		printf( "|    'T':                         toggle textures backward (multi-mode only) |\n" );
		printf( "|    'y':                         prescribe interpolation weight             |\n" );
		printf( "+----------------------------------------------------------------------------+\n" );
	}
	if( Double.set ) _main< double >( argc , argv );
	else             _main< float  >( argc , argv );
	return 0;
}
