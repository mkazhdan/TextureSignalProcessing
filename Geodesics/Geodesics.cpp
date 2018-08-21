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

#include <Src/PrecisionType.inl>
#include <Misha/CmdLineParser.h> 
#include <Misha/Miscellany.h>
#include <Src/SimpleMesh.h>
#include <Src/Basis.h>
#include <Misha/FEM.h>
#include <Src/Solver.h>
#include <Src/Hierarchy.h>
#include <Src/VectorFieldIntergration.inl>
#include <Src/MassAndStiffness.h>
#include <Src/Padding.h>
#include <Src/TexturedMeshVisualization.h>

cmdLineParameter< char* > Input( "in" );
cmdLineParameter< int   > Width( "width" , 1024 );
cmdLineParameter< int   > Height( "height" , 1024 );
cmdLineParameter< float > DiffusionInterpolationWeight( "interpolation" , 1e3 );
cmdLineParameter< int   > Levels( "levels" , 4 );
cmdLineParameter< char* > CameraConfig( "camera" );
cmdLineParameter< int   > Threads( "threads" , omp_get_num_procs() );
cmdLineParameter< int   > DisplayMode( "display" , TWO_REGION_DISPLAY );
cmdLineParameter< int   > MatrixQuadrature( "mQuadrature" , 6 );
cmdLineParameter< int   > VectorFieldQuadrature( "vfQuadrature" , 6 );

cmdLineParameter< int   > MultigridBlockHeight( "mBlockH" , 16 );
cmdLineParameter< int   > MultigridBlockWidth( "mBlockW" , 128 );
cmdLineParameter< int   > MultigridPaddedHeight( "mPadH" , 0 );
cmdLineParameter< int   > MultigridPaddedWidth( "mPadW" , 2 );
cmdLineParameter< int   > MultigridUpdateVcycles( "mVCycles" , 6 );

cmdLineReadable RandomJitter( "jitter" );
cmdLineReadable Verbose( "verbose" );
cmdLineReadable DetailVerbose( "detail" );
cmdLineReadable UseDirectSolver( "useDirectSolver" );
cmdLineReadable Double( "double" );

cmdLineReadable* params[] =
{
	&Input , &Width , &Height , &DiffusionInterpolationWeight , &CameraConfig , &Levels , &UseDirectSolver , &Threads , &DisplayMode , &MultigridBlockHeight , &MultigridBlockWidth , &MultigridPaddedHeight , &MultigridPaddedWidth , &MultigridUpdateVcycles ,
	&Verbose , &DetailVerbose ,
	&RandomJitter ,
	&Double ,
	&MatrixQuadrature , &VectorFieldQuadrature ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n", ex );
	printf( "\t --%s <input mesh>\n" , Input.name );
	printf( "\t[--%s <texture width>=%d]\n" , Width.name , Width.value );
	printf( "\t[--%s <texture height>=%d]\n" , Height.name , Height.value );
	printf( "\t[--%s <diffusion interpolation weight>=%f]\n" , DiffusionInterpolationWeight.name , DiffusionInterpolationWeight.value );
	printf( "\t[--%s <system matrix quadrature points per triangle>=%d]\n" , MatrixQuadrature.name , MatrixQuadrature.value );
	printf( "\t[--%s <normalized vector field quadrature points per triangle>=%d]\n" , VectorFieldQuadrature.name , VectorFieldQuadrature.value );
	printf( "\t[--%s]\n" , UseDirectSolver.name );
	printf( "\t[--%s]\n" , RandomJitter.name );
	printf( "\t[--%s]\n" , Verbose.name );

	printf( "\t[--%s <camera configuration file>]\n" , CameraConfig.name );
	printf( "\t[--%s <hierarchy levels>=%d]\n" , Levels.name , Levels.value );
	printf( "\t[--%s <threads>=%d]\n" , Threads.name , Threads.value );
	printf( "\t[--%s]\n" , DetailVerbose.name );
	printf( "\t[--%s <display mode>=%d]\n" , DisplayMode.name , DisplayMode.value );
	printf( "\t\t%d] One Region \n" , ONE_REGION_DISPLAY );
	printf( "\t\t%d] Two Region \n" , TWO_REGION_DISPLAY );
	printf( "\t[--%s <multigrid block width>=%d]\n" , MultigridBlockWidth.name , MultigridBlockWidth.value );
	printf( "\t[--%s <multigrid block height>=%d]\n" , MultigridBlockHeight.name , MultigridBlockHeight.value );
	printf( "\t[--%s <multigrid padded width>=%d]\n" , MultigridPaddedWidth.name , MultigridPaddedWidth.value );
	printf( "\t[--%s <multigrid padded height>=%d]\n" , MultigridPaddedHeight.name , MultigridPaddedHeight.value );
	printf( "\t[--%s <multigrid update VCycles>=%d]\n" , MultigridUpdateVcycles.name , MultigridUpdateVcycles.value );
}

template< class Real >
class Geodesics
{
public:
	static TexturedMesh mesh;
	static int textureWidth;
	static int textureHeight;
	static double diffusionInterpolationWeight;
	static double geodesicInterpolationWeight;
	static int levels;
	
	static Padding padding;
	
	static std::vector<Point3D<float>> textureNodePositions;

	static HierarchicalSystem hierarchy;

	static std::vector<CellIndex> cellIndices;

	static std::vector<TextureNodeInfo> textureNodes;
	static Image<int> nodeIndex;

	static SparseMatrix<double, int> mass;
	static SparseMatrix<double, int> stiffness;

	static SparseMatrix< double, int > smoothImpulseMatrix;
	static SparseMatrix< double, int > geodesicDistanceMatrix;

	static int impulseTexel;

	static std::vector<AtlasChart> atlasCharts;
	static std::vector<std::vector<SquareMatrix<double, 2>>> parameterMetric;

	static double smoothImpulseRange;
	static double geodesicDistanceRange;

	static std::vector<Real> exactSmoothImpulseSolution;
	static std::vector<Real> exactGeodesicDistanceSolution;

	static std::vector<Real> smoothedImpulseProlongation;
	static std::vector<Real> integratedVFProlongation;

	//Impulse Smoothing
	static std::vector<MultigridLevelCoefficients<Real>> multigridSmoothImpulseCoefficients;
	static std::vector<MultigridLevelVariables<Real>> multigridSmoothImpulseVariables;
	
	//Geodesic Distance
	static std::vector<MultigridLevelCoefficients<Real>> multigridGeodesicDistanceCoefficients;
	static std::vector<MultigridLevelVariables<Real>> multigridGeodesicDistanceVariables;

#if USE_CHOLMOD
	typedef  std::vector<CholmodCholeskySolver1<Real>> BoundarySolverType;
	typedef  CholmodCholeskySolver1<Real>  CoarseSolverType;
	typedef  CholmodCholeskySolver1<Real> DirectSolverType;
#elif USE_EIGEN_SIMPLICIAL
	typedef  std::vector<EigenCholeskySolver1<Real>> BoundarySolverType;
	typedef  EigenCholeskySolver1<Real>  CoarseSolverType;
	typedef  EigenCholeskySolver1<Real> DirectSolverType;
#else
	typedef  std::vector<EigenPardisoSolver1<Real>> BoundarySolverType;
	typedef  EigenPardisoSolver1<Real> CoarseSolverType;
	typedef  EigenPardisoSolver1<Real> DirectSolverType;
#endif

	static BoundarySolverType	boundarySmoothImpulseSolver;
	static BoundarySolverType	boundaryGeodesicDistanceSolver;

	static CoarseSolverType coarseSmoothImpulseSolver;
	static CoarseSolverType coarseGeodesicDistanceSolver;

	static DirectSolverType fineSmoothImpulseSolver;
	static DirectSolverType fineGeodesicDistanceSolver;

	static std::vector<MultigridLevelIndices<Real>> multigridIndices;

	static SparseMatrix<Real, int> coarseBoundaryFineBoundaryProlongation;
	static SparseMatrix<Real, int> fineBoundaryCoarseBoundaryRestriction;
	static std::vector<Real> coarseBoundaryValues;
	static std::vector<Real> coarseBoundaryRHS;
	static std::vector<Real> fineBoundaryValues;
	static std::vector<Real> fineBoundaryRHS;

	//Samples
	static std::vector< TriangleElementSampleInfo< Real > > triangleElementSamples;
	static std::vector< std::vector< SquareElementLineSampleInfo< Real > > > squareElementLineSamples;
	static std::vector<InteriorCellLine> interiorCellLines;
	static std::vector<std::pair<int, int>> interiorCellLineIndex;

	static int updateVCycles;
	static int cycleCount;

	static unsigned char * outputBuffer;

	//Visulization
	static TexturedMeshVisualization visualization;
	static int mouseX, mouseY;
	static bool mouseSelectionActive;

	static void UpdateOutputBuffer(const std::vector<Real> & solution);

	static void ExportTextureCallBack(Visualization* v, const char* prompt);
	static int Init();
	static void InitializeVisualization(const int width, const int height);
	static void ComputeExactSolution(bool verbose = false);
	static int UpdateSolution(bool verbose = false, bool detailVerbose = false);
	static int InitializeSystem(const int width, const int height);

	static void Display(void) { visualization.Display(); }
	static void MouseFunc(int button, int state, int x, int y);
	static void MotionFunc(int x, int y);
	static void Reshape(int w, int h) { visualization.Reshape(w, h); }
	static void KeyboardFunc(unsigned char key, int x, int y) { visualization.KeyboardFunc(key, x, y); }
	static void Idle();
};

template<class Real> TexturedMesh												Geodesics<Real>::mesh;
template<class Real> int														Geodesics<Real>::textureWidth;
template<class Real> int														Geodesics<Real>::textureHeight;
template<class Real> TexturedMeshVisualization									Geodesics<Real>::visualization;
template<class Real> int														Geodesics<Real>::mouseX = -1;
template<class Real> int														Geodesics<Real>::mouseY = -1;
template<class Real> bool														Geodesics<Real>::mouseSelectionActive = false;
template<class Real> Padding													Geodesics<Real>::padding;

template<class Real> std::vector<AtlasChart>									Geodesics<Real>::atlasCharts;
template<class Real> std::vector<std::vector<SquareMatrix<double, 2>>>			Geodesics<Real>::parameterMetric;
template<class Real> SparseMatrix<double, int>									Geodesics<Real>::mass;
template<class Real> SparseMatrix<double, int>									Geodesics<Real>::stiffness;

template<class Real> SparseMatrix<double, int>									Geodesics<Real>::smoothImpulseMatrix;
template<class Real> SparseMatrix<double, int>									Geodesics<Real>::geodesicDistanceMatrix;

template<class Real> double														Geodesics<Real>::diffusionInterpolationWeight;
template<class Real> double														Geodesics<Real>::geodesicInterpolationWeight = 1e-6;


template<class Real> std::vector<TextureNodeInfo>								Geodesics<Real>::textureNodes;
template<class Real> Image<int>													Geodesics<Real>::nodeIndex;
template<class Real> std::vector<CellIndex>										Geodesics<Real>::cellIndices;

template<class Real> int														Geodesics<Real>::levels;
template<class Real> HierarchicalSystem											Geodesics<Real>::hierarchy;

template<class Real> unsigned char *											Geodesics<Real>::outputBuffer;
template<class Real> std::vector<MultigridLevelIndices<Real>>					Geodesics<Real>::multigridIndices;

//Impulse Smoothing
template<class Real> std::vector<MultigridLevelCoefficients<Real>>						Geodesics<Real>::multigridSmoothImpulseCoefficients;
template<class Real> std::vector<MultigridLevelVariables<Real>>							Geodesics<Real>::multigridSmoothImpulseVariables;
template<class Real> typename Geodesics<Real>::CoarseSolverType						Geodesics<Real>::coarseSmoothImpulseSolver;

//Geodesic Distance
template<class Real> std::vector<MultigridLevelCoefficients<Real>>						Geodesics<Real>::multigridGeodesicDistanceCoefficients;
template<class Real> std::vector<MultigridLevelVariables<Real>>							Geodesics<Real>::multigridGeodesicDistanceVariables;
template<class Real> typename Geodesics<Real>::CoarseSolverType						Geodesics<Real>::coarseGeodesicDistanceSolver;

template<class Real> typename Geodesics<Real>::DirectSolverType						Geodesics<Real>::fineSmoothImpulseSolver;
template<class Real> typename Geodesics<Real>::DirectSolverType						Geodesics<Real>::fineGeodesicDistanceSolver;

template<class Real>  typename Geodesics<Real>::BoundarySolverType					Geodesics<Real>::boundarySmoothImpulseSolver;
template<class Real>  typename Geodesics<Real>::BoundarySolverType					Geodesics<Real>::boundaryGeodesicDistanceSolver;



template<class Real> std::vector<Real>												Geodesics<Real>::exactSmoothImpulseSolution;
template<class Real> std::vector<Real>												Geodesics<Real>::exactGeodesicDistanceSolution;

template<class Real> std::vector<Real>												Geodesics<Real>::smoothedImpulseProlongation;
template<class Real> std::vector<Real>												Geodesics<Real>::integratedVFProlongation;

//Samples
template<class Real> std::vector< TriangleElementSampleInfo< Real > >				Geodesics<Real>::triangleElementSamples;
template<class Real> std::vector< std::vector< SquareElementLineSampleInfo< Real > > >	Geodesics<Real>::squareElementLineSamples;
template<class Real> std::vector<InteriorCellLine>									Geodesics<Real>::interiorCellLines;
template<class Real> std::vector<std::pair<int, int>>								Geodesics<Real>::interiorCellLineIndex;

template<class Real> int															Geodesics<Real>::impulseTexel = -1;
template<class Real> std::vector<Point3D<float>>									Geodesics<Real>::textureNodePositions;

template<class Real> double															Geodesics<Real>::smoothImpulseRange;
template<class Real> double															Geodesics<Real>::geodesicDistanceRange;

template<class Real> SparseMatrix<Real, int>										Geodesics<Real>::coarseBoundaryFineBoundaryProlongation;
template<class Real> SparseMatrix<Real, int>										Geodesics<Real>::fineBoundaryCoarseBoundaryRestriction;

template<class Real> std::vector<Real>												Geodesics<Real>::coarseBoundaryValues;
template<class Real> std::vector<Real>												Geodesics<Real>::coarseBoundaryRHS;
template<class Real> std::vector<Real>												Geodesics<Real>::fineBoundaryValues;
template<class Real> std::vector<Real>												Geodesics<Real>::fineBoundaryRHS;

template<class Real> int															Geodesics<Real>::updateVCycles;
template<class Real> int															Geodesics<Real>::cycleCount;

template<class Real>
void Geodesics<Real>::ComputeExactSolution( bool verbose )
{
	clock_t begin;

	//(1) Smoothing impulse	
	if (verbose) begin = clock();
	solve( fineSmoothImpulseSolver , exactSmoothImpulseSolution , multigridSmoothImpulseVariables[0].rhs );
	if (verbose) printf("Smoothing impulse %.4f \n", double(clock() - begin) / CLOCKS_PER_SEC);

	//(1) Integrating vector field	
	const std::vector<int> & boundaryGlobalIndex = hierarchy.gridAtlases[0].boundaryGlobalIndex;

	if (verbose) begin = clock();
#pragma omp parallel for
	for (int i = 0; i < boundaryGlobalIndex.size(); i++) coarseBoundaryValues[i] = exactSmoothImpulseSolution[boundaryGlobalIndex[i]];
	coarseBoundaryFineBoundaryProlongation.Multiply(&coarseBoundaryValues[0], &fineBoundaryValues[0]);
	if (verbose) printf("Coarse to fine %.4f \n", double(clock() - begin) / CLOCKS_PER_SEC);


	std::vector< Real > fineGeodesicDistanceRHS( exactGeodesicDistanceSolution.size() );
	if( verbose ) begin = clock();
	auto VectorFunction = []( Point2D< Real > v , SquareMatrix< Real , 2 > tensor )
	{
		Point2D< Real > _v = tensor * v;
		double len2 = Point2D< Real >::Dot( v , _v );
		if( len2>0 ) return -v / (Real)sqrt( len2 );
		else         return -v;
	};
	if( !IntegrateVectorField< Real >( interiorCellLines , squareElementLineSamples , triangleElementSamples , exactSmoothImpulseSolution , VectorFunction , fineGeodesicDistanceRHS , fineBoundaryValues , fineBoundaryRHS , true ) )
	{
		printf( "[ERROR] Unable to integrate normalized vector field!\n" );
	}
	if( verbose ) printf( "Integrating normalized vector field %.4f\n" , double( clock()-begin ) / CLOCKS_PER_SEC);

	if( verbose ) begin = clock();
	fineBoundaryCoarseBoundaryRestriction.Multiply( &fineBoundaryRHS[0] , &coarseBoundaryRHS[0] );
#pragma omp parallel for
	for( int i=0 ; i<boundaryGlobalIndex.size() ; i++ ) fineGeodesicDistanceRHS[ boundaryGlobalIndex[i] ] += coarseBoundaryRHS[i];
	if( verbose ) printf( "Fine to coarse %.4f\n" , double( clock()-begin ) / CLOCKS_PER_SEC );


	//(3) Update geodesic distance solution	
	if( verbose ) begin = clock();
	solve( fineGeodesicDistanceSolver , exactGeodesicDistanceSolution , fineGeodesicDistanceRHS );
	if( verbose ) printf( "Computing geodesic distance %.4f\n" , double( clock()-begin ) / CLOCKS_PER_SEC );

	Real expectedMinDistance = exactGeodesicDistanceSolution[impulseTexel];

#pragma omp parallel for
	for( int i=0 ; i<exactGeodesicDistanceSolution.size() ; i++ ) exactGeodesicDistanceSolution[i] -= expectedMinDistance;
}


template<class Real>
void Geodesics<Real>::UpdateOutputBuffer(const std::vector<Real> & solution) {
	//Real expectedMaxDistance = solution[impulseTexel];

	double attenuationRadius = 0.2; //value within 0 and 0.5
	double attenuationStart = 0.8 - attenuationRadius;
#pragma omp parallel for
	for (int i = 0; i < textureNodes.size(); i++) {
		int ci = textureNodes[i].ci;
		int cj = textureNodes[i].cj;
		int offset = 3 * (textureWidth*cj + ci);
		Real value = solution[i];
		Point3D< float > red(255, 64, 64), blue(64, 64, 255), green(64, 255, 64); ;
		Point3D< float > color = value < 1.0 ? red * (1 - value) + blue * value : blue * (2.0 - value) + green * (value - 1.0);
		double scaledDistance = value * 60.0;
		int closestInteger = floor(scaledDistance);

		double residual = scaledDistance - double(closestInteger);
		if (closestInteger % 2 != 0) residual = 1.0 - residual;
		double attenuationWeight = (residual - attenuationStart) / (2.0 * attenuationRadius);
		attenuationWeight = std::min<double>(std::max<double>(0, attenuationWeight), 1.0);
		attenuationWeight = attenuationWeight*attenuationWeight*(3.0 - 2.0 * attenuationWeight);
		color = color*(1.0 - attenuationWeight) + Point3D< float >(0, 0, 0)*attenuationWeight;

		outputBuffer[offset + 0] = (unsigned char)(color[0]);
		outputBuffer[offset + 1] = (unsigned char)(color[1]);
		outputBuffer[offset + 2] = (unsigned char)(color[2]);

	}

	glBindTexture(GL_TEXTURE_2D, visualization.textureBuffer);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, textureWidth, textureHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)&outputBuffer[0]);
	glBindTexture(GL_TEXTURE_2D, 0);

	glutPostRedisplay();
}

template<class Real>
void Geodesics<Real>::Idle() {
	int selectedTexel = -1;
	if (mouseSelectionActive){
		Point3D<float> selectedPoint;
		bool validSelection = false;
		if (visualization.showMesh) {
			visualization.select(mouseX, mouseY, selectedPoint);
			float minDistance = FLT_MAX;
			for (int i = 0; i < textureNodePositions.size(); i++) {
				float squaredDistance = Point3D<float>::SquareNorm(textureNodePositions[i] - selectedPoint);
				if (squaredDistance < minDistance) {
					minDistance = squaredDistance;
					selectedTexel = i;
				}
			}
		}
		else {
			Point2D<float> ip = visualization.selectImagePos(mouseX, mouseY);
			//printf("Texture Coord %f %f \n", ip[0], ip[1]);
			int i = floor(ip[0] * float(nodeIndex.width()) - 0.5f);
			int j = floor((1.0 - ip[1]) * float(nodeIndex.height()) - 0.5f);
			//printf("Image pos %d %d \n", i, j);
			if (i >= 0 && i < nodeIndex.width() && j >= 0 && j < nodeIndex.height()) {
				selectedTexel = nodeIndex(i, j);
			}
		}

		if (impulseTexel != selectedTexel && selectedTexel != -1){
			impulseTexel = selectedTexel;
			memset(&multigridSmoothImpulseVariables[0].rhs[0], 0, multigridSmoothImpulseVariables[0].rhs.size() * sizeof(Real));
			multigridSmoothImpulseVariables[0].rhs[impulseTexel] = 1.0;
			memset(&multigridSmoothImpulseVariables[0].x[0], 0, multigridSmoothImpulseVariables[0].x.size() * sizeof(Real));
			memset(&multigridGeodesicDistanceVariables[0].x[0], 0, multigridGeodesicDistanceVariables[0].x.size() * sizeof(Real));

			cycleCount = 0;
		}

	}

	if (impulseTexel != -1) {
		if (cycleCount < updateVCycles) {
			if (!UpdateSolution()) {
				printf("Updated solution failed! \n");
			}
			UpdateOutputBuffer(multigridGeodesicDistanceVariables[0].x);
			cycleCount++;
		}
	}
}

template<class Real>
void Geodesics<Real>::MouseFunc(int button, int state, int x, int y) {

	visualization.newX = x; visualization.newY = y;
	visualization.rotating = visualization.scaling = visualization.panning = false;

	if (state == GLUT_DOWN && glutGetModifiers() & GLUT_ACTIVE_SHIFT) {

		int selectedTexel = -1;
		if (visualization.showMesh) {
			Point3D<float> selectedPoint;
			if (visualization.select(x, y, selectedPoint)) {
				mouseSelectionActive = true;
				mouseX = x;
				mouseY = y;
				float minDistance = FLT_MAX;
				for (int i = 0; i < textureNodePositions.size(); i++) {
					float squaredDistance = Point3D<float>::SquareNorm(textureNodePositions[i] - selectedPoint);
					if (squaredDistance < minDistance) {
						minDistance = squaredDistance;
						selectedTexel = i;
					}
				}
			}
		}
		else {
			Point2D<float> ip = visualization.selectImagePos(x, y);
			int i = floor(ip[0] * float(nodeIndex.width()) - 0.5f);
			int j = floor((1.0 - ip[1])*float(nodeIndex.height()) - 0.5f);
			if (i >= 0 && i < nodeIndex.width() && j >= 0 && j < nodeIndex.height()) {
				mouseSelectionActive = true;
				selectedTexel = nodeIndex(i, j);
			}
		}
		if (selectedTexel != -1 && selectedTexel != impulseTexel) {
			impulseTexel = selectedTexel;
			memset(&multigridSmoothImpulseVariables[0].rhs[0], 0, multigridSmoothImpulseVariables[0].rhs.size() * sizeof(Real));
			multigridSmoothImpulseVariables[0].rhs[impulseTexel] = 1.0;
			memset(&multigridSmoothImpulseVariables[0].x[0], 0, multigridSmoothImpulseVariables[0].x.size() * sizeof(Real));
			memset(&multigridGeodesicDistanceVariables[0].x[0], 0, multigridGeodesicDistanceVariables[0].x.size() * sizeof(Real));
		}
		if (impulseTexel != -1) {
			if (UseDirectSolver.set){
				clock_t begin;
				if( DetailVerbose.set ) begin = clock();
				ComputeExactSolution(Verbose.set);
				if( DetailVerbose.set ) printf("Exact solution %.4f \n", double(clock() - begin) / CLOCKS_PER_SEC);

				UpdateOutputBuffer(exactGeodesicDistanceSolution);
			}
			else {
				UpdateSolution();
				UpdateOutputBuffer(multigridGeodesicDistanceVariables[0].x);
			}
			glutPostRedisplay();
		}
	}
	else {
		mouseSelectionActive = false;
		if (button == GLUT_LEFT_BUTTON) {
			if (glutGetModifiers() & GLUT_ACTIVE_CTRL) visualization.panning = true;
			else                                        visualization.rotating = true;
		}
		else if (button == GLUT_RIGHT_BUTTON) visualization.scaling = true;
	}
}

template<class Real>
void Geodesics<Real>::MotionFunc(int x, int y) {
	if (mouseSelectionActive) {
		mouseX = x;
		mouseY = y;
		glutPostRedisplay();
	}
	else {
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

#ifdef GLM_FORCE_RADIANS
			if      ( visualization.rotating ) visualization.camera.rotateUp( visualization.center , -rUp ) , visualization.camera.rotateRight( visualization.center , rRight );
#else // !GLM_FORCE_RADIANS
			if      ( visualization.rotating ) visualization.camera.rotateUp( visualization.center , rUp ) , visualization.camera.rotateRight( visualization.center , rRight );
#endif // GLM_FORCE_RADIANS
			else if( visualization.scaling   ) visualization.camera.moveForward( pForward);
			else if( visualization.panning   ) visualization.camera.moveRight( pRight ), visualization.camera.moveUp( pUp );
		}
		glutPostRedisplay();
	}
}

template<class Real>
void Geodesics<Real>::ExportTextureCallBack(Visualization* v, const char* prompt) {

	Image<Point3D<float>> outputImage;
	outputImage.resize(textureWidth, textureHeight);
	for (int i = 0; i < outputImage.size(); i++) outputImage[i] = Point3D<float>(outputBuffer[3 * i], outputBuffer[3 * i + 1], outputBuffer[3 * i + 2]) / float(255.0);
	if( padding.nonTrivial ) UnpadImage( padding , outputImage );
	outputImage.write(prompt);
}

template<class Real>
int Geodesics<Real>::UpdateSolution(bool verbose, bool detailVerbose) {
	clock_t begin;

	//(1)Update smoothed input solution
	if( verbose ) begin = clock();
	VCycle( multigridSmoothImpulseVariables , multigridSmoothImpulseCoefficients , multigridIndices , boundarySmoothImpulseSolver , coarseSmoothImpulseSolver , detailVerbose , detailVerbose );
	if( verbose ) printf( "Smoothing impulse %.4f\n" , double(clock() - begin) / CLOCKS_PER_SEC );

	//(2) Integrate normalized vector field
	const std::vector<int> & boundaryGlobalIndex = hierarchy.gridAtlases[0].boundaryGlobalIndex;

	if (verbose) begin = clock();
#pragma omp parallel for
	for (int i = 0; i < boundaryGlobalIndex.size(); i++) coarseBoundaryValues[i] = multigridSmoothImpulseVariables[0].x[boundaryGlobalIndex[i]];
	coarseBoundaryFineBoundaryProlongation.Multiply(&coarseBoundaryValues[0], &fineBoundaryValues[0]);
	if (verbose) printf("Coarse to fine %.4f \n", double(clock() - begin) / CLOCKS_PER_SEC);

	if (verbose) begin = clock();
	auto VectorFunction = []( Point2D< Real > v , SquareMatrix< Real , 2 > tensor )
	{
		Point2D< Real > _v = tensor * v;
		double len2 = Point2D<Real>::Dot( v , _v );
		if( len2>0 ) return -v / (Real)sqrt( len2 );
		else         return -v;
	};
	if( !IntegrateVectorField< Real >( interiorCellLines , squareElementLineSamples , triangleElementSamples , multigridSmoothImpulseVariables[0].x , VectorFunction , multigridGeodesicDistanceVariables[0].rhs , fineBoundaryValues , fineBoundaryRHS ) )
	{
		printf( "[ERROR] Unable to integrate normalized vector field!\n" );
	}
	if (verbose) printf("Integrating normalized vector field %.4f \n", double(clock() - begin) / CLOCKS_PER_SEC);

	if (verbose) begin = clock();
	fineBoundaryCoarseBoundaryRestriction.Multiply(&fineBoundaryRHS[0], &coarseBoundaryRHS[0]);
#pragma omp parallel for
	for (int i = 0; i < boundaryGlobalIndex.size(); i++) multigridGeodesicDistanceVariables[0].rhs[boundaryGlobalIndex[i]] += coarseBoundaryRHS[i];
	if (verbose) printf("Fine to coarse %.4f \n", double(clock() - begin) / CLOCKS_PER_SEC);


	//(3) Update geodesic distance solution	
	if (verbose) begin = clock();
	VCycle(multigridGeodesicDistanceVariables, multigridGeodesicDistanceCoefficients, multigridIndices, boundaryGeodesicDistanceSolver, coarseGeodesicDistanceSolver, detailVerbose, detailVerbose);
	if (verbose) printf("Solving geodesic distance %.4f \n", double(clock() - begin) / CLOCKS_PER_SEC);

	Real expectedMinDistance = multigridGeodesicDistanceVariables[0].x[impulseTexel];

#pragma omp parallel for
	for (int i = 0; i < multigridGeodesicDistanceVariables[0].x.size(); i++) multigridGeodesicDistanceVariables[0].x[i] -= expectedMinDistance;

	return 1;
}


template<class Real>
int Geodesics<Real>::InitializeSystem( const int width , const int height )
{
	clock_t t_begin;
	
	t_begin = clock();
	MultigridBlockInfo multigridBlockInfo(MultigridBlockWidth.value, MultigridBlockHeight.value, MultigridPaddedWidth.value, MultigridPaddedHeight.value, 0);
	if( !InitializeHierarchy( mesh , width , height , levels , textureNodes , cellIndices , hierarchy , atlasCharts , multigridBlockInfo , true , DetailVerbose.set ) )
	{
		printf("ERROR : Failed intialization! \n");
		return 0;
	}
	if( Verbose.set ) printf( "\tInitialized hierarchy: %.2f(s)\n" , double(clock() - t_begin) / CLOCKS_PER_SEC);

	//Initialize node index
	nodeIndex.resize(width, height);
	for (int i = 0; i < nodeIndex.size(); i++)nodeIndex[i] = -1;
	for (int i = 0; i < textureNodes.size(); i++) {
		if (nodeIndex(textureNodes[i].ci, textureNodes[i].cj) != -1) {
			if(0)printf("WARNING: Multiple nodes mapped to pixel %d %d!\n", textureNodes[i].ci, textureNodes[i].cj);
		}
		nodeIndex(textureNodes[i].ci, textureNodes[i].cj) = i;
	}

	BoundaryProlongationData boundaryProlongation;
	if (!InitializeBoundaryProlongationData(hierarchy.gridAtlases[0], boundaryProlongation)) {
		printf("ERROR : Failed boundary prolongation! \n");
		return 0;
	}


	std::vector<double> deepMassCoefficients;
	std::vector<double> deepStiffnessCoefficients;
	SparseMatrix<double, int> boundaryBoundaryMassMatrix;
	SparseMatrix<double, int> boundaryBoundaryStiffnessMatrix;
	SparseMatrix<double, int> boundaryDeepMassMatrix;
	SparseMatrix<double, int> boundaryDeepStiffnessMatrix;


	std::vector<Point3D<double>> __inputSignal;
	std::vector<double> __texelToCellCoeffs;
	SparseMatrix<double, int> __boundaryCellBasedStiffnessRHSMatrix[3];

	if (!InitializeMetric(mesh,EMBEDDING_METRIC,atlasCharts, parameterMetric)) {
		printf("ERROR: Unable to initialize metric \n");
		return 0;
	}

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
		case 32:
			ret = InitializeMassAndStiffness<32>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix );
			break;
		default:
			fprintf( stderr , "[ERROR] Only 1-, 3-, 6-, and 32-point quadrature supported for triangles\n" );
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
		FullMatrixConstruction(hierarchy.gridAtlases[0], deepMassCoefficients, boundaryBoundaryMassMatrix, boundaryDeepMassMatrix, mass);
		FullMatrixConstruction(hierarchy.gridAtlases[0], deepStiffnessCoefficients, boundaryBoundaryStiffnessMatrix, boundaryDeepStiffnessMatrix, stiffness);
		smoothImpulseMatrix = mass * diffusionInterpolationWeight + stiffness;
		geodesicDistanceMatrix = mass * geodesicInterpolationWeight + stiffness;
		printf("Assembling matrices =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
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

	t_begin = clock();
	if (!UpdateLinearSystem( diffusionInterpolationWeight , 1.0 , hierarchy , multigridSmoothImpulseCoefficients ,
		deepMassCoefficients, deepStiffnessCoefficients,
		boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix,
		boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix,
		coarseSmoothImpulseSolver, boundarySmoothImpulseSolver, fineSmoothImpulseSolver,
		smoothImpulseMatrix, DetailVerbose.set, true, UseDirectSolver.set)){
		printf("ERROR : Failed system update! \n");
		return 0;
	}

	if (!UpdateLinearSystem( geodesicInterpolationWeight , 1. , hierarchy , multigridGeodesicDistanceCoefficients ,
		deepMassCoefficients, deepStiffnessCoefficients,
		boundaryBoundaryMassMatrix, boundaryBoundaryStiffnessMatrix,
		boundaryDeepMassMatrix, boundaryDeepStiffnessMatrix,
		coarseGeodesicDistanceSolver, boundaryGeodesicDistanceSolver, fineGeodesicDistanceSolver,
		geodesicDistanceMatrix, DetailVerbose.set, true, UseDirectSolver.set)) {
		printf("ERROR : Failed system update! \n");
		return 0;
	}
	if( Verbose.set ) printf( "\tInitialized multigrid coefficients: %.2f(s)\n" , double(clock() - t_begin) / CLOCKS_PER_SEC);

//////////////////////////////////// Initialize multigrid variables

	multigridSmoothImpulseVariables.resize(levels);
	for (int i = 0; i < levels; i++) {
		MultigridLevelVariables<Real> & variables = multigridSmoothImpulseVariables[i];
		variables.x.resize(hierarchy.gridAtlases[i].numTexels);
		variables.rhs.resize(hierarchy.gridAtlases[i].numTexels);
		variables.residual.resize(hierarchy.gridAtlases[i].numTexels);
		variables.boundary_rhs.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
		variables.boundary_value.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
		variables.variable_boundary_value.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
	}

	multigridGeodesicDistanceVariables.resize(levels);
	for (int i = 0; i < levels; i++) {
		MultigridLevelVariables<Real> & variables = multigridGeodesicDistanceVariables[i];
		variables.x.resize(hierarchy.gridAtlases[i].numTexels);
		variables.rhs.resize(hierarchy.gridAtlases[i].numTexels);
		variables.residual.resize(hierarchy.gridAtlases[i].numTexels);
		variables.boundary_rhs.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
		variables.boundary_value.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
		variables.variable_boundary_value.resize(hierarchy.gridAtlases[i].boundaryGlobalIndex.size());
	}

//////////////////////////////////// Initialize cell samples

	if (!InitializeGridAtlasInteriorCellLines(atlasCharts, hierarchy.gridAtlases[0].gridCharts, interiorCellLines, interiorCellLineIndex)) {
		printf("Unable to initialize interior cell lines! \n");
		return 0;
	}
	if (interiorCellLineIndex.size() != hierarchy.gridAtlases[0].numInteriorCells) {
		printf("ERROR: Inconsistent number of interior cells!. Expected %d . Result %d. \n", hierarchy.gridAtlases[0].numInteriorCells, (int)interiorCellLineIndex.size());
		return 0;
	}

	coarseBoundaryFineBoundaryProlongation = boundaryProlongation.coarseBoundaryFineBoundaryProlongation;
	fineBoundaryCoarseBoundaryRestriction = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction;
	std::vector<int> fineBoundaryIndex = boundaryProlongation.fineBoundaryIndex;
	int numFineBoundarNodes = boundaryProlongation.numFineBoundarNodes;

	squareElementLineSamples.resize( interiorCellLines.size() );

	t_begin = clock();
	{
		int ret = 0;
		switch( VectorFieldQuadrature.value )
		{
		case 1:
			ret = InitializeVectorFieldIntegration< 1>( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , squareElementLineSamples , triangleElementSamples );
			break;
		case 3:
			ret = InitializeVectorFieldIntegration< 3>( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , squareElementLineSamples , triangleElementSamples );
			break;
		case 6:
			ret = InitializeVectorFieldIntegration< 6>( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , squareElementLineSamples , triangleElementSamples );
			break;
		case 32:
			ret = InitializeVectorFieldIntegration<32>( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , squareElementLineSamples , triangleElementSamples );
			break;
		default:
			fprintf( stderr , "[ERROR] Only 1-, 3-, 6-, and 32-point quadrature supported for triangles\n" );
		}
		if( !ret )
		{
			fprintf( stderr , "[ERROR] Unable to initialize vector field integration samples!\n" );
			return 0;
		}
	}
	if( Verbose.set) printf( "\tInitialized vector field integration: %.2f(s)\n" , double(clock() - t_begin) / CLOCKS_PER_SEC);

	coarseBoundaryValues.resize(hierarchy.gridAtlases[0].numTexels - hierarchy.gridAtlases[0].numDeepTexels);
	coarseBoundaryRHS.resize(hierarchy.gridAtlases[0].numTexels - hierarchy.gridAtlases[0].numDeepTexels);
	fineBoundaryValues.resize(numFineBoundarNodes);
	fineBoundaryRHS.resize(numFineBoundarNodes);

	for( int i=0 ; i<squareElementLineSamples.size() ; i++ )  std::sort( squareElementLineSamples[i].begin() , squareElementLineSamples[i].end() , SquareElementLineSampleInfo< Real >::Compare );

	int numTexels = hierarchy.gridAtlases[0].numTexels;
	int numFineNodes = hierarchy.gridAtlases[0].numFineNodes;

	exactSmoothImpulseSolution.resize(numTexels);
	exactGeodesicDistanceSolution.resize(numTexels);
	smoothedImpulseProlongation.resize(numFineNodes);
	integratedVFProlongation.resize(numFineNodes);

	return 1;

}

template<class Real>
void Geodesics<Real>::InitializeVisualization( const int width , const int height )
{
	outputBuffer = new unsigned char[ height*width* 3];
	memset( outputBuffer , 204 , height * width * 3 * sizeof(unsigned char) );

	int tCount = (int)mesh.triangles.size();

	visualization.triangles.resize( tCount );
	visualization.vertices.resize( 3*tCount );
	visualization.colors.resize( 3*tCount , Point3D< double >( 0.75 , 0.75 , 0.75 ) );
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
	if( !InitializeBoundaryEdges( mesh , boundaryEdges ) ) fprintf( stderr , "[WARNING] Unable to initialize boundary edges!\n" );

	for( int e=0 ; e<boundaryEdges.size() ; e++ )
	{
		int tIndex = boundaryEdges[e] / 3;
		int kIndex = boundaryEdges[e] % 3;
		for( int c=0 ; c<2 ; c++ )
		{
			Point3D< double > v = mesh.vertices[ mesh.triangles[tIndex][ (kIndex+c)%3 ] ];
			visualization.boundaryEdgeVertices.push_back(v);
		}
	}

	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization,  's' , "export texture" , "Output Texture" , ExportTextureCallBack ) );

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
}

template< class Real >
int Geodesics<Real>::Init( void ){
	levels = Levels.value;
	diffusionInterpolationWeight = DiffusionInterpolationWeight.value;
	textureWidth = Width.value;
	textureHeight = Height.value;
	updateVCycles = MultigridUpdateVcycles.value;
	cycleCount = updateVCycles;

	if (!ReadTexturedMesh( mesh , Input.value , NULL , DetailVerbose.set ) )
	{
		printf("Unable to read mesh data\n");
		return 0;
	}
	if (1) for (int i = 0; i < mesh.textureCoordinates.size(); i++)mesh.textureCoordinates[i][1] = 1.0 - mesh.textureCoordinates[i][1];

	if( RandomJitter.set )
	{
		srand(time(NULL));
		std::vector<Point2D < double >>randomOffset(mesh.vertices.size());
		double jitterScale = 1e-3 / double(std::max<int>(textureWidth, textureHeight));
		for (int i = 0; i < randomOffset.size(); i++) randomOffset[i] = Point2D < double >(1.0 - 2.0 * double(rand()) / double(RAND_MAX), 1.0 - 2.0 *  double(rand()) / double(RAND_MAX))*jitterScale;
		for (int i = 0; i < mesh.triangles.size(); i++) for (int k = 0; k < 3; k++)mesh.textureCoordinates[3 * i + k] += randomOffset[mesh.triangles[i][k]];
	}

	ComputePadding( padding , textureWidth , textureHeight , mesh.textureCoordinates , DetailVerbose.set );
	if (padding.nonTrivial) {
		PaddTextureCoordinates(padding, textureWidth, textureHeight, mesh.textureCoordinates);
		textureWidth += (padding.left + padding.right);
		textureHeight += (padding.bottom + padding.top);
	}

	//Define centroid and scale for visualization
	Point3D< double > centroid(0.f, 0.f, 0.f);
	for( int i=0 ; i<mesh.vertices.size() ; i++ ) centroid += mesh.vertices[i];
	centroid /= (double)mesh.vertices.size();
	double radius = 0;
	for( int i=0 ; i<mesh.vertices.size() ; i++ ) radius = std::max< double >( radius , Point3D< double >::Length( mesh.vertices[i]-centroid) );
	for( int i=0 ; i<mesh.vertices.size() ; i++ ) mesh.vertices[i] = ( mesh.vertices[i]-centroid ) / radius;


	clock_t t = clock();
	if( !InitializeSystem( textureWidth , textureHeight ) ){ fprintf( stderr , "[ERROR] Unable to initialize system\n") ; return 0; }

	if( Verbose.set )
	{
		printf( "Resolution: %d / %d x %d\n" , (int)textureNodes.size() , textureWidth , textureHeight );
		printf( "Initialized system %.2f(s)\n" , double(clock()-t) / CLOCKS_PER_SEC );
		printf( "Peak Memory (MB): %d\n" , Miscellany::MemoryInfo::PeakMemoryUsageMB() );
	}

	//Assign position to exterior nodes using baricentric-exponential map
	FEM::RiemannianMesh< double > * Rmesh = new FEM::RiemannianMesh< double >(GetPointer(mesh.triangles), mesh.triangles.size());
	Rmesh->setMetricFromEmbedding(GetPointer(mesh.vertices));
	Rmesh->makeUnitArea();
	Pointer(FEM::CoordinateXForm< double >) xForms = Rmesh->getCoordinateXForms();

	for (int i = 0; i<textureNodes.size(); i++) {
		if (textureNodes[i].tId != -1 && !textureNodes[i].isInterior) {
			FEM::HermiteSamplePoint< double > _p;
			_p.tIdx = textureNodes[i].tId;
			_p.p = Point2D< double >((double)1. / 3, (double)1. / 3);
			_p.v = textureNodes[i].baricentricCoords - _p.p;

			Rmesh->exp(xForms, _p);

			textureNodes[i].tId = _p.tIdx;
			textureNodes[i].baricentricCoords = _p.p;
		}
	}
	
	textureNodePositions.resize( textureNodes.size() );
	for( int i=0 ; i<textureNodePositions.size() ; i++ )
	{
		Point2D< double > barincetricCoords = textureNodes[i].baricentricCoords;
		int tId = textureNodes[i].tId;
		Point3D<float> surfacePosition =
			mesh.vertices[ mesh.triangles[tId][0] ] * ( 1.0-barincetricCoords[0]-barincetricCoords[1] ) +
			mesh.vertices[ mesh.triangles[tId][1] ] * barincetricCoords[0] +
			mesh.vertices[ mesh.triangles[tId][2] ] * barincetricCoords[1];
		textureNodePositions[i] = surfacePosition;
	}
	return 1;
}

template<class Real>
int _main(int argc, char* argv[])
{
	if( !Geodesics<Real>::Init() ) return 0;

	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	Geodesics<Real>::visualization.displayMode = DisplayMode.value;
	if     ( DisplayMode.value==ONE_REGION_DISPLAY ) Geodesics<Real>::visualization.screenWidth =  800 , Geodesics<Real>::visualization.screenHeight = 800;
	else if( DisplayMode.value==TWO_REGION_DISPLAY ) Geodesics<Real>::visualization.screenWidth = 1440 , Geodesics<Real>::visualization.screenHeight = 720;

	Geodesics< Real >::visualization.UpdateMainFrameSize();
	glutInitWindowSize( Geodesics< Real >::visualization.screenWidth , Geodesics< Real >::visualization.screenHeight );

	glutInit(&argc, argv);
	char windowName[1024];
	sprintf(windowName, "Goedsics");
	glutCreateWindow(windowName);
	if (glewInit() != GLEW_OK) fprintf(stderr, "[ERROR] glewInit failed\n"), exit(0);
	glutDisplayFunc ( Geodesics< Real >::Display);
	glutReshapeFunc ( Geodesics< Real >::Reshape);
	glutMouseFunc   ( Geodesics< Real >::MouseFunc);
	glutMotionFunc  ( Geodesics< Real >::MotionFunc);
	glutKeyboardFunc( Geodesics< Real >::KeyboardFunc);
	if( !UseDirectSolver.set ) glutIdleFunc( Geodesics< Real >::Idle );
	if( CameraConfig.set ) Geodesics< Real >::visualization.ReadSceneConfigurationCallBack( &Geodesics<Real>::visualization , CameraConfig.value );
	Geodesics< Real >::InitializeVisualization( Geodesics< Real >::textureWidth , Geodesics<Real>::textureHeight );
	glutMainLoop();

	return 0;

}

int main( int argc , char* argv[] )
{
	cmdLineParse(argc - 1, argv + 1, params);
	if( !Input.set ){ ShowUsage( argv[0] ) ; return EXIT_FAILURE; }
	omp_set_num_threads( Threads.value );
	if( Double.set ) _main< double >( argc , argv );
	else             _main< float  >( argc , argv );
	return 0;
}
