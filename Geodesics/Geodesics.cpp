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
#include <Misha/FEM.h>
#include <Src/Hierarchy.h>
#include <Src/SimpleTriangleMesh.h>
#include <Src/Basis.h>
#include <Src/Solver.h>
#include <Src/QuadratureIntergration.inl>
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

cmdLineParameter< int   > MultigridBlockHeight ( "mBlockH" ,  16 );
cmdLineParameter< int   > MultigridBlockWidth  ( "mBlockW" , 128 );
cmdLineParameter< int   > MultigridPaddedHeight( "mPadH"   ,   0 );
cmdLineParameter< int   > MultigridPaddedWidth ( "mPadW"   ,   2 );

cmdLineParameter< int   > RandomJitter( "jitter" , 0 );
cmdLineReadable Verbose( "verbose" );
cmdLineReadable NoHelp( "noHelp" );
cmdLineReadable DetailVerbose( "detail" );
cmdLineReadable UseDirectSolver( "useDirectSolver" );
cmdLineReadable Double( "double" );
cmdLineReadable PreciseIntegration( "preciseIntegration" );

cmdLineReadable* params[] =
{
	&Input , &Width , &Height , &DiffusionInterpolationWeight , &CameraConfig , &Levels , &UseDirectSolver , &Threads , &DisplayMode , &MultigridBlockHeight , &MultigridBlockWidth , &MultigridPaddedHeight , &MultigridPaddedWidth ,
	&Verbose , &DetailVerbose ,
	&RandomJitter ,
	&Double ,
	&MatrixQuadrature , &VectorFieldQuadrature ,
	&PreciseIntegration ,
	&NoHelp ,
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
	printf( "\t[--%s]\n" , PreciseIntegration.name );
	printf( "\t[--%s]\n" , UseDirectSolver.name );
	printf( "\t[--%s <jittering seed>]\n" , RandomJitter.name );
	printf( "\t[--%s]\n" , Verbose.name );

	printf( "\t[--%s <camera configuration file>]\n" , CameraConfig.name );
	printf( "\t[--%s <hierarchy levels>=%d]\n" , Levels.name , Levels.value );
	printf( "\t[--%s <threads>=%d]\n" , Threads.name , Threads.value );
	printf( "\t[--%s]\n" , DetailVerbose.name );
	printf( "\t[--%s <display mode>=%d]\n" , DisplayMode.name , DisplayMode.value );
	printf( "\t\t%d] One Region \n" , ONE_REGION_DISPLAY );
	printf( "\t\t%d] Two Region \n" , TWO_REGION_DISPLAY );
	printf( "\t[--%s <multigrid block width>=%d]\n"   , MultigridBlockWidth.name   , MultigridBlockWidth.value   );
	printf( "\t[--%s <multigrid block height>=%d]\n"  , MultigridBlockHeight.name  , MultigridBlockHeight.value  );
	printf( "\t[--%s <multigrid padded width>=%d]\n"  , MultigridPaddedWidth.name  , MultigridPaddedWidth.value  );
	printf( "\t[--%s <multigrid padded height>=%d]\n" , MultigridPaddedHeight.name , MultigridPaddedHeight.value );
	printf( "\t[--%s]\n" , NoHelp.name );
}

template< typename PreReal , typename Real >
class Geodesics
{
public:
	static OrientedTexturedTriangleMesh< PreReal > mesh;
	static int textureWidth;
	static int textureHeight;
	static Real diffusionInterpolationWeight;
	static Real geodesicInterpolationWeight;
	static int levels;
	static int updateCount;
	
	static int steps;
	static char stepsString[];

	static Padding padding;
	
	static std::vector< Point3D< float > > textureNodePositions;

	static HierarchicalSystem< PreReal , Real > hierarchy;

	static std::vector< BilinearElementIndex > bilinearElementIndices;

	static std::vector< TextureNodeInfo< PreReal > > textureNodes;
	static Image<int> nodeIndex;

	static SparseMatrix< Real , int > mass;
	static SparseMatrix< Real , int > stiffness;

	static SparseMatrix< Real , int > smoothImpulseMatrix;
	static SparseMatrix< Real , int > geodesicDistanceMatrix;

	static int impulseTexel;

	static std::vector< AtlasChart< PreReal > > atlasCharts;
	static std::vector< std::vector< SquareMatrix< PreReal , 2 > > > parameterMetric;

	static Real smoothImpulseRange;
	static Real geodesicDistanceRange;

	//Impulse Smoothing
	static std::vector< SystemCoefficients< Real > > multigridSmoothImpulseCoefficients;
	static std::vector<MultigridLevelVariables<Real>> multigridSmoothImpulseVariables;
	
	//Geodesic Distance
	static std::vector< SystemCoefficients< Real > > multigridGeodesicDistanceCoefficients;
	static std::vector<MultigridLevelVariables<Real>> multigridGeodesicDistanceVariables;

#if defined( USE_CHOLMOD )
	typedef CholmodCholeskySolver< Real , 1 > DirectSolver;
#elif defined( USE_EIGEN_SIMPLICIAL )
	typedef EigenCholeskySolver< Real , 1 > DirectSolver;
#elif defined( USE_EIGEN_PARDISO )
	typedef EigenPardisoSolver< Real , 1 > DirectSolver;
#else
#error "[ERROR] No solver defined!"
#endif

	static VCycleSolvers< DirectSolver > smoothImpulseSolvers;
	static VCycleSolvers< DirectSolver > geodesicDistanceSolvers;
	static DirectSolver fineSmoothImpulseSolver;
	static DirectSolver fineGeodesicDistanceSolver;

	static std::vector<MultigridLevelIndices<Real>> multigridIndices;

	static SparseMatrix<Real, int> coarseBoundaryFineBoundaryProlongation;
	static SparseMatrix<Real, int> fineBoundaryCoarseBoundaryRestriction;
	static std::vector<Real> coarseBoundaryValues;
	static std::vector<Real> coarseBoundaryRHS;
	static std::vector<Real> fineBoundaryValues;
	static std::vector<Real> fineBoundaryRHS;

	//Samples
	static GradientElementSamples< Real > gradientSamples;
	static std::vector<InteriorCellLine> interiorCellLines;
	static std::vector<std::pair<int, int>> interiorCellLineIndex;

	static unsigned char * outputBuffer;

	//Visulization
	static TexturedMeshVisualization visualization;
	static int mouseX, mouseY;
	static bool mouseSelectionActive;

	static void UpdateOutputBuffer( const std::vector< Real >& solution );

	static void ToggleUpdateCallBack( Visualization* v , const char* prompt );
	static void IncrementUpdateCallBack( Visualization* v , const char* prompt );
	static void ExportTextureCallBack( Visualization* v , const char* prompt );
	static void Init( void );
	static void InitializeVisualization( int width , int height );
	static void ComputeExactSolution( bool verbose=false );
	static void UpdateSolution( bool verbose=false , bool detailVerbose=false );
	static void InitializeSystem( int width , int height );

	static void Display( void ){ visualization.Display(); }
	static void MouseFunc( int button , int state , int x , int y );
	static void MotionFunc( int x , int y );
	static void Reshape( int w , int h ){ visualization.Reshape(w,h); }
	static void KeyboardFunc( unsigned char key , int x , int y ){ visualization.KeyboardFunc( key , x , y ); }
	static void Idle( void );
};

template< typename PreReal , typename Real > OrientedTexturedTriangleMesh< PreReal >						Geodesics< PreReal , Real >::mesh;
template< typename PreReal , typename Real > int															Geodesics< PreReal , Real >::textureWidth;
template< typename PreReal , typename Real > int															Geodesics< PreReal , Real >::textureHeight;
template< typename PreReal , typename Real > TexturedMeshVisualization										Geodesics< PreReal , Real >::visualization;
template< typename PreReal , typename Real > int															Geodesics< PreReal , Real >::mouseX = -1;
template< typename PreReal , typename Real > int															Geodesics< PreReal , Real >::mouseY = -1;
template< typename PreReal , typename Real > bool															Geodesics< PreReal , Real >::mouseSelectionActive = false;
template< typename PreReal , typename Real > Padding														Geodesics< PreReal , Real >::padding;

template< typename PreReal , typename Real > std::vector< AtlasChart< PreReal > >							Geodesics< PreReal , Real >::atlasCharts;
template< typename PreReal , typename Real > std::vector< std::vector< SquareMatrix< PreReal , 2 > > >		Geodesics< PreReal , Real >::parameterMetric;

template< typename PreReal , typename Real > SparseMatrix< Real , int >										Geodesics< PreReal , Real >::mass;
template< typename PreReal , typename Real > SparseMatrix< Real , int >										Geodesics< PreReal , Real >::stiffness;

template< typename PreReal , typename Real > SparseMatrix< Real , int >										Geodesics< PreReal , Real >::smoothImpulseMatrix;
template< typename PreReal , typename Real > SparseMatrix< Real , int >										Geodesics< PreReal , Real >::geodesicDistanceMatrix;

template< typename PreReal , typename Real > Real															Geodesics< PreReal , Real >::diffusionInterpolationWeight;
template< typename PreReal , typename Real > Real															Geodesics< PreReal , Real >::geodesicInterpolationWeight = 1e-6;

template< typename PreReal , typename Real > std::vector< TextureNodeInfo< PreReal > >						Geodesics< PreReal , Real >::textureNodes;
template< typename PreReal , typename Real > Image<int>														Geodesics< PreReal , Real >::nodeIndex;
template< typename PreReal , typename Real > std::vector< BilinearElementIndex >							Geodesics< PreReal , Real >::bilinearElementIndices;

template< typename PreReal , typename Real > int															Geodesics< PreReal , Real >::steps;
template< typename PreReal , typename Real > char															Geodesics< PreReal , Real >::stepsString[1024];

template< typename PreReal , typename Real > int															Geodesics< PreReal , Real >::levels;
template< typename PreReal , typename Real > HierarchicalSystem< PreReal , Real >							Geodesics< PreReal , Real >::hierarchy;

template< typename PreReal , typename Real > unsigned char *												Geodesics< PreReal , Real >::outputBuffer;
template< typename PreReal , typename Real > std::vector<MultigridLevelIndices<Real>>						Geodesics< PreReal , Real >::multigridIndices;

//Impulse Smoothing
template< typename PreReal , typename Real > std::vector< SystemCoefficients< Real > >						Geodesics< PreReal , Real >::multigridSmoothImpulseCoefficients;
template< typename PreReal , typename Real > std::vector<MultigridLevelVariables<Real>>						Geodesics< PreReal , Real >::multigridSmoothImpulseVariables;
template< typename PreReal , typename Real > VCycleSolvers< typename Geodesics< PreReal , Real >::DirectSolver > Geodesics< PreReal , Real >::smoothImpulseSolvers;

//Geodesic Distance
template< typename PreReal , typename Real > std::vector< SystemCoefficients< Real > >						Geodesics< PreReal , Real >::multigridGeodesicDistanceCoefficients;
template< typename PreReal , typename Real > std::vector<MultigridLevelVariables<Real>>						Geodesics< PreReal , Real >::multigridGeodesicDistanceVariables;
template< typename PreReal , typename Real > VCycleSolvers< typename Geodesics< PreReal , Real >::DirectSolver > Geodesics< PreReal , Real >::geodesicDistanceSolvers;

template< typename PreReal , typename Real > typename Geodesics< PreReal , Real >::DirectSolver				Geodesics< PreReal , Real >::fineSmoothImpulseSolver;
template< typename PreReal , typename Real > typename Geodesics< PreReal , Real >::DirectSolver				Geodesics< PreReal , Real >::fineGeodesicDistanceSolver;

//Samples
template< typename PreReal , typename Real > GradientElementSamples< Real >									Geodesics< PreReal , Real >::gradientSamples;
template< typename PreReal , typename Real > std::vector<InteriorCellLine>									Geodesics< PreReal , Real >::interiorCellLines;
template< typename PreReal , typename Real > std::vector<std::pair<int, int>>								Geodesics< PreReal , Real >::interiorCellLineIndex;

template< typename PreReal , typename Real > int															Geodesics< PreReal , Real >::impulseTexel = -1;
template< typename PreReal , typename Real > std::vector<Point3D< float > >									Geodesics< PreReal , Real >::textureNodePositions;

template< typename PreReal , typename Real > Real															Geodesics< PreReal , Real >::smoothImpulseRange;
template< typename PreReal , typename Real > Real															Geodesics< PreReal , Real >::geodesicDistanceRange;

template< typename PreReal , typename Real > SparseMatrix<Real, int>										Geodesics< PreReal , Real >::coarseBoundaryFineBoundaryProlongation;
template< typename PreReal , typename Real > SparseMatrix<Real, int>										Geodesics< PreReal , Real >::fineBoundaryCoarseBoundaryRestriction;

template< typename PreReal , typename Real > std::vector<Real>												Geodesics< PreReal , Real >::coarseBoundaryValues;
template< typename PreReal , typename Real > std::vector<Real>												Geodesics< PreReal , Real >::coarseBoundaryRHS;
template< typename PreReal , typename Real > std::vector<Real>												Geodesics< PreReal , Real >::fineBoundaryValues;
template< typename PreReal , typename Real > std::vector<Real>												Geodesics< PreReal , Real >::fineBoundaryRHS;

template< typename PreReal , typename Real > int															Geodesics< PreReal , Real >::updateCount = -1;


template< typename PreReal , typename Real >
void Geodesics< PreReal , Real >::ComputeExactSolution( bool verbose )
{
	Miscellany::Timer timer;

	//(1) Smoothing impulse	
	timer.reset();
	solve( fineSmoothImpulseSolver , multigridSmoothImpulseVariables[0].x , multigridSmoothImpulseVariables[0].rhs );
	if( verbose ) printf( "Smoothing impulse %.4f\n" , timer.elapsed() );

	//(1) Integrating vector field	
	const std::vector<int> & boundaryGlobalIndex = hierarchy.gridAtlases[0].boundaryGlobalIndex;

	if( verbose ) timer.reset();
#pragma omp parallel for
	for (int i = 0; i < boundaryGlobalIndex.size(); i++) coarseBoundaryValues[i] = multigridSmoothImpulseVariables[0].x[boundaryGlobalIndex[i]];
	coarseBoundaryFineBoundaryProlongation.Multiply(&coarseBoundaryValues[0], &fineBoundaryValues[0]);
	if( verbose ) printf("Coarse to fine %.4f \n" , timer.elapsed() );


	std::vector< Real >& fineGeodesicDistanceRHS = multigridGeodesicDistanceVariables[0].rhs;
	if( verbose ) timer.reset();
	auto VectorFunction = []( Point2D< Real > v , SquareMatrix< Real , 2 > tensor )
	{
		Point2D< Real > _v = tensor * v;
		Real len2 = Point2D< Real >::Dot( v , _v );
		if( len2>0 ) return -v / (Real)sqrt( len2 );
		else         return -v;
	};
	memset( &multigridGeodesicDistanceVariables[0].rhs[0] , 0 , multigridGeodesicDistanceVariables[0].rhs.size() * sizeof(Real) );
	memset( &fineBoundaryRHS[0] , 0 , fineBoundaryRHS.size() * sizeof(Real) );
	Integrate< Real >( interiorCellLines , gradientSamples , multigridSmoothImpulseVariables[0].x , fineBoundaryValues , VectorFunction , fineGeodesicDistanceRHS , fineBoundaryRHS );
	if( verbose ) printf( "Integrating normalized vector field %.4f\n" , timer.elapsed() );

	if( verbose ) timer.reset();
	fineBoundaryCoarseBoundaryRestriction.Multiply( &fineBoundaryRHS[0] , &coarseBoundaryRHS[0] );
#pragma omp parallel for
	for( int i=0 ; i<boundaryGlobalIndex.size() ; i++ ) fineGeodesicDistanceRHS[ boundaryGlobalIndex[i] ] += coarseBoundaryRHS[i];
	if( verbose ) printf( "Fine to coarse %.4f\n" , timer.elapsed() );


	//(3) Update geodesic distance solution	
	if( verbose ) timer.reset();
	solve( fineGeodesicDistanceSolver , multigridGeodesicDistanceVariables[0].x , fineGeodesicDistanceRHS );
	if( verbose ) printf( "Computing geodesic distance %.4f\n" , timer.elapsed() );

	Real expectedMinDistance = multigridGeodesicDistanceVariables[0].x[impulseTexel];

#pragma omp parallel for
	for( int i=0 ; i<multigridGeodesicDistanceVariables[0].x.size() ; i++ ) multigridGeodesicDistanceVariables[0].x[i] -= expectedMinDistance;
}

template< typename PreReal , typename Real >
void Geodesics< PreReal , Real >::UpdateOutputBuffer( const std::vector< Real > &solution )
{
	//Real expectedMaxDistance = solution[impulseTexel];

	Real attenuationRadius = (Real)0.2; //value within 0 and 0.5
	Real attenuationStart  = (Real)0.8 - attenuationRadius;
#pragma omp parallel for
	for (int i = 0; i < textureNodes.size(); i++) {
		int ci = textureNodes[i].ci;
		int cj = textureNodes[i].cj;
		int offset = 3 * (textureWidth*cj + ci);
		Real value = solution[i];
		Point3D< float > red( 255.f , 64.f , 64.f ) , blue( 64.f , 64.f , 255.f ) , green( 64.f , 255.f , 64.f ); ;
		Point3D< float > color = value < 1.f ? red * (1.f - value) + blue * value : blue * (2.f - value) + green * ( value-1.f );
		float scaledDistance = value * 60.f;
		int closestInteger = floor(scaledDistance);

		float residual = scaledDistance - float(closestInteger);
		if( closestInteger%2!=0 ) residual = 1.f - residual;
		float attenuationWeight = ( residual-attenuationStart ) / ( 2.f * attenuationRadius );
		attenuationWeight = std::min< float >( std::max< float >( 0 , attenuationWeight ), 1.f );
		attenuationWeight = attenuationWeight*attenuationWeight*(3.f - 2.f * attenuationWeight );
		color = color * ( 1.f-attenuationWeight );

		outputBuffer[offset + 0] = (unsigned char)(color[0]);
		outputBuffer[offset + 1] = (unsigned char)(color[1]);
		outputBuffer[offset + 2] = (unsigned char)(color[2]);
	}

	glBindTexture(GL_TEXTURE_2D, visualization.textureBuffer);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, textureWidth, textureHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)&outputBuffer[0]);
	glBindTexture(GL_TEXTURE_2D, 0);

	glutPostRedisplay();
}

template< typename PreReal , typename Real >
void Geodesics< PreReal , Real >::Idle( void )
{
	visualization.Idle();
	int selectedTexel = -1;
	if( mouseSelectionActive )
	{
		Point3D< float > selectedPoint;
		if( visualization.showMesh )
		{
			visualization.select(mouseX, mouseY, selectedPoint);
			float minDistance = FLT_MAX;
			for( int i=0 ; i<textureNodePositions.size() ; i++ )
			{
				float squaredDistance = Point3D< float >::SquareNorm( textureNodePositions[i]-selectedPoint );
				if( squaredDistance<minDistance ) minDistance = squaredDistance , selectedTexel = i;
			}
		}
		else
		{
			Point2D< float > ip = visualization.selectImagePos(mouseX, mouseY);
			//printf("Texture Coord %f %f \n", ip[0], ip[1]);
			int i = floor(ip[0] * float(nodeIndex.width()) - 0.5f);
			int j = floor((1.0 - ip[1]) * float(nodeIndex.height()) - 0.5f);
			//printf("Image pos %d %d \n", i, j);
			if (i >= 0 && i < nodeIndex.width() && j >= 0 && j < nodeIndex.height()) {
				selectedTexel = nodeIndex(i, j);
			}
		}

		if( impulseTexel!=selectedTexel && selectedTexel!=-1 )
		{
			impulseTexel = selectedTexel;
			steps = 0;
			memset( &multigridSmoothImpulseVariables[0].rhs[0] , 0 , multigridSmoothImpulseVariables[0].rhs.size() * sizeof(Real) );
			multigridSmoothImpulseVariables[0].rhs[impulseTexel] = 1.0;
			memset( &multigridSmoothImpulseVariables[0].x[0] , 0 , multigridSmoothImpulseVariables[0].x.size() * sizeof(Real) );
			memset( &multigridGeodesicDistanceVariables[0].x[0] , 0 , multigridGeodesicDistanceVariables[0].x.size() * sizeof(Real) );
		}
	}

	if( impulseTexel!=-1 )
	{
		if( updateCount && !visualization.promptCallBack )
		{
			UpdateSolution();
			if( updateCount>0 ) updateCount--;
			steps++;
			sprintf( stepsString , "Steps: %d" , steps );
			UpdateOutputBuffer( multigridGeodesicDistanceVariables[0].x );
		}
	}
	if( !visualization.promptCallBack ) glutPostRedisplay();	
}

template< typename PreReal , typename Real>
void Geodesics< PreReal , Real >::MouseFunc( int button , int state , int x , int y )
{
	visualization.newX = x; visualization.newY = y;
	visualization.rotating = visualization.scaling = visualization.panning = false;

	if (state == GLUT_DOWN && glutGetModifiers() & GLUT_ACTIVE_SHIFT) {

		int selectedTexel = -1;
		if (visualization.showMesh)
		{
			Point3D< float > selectedPoint;
			if (visualization.select(x, y, selectedPoint))
			{
				mouseSelectionActive = true;
				mouseX = x;
				mouseY = y;
				float minDistance = FLT_MAX;
				for( int i=0 ; i<textureNodePositions.size() ; i++ )
				{
					float squaredDistance = Point3D< float >::SquareNorm(textureNodePositions[i] - selectedPoint);
					if( squaredDistance<minDistance ) minDistance = squaredDistance , selectedTexel = i;
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
		if( selectedTexel!=-1 && selectedTexel!=impulseTexel )
		{
			impulseTexel = selectedTexel;
			steps = 0;
			memset(&multigridSmoothImpulseVariables[0].rhs[0], 0, multigridSmoothImpulseVariables[0].rhs.size() * sizeof(Real));
			multigridSmoothImpulseVariables[0].rhs[impulseTexel] = 1.0;
			memset(&multigridSmoothImpulseVariables[0].x[0], 0, multigridSmoothImpulseVariables[0].x.size() * sizeof(Real));
			memset(&multigridGeodesicDistanceVariables[0].x[0], 0, multigridGeodesicDistanceVariables[0].x.size() * sizeof(Real));
		}
		if( impulseTexel!=-1 )
		{
			if( UseDirectSolver.set )
			{
				Miscellany::Timer timer;
				ComputeExactSolution( DetailVerbose.set );
				if( DetailVerbose.set ) printf( "Exact solution %.4f \n" , timer.elapsed() );

				UpdateOutputBuffer(multigridGeodesicDistanceVariables[0].x);
			}
			glutPostRedisplay();
		}
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
void Geodesics< PreReal , Real >::MotionFunc( int x , int y )
{
	if (mouseSelectionActive) {
		mouseX = x;
		mouseY = y;
		glutPostRedisplay();
	}
	else
	{
		if( !visualization.showMesh )
		{
			visualization.oldX = visualization.newX, visualization.oldY = visualization.newY, visualization.newX = x, visualization.newY = y;

			if (visualization.panning) visualization.xForm.offset[0] -= (visualization.newX - visualization.oldX) / visualization.imageToScreenScale(), visualization.xForm.offset[1] += (visualization.newY - visualization.oldY) / visualization.imageToScreenScale();
			else
			{
				float dz = (float)pow( 1.1 , (float)( visualization.newY-visualization.oldY ) / 8.f );
				visualization.xForm.zoom *= dz;
			}

		}
		else
		{
			visualization.oldX = visualization.newX, visualization.oldY = visualization.newY, visualization.newX = x, visualization.newY = y;
			int screenSize = std::min< int >(visualization.screenWidth, visualization.screenHeight);
			float rel_x = (visualization.newX - visualization.oldX) / (float)screenSize * 2;
			float rel_y = (visualization.newY - visualization.oldY) / (float)screenSize * 2;

			float pRight = rel_x * visualization.zoom, pUp = -rel_y * visualization.zoom;
			float pForward = rel_y * visualization.zoom;
			float rRight = -rel_y, rUp = -rel_x;

			if     ( visualization.rotating ) visualization.camera.rotateUp( -rUp ) , visualization.camera.rotateRight( -rRight );
			else if( visualization.scaling  ) visualization.camera.translate( visualization.camera.forward*pForward );
			else if( visualization.panning  ) visualization.camera.translate( -( visualization.camera.right*pRight + visualization.camera.up*pUp ) );
		}
		glutPostRedisplay();
	}
}

template< typename PreReal , typename Real>
void Geodesics< PreReal , Real >::ExportTextureCallBack( Visualization * /*v*/ , const char* prompt )
{

	Image< Point3D< float > > outputImage;
	outputImage.resize( textureWidth , textureHeight );
	for( int i=0 ; i<outputImage.size() ; i++ ) outputImage[i] = Point3D< float >( outputBuffer[3*i] , outputBuffer[3*i+1] , outputBuffer[3*i+2]) / float(255.0);
	padding.unpad( outputImage );
#ifdef VARIABLE_SIZED_IMAGE_CHANNEL
	outputImage.template write< 8 >(prompt);
#else // !VARIABLE_SIZED_IMAGE_CHANNEL
	outputImage.write(prompt);
#endif // VARIABLE_SIZED_IMAGE_CHANNEL
}

template< typename PreReal , typename Real >
void Geodesics< PreReal , Real >::UpdateSolution( bool verbose , bool detailVerbose )
{
	Miscellany::Timer timer;

	// (1) Update smoothed input solution
	if( verbose ) timer.reset();
	VCycle( multigridSmoothImpulseVariables , multigridSmoothImpulseCoefficients , multigridIndices , smoothImpulseSolvers , detailVerbose , detailVerbose );
	if( verbose ) printf( "Smoothing impulse %.4f\n" , timer.elapsed() );

	// (2) Integrate normalized vector field
	const std::vector<int> & boundaryGlobalIndex = hierarchy.gridAtlases[0].boundaryGlobalIndex;

	if( verbose ) timer.reset();
#pragma omp parallel for
	for (int i = 0; i < boundaryGlobalIndex.size(); i++) coarseBoundaryValues[i] = multigridSmoothImpulseVariables[0].x[boundaryGlobalIndex[i]];
	coarseBoundaryFineBoundaryProlongation.Multiply(&coarseBoundaryValues[0], &fineBoundaryValues[0]);
	if( verbose ) printf( "Coarse to fine %.4f\n" , timer.elapsed() );

	if( verbose ) timer.reset();
	auto VectorFunction = []( Point2D< Real > v , SquareMatrix< Real , 2 > tensor )
	{
		Point2D< Real > _v = tensor * v;
		Real len2 = Point2D<Real>::Dot( v , _v );
		if( len2>0 ) return -v / (Real)sqrt( len2 );
		else         return -v;
	};

	memset( &multigridGeodesicDistanceVariables[0].rhs[0] , 0 , multigridGeodesicDistanceVariables[0].rhs.size() * sizeof(Real) );
	memset( &fineBoundaryRHS[0] , 0 , fineBoundaryRHS.size() * sizeof(Real) );
	Integrate< Real >( interiorCellLines , gradientSamples , multigridSmoothImpulseVariables[0].x , fineBoundaryValues , VectorFunction , multigridGeodesicDistanceVariables[0].rhs , fineBoundaryRHS );
	if( verbose ) printf( "Integrating normalized vector field %.4f \n" , timer.elapsed() );

	if( verbose ) timer.reset();
	fineBoundaryCoarseBoundaryRestriction.Multiply(&fineBoundaryRHS[0], &coarseBoundaryRHS[0]);
#pragma omp parallel for
	for (int i = 0; i < boundaryGlobalIndex.size(); i++) multigridGeodesicDistanceVariables[0].rhs[boundaryGlobalIndex[i]] += coarseBoundaryRHS[i];
	if( verbose ) printf( "Fine to coarse %.4f\n" , timer.elapsed() );


	// (3) Update geodesic distance solution	
	if( verbose ) timer.reset();
	VCycle( multigridGeodesicDistanceVariables , multigridGeodesicDistanceCoefficients , multigridIndices, geodesicDistanceSolvers , detailVerbose , detailVerbose );
	if( verbose ) printf( "Solving geodesic distance %.4f \n" , timer.elapsed() );

	Real expectedMinDistance = multigridGeodesicDistanceVariables[0].x[impulseTexel];

#pragma omp parallel for
	for (int i = 0; i < multigridGeodesicDistanceVariables[0].x.size(); i++) multigridGeodesicDistanceVariables[0].x[i] -= expectedMinDistance;
}

template< typename PreReal , typename Real >
void Geodesics< PreReal , Real >::InitializeSystem( int width , int height )
{
	Miscellany::Timer timer;
	MultigridBlockInfo multigridBlockInfo(MultigridBlockWidth.value, MultigridBlockHeight.value, MultigridPaddedWidth.value, MultigridPaddedHeight.value, 0);
	InitializeHierarchy( mesh , width , height , levels , textureNodes , bilinearElementIndices , hierarchy , atlasCharts , multigridBlockInfo , true , DetailVerbose.set );
	if( Verbose.set ) printf( "\tInitialized hierarchy: %.2f(s)\n" , timer.elapsed() );

	//Initialize node index
	nodeIndex.resize(width, height);
	for( int i=0 ; i<nodeIndex.size() ; i++ ) nodeIndex[i] = -1;
	for( int i=0 ; i<textureNodes.size() ; i++ ) nodeIndex( textureNodes[i].ci , textureNodes[i].cj ) = i;

	BoundaryProlongationData< Real > boundaryProlongation;
	InitializeBoundaryProlongationData( hierarchy.gridAtlases[0] , boundaryProlongation );

	SystemCoefficients< Real > massCoefficients;
	SystemCoefficients< Real > stiffnessCoefficients;

	std::vector< Point3D< Real > > __inputSignal;
	std::vector< Real > __texelToCellCoeffs;
	SparseMatrix< Real , int > __boundaryCellBasedStiffnessRHSMatrix[3];

	InitializeMetric( mesh , EMBEDDING_METRIC , atlasCharts , parameterMetric );

	timer.reset();
	{
		switch( MatrixQuadrature.value )
		{
		case  1: InitializeMassAndStiffness< 1>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
		case  3: InitializeMassAndStiffness< 3>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
		case  6: InitializeMassAndStiffness< 6>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
		case 12: InitializeMassAndStiffness<12>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
		case 24: InitializeMassAndStiffness<24>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
		case 32: InitializeMassAndStiffness<32>( massCoefficients , stiffnessCoefficients , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
		default: Miscellany::Throw( "Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles" );
		}
	}
	if( Verbose.set ) printf( "\tInitialized mass and stiffness: %.2f(s)\n" , timer.elapsed() );

	if( UseDirectSolver.set )
	{
		Miscellany::Timer tmr;
		FullMatrixConstruction( hierarchy.gridAtlases[0] , massCoefficients , mass );
		FullMatrixConstruction( hierarchy.gridAtlases[0] , stiffnessCoefficients , stiffness );
		smoothImpulseMatrix = mass * diffusionInterpolationWeight + stiffness;
		geodesicDistanceMatrix = mass * geodesicInterpolationWeight + stiffness;
		printf( "\tAssembled matrices: %.2f(s) \n" , tmr.elapsed() );
	}

//////////////////////////////////// Initialize multigrid indices

	multigridIndices.resize(levels);
	for (int i = 0; i < levels; i++)
	{
		const GridAtlas< PreReal , Real > &gridAtlas = hierarchy.gridAtlases[i];
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

	timer.reset();
	UpdateLinearSystem( diffusionInterpolationWeight , (Real)1. , hierarchy , multigridSmoothImpulseCoefficients , massCoefficients , stiffnessCoefficients , smoothImpulseSolvers , fineSmoothImpulseSolver , smoothImpulseMatrix , DetailVerbose.set , true , UseDirectSolver.set );
	UpdateLinearSystem(  geodesicInterpolationWeight , (Real)1. , hierarchy , multigridGeodesicDistanceCoefficients , massCoefficients , stiffnessCoefficients , geodesicDistanceSolvers , fineGeodesicDistanceSolver , geodesicDistanceMatrix , DetailVerbose.set , true , UseDirectSolver.set );
	if( Verbose.set ) printf( "\tInitialized multigrid coefficients: %.2f(s)\n" , timer.elapsed() );

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

	InitializeGridAtlasInteriorCellLines( hierarchy.gridAtlases[0].gridCharts , interiorCellLines , interiorCellLineIndex );
	if( interiorCellLineIndex.size()!=hierarchy.gridAtlases[0].numInteriorCells )
		Miscellany::Throw( "Inconsistent number of interior cells: %d!=%d" , hierarchy.gridAtlases[0].numInteriorCells , (int)interiorCellLineIndex.size() );

	coarseBoundaryFineBoundaryProlongation = boundaryProlongation.coarseBoundaryFineBoundaryProlongation;
	fineBoundaryCoarseBoundaryRestriction = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction;
	std::vector<int> fineBoundaryIndex = boundaryProlongation.fineBoundaryIndex;
	int numFineBoundarNodes = boundaryProlongation.numFineBoundarNodes;

	gradientSamples.resize( interiorCellLines.size() );

	timer.reset();
	{
		switch( VectorFieldQuadrature.value )
		{
		case  1: InitializeIntegration<  1 >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , gradientSamples , !PreciseIntegration.set ) ; break;
		case  3: InitializeIntegration<  3 >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , gradientSamples , !PreciseIntegration.set ) ; break;
		case  6: InitializeIntegration<  6 >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , gradientSamples , !PreciseIntegration.set ) ; break;
		case 12: InitializeIntegration< 12 >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , gradientSamples , !PreciseIntegration.set ) ; break;
		case 24: InitializeIntegration< 24 >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , gradientSamples , !PreciseIntegration.set ) ; break;
		case 32: InitializeIntegration< 32 >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , gradientSamples , !PreciseIntegration.set ) ; break;
		Miscellany::Throw( "Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles" );
		}
	}
	if( Verbose.set ) printf( "\tInitialized vector field integration: %.2f(s)\n" , timer.elapsed() );
	coarseBoundaryValues.resize(hierarchy.gridAtlases[0].numTexels - hierarchy.gridAtlases[0].numDeepTexels);
	coarseBoundaryRHS.resize(hierarchy.gridAtlases[0].numTexels - hierarchy.gridAtlases[0].numDeepTexels);
	fineBoundaryValues.resize(numFineBoundarNodes);
	fineBoundaryRHS.resize(numFineBoundarNodes);

	gradientSamples.sort();
}

 template< typename PreReal , typename Real>
void Geodesics< PreReal , Real >::InitializeVisualization( int width , int height )
{
	outputBuffer = new unsigned char[ height*width* 3];
	memset( outputBuffer , 204 , height * width * 3 * sizeof(unsigned char) );

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
	visualization.info.push_back( stepsString );

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

template< typename PreReal , typename Real >
void Geodesics< PreReal , Real >::ToggleUpdateCallBack( Visualization * /*v*/ , const char * /*prompt*/ )
{
	if( updateCount ) updateCount = 0;
	else              updateCount = -1;
}

template< typename PreReal , typename Real >
void Geodesics< PreReal , Real >::IncrementUpdateCallBack( Visualization * /*v*/ , const char * /*prompt*/ )
{
	if( updateCount<0 ) updateCount = 1;
	else updateCount++;
}

template< typename PreReal , typename Real >
void Geodesics< PreReal , Real >::Init( void )
{
	sprintf( stepsString , "Steps: 0" );
	levels = Levels.value;
	diffusionInterpolationWeight = DiffusionInterpolationWeight.value;
	textureWidth = Width.value;
	textureHeight = Height.value;

	mesh.read( Input.value , DetailVerbose.set );

	if( RandomJitter.set )
	{
		if( RandomJitter.value ) srand( RandomJitter.value );
		else                     srand( time(NULL) );
		std::vector< Point2D< PreReal > > randomOffset( mesh.vertices.size() );
		PreReal jitterScale = (PreReal)1e-3 / std::max< int >( textureWidth , textureHeight );
		for( int i=0 ; i<randomOffset.size() ; i++ ) randomOffset[i] = Point2D< PreReal >( (PreReal)1. - Random< PreReal >()*2 , (PreReal)1. - Random<PreReal>()*2 ) * jitterScale;
		for( int i=0 ; i<mesh.triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ ) mesh.textureCoordinates[ 3*i+k ] += randomOffset[ mesh.triangles[i][k] ];
	}

	{
		padding = Padding::Init( textureWidth , textureHeight , mesh.textureCoordinates , DetailVerbose.set );
		padding.pad( textureWidth , textureHeight , mesh.textureCoordinates );
		textureWidth  += padding.width();
		textureHeight += padding.height();
	}

	//Define centroid and scale for visualization
	Point3D< PreReal > centroid;
	for( int i=0 ; i<mesh.vertices.size() ; i++ ) centroid += mesh.vertices[i];
	centroid /= (int)mesh.vertices.size();
	PreReal radius = 0;
	for( int i=0 ; i<mesh.vertices.size() ; i++ ) radius = std::max< PreReal >( radius , Point3D< PreReal >::Length( mesh.vertices[i]-centroid) );
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

	textureNodePositions.resize( textureNodes.size() );
	for( int i=0 ; i<textureNodePositions.size() ; i++ )
	{
		Point2D< PreReal > baryncetricCoords = textureNodes[i].barycentricCoords;
		int tID = textureNodes[i].tID;
		Point3D< PreReal > p =
			mesh.vertices[ mesh.triangles[tID][0] ] * ( (PreReal)1.-baryncetricCoords[0]-baryncetricCoords[1] ) +
			mesh.vertices[ mesh.triangles[tID][1] ] *               baryncetricCoords[0]                        +
			mesh.vertices[ mesh.triangles[tID][2] ] *                                    baryncetricCoords[1]   ;
		textureNodePositions[i] = Point3D< float >( p );
	}
}

template< typename PreReal , typename Real >
void _main( int argc, char *argv[] )
{
	Geodesics< PreReal , Real >::Init();

	glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
	Geodesics< PreReal , Real >::visualization.displayMode = DisplayMode.value;
	if     ( DisplayMode.value==ONE_REGION_DISPLAY ) Geodesics< PreReal , Real >::visualization.screenWidth =  800 , Geodesics< PreReal , Real >::visualization.screenHeight = 800;
	else if( DisplayMode.value==TWO_REGION_DISPLAY ) Geodesics< PreReal , Real >::visualization.screenWidth = 1440 , Geodesics< PreReal , Real >::visualization.screenHeight = 720;

	glutInitWindowSize( Geodesics< PreReal , Real >::visualization.screenWidth , Geodesics< PreReal , Real >::visualization.screenHeight );

	glutInit( &argc , argv );
	char windowName[1024];
	sprintf( windowName , "Goedsics" );
	glutCreateWindow( windowName );
	if( glewInit()!=GLEW_OK ) Miscellany::Throw( "glewInit failed" );
	glutDisplayFunc ( Geodesics< PreReal , Real >::Display );
	glutReshapeFunc ( Geodesics< PreReal , Real >::Reshape );
	glutMouseFunc   ( Geodesics< PreReal , Real >::MouseFunc );
	glutMotionFunc  ( Geodesics< PreReal , Real >::MotionFunc );
	glutKeyboardFunc( Geodesics< PreReal , Real >::KeyboardFunc );
	if( !UseDirectSolver.set ) glutIdleFunc( Geodesics< PreReal , Real >::Idle );
	if( CameraConfig.set ) Geodesics< PreReal , Real >::visualization.ReadSceneConfigurationCallBack( &Geodesics< PreReal , Real >::visualization , CameraConfig.value );
	Geodesics< PreReal , Real >::InitializeVisualization( Geodesics< PreReal , Real >::textureWidth , Geodesics< PreReal , Real >::textureHeight );
	glutMainLoop();
}

int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , argv+1 , params );
	if( !Input.set ){ ShowUsage( argv[0] ) ; return EXIT_FAILURE; }
	omp_set_num_threads( Threads.value );
	if( !NoHelp.set )
	{
		printf( "+---------------------------------------------+\n" );
		printf( "| Interface Controls:                         |\n" );
		printf( "|    [Left Mouse]:                 rotate     |\n" );
		printf( "|    [Right Mouse]:                zoom       |\n" );
		printf( "|    [Left/Right Mouse] + [CTRL]:  pan        |\n" );
		printf( "|    [Left/Right Mouse] + [SHIFT]: set source |\n" );
		printf( "+---------------------------------------------+\n" );
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
