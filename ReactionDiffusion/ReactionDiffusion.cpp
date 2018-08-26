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
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDINDetailVerboseG NEGLIGENCE OR OTHERWISE) ARISING IN
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
#include <Src/QuadratureIntergration.inl>
#include <Src/MassAndStiffness.h>
#include <Src/Padding.h>
#include <Src/TexturedMeshVisualization.h>

const float StripeRates[] = { 0.062f , 0.062f };
const float    DotRates[] = { 0.0367f , 0.0649f };
const float DefaultStripesSamplesFraction = 0.01f;
const float DefaultDotsSamplesFraction = 0.001f;

cmdLineParameter< char* > Input( "in" );
cmdLineParameter< char* > Output( "out" );
cmdLineParameter< int   > OutputSteps( "outSteps" , 1000 );
cmdLineParameter< int   > Width( "width" , 512 );
cmdLineParameter< int   > Height( "height" , 512 );
cmdLineParameter< float > Speed( "speed" , 10.f );
cmdLineParameterArray< float , 2 > FeedKillRates( "fk" , StripeRates );
cmdLineParameter< float > DiffusionScale( "diff" , 1.f );
cmdLineParameter< float > SamplesFraction( "samples" );
cmdLineParameter< float > AreaScale( "areaScale" , 1.f );
cmdLineParameter< int   > Levels( "levels" , 4 );
cmdLineParameter< char* > CameraConfig( "camera" );
cmdLineParameter< int   > Threads( "threads" , omp_get_num_procs() );
cmdLineParameter< int   > DisplayMode( "display" , TWO_REGION_DISPLAY );
cmdLineParameter< int   > MatrixQuadrature( "mQuadrature" , 6 );
cmdLineParameter< int   > RHSQuadrature( "rhsQuadrature" , 3 );

cmdLineParameter< int   > MultigridBlockHeight ( "mBlockH" ,  16 );
cmdLineParameter< int   > MultigridBlockWidth  ( "mBlockW" , 128 );
cmdLineParameter< int   > MultigridPaddedHeight( "mPadH"   ,   0 );
cmdLineParameter< int   > MultigridPaddedWidth ( "mPadW"   ,   2 );

cmdLineReadable RandomJitter( "jitter" );
cmdLineReadable Verbose( "verbose" );
cmdLineReadable DetailVerbose( "detail" );
cmdLineReadable UseDirectSolver( "useDirectSolver" );
cmdLineReadable Double( "double" );
cmdLineReadable ApproximateIntegration( "approximateIntegration" );
cmdLineReadable Dots( "dots" );

cmdLineReadable* params[] =
{
	&Input , &Output , &OutputSteps , &Width , &Height , &Speed , &FeedKillRates , &DiffusionScale , &SamplesFraction , &CameraConfig , &Levels , &UseDirectSolver , &Threads , &DisplayMode , &MultigridBlockHeight , &MultigridBlockWidth , &MultigridPaddedHeight , &MultigridPaddedWidth ,
	&Verbose , &DetailVerbose , &AreaScale ,
	&RandomJitter ,
	&Double ,
	&MatrixQuadrature , &RHSQuadrature ,
	&ApproximateIntegration , &Dots ,
	NULL
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n", ex );
	printf( "\t --%s <input mesh>\n" , Input.name );
	printf( "\t[--%s <output texture>\n" , Output.name );
	printf( "\t[--%s <output steps>=%d]\n" , OutputSteps.name , OutputSteps.value );
	printf( "\t[--%s <texture width>=%d]\n" , Width.name , Width.value );
	printf( "\t[--%s <texture height>=%d]\n" , Height.name , Height.value );
	printf( "\t[--%s <time-step>=%f]\n" , Speed.name , Speed.value );
	printf( "\t[--%s <feed/kill rates>=%f %f]\n" , FeedKillRates.name , FeedKillRates.values[0] , FeedKillRates.values[1] );
	printf( "\t[--%s <diffusion scale>=%f]\n" , DiffusionScale.name , DiffusionScale.value );
	printf( "\t[--%s <samples fraction>=%f / %f]\n" , SamplesFraction.name , DefaultStripesSamplesFraction , DefaultDotsSamplesFraction );
	printf( "\t[--%s <area scale factor>=%f]\n" , AreaScale.name , AreaScale.value );
	printf( "\t[--%s <system matrix quadrature points per triangle>=%d]\n" , MatrixQuadrature.name , MatrixQuadrature.value );
	printf( "\t[--%s <right-hand-side quadrature points per triangle>=%d]\n" , RHSQuadrature.name , RHSQuadrature.value );
	printf( "\t[--%s]\n" , ApproximateIntegration.name );
	printf( "\t[--%s]\n" , UseDirectSolver.name );
	printf( "\t[--%s]\n" , RandomJitter.name );
	printf( "\t[--%s]\n" , Dots.name );
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
}

template< class Real >
class GrayScottReactionDiffusion
{
public:
	static TexturedMesh mesh;
	static int textureWidth;
	static int textureHeight;
	static double diffusionRates[2];
	static double kill;
	static double feed;
	static double speed;
	static double samplesFraction;
	static int levels;
	static int steps;
	static char stepsString[];
	static int whichConcentration;

	static Padding padding;

	static std::vector< Point3D< float > > textureNodePositions;

	static HierarchicalSystem hierarchy;

	static std::vector< BilinearElementIndex > bilinearElementIndices;

	static std::vector< TextureNodeInfo > textureNodes;
	static Image< int > nodeIndex;

	static SparseMatrix< double , int > mass;
	static SparseMatrix< double , int > stiffness;

	static SparseMatrix< double , int > systemMatrices[2];

	static int seedTexel;

	static std::vector< AtlasChart > atlasCharts;
	static std::vector< std::vector< SquareMatrix< double , 2 > > > parameterMetric;

	static std::vector< MultigridLevelCoefficients< Real > > multigridCoefficients[2];
	static std::vector< MultigridLevelVariables   < Real > > multigridVariables[2];

#if defined( USE_CHOLMOD )
	typedef  std::vector< CholmodCholeskySolver< Real , 1 > > BoundarySolverType;
	typedef  CholmodCholeskySolver< Real , 1 >  CoarseSolverType;
	typedef  CholmodCholeskySolver< Real , 1 > DirectSolverType;
#elif defined( USE_EIGEN_SIMPLICIAL )
	typedef  std::vector< EigenCholeskySolver< Real , 1 > > BoundarySolverType;
	typedef  EigenCholeskySolver< Real , 1 >  CoarseSolverType;
	typedef  EigenCholeskySolver< Real , 1 > DirectSolverType;
#elif defined( USE_EIGEN_PARDISO )
	typedef  std::vector< EigenPardisoSolver< Real , 1 > > BoundarySolverType;
	typedef  EigenPardisoSolver< Real , 1 > CoarseSolverType;
	typedef  EigenPardisoSolver< Real , 1 > DirectSolverType;
#else
#error "[ERROR] No solver defined!"
#endif

	static BoundarySolverType boundarySolvers[2];
	static CoarseSolverType coarseSolvers[2];
	static DirectSolverType fineSolvers[2];

	static std::vector< MultigridLevelIndices< Real > > multigridIndices;

	static SparseMatrix< Real , int > coarseBoundaryFineBoundaryProlongation;
	static SparseMatrix< Real , int > fineBoundaryCoarseBoundaryRestriction;
	static std::vector< Point2D< Real > > coarseBoundaryValues;
	static std::vector< Point2D< Real > > coarseBoundaryRHS;
	static std::vector< Point2D< Real > > fineBoundaryValues;
	static std::vector< Point2D< Real > > fineBoundaryRHS;

	//Samples
	static std::vector< QuadraticElementScalarSample< Real > > quadraticElementScalarSamples;
	static std::vector< std::vector< BilinearElementScalarSample< Real > > > bilinearElementScalarSamples;
	static std::vector< InteriorCellLine > interiorCellLines;
	static std::vector< std::pair< int , int > > interiorCellLineIndex;

	//Linear Operators
	static std::vector< double > deepMassCoefficients;
	static std::vector< double > deepStiffnessCoefficients;
	static SparseMatrix< double , int > boundaryBoundaryMassMatrix;
	static SparseMatrix< double , int > boundaryBoundaryStiffnessMatrix;
	static SparseMatrix< double , int > boundaryDeepMassMatrix;
	static SparseMatrix< double , int > boundaryDeepStiffnessMatrix;

	static unsigned char * outputBuffer;

	//Visulization
	static TexturedMeshVisualization visualization;
	static int mouseX , mouseY;
	static bool mouseSelectionActive;

	static void SetOutputBuffer( const std::vector< Real >& solution );
	static void UpdateOutputBuffer( const std::vector< Real >& solution );

	static int updateCount;

	static void SetConcentration1CallBack( Visualization* v , const char* prompt );
	static void SetConcentration2CallBack( Visualization* v , const char* prompt );
	static void ToggleUpdateCallBack( Visualization* v , const char* prompt );
	static void IncrementUpdateCallBack( Visualization* v , const char* prompt );
	static void ExportTextureCallBack( Visualization* v , const char* prompt );
	static int Init( void );
	static void InitializeVisualization( const int width , const int height );
	static int SetRightHandSide( bool verbose=false );
	static int UpdateExactSolution( bool verbose=false );
	static int UpdateApproximateSolution( bool verbose=false , bool detailVerbose=false );
	static int InitializeSystem( const int width , const int height , float areaScale );
	static int InitializeConcentrations( void );

	static void Display( void ){ visualization.Display(); }
	static void MouseFunc( int button , int state , int x , int y );
	static void MotionFunc( int x , int y );
	static void Reshape( int w , int h ){ visualization.Reshape(w,h); }
	static void KeyboardFunc( unsigned char key , int x , int y ){ visualization.KeyboardFunc( key , x , y ); }
	static void Idle( void );
};

template< class Real > TexturedMesh													GrayScottReactionDiffusion< Real >::mesh;
template< class Real > int															GrayScottReactionDiffusion< Real >::textureWidth;
template< class Real > int															GrayScottReactionDiffusion< Real >::textureHeight;
template< class Real > TexturedMeshVisualization									GrayScottReactionDiffusion< Real >::visualization;
template< class Real > int															GrayScottReactionDiffusion< Real >::mouseX = -1;
template< class Real > int															GrayScottReactionDiffusion< Real >::mouseY = -1;
template< class Real > bool															GrayScottReactionDiffusion< Real >::mouseSelectionActive = false;
template< class Real > Padding														GrayScottReactionDiffusion< Real >::padding;

template< class Real > std::vector< AtlasChart >									GrayScottReactionDiffusion< Real> ::atlasCharts;
template< class Real > std::vector< std::vector< SquareMatrix< double , 2 > > >		GrayScottReactionDiffusion< Real >::parameterMetric;
template< class Real > SparseMatrix< double , int >									GrayScottReactionDiffusion< Real >::mass;
template< class Real > SparseMatrix< double , int >									GrayScottReactionDiffusion< Real >::stiffness;

template< class Real > SparseMatrix< double , int >									GrayScottReactionDiffusion< Real >::systemMatrices[2];

template< class Real > double														GrayScottReactionDiffusion< Real >::diffusionRates[2];
template< class Real > double														GrayScottReactionDiffusion< Real >::speed;
template< class Real > double														GrayScottReactionDiffusion< Real >::kill;
template< class Real > double														GrayScottReactionDiffusion< Real >::feed;
template< class Real > double														GrayScottReactionDiffusion< Real >::samplesFraction;


template< class Real > std::vector< TextureNodeInfo >								GrayScottReactionDiffusion< Real >::textureNodes;
template< class Real > Image< int >													GrayScottReactionDiffusion< Real >::nodeIndex;
template< class Real > std::vector< BilinearElementIndex >							GrayScottReactionDiffusion< Real >::bilinearElementIndices;

template< class Real > int															GrayScottReactionDiffusion< Real >::steps;
template< class Real > char															GrayScottReactionDiffusion< Real >::stepsString[1024];
template< class Real > int															GrayScottReactionDiffusion< Real >::levels;
template< class Real > HierarchicalSystem											GrayScottReactionDiffusion< Real >::hierarchy;

template< class Real > unsigned char *												GrayScottReactionDiffusion< Real >::outputBuffer;
template< class Real > std::vector< MultigridLevelIndices< Real > >					GrayScottReactionDiffusion< Real >::multigridIndices;

template< class Real > std::vector< MultigridLevelCoefficients< Real > >			GrayScottReactionDiffusion< Real >::multigridCoefficients[2];
template< class Real > std::vector< MultigridLevelVariables< Real > >				GrayScottReactionDiffusion< Real >::multigridVariables[2];
template< class Real > typename GrayScottReactionDiffusion< Real >::CoarseSolverType		GrayScottReactionDiffusion< Real >::coarseSolvers[2];

template< class Real > typename GrayScottReactionDiffusion< Real >::DirectSolverType		GrayScottReactionDiffusion< Real >::fineSolvers[2];
template< class Real > typename GrayScottReactionDiffusion< Real >::BoundarySolverType		GrayScottReactionDiffusion< Real >::boundarySolvers[2];


//Samples
template< class Real > std::vector< QuadraticElementScalarSample< Real > >			GrayScottReactionDiffusion< Real >::quadraticElementScalarSamples;
template< class Real > std::vector< std::vector< BilinearElementScalarSample< Real > > >	GrayScottReactionDiffusion< Real >::bilinearElementScalarSamples;
template< class Real > std::vector< InteriorCellLine >								GrayScottReactionDiffusion< Real >::interiorCellLines;
template< class Real > std::vector< std::pair< int , int > >						GrayScottReactionDiffusion< Real >::interiorCellLineIndex;

template< class Real > int															GrayScottReactionDiffusion< Real >::seedTexel = -1;
template< class Real > std::vector< Point3D< float > >								GrayScottReactionDiffusion< Real >::textureNodePositions;

template< class Real > SparseMatrix< Real , int >									GrayScottReactionDiffusion< Real >::coarseBoundaryFineBoundaryProlongation;
template< class Real > SparseMatrix< Real , int >									GrayScottReactionDiffusion< Real >::fineBoundaryCoarseBoundaryRestriction;

template< class Real > std::vector< Point2D< Real > >								GrayScottReactionDiffusion< Real >::coarseBoundaryValues;
template< class Real > std::vector< Point2D< Real > >								GrayScottReactionDiffusion< Real >::coarseBoundaryRHS;
template< class Real > std::vector< Point2D< Real > >								GrayScottReactionDiffusion< Real >::fineBoundaryValues;
template< class Real > std::vector< Point2D< Real > >								GrayScottReactionDiffusion< Real >::fineBoundaryRHS;

template< class Real > std::vector< double >										GrayScottReactionDiffusion< Real >::deepMassCoefficients;
template< class Real > std::vector< double >										GrayScottReactionDiffusion< Real >::deepStiffnessCoefficients;
template< class Real > SparseMatrix< double , int >									GrayScottReactionDiffusion< Real >::boundaryBoundaryMassMatrix;
template< class Real > SparseMatrix< double , int >									GrayScottReactionDiffusion< Real >::boundaryBoundaryStiffnessMatrix;
template< class Real > SparseMatrix< double , int >									GrayScottReactionDiffusion< Real >::boundaryDeepMassMatrix;
template< class Real > SparseMatrix< double , int >									GrayScottReactionDiffusion< Real >::boundaryDeepStiffnessMatrix;

template< class Real > int															GrayScottReactionDiffusion< Real >::whichConcentration = 1;
template< class Real > int															GrayScottReactionDiffusion< Real >::updateCount = 0;

template< class Real >
int GrayScottReactionDiffusion< Real >::SetRightHandSide( bool verbose )
{
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// ( D + d * diff1 * S ) a[t+d] = D( a[t] + d * ( - a[t] * b[t] * b[t] + feed * ( 1 - a[t] ) ) )    //
	// ( D + d * diff2 * S ) b[t+d] = D( b[t] + d * (   a[t] * b[t] * b[t] - ( kill + feed ) * b[t] ) ) //
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	// [MK] Should pull this out and initialize it once
	static std::vector< Point2D< Real > > ab_x( multigridVariables[0][0].x.size() );
	static std::vector< Point2D< Real > > ab_rhs( multigridVariables[0][0].rhs.size() );
	for( int ab=0 ; ab<2 ; ab++ )
#pragma omp parallel for
		for( int i=0 ; i<ab_x.size() ; i++ ) ab_x[i][ab] = multigridVariables[ab][0].x[i];
	const std::vector< int >& boundaryGlobalIndex = hierarchy.gridAtlases[0].boundaryGlobalIndex;
#pragma omp parallel for
	for( int i=0 ; i<boundaryGlobalIndex.size() ; i++ ) coarseBoundaryValues[i] = ab_x[ boundaryGlobalIndex[i] ];
	coarseBoundaryFineBoundaryProlongation.Multiply( &coarseBoundaryValues[0] , &fineBoundaryValues[0] );

	for( int ab=0 ; ab<2 ; ab++ ) MultiplyBySystemMatrix_NoReciprocals( deepMassCoefficients , boundaryDeepMassMatrix , boundaryBoundaryMassMatrix , hierarchy.gridAtlases[0].boundaryGlobalIndex , hierarchy.gridAtlases[0].rasterLines , multigridVariables[ab][0].x , multigridVariables[ab][0].rhs );
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
	if( !Integrate< Real >( interiorCellLines , bilinearElementScalarSamples , quadraticElementScalarSamples , ab_x , fineBoundaryValues , ABFunction , ab_rhs , fineBoundaryRHS ) )
	{
		fprintf( stderr , "[ERROR] Unable to integrate concentrations!\n" );
	}
	fineBoundaryCoarseBoundaryRestriction.Multiply( &fineBoundaryRHS[0] , &coarseBoundaryRHS[0] );
	for( int ab=0 ; ab<2 ; ab++ )
	{
#pragma omp parallel for
		for( int i=0 ; i<ab_rhs.size() ; i++ ) multigridVariables[ab][0].rhs[i] += ab_rhs[i][ab];
#pragma omp parallel for
		for( int i=0 ; i<boundaryGlobalIndex.size() ; i++ ) multigridVariables[ab][0].rhs[ boundaryGlobalIndex[i] ] += coarseBoundaryRHS[i][ab];
	}
	return 1;
}

template< class Real >
int GrayScottReactionDiffusion< Real >::UpdateExactSolution( bool verbose )
{
	// (1) Compute the right-hand-sides
	{
		clock_t begin = clock();
		SetRightHandSide( verbose );
		if( verbose ) printf( "Integrating normalized vector field %.4f\n" , double(clock() - begin) / CLOCKS_PER_SEC );
	}

	// (2) Solve the linear systems
	{
		clock_t begin = clock();
		for( int ab=0 ; ab<2 ; ab++ )
		{
			solve( fineSolvers[ab] , multigridVariables[ab][0].x , multigridVariables[ab][0].rhs );
			for( int i=0 ; i<multigridVariables[ab][0].x.size() ; i++ ) multigridVariables[ab][0].x[i] = std::max< Real >( multigridVariables[ab][0].x[i] , 0 );
		}
		if( verbose ) printf( "Performed direct solve %.4f\n" , double(clock() - begin) / CLOCKS_PER_SEC );
	}
	return 1;
}
template< class Real >
int GrayScottReactionDiffusion< Real >::UpdateApproximateSolution( bool verbose , bool detailVerbose )
{
	// (1) Compute the right-hand-sides
	{
		clock_t begin = clock();
		SetRightHandSide( verbose );
		if( verbose ) printf( "Integrated constraints %.3f(s)\n" , double(clock() - begin) / CLOCKS_PER_SEC );
	}

	// (2) Solve the linear systems
	{
		clock_t begin = clock();
		for( int ab=0 ; ab<2 ; ab++ )
		{
			VCycle( multigridVariables[ab] , multigridCoefficients[ab] , multigridIndices , boundarySolvers[ab] , coarseSolvers[ab] , detailVerbose , detailVerbose );
			for( int i=0 ; i<multigridVariables[ab][0].x.size() ; i++ ) multigridVariables[ab][0].x[i] = std::max< Real >( multigridVariables[ab][0].x[i] , 0 );
		}
		if( verbose ) printf( "Performed v-cycle: %.3f(s)\n" , double(clock() - begin) / CLOCKS_PER_SEC );
	}

	return 1;
}

template< class Real >
void GrayScottReactionDiffusion< Real >::SetOutputBuffer( const std::vector< Real > & solution )
{
#pragma omp parallel for
	for( int i=0 ; i<textureNodes.size() ; i++ )
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
}
template< class Real >
void GrayScottReactionDiffusion< Real >::UpdateOutputBuffer( const std::vector< Real > & solution )
{
	SetOutputBuffer( solution );
	glBindTexture( GL_TEXTURE_2D , visualization.textureBuffer );
	glTexImage2D( GL_TEXTURE_2D , 0 , GL_RGBA , textureWidth , textureHeight , 0 , GL_LUMINANCE , GL_UNSIGNED_BYTE , (GLvoid*)&outputBuffer[0] );
	glBindTexture( GL_TEXTURE_2D , 0 );

	glutPostRedisplay();
}

template<class Real>
void GrayScottReactionDiffusion< Real >::Idle( void )
{
	if( updateCount )
	{
		if( UseDirectSolver.set ){ if( !UpdateExactSolution()       ) fprintf( stderr , "[WARNING] Exact update failed!\n" ); }
		else                     { if( !UpdateApproximateSolution() ) fprintf( stderr , "[WARNING] Approximate update failed!\n" ); }
		if( updateCount>0 ) updateCount--;
		steps++;
		sprintf( stepsString , "Steps: %d" , steps );
	}
	UpdateOutputBuffer( multigridVariables[whichConcentration][0].x );

}

template< class Real >
void GrayScottReactionDiffusion< Real >::MouseFunc( int button , int state , int x , int y )
{
	visualization.newX = x , visualization.newY = y;
	visualization.rotating = visualization.scaling = visualization.panning = false;

	if( state==GLUT_DOWN && glutGetModifiers() & GLUT_ACTIVE_SHIFT )
	{
		int selectedTexel = -1;
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
					float squaredDistance = Point3D< float >::SquareNorm( textureNodePositions[i] - selectedPoint );
					if( squaredDistance<minDistance ) minDistance = squaredDistance , selectedTexel = i;
				}
			}
		}
		else
		{
			Point2D< float > ip = visualization.selectImagePos( x , y );
			int i = floor( ip[0] * float( nodeIndex.width()) - 0.5f );
			int j = floor( (1.0-ip[1])*float( nodeIndex.height() ) - 0.5f );
			if( i>=0 && i<nodeIndex.width() && j>=0 && j<nodeIndex.height() ) mouseSelectionActive = true , selectedTexel = nodeIndex(i,j);
		}
		if( selectedTexel!=-1 && selectedTexel!=seedTexel )
		{
			seedTexel = selectedTexel;
			InitializeConcentrations();
		}
	}
	else
	{
		mouseSelectionActive = false;
		if( button==GLUT_LEFT_BUTTON )
		{
			if( glutGetModifiers() & GLUT_ACTIVE_CTRL ) visualization.panning = true;
			else                                        visualization.rotating = true;
		}
		else if( button==GLUT_RIGHT_BUTTON ) visualization.scaling = true;
	}
}

template< class Real >
void GrayScottReactionDiffusion< Real >::MotionFunc( int x , int y )
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

			int imageSize = std::min< int >( visualization.screenWidth , visualization.screenHeight );
			if( visualization.panning ) visualization.xForm.offset[0] -= ( visualization.newX - visualization.oldX ) / visualization.imageToScreenScale() , visualization.xForm.offset[1] += ( visualization.newY - visualization.oldY ) / visualization.imageToScreenScale();
			else
			{
				float dz = float( pow( 1.1 , double( visualization.newY - visualization.oldY ) / 8 ) );
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

#ifdef GLM_FORCE_RADIANS
			if      ( visualization.rotating ) visualization.camera.rotateUp( visualization.center , -rUp ) , visualization.camera.rotateRight( visualization.center , rRight );
#else // !GLM_FORCE_RADIANS
			if      ( visualization.rotating ) visualization.camera.rotateUp( visualization.center , rUp ) , visualization.camera.rotateRight( visualization.center , rRight );
#endif // GLM_FORCE_RADIANS
			else if( visualization.scaling   ) visualization.camera.moveForward( pForward);
			else if( visualization.panning   ) visualization.camera.moveRight( -pRight ) , visualization.camera.moveUp( -pUp );
		}
		glutPostRedisplay();
	}
}

template< class Real >
void GrayScottReactionDiffusion< Real >::ExportTextureCallBack( Visualization* v , const char* prompt )
{
	Image< Point3D< float > > outputImage;
	outputImage.resize( textureWidth , textureHeight );
	for( int i=0 ; i<outputImage.size() ; i++ ) outputImage[i] = Point3D< float >( outputBuffer[i] , outputBuffer[i] , outputBuffer[i] ) / 255.f;
	if( padding.nonTrivial ) UnpadImage( padding , outputImage );
	outputImage.write( prompt );
}

template< class Real >
int GrayScottReactionDiffusion< Real >::InitializeConcentrations( void )
{
	steps = 0;
	for( int i=0 ; i<multigridVariables[0][0].x.size() ; i++ ) multigridVariables[0][0].x[i] = 1;
	for( int i=0 ; i<multigridVariables[1][0].x.size() ; i++ ) multigridVariables[1][0].x[i] = 0;
	if( seedTexel!=-1 ) multigridVariables[1][0].x[seedTexel] = 1;
	else
	{
		const int SAMPLES = (int)( multigridVariables[1][0].x.size() * samplesFraction );
		for( int i=0 ; i<SAMPLES ; i++ ) multigridVariables[1][0].x[ ( multigridVariables[1][0].x.size() * i ) / SAMPLES ] = 1;
	}
	return 1;
}

template< class Real >
int GrayScottReactionDiffusion< Real >::InitializeSystem( const int width , const int height , float areaScale )
{
	clock_t t_begin;

	t_begin = clock();
	MultigridBlockInfo multigridBlockInfo( MultigridBlockWidth.value , MultigridBlockHeight.value , MultigridPaddedWidth.value , MultigridPaddedHeight.value , 0 );
	if( !InitializeHierarchy( mesh , width , height , levels , textureNodes , bilinearElementIndices , hierarchy , atlasCharts , multigridBlockInfo , true , DetailVerbose.set ) ){ fprintf( stderr , "[ERROR] Failed intialization!\n" ) ; return 0; }
	if( Verbose.set ) printf( "\tInitialized hierarchy: %.2f(s)\n" , double(clock() - t_begin) / CLOCKS_PER_SEC );

	//Initialize node index
	nodeIndex.resize( width , height );
	for( int i=0 ; i<nodeIndex.size() ; i++ ) nodeIndex[i] = -1;
	for( int i=0 ; i<textureNodes.size() ; i++ )
	{
		if( nodeIndex( textureNodes[i].ci , textureNodes[i].cj )!=-1 ) if( false ) fprintf( stderr , "[WARNING] Multiple nodes mapped to pixel %d %d!\n" , textureNodes[i].ci , textureNodes[i].cj );
		nodeIndex( textureNodes[i].ci , textureNodes[i].cj ) = i;
	}

	BoundaryProlongationData boundaryProlongation;
	if( !InitializeBoundaryProlongationData( hierarchy.gridAtlases[0] , boundaryProlongation ) ) { fprintf( stderr , "[ERROR] Failed boundary prolongation!\n" ) ; return 0; }

	std::vector< Point3D< double > > __inputSignal;
	std::vector< double > __texelToCellCoeffs;
	SparseMatrix< double , int > __boundaryCellBasedStiffnessRHSMatrix[3];

	if( !InitializeMetric( mesh , EMBEDDING_METRIC , atlasCharts , parameterMetric ) ){ fprintf( stderr , "[ERROR] Unable to initialize metric\n") ; return 0; }

	// Scale the metric so that the area is equal to the resolution
	for( int i=0 ; i<parameterMetric.size() ; i++ ) for( int j=0 ; j<parameterMetric[i].size() ; j++ ) parameterMetric[i][j] *= areaScale * textureNodes.size() / 2;

	t_begin = clock();
	{
		int ret = 0;
		switch( MatrixQuadrature.value )
		{
			case 1:  ret = InitializeMassAndStiffness< 1>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case 3:  ret = InitializeMassAndStiffness< 3>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case 6:  ret = InitializeMassAndStiffness< 6>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case 12: ret = InitializeMassAndStiffness<12>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case 24: ret = InitializeMassAndStiffness<24>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			case 32: ret = InitializeMassAndStiffness<32>( deepMassCoefficients , deepStiffnessCoefficients , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix , hierarchy , parameterMetric , atlasCharts , boundaryProlongation , false , __inputSignal , __texelToCellCoeffs , __boundaryCellBasedStiffnessRHSMatrix ) ; break;
			default: fprintf( stderr , "[ERROR] Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles\n" );
		}
		if( !ret ){ fprintf( stderr , "[ERROR] Failed intialization!\n" ) ; return 0; }
	}
	if( Verbose.set ) printf( "\tInitialized mass and stiffness: %.2f(s)\n" , double(clock() - t_begin) / CLOCKS_PER_SEC );

	if( UseDirectSolver.set )
	{
		clock_t t_begin = clock();
		FullMatrixConstruction( hierarchy.gridAtlases[0] , deepMassCoefficients , boundaryBoundaryMassMatrix , boundaryDeepMassMatrix , mass );
		FullMatrixConstruction( hierarchy.gridAtlases[0] , deepStiffnessCoefficients , boundaryBoundaryStiffnessMatrix , boundaryDeepStiffnessMatrix , stiffness );
		systemMatrices[0] = mass + stiffness * diffusionRates[0] * speed;
		systemMatrices[1] = mass + stiffness * diffusionRates[1] * speed;
		if( Verbose.set ) printf( "\tAssembled matrices: %.2f(s)\n" , double(clock() - t_begin) / CLOCKS_PER_SEC );
	}

	//////////////////////////////////// Initialize multigrid indices

	multigridIndices.resize( levels );
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

	//////////////////////////////////// Initialize multigrid coefficients

	t_begin = clock();
	for( int ab=0 ; ab<2 ; ab++ )
		if( !UpdateLinearSystem( 1.0 , diffusionRates[ab] * speed , hierarchy , multigridCoefficients[ab] ,
			deepMassCoefficients , deepStiffnessCoefficients ,
			boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix ,
			boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix ,
			coarseSolvers[ab] , boundarySolvers[ab] , fineSolvers[ab] ,
			systemMatrices[ab] , DetailVerbose.set , true , UseDirectSolver.set ) )
		{ fprintf( stderr , "[ERROR] Failed system update!\n" ) ; return 0; }
	if( Verbose.set ) printf( "\tInitialized multigrid coefficients: %.2f(s)\n" , double(clock() - t_begin) / CLOCKS_PER_SEC );

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

	if( !InitializeGridAtlasInteriorCellLines( atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLines , interiorCellLineIndex ) ){ fprintf( stderr , "[ERROR] Unable to initialize interior cell lines!\n" ) ; return 0; }
	if( interiorCellLineIndex.size()!=hierarchy.gridAtlases[0].numInteriorCells ){ fprintf( stderr , "[ERROR] Inconsistent number of interior cells! Expected %d . Result %d.\n" , hierarchy.gridAtlases[0].numInteriorCells , (int)interiorCellLineIndex.size() ) ; return 0; }

	coarseBoundaryFineBoundaryProlongation = boundaryProlongation.coarseBoundaryFineBoundaryProlongation;
	fineBoundaryCoarseBoundaryRestriction = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction;
	std::vector< int > fineBoundaryIndex = boundaryProlongation.fineBoundaryIndex;
	int numFineBoundarNodes = boundaryProlongation.numFineBoundarNodes;

	bilinearElementScalarSamples.resize( interiorCellLines.size() );

	t_begin = clock();
	{
		int ret = 0;
		switch( RHSQuadrature.value )
		{
		case  1: ret = InitializeIntegration<  1 , Real >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , bilinearElementScalarSamples , quadraticElementScalarSamples , ApproximateIntegration.set ) ; break;
		case  3: ret = InitializeIntegration<  3 , Real >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , bilinearElementScalarSamples , quadraticElementScalarSamples , ApproximateIntegration.set ) ; break;
		case  6: ret = InitializeIntegration<  6 , Real >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , bilinearElementScalarSamples , quadraticElementScalarSamples , ApproximateIntegration.set ) ; break;
		case 12: ret = InitializeIntegration< 12 , Real >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , bilinearElementScalarSamples , quadraticElementScalarSamples , ApproximateIntegration.set ) ; break;
		case 24: ret = InitializeIntegration< 24 , Real >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , bilinearElementScalarSamples , quadraticElementScalarSamples , ApproximateIntegration.set ) ; break;
		case 32: ret = InitializeIntegration< 32 , Real >( parameterMetric , atlasCharts , hierarchy.gridAtlases[0].gridCharts , interiorCellLineIndex , fineBoundaryIndex , bilinearElementScalarSamples , quadraticElementScalarSamples , ApproximateIntegration.set ) ; break;
		default: fprintf( stderr , "[ERROR] Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles\n" );
		}
		if( !ret ){ fprintf( stderr , "[ERROR] Unable to initialize vector field integration samples!\n" ) ; return 0; }
	}
	if( Verbose.set ) printf( "\tInitialized vector field integration: %.2f(s)\n" , double(clock() - t_begin) / CLOCKS_PER_SEC );
	coarseBoundaryValues.resize( hierarchy.gridAtlases[0].numTexels - hierarchy.gridAtlases[0].numDeepTexels );
	coarseBoundaryRHS.resize   ( hierarchy.gridAtlases[0].numTexels - hierarchy.gridAtlases[0].numDeepTexels );
	fineBoundaryValues.resize( numFineBoundarNodes );
	fineBoundaryRHS.resize   ( numFineBoundarNodes );

	for( int i=0 ; i<bilinearElementScalarSamples.size() ; i++ ) std::sort( bilinearElementScalarSamples[i].begin() , bilinearElementScalarSamples[i].end() , BilinearElementScalarSample< Real >::Compare );

	int numTexels = hierarchy.gridAtlases[0].numTexels;
	int numFineNodes = hierarchy.gridAtlases[0].numFineNodes;

	return 1;
}

template< class Real >
void GrayScottReactionDiffusion< Real >::InitializeVisualization( const int width , const int height )
{
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
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , ' ' , "toggle update" , ToggleUpdateCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , '+' , "increment update" , IncrementUpdateCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , '2' , "show concentration 2" , SetConcentration2CallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , '1' , "show concentration 1" , SetConcentration1CallBack ) );

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
template< class Real > void GrayScottReactionDiffusion< Real >::SetConcentration1CallBack( Visualization* , const char* ){ whichConcentration = 0; }
template< class Real > void GrayScottReactionDiffusion< Real >::SetConcentration2CallBack( Visualization* , const char* ){ whichConcentration = 1; }
template< class Real > void GrayScottReactionDiffusion< Real >::ToggleUpdateCallBack( Visualization* , const char* )
{
	if( updateCount ) updateCount = 0;
	else              updateCount = -1;
}
template< class Real > void GrayScottReactionDiffusion< Real >::IncrementUpdateCallBack( Visualization* , const char* )
{
	if( updateCount<0 ) updateCount = 1;
	else updateCount++;
}

template< class Real >
int GrayScottReactionDiffusion< Real >::Init( void )
{
	sprintf( stepsString , "Steps: 0" );
	levels = Levels.value;
	feed = FeedKillRates.values[0];
	kill = FeedKillRates.values[1];
	speed = Speed.value;
	samplesFraction = SamplesFraction.value;
	diffusionRates[0] = 1.0;
	diffusionRates[1] = 0.5;
	for( int ab=0 ; ab<2 ; ab++ ) diffusionRates[ab] *= DiffusionScale.value / 10.;
	textureWidth = Width.value;
	textureHeight = Height.value;

	if( !ReadTexturedMesh( mesh , Input.value , NULL , DetailVerbose.set ) ){ fprintf( stderr , "[ERROR] Unable to read mesh data\n" ) ; return 0; }

	if( true ) for( int i=0 ; i<mesh.textureCoordinates.size() ; i++ ) mesh.textureCoordinates[i][1] = 1.0 - mesh.textureCoordinates[i][1];

	if( RandomJitter.set )
	{
		srand( time( NULL ) );
		std::vector< Point2D< double > > randomOffset( mesh.vertices.size() );
		double jitterScale = 1e-3 / double( std::max< int >( textureWidth , textureHeight ) );
		for( int i=0 ; i<randomOffset.size() ; i++ ) randomOffset[i] = Point2D< double >(1.0 - 2.0 * double( rand() ) / double(RAND_MAX) , 1.0 - 2.0 *  double( rand() ) / double(RAND_MAX) )*jitterScale;
		for( int i=0 ; i<mesh.triangles.size() ; i++ ) for( int k=0 ; k<3 ; k++ ) mesh.textureCoordinates[3*i+k] += randomOffset[ mesh.triangles[i][k] ];
	}

	ComputePadding( padding , textureWidth , textureHeight , mesh.textureCoordinates , DetailVerbose.set );
	if( padding.nonTrivial )
	{
		PaddTextureCoordinates( padding , textureWidth , textureHeight , mesh.textureCoordinates );
		textureWidth  += ( padding.left   + padding.right );
		textureHeight += ( padding.bottom + padding.top   );
	}

	// Define centroid and scale for visualization
	Point3D< double > centroid;
	for( int i=0 ; i<mesh.vertices.size() ; i++ ) centroid += mesh.vertices[i];
	centroid /= (double)mesh.vertices.size();
	double radius = 0;
	for( int i=0 ; i<mesh.vertices.size() ; i++ ) radius = std::max< double >( radius , Point3D< double >::Length( mesh.vertices[i]-centroid ) );
	for( int i=0 ; i<mesh.vertices.size() ; i++ ) mesh.vertices[i] = ( mesh.vertices[i]-centroid ) / radius;


	clock_t t = clock();
	if( !InitializeSystem( textureWidth , textureHeight , AreaScale.value ) ){ fprintf( stderr , "[ERROR] Failed to initialize system\n") ; return 0; }
	if( !InitializeConcentrations() ){ fprintf( stderr , "[ERROR] Failed to initialize concentrations\n") ; return 0; }

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

	textureNodePositions.resize( textureNodes.size() );
	for( int i=0 ; i<textureNodePositions.size() ; i++ )
	{
		Point2D< double > barycentricCoords = textureNodes[i].barycentricCoords;
		int tId = textureNodes[i].tId;
		Point3D<float> surfacePosition =
			mesh.vertices[ mesh.triangles[tId][0] ] * ( 1.0-barycentricCoords[0]-barycentricCoords[1] ) +
			mesh.vertices[ mesh.triangles[tId][1] ] * barycentricCoords[0] +
			mesh.vertices[ mesh.triangles[tId][2] ] * barycentricCoords[1];
		textureNodePositions[i] = surfacePosition;
	}

	outputBuffer = new unsigned char[ textureHeight*textureWidth];
	memset( outputBuffer , 204 , textureHeight * textureWidth * sizeof(unsigned char) );

	return 1;
}

template< class Real >
int _main( int argc , char* argv[] )
{
	if( !GrayScottReactionDiffusion< Real >::Init() ) return 0;

	if( !Output.set )
	{
		glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
		GrayScottReactionDiffusion< Real >::visualization.displayMode = DisplayMode.value;
		if     ( DisplayMode.value==ONE_REGION_DISPLAY ) GrayScottReactionDiffusion< Real >::visualization.screenWidth =  800 , GrayScottReactionDiffusion< Real >::visualization.screenHeight = 800;
		else if( DisplayMode.value==TWO_REGION_DISPLAY ) GrayScottReactionDiffusion< Real >::visualization.screenWidth = 1440 , GrayScottReactionDiffusion< Real >::visualization.screenHeight = 720;

		GrayScottReactionDiffusion< Real >::visualization.UpdateMainFrameSize();
		glutInitWindowSize( GrayScottReactionDiffusion< Real >::visualization.screenWidth , GrayScottReactionDiffusion< Real >::visualization.screenHeight );

		glutInit(&argc, argv);
		char windowName[1024];
		sprintf( windowName , "Gray-Scott Reaction Diffusion");
		glutCreateWindow( windowName );
		if( glewInit()!=GLEW_OK ) fprintf(stderr, "[ERROR] glewInit failed\n" ) , exit(0);
		glutDisplayFunc ( GrayScottReactionDiffusion< Real >::Display );
		glutReshapeFunc ( GrayScottReactionDiffusion< Real >::Reshape );
		glutMouseFunc   ( GrayScottReactionDiffusion< Real >::MouseFunc );
		glutMotionFunc  ( GrayScottReactionDiffusion< Real >::MotionFunc );
		glutKeyboardFunc( GrayScottReactionDiffusion< Real >::KeyboardFunc );
		glutIdleFunc    ( GrayScottReactionDiffusion< Real >::Idle );
		if( CameraConfig.set ) GrayScottReactionDiffusion< Real >::visualization.ReadSceneConfigurationCallBack( &GrayScottReactionDiffusion< Real >::visualization , CameraConfig.value );
		GrayScottReactionDiffusion< Real >::InitializeVisualization( GrayScottReactionDiffusion< Real >::textureWidth , GrayScottReactionDiffusion<Real>::textureHeight );
		glutMainLoop();
	}
	else
	{
		clock_t t = clock();
		for( int i=0 ; i<OutputSteps.value ; i++ )
		{
			if( Verbose.set ) printf( "%d / %d \r" , i+1 , OutputSteps.value );
			if( UseDirectSolver.set ) GrayScottReactionDiffusion< Real >::UpdateExactSolution();
			else                      GrayScottReactionDiffusion< Real >::UpdateApproximateSolution();
		}
		if( Verbose.set )
		{
			printf( "\n" );
			double total_time = double( clock()-t ) / CLOCKS_PER_SEC;
			printf( "Reaction-diffusion total time / time per iteration: %.2f(s) / %.4f(s)\n" , total_time , total_time / OutputSteps.value );
		}
		GrayScottReactionDiffusion< Real >::SetOutputBuffer( GrayScottReactionDiffusion< Real >::multigridVariables[1][0].x );
		GrayScottReactionDiffusion< Real >::ExportTextureCallBack( &GrayScottReactionDiffusion<Real>::visualization , Output.value );
	}
	return 0;

}

int main( int argc , char* argv[] )
{
	cmdLineParse( argc-1 , argv+1 , params );
	if( !Input.set ){ ShowUsage( argv[0] ) ; return EXIT_FAILURE; }
	if( !SamplesFraction.set )
	{
		if( Dots.set ) SamplesFraction.value = DefaultDotsSamplesFraction;
		else           SamplesFraction.value = DefaultStripesSamplesFraction;
	}
	if( Dots.set && !FeedKillRates.set ) FeedKillRates.values[0] = DotRates[0] , FeedKillRates.values[1] = DotRates[1];
	omp_set_num_threads( Threads.value );
	if( Double.set ) _main< double >( argc , argv );
	else             _main< float  >( argc , argv );
	return 0;
}