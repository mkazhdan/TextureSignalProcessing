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

#ifdef USE_EIGEN_PARDISO
#include <Eigen/PardisoSupport>
#endif // USE_EIGEN_PARDISO

#include <Misha/CmdLineParser.h> 
#include <Misha/Miscellany.h>
#include <Misha/Exceptions.h>
#include <Misha/FEM.h>
#include <Misha/MultiThreading.h>
#include <Misha/AutoDiff/EquationParser.h>
#include <Misha/SquaredEDT.h>
#include <Src/Hierarchy.h>
#include <Src/SimpleTriangleMesh.h>
#include <Src/Basis.h>
#include <Src/EigenSolverWrapper.h>
#include <Src/Operators.h>
#include <Src/Padding.h>
#include <Src/ImageIO.h>
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
#include <Src/StitchingVisualization.h>
#endif // NO_OPEN_GL_VISUALIZATION

using namespace MishaK;
using namespace MishaK::TSP;

enum struct WeightingType
{
	Metric ,
	Confidence ,
	GradientScale
};

const std::vector< std::string > WeightingTypeNames = { "metric" , "confidence" , "gradient scale" };

CmdLineParameterArray< std::string , 2 >
	In( "in" );

CmdLineParameter< std::string >
	InMask( "mask" ),
	InputLowFrequency( "inLow" ) ,
	ConfidenceFilter( "confFilter" ) ,
	Output( "out" );

CmdLineParameter< unsigned int >
	WeightType( "wType" , static_cast< unsigned int >( WeightingType::GradientScale ) ) ,
	Levels( "levels" , 4 ) ,
	OutputVCycles( "outVCycles" , 6 ) ,
	MatrixQuadrature( "mQuadrature" , 6 ) ,
	VectorQuadrature( "vQuadrature" , 6 ) ,
	MultigridBlockHeight ( "mBlockH" ,  16 ) ,
	MultigridBlockWidth  ( "mBlockW" , 128 ) ,
	MultigridPaddedHeight( "mPadH"   ,   0 ) ,
	MultigridPaddedWidth ( "mPadW"   ,   2 ) ,
	ChartMaskErode( "erode" , 0 ) ,
	RandomJitter( "jitter" , 0 );

CmdLineParameter< double >
	GradientFittingWeight( "gWeight" , 1e-2 ) ,
	MinimumConfidence( "minConf" , 0. ) ,
	ChartMaskErosionEps( "erodeEps" , 0. ) ,
	ConfidenceRegularizationWeight( "confRegWeight" , 1e-5 ) ,
	CollapseEpsilon( "collapse" , 0 );


#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
CmdLineParameter< std::string >
	CameraConfig( "camera" );
#endif // NO_OPEN_GL_VISUALIZATION

CmdLineReadable
	UseDirectSolver( "useDirectSolver" ) ,
	SanityCheck( "sanityCheck" ) ,
	NoHelp( "noHelp" ) ,
	DetailVerbose( "detail" ) ,
	Double( "double" ) ,
	MultiInput( "multi" ) ,
	Serial( "serial" ) ,
	Run( "run" ) ,
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	Nearest( "nearest" ) ,
#endif // NO_OPEN_GL_VISUALIZATION
	SingleView( "single" ) ,
	Verbose( "verbose" );

std::vector< CmdLineReadable * > params =
{
	&In ,
	&InMask ,
	&Output ,
	&GradientFittingWeight ,
	&Levels ,
	&UseDirectSolver ,
	&Serial,
	&Verbose ,
	&InputLowFrequency ,
	&DetailVerbose ,
	&MultigridBlockHeight ,
	&MultigridBlockWidth ,
	&MultigridPaddedHeight ,
	&MultigridPaddedWidth ,
	&RandomJitter ,
	&Double ,
	&MatrixQuadrature ,
	&VectorQuadrature ,
	&WeightType ,
	&MinimumConfidence ,
	&OutputVCycles ,
	&NoHelp ,
	&MultiInput ,
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	&CameraConfig ,
#endif // NO_OPEN_GL_VISUALIZATION
	&ChartMaskErode ,
#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	&Nearest ,
#endif // NO_OPEN_GL_VISUALIZATION
	&SingleView ,
	&CollapseEpsilon ,
	&SanityCheck ,
	&ChartMaskErosionEps ,
	&ConfidenceRegularizationWeight ,
	&ConfidenceFilter ,
	&Run
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
	printf( "\t[--%s <chart mask erosion epsilon>=%f]\n" , ChartMaskErosionEps.name.c_str() , ChartMaskErosionEps.value );

	printf( "\t[--%s <output v-cycles>=%d]\n" , OutputVCycles.name.c_str() , OutputVCycles.value );
	printf( "\t[--%s <gradient fitting weight>=%f]\n" , GradientFittingWeight.name.c_str() , GradientFittingWeight.value );
	printf( "\t[--%s <system matrix quadrature points per triangle>=%d]\n" , MatrixQuadrature.name.c_str(), MatrixQuadrature.value );
	printf( "\t[--%s <system vector quadrature points per triangle>=%d]\n" , VectorQuadrature.name.c_str(), VectorQuadrature.value );
	printf( "\t[--%s <weight type>=%d]\n" , WeightType.name.c_str() , WeightType.value );
	for( unsigned int i=0 ; i<WeightingTypeNames.size() ; i++ ) printf( "\t\t%d] %s\n" , i , WeightingTypeNames[i].c_str() );
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
	printf( "\t[--%s <hierarchy levels>=%d]\n" , Levels.name.c_str(), static_cast< int >(Levels.value) );
	printf( "\t[--%s <multigrid block width>=%d]\n" , MultigridBlockWidth.name.c_str() , MultigridBlockWidth.value );
	printf( "\t[--%s <multigrid block height>=%d]\n" , MultigridBlockHeight.name.c_str() , MultigridBlockHeight.value );
	printf( "\t[--%s <multigrid padded width>=%d]\n" , MultigridPaddedWidth.name.c_str() , MultigridPaddedWidth.value );
	printf( "\t[--%s <multigrid padded height>=%d]\n" , MultigridPaddedHeight.name.c_str() , MultigridPaddedHeight.value );
	printf( "\t[--%s <collapse epsilon>=%g]\n" , CollapseEpsilon.name.c_str() , CollapseEpsilon.value );
	printf( "\t[--%s <minimum confidence>=%g]\n" , MinimumConfidence.name.c_str() , MinimumConfidence.value );
	printf( "\t[--%s <confidence regularization weight>=%g]\n" , ConfidenceRegularizationWeight.name.c_str() , ConfidenceRegularizationWeight.value );
	printf( "\t[--%s <confidence filter equation (as a function in variable 's')>]\n" , ConfidenceFilter.name.c_str() );
	printf( "\t[--%s]\n", SingleView.name.c_str() );
	printf( "\t[--%s]\n", Serial.name.c_str() );
	printf( "\t[--%s]\n", Run.name.c_str() );
	printf( "\t[--%s]\n", DetailVerbose.name.c_str() );
	printf( "\t[--%s]\n" , SanityCheck.name.c_str() );
	printf( "\t[--%s]\n" , NoHelp.name.c_str() );
}

template< typename ActiveTexelFunctor /* = std::function< bool ( typename RegularGrid< 2 >::Index ) > */ >
RegularGrid< 2 , unsigned int > SquareDistanceToTextureBoundary( unsigned int width , unsigned int height , ActiveTexelFunctor && activeF )
{
	static_assert( std::is_convertible_v< ActiveTexelFunctor , std::function< bool ( typename RegularGrid< 2 >::Index ) > > , "[ERROR] ActiveTexelFunctor poorly formed" );
	using Index = typename RegularGrid< 2 >::Index;
	using Range = typename RegularGrid< 2 >::Range; 


	unsigned int res[] = { width , height };
	RegularGrid< 2 , unsigned int > d2( res );
	RegularGrid< 2 , unsigned char > raster( res ) , boundary( res );

	Range range;
	for( unsigned int d=0 ; d<2 ; d++ ) range.second[d] = res[d];

	// Identify all active texels
	range.processParallel( [&]( unsigned int , Index I ){ raster( I ) = activeF( I ) ? 1 : 0; } );

	// Identify the boundary texels
	{
		auto Kernel = [&]( unsigned int , Index I )
		{
			if( raster( I ) ) boundary( I ) = 0;
			else
			{
				bool hasActiveNeighbors = false;
				Range::Intersect( range , Range(I).dilate(1) ).process( [&]( Index I ){ hasActiveNeighbors |= raster(I)!=0; } );
				boundary( I ) = hasActiveNeighbors ? 1 : 0;
			}
		};
		range.processParallel( Kernel );
	}

	d2 = SquaredEDT< double , 2 >::Saito( boundary );

	ThreadPool::ParallelFor( 0 , d2.size() , [&]( size_t i ){ if( !raster[i] ) d2[i] = 0; } );

	return d2;
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
class Stitching
{
public:
	static int inputMode;
	static TexturedTriangleMesh< PreReal > mesh;
	static int textureWidth;
	static int textureHeight;
	static Real gradientFittingWeight;
	static unsigned int levels;

	static HierarchicalSystem< PreReal , Real > hierarchy;
	static bool rhsNeedsUpdating;
	static bool positiveModulation;

	// Single input mode
	static RegularGrid< 2 , int > inputMask;
	static RegularGrid< 2 , Point3D< Real > > lowFrequencyTexture;
	static RegularGrid< 2 , Point3D< Real > > inputComposition;
	static RegularGrid< 2 , Point3D< Real > > inputColorMask;

	// Multiple input mode
	static int numTextures;
	static std::vector< RegularGrid< 2 , Real > > inputTexelConfidence;
	static std::vector< RegularGrid< 2 , Real > > inputCellConfidence;
	static std::vector< RegularGrid< 2 , Point3D< Real > > > inputTextures;
	static std::vector< std::vector< Point3D< Real > > > partialTexelValues;
	static std::vector< std::vector< Point3D< Real > > > partialEdgeValues;

	static RegularGrid< 2 , Point3D< Real > > filteredTexture;

#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
	// UI
	static char gradientFittingStr[1024];
	static char referenceTextureStr[1024];
#endif // NO_OPEN_GL_VISUALIZATION

	static std::vector< Point3D< float > >textureNodePositions;
	static std::vector< Point3D< float > >textureEdgePositions;

	static std::vector< TextureNodeInfo< PreReal > > textureNodes;

	static SparseMatrix< Real , int > mass;
	static SparseMatrix< Real , int > stiffness;
	static SparseMatrix< Real , int > stitchingMatrix;

	static std::vector< Point3D< Real > > texelMass;
	static std::vector< Point3D< Real > > texelDivergence;

	static std::vector< SystemCoefficients< Real > > multigridStitchingCoefficients;
	static std::vector< MultigridLevelVariables< Point3D< Real > > > multigridStitchingVariables;
	static std::vector< MultigridLevelIndices< Real > > multigridIndices;

#ifdef USE_EIGEN_PARDISO
	using EigenSolver = Eigen::PardisoLDLT< Eigen::SparseMatrix< double > >;
#else // !USE_EIGEN_PARDISO
	using EigenSolver = Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > >;
#endif // USE_EIGEN_PARDISO
	static VCycleSolvers< EigenSolver > vCycleSolvers;
	static EigenSolverWrapper< EigenSolver > directSolver;


	static DivergenceOperator< Real > divergenceOperator;

	static std::vector< bool > unobservedTexel;
	static std::vector< Point3D< Real > > texelValues;
	static std::vector< Point3D< Real > > edgeValues;

	// Linear Operators
	static MassAndStiffnessOperators< Real > massAndStiffnessOperators;
	static ScalarIntegrator< Real > scalarIntegrator;
	static GradientIntegrator< Real > gradientIntegrator;

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

	static void ToggleInputOutputCallBack             ( Visualization *v , const char *prompt );
	static void ToggleForwardReferenceTextureCallBack ( Visualization *v , const char *prompt );
	static void ToggleBackwardReferenceTextureCallBack( Visualization *v , const char *prompt );
	static void ToggleMaskCallBack                    ( Visualization *v , const char *prompt );
	static void ToggleUpdateCallBack                  ( Visualization *v , const char *prompt );
	static void IncrementUpdateCallBack               ( Visualization *v , const char *prompt );
	static void ExportTextureCallBack                 ( Visualization *v , const char *prompt );
	static void GradientFittingWeightCallBack         ( Visualization *v , const char *prompt );
#endif // NO_OPEN_GL_VISUALIZATION

	static RegularGrid< 2 , Point3D< unsigned char > > GetChartMask( void );
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
template< typename PreReal , typename Real , unsigned int TextureBitDepth > char															Stitching< PreReal , Real , TextureBitDepth >::gradientFittingStr[1024];
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

template< typename PreReal , typename Real , unsigned int TextureBitDepth > Real															Stitching< PreReal , Real , TextureBitDepth >::gradientFittingWeight;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< TextureNodeInfo< PreReal > >						Stitching< PreReal , Real , TextureBitDepth >::textureNodes;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > int																Stitching< PreReal , Real , TextureBitDepth >::steps;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > char															Stitching< PreReal , Real , TextureBitDepth >::stepsString[1024];
template< typename PreReal , typename Real , unsigned int TextureBitDepth > unsigned int													Stitching< PreReal , Real , TextureBitDepth >::levels;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > HierarchicalSystem< PreReal , Real >							Stitching< PreReal , Real , TextureBitDepth >::hierarchy;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > bool															Stitching< PreReal , Real , TextureBitDepth >::rhsNeedsUpdating = false;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > bool															Stitching< PreReal , Real , TextureBitDepth >::positiveModulation = true;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > RegularGrid< 2 , Point3D< Real > >								Stitching< PreReal , Real , TextureBitDepth >::filteredTexture;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > RegularGrid< 2 , int >										    Stitching< PreReal , Real , TextureBitDepth >::inputMask;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > RegularGrid< 2 , Point3D< Real > >								Stitching< PreReal , Real , TextureBitDepth >::lowFrequencyTexture;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > RegularGrid< 2 , Point3D< Real > >								Stitching< PreReal , Real , TextureBitDepth >::inputComposition;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > RegularGrid< 2 , Point3D< Real > >								Stitching< PreReal , Real , TextureBitDepth >::inputColorMask;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > int																Stitching< PreReal , Real , TextureBitDepth >::numTextures;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< RegularGrid< 2 , Real > >							Stitching< PreReal , Real , TextureBitDepth >::inputTexelConfidence;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< RegularGrid< 2 , Real > >							Stitching< PreReal , Real , TextureBitDepth >::inputCellConfidence;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< RegularGrid< 2 , Point3D< Real > > >				Stitching< PreReal , Real , TextureBitDepth >::inputTextures;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Point3D< Real > >									Stitching< PreReal , Real , TextureBitDepth >::texelMass;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Point3D< Real > >									Stitching< PreReal , Real , TextureBitDepth >::texelDivergence;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< SystemCoefficients< Real > >						Stitching< PreReal , Real , TextureBitDepth >::multigridStitchingCoefficients;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< MultigridLevelVariables< Point3D< Real > > >		Stitching< PreReal , Real , TextureBitDepth >::multigridStitchingVariables;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< MultigridLevelIndices< Real > >					Stitching< PreReal , Real , TextureBitDepth >::multigridIndices;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > VCycleSolvers< typename Stitching< PreReal , Real , TextureBitDepth >::EigenSolver >		Stitching< PreReal , Real , TextureBitDepth >::vCycleSolvers;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > EigenSolverWrapper< typename Stitching< PreReal , Real , TextureBitDepth >::EigenSolver >	Stitching< PreReal , Real , TextureBitDepth >::directSolver;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Point3D< float > >									Stitching< PreReal , Real , TextureBitDepth >::textureNodePositions;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > std::vector< Point3D< float > >									Stitching< PreReal , Real , TextureBitDepth >::textureEdgePositions;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > Padding															Stitching< PreReal , Real , TextureBitDepth >::padding;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > MassAndStiffnessOperators< Real >								Stitching< PreReal , Real , TextureBitDepth >::massAndStiffnessOperators;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > ScalarIntegrator< Real >										Stitching< PreReal , Real , TextureBitDepth >::scalarIntegrator;
template< typename PreReal , typename Real , unsigned int TextureBitDepth > GradientIntegrator< Real >										Stitching< PreReal , Real , TextureBitDepth >::gradientIntegrator;

template< typename PreReal , typename Real , unsigned int TextureBitDepth > unsigned int													Stitching< PreReal , Real , TextureBitDepth >::updateCount = static_cast< unsigned int >(-1);

template< typename PreReal , typename Real , unsigned int TextureBitDepth >  DivergenceOperator< Real >										Stitching< PreReal , Real , TextureBitDepth >::divergenceOperator;

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
			[&]( size_t i )
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
			[&]( size_t i )
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
	RegularGrid< 2 , Point3D< Real > > outputTexture = filteredTexture;
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
				ThreadPool::ParallelFor( 0 , textureNodePositions.size() , [&]( size_t i ){ if( Point3D< float >::SquareNorm( textureNodePositions[i]-selectedPoint )<radiusSquared ) texelValues[i] = partialTexelValues[ textureIndex ][i]; } );
				ThreadPool::ParallelFor( 0 , textureNodePositions.size() , [&]( size_t i ){ if( Point3D< float >::SquareNorm( textureEdgePositions[i]-selectedPoint )<radiusSquared )  edgeValues[i] =  partialEdgeValues[ textureIndex ][i]; } );
				rhsNeedsUpdating = true;
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

template< typename PreReal , typename Real , unsigned int TextureBitDepth > void Stitching< PreReal , Real , TextureBitDepth >::ToggleInputOutputCallBack( Visualization * /*v*/ , const char* )
{
	visualization.showSource = !visualization.showSource;
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
	RegularGrid< 2 , Point3D< Real > > outputTexture = filteredTexture;
	padding.unpad( outputTexture );
	WriteImage< TextureBitDepth >( outputTexture , prompt );
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void  Stitching< PreReal , Real , TextureBitDepth >::GradientFittingWeightCallBack( Visualization * /*v*/ , const char *prompt )
{
	Miscellany::PerformanceMeter pMeter( '.' );

	gradientFittingWeight = atof(prompt);
	if( UseDirectSolver.set )
	{
		stitchingMatrix = mass + stiffness * gradientFittingWeight;
		if( Verbose.set ) std::cout << pMeter( "Stitching matrix" ) << std::endl;
	}

	UpdateLinearSystem( static_cast< Real >( 1. ) , gradientFittingWeight , hierarchy , multigridStitchingCoefficients , massAndStiffnessOperators , vCycleSolvers , directSolver , stitchingMatrix , DetailVerbose.set , false , UseDirectSolver.set );
	if( Verbose.set ) std::cout << pMeter( "Initialized MG" ) << std::endl;

	ThreadPool::ParallelFor( 0 , multigridStitchingVariables[0].rhs.size() , [&]( size_t i ){ multigridStitchingVariables[0].rhs[i] = texelMass[i] + texelDivergence[i] * gradientFittingWeight; } );

	if( UseDirectSolver.set ) ComputeExactSolution( Verbose.set );
	else for( unsigned int i=0 ; i<OutputVCycles.value ; i++ ) VCycle( multigridStitchingVariables , multigridStitchingCoefficients , multigridIndices , vCycleSolvers , 2 , false , false );

	UpdateFilteredColorTexture( multigridStitchingVariables[0].x );
	visualization.UpdateColorTextureBuffer();
	sprintf( gradientFittingStr , "Gradient fitting weight: %e\n" , gradientFittingWeight );
}
#endif // NO_OPEN_GL_VISUALIZATION


template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::ComputeExactSolution( bool verbose )
{
	Miscellany::PerformanceMeter pMeter( '.' );
	directSolver.solve( multigridStitchingVariables[0].x , multigridStitchingVariables[0].rhs );
	if( verbose ) std::cout << pMeter( "Solved" ) << std::endl;
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::UpdateSolution( bool verbose , bool detailVerbose )
{
	Miscellany::PerformanceMeter pMeter( '.' );
	if( rhsNeedsUpdating )
	{
		if( inputMode==MULTIPLE_INPUT_MODE && inputCellConfidence.size() )
		{
			if( WeightType.value==static_cast< unsigned int >( WeightingType::Metric ) || WeightType.value==static_cast< unsigned int >( WeightingType::Confidence ) )
			{
				std::vector< Point3D< Real > > _texelMass( textureNodes.size() );
				for( unsigned int i=0 ; i<partialTexelValues.size() ; i++ )
				{
					scalarIntegrator( partialTexelValues[i] , []( Point3D< Real > v , SquareMatrix< Real , 2 >  ){ return v; } , _texelMass , [&]( typename RegularGrid< 2 >::Index I ){ return inputCellConfidence[i]( I ); } );
					for( unsigned int i=0 ; i<_texelMass.size() ; i++ ) texelMass[i] += _texelMass[i];
				}
			}
			else massAndStiffnessOperators.mass( texelValues , texelMass );
			if( WeightType.value==static_cast< unsigned int >( WeightingType::Confidence ) )
			{
				std::vector< Point3D< Real > > _texelDivergence( textureNodes.size() );
				for( unsigned int i=0 ; i<partialTexelValues.size() ; i++ )
				{
					gradientIntegrator( partialTexelValues[i] , []( Point2D< Point3D< Real > > v , SquareMatrix< Real , 2 > tensor ){ return v; } , _texelDivergence , [&]( typename RegularGrid< 2 >::Index I ){ return inputCellConfidence[i]( I ); } );
					for( unsigned int j=0 ; j<_texelDivergence.size() ; j++ ) texelDivergence[j] += _texelDivergence[j];
				}
			}
			else divergenceOperator( edgeValues , texelDivergence );
		}
		else
		{
			massAndStiffnessOperators.mass( texelValues , texelMass );
			divergenceOperator( edgeValues , texelDivergence );
		}

		ThreadPool::ParallelFor( 0 , textureNodes.size() , [&]( size_t i ){ multigridStitchingVariables[0].rhs[i] = texelMass[i] + texelDivergence[i] * gradientFittingWeight; } );

		if( verbose ) std::cout << pMeter( "RHS" ) << std::endl;
		rhsNeedsUpdating = false;
	}

	VCycle( multigridStitchingVariables , multigridStitchingCoefficients , multigridIndices , vCycleSolvers , 2 , verbose , detailVerbose );
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::InitializeSystem( int width , int height )
{
	Miscellany::PerformanceMeter pMeter( '.' );

	ExplicitIndexVector< ChartIndex , AtlasChart< PreReal > > atlasCharts;

	MultigridBlockInfo multigridBlockInfo( MultigridBlockWidth.value , MultigridBlockHeight.value , MultigridPaddedWidth.value , MultigridPaddedHeight.value );
	InitializeHierarchy( mesh , width , height , levels , textureNodes , hierarchy , atlasCharts , multigridBlockInfo , SanityCheck.set );

	if( Verbose.set ) std::cout << pMeter( "Hierarchy" ) << std::endl;

	ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< PreReal , 2 > > > parameterMetric;
	InitializeMetric( mesh , EMBEDDING_METRIC , atlasCharts , parameterMetric );

	pMeter.reset();
	if( inputMode==MULTIPLE_INPUT_MODE && inputCellConfidence.size() )
	{
		RegularGrid< 2 , Real > cellConfidence( inputCellConfidence[0].res() );
		for( unsigned int j=0 ; j<cellConfidence.size() ; j++ ) cellConfidence[j] = 0;
		for( unsigned int i=0 ; i<inputCellConfidence.size() ; i++ ) for( size_t j=0 ; j<cellConfidence.size() ; j++ ) cellConfidence[j] += inputCellConfidence[i][j];
		switch( WeightType.value )
		{
		case static_cast< unsigned int >( WeightingType::Metric ):
			OperatorInitializer::Initialize( MatrixQuadrature.value , massAndStiffnessOperators , hierarchy.gridAtlases[0] , parameterMetric , atlasCharts , divergenceOperator , std::tie( scalarIntegrator , gradientIntegrator ) , VectorQuadrature.value , false , SanityCheck.set , [&]( typename RegularGrid< 2 >::Index I ){ return cellConfidence(I); } );
			break;
		case static_cast< unsigned int >( WeightingType::Confidence ):
			OperatorInitializer::Initialize( MatrixQuadrature.value , massAndStiffnessOperators , hierarchy.gridAtlases[0] , parameterMetric , atlasCharts , divergenceOperator , std::tie( scalarIntegrator , gradientIntegrator ) , VectorQuadrature.value , false , SanityCheck.set , [&]( typename RegularGrid< 2 >::Index I ){ return cellConfidence(I); } , [&]( typename RegularGrid< 2 >::Index I ){ return cellConfidence(I) + ConfidenceRegularizationWeight.value; } );
			break;
		default:
			OperatorInitializer::Initialize( MatrixQuadrature.value , massAndStiffnessOperators , hierarchy.gridAtlases[0] , parameterMetric , atlasCharts , divergenceOperator , std::tie( scalarIntegrator , gradientIntegrator ) , VectorQuadrature.value , false , SanityCheck.set );
		}
	}
	else OperatorInitializer::Initialize( MatrixQuadrature.value , massAndStiffnessOperators , hierarchy.gridAtlases[0] , parameterMetric , atlasCharts , divergenceOperator , SanityCheck.set );
	if( Verbose.set ) std::cout << pMeter( "Mass and stiffness" ) << std::endl;

	texelMass.resize( textureNodes.size() );
	texelDivergence.resize( textureNodes.size() );

	if( UseDirectSolver.set )
	{
		FullMatrixConstruction( hierarchy.gridAtlases[0] , massAndStiffnessOperators.massCoefficients , mass );
		FullMatrixConstruction( hierarchy.gridAtlases[0] , massAndStiffnessOperators.stiffnessCoefficients , stiffness );
		stitchingMatrix = mass + stiffness * gradientFittingWeight;
	}

	multigridIndices.resize( levels );
	for( unsigned int i=0 ; i<levels ; i++ )
	{
		const typename GridAtlas<>::IndexConverter & indexConverter = hierarchy.gridAtlases[i].indexConverter;
		const GridAtlas< PreReal , Real > &gridAtlas = hierarchy.gridAtlases[i];
		multigridIndices[i].threadTasks = gridAtlas.threadTasks;
		multigridIndices[i].boundaryToCombined = indexConverter.boundaryToCombined();
		multigridIndices[i].segmentedLines = gridAtlas.segmentedLines;
		multigridIndices[i].rasterLines = gridAtlas.rasterLines;
		multigridIndices[i].restrictionLines = gridAtlas.restrictionLines;
		multigridIndices[i].prolongationLines = gridAtlas.prolongationLines;
		if( i<levels-1 ) multigridIndices[i].boundaryRestriction = hierarchy.boundaryRestriction[i];
	}

	pMeter.reset();
	UpdateLinearSystem( static_cast< Real >( 1. ) , gradientFittingWeight , hierarchy , multigridStitchingCoefficients , massAndStiffnessOperators , vCycleSolvers , directSolver , stitchingMatrix , DetailVerbose.set , true , UseDirectSolver.set );
	if( Verbose.set ) std::cout << pMeter( "Initialize MG" ) << std::endl;

	multigridStitchingVariables.resize( levels );
	for( unsigned int i=0 ; i<levels ; i++ )
	{
		const typename GridAtlas<>::IndexConverter & indexConverter = hierarchy.gridAtlases[i].indexConverter;
		MultigridLevelVariables< Point3D< Real > >& variables = multigridStitchingVariables[i];
		variables.x.resize( indexConverter.numCombined() );
		variables.rhs.resize( indexConverter.numCombined() );
		variables.residual.resize( indexConverter.numCombined() );
		variables.boundary_rhs.resize( indexConverter.numBoundary() );
		variables.boundary_value.resize( indexConverter.numBoundary() );
		variables.variable_boundary_value.resize( indexConverter.numBoundary() );
	}
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::SetUpSystem( void )
{
	texelMass.resize( textureNodes.size() );
	texelDivergence.resize( textureNodes.size() );
	for( unsigned int i=0 ; i<texelMass.size() ; i++ ) texelMass[i] = texelDivergence[i] = Point3D< Real >();

	if( inputMode==MULTIPLE_INPUT_MODE && inputCellConfidence.size() )
	{
		if( WeightType.value==static_cast< unsigned int >( WeightingType::Metric ) || WeightType.value==static_cast< unsigned int >( WeightingType::Confidence ) )
		{
			std::vector< Point3D< Real > > _texelMass( textureNodes.size() );
			for( unsigned int i=0 ; i<partialTexelValues.size() ; i++ )
			{
				scalarIntegrator( partialTexelValues[i] , []( Point3D< Real > v , SquareMatrix< Real , 2 >  ){ return v; } , _texelMass , [&]( typename RegularGrid< 2 >::Index I ){ return inputCellConfidence[i]( I ); } );
				for( unsigned int j=0 ; j<_texelMass.size() ; j++ ) texelMass[j] += _texelMass[j];
			}
		}
		else massAndStiffnessOperators.mass( texelValues , texelMass );

		if( WeightType.value==static_cast< unsigned int >( WeightingType::Confidence ) )
		{
			std::vector< Point3D< Real > > _texelDivergence( textureNodes.size() );
			for( unsigned int i=0 ; i<partialTexelValues.size() ; i++ )
			{
				gradientIntegrator( partialTexelValues[i] , []( Point2D< Point3D< Real > > v , SquareMatrix< Real , 2 > tensor ){ return v; } , _texelDivergence , [&]( typename RegularGrid< 2 >::Index I ){ return inputCellConfidence[i]( I ); } );
				for( unsigned int j=0 ; j<_texelDivergence.size() ; j++ ) texelDivergence[j] += _texelDivergence[j];
			}
		}
		else divergenceOperator( edgeValues , texelDivergence );
	}
	else
	{
		massAndStiffnessOperators.mass( texelValues , texelMass );
		divergenceOperator( edgeValues , texelDivergence );
	}

	ThreadPool::ParallelFor
		(
			0 , texelValues.size() ,
			[&]( size_t i )
			{
				multigridStitchingVariables[0].x[i] = texelValues[i];
				multigridStitchingVariables[0].rhs[i] = texelMass[i] + texelDivergence[i] * gradientFittingWeight;
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
	else for( unsigned int i=0 ; i<OutputVCycles.value ; i++ ) VCycle( multigridStitchingVariables , multigridStitchingCoefficients , multigridIndices , vCycleSolvers , 2 , false , false );
	UpdateFilteredTexture( multigridStitchingVariables[0].x );
}

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
RegularGrid< 2 , Point3D< unsigned char > > Stitching< PreReal , Real , TextureBitDepth >::GetChartMask( void )
{
	auto IsBlack = []( Point3D< unsigned char > c ){ return !c[0] && !c[1] && !c[2]; };

	RegularGrid< 2 , Point3D< unsigned char > > chartMask;
	chartMask.resize( textureWidth , textureHeight );
	for( unsigned int i=0 ; i<(unsigned int)textureWidth ; i++ ) for( unsigned int j=0 ; j<(unsigned int)textureHeight ; j++ ) chartMask(i,j) = Point3D< unsigned char >(0,0,0);

	std::map< ChartIndex , Point3D< unsigned char > > chartColors;
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
	for( unsigned int e=0 ; e<ChartMaskErode.value ; e++ )
	{
		RegularGrid< 2 , Point3D< unsigned char > > _chartMask = chartMask;
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
				[&]( size_t i )
				{
					char textureName[256];
					sprintf( textureName , In.values[1].c_str() , i );
					ReadImage< TextureBitDepth >( inputTextures[i] , textureName );
				}
			);

		textureWidth = inputTextures[0].res(0);
		textureHeight = inputTextures[0].res(1);
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
		RegularGrid< 2 , unsigned char > mask;
		if( ChartMaskErode.value )
		{
			static const bool Nearest = false;
			static const bool NodeAtCellCenter = true;
			using Index = unsigned int;
			using TexelInfo = typename Texels< NodeAtCellCenter , Index >::template TexelInfo< 2 >;

			mask.resize( textureWidth , textureHeight );
			for( size_t i=0 ; i<mask.size() ; i++ ) mask[i] = 0;
			TexturedTriangleMesh< PreReal > mesh;
			mesh.read( In.values[0] , false , CollapseEpsilon.value );
			auto TextureSimplexFunctor = [&]( size_t sIdx )
			{
				Simplex< double , 2 , 2 > simplex;
				for( unsigned int k=0 ; k<=2 ; k++ ) simplex[k] = mesh.texture.vertices[ mesh.texture.triangles[sIdx][k] ];
				return simplex;
			};

			RegularGrid< 2 , TexelInfo > activeTexels = Texels< NodeAtCellCenter , Index >::template GetSupportedTexelInfo< Nearest >( mesh.numTriangles() , TextureSimplexFunctor , mask.res() , 0 , false );
			for( size_t i=0 ; i<activeTexels.size() ; i++ ) if( activeTexels[i].sIdx!=-1 ) mask[i] = 1;
		}
		inputTexelConfidence.resize( numTextures );
		inputCellConfidence.resize( numTextures );
		ThreadPool::ParallelFor
			(
				0 , numTextures ,
				[&]( size_t i )
				{
					char confidenceName[256];
					sprintf( confidenceName , InMask.value.c_str() , i );
					RegularGrid< 2 , Point3D< Real > > confidence;
					ReadImage< 8 >( confidence , confidenceName );
					inputTexelConfidence[i].resize( textureWidth , textureHeight );
					for( size_t p=0 ; p<confidence.size() ; p++ )
						inputTexelConfidence[i][p] = Point3D< Real >::Dot( confidence[p] , Point3D< Real >( (Real)1./3 , (Real)1./3 , (Real)1./3 ) );
					if( ChartMaskErode.value )
					{
						auto Erode = [&mask]( RegularGrid< 2 , Real > & conf , int radius , Real epsilon )
						{
							unsigned int radius2 = radius * radius;
							RegularGrid< 2 , unsigned int > d2 = SquareDistanceToTextureBoundary( conf.res(0) , conf.res(1) , [&]( typename RegularGrid< 2 >::Index I ){ return conf( I )>epsilon || !mask( I ); } );
							if( radius>0 ) ThreadPool::ParallelFor( 0 , conf.size() , [&]( size_t i ){ if( d2[i]<radius2 ) conf[i] = 0; } );
							else           ThreadPool::ParallelFor( 0 , conf.size() , [&]( size_t i ){ if( d2[i]<radius2 ) conf[i] *= static_cast< Real >( sqrt( static_cast< double >( d2[i] )/radius2 ) ); } );
						};
						Erode( inputTexelConfidence[i] , ChartMaskErode.value , ChartMaskErosionEps.value );
					}
					for( size_t p=0 ; p<confidence.size() ; p++ ) inputTexelConfidence[i][p] = std::max< Real >( inputTexelConfidence[i][p] , MinimumConfidence.value );
				}
			);
		if( ConfidenceFilter.set )
		{
			try
			{
				EquationParser::Node filter( ConfidenceFilter.value , { "s" } );
				ThreadPool::ParallelFor
				(
					0 , inputTexelConfidence.size() ,
					[&]( size_t i )
					{
						for( size_t j=0 ; j<inputTexelConfidence[i].size() ; j++ )
						{
							double value = std::max< Real >( 0 , std::min< Real >( 1 , inputTexelConfidence[i][j] ) );
							inputTexelConfidence[i][j] = filter( &value );
						}
					}
				);
			}
			catch( const Exception & exception )
			{
				std::cout << exception.what() << std::endl;
				MK_WARN( "Could not parse confidence filter: " , ConfidenceFilter.value );
			}
		}

#ifdef NORMALIZE_CONFIDENCE
		// Normalizing so that the sum of probabilities is in [0,1]
		ThreadPool::ParallelFor
		(
			0 , inputTexelConfidence[0].size() ,
			[&]( size_t j )
			{
				Real sum = 0 , prob = 1;
				for( unsigned int i=0 ; i<inputTexelConfidence.size() ; i++ )
				{
					inputTexelConfidence[i][j] = std::max< Real >( 0 , std::min< Real >( 1 , inputTexelConfidence[i][j] ) );
					prob *= ( 1 - inputTexelConfidence[i][j] );
					sum += inputTexelConfidence[i][j];
				}
				prob = 1-prob;
				if( sum>0 && prob>0 ) for( unsigned int i=0 ; i<inputTexelConfidence.size() ; i++ ) inputTexelConfidence[i][j] *= prob / sum;
				else                  for( unsigned int i=0 ; i<inputTexelConfidence.size() ; i++ ) inputTexelConfidence[i][j] *= 0;
			}
		);
#endif // NORMALIZE_CONFIDENCE
		// Turn per-texel confidences into per-cell confidences by computing the geometric mean
		ThreadPool::ParallelFor
			(
				0 , numTextures ,
				[&]( size_t i )
				{
					inputCellConfidence[i].resize( textureWidth-1 , textureHeight-1 );
					for( unsigned int x=0 ; x<static_cast< unsigned int >( textureWidth-1 ) ; x++ ) for( unsigned int y=0 ; y<static_cast< unsigned int >( textureHeight-1 ) ; y++ )
						inputCellConfidence[i](x,y) = pow( fabs( inputTexelConfidence[i](x,y) * inputTexelConfidence[i](x+1,y) * inputTexelConfidence[i](x+1,y+1) * inputTexelConfidence[i](x,y+1) ) , 0.25 );
				}
			);
	}
	else
	{
		RegularGrid< 2 , Point3D< unsigned char > > textureConfidence;
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
	unsigned int numNodes = (unsigned int)textureNodes.size();
	unsigned int numEdges = static_cast< unsigned int >( divergenceOperator.edges.size() );
	unobservedTexel.resize( numNodes );
	texelValues.resize( numNodes );
 	edgeValues.resize( numEdges );

	if( inputMode==MULTIPLE_INPUT_MODE )
	{
		std::vector< Real > texelWeight( numNodes , 0 );
#ifdef NORMALIZE_CONFIDENCE
#else // !NORMALIZE_CONFIDENCE
		std::vector< Real > edgeWeight( numEdges , 0 );
#endif // NORMALIZE_CONFIDENCE

		partialTexelValues.resize( numTextures );
		for( int i=0 ; i<numTextures ; i++ ) partialTexelValues[i].resize( numNodes );
		partialEdgeValues.resize( numTextures );
		for( int i=0 ; i<numTextures ; i++ ) partialEdgeValues[i].resize( divergenceOperator.edges.size() );

		for( unsigned int i=0 ; i<textureNodes.size() ; i++ ) for( int textureIter=0 ; textureIter<numTextures ; textureIter++ )
		{
			const RegularGrid< 2 , Point3D< Real > > & textureValues = InputLowFrequency.set ? lowFrequencyTexture : inputTextures[textureIter];
			const RegularGrid< 2 , Real > & textureConfidence = inputTexelConfidence[textureIter];

			Real weight = textureConfidence( textureNodes[i].ci , textureNodes[i].cj );
			Point3D< Real > value = textureValues( textureNodes[i].ci , textureNodes[i].cj );
			partialTexelValues[textureIter][i] = value;
			texelValues[i] += value * weight;
			texelWeight[i] += weight;
		}

		for( unsigned int e=0 ; e<divergenceOperator.edges.size() ; e++ )
		{
#ifdef NORMALIZE_CONFIDENCE
			for( int textureIter=0 ; textureIter<numTextures ; textureIter++ )
			{
				const RegularGrid< 2 , Point3D< Real > > & textureValues = InputLowFrequency.set ? lowFrequencyTexture : inputTextures[textureIter];
				const RegularGrid< 2 , Real > & textureConfidence = inputTexelConfidence[textureIter];

				Real weight[2];
				Point3D< Real > value[2];
				for( unsigned int k=0 ; k<2 ; k++ )
				{
					AtlasTexelIndex i = divergenceOperator.edges[e][k];
					weight[k] = textureConfidence( textureNodes[ static_cast< unsigned int >(i) ].ci , textureNodes[ static_cast< unsigned int >(i) ].cj );
					value[k] = textureValues( textureNodes[ static_cast< unsigned int >(i) ].ci , textureNodes[ static_cast< unsigned int >(i) ].cj );
				}
				Real eWeight = sqrt( weight[0] * weight[1] );
				Point3D< Real > eValue = value[1] - value[0];
				partialEdgeValues[textureIter][e] = eWeight > 0 ? eValue : Point3D< Real >();

				edgeValues[e] += eValue * eWeight;
			}
#else // !NORMALIZE_CONFIDENCE
			Real eWeight = 0 , maxEWeight = 0;
			Point3D< Real > eValue;

			for( int textureIter=0 ; textureIter<numTextures ; textureIter++ )
			{
				const RegularGrid< 2 , Point3D< Real > > & textureValues = InputLowFrequency.set ? lowFrequencyTexture : inputTextures[textureIter];
				const RegularGrid< 2 , Real > & textureConfidence = inputTexelConfidence[textureIter];

				Real weight[2];
				Point3D< Real > value[2];
				for( unsigned int k=0 ; k<2 ; k++ )
				{
					AtlasTexelIndex i = divergenceOperator.edges[e][k];
					weight[k] = textureConfidence( textureNodes[ static_cast< unsigned int >(i) ].ci , textureNodes[ static_cast< unsigned int >(i) ].cj );
					value[k] = textureValues( textureNodes[ static_cast< unsigned int >(i) ].ci , textureNodes[ static_cast< unsigned int >(i) ].cj );
				}
				Real _eWeight = weight[0] * weight[1];
				Point3D< Real > _eValue = value[1] - value[0];
				partialEdgeValues[textureIter][e] = _eWeight > 0 ? _eValue : Point3D< Real >();

				maxEWeight = std::max< Real >( maxEWeight , _eWeight );
				eValue += _eValue * _eWeight;
				eWeight += _eWeight;
			}
			edgeValues[e] = eValue;
			edgeWeight[e] = eWeight;
#endif // NORMALIZE_CONFIDENCE
		}

		for( unsigned int i=0 ; i<numNodes ; i++ )
		{
			unobservedTexel[i] = texelWeight[i]==0;
			if( texelWeight[i]>0 ) texelValues[i] /= texelWeight[i];
		}

#ifdef NORMALIZE_CONFIDENCE
#else // !NORMALIZE_CONFIDENCE
		if( WeightType.value==static_cast< unsigned int >( WeightingType::GradientScale ) ) for( unsigned int i=0 ; i<numEdges ; i++ ) if( edgeWeight[i]>0 ) edgeValues[i] /= edgeWeight[i];
#endif // NORMALIZE_CONFIDENCE
	}
	else
	{
		for( int i=0 ; i<textureNodes.size() ; i++ )
		{
			unobservedTexel[i] = inputMask( textureNodes[i].ci , textureNodes[i].cj )==-1;
			texelValues[i] = InputLowFrequency.set ? lowFrequencyTexture( textureNodes[i].ci , textureNodes[i].cj ) : inputComposition( textureNodes[i].ci , textureNodes[i].cj );
		}
		for( unsigned int e=0 ; e<divergenceOperator.edges.size() ; e++ )
		{
			const SimplexIndex< 1 , AtlasTexelIndex > &edgeCorners = divergenceOperator.edges[e];
			unsigned int ci[] = { textureNodes[ static_cast< unsigned int >(edgeCorners[0]) ].ci , textureNodes[ static_cast< unsigned int >(edgeCorners[1]) ].ci };
			unsigned int cj[] = { textureNodes[ static_cast< unsigned int >(edgeCorners[0]) ].cj , textureNodes[ static_cast< unsigned int >(edgeCorners[1]) ].cj };
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

	ThreadPool::ParallelFor( 0 , textureNodes.size() , [&]( size_t i ){ multigridStitchingVariables[0].x[i] = texelValues[i]; } );
}

#ifdef NO_OPEN_GL_VISUALIZATION
#else // !NO_OPEN_GL_VISUALIZATION
template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::InitializeVisualization( void )
{
	sprintf( referenceTextureStr , "Reference Texture: %02d of %02d\n" , textureIndex,numTextures );
	sprintf( gradientFittingStr , "Gradient fitting weight: %.2e\n" , gradientFittingWeight );

	visualization.textureWidth = textureWidth;
	visualization.textureHeight = textureHeight;

	visualization.colorTextureBuffer = new unsigned char[textureHeight*textureWidth * 3];
	memset( visualization.colorTextureBuffer , 204 , textureHeight * textureWidth * 3 * sizeof(unsigned char) );

	unsigned int tCount = (unsigned int)mesh.numTriangles();

	visualization.triangles.resize( tCount );
	visualization.vertices.resize( 3*tCount );
	visualization.colors.resize( 3*tCount , Point3D< float >( 0.75f , 0.75f , 0.75f ) );
	visualization.textureCoordinates.resize( 3*tCount );
	visualization.normals.resize( 3*tCount );


	for( unsigned int t=0 , idx=0 ; t<tCount ; t++ )
	{
		Point3D< float > n = static_cast< Point3D< float > >( mesh.surfaceTriangle( t ).normal() );
		n /= Point3D< float >::Length( n );

		for( int k=0 ; k<3 ; k++ , idx++ )
		{
			visualization.triangles[t][k] = idx;
			visualization.vertices[idx] = static_cast< Point3D< float > >( mesh.surface.vertices[ mesh.surface.triangles[t][k] ] );
			visualization.normals[idx] = n;
			visualization.textureCoordinates[idx] = static_cast< Point2D< float > >( mesh.texture.vertices[ mesh.texture.triangles[t][k] ] );
		}
	}

	std::vector< unsigned int > boundaryHalfEdges = mesh.texture.boundaryHalfEdges();

	for( int e=0 ; e<boundaryHalfEdges.size() ; e++ )
	{
		SimplexIndex< 1 > eIndex = mesh.surface.edgeIndex( boundaryHalfEdges[e] );
		for( int i=0 ; i<2 ; i++ ) visualization.chartBoundaryVertices.push_back( static_cast< Point3D< float > >( mesh.surface.vertices[ eIndex[i] ] ) );
	}

	if( inputMode==SINGLE_INPUT_MODE ) visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'M' , "toggle mask" , ToggleMaskCallBack ) );
	else                               visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'M' , "toggle weights" , ToggleMaskCallBack ) );
	if( SingleView.set ) visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'S' , "toggle input/output" , ToggleInputOutputCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 's' , "export texture" , "Output Texture" , ExportTextureCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'y' , "gradient fitting weight" , "Gradient Fitting Weight" , GradientFittingWeightCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , ' ' , "toggle update" , ToggleUpdateCallBack ) );
	visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , '+' , "increment update" , IncrementUpdateCallBack ) );
	
	if( inputMode==MULTIPLE_INPUT_MODE )
	{
		visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 't' , "toggle reference" , ToggleForwardReferenceTextureCallBack ) );
		visualization.callBacks.push_back( Visualization::KeyboardCallBack( &visualization , 'T' , "toggle reference" , ToggleBackwardReferenceTextureCallBack ) );
	}
	
	visualization.info.push_back( stepsString );
	
	visualization.info.push_back( gradientFittingStr );
	if( inputMode==MULTIPLE_INPUT_MODE ) visualization.info.push_back( referenceTextureStr );

	visualization.UpdateVertexBuffer();
	visualization.UpdateFaceBuffer();
	visualization.UpdateTextureBuffer( filteredTexture );

	if( inputMode==MULTIPLE_INPUT_MODE ) visualization.UpdateReferenceTextureBuffers( inputTextures );
	if( inputMode==MULTIPLE_INPUT_MODE ) visualization.UpdateReferenceConfidenceBuffers( inputTexelConfidence );
	if( inputMode==SINGLE_INPUT_MODE )   visualization.UpdateCompositeTextureBuffer( inputComposition );
	if( inputMode==SINGLE_INPUT_MODE )   visualization.UpdateMaskTextureBuffer( inputColorMask );

}
#endif // NO_OPEN_GL_VISUALIZATION

template< typename PreReal , typename Real , unsigned int TextureBitDepth >
void Stitching< PreReal , Real , TextureBitDepth >::Init( void )
{
	sprintf( stepsString , "Steps: 0" );
	levels = std::max< unsigned int >( Levels.value , 1 );
	gradientFittingWeight = GradientFittingWeight.value;

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
		if( inputMode==MULTIPLE_INPUT_MODE ) for( int i=0 ; i<numTextures ; i++ ) padding.pad( inputTextures[i] ) , padding.pad( inputTexelConfidence[i] );
		else padding.pad( inputComposition ) , padding.pad( inputMask );
		textureWidth  += padding.width();
		textureHeight += padding.height();
	}

	Miscellany::PerformanceMeter pMeter( '.' );
	InitializeSystem( textureWidth , textureHeight );
	if( Verbose.set ) std::cout << pMeter( "Initialized" ) << std::endl;
	if( Verbose.set ) printf( "Resolution: %d / %d x %d\n" , (int)textureNodes.size() , textureWidth , textureHeight );

	// Assign position to exterior nodes using barycentric-exponential map
	{
		FEM::RiemannianMesh< PreReal , unsigned int > rMesh( GetPointer( mesh.surface.triangles) , mesh.surface.triangles.size() );
		rMesh.setMetricFromEmbedding( GetPointer( mesh.surface.vertices ) );
		rMesh.makeUnitArea();
		Pointer( FEM::CoordinateXForm< PreReal > ) xForms = rMesh.getCoordinateXForms();

		for( unsigned int i=0 ; i<textureNodes.size() ; i++ ) if( textureNodes[i].tID!=AtlasMeshTriangleIndex(-1) && !textureNodes[i].isInterior )
		{
			FEM::HermiteSamplePoint< PreReal > _p;
			_p.tIdx = static_cast< unsigned int >( textureNodes[i].tID );
			_p.p = Point2D< PreReal >( (PreReal)1./3 , (PreReal)1./3 );
			_p.v = textureNodes[i].barycentricCoords - _p.p;

			if( SanityCheck.set ) rMesh.exp( xForms , _p , 0 , false );
			else rMesh.exp( xForms , _p );

			textureNodes[i].tID = AtlasMeshTriangleIndex( _p.tIdx );
			textureNodes[i].barycentricCoords = _p.p;
		}
	}

	textureNodePositions.resize( textureNodes.size() );
	for( int i=0 ; i<textureNodePositions.size() ; i++ ) textureNodePositions[i] = static_cast< Point3D< float > >( mesh.surface( textureNodes[i] ) );

	textureEdgePositions.resize( divergenceOperator.edges.size() );
	for( int i=0 ; i<divergenceOperator.edges.size() ; i++ ) textureEdgePositions[i] = ( textureNodePositions[ static_cast< unsigned int >(divergenceOperator.edges[i][0]) ] + textureNodePositions[ static_cast< unsigned int >(divergenceOperator.edges[i][1]) ] ) / 2;

	{
		unsigned int multiChartTexelCount = 0;
		RegularGrid< 2 , int > texelId;
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
	Stitching< PreReal , Real , TextureBitDepth >::updateCount = Run.set ? static_cast< unsigned int >(-1) : 0;

	Stitching< PreReal , Real , TextureBitDepth >::LoadTextures();
	Stitching< PreReal , Real , TextureBitDepth >::LoadMasks();
	Stitching< PreReal , Real , TextureBitDepth >::Init();
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
		Stitching< PreReal , Real , TextureBitDepth >::visualization.displayMode = SingleView.set ? ONE_REGION_DISPLAY : TWO_REGION_DISPLAY;
		Stitching< PreReal , Real , TextureBitDepth >::visualization.screenWidth = SingleView.set ? 800 : 1600;
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
#ifdef NORMALIZE_CONFIDENCE
#else // !NORMALIZE_CONFIDENCE
	if( MultiInput.set && WeightType.value==static_cast< unsigned int >( WeightingType::Confidence ) && MinimumConfidence.value<=0 )
		MK_WARN( "Likely need a positive minimum confidence for confidence weighting" );
#endif // NORMALIZE_CONFIDENCE

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
					GetImageInfo( textureName , _width , _height , _channels , _bitDepth );
					if( !numTextures ) width = _width , height = _height , channels = _channels , bitDepth = _bitDepth;
					else if( width!=_width || height!=_height || channels!=_channels || bitDepth!=_bitDepth )
						MK_THROW( "Image properties don't match: (" , width , " " , height , " " , channels , " " , bitDepth , ") != (" , _width , " " , _height , " " , _channels , " " , _bitDepth , ")" );
					numTextures++;
				}
				else break;
			}
		}
		else GetImageInfo( In.values[1].c_str() , width , height , channels , bitDepth );
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
