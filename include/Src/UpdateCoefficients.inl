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

#include <Misha/Miscellany.h>


template< class Real >
void DeepCoefficientRestriction( const std::vector< Real >& fineDeepCoefficients , std::vector< Real >& coarseDeepCoefficients , const std::vector< DeepLine > & deepLines )
{
	auto UpdateRow = [&]( int r )
	{
		const Real * prevLineCoeff = fineDeepCoefficients.data() + deepLines[r].finePrevLineIndex * 10;
		const Real * currentLineCoeff = fineDeepCoefficients.data() + deepLines[r].fineCurrentLineIndex * 10;
		const Real * nextLineCoeff = fineDeepCoefficients.data() + deepLines[r].fineNextLineIndex * 10;
		Real * coarseLineCoeff = coarseDeepCoefficients.data() + deepLines[r].coarseLineStartIndex * 10;

		int lineLength = deepLines[r].coarseLineEndIndex - deepLines[r].coarseLineStartIndex + 1;

		auto UpdateNode = [&]( int i )
		{
			const Real * prevCoeff = prevLineCoeff + 20 * i;
			const Real * currentCoeff = currentLineCoeff + 20 * i;
			const Real * nextCoeff = nextLineCoeff + 20 * i;
			Real * coarseCoeff = coarseLineCoeff + 10 * i;
			Real cumValues[] =
			{
				0,0,0,
				0,0,0,
				0,0,0
			};

			static const Real HALF = (Real)0.5;
			static const Real QUARTER = (Real)0.25;
			//Procces middle cell node
			for( int k=0 ; k<2 ; k++ ) for( int l=0 ; l<2 ; l++ )
			{
				const Real * nodeCoeff = (l == 0) ? prevCoeff + 10 * (2 * k - 1) : nextCoeff + 10 * (2 * k - 1);
				cumValues[3 * l + k] += (nodeCoeff[0] + nodeCoeff[1] * HALF + nodeCoeff[3] * HALF + nodeCoeff[4] * QUARTER ) * QUARTER;
				cumValues[3 * l + k + 1] += (nodeCoeff[1] * HALF + nodeCoeff[2] + nodeCoeff[4] * QUARTER + nodeCoeff[5] * HALF) * QUARTER;
				cumValues[3 * (l + 1) + k] += (nodeCoeff[3] * HALF + nodeCoeff[4] * QUARTER + nodeCoeff[6] + nodeCoeff[7] * HALF) * QUARTER;
				cumValues[3 * (l + 1) + k + 1] += (nodeCoeff[4] * QUARTER + nodeCoeff[5] * HALF + nodeCoeff[7] * HALF + nodeCoeff[8]) * QUARTER;
			}
			//Procces vertical edge node
			for (int l = 0; l < 2; l++) {
				const Real * nodeCoeff = (l == 0) ? prevCoeff : nextCoeff;
				cumValues[3 * l] += (nodeCoeff[0] * HALF + nodeCoeff[3] * QUARTER) * HALF;
				cumValues[3 * l + 1] += (nodeCoeff[0] * HALF + nodeCoeff[1] + nodeCoeff[2] * HALF + nodeCoeff[3] * QUARTER + nodeCoeff[4] * HALF + nodeCoeff[5] * QUARTER) * HALF;
				cumValues[3 * l + 2] += (nodeCoeff[2] * HALF + nodeCoeff[5] * QUARTER) * HALF;

				cumValues[3 * (l + 1)] += (nodeCoeff[3] * QUARTER + nodeCoeff[6] * HALF) * HALF;
				cumValues[3 * (l + 1) + 1] += (nodeCoeff[3] * QUARTER + nodeCoeff[4] * HALF + nodeCoeff[5] * QUARTER + nodeCoeff[6] * HALF + nodeCoeff[7] + nodeCoeff[8] * HALF)*HALF;
				cumValues[3 * (l + 1) + 2] += (nodeCoeff[5] * QUARTER + nodeCoeff[8] * HALF) * HALF;
			}
			//Procces horizontal edge node
			for (int k = 0; k < 2; k++) {
				const Real * nodeCoeff = (k == 0) ? currentCoeff - 10 : currentCoeff + 10;
				cumValues[k] += (nodeCoeff[0] * HALF + nodeCoeff[1] * QUARTER)*HALF;
				cumValues[3 + k] += (nodeCoeff[0] * HALF + nodeCoeff[1] * QUARTER + nodeCoeff[3] + nodeCoeff[4] * HALF + nodeCoeff[6] * HALF + nodeCoeff[7] * QUARTER)*HALF;
				cumValues[6 + k] += (nodeCoeff[6] * HALF + nodeCoeff[7] * QUARTER)*HALF;

				cumValues[k + 1] += (nodeCoeff[1] * QUARTER + nodeCoeff[2] * HALF) * HALF;
				cumValues[3 + k + 1] += (nodeCoeff[1] * QUARTER + nodeCoeff[2] * HALF + nodeCoeff[4] * HALF + nodeCoeff[5] + nodeCoeff[7] * QUARTER + nodeCoeff[8] * HALF)*HALF;
				cumValues[6 + k + 1] += (nodeCoeff[7] * QUARTER + nodeCoeff[8] * HALF) * HALF;
			}
			//Procces center node
			{
				const Real * nodeCoeff = currentCoeff;
				cumValues[0] += (nodeCoeff[0] * QUARTER);
				cumValues[1] += (nodeCoeff[0] * QUARTER + nodeCoeff[1] * HALF + nodeCoeff[2] * QUARTER);
				cumValues[2] += (nodeCoeff[2] * QUARTER);
				cumValues[3] += (nodeCoeff[0] * QUARTER + nodeCoeff[3] * HALF + nodeCoeff[6] * QUARTER);

				cumValues[4] += (nodeCoeff[0] * QUARTER + nodeCoeff[1] * HALF + nodeCoeff[2] * QUARTER + nodeCoeff[3] * HALF + nodeCoeff[4] + nodeCoeff[5] * HALF +
					nodeCoeff[6] * QUARTER + nodeCoeff[7] * HALF + nodeCoeff[8] * QUARTER);

				cumValues[5] += (nodeCoeff[2] * QUARTER + nodeCoeff[5] * HALF + nodeCoeff[8] * QUARTER);
				cumValues[6] += (nodeCoeff[6] * QUARTER);
				cumValues[7] += (nodeCoeff[6] * QUARTER + nodeCoeff[7] * HALF + nodeCoeff[8] * QUARTER);
				cumValues[8] += nodeCoeff[8] * QUARTER;
			}
			for( int k=0 ; k<3 ; k++ ) for( int l=0 ; l<3 ; l++ ) coarseCoeff[3*l+k] = cumValues[3*l+k];
		};
		for (int i = 0; i < lineLength;i++)UpdateNode(i);
	};
	int threads = omp_get_max_threads();
#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		int tID = omp_get_thread_num();
		int firstLine = ((int)deepLines.size() * tID) / threads;
		int lastLine = ((int)deepLines.size() * (tID + 1)) / threads;
		for (int r = firstLine; r < lastLine; r++) UpdateRow(r);
	}
}

template< typename MatrixReal >
void BoundaryDeepMatrixConstruction( int numBoundayTexels , int numTexels , const std::vector< MatrixReal > &deepCoefficients , const std::vector< BoundaryDeepIndex > &boundaryDeepIndices , SparseMatrix< MatrixReal , int > &boundaryDeepMatrix )
{
	std::vector< Eigen::Triplet< MatrixReal > > boundaryDeepTriplets;
	for( int i=0 ; i<boundaryDeepIndices.size() ; i++ )
	{
		boundaryDeepTriplets.push_back( Eigen::Triplet< MatrixReal >( boundaryDeepIndices[i].boundaryIndex , boundaryDeepIndices[i].deepGlobalIndex , deepCoefficients[ 10*boundaryDeepIndices[i].deepIndex + boundaryDeepIndices[i].offset ] ) );
		//NOTE: Deep coefficient is not multiplied by reciprocal yet.
	}
	boundaryDeepMatrix = SetSparseMatrix( boundaryDeepTriplets , numBoundayTexels , numTexels , false );
}

template< typename MatrixReal >
void BoundaryBoundaryMatrixConstruction( int numBoundayTexels , const std::vector< MatrixReal > &fineDeepCoefficients , const SparseMatrix< MatrixReal , int > &fineBoundaryBoundaryMatrix , const SparseMatrix< MatrixReal , int > &boundaryCoarseFineProlongation , const SparseMatrix< MatrixReal , int > &boundaryFineCoarseRestriction , const std::vector< BoundaryBoundaryIndex< MatrixReal > > &boundaryBoundaryIndices , SparseMatrix< MatrixReal , int > &coarseBoundaryBoundaryMatrix , bool verbose=false )
{
	Miscellany::Timer timer;

	SparseMatrix< MatrixReal , int > partialRestriction = fineBoundaryBoundaryMatrix * boundaryCoarseFineProlongation;
	partialRestriction = boundaryFineCoarseRestriction *partialRestriction;

	if( verbose ) printf( "\t Multiplication  =  %.4f\n" , timer.elapsed() );

	if( verbose ) timer.reset();
	std::vector< Eigen::Triplet< MatrixReal > > boundaryBoundaryTriplets;
	for( int i=0 ; i<boundaryBoundaryIndices.size() ; i++ )
	{
		boundaryBoundaryTriplets.push_back( Eigen::Triplet< MatrixReal >( boundaryBoundaryIndices[i].coarsePrincipalBoundaryIndex , boundaryBoundaryIndices[i].coarseSecondaryBoundaryIndex , fineDeepCoefficients[ boundaryBoundaryIndices[i].fineDeepIndex*10 + boundaryBoundaryIndices[i].offset ] * boundaryBoundaryIndices[i].weight ) );
		//NOTE: Deep coefficient is not multiplied by reciprocal yet.
	}
	if( verbose ) printf( "\t Triplets  =  %.4f\n" , timer.elapsed() );

	if( verbose ) timer.reset();
	coarseBoundaryBoundaryMatrix = SetSparseMatrix( boundaryBoundaryTriplets , numBoundayTexels , numBoundayTexels , false );
	if( verbose ) printf( "\t Parsing =  %.4f \n" , timer.elapsed() );

	if( verbose ) timer.reset();
	coarseBoundaryBoundaryMatrix += partialRestriction;
	if( verbose ) printf( "\t Adding =  %.4f\n" , timer.elapsed() );
}

template< typename GeometryReal , typename MatrixReal >
void FullMatrixConstruction( const GridAtlas< GeometryReal , MatrixReal > &gridAtlas , const SystemCoefficients< MatrixReal > &systemCoefficients , SparseMatrix< MatrixReal , int > &fullMatrix )
{
	int numTexels = gridAtlas.numTexels;
	fullMatrix.resize( numTexels );
	std::vector< int > deepGlobalIndex = gridAtlas.deepGlobalIndex;

	//Add deep triplets
	//Initalization deep
	for (int i = 0; i < deepGlobalIndex.size(); i++) fullMatrix.SetRowSize(deepGlobalIndex[i], 9);

	const std::vector<RasterLine> & rasterLines = gridAtlas.rasterLines;
	for (int r = 0; r < rasterLines.size(); r++) {
		const RasterLine & currentLine = rasterLines[r];
		int lineLength = currentLine.lineEndIndex - currentLine.lineStartIndex + 1;
		int previousLineStart = currentLine.prevLineIndex;
		int currentLineStart = currentLine.lineStartIndex;
		int nextLineStart = currentLine.nextLineIndex;
		for (int i = 0; i < lineLength; i++) {
			fullMatrix[currentLineStart][0].N = previousLineStart - 1;
			fullMatrix[currentLineStart][1].N = previousLineStart;
			fullMatrix[currentLineStart][2].N = previousLineStart + 1;

			fullMatrix[currentLineStart][3].N = currentLineStart - 1;
			fullMatrix[currentLineStart][4].N = currentLineStart;
			fullMatrix[currentLineStart][5].N = currentLineStart + 1;

			fullMatrix[currentLineStart][6].N = nextLineStart - 1;
			fullMatrix[currentLineStart][7].N = nextLineStart;
			fullMatrix[currentLineStart][8].N = nextLineStart + 1;

			previousLineStart++;
			currentLineStart++;
			nextLineStart++;
		}
	}

	//Update deep
	for (int r = 0; r < rasterLines.size(); r++) {
		const RasterLine & currentLine = rasterLines[r];
		int deepOffset = currentLine.coeffStartIndex;
		int lineLength = currentLine.lineEndIndex - currentLine.lineStartIndex + 1;
		int previousLineStart = currentLine.prevLineIndex;
		int currentLineStart = currentLine.lineStartIndex;
		int nextLineStart = currentLine.nextLineIndex;
		const MatrixReal *coefficients = systemCoefficients.deepCoefficients.data() + deepOffset * 10;
		for (int i = 0; i < lineLength; i++) {
			for (int k = 0; k < 9; k++) fullMatrix[currentLineStart][k].Value = coefficients[k];

			previousLineStart++;
			currentLineStart++;
			nextLineStart++;
			coefficients += 10;
		}
	}

	std::vector< Eigen::Triplet< MatrixReal > > boundaryTriplets;
	std::vector<int> boundaryGlobalIndex = gridAtlas.boundaryGlobalIndex;
	for (int i = 0; i < boundaryGlobalIndex.size(); i++) {
		int globalIndex = boundaryGlobalIndex[i];
		for( int j=0 ; j<systemCoefficients.boundaryDeepMatrix.RowSize(i) ; j++ )
			boundaryTriplets.push_back( Eigen::Triplet< MatrixReal >( globalIndex , systemCoefficients.boundaryDeepMatrix[i][j].N , systemCoefficients.boundaryDeepMatrix[i][j].Value ) );
		for( int j=0 ; j<systemCoefficients.boundaryBoundaryMatrix.RowSize(i) ; j++ )
		{
			int neighbourGlobalIndex = boundaryGlobalIndex[ systemCoefficients.boundaryBoundaryMatrix[i][j].N ];
			boundaryTriplets.push_back( Eigen::Triplet< MatrixReal >( globalIndex , neighbourGlobalIndex , systemCoefficients.boundaryBoundaryMatrix[i][j].Value ) );
		}
	}

	SparseMatrix< MatrixReal , int > boundaryMatrix = SetSparseMatrix( boundaryTriplets , numTexels , numTexels , false );

	fullMatrix += boundaryMatrix;

	SparseMatrix< MatrixReal , int > _fullMatrix;
	CompressSparseMatrix( fullMatrix , _fullMatrix );
	fullMatrix = _fullMatrix;
}

template< typename GeometryReal , typename MatrixReal >
void UpdateMultigridCoefficients( const HierarchicalSystem< GeometryReal , MatrixReal > &hierarchy , std::vector< SystemCoefficients< MatrixReal > > &multigridCoefficients , const SystemCoefficients< MatrixReal > &inputCoefficients , SparseMatrix< MatrixReal , int > &coarseSystemMatrix , bool useDeepReciprocals, bool verbose=false)
{
	int levels = (int)hierarchy.gridAtlases.size();

	multigridCoefficients.resize( levels );

	multigridCoefficients[0].deepCoefficients = inputCoefficients.deepCoefficients;
	multigridCoefficients[0].boundaryBoundaryMatrix = inputCoefficients.boundaryBoundaryMatrix;
	multigridCoefficients[0].boundaryDeepMatrix = inputCoefficients.boundaryDeepMatrix;

	SparseMatrix< MatrixReal , int > currentSystemMatrix;
	Miscellany::Timer timer;

	for (int i = 1; i < levels; i++)
	{
		if (verbose) printf("Level %d \n", i);

		const GridAtlas< GeometryReal , MatrixReal > &gridAtlas = hierarchy.gridAtlases[i];
		int numTexels = gridAtlas.numTexels;
		int numDeepTexels = gridAtlas.numDeepTexels;
		int numBoundaryTexels = gridAtlas.numBoundaryTexels;

		const SystemCoefficients< MatrixReal > &fineCoefficients = multigridCoefficients[i-1];
		SystemCoefficients< MatrixReal > &coarseCoefficients = multigridCoefficients[i];

		// Deep restriction
		coarseCoefficients.deepCoefficients.resize( numDeepTexels * 10 , 0 );
		if( verbose ) timer.reset();
		DeepCoefficientRestriction( fineCoefficients.deepCoefficients , coarseCoefficients.deepCoefficients , gridAtlas.deepLines );
		if( verbose ) printf( "Deep Deep Restriction %d  =  %.4f\n" , i , timer.elapsed() );

		// Boundary Deep Matrix
		if( verbose ) timer.reset();
		BoundaryDeepMatrixConstruction( numBoundaryTexels , numTexels , coarseCoefficients.deepCoefficients , hierarchy.boundaryDeepIndices[i] , coarseCoefficients.boundaryDeepMatrix );
		if( verbose ) printf( "Boundary Deep Restriction %d  =  %.4f \n" , i , timer.elapsed() );

		// Boundary Boundary Matrix
		if( verbose ) timer.reset();
		BoundaryBoundaryMatrixConstruction( numBoundaryTexels , fineCoefficients.deepCoefficients , fineCoefficients.boundaryBoundaryMatrix , hierarchy.boundaryCoarseFineProlongation[i] , hierarchy.boundaryFineCoarseRestriction[i-1] , hierarchy.boundaryBoundaryIndices[i] , coarseCoefficients.boundaryBoundaryMatrix );
		if( verbose ) printf( "Boundary Boundary Restriction %d  =  %.4f \n" , i , timer.elapsed() );

		if( i==(levels-1) )	FullMatrixConstruction( gridAtlas , coarseCoefficients , coarseSystemMatrix );
	}

	if( useDeepReciprocals )
	{
		if( verbose ) timer.reset();
		//Set Deep coefficients to reciprocal
		for (int i = 0; i < levels - 1; i++)
		{
			std::vector< MatrixReal > &deepCoefficients = multigridCoefficients[i].deepCoefficients;
			int numDeepTexels = hierarchy.gridAtlases[i].numDeepTexels;
#pragma omp parallel for
			for (int ii=0 ; ii<numDeepTexels ; ii++ )
			{
				MatrixReal centerValue = deepCoefficients[10 * ii + 4];
				MatrixReal reciprocal = (MatrixReal)1. / centerValue;
				deepCoefficients[10 * ii + 0] *= reciprocal;
				deepCoefficients[10 * ii + 1] *= reciprocal;
				deepCoefficients[10 * ii + 2] *= reciprocal;
				deepCoefficients[10 * ii + 3] *= reciprocal;
				deepCoefficients[10 * ii + 4]  = reciprocal;
				deepCoefficients[10 * ii + 5] *= reciprocal;
				deepCoefficients[10 * ii + 6] *= reciprocal;
				deepCoefficients[10 * ii + 7] *= reciprocal;
				deepCoefficients[10 * ii + 8] *= reciprocal;
				deepCoefficients[10 * ii + 9]  = centerValue;
			}
		}
		if( verbose ) printf( "Deep coeff inversion =  %.4f\n" , timer.elapsed() );
	}
}


template< class DirectSolver , typename GeometryReal , typename MatrixReal >
void UpdateMultigridSolvers( const HierarchicalSystem< GeometryReal , MatrixReal > &hierarchy , const SparseMatrix< MatrixReal , int > &coarseSystemMatrix , const std::vector< SystemCoefficients< MatrixReal > > &multigridCoefficients , VCycleSolvers< DirectSolver > &vCycleSolvers , bool verbose=false , bool initSolvers=true )
{
	Miscellany::Timer timer;
	double cumInitTime = 0;
	double initTime;

	double cumUpdateTime = 0;
	double updateTime;

	timer.reset();
	if( initSolvers ) vCycleSolvers.coarse.init( coarseSystemMatrix );
	initTime = timer.elapsed();
	cumInitTime += initTime;

	timer.reset();
	vCycleSolvers.coarse.update( coarseSystemMatrix );
	updateTime = timer.elapsed();
	cumUpdateTime += updateTime;

	int levels = (int)hierarchy.gridAtlases.size();
	vCycleSolvers.boundary.resize( levels-1 );
	for( int i=0 ; i<levels-1 ; i++ )
	{
		timer.reset();
		if( initSolvers ) vCycleSolvers.boundary[i].init( multigridCoefficients[i].boundaryBoundaryMatrix );
		initTime = timer.elapsed();
		cumInitTime += initTime;

		timer.reset();
		vCycleSolvers.boundary[i].update( multigridCoefficients[i].boundaryBoundaryMatrix );
		updateTime = timer.elapsed();
		cumUpdateTime += updateTime;
	}

	if( verbose ) printf( "Hierarchy solvers init time %.4f\n" , cumInitTime );
	if( verbose ) printf( "Hierarchy solvers update time %.4f\n" , cumUpdateTime );
}


template< class DirectSolver , typename GeometryReal , typename MatrixReal >
void UpdateMultigridCoefficientsAndSolvers( const HierarchicalSystem< GeometryReal , MatrixReal > &hierarchy , std::vector< SystemCoefficients< MatrixReal > > &multigridCoefficients , const SystemCoefficients< MatrixReal > & systemCoefficients , VCycleSolvers< DirectSolver > &vCycleSolvers , bool detailVerbose=false , bool initSolvers=true )
{
	SparseMatrix< MatrixReal , int > coarseSystemMatrix;
	Miscellany::Timer timer;
	UpdateMultigridCoefficients( hierarchy , multigridCoefficients , systemCoefficients , coarseSystemMatrix , true , detailVerbose );
	if( detailVerbose ) printf( "Hierarchy coefficients update =  %.4f \n" , timer.elapsed() );

	timer.reset();
	UpdateMultigridSolvers( hierarchy , coarseSystemMatrix , multigridCoefficients , vCycleSolvers , detailVerbose , initSolvers );
	if( detailVerbose ) printf( "Hierarchy solvers =  %.4f \n" , timer.elapsed() );
}

template< typename GeometryReal , typename MatrixReal  , class DirectSolver >
void UpdateLinearSystem
(
	MatrixReal screenWeight , MatrixReal stiffnessWeight , const HierarchicalSystem< GeometryReal , MatrixReal > &hierarchy , std::vector< SystemCoefficients< MatrixReal > > &multigridCoefficients ,
	const SystemCoefficients< MatrixReal > &mass , const SystemCoefficients< MatrixReal > &stiffness ,
	VCycleSolvers< DirectSolver > &vCycleSolvers , DirectSolver &fineSolver ,
	const SparseMatrix< MatrixReal , int > &fineSystemMatrix ,
	bool detailVerbose=false , bool initSolvers=true , bool updateFineSolver=false
)
{
	SystemCoefficients< MatrixReal > fineCoefficients;
	int numDeepCoefficients = (int)mass.deepCoefficients.size();
	fineCoefficients.deepCoefficients.resize( numDeepCoefficients );

#pragma omp parallel for
	for( int i=0 ; i<numDeepCoefficients ; i++ ) fineCoefficients.deepCoefficients[i] = mass.deepCoefficients[i] * screenWeight + stiffness.deepCoefficients[i] * stiffnessWeight;

	const SparseMatrix< MatrixReal , int >* in[][2] = { { &mass.boundaryBoundaryMatrix , &stiffness.boundaryBoundaryMatrix } , { &mass.boundaryDeepMatrix , &stiffness.boundaryDeepMatrix } };
	SparseMatrix< MatrixReal , int >* out[] = { &fineCoefficients.boundaryBoundaryMatrix , &fineCoefficients.boundaryDeepMatrix };
	for( int ii=0 ; ii<2 ; ii++ )
	{
		const SparseMatrix< MatrixReal , int >& _in1 = *(in[ii][0]);
		const SparseMatrix< MatrixReal , int >& _in2 = *(in[ii][1]);
		SparseMatrix< MatrixReal , int >& _out = *(out[ii]);
		_out.resize( _in1.rows );
#pragma omp parallel for
		for( int i=0 ; i<_out.rows ; i++ )
		{
			_out.SetRowSize( i , _in1.rowSizes[i] );
			for( int j=0 ; j<_out.rowSizes[i] ; j++ ) _out[i][j] = MatrixEntry< MatrixReal , int >( _in1[i][j].N , _in1[i][j].Value * screenWeight + _in2[i][j].Value * stiffnessWeight );
		}
	}
	UpdateMultigridCoefficientsAndSolvers( hierarchy , multigridCoefficients , fineCoefficients , vCycleSolvers , detailVerbose , initSolvers );
	if( updateFineSolver )
	{
		Miscellany::Timer timer;
		if( initSolvers ) fineSolver.init(fineSystemMatrix);
		if( detailVerbose ) printf( "Analyze =  %.4f\n" , timer.elapsed() );

		if( detailVerbose ) timer.reset();
		fineSolver.update(fineSystemMatrix);
		if( detailVerbose ) printf( "Factorize =  %.4f\n" , timer.elapsed() );
	}
}
