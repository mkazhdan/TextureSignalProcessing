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

template< class Real >
int DeepCoefficientRestriction( const std::vector< Real >& fineDeepCoefficients , std::vector< Real >& coarseDeepCoefficients , const std::vector< DeepLine > & deepLines )
{

	int numCoarseDeepCoeffs = (int)coarseDeepCoefficients.size();

	auto UpdateRow = [&](int r) {
		const Real * prevLineCoeff = fineDeepCoefficients.data() + deepLines[r].finePrevLineIndex * 10;
		const Real * currentLineCoeff = fineDeepCoefficients.data() + deepLines[r].fineCurrentLineIndex * 10;
		const Real * nextLineCoeff = fineDeepCoefficients.data() + deepLines[r].fineNextLineIndex * 10;
		Real * coarseLineCoeff = coarseDeepCoefficients.data() + deepLines[r].coarseLineStartIndex * 10;

		const int lineLength = deepLines[r].coarseLineEndIndex - deepLines[r].coarseLineStartIndex + 1;

		auto UpdateNode = [&](int i) {
			const Real * prevCoeff = prevLineCoeff + 20 * i;
			const Real * currentCoeff = currentLineCoeff + 20 * i;
			const Real * nextCoeff = nextLineCoeff + 20 * i;
			Real * coarseCoeff = coarseLineCoeff + 10 * i;
			double cumValues[] =
			{
				0,0,0,
				0,0,0,
				0,0,0
			};

			//Procces middle cell node
			for (int k = 0; k < 2; k++) for (int l = 0; l < 2; l++) {
				const Real * nodeCoeff = (l == 0) ? prevCoeff + 10 * (2 * k - 1) : nextCoeff + 10 * (2 * k - 1);
				cumValues[3 * l + k] += (nodeCoeff[0] + nodeCoeff[1] * 0.5 + nodeCoeff[3] * 0.5 + nodeCoeff[4] * 0.25) * 0.25;
				cumValues[3 * l + k + 1] += (nodeCoeff[1] * 0.5 + nodeCoeff[2] + nodeCoeff[4] * 0.25 + nodeCoeff[5] * 0.5) * 0.25;
				cumValues[3 * (l + 1) + k] += (nodeCoeff[3] * 0.5 + nodeCoeff[4] * 0.25 + nodeCoeff[6] + nodeCoeff[7] * 0.5) * 0.25;
				cumValues[3 * (l + 1) + k + 1] += (nodeCoeff[4] * 0.25 + nodeCoeff[5] * 0.5 + nodeCoeff[7] * 0.5 + nodeCoeff[8]) * 0.25;
			}
			//Procces vertical edge node
			for (int l = 0; l < 2; l++) {
				const Real * nodeCoeff = (l == 0) ? prevCoeff : nextCoeff;
				cumValues[3 * l] += (nodeCoeff[0] * 0.5 + nodeCoeff[3] * 0.25)*0.5;
				cumValues[3 * l + 1] += (nodeCoeff[0] * 0.5 + nodeCoeff[1] + nodeCoeff[2] * 0.5 + nodeCoeff[3] * 0.25 + nodeCoeff[4] * 0.5 + nodeCoeff[5] * 0.25)*0.5;
				cumValues[3 * l + 2] += (nodeCoeff[2] * 0.5 + nodeCoeff[5] * 0.25)*0.5;

				cumValues[3 * (l + 1)] += (nodeCoeff[3] * 0.25 + nodeCoeff[6] * 0.5) * 0.5;
				cumValues[3 * (l + 1) + 1] += (nodeCoeff[3] * 0.25 + nodeCoeff[4] * 0.5 + nodeCoeff[5] * 0.25 + nodeCoeff[6] * 0.5 + nodeCoeff[7] + nodeCoeff[8] * 0.5)*0.5;
				cumValues[3 * (l + 1) + 2] += (nodeCoeff[5] * 0.25 + nodeCoeff[8] * 0.5) * 0.5;
			}
			//Procces horizontal edge node
			for (int k = 0; k < 2; k++) {
				const Real * nodeCoeff = (k == 0) ? currentCoeff - 10 : currentCoeff + 10;
				cumValues[k] += (nodeCoeff[0] * 0.5 + nodeCoeff[1] * 0.25)*0.5;
				cumValues[3 + k] += (nodeCoeff[0] * 0.5 + nodeCoeff[1] * 0.25 + nodeCoeff[3] + nodeCoeff[4] * 0.5 + nodeCoeff[6] * 0.5 + nodeCoeff[7] * 0.25)*0.5;
				cumValues[6 + k] += (nodeCoeff[6] * 0.5 + nodeCoeff[7] * 0.25)*0.5;

				cumValues[k + 1] += (nodeCoeff[1] * 0.25 + nodeCoeff[2] * 0.5) * 0.5;
				cumValues[3 + k + 1] += (nodeCoeff[1] * 0.25 + nodeCoeff[2] * 0.5 + nodeCoeff[4] * 0.5 + nodeCoeff[5] + nodeCoeff[7] * 0.25 + nodeCoeff[8] * 0.5)*0.5;
				cumValues[6 + k + 1] += (nodeCoeff[7] * 0.25 + nodeCoeff[8] * 0.5) * 0.5;
			}
			//Procces center node
			{
				const Real * nodeCoeff = currentCoeff;
				cumValues[0] += (nodeCoeff[0] * 0.25);
				cumValues[1] += (nodeCoeff[0] * 0.25 + nodeCoeff[1] * 0.5 + nodeCoeff[2] * 0.25);
				cumValues[2] += (nodeCoeff[2] * 0.25);
				cumValues[3] += (nodeCoeff[0] * 0.25 + nodeCoeff[3] * 0.5 + nodeCoeff[6] * 0.25);

				cumValues[4] += (nodeCoeff[0] * 0.25 + nodeCoeff[1] * 0.5 + nodeCoeff[2] * 0.25 + nodeCoeff[3] * 0.5 + nodeCoeff[4] + nodeCoeff[5] * 0.5 +
					nodeCoeff[6] * 0.25 + nodeCoeff[7] * 0.5 + nodeCoeff[8] * 0.25);

				cumValues[5] += (nodeCoeff[2] * 0.25 + nodeCoeff[5] * 0.5 + nodeCoeff[8] * 0.25);
				cumValues[6] += (nodeCoeff[6] * 0.25);
				cumValues[7] += (nodeCoeff[6] * 0.25 + nodeCoeff[7] * 0.5 + nodeCoeff[8] * 0.25);
				cumValues[8] += nodeCoeff[8] * 0.25;
			}
			for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) {
				coarseCoeff[3 * l + k] = cumValues[3 * l + k];
			}
		};
		for (int i = 0; i < lineLength;i++)UpdateNode(i);
	};
	int threads = omp_get_max_threads();
#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		const int tId = omp_get_thread_num();
		const int firstLine = ((int)deepLines.size() * tId) / threads;
		const int lastLine = ((int)deepLines.size() * (tId + 1)) / threads;
		for (int r = firstLine; r < lastLine; r++) UpdateRow(r);
	}

	//for (int r = 0; r < deepLines.size(); r++) UpdateRow(r);

	return 1;
}

int BoundaryDeepMatrixConstruction(const int numBoundayTexels, const int numTexels, const std::vector<double> & deepCoefficients, const std::vector<BoundaryDeepIndex> & boundaryDeepIndices, SparseMatrix<double, int> & boundaryDeepMatrix) {
	std::vector<Eigen::Triplet<double>> boundaryDeepTriplets;
	for (int i = 0; i < boundaryDeepIndices.size(); i++) {
		boundaryDeepTriplets.push_back(Eigen::Triplet<double>(boundaryDeepIndices[i].boundaryIndex, boundaryDeepIndices[i].deepGlobalIndex, deepCoefficients[10 * boundaryDeepIndices[i].deepIndex + boundaryDeepIndices[i].offset])); //NOTE: Deep coefficient is not multiplied by reciprocal yet.
	}
	boundaryDeepMatrix = SetSparseMatrix( boundaryDeepTriplets , numBoundayTexels , numTexels , false );

	return 1;
}


int BoundaryBoundaryMatrixConstruction(const std::vector<int> & boundaryGlobalIndex, const int numBoundayTexels, const std::vector<double> & fineDeepCoefficients, const SparseMatrix<double, int> & fineBoundaryBoundaryMatrix, const SparseMatrix<double, int> & boundaryCoarseFineProlongation, const SparseMatrix<double, int> & boundaryFineCoarseRestriction, const std::vector<BoundaryBoundaryIndex> & boundaryBoundaryIndices, SparseMatrix<double, int> & coarseBoundaryBoundaryMatrix, bool verbose = false) {

	clock_t t_begin;

	if (verbose)t_begin = clock();

	SparseMatrix<double, int> partialRestriction = fineBoundaryBoundaryMatrix * boundaryCoarseFineProlongation;
	partialRestriction = boundaryFineCoarseRestriction *partialRestriction;

	if (verbose) printf("\t Multiplication  =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);

	if (verbose) t_begin = clock();
	std::vector<Eigen::Triplet<double>> boundaryBoundaryTriplets;
	for (int i = 0; i < boundaryBoundaryIndices.size(); i++) {
		boundaryBoundaryTriplets.push_back(Eigen::Triplet<double>(boundaryBoundaryIndices[i].coarsePrincipalBoundaryIndex, boundaryBoundaryIndices[i].coarseSecondaryBoundaryIndex, fineDeepCoefficients[boundaryBoundaryIndices[i].fineDeepIndex * 10 + boundaryBoundaryIndices[i].offset] * boundaryBoundaryIndices[i].weight)); //NOTE: Deep coefficient is not multiplied by reciprocal yet.
	}
	if (verbose) printf("\t Triplets  =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);

	if (verbose) t_begin = clock();
	coarseBoundaryBoundaryMatrix = SetSparseMatrix( boundaryBoundaryTriplets , numBoundayTexels , numBoundayTexels , false );
	if (verbose) printf("\t Parsing =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);

	if (verbose) t_begin = clock();
	coarseBoundaryBoundaryMatrix += partialRestriction;
	if (verbose) printf("\t Adding =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);

	return 1;
}

int FullMatrixConstruction(const GridAtlas & gridAtlas, const std::vector<double> & deepCoefficients, const SparseMatrix<double, int> & boundaryBoundaryMatrix, const SparseMatrix<double, int> & boundaryDeepMatrix, SparseMatrix<double, int> & fullMatrix) {
	int numTexels = gridAtlas.numTexels;
	fullMatrix.resize(numTexels);
	std::vector<int> deepGlobalIndex = gridAtlas.deepGlobalIndex;

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
		const double * coefficients = deepCoefficients.data() + deepOffset * 10;
		for (int i = 0; i < lineLength; i++) {
			for (int k = 0; k < 9; k++) fullMatrix[currentLineStart][k].Value = coefficients[k];

			previousLineStart++;
			currentLineStart++;
			nextLineStart++;
			coefficients += 10;
		}
	}

	std::vector<Eigen::Triplet<double>> boundaryTriplets;
	std::vector<int> boundaryGlobalIndex = gridAtlas.boundaryGlobalIndex;
	for (int i = 0; i < boundaryGlobalIndex.size(); i++) {
		int globalIndex = boundaryGlobalIndex[i];
		for (int j = 0; j < boundaryDeepMatrix.RowSize(i); j++) {
			boundaryTriplets.push_back(Eigen::Triplet<double>(globalIndex, boundaryDeepMatrix[i][j].N, boundaryDeepMatrix[i][j].Value));
		}
		for (int j = 0; j < boundaryBoundaryMatrix.RowSize(i); j++) {
			int neighbourGlobalIndex = boundaryGlobalIndex[boundaryBoundaryMatrix[i][j].N];
			boundaryTriplets.push_back(Eigen::Triplet<double>(globalIndex, neighbourGlobalIndex, boundaryBoundaryMatrix[i][j].Value));
		}
	}

	SparseMatrix< double , int > boundaryMatrix = SetSparseMatrix( boundaryTriplets , numTexels , numTexels , false );

	fullMatrix += boundaryMatrix;

	SparseMatrix<double, int> _fullMatrix;
	CompressSparseMatrix(fullMatrix, _fullMatrix);
	fullMatrix = _fullMatrix;
	return 1;
}


template<class Real>
int CompareMatrices(const SparseMatrix<Real, int> & M1, const SparseMatrix<Real, int> & M2, const std::vector<GridNodeInfo> & nodeInfo, const std::vector<int> & boundaryAndDeepIndex) {
	double threshold = 1e-10;

	if (M1.Rows() != M2.Rows()) {
		printf("ERROR: Different num rows! %d %d \n", (int)M1.Rows(), (int)M2.Rows());
		return 0;
	}
	for (int i = 0; i < M1.Rows(); i++) {
		if (M1.RowSize(i) != M2.RowSize(i)) {
			printf("ERROR: Different row sizes! %d : %d %d \n", i, (int)M1.RowSize(i), (int)M2.RowSize(i));
			printf("Node Type %d.Chart Id %d.Pos(%d, %d) \n", nodeInfo[i].nodeType, nodeInfo[i].chartId, nodeInfo[i].ci, nodeInfo[i].cj);
			printf("Boundary and deep index %d \n", boundaryAndDeepIndex[i]);
			printf("Matrix 1 \n");
			for (int j1 = 0; j1 < M1.RowSize(i); j1++) {
				printf("(%d,%g) \n", M1[i][j1].N, M1[i][j1].Value);
				GridNodeInfo _nodeInfo = nodeInfo[M1[i][j1].N];
				printf("Node Type %d. Chart Id %d. Pos (%d ,%d) \n", _nodeInfo.nodeType, _nodeInfo.chartId, _nodeInfo.ci, _nodeInfo.cj);
			}
			printf("Matrix 2 \n");
			for (int j2 = 0; j2 < M2.RowSize(i); j2++) {
				printf("(%d,%g) \n", M2[i][j2].N, M1[i][j2].Value);
				GridNodeInfo _nodeInfo = nodeInfo[M2[i][j2].N];
				printf("Node Type %d. Chart Id %d. Pos (%d ,%d) \n", _nodeInfo.nodeType, _nodeInfo.chartId, _nodeInfo.ci, _nodeInfo.cj);
			}
			return 0;
		}

		if (0) {
			printf("\t \t Row %d \n", i);
			for (int j = 0; j < M1.RowSize(i); j++) {
				printf("(%d %g) (%d %g) \n", M1[i][j].N, M1[i][j].Value, M2[i][j].N, M2[i][j].Value);
				GridNodeInfo _nodeInfo = nodeInfo[M1[i][j].N];
				printf("Neighbour Node type %d Chart Id %d Pos (%d ,%d) \n", _nodeInfo.nodeType, _nodeInfo.chartId, _nodeInfo.ci, _nodeInfo.cj);
			}
		}

		for (int j = 0; j < M1.RowSize(i); j++) {
			if (M1[i][j].N != M2[i][j].N) {
				printf("ERROR: Different column values! %d %d : %d %d \n", i, j, M1[i][j].N, M2[i][j].N);
				printf("Matrix 1 \n");
				for (int j1 = 0; j1 < M1.RowSize(i); j1++) printf("(%d,%g \n", M1[i][j1].N, M1[i][j1].Value);
				printf("Matrix 2 \n");
				for (int j2 = 0; j2 < M2.RowSize(i); j2++) printf("(%d,%g \n", M2[i][j2].N, M2[i][j2].Value);

				return 0;
			}
			if (fabs(M1[i][j].Value - M2[i][j].Value) > threshold) {
				printf("ERROR: Numerical imprecision! %d %d : %g %g \n", i, M1[i][j].N, M1[i][j].Value, M2[i][j].Value);
				GridNodeInfo _nodeInfo = nodeInfo[i];
				printf("Current Node type %d Chart Id %d Pos (%d ,%d) \n", _nodeInfo.nodeType, _nodeInfo.chartId, _nodeInfo.ci, _nodeInfo.cj);
				_nodeInfo = nodeInfo[M1[i][j].N];
				printf("Neighbour Node type %d Chart Id %d Pos (%d ,%d) \n", _nodeInfo.nodeType, _nodeInfo.chartId, _nodeInfo.ci, _nodeInfo.cj);
				return 0;
			}
		}
	}

	printf("Succesful comparison! \n");
	return 1;

}

int UpdateMultigridCoefficients(const HierarchicalSystem & hierarchy, std::vector<MultigridLevelCoefficients<double>> & multigridCoefficients, const std::vector<double> & inputDeepCoefficients, const SparseMatrix<double, int> & inputBoundaryBoundaryMatrix, const SparseMatrix<double, int> & inputBoundaryDeepMatrix, const SparseMatrix<double, int> & systemMatrix, SparseMatrix<double, int> & coarseSystemMatrix, bool useDeepReciprocals, bool verbose = false) {

	int levels = (int)hierarchy.gridAtlases.size();

	multigridCoefficients.resize(levels);

	multigridCoefficients[0].deepCoefficients = inputDeepCoefficients;
	multigridCoefficients[0].boundaryBoundaryMatrix = inputBoundaryBoundaryMatrix;
	multigridCoefficients[0].boundaryDeepMatrix = inputBoundaryDeepMatrix;

	if (0) {//For debugging
		printf("Enable for debugging! \n");
		SparseMatrix<double, int> fullMatrix;
		FullMatrixConstruction(hierarchy.gridAtlases[0], multigridCoefficients[0].deepCoefficients, multigridCoefficients[0].boundaryBoundaryMatrix, multigridCoefficients[0].boundaryDeepMatrix, fullMatrix);

		SparseMatrix<double, int> compressedSystemMatrix;
		CompressSparseMatrix(systemMatrix, compressedSystemMatrix);
		if (!CompareMatrices(fullMatrix, compressedSystemMatrix, hierarchy.gridAtlases[0].nodeInfo, hierarchy.gridAtlases[0].boundaryAndDeepIndex)) {
			return 0;
		}
	}

	SparseMatrix<double, int> currentSystemMatrix;
	if (0) {//For debugging
		printf("Enable for debugging! \n");
		currentSystemMatrix = systemMatrix;
	}
	clock_t t_begin;

	for (int i = 1; i < levels; i++) {
		if (verbose) printf("Level %d \n", i);

		const GridAtlas & gridAtlas = hierarchy.gridAtlases[i];
		int numTexels = gridAtlas.numTexels;
		int numDeepTexels = gridAtlas.numDeepTexels;
		int numBoundaryTexels = gridAtlas.numBoundaryTexels;

		//Deep restriction
		const std::vector<double> & fineDeepCoefficients = multigridCoefficients[i - 1].deepCoefficients;
		std::vector<double> & coarseDeepCoefficients = multigridCoefficients[i].deepCoefficients;
		coarseDeepCoefficients.resize(numDeepTexels * 10, 0);
		if (verbose) t_begin = clock();
		DeepCoefficientRestriction(fineDeepCoefficients, coarseDeepCoefficients, gridAtlas.deepLines);
		if (verbose) printf("Deep Deep Restriction %d  =  %.4f \n", i, double(clock() - t_begin) / CLOCKS_PER_SEC);

		//Boundary Deep Matrix
		SparseMatrix<double, int> & coarseBoundaryDeepMatrix = multigridCoefficients[i].boundaryDeepMatrix;
		if (verbose) t_begin = clock();
		BoundaryDeepMatrixConstruction(numBoundaryTexels, numTexels, coarseDeepCoefficients, hierarchy.boundaryDeepIndices[i], coarseBoundaryDeepMatrix);
		if (verbose) printf("Boundary Deep Restriction %d  =  %.4f \n", i, double(clock() - t_begin) / CLOCKS_PER_SEC);

		//Boundary Boundary Matrix
		const SparseMatrix<double, int> & fineBoundaryBoundaryMatrix = multigridCoefficients[i - 1].boundaryBoundaryMatrix;
		SparseMatrix<double, int> & coarseBoundaryBoundaryMatrix = multigridCoefficients[i].boundaryBoundaryMatrix;
		if (verbose) t_begin = clock();
		BoundaryBoundaryMatrixConstruction(gridAtlas.boundaryGlobalIndex, numBoundaryTexels, fineDeepCoefficients, fineBoundaryBoundaryMatrix, hierarchy.boundaryCoarseFineProlongation[i], hierarchy.boundaryFineCoarseRestriction[i - 1], hierarchy.boundaryBoundaryIndices[i], coarseBoundaryBoundaryMatrix);
		if (verbose) printf("Boundary Boundary Restriction %d  =  %.4f \n", i, double(clock() - t_begin) / CLOCKS_PER_SEC);


		if (i == (levels - 1)) {
			FullMatrixConstruction(gridAtlas, coarseDeepCoefficients, coarseBoundaryBoundaryMatrix, coarseBoundaryDeepMatrix, coarseSystemMatrix);
		}

		if (0) {//For debugging
			printf("Enable for debugging! \n");
			SparseMatrix<double, int> fullMatrix;
			FullMatrixConstruction(gridAtlas, coarseDeepCoefficients, coarseBoundaryBoundaryMatrix, coarseBoundaryDeepMatrix, fullMatrix);
			std::vector<double> ones(numTexels, 1.0);
			std::vector<double> onesProduct(numTexels, 0);
			fullMatrix.Multiply(&ones[0], &onesProduct[0]);
			double cumMass = 0;
			for (int k = 0; k < numTexels; k++) cumMass += onesProduct[k];
			if (verbose) printf("Cum Mass Full Matrix %g \n", cumMass);

			currentSystemMatrix = hierarchy.prolongation[i - 1].transpose() * currentSystemMatrix * hierarchy.prolongation[i - 1];
			SparseMatrix<double, int> compressedSystemMatrix;
			CompressSparseMatrix(currentSystemMatrix, compressedSystemMatrix);


			compressedSystemMatrix.Multiply(&ones[0], &onesProduct[0]);
			cumMass = 0;
			for (int k = 0; k < numTexels; k++)cumMass += onesProduct[k];
			if (verbose) printf("Cum Mass System Matrix %g \n", cumMass);
			if (!CompareMatrices(fullMatrix, compressedSystemMatrix, hierarchy.gridAtlases[i].nodeInfo, hierarchy.gridAtlases[i].boundaryAndDeepIndex)) {
				return 0;
			}
		}
	}

	if (useDeepReciprocals) {
		if (verbose) t_begin = clock();
		//Set Deep coefficients to reciprocal
		for (int i = 0; i < levels - 1; i++) {
			std::vector<double> & deepCoefficients = multigridCoefficients[i].deepCoefficients;
			int numDeepTexels = hierarchy.gridAtlases[i].numDeepTexels;
#pragma omp parallel for
			for (int i = 0; i < numDeepTexels; i++) {
				double centerValue = deepCoefficients[10 * i + 4];
				double reciprocal = 1.0 / centerValue;
				deepCoefficients[10 * i] *= reciprocal;
				deepCoefficients[10 * i + 1] *= reciprocal;
				deepCoefficients[10 * i + 2] *= reciprocal;
				deepCoefficients[10 * i + 3] *= reciprocal;
				deepCoefficients[10 * i + 4] = reciprocal;
				deepCoefficients[10 * i + 5] *= reciprocal;
				deepCoefficients[10 * i + 6] *= reciprocal;
				deepCoefficients[10 * i + 7] *= reciprocal;
				deepCoefficients[10 * i + 8] *= reciprocal;
				deepCoefficients[10 * i + 9] = centerValue;
			}
		}
		if (verbose) printf("Deep coeff inversion =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
	}
	return 1;
}


template<class CoarseSolverType, class BoundarySolverType>
int UpdateMultigridSolvers(const HierarchicalSystem & hierarchy, const SparseMatrix<double, int> & coarseSystemMatrix, const std::vector<MultigridLevelCoefficients<double>> & multigridCoefficients, CoarseSolverType & coarseSolver, std::vector<BoundarySolverType> & boundarySolver, bool verbose = false, bool initSolvers = true) {

	clock_t p_begin;
	double cumInitTime = 0;
	double initTime;

	double cumUpdateTime = 0;
	double updateTime;

	p_begin = clock();
	if(initSolvers)coarseSolver.init(coarseSystemMatrix);
	initTime = double(clock() - p_begin) / CLOCKS_PER_SEC;
	cumInitTime += initTime;

	p_begin = clock();
	coarseSolver.update(coarseSystemMatrix);
	updateTime = double(clock() - p_begin) / CLOCKS_PER_SEC;
	cumUpdateTime += updateTime;

	int levels = (int)hierarchy.gridAtlases.size();
	boundarySolver.resize(levels - 1);
	for (int i = 0; i < levels - 1; i++) {
		p_begin = clock();
		if (initSolvers) boundarySolver[i].init(multigridCoefficients[i].boundaryBoundaryMatrix);
		initTime = double(clock() - p_begin) / CLOCKS_PER_SEC;
		cumInitTime += initTime;

		p_begin = clock();
		boundarySolver[i].update(multigridCoefficients[i].boundaryBoundaryMatrix);
		updateTime = double(clock() - p_begin) / CLOCKS_PER_SEC;
		cumUpdateTime += updateTime;
	}

	if(verbose) printf("Hierarchy solvers init time %.4f \n", cumInitTime);
	if(verbose) printf("Hierarchy solvers update time %.4f \n", cumUpdateTime);
	return 1;
}


template <class CoarseSolverType, class BoundarySolverType>
int UpdateMultigridCoefficientsAndSolvers(const HierarchicalSystem & hierarchy, std::vector<MultigridLevelCoefficients<double>> & multigridCoefficients, const std::vector<double> & deepCoefficients, const SparseMatrix<double, int> & boundaryBoundaryMatrix, const SparseMatrix<double, int> & boundaryDeepMatrix, CoarseSolverType & coarseSolver, std::vector<BoundarySolverType> & boundarySolver, const SparseMatrix<double, int> & fineSystemMatrix, bool detailVerbose = false, bool initSolvers = true) {
	SparseMatrix<double, int> coarseSystemMatrix;
	clock_t t_begin;
	if (detailVerbose)t_begin = clock();
	if (!UpdateMultigridCoefficients(hierarchy, multigridCoefficients, deepCoefficients, boundaryBoundaryMatrix, boundaryDeepMatrix, fineSystemMatrix, coarseSystemMatrix, true, detailVerbose)) {
		printf("ERROR: Unable to update multigrid coefficients! \n");
		return 0;
	}
	if(detailVerbose)printf("Hierarchy coefficients update =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);

	t_begin = clock();
	if (!UpdateMultigridSolvers(hierarchy, coarseSystemMatrix, multigridCoefficients, coarseSolver, boundarySolver, detailVerbose, initSolvers)) {
		printf("ERROR: Unable to update multigrid solvers! \n");
		return 0;
	}
	if (detailVerbose)printf("Hierarchy solvers =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);

	return 1;
}

template< class Real , class CoarseSolverType , class BoundarySolverType , class DirectSolverType >
typename std::enable_if< std::is_same< Real , float >::value , int >::type UpdateLinearSystem( double screenWeight , double stiffnessWeight , const HierarchicalSystem & hierarchy, std::vector<MultigridLevelCoefficients<Real>> & multigridCoefficients,
	const std::vector<double> & deepMassCoefficients, const std::vector<double> & deepStiffnessCoefficients,
	const SparseMatrix<double, int> & boundaryBoundaryMassMatrix,const SparseMatrix<double, int> & boundaryBoundaryStiffnessMatrix,
	const SparseMatrix<double, int> & boundaryDeepMassMatrix, const SparseMatrix<double, int> & boundaryDeepStiffnessMatrix,
	CoarseSolverType & coarseSolver, std::vector<BoundarySolverType> & boundarySolver, DirectSolverType & fineSolver,
	const SparseMatrix<double, int> & fineSystemMatrix, bool detailVerbose = false, bool initSolvers = true, bool updateFineSolver = false){
	std::vector<double> fineDeepCoefficients;
	int numDeepCoefficients = (int)deepMassCoefficients.size();
	fineDeepCoefficients.resize(numDeepCoefficients);

#pragma omp parallel for
	for (int i = 0; i < numDeepCoefficients; i++) fineDeepCoefficients[i] = deepMassCoefficients[i] * screenWeight + deepStiffnessCoefficients[i] * stiffnessWeight;


	SparseMatrix< double , int > fineBoundaryBoundaryMatrix , fineBoundaryDeepMatrix;
	const SparseMatrix< double , int >* in[][2] = { { &boundaryBoundaryMassMatrix , &boundaryBoundaryStiffnessMatrix } , { &boundaryDeepMassMatrix , &boundaryDeepStiffnessMatrix } };
	SparseMatrix< double , int >* out[] = { &fineBoundaryBoundaryMatrix , &fineBoundaryDeepMatrix };
	for( int ii=0 ; ii<2 ; ii++ )
	{
		const SparseMatrix< double , int >& _in1 = *(in[ii][0]);
		const SparseMatrix< double , int >& _in2 = *(in[ii][1]);
		SparseMatrix< double , int >& _out = *(out[ii]);
		_out.resize( _in1.rows );
#pragma omp parallel for
		for( int i=0 ; i<_out.rows ; i++ )
		{
			_out.SetRowSize( i , _in1.rowSizes[i] );
			for( int j=0 ; j<_out.rowSizes[i] ; j++ ) _out[i][j] = MatrixEntry< double , int >( _in1[i][j].N , _in1[i][j].Value * screenWeight + _in2[i][j].Value * stiffnessWeight );
		}
	}
	std::vector<MultigridLevelCoefficients<double>> _multigridCoefficients;
	if (!UpdateMultigridCoefficientsAndSolvers(hierarchy, _multigridCoefficients, fineDeepCoefficients, fineBoundaryBoundaryMatrix, fineBoundaryDeepMatrix, coarseSolver, boundarySolver, fineSystemMatrix, detailVerbose,initSolvers)) {
		printf("ERROR: Unable to initialize multigrid coefficients and solver! \n");
		return 0;
	}
	int levels = (int)hierarchy.gridAtlases.size();
	multigridCoefficients.resize(levels);
	for (int l = 0; l < levels; l++) {
		multigridCoefficients[l].boundaryDeepMatrix = _multigridCoefficients[l].boundaryDeepMatrix;
		multigridCoefficients[l].boundaryBoundaryMatrix = _multigridCoefficients[l].boundaryBoundaryMatrix;
		const std::vector<double> & _deepCoefficients = _multigridCoefficients[l].deepCoefficients;
		multigridCoefficients[l].deepCoefficients.resize(_deepCoefficients.size());
		for (int i = 0; i < (int)_deepCoefficients.size(); i++)multigridCoefficients[l].deepCoefficients[i] = (Real)_deepCoefficients[i];
	}

	if (updateFineSolver){
		clock_t t_begin;

		if (detailVerbose) t_begin = clock();
		if(initSolvers) fineSolver.init(fineSystemMatrix);
		if(detailVerbose)printf("Analyze =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);

		if (detailVerbose) t_begin = clock();
		fineSolver.update(fineSystemMatrix);
		if(detailVerbose)printf("Factorize =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
	}

	return 1;
}
template< class Real , class CoarseSolverType , class BoundarySolverType , class DirectSolverType >
typename std::enable_if< std::is_same< Real , double >::value , int >::type UpdateLinearSystem( double screenWeight , double stiffnessWeight , const HierarchicalSystem & hierarchy, std::vector<MultigridLevelCoefficients<Real>> & multigridCoefficients,
	const std::vector<double> & deepMassCoefficients, const std::vector<double> & deepStiffnessCoefficients,
	const SparseMatrix<double, int> & boundaryBoundaryMassMatrix,const SparseMatrix<double, int> & boundaryBoundaryStiffnessMatrix,
	const SparseMatrix<double, int> & boundaryDeepMassMatrix, const SparseMatrix<double, int> & boundaryDeepStiffnessMatrix,
	CoarseSolverType & coarseSolver, std::vector<BoundarySolverType> & boundarySolver, DirectSolverType & fineSolver,
	const SparseMatrix<double, int> & fineSystemMatrix, bool detailVerbose = false, bool initSolvers = true, bool updateFineSolver = false){
	std::vector<double> fineDeepCoefficients;
	int numDeepCoefficients = (int)deepMassCoefficients.size();
	fineDeepCoefficients.resize(numDeepCoefficients);

#pragma omp parallel for
	for (int i = 0; i < numDeepCoefficients; i++) fineDeepCoefficients[i] = deepMassCoefficients[i] * screenWeight + deepStiffnessCoefficients[i] * stiffnessWeight;

	SparseMatrix< double , int > fineBoundaryBoundaryMatrix , fineBoundaryDeepMatrix;
	const SparseMatrix< double , int >* in[][2] = { { &boundaryBoundaryMassMatrix , &boundaryBoundaryStiffnessMatrix } , { &boundaryDeepMassMatrix , &boundaryDeepStiffnessMatrix } };
	SparseMatrix< double , int >* out[] = { &fineBoundaryBoundaryMatrix , &fineBoundaryDeepMatrix };
	for( int ii=0 ; ii<2 ; ii++ )
	{
		const SparseMatrix< double , int >& _in1 = *(in[ii][0]);
		const SparseMatrix< double , int >& _in2 = *(in[ii][1]);
		SparseMatrix< double , int >& _out = *(out[ii]);
		_out.resize( _in1.rows );
#pragma omp parallel for
		for( int i=0 ; i<_out.rows ; i++ )
		{
			_out.SetRowSize( i , _in1.rowSizes[i] );
			for( int j=0 ; j<_out.rowSizes[i] ; j++ ) _out[i][j] = MatrixEntry< double , int >( _in1[i][j].N , _in1[i][j].Value * screenWeight + _in2[i][j].Value * stiffnessWeight );
		}
	}
	if (!UpdateMultigridCoefficientsAndSolvers(hierarchy, multigridCoefficients, fineDeepCoefficients, fineBoundaryBoundaryMatrix, fineBoundaryDeepMatrix, coarseSolver, boundarySolver, fineSystemMatrix, detailVerbose, initSolvers)) {
		printf("ERROR: Unable to initialize multigrid coefficients and solver! \n");
		return 0;
	}

	if (updateFineSolver){
		clock_t t_begin;

		if (detailVerbose) t_begin = clock();
		if(initSolvers) fineSolver.init(fineSystemMatrix);
		if(detailVerbose)printf("Analyze =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);

		if (detailVerbose) t_begin = clock();
		fineSolver.update(fineSystemMatrix);
		if(detailVerbose)printf("Factorize =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
	}

	return 1;
}
