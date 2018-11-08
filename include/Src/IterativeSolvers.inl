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

#include<Eigen/SparseLU>
#include "SparseMatrixParser.h"

template<class Real, class Data, class Solver>
int Relaxation(const std::vector<Real> & deepCoefficients, const SparseMatrix< Real, int> & boundaryDeepMatrix, Solver & boundarySolver, const std::vector<int> & boundaryGlobalIndex, const std::vector<SegmentedRasterLine> & segmentedLines, const std::vector<Data> & rhs, std::vector<Data> & x0, std::vector<Data> & boundaryRHS, std::vector<Data> & boundarySolution, std::vector<Data> & variableBoundaryRHS, const int numIterations = 2, bool boundaryFirst = true, bool verbose = false) {


	int numBoundaryVariables = boundaryGlobalIndex.size();
	
#pragma omp parallel for
	for (int i = 0; i < numBoundaryVariables; i++) {
		boundaryRHS[i] = rhs[boundaryGlobalIndex[i]];
	}

	clock_t t_begin;

	int it_offset = boundaryFirst ? 0 : 1;
	for (int it = 0; it < numIterations; it++){
		if ((it + it_offset) % 2 == 0){//Update Boundary;

			if (verbose) t_begin = clock();

#pragma omp parallel for 
			for (int i = 0; i < numBoundaryVariables; i++){
				boundarySolution[i] = x0[boundaryGlobalIndex[i]];
				variableBoundaryRHS[i] = boundaryRHS[i];
			}

			boundaryDeepMatrix.Multiply((Data*)&x0[0], (Data*)&variableBoundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);

			if (verbose) printf("\t Boundary initialization =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);

			if (verbose) t_begin = clock();
			solve(boundarySolver,boundarySolution, variableBoundaryRHS);
			if (verbose)printf("\t Boundary update =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);

#pragma omp parallel for
			for (int i = 0; i < numBoundaryVariables; i++) x0[boundaryGlobalIndex[i]] = boundarySolution[i];
		}


		if ((it + it_offset) % 2 == 1) {//Update Interior


			auto UpdateRow = [&](int k)
			{
				const std::vector<RasterLine> & rasterLines = segmentedLines[k].segments;
				for (int r = 0; r < rasterLines.size(); r++) {

					Data* _xCurrent = (Data*)x0.data() + rasterLines[r].lineStartIndex;
					const Data* _xPrevious = (Data*)x0.data() + rasterLines[r].prevLineIndex;
					const Data* _xNext = (Data*)x0.data() + rasterLines[r].nextLineIndex;
					const Data* _rhs = (Data*)rhs.data() + rasterLines[r].lineStartIndex;


					int lineLength = (rasterLines[r].lineEndIndex - rasterLines[r].lineStartIndex + 1);
					int lineDeepStart = rasterLines[r].coeffStartIndex;
					const Real* _deepCoefficients = &deepCoefficients[10 * lineDeepStart];
					for (int i = 0; i < lineLength; _deepCoefficients += 10, i++)
					{
						_xCurrent[i] =
							_rhs[i] * _deepCoefficients[4] -
							(
								_xPrevious[i - 1] * _deepCoefficients[0] +
								_xPrevious[i + 0] * _deepCoefficients[1] +
								_xPrevious[i + 1] * _deepCoefficients[2] +
								_xCurrent[i - 1] * _deepCoefficients[3] +
								_xCurrent[i + 1] * _deepCoefficients[5] +
								_xNext[i - 1] * _deepCoefficients[6] +
								_xNext[i + 0] * _deepCoefficients[7] +
								_xNext[i + 1] * _deepCoefficients[8]
								);
					}
				}
			};


			auto UpdateBlock = [&](const int firstLine, const int lastLine, const int passes){
				for (int r = firstLine; r < lastLine + passes - 1; r++) for (int p = 0; p < passes; p++) if (r - p >= firstLine && r - p < lastLine) UpdateRow(r - p);
			};

			if (verbose) t_begin = clock();

			const int blockSize = 16;
			const int totalLines = segmentedLines.size();
			const int numBlocks = ((totalLines - 1) / blockSize) + 1;

#pragma omp parallel for
			for (int b = 0; b < numBlocks; b+=2){
				const int firstLine = (totalLines * b)/ numBlocks;
				const int lastLine = (totalLines * (b + 1)) / numBlocks;
				UpdateBlock(firstLine, lastLine, 3);
			}

#pragma omp parallel for
			for (int b = 1; b < numBlocks; b += 2){
				const int firstLine = (totalLines * b) / numBlocks;
				const int lastLine = (totalLines * (b + 1)) / numBlocks;
				UpdateBlock(firstLine, lastLine, 3);
			}
			if (verbose) printf("\t GS Deep update =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
		}
	}

	return 1;
}




template<class Real, class Data, class Solver>
int Relaxation(const std::vector<Real> & deepCoefficients, const SparseMatrix< Real, int> & boundaryDeepMatrix, Solver & boundarySolver, const std::vector<int> & boundaryGlobalIndex, const std::vector<ThreadTask> & threadTasks, const std::vector<Data> & rhs, std::vector<Data> & x0, std::vector<Data> & boundaryRHS, std::vector<Data> & boundarySolution, std::vector<Data> & variableBoundaryRHS, const int numIterations = 2, bool boundaryFirst = true, bool verbose = false) {


	int numBoundaryVariables = (int)boundaryGlobalIndex.size();

#pragma omp parallel for
	for (int i = 0; i < numBoundaryVariables; i++) {
		boundaryRHS[i] = rhs[boundaryGlobalIndex[i]];
	}

	clock_t t_begin;

	int it_offset = boundaryFirst ? 0 : 1;
	for (int it = 0; it < numIterations; it++) {
		if ((it + it_offset) % 2 == 0) {//Update Boundary;

			if (verbose) t_begin = clock();

#pragma omp parallel for 
			for (int i = 0; i < numBoundaryVariables; i++) {
				boundarySolution[i] = x0[boundaryGlobalIndex[i]];
				variableBoundaryRHS[i] = boundaryRHS[i];
			}

			boundaryDeepMatrix.Multiply((Data*)&x0[0], (Data*)&variableBoundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);

			if (verbose) printf("\t Boundary initialization =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);

			if (verbose) t_begin = clock();
			solve(boundarySolver, boundarySolution, variableBoundaryRHS);
			if (verbose)printf("\t Boundary update =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);

#pragma omp parallel for
			for (int i = 0; i < numBoundaryVariables; i++) x0[boundaryGlobalIndex[i]] = boundarySolution[i];
		}


		if ((it + it_offset) % 2 == 1) {//Update Interior
			
			if (verbose) t_begin = clock();
			auto UpdateRow = [&](const BlockDeepSegmentedLine & deepSegmentedLine){
				for (int s = 0; s < deepSegmentedLine.blockDeepSegments.size(); s++) {
					const BlockDeepSegment & deepSegment = deepSegmentedLine.blockDeepSegments[s];
					int length = deepSegment.currentEnd - deepSegment.currentStart + 1;

					Data* _xCurrent = (Data*)x0.data() + deepSegment.currentStart;
					const Data* _xPrevious = (Data*)x0.data() + deepSegment.previousStart;
					const Data* _xNext = (Data*)x0.data() + deepSegment.nextStart;
					const Data* _rhs = (Data*)rhs.data() + deepSegment.currentStart;
					const Real* _deepCoefficients = &deepCoefficients[10 * deepSegment.deepStart];

					for (int i = 0; i < length; _deepCoefficients += 10, i++)
					{
						_xCurrent[i] =
							_rhs[i] * _deepCoefficients[4] -
							(
								_xPrevious[i - 1] * _deepCoefficients[0] +
								_xPrevious[i + 0] * _deepCoefficients[1] +
								_xPrevious[i + 1] * _deepCoefficients[2] +
								_xCurrent[i - 1] * _deepCoefficients[3] +
								_xCurrent[i + 1] * _deepCoefficients[5] +
								_xNext[i - 1] * _deepCoefficients[6] +
								_xNext[i + 0] * _deepCoefficients[7] +
								_xNext[i + 1] * _deepCoefficients[8]
								);
					}
				}
			};

			auto UpdateBlock = [&](int t, int b, int gsIters)
			{
				const std::vector<BlockDeepSegmentedLine> & blockDeepSegmentedLines = threadTasks[t].blockTasks[b].blockDeepSegmentedLines;
				for (int j = 0; j < blockDeepSegmentedLines.size() + gsIters; j++){
					for (int gsIter = 0; gsIter < gsIters; gsIter++) {
						if ((j - gsIter) >= 0 && (j - gsIter) < blockDeepSegmentedLines.size()) UpdateRow(blockDeepSegmentedLines[j - gsIter]);
					}
				}
			};

			//Update even blocks
#pragma omp parallel for
			for (int t = 0; t < (threadTasks.size() + 1) / 2; t++){
				for (int b = 0; b < threadTasks[2 * t].blockTasks.size(); b++) UpdateBlock(2 * t, b, 3);
			}

			//Update odd blocks
#pragma omp parallel for
			for (int t = 0; t < threadTasks.size() / 2; t++){
				for (int b = 0; b < threadTasks[2 * t + 1].blockTasks.size(); b++) UpdateBlock(2 * t + 1, b ,3);
			}
			if (verbose) printf("\t GS Deep update =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
		}
	}

	return 1;
}


template<class Real, class Data, class Solver>
int RelaxationAndResidual(const std::vector<Real> & deepCoefficients, const SparseMatrix< Real, int> & boundaryDeepMatrix, Solver & boundarySolver, const std::vector<int> & boundaryGlobalIndex, const std::vector<ThreadTask> & threadTasks, const std::vector<Data> & rhs, std::vector<Data> & x0, std::vector<Data> & boundaryRHS, std::vector<Data> & boundaryValue, std::vector<Data> & variableBoundaryRHS, const SparseMatrix< Real, int> & boundaryBoundaryMatrix, std::vector<Data> & residual, const int numIterations = 2, bool verbose = false) {

	int numBoundaryVariables = (int)boundaryGlobalIndex.size();

#pragma omp parallel for
	for (int i = 0; i < numBoundaryVariables; i++) {
		boundaryRHS[i] = rhs[boundaryGlobalIndex[i]];
	}

	clock_t t_begin;
	for (int it = 0; it < numIterations; it++) {
		if (it % 2 == 0) {//Update Boundary

			if (verbose) t_begin = clock();

#pragma omp parallel for
			for (int i = 0; i < numBoundaryVariables; i++) {
				boundaryValue[i] = x0[boundaryGlobalIndex[i]];
				variableBoundaryRHS[i] = boundaryRHS[i];
			}

			boundaryDeepMatrix.Multiply((Data*)&x0[0], (Data*)&variableBoundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);
			if (verbose) printf("\t Boundary initialization =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);


			if (verbose) t_begin = clock();
			solve(boundarySolver, boundaryValue, variableBoundaryRHS);

			if (verbose) printf("\t Boundary update =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
#pragma omp parallel for
			for (int i = 0; i < numBoundaryVariables; i++) x0[boundaryGlobalIndex[i]] = boundaryValue[i];
		}

		if (it % 2 == 1) {//Update Interior


			if (verbose) t_begin = clock();


			auto UpdateRow = [&](const BlockDeepSegmentedLine & deepSegmentedLine) {
				for (int s = 0; s < deepSegmentedLine.blockDeepSegments.size(); s++) {
					const BlockDeepSegment & deepSegment = deepSegmentedLine.blockDeepSegments[s];
					int length = deepSegment.currentEnd - deepSegment.currentStart + 1;

					Data* _xCurrent = (Data*)x0.data() + deepSegment.currentStart;
					const Data* _xPrevious = (Data*)x0.data() + deepSegment.previousStart;
					const Data* _xNext = (Data*)x0.data() + deepSegment.nextStart;
					const Data* _rhs = (Data*)rhs.data() + deepSegment.currentStart;
					const Real* _deepCoefficients = &deepCoefficients[10 * deepSegment.deepStart];

					for (int i = 0; i < length; _deepCoefficients += 10, i++)
					{
						_xCurrent[i] =
							_rhs[i] * _deepCoefficients[4] -
							(
								_xPrevious[i - 1] * _deepCoefficients[0] +
								_xPrevious[i + 0] * _deepCoefficients[1] +
								_xPrevious[i + 1] * _deepCoefficients[2] +
								_xCurrent[i - 1] * _deepCoefficients[3] +
								_xCurrent[i + 1] * _deepCoefficients[5] +
								_xNext[i - 1] * _deepCoefficients[6] +
								_xNext[i + 0] * _deepCoefficients[7] +
								_xNext[i + 1] * _deepCoefficients[8]
								);
					}
				}
			};

			auto UpdateResidual = [&](const BlockDeepSegmentedLine & deepSegmentedLine) {
				for (int s = 0; s < deepSegmentedLine.blockDeepSegments.size(); s++) {
					const BlockDeepSegment & deepSegment = deepSegmentedLine.blockDeepSegments[s];
					int length = deepSegment.currentEnd - deepSegment.currentStart + 1;

					Data* _residual = (Data*)residual.data() + deepSegment.currentStart;
					const Data* _xCurrent = (Data*)x0.data() + deepSegment.currentStart;
					const Data* _xPrevious = (Data*)x0.data() + deepSegment.previousStart;
					const Data* _xNext = (Data*)x0.data() + deepSegment.nextStart;
					const Data* _rhs = (Data*)rhs.data() + deepSegment.currentStart;
					const Real* _deepCoefficients = &deepCoefficients[10 * deepSegment.deepStart];

					for (int i = 0; i < length; _deepCoefficients += 10, i++)
					{
						_residual[i] = _rhs[i] -
							(
								_xPrevious[i - 1] * _deepCoefficients[0] +
								_xPrevious[i + 0] * _deepCoefficients[1] +
								_xPrevious[i + 1] * _deepCoefficients[2] +
								_xCurrent[i - 1] * _deepCoefficients[3] +
								_xCurrent[i] +
								_xCurrent[i + 1] * _deepCoefficients[5] +
								_xNext[i - 1] * _deepCoefficients[6] +
								_xNext[i + 0] * _deepCoefficients[7] +
								_xNext[i + 1] * _deepCoefficients[8]
								) *_deepCoefficients[9];
					}
				}
			};

			auto UpdateBlock = [&](int t, int b, int gsIters)
			{
				const std::vector<BlockDeepSegmentedLine> & blockDeepSegmentedLines = threadTasks[t].blockTasks[b].blockDeepSegmentedLines;
				const std::vector<BlockDeepSegmentedLine> & blockPaddedSegmentedLines = threadTasks[t].blockTasks[b].blockPaddedSegmentedLines;

				for (int j = 0; j < blockDeepSegmentedLines.size() + gsIters; j++) {
					for (int gsIter = 0; gsIter < gsIters; gsIter++) {
						if ((j - gsIter) >= 0 && (j - gsIter) < blockDeepSegmentedLines.size()) UpdateRow(blockDeepSegmentedLines[j - gsIter]);
					}
					if ((j - gsIters) >= 0 && (j - gsIters) < blockPaddedSegmentedLines.size()) UpdateResidual(blockPaddedSegmentedLines[j - gsIters]);
				}
				for (int j = (int)blockDeepSegmentedLines.size(); j < blockPaddedSegmentedLines.size(); j++)UpdateResidual(blockPaddedSegmentedLines[j]);
			};

			//Update even blocks
#pragma omp parallel for
			for (int t = 0; t < (threadTasks.size() + 1) / 2; t++) {
				for (int b = 0; b < threadTasks[2 * t].blockTasks.size(); b++) UpdateBlock(2 * t, b, 3);
			}
			//Update odd blocks
#pragma omp parallel for
			for (int t = 0; t < threadTasks.size() / 2; t++) {
				for (int b = 0; b < threadTasks[2 * t + 1].blockTasks.size(); b++) UpdateBlock(2 * t + 1, b, 3);
			}

			if (verbose) printf("\t GS Deep update =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
		}
	}

	{//Compute boundary residual
		if (verbose) t_begin = clock();
		boundaryDeepMatrix.Multiply((Data*)&x0[0], (Data*)&boundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);
		boundaryBoundaryMatrix.Multiply((Data*)&boundaryValue[0], (Data*)&boundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);
		for (int i = 0; i < numBoundaryVariables; i++)residual[boundaryGlobalIndex[i]] = boundaryRHS[i];
		if (verbose) printf("\t Boundary residual =  %.4f \n", double(clock() - t_begin) / (CLOCKS_PER_SEC));
	}

	return 1;
}


template<class Real, class Data, class Solver>
int RelaxationAndResidual(const std::vector<Real> & deepCoefficients, const SparseMatrix< Real, int> & boundaryDeepMatrix, Solver & boundarySolver, const std::vector<int> & boundaryGlobalIndex, const std::vector<SegmentedRasterLine> & segmentedLines, const std::vector<Data> & rhs, std::vector<Data> & x0, std::vector<Data> & boundaryRHS,std::vector<Data> & boundaryValue, std::vector<Data> & variableBoundaryRHS, const SparseMatrix< Real, int> & boundaryBoundaryMatrix, std::vector<Data> & residual, const int numIterations = 2, bool verbose = false) {

	int numBoundaryVariables = boundaryGlobalIndex.size();

#pragma omp parallel for
	for (int i = 0; i < numBoundaryVariables; i++){
		boundaryRHS[i] = rhs[boundaryGlobalIndex[i]];
	}

	clock_t t_begin;
	for (int it = 0; it < numIterations; it++) {
		if (it % 2 == 0) {//Update Boundary

			if (verbose) t_begin = clock();

#pragma omp parallel for
			for (int i = 0; i < numBoundaryVariables; i++){
				boundaryValue[i] = x0[boundaryGlobalIndex[i]];
				variableBoundaryRHS[i] = boundaryRHS[i];
			}

			boundaryDeepMatrix.Multiply((Data*)&x0[0], (Data*)&variableBoundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);
			if (verbose) printf("\t Boundary initialization =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);


			if (verbose) t_begin = clock();
			solve(boundarySolver, boundaryValue, variableBoundaryRHS);

			if (verbose) printf("\t Boundary update =  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
#pragma omp parallel for
			for (int i = 0; i < numBoundaryVariables; i++) x0[boundaryGlobalIndex[i]] = boundaryValue[i];
		}

		if (it % 2 == 1) {//Update Interior

										//static const int GSIterations = 3;
										//for (int r = 0; r < rasterLines.size() + GSIterations; r++)
										//	for (int j = 0; j < GSIterations; j++)
										//		if (r - j >= 0 && r - j < rasterLines.size()) UpdateRow(r - j);

			if (verbose) t_begin = clock();


			auto UpdateResidual = [&](int k)
			{
				const std::vector<RasterLine> & rasterLines = segmentedLines[k].segments;
				for (int r = 0; r < rasterLines.size(); r++) {

					Data* _residual = (Data*)residual.data() + rasterLines[r].lineStartIndex;
					const Data* _xCurrent = (Data*)x0.data() + rasterLines[r].lineStartIndex;
					const Data* _xPrevious = (Data*)x0.data() + rasterLines[r].prevLineIndex;
					const Data* _xNext = (Data*)x0.data() + rasterLines[r].nextLineIndex;
					const Data* _rhs = (Data*)rhs.data() + rasterLines[r].lineStartIndex;


					int lineLength = (rasterLines[r].lineEndIndex - rasterLines[r].lineStartIndex + 1);
					int lineDeepStart = rasterLines[r].coeffStartIndex;
					const Real* _deepCoefficients = &deepCoefficients[10 * lineDeepStart];
					for (int i = 0; i < lineLength; _deepCoefficients += 10, i++)
					{
						_residual[i] = _rhs[i] -
							(
								_xPrevious[i - 1] * _deepCoefficients[0] +
								_xPrevious[i + 0] * _deepCoefficients[1] +
								_xPrevious[i + 1] * _deepCoefficients[2] +
								_xCurrent[i - 1] * _deepCoefficients[3] +
								_xCurrent[i] +
								_xCurrent[i + 1] * _deepCoefficients[5] +
								_xNext[i - 1] * _deepCoefficients[6] +
								_xNext[i + 0] * _deepCoefficients[7] +
								_xNext[i + 1] * _deepCoefficients[8]
								) *_deepCoefficients[9];
					}
				}
			};


			auto UpdateRow = [&](int k)
			{
				const std::vector<RasterLine> & rasterLines = segmentedLines[k].segments;
				for (int r = 0; r < rasterLines.size(); r++) {

					Data* _xCurrent = (Data*)x0.data() + rasterLines[r].lineStartIndex;
					const Data* _xPrevious = (Data*)x0.data() + rasterLines[r].prevLineIndex;
					const Data* _xNext = (Data*)x0.data() + rasterLines[r].nextLineIndex;
					const Data* _rhs = (Data*)rhs.data() + rasterLines[r].lineStartIndex;


					int lineLength = (rasterLines[r].lineEndIndex - rasterLines[r].lineStartIndex + 1);
					int lineDeepStart = rasterLines[r].coeffStartIndex;
					const Real* _deepCoefficients = &deepCoefficients[10 * lineDeepStart];
					for (int i = 0; i < lineLength; _deepCoefficients += 10, i++)
					{
						_xCurrent[i] =
							_rhs[i] * _deepCoefficients[4] -
							(
								_xPrevious[i - 1] * _deepCoefficients[0] +
								_xPrevious[i + 0] * _deepCoefficients[1] +
								_xPrevious[i + 1] * _deepCoefficients[2] +
								_xCurrent[i - 1] * _deepCoefficients[3] +
								_xCurrent[i + 1] * _deepCoefficients[5] +
								_xNext[i - 1] * _deepCoefficients[6] +
								_xNext[i + 0] * _deepCoefficients[7] +
								_xNext[i + 1] * _deepCoefficients[8]
								);
					}
				}
			};

			auto UpdateHeaderBlock = [&](const int firstLine, const int lastLine, const int passes, const bool firstResidual,const bool prevResidual){
				for (int r = firstLine; r < firstLine + passes + 1; r++){
					for (int p = 0; p < passes; p++) if (r - p >= firstLine && r - p < lastLine) UpdateRow(r - p);
				}
				if(firstResidual)UpdateResidual(firstLine);
				if(prevResidual)UpdateResidual(firstLine - 1);
			};

			auto UpdateMiddleBlock = [&](const int firstLine, const int lastLine, const int passes) {
				for (int r = firstLine + passes + 1; r < lastLine + passes - 1; r++) {
					for (int p = 0; p < passes; p++) if (r - p > firstLine && r - p < lastLine) UpdateRow(r - p);
					if (r - passes> firstLine && r - passes < lastLine) UpdateResidual(r - passes);
				}
			};

			auto UpdateSuffixBlock = [&](const int lastLine,const bool lastResidual, const bool nextResidual) {
				if (lastResidual)UpdateResidual(lastLine - 1);
				if (nextResidual)UpdateResidual(lastLine);
			};


			const int blockSize = 16;
			const int totalLines = segmentedLines.size();
			const int numBlocks = ((totalLines - 1) / blockSize) + 1;

			//Even blocks
#pragma omp parallel for
			for (int b = 0; b < numBlocks; b+=2){
				const int firstLine = (totalLines * b)/ numBlocks;
				const int lastLine = (totalLines * (b + 1)) / numBlocks;
				UpdateHeaderBlock(firstLine, lastLine, 3, firstLine == 0, false);
				UpdateMiddleBlock(firstLine, lastLine, 3);
				UpdateSuffixBlock(lastLine, lastLine == totalLines, false);
			}

			//Odd blocks
#pragma omp parallel for
			for (int b = 1; b < numBlocks; b += 2) {
				const int firstLine = (totalLines * b) / numBlocks;
				const int lastLine = (totalLines * (b + 1)) / numBlocks;

				UpdateHeaderBlock(firstLine, lastLine, 3, true, firstLine != 0);
				UpdateMiddleBlock(firstLine, lastLine, 3);
				UpdateSuffixBlock(lastLine, true, lastLine != totalLines);
			}
			if (verbose) printf("\t GS Deep update =  %.4f \n", double(clock() - t_begin) / (CLOCKS_PER_SEC));
		}
	}

	{//Compute boundary residual
		if (verbose) t_begin = clock();
		boundaryDeepMatrix.Multiply((Data*)&x0[0], (Data*)&boundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);
		boundaryBoundaryMatrix.Multiply((Data*)&boundaryValue[0], (Data*)&boundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);
		for (int i = 0; i < numBoundaryVariables; i++)residual[boundaryGlobalIndex[i]] = boundaryRHS[i];
		if (verbose) printf("\t Boundary residual =  %.4f \n", double(clock() - t_begin) / (CLOCKS_PER_SEC));
	}

	return 1;
}

template< class Real , class Data >
int MultiplyBySystemMatrix( const SystemCoefficients< double > &systemCoefficients , const std::vector< int > &boundaryGlobalIndex , const std::vector< RasterLine > &rasterLines , const std::vector< Data > &in , std::vector< Data > &out , bool verbose=false )
{

	clock_t t_begin;

	int numBoundaryVariables = boundaryGlobalIndex.size();

	std::vector<Data> inBoundaryValues;
	inBoundaryValues.resize(numBoundaryVariables);
	for (int i = 0; i < numBoundaryVariables; i++) inBoundaryValues[i] = in[boundaryGlobalIndex[i]];

	std::vector<Data> outBoundaryValues;
	outBoundaryValues.resize(numBoundaryVariables);

	if (verbose) t_begin = clock();
	systemCoefficients.boundaryBoundaryMatrix.Multiply(&inBoundaryValues[0], &outBoundaryValues[0]);
	systemCoefficients.boundaryDeepMatrix.Multiply(&in[0], &outBoundaryValues[0], MULTIPLY_ADD);
	if (verbose) printf("\t Multiply boundary =  %.4f \n", double(clock() - t_begin) / (CLOCKS_PER_SEC));

#pragma omp parallel for
	for (int i = 0; i < numBoundaryVariables; i++) out[boundaryGlobalIndex[i]] = outBoundaryValues[i];

	auto UpdateRow = [&](int r)
	{
		Data* _out = (Data*)out.data() + rasterLines[r].lineStartIndex;
		const Data* _inCurrent = (Data*)in.data() + rasterLines[r].lineStartIndex;
		const Data* _inPrevious = (Data*)in.data() + rasterLines[r].prevLineIndex;
		const Data* _inNext = (Data*)in.data() + rasterLines[r].nextLineIndex;

		int lineLength = (rasterLines[r].lineEndIndex - rasterLines[r].lineStartIndex + 1);
		int lineDeepStart = rasterLines[r].coeffStartIndex;
		const Real* _deepCoefficients = &systemCoefficients.deepCoefficients[10 * lineDeepStart];
		for (int i = 0; i < lineLength; _deepCoefficients += 10, i++)
		{
			_out[i] =
				(_inPrevious[i - 1] * _deepCoefficients[0] +
				_inPrevious[i] * _deepCoefficients[1] +
				_inPrevious[i + 1] * _deepCoefficients[2] +
				_inCurrent[i - 1] * _deepCoefficients[3] +
				_inCurrent[i] +
				_inCurrent[i + 1] * _deepCoefficients[5] +
				_inNext[i - 1] * _deepCoefficients[6] +
				_inNext[i] * _deepCoefficients[7] +
				_inNext[i + 1] * _deepCoefficients[8])*_deepCoefficients[9];
		}
	};

	if (verbose) t_begin = clock();

	int threads = omp_get_max_threads();
	std::vector<int> lineRange(threads + 1);
	int blockSize = rasterLines.size() / threads;
	for (int t = 0; t < threads; t++) lineRange[t] = t*blockSize;
	lineRange[threads] = rasterLines.size();
#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		const int tId = omp_get_thread_num();
		const int firstLine = lineRange[tId];
		const int lastLine = lineRange[tId + 1];
		for (int r = firstLine; r < lastLine; r++) UpdateRow(r);
	}
	if (verbose) printf("\t Multiply deep %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
	return 1;
}

template< class Real , class Data >
int MultiplyBySystemMatrix_NoReciprocals
(
	const SystemCoefficients< Real > &systemCoefficients ,
	const std::vector< int >& boundaryGlobalIndex ,
	const std::vector< RasterLine >& rasterLines ,
	const std::vector< Data >& in ,
	std::vector< Data > & out ,
	bool verbose=false
)
{
	clock_t t_begin;

	int numBoundaryVariables = (int)boundaryGlobalIndex.size();

	// Pull out the boundary values from the input array
	std::vector< Data > outBoundaryValues( numBoundaryVariables );
	std::vector< Data >  inBoundaryValues( numBoundaryVariables );
	for( int i=0 ; i<numBoundaryVariables ; i++ ) inBoundaryValues[i] = in[ boundaryGlobalIndex[i] ];

	t_begin = clock();
	//  Perform the boundary -> boundary multiplication
	systemCoefficients.boundaryBoundaryMatrix.Multiply( &inBoundaryValues[0] , &outBoundaryValues[0] );
	// Perform the interior -> boundary multiplication
	systemCoefficients.boundaryDeepMatrix.Multiply( &in[0] , &outBoundaryValues[0] , MULTIPLY_ADD );
	if( verbose ) printf( "\tMultiply boundary = %.4f\n" , double(clock() - t_begin) / (CLOCKS_PER_SEC) );

	// Write the boundary values back into the output array
#pragma omp parallel for
	for( int i=0 ; i<numBoundaryVariables ; i++ ) out[ boundaryGlobalIndex[i] ] = outBoundaryValues[i];

	// Perform the interior -> interior multiplication
	auto UpdateRow = [&]( int r )
	{

		Data* _out = (Data*)out.data() + rasterLines[r].lineStartIndex;
		const Data* _inCurrent  = (Data*)in.data() + rasterLines[r].lineStartIndex;
		const Data* _inPrevious = (Data*)in.data() + rasterLines[r].prevLineIndex;
		const Data* _inNext     = (Data*)in.data() + rasterLines[r].nextLineIndex;

		int lineLength = ( rasterLines[r].lineEndIndex - rasterLines[r].lineStartIndex + 1 );
		int lineDeepStart = rasterLines[r].coeffStartIndex;
		const Real* _deepCoefficients = &systemCoefficients.deepCoefficients[10*lineDeepStart];
		for( int i=0 ; i<lineLength ; _deepCoefficients+=10 , i++)
		{
			_out[i] =(Data)
				(
					_inPrevious[i-1] * _deepCoefficients[0] +
					_inPrevious[i+0] * _deepCoefficients[1] +
					_inPrevious[i+1] * _deepCoefficients[2] +
					_inCurrent [i-1] * _deepCoefficients[3] +
					_inCurrent [i+0] * _deepCoefficients[4] +
					_inCurrent [i+1] * _deepCoefficients[5] +
					_inNext    [i-1] * _deepCoefficients[6] +
					_inNext    [i+0] * _deepCoefficients[7] +
					_inNext    [i+1] * _deepCoefficients[8]
				);
		}
	};


	t_begin = clock();

	int threads = omp_get_max_threads();
	std::vector< int > lineRange( threads+1 );
	int blockSize = (int)rasterLines.size() / threads;
	for( int t=0 ; t<threads ; t++ ) lineRange[t] = t*blockSize;
	lineRange[threads] = (int)rasterLines.size();
#pragma omp parallel for
	for( int t=0 ; t<threads ; t++ )
	{
		const int tId = omp_get_thread_num();
		const int firstLine = lineRange[tId];
		const int lastLine = lineRange[tId+1];
		for( int r=firstLine ; r<lastLine ; r++ ) UpdateRow(r);
	}
	if( verbose ) printf( "\tMultiply deep = %.4f \n" , double(clock() - t_begin) / CLOCKS_PER_SEC );
	return 1;
}


template<class Real, class Data>
int ComputeSystemResdiual(const std::vector<Real> & deepCoefficients, const SparseMatrix< Real, int> & boundaryDeepMatrix, const SparseMatrix< Real, int> & boundaryBoundaryMatrix, const std::vector<int> & boundaryGlobalIndex, const std::vector<RasterLine> & rasterLines, const std::vector<Data> & boundaryRHS, const std::vector<Data> & boundaryValues, const std::vector<Data> & in, const std::vector<Data> & rhs, std::vector<Data> & out, bool verbose = false) {

	clock_t t_begin;
	//Boundary residual
	{
		int numBoundaryVariables = boundaryGlobalIndex.size();
		if (verbose) t_begin = clock();
		boundaryDeepMatrix.Multiply((Data*)&in[0], (Data*)&boundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);
		boundaryBoundaryMatrix.Multiply((Data*)&boundaryValues[0], (Data*)&boundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);
		for (int i = 0; i < numBoundaryVariables; i++)out[boundaryGlobalIndex[i]] = boundaryRHS[i];
		if (verbose) printf("\t Residual Boundary =  %.4f \n", double(clock() - t_begin) / (CLOCKS_PER_SEC));
	}

	auto UpdateResidual = [&](int r)
	{

		Data* _out = (Data*)out.data() + rasterLines[r].lineStartIndex;
		const Data* _inCurrent = (Data*)in.data() + rasterLines[r].lineStartIndex;
		const Data* _inPrevious = (Data*)in.data() + rasterLines[r].prevLineIndex;
		const Data* _inNext = (Data*)in.data() + rasterLines[r].nextLineIndex;
		const Data* _rhs = (Data*)rhs.data() + rasterLines[r].lineStartIndex;

		int lineLength = (rasterLines[r].lineEndIndex - rasterLines[r].lineStartIndex + 1);
		int lineDeepStart = rasterLines[r].coeffStartIndex;
		const Real* _deepCoefficients = &deepCoefficients[10 * lineDeepStart];
		for (int i = 0; i < lineLength; _deepCoefficients += 10, i++)
		{
			_out[i] = _rhs[i] -
				(_inPrevious[i - 1] * _deepCoefficients[0] +
					_inPrevious[i] * _deepCoefficients[1] +
					_inPrevious[i + 1] * _deepCoefficients[2] +
					_inCurrent[i - 1] * _deepCoefficients[3] +
					_inCurrent[i] +
					_inCurrent[i + 1] * _deepCoefficients[5] +
					_inNext[i - 1] * _deepCoefficients[6] +
					_inNext[i] * _deepCoefficients[7] +
					_inNext[i + 1] * _deepCoefficients[8])*_deepCoefficients[9];
		}
	};


	if (verbose) t_begin = clock();

	int threads = omp_get_max_threads();
#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		const int tId = omp_get_thread_num();
		const int firstLine = (rasterLines.size() *tId) / threads;
		const int lastLine = (rasterLines.size() * (tId + 1)) / threads;
		for (int r = firstLine; r < lastLine; r++) UpdateResidual(r);
	}
	if (verbose) printf("\t Residual Deep %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
	return 1;
}

template<class Real, class Data>
int MultiplyByRestriction(const SparseMatrix< Real, int> & __boundaryRestrictionMatrix, const std::vector<int> & boundaryGlobalIndex, std::vector<Data> & boundaryValue, const std::vector<RasterLine> & restrictionLines, const std::vector<Data> & in, std::vector<Data> & out, bool verbose = false) {

	clock_t t_begin;

	int numBoundaryVariables = (int)boundaryGlobalIndex.size();

	//std::vector<Data> boundaryValues;
	//boundaryValues.resize(numBoundaryVariables);

	//ConstPointer(Data) _in = (ConstPointer(Data))GetPointer(in);
	//Pointer(Data) _out = (Pointer(Data))GetPointer(out);
	//Pointer(Data) _boundaryValues = (Pointer(Data))GetPointer(boundaryValues);

	//if (verbose) t_begin = clock();
	//__boundaryRestrictionMatrix.Multiply(_in, _boundaryValues);
	//if (verbose) printf("\t Restriction boundary  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
//#pragma omp parallel for
	//for (int i = 0; i < numBoundaryVariables; i++) _out[boundaryGlobalIndex[i]] = _boundaryValues[i];

	if (verbose) t_begin = clock();
	__boundaryRestrictionMatrix.Multiply((Data*)&in[0], (Data*)&boundaryValue[0]);
	if (verbose) printf("\t Restriction boundary  %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);

#pragma omp parallel for
	for (int i = 0; i < numBoundaryVariables; i++) out[boundaryGlobalIndex[i]] = boundaryValue[i];


	auto UpdateRow = [&](int r)
	{
		Data* _out = (Data*)out.data() + restrictionLines[r].coeffStartIndex;
		const Data* _inCurrent = (Data*)in.data() + restrictionLines[r].lineStartIndex;
		const Data* _inPrevious = (Data*)in.data() + restrictionLines[r].prevLineIndex;
		const Data* _inNext = (Data*)in.data() + restrictionLines[r].nextLineIndex;

		int lineLength = ((restrictionLines[r].lineEndIndex - restrictionLines[r].lineStartIndex) / 2 + 1);
		int j = 0;
		for (int i = 0; i < lineLength; i++, j += 2)
		{
			_out[i] =
				_inCurrent[j] +
				(_inPrevious[j] + _inCurrent[j - 1] + _inCurrent[j + 1] + _inNext[j])* (Real)0.5 +
				(_inPrevious[j - 1] + _inPrevious[j + 1] + _inNext[j - 1] + _inNext[j + 1])* (Real)0.25;
		}
	};

	if (verbose) t_begin = clock();
	int threads = omp_get_max_threads();
	std::vector<int> lineRange(threads + 1);
	int blockSize = (int)restrictionLines.size() / threads;
	for (int t = 0; t < threads; t++) lineRange[t] = t*blockSize;
	lineRange[threads] = (int)restrictionLines.size();
#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		const int tId = omp_get_thread_num();
		const int firstLine = lineRange[tId];
		const int lastLine = lineRange[tId + 1];
		for (int r = firstLine; r < lastLine; r++) UpdateRow(r);
	}

	if (verbose)printf("\t Restriction deep %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
	return 1;
}

template<class Real, class Data>
int MultiplyByProlongation(const std::vector<ProlongationLine> & prolongationLines, const std::vector<Data> & in, std::vector<Data> & out, bool verbose = false) {

	clock_t t_begin;

	
	auto UpdateRow = [&](int r)
	{
		Data* _out = out.data() + prolongationLines[r].startIndex;
		int lineLenght = prolongationLines[r].length;
		int offset = prolongationLines[r].alignedStart ? 0 : 1;

		if (prolongationLines[r].centerLineIndex != -1) {
			const Data* _inCurrent = in.data() + prolongationLines[r].centerLineIndex;
			for (int i = 0; i < lineLenght; i++) {
				int prev = ((i + offset) / 2);
				int next = ((i + offset + 1) / 2);
				_out[i] = (_inCurrent[prev] + _inCurrent[next]) * Real(0.5);
			}
		}
		else {
			const Data* _inPrevious = in.data() + prolongationLines[r].prevLineIndex;
			const Data* _inNext = in.data() + prolongationLines[r].nextLineIndex;
			for (int i = 0; i < lineLenght; i++) {
				int prev =  ((i + offset) / 2);
				int next =  ((i + offset + 1) / 2);
				_out[i] = (_inPrevious[prev]  + _inPrevious[next] + _inNext[prev] + _inNext[next]) * Real(0.25);
			}
		}
	};

	if (verbose) t_begin = clock();
	int threads = omp_get_max_threads();
	std::vector<int> lineRange(threads + 1);
	int blockSize = prolongationLines.size() / threads;
	for (int t = 0; t < threads; t++) lineRange[t] = t*blockSize;
	lineRange[threads] = prolongationLines.size();
#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		const int tId = omp_get_thread_num();
		const int firstLine = lineRange[tId];
		const int lastLine = lineRange[tId + 1];
		for (int r = firstLine; r < lastLine; r++) UpdateRow(r);
	}
	if (verbose) printf("\t Prolongation deep %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC); 
	return 1;
}


template<class Real, class Data>
int AddProlongation(const std::vector<ProlongationLine> & prolongationLines, const std::vector<Data> & in, std::vector<Data> & out, bool verbose = false) {

	clock_t t_begin;


	auto UpdateRow = [&](int r)
	{
		Data* _out = out.data() + prolongationLines[r].startIndex;
		int lineLenght = prolongationLines[r].length;
		int offset = prolongationLines[r].alignedStart ? 0 : 1;

		if (prolongationLines[r].centerLineIndex != -1) {
			const Data* _inCurrent = in.data() + prolongationLines[r].centerLineIndex;
			for (int i = 0; i < lineLenght; i++) {
				int prev = ((i + offset) / 2);
				int next = ((i + offset + 1) / 2);
				_out[i] += (_inCurrent[prev] + _inCurrent[next]) * Real(0.5);
			}
		}
		else {
			const Data* _inPrevious = in.data() + prolongationLines[r].prevLineIndex;
			const Data* _inNext = in.data() + prolongationLines[r].nextLineIndex;
			for (int i = 0; i < lineLenght; i++) {
				int prev = ((i + offset) / 2);
				int next = ((i + offset + 1) / 2);
				_out[i] += (_inPrevious[prev] + _inPrevious[next] + _inNext[prev] + _inNext[next]) * Real(0.25);
			}
		}
	};

	if (verbose) t_begin = clock();
	int threads = omp_get_max_threads();
	std::vector<int> lineRange(threads + 1);
	int blockSize = (int)prolongationLines.size() / threads;
	for (int t = 0; t < threads; t++) lineRange[t] = t*blockSize;
	lineRange[threads] = (int)prolongationLines.size();
#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		const int tId = omp_get_thread_num();
		const int firstLine = lineRange[tId];
		const int lastLine = lineRange[tId + 1];
		for (int r = firstLine; r < lastLine; r++) UpdateRow(r);
	}
	if (verbose) printf("\t Prolongation deep %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
	return 1;
}

template< class Real , unsigned int Channels , class Data >
int CellStiffnessToTexelStiffness
(
	const std::vector< Real >& cellSharpenningMask ,
	const std::vector< InteriorTexelToCellLine >& interiorTexelToCellLines ,
	const std::vector< Data >& interiorTexelToCellCoeffs ,
	SparseMatrix< Real , int > boundaryCellBasedStiffnessRHSMatrix[Channels] ,
	std::vector< Real > boundaryTexelValues[Channels] ,
	const std::vector< int > & boundaryGlobalIndex ,
	std::vector< Data >& texelModulatedStiffness ,
	bool verbose=false
)
{
	//Update Boundary Texels
	clock_t t_begin; 
	if( verbose ) t_begin = clock();
#pragma omp parallel for
	for( int c=0 ; c<Channels ; c++ ) boundaryCellBasedStiffnessRHSMatrix[c].Multiply( &cellSharpenningMask[0] , &boundaryTexelValues[c][0] );
#pragma omp parallel for
	for( int i=0 ; i<boundaryGlobalIndex.size() ; i++ ) for( int c=0 ; c<Channels ; c++ ) texelModulatedStiffness[ boundaryGlobalIndex[i] ][c] = boundaryTexelValues[c][i];
	if( verbose ) printf( "\tBoundary updated: %.3f(s) \n" , double(clock() - t_begin) / CLOCKS_PER_SEC );
	//Update Interior Texels

	auto UpdateRow = [&](int r)
	{
		Data* out = texelModulatedStiffness.data() + interiorTexelToCellLines[r].texelStartIndex;
		const Real* previousCellRow = cellSharpenningMask.data() + interiorTexelToCellLines[r].previousCellStartIndex;
		const Real* nextCellRow = cellSharpenningMask.data() + interiorTexelToCellLines[r].nextCellStartIndex;
		const Data * coeff = interiorTexelToCellCoeffs.data() + interiorTexelToCellLines[r].coeffOffset * 4;
		int lineLenght = interiorTexelToCellLines[r].texelEndIndex - interiorTexelToCellLines[r].texelStartIndex + 1;
		for (int i = 0; i < lineLenght; i++) {
			out[i] = coeff[4 * i + 0] * previousCellRow[i] + coeff[4 * i + 1] * previousCellRow[i + 1] +
				coeff[4 * i + 2] * nextCellRow[i] + coeff[4 * i + 3] * nextCellRow[i + 1];
		}
	};

	if (verbose) t_begin = clock();
	int threads = omp_get_max_threads();
#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		const int tId = omp_get_thread_num();
		const int firstLine = (int)(interiorTexelToCellLines.size() *tId) / threads;
		const int lastLine = (int)(interiorTexelToCellLines.size() * (tId + 1)) / threads;
		for (int r = firstLine; r < lastLine; r++) UpdateRow(r);
	}
	if (verbose) printf("\t Interior update %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
	return 1;
}
template< class Real >
int CellStiffnessToTexelStiffness
(
	const std::vector< Real >& cellSharpenningMask ,
	const std::vector< InteriorTexelToCellLine >& interiorTexelToCellLines ,
	const std::vector< Real >& interiorTexelToCellCoeffs ,
	SparseMatrix< Real , int > boundaryCellBasedStiffnessRHSMatrix ,
	std::vector< Real > boundaryTexelValues ,
	const std::vector< int > & boundaryGlobalIndex ,
	std::vector< Real >& texelModulatedStiffness ,
	bool verbose=false
)
{
	//Update Boundary Texels
	clock_t t_begin; 
	if( verbose ) t_begin = clock();
	boundaryCellBasedStiffnessRHSMatrix.Multiply( &cellSharpenningMask[0] , &boundaryTexelValues[0] );
#pragma omp parallel for
	for( int i=0 ; i<boundaryGlobalIndex.size() ; i++ ) texelModulatedStiffness[ boundaryGlobalIndex[i] ] = boundaryTexelValues[i];
	if( verbose ) printf( "\tBoundary updated: %.3f(s) \n" , double(clock() - t_begin) / CLOCKS_PER_SEC );
	//Update Interior Texels

	auto UpdateRow = [&](int r)
	{
		Real* out = texelModulatedStiffness.data() + interiorTexelToCellLines[r].texelStartIndex;
		const Real* previousCellRow = cellSharpenningMask.data() + interiorTexelToCellLines[r].previousCellStartIndex;
		const Real* nextCellRow = cellSharpenningMask.data() + interiorTexelToCellLines[r].nextCellStartIndex;
		const Real * coeff = interiorTexelToCellCoeffs.data() + interiorTexelToCellLines[r].coeffOffset * 4;
		int lineLenght = interiorTexelToCellLines[r].texelEndIndex - interiorTexelToCellLines[r].texelStartIndex + 1;
		for (int i = 0; i < lineLenght; i++) {
			out[i] = coeff[4 * i + 0] * previousCellRow[i] + coeff[4 * i + 1] * previousCellRow[i + 1] +
				coeff[4 * i + 2] * nextCellRow[i] + coeff[4 * i + 3] * nextCellRow[i + 1];
		}
	};

	if (verbose) t_begin = clock();
	int threads = omp_get_max_threads();
#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		const int tId = omp_get_thread_num();
		const int firstLine = (int)(interiorTexelToCellLines.size() *tId) / threads;
		const int lastLine = (int)(interiorTexelToCellLines.size() * (tId + 1)) / threads;
		for (int r = firstLine; r < lastLine; r++) UpdateRow(r);
	}
	if (verbose) printf("\t Interior update %.4f \n", double(clock() - t_begin) / CLOCKS_PER_SEC);
	return 1;
}

template< class Real , class DataType , class BoundarySolver , class CoarseSolver >
int VCycle(	std::vector< MultigridLevelVariables< DataType > > &variables , const std::vector< SystemCoefficients< Real > > &coefficients , const std::vector< MultigridLevelIndices< Real > > &indices , VCycleSolvers< BoundarySolver , CoarseSolver > &vCycleSolvers , bool verbose , bool detailVerbose )
{
	clock_t p_begin;

	int levels = (int)variables.size();

	//Set x0[i = 1 : levels -1] = 0
	if (verbose) p_begin = clock();
	for (int i = 1; i < levels; i++){
		memset(&variables[i].x[0], 0, variables[i].x.size() * sizeof(DataType));
	}
	if (verbose) printf("Zero arrays %.4f \n", double(clock() - p_begin) / CLOCKS_PER_SEC);

	//Reduction phase
	for( int i=0 ; i<levels-1 ; i++ )
	{
		const SystemCoefficients< Real > & _coefficients = coefficients[i];
		const MultigridLevelIndices<Real> & _indices = indices[i];
		MultigridLevelVariables<DataType> & _variables = variables[i];

		const MultigridLevelIndices<Real> & nextLevelIndices = indices[i + 1];
		MultigridLevelVariables<DataType> & nextLevelVariables = variables[i + 1];

		if (verbose) printf("Level %d \n", i);
		clock_t p_begin;

		if (verbose) p_begin = clock();
		if( !RelaxationAndResidual( _coefficients.deepCoefficients , _coefficients.boundaryDeepMatrix , vCycleSolvers.boundary[i], _indices.boundaryGlobalIndex, _indices.threadTasks, _variables.rhs, _variables.x, _variables.boundary_rhs, _variables.boundary_value, _variables.variable_boundary_value, _coefficients.boundaryBoundaryMatrix, _variables.residual, 2, detailVerbose)) {
			printf("ERROR : Failed hybrid gauss seidel solver! \n");
			return 0;
		}
		if (verbose) printf("Relaxation  + Residual %.4f  \n", double(clock() - p_begin) / CLOCKS_PER_SEC);

		if (verbose) p_begin = clock();
		MultiplyByRestriction(_indices.boundaryRestriction, nextLevelIndices.boundaryGlobalIndex, nextLevelVariables.boundary_value, nextLevelIndices.restrictionLines, _variables.residual, nextLevelVariables.rhs, detailVerbose);
		if (verbose) printf("Restriction %.4f  \n", double(clock() - p_begin) / CLOCKS_PER_SEC);
	}

	//Prolongation phase
	for (int i = levels - 1; i >= 0; i--) {
		if (verbose) printf("Level %d \n", i);
		if (i < levels - 1){
		//if (i < levels - 1 && i > 0){
			clock_t p_begin;

			const SystemCoefficients< Real > & _coefficients = coefficients[i];
			const MultigridLevelIndices<Real> & _indices = indices[i];
			MultigridLevelVariables<DataType> & _variables = variables[i];

			if (verbose) p_begin = clock();
			if( !Relaxation( _coefficients.deepCoefficients , _coefficients.boundaryDeepMatrix , vCycleSolvers.boundary[i] , _indices.boundaryGlobalIndex , _indices.threadTasks , _variables.rhs , _variables.x , _variables.boundary_rhs , _variables.boundary_value , _variables.variable_boundary_value , 2 , true , detailVerbose ) )
			{
				fprintf( stderr , "[ERROR] Failed hybrid gauss seidel solver!\n" );
				return 0;
			}
			if (verbose) printf("Gauss Seidel %.4f  \n", double(clock() - p_begin) / CLOCKS_PER_SEC);
		}
		else if( i==levels-1 )
		{
			MultigridLevelVariables<DataType> & _variables = variables[i];
			clock_t p_begin;
			if (verbose) p_begin = clock();
			solve( vCycleSolvers.coarse , _variables.x , _variables.rhs );
			if (verbose) printf("Direct solver %.4f  \n", double(clock() - p_begin) / CLOCKS_PER_SEC);
		}

		if (i > 0) {
			MultigridLevelVariables<DataType> & _variables = variables[i];
			const MultigridLevelIndices<Real> & previousLevelIndices = indices[i - 1];
			MultigridLevelVariables<DataType> & previousLevelVariables = variables[i - 1];

			clock_t p_begin;
			if (verbose) p_begin = clock();
			AddProlongation<Real>(previousLevelIndices.prolongationLines, _variables.x, previousLevelVariables.x, detailVerbose);
			if (verbose) printf("Prolongation %.4f \n", double(clock() - p_begin) / CLOCKS_PER_SEC);
		}
	}

	return 1;
}
