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

template< class Real , class Data , class Solver >
void Relaxation
(
	const std::vector< Real > &deepCoefficients , const SparseMatrix< Real , int > &boundaryDeepMatrix , Solver &boundarySolver , const std::vector< int > &boundaryToSupported ,
	const std::vector< SegmentedRasterLine > &segmentedLines ,
	const std::vector< Data > &rhs , std::vector< Data > &x0 , std::vector< Data > &boundaryRHS , std::vector< Data > &boundarySolution , std::vector< Data > &variableBoundaryRHS ,
	int numIterations=2 , bool boundaryFirst=true , bool verbose=false
)
{
	int numBoundaryVariables = boundaryToSupported.size();
	
	ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( unsigned int , size_t i ){ boundaryRHS[i] = rhs[ boundaryToSupported[i] ]; } );

	Miscellany::Timer timer;

	int it_offset = boundaryFirst ? 0 : 1;
	for (int it = 0; it < numIterations; it++){
		if ((it + it_offset) % 2 == 0){//Update Boundary;

			if( verbose ) timer.reset();

			ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( unsigned int , size_t i ){ boundarySolution[i] = x0[ boundaryToSupported[i] ] ; variableBoundaryRHS[i] = boundaryRHS[i]; } );

			boundaryDeepMatrix.Multiply((Data*)&x0[0], (Data*)&variableBoundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);

			if( verbose ) printf("\t Boundary initialization =  %.4f \n" , timer.elapsed() );

			if( verbose ) timer.reset();
			solve(boundarySolver,boundarySolution, variableBoundaryRHS);
			if( verbose ) printf("\t Boundary update =  %.4f \n" , timer.elapsed() );

			ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( unsigned int , size_t i ){ x0[ boundaryToSupported[i] ] = boundarySolution[i]; } );
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


			auto UpdateBlock = [&](int firstLine, int lastLine, int passes){
				for (int r = firstLine; r < lastLine + passes - 1; r++) for (int p = 0; p < passes; p++) if (r - p >= firstLine && r - p < lastLine) UpdateRow(r - p);
			};

			if( verbose ) timer.reset();

			int blockSize = 16;
			int totalLines = segmentedLines.size();
			int numBlocks = ((totalLines - 1) / blockSize) + 1;

			ThreadPool::ParallelFor
				(
					0 , (numBlocks+1)/2 ,
					[&]( unsigned int , size_t _b )
					{
						size_t b = _b * 2;
						int firstLine = (totalLines * b)/ numBlocks;
						int lastLine = (totalLines * (b + 1)) / numBlocks;
						UpdateBlock(firstLine, lastLine, 3);
					}
				);
			ThreadPool::ParallelFor
			(
				0 , numBlocks/2 ,
				[&]( unsigned int , size_t _b )
				{
					size_t b = _b * 2 + 1;
					int firstLine = (totalLines * b) / numBlocks;
					int lastLine = (totalLines * (b + 1)) / numBlocks;
					UpdateBlock(firstLine, lastLine, 3);
				}
			);
			if( verbose ) printf( "\t GS Deep update =  %.4f\n" , timer.elapsed() );
		}
	}
}




template< class Real , class Data , class Solver >
void Relaxation
(
	const std::vector< Real > &deepCoefficients ,
	const SparseMatrix< Real , int > &boundaryDeepMatrix ,
	Solver &boundarySolver ,
	const std::vector< unsigned int > &boundaryToSupported ,
	const std::vector< ThreadTask > &threadTasks ,
	const std::vector< Data > &rhs ,
	std::vector< Data > &x0 ,
	std::vector< Data > &boundaryRHS ,
	std::vector< Data > &boundarySolution ,
	std::vector< Data > &variableBoundaryRHS ,
	unsigned int numIterations=2 ,
	bool boundaryFirst=true ,
	bool verbose=false
	)
{


	int numBoundaryVariables = (int)boundaryToSupported.size();

	ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( unsigned int , size_t i ){ boundaryRHS[i] = rhs[ boundaryToSupported[i] ]; } );

	Miscellany::Timer timer;

	int it_offset = boundaryFirst ? 0 : 1;
	for( unsigned int it=0 ; it<numIterations ; it++ )
	{
		if ((it + it_offset) % 2 == 0) {//Update Boundary;

			if( verbose ) timer.reset();

			ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( unsigned int , size_t i ){ boundarySolution[i] = x0[ boundaryToSupported[i] ] ; variableBoundaryRHS[i] = boundaryRHS[i]; } );

			boundaryDeepMatrix.Multiply((Data*)&x0[0], (Data*)&variableBoundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);

			if( verbose ) printf( "\t Boundary initialization =  %.4f \n" , timer.elapsed() );

			if( verbose ) timer.reset();
			solve(boundarySolver, boundarySolution, variableBoundaryRHS);
			if( verbose ) printf( "\t Boundary update =  %.4f\n" , timer.elapsed() );

			ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( unsigned int , size_t i ){ x0[ boundaryToSupported[i] ] = boundarySolution[i]; } );
		}


		if ((it + it_offset) % 2 == 1) {//Update Interior
			
			if( verbose ) timer.reset();
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

			// Update even blocks
			ThreadPool::ParallelFor
				(
					0 , (threadTasks.size() + 1) / 2 ,
					[&]( unsigned int , size_t t ){ for( int b=0 ; b<threadTasks[2*t+0].blockTasks.size() ; b++ ) UpdateBlock( (int)(2*t+0) , b , 3 ); }
				);
			ThreadPool::ParallelFor
				(
					0 , threadTasks.size() / 2 ,
					[&]( unsigned int , size_t t ){ for( int b=0 ; b<threadTasks[2*t+1].blockTasks.size() ; b++ ) UpdateBlock( (int)(2*t+1) , b , 3 ); }
				);
			if( verbose ) printf( "\t GS Deep update =  %.4f\n" , timer.elapsed() );
		}
	}
}


template< class Real , class Data , class Solver >
void RelaxationAndResidual
(
	const std::vector< Real > &deepCoefficients ,
	const SparseMatrix< Real , int > &boundaryDeepMatrix ,
	Solver &boundarySolver ,
	const std::vector< unsigned int > &boundaryToSupported ,
	const std::vector< ThreadTask > &threadTasks ,
	const std::vector< Data > &rhs ,
	std::vector< Data > &x0 ,
	std::vector< Data > &boundaryRHS ,
	std::vector< Data > &boundaryValue ,
	std::vector< Data > &variableBoundaryRHS ,
	const SparseMatrix< Real , int > & boundaryBoundaryMatrix ,
	std::vector< Data > & residual ,
	unsigned int numIterations=2 ,
	bool verbose=false
)
{
	int numBoundaryVariables = (int)boundaryToSupported.size();

	ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( unsigned int , size_t i ){ boundaryRHS[i] = rhs[ boundaryToSupported[i] ]; } );

	Miscellany::Timer timer;
	for( unsigned int it=0 ; it<numIterations ; it++ )
	{
		if (it % 2 == 0) {//Update Boundary

			if( verbose ) timer.reset();

			ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( unsigned int , size_t i ){ boundaryValue[i] = x0[ boundaryToSupported[i] ] ; variableBoundaryRHS[i] = boundaryRHS[i]; } );

			boundaryDeepMatrix.Multiply((Data*)&x0[0], (Data*)&variableBoundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);
			if( verbose ) printf( "\t Boundary initialization =  %.4f\n" , timer.elapsed() );


			if( verbose ) timer.reset();
			solve(boundarySolver, boundaryValue, variableBoundaryRHS);

			if( verbose ) printf( "\t Boundary update =  %.4f\n" , timer.elapsed() );
			ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( unsigned int , size_t i ){ x0[ boundaryToSupported[i] ] = boundaryValue[i]; } );
		}

		if (it % 2 == 1) {//Update Interior


			if( verbose ) timer.reset();


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
			ThreadPool::ParallelFor
				(
					0 , (threadTasks.size() + 1) / 2 ,
					[&]( unsigned int , size_t t ){ for( int b=0 ; b<threadTasks[2*t+0].blockTasks.size() ; b++ ) UpdateBlock( (int)(2*t+0) , b , 3 ); }
				);
			ThreadPool::ParallelFor
				(
					0 , threadTasks.size() / 2 ,
					[&]( unsigned int , size_t t ){ for( int b=0 ; b<threadTasks[2*t+1].blockTasks.size() ; b++ ) UpdateBlock( (int)(2*t+1) , b , 3 ); }
				);

			if( verbose ) printf( "\t GS Deep update =  %.4f\n" , timer.elapsed() );
		}
	}

	{//Compute boundary residual
		if( verbose ) timer.reset();
		boundaryDeepMatrix.Multiply((Data*)&x0[0], (Data*)&boundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);
		boundaryBoundaryMatrix.Multiply((Data*)&boundaryValue[0], (Data*)&boundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);
		for (int i = 0; i < numBoundaryVariables; i++)residual[ boundaryToSupported[i] ] = boundaryRHS[i];
		if( verbose ) printf( "\t Boundary residual =  %.4f\n" , timer.elapsed() );
	}
}


template< class Real , class Data , class Solver >
void RelaxationAndResidual
(
	const std::vector< Real > &deepCoefficients , const SparseMatrix< Real , int > &boundaryDeepMatrix , Solver &boundarySolver , const std::vector< int > &boundaryToSupported ,
	const std::vector< SegmentedRasterLine > &segmentedLines ,
	const std::vector< Data > &rhs , std::vector< Data > &x0 , std::vector< Data > &boundaryRHS , std::vector< Data > &boundaryValue , std::vector< Data > &variableBoundaryRHS ,
	const SparseMatrix< Real , int > &boundaryBoundaryMatrix , std::vector< Data > &residual ,
	int numIterations=2 , bool verbose=false
)
{
	int numBoundaryVariables = boundaryToSupported.size();

	ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( unsigned int , size_t i ){ boundaryRHS[i] = rhs[ boundaryToSupported[i] ]; } );

	Miscellany::Timer timer;
	for (int it = 0; it < numIterations; it++) {
		if (it % 2 == 0) {//Update Boundary

			if( verbose ) timer.reset();

			ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( unsigned int , size_t i ){ boundaryValue[i] = x0[ boundaryToSupported[i] ] ; variableBoundaryRHS[i] = boundaryRHS[i] ; } );

			boundaryDeepMatrix.Multiply((Data*)&x0[0], (Data*)&variableBoundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);
			if( verbose ) printf( "\t Boundary initialization =  %.4f\n" , timer.elapsed() );

			if( verbose ) timer.reset();
			solve(boundarySolver, boundaryValue, variableBoundaryRHS);

			if( verbose ) printf( "\t Boundary update =  %.4f\n" , timer.elapsed() );
			ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( unsigned int , size_t i ){ x0[ boundaryToSupported[i] ] = boundaryValue[i]; } );
		}

		if (it % 2 == 1) {//Update Interior

										//static int GSIterations = 3;
										//for (int r = 0; r < rasterLines.size() + GSIterations; r++)
										//	for (int j = 0; j < GSIterations; j++)
										//		if (r - j >= 0 && r - j < rasterLines.size()) UpdateRow(r - j);

			if( verbose ) timer.reset();


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

			auto UpdateHeaderBlock = [&](int firstLine, int lastLine, int passes, const bool firstResidual,const bool prevResidual){
				for (int r = firstLine; r < firstLine + passes + 1; r++){
					for (int p = 0; p < passes; p++) if (r - p >= firstLine && r - p < lastLine) UpdateRow(r - p);
				}
				if(firstResidual)UpdateResidual(firstLine);
				if(prevResidual)UpdateResidual(firstLine - 1);
			};

			auto UpdateMiddleBlock = [&](int firstLine, int lastLine, int passes) {
				for (int r = firstLine + passes + 1; r < lastLine + passes - 1; r++) {
					for (int p = 0; p < passes; p++) if (r - p > firstLine && r - p < lastLine) UpdateRow(r - p);
					if (r - passes> firstLine && r - passes < lastLine) UpdateResidual(r - passes);
				}
			};

			auto UpdateSuffixBlock = [&](int lastLine,const bool lastResidual, const bool nextResidual) {
				if (lastResidual)UpdateResidual(lastLine - 1);
				if (nextResidual)UpdateResidual(lastLine);
			};


			int blockSize = 16;
			int totalLines = segmentedLines.size();
			int numBlocks = ((totalLines - 1) / blockSize) + 1;

			//Even blocks
			ThreadPool::ParallelFor
				(
					0 , (numBlocks+1)/2 ,
					[&]( unsigned int , size_t _b )
					{
						size_t b = 2*_b;
						int firstLine = (totalLines * b)/ numBlocks;
						int lastLine = (totalLines * (b + 1)) / numBlocks;
						UpdateHeaderBlock(firstLine, lastLine, 3, firstLine == 0, false);
						UpdateMiddleBlock(firstLine, lastLine, 3);
						UpdateSuffixBlock(lastLine, lastLine == totalLines, false);
					}
				);

			//Odd blocks
			ThreadPool::ParallelFor
			(
				0 , numBlocks/2 ,
				[&]( unsigned int , size_t _b )
				{
					size_t b = 2*_b+1;
					int firstLine = (totalLines * b) / numBlocks;
					int lastLine = (totalLines * (b + 1)) / numBlocks;
					UpdateHeaderBlock(firstLine, lastLine, 3, true, firstLine != 0);
					UpdateMiddleBlock(firstLine, lastLine, 3);
					UpdateSuffixBlock(lastLine, true, lastLine != totalLines);
				}
			);

			if( verbose ) printf( "\t GS Deep update =  %.4f\n" , timer.elapsed() );
		}
	}

	{//Compute boundary residual
		if( verbose ) timer.reset();
		boundaryDeepMatrix.Multiply((Data*)&x0[0], (Data*)&boundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);
		boundaryBoundaryMatrix.Multiply((Data*)&boundaryValue[0], (Data*)&boundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);
		for (int i = 0; i < numBoundaryVariables; i++)residual[ boundaryToSupported[i] ] = boundaryRHS[i];
		if( verbose ) printf( "\t Boundary residual =  %.4f\n" , timer.elapsed() );
	}
}

template< class Real , class Data , class DataReal=Real >
void MultiplyBySystemMatrix( const SystemCoefficients< Real > &systemCoefficients , const std::vector< int > &boundaryToSupported , const std::vector< RasterLine > &rasterLines , const std::vector< Data > &in , std::vector< Data > &out , bool verbose=false )
{
	Miscellany::Timer timer;

	int numBoundaryVariables = boundaryToSupported.size();

	std::vector<Data> inBoundaryValues;
	inBoundaryValues.resize(numBoundaryVariables);
	for (int i = 0; i < numBoundaryVariables; i++) inBoundaryValues[i] = in[ boundaryToSupported[i] ];

	std::vector<Data> outBoundaryValues;
	outBoundaryValues.resize(numBoundaryVariables);

	if( verbose ) timer.reset();
	systemCoefficients.boundaryBoundaryMatrix.Multiply(&inBoundaryValues[0], &outBoundaryValues[0]);
	systemCoefficients.boundaryDeepMatrix.Multiply(&in[0], &outBoundaryValues[0], MULTIPLY_ADD);
	if( verbose ) printf( "\t Multiply boundary =  %.4f\n" , timer.elapsed() );

	ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( unsigned int , size_t i ){ out[ boundaryToSupported[i] ] = outBoundaryValues[i]; } );

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
				(
					_inPrevious[i-1] * (DataReal)_deepCoefficients[0] +
					_inPrevious[i+0] * (DataReal)_deepCoefficients[1] +
					_inPrevious[i+1] * (DataReal)_deepCoefficients[2] +
					_inCurrent [i-1] * (DataReal)_deepCoefficients[3] +
					_inCurrent [i+0] +
					_inCurrent [i+1] * (DataReal)_deepCoefficients[5] +
					_inNext    [i-1] * (DataReal)_deepCoefficients[6] +
					_inNext    [i+0] * (DataReal)_deepCoefficients[7] +
					_inNext    [i+1] * (DataReal)_deepCoefficients[8]
				) * (DataReal)_deepCoefficients[9];
		}
	};

	if( verbose ) timer.reset();

	unsigned int threads = ThreadPool::NumThreads();
	std::vector<int> lineRange(threads + 1);
	int blockSize = rasterLines.size() / threads;
	for (int t = 0; t < threads; t++) lineRange[t] = t*blockSize;
	lineRange[threads] = rasterLines.size();
	ThreadPool::ParallelFor
		(
			0 , threads ,
			[&]( unsigned int , size_t t )
			{
				int firstLine = lineRange[t];
				int lastLine = lineRange[t + 1];
				for( int r=firstLine ; r<lastLine ; r++ ) UpdateRow(r);
			}
		);

	if( verbose ) printf( "\t Multiply deep %.4f\n" , timer.elapsed() );
}

template< class Real , class Data , class DataReal=Real >
void MultiplyBySystemMatrix_NoReciprocals
(
	const SystemCoefficients< Real > &systemCoefficients ,
	const typename GridAtlas<>::IndexConverter & indexConverter ,
	const std::vector< RasterLine >& rasterLines ,
	const std::vector< Data >& in ,
	std::vector< Data > & out ,
	bool verbose=false
)
{
	Miscellany::Timer timer;

	unsigned int numBoundaryVariables = (unsigned int)indexConverter.numBoundary();

	// Pull out the boundary values from the input array
	std::vector< Data > outBoundaryValues( numBoundaryVariables );
	std::vector< Data >  inBoundaryValues( numBoundaryVariables );
	for( unsigned int i=0 ; i<numBoundaryVariables ; i++ ) inBoundaryValues[i] = in[ indexConverter.boundaryToSupported(i) ];

	timer.reset();
	//  Perform the boundary -> boundary multiplication
	systemCoefficients.boundaryBoundaryMatrix.template Multiply< Data , DataReal >( &inBoundaryValues[0] , &outBoundaryValues[0] );
	// Perform the interior -> boundary multiplication
	systemCoefficients.boundaryDeepMatrix.template Multiply< Data , DataReal >( &in[0] , &outBoundaryValues[0] , MULTIPLY_ADD );
	if( verbose ) printf( "\tMultiply boundary = %.4f\n" , timer.elapsed() );

	// Write the boundary values back into the output array
	ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( unsigned int , size_t i ){ out[ indexConverter.boundaryToSupported((unsigned int)i) ] = outBoundaryValues[i]; } );

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
			_out[i] =
				(
					_inPrevious[i-1] * (DataReal)_deepCoefficients[0] +
					_inPrevious[i+0] * (DataReal)_deepCoefficients[1] +
					_inPrevious[i+1] * (DataReal)_deepCoefficients[2] +
					_inCurrent [i-1] * (DataReal)_deepCoefficients[3] +
					_inCurrent [i+0] * (DataReal)_deepCoefficients[4] +
					_inCurrent [i+1] * (DataReal)_deepCoefficients[5] +
					_inNext    [i-1] * (DataReal)_deepCoefficients[6] +
					_inNext    [i+0] * (DataReal)_deepCoefficients[7] +
					_inNext    [i+1] * (DataReal)_deepCoefficients[8]
				);
		}
	};



	timer.reset();

	unsigned int threads = ThreadPool::NumThreads();
	std::vector< int > lineRange( threads+1 );
	int blockSize = (int)rasterLines.size() / threads;
	for( unsigned int t=0 ; t<threads ; t++ ) lineRange[t] = t*blockSize;
	lineRange[threads] = (int)rasterLines.size();
	ThreadPool::ParallelFor
		(
			0 , threads ,
			[&]( unsigned int , size_t t )
			{
				int firstLine = lineRange[t];
				int lastLine = lineRange[t+1];
				for( int r=firstLine ; r<lastLine ; r++ ) UpdateRow(r);
			}
		);

	if( verbose ) printf( "\tMultiply deep = %.4f \n" , timer.elapsed() );
}


template< class Real , class Data >
void ComputeSystemResidual
(
	const std::vector< Real > &deepCoefficients , const SparseMatrix< Real , int > &boundaryDeepMatrix , const SparseMatrix< Real , int > &boundaryBoundaryMatrix ,
	const std::vector< int > &boundaryToSupported , const std::vector< RasterLine > &rasterLines ,
	const std::vector< Data > &boundaryRHS , const std::vector< Data > &boundaryValues , const std::vector< Data > &in , const std::vector< Data > &rhs , std::vector< Data > &out ,
	bool verbose=false
)
{
	Miscellany::Timer timer;
	//Boundary residual
	{
		int numBoundaryVariables = boundaryToSupported.size();
		if( verbose ) timer.reset();
		boundaryDeepMatrix.Multiply((Data*)&in[0], (Data*)&boundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);
		boundaryBoundaryMatrix.Multiply((Data*)&boundaryValues[0], (Data*)&boundaryRHS[0], MULTIPLY_ADD | MULTIPLY_NEGATE);
		for (int i = 0; i < numBoundaryVariables; i++) out[ boundaryToSupported[i] ] = boundaryRHS[i];
		if( verbose ) printf( "\t Residual Boundary =  %.4f \n" , timer.elapsed() );
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


	if( verbose ) timer.reset();

	unsigned int threads = ThreadPool::NumThreads();
	ThreadPool::ParallelFor
		(
			0 , threads ,
			[&]( unsigned int , size_t t )
			{
				int firstLine = (rasterLines.size() * t) / threads;
				int lastLine = (rasterLines.size() * (t+1)) / threads;
				for( int r=firstLine ; r<lastLine ; r++ ) UpdateResidual(r);
			}
		);

	if( verbose ) printf( "\t Residual Deep %.4f\n" , timer.elapsed() );
}

template< class Real , class Data >
void MultiplyByRestriction
(
	const SparseMatrix< Real , int > & __boundaryRestrictionMatrix ,
	const std::vector< unsigned int > &boundaryToSupported ,
	std::vector< Data > &boundaryValue ,
	const std::vector< RasterLine > &restrictionLines ,
	const std::vector< Data > &in ,
	std::vector< Data > &out ,
	bool verbose=false
)
{
	Miscellany::Timer timer;

	unsigned int numBoundaryVariables = (unsigned int)boundaryToSupported.size();

	if( verbose ) timer.reset();
	__boundaryRestrictionMatrix.Multiply((Data*)&in[0], (Data*)&boundaryValue[0]);
	if( verbose ) printf( "\t Restriction boundary  %.4f \n" , timer.elapsed() );
	
	ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( unsigned int , size_t i ){ out[ boundaryToSupported[i] ] = boundaryValue[i]; } );


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

	if( verbose ) timer.reset();
	unsigned int threads = ThreadPool::NumThreads();
	std::vector<int> lineRange(threads + 1);
	int blockSize = (int)restrictionLines.size() / threads;
	for( unsigned int t=0 ; t<threads ; t++ ) lineRange[t] = t*blockSize;
	lineRange[threads] = (int)restrictionLines.size();
	ThreadPool::ParallelFor
		(
			0 , threads ,
			[&]( unsigned int , size_t t )
			{
				int firstLine = lineRange[t];
				int lastLine = lineRange[t+1];
				for( int r=firstLine ; r<lastLine ; r++ ) UpdateRow(r);
			}
		);

	if( verbose ) printf( "\t Restriction deep %.4f \n" , timer.elapsed() );
}

template< class Real , class Data >
void MultiplyByProlongation
(
	const std::vector< ProlongationLine > &prolongationLines ,
	const std::vector< Data > &in ,
	std::vector< Data > &out ,
	bool verbose=false
)
{
	Miscellany::Timer timer;

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

	if( verbose ) timer.reset();
	unsigned int threads = ThreadPool::NumThreads();
	std::vector<int> lineRange(threads + 1);
	int blockSize = prolongationLines.size() / threads;
	for (int t = 0; t < threads; t++) lineRange[t] = t*blockSize;
	lineRange[threads] = prolongationLines.size();
	ThreadPool::ParallelFor
		(
			0 , threads ,
			[&]( unsigned int , size_t t )
			{
				int firstLine = lineRange[t];
				int lastLine = lineRange[t + 1];
				for( int r=firstLine ; r<lastLine ; r++ ) UpdateRow(r);
			}
		);
	if( verbose ) printf( "\t Prolongation deep %.4f\n" , timer.elapsed() ); 
}


template< class Real , class Data >
void AddProlongation( const std::vector< ProlongationLine > &prolongationLines , const std::vector< Data > &in , std::vector< Data > &out , bool verbose=false )
{
	Miscellany::Timer timer;

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

	if( verbose ) timer.reset();
	unsigned int threads = ThreadPool::NumThreads();
	std::vector<int> lineRange(threads + 1);
	int blockSize = (int)prolongationLines.size() / threads;
	for( unsigned int t=0 ; t<threads ; t++ ) lineRange[t] = t*blockSize;
	lineRange[threads] = (int)prolongationLines.size();
	ThreadPool::ParallelFor
		(
			0 , threads ,
			[&]( unsigned int , size_t t )
			{
				int firstLine = lineRange[t];
				int lastLine = lineRange[t+1];
				for( int r=firstLine ; r<lastLine ; r++ ) UpdateRow(r);
			}
		);

	if( verbose ) printf( "\t Prolongation deep %.4f \n" , timer.elapsed() );
}

template< class Real , unsigned int Channels , class Data >
void CellStiffnessToTexelStiffness
(
	const std::vector< Real >& cellSharpenningMask ,
	const std::vector< InteriorTexelToCellLine >& interiorTexelToCellLines ,
	const std::vector< Data >& interiorTexelToCellCoeffs ,
	SparseMatrix< Real , int > boundaryCellBasedStiffnessRHSMatrix[Channels] ,
	std::vector< Real > boundaryTexelValues[Channels] ,
	const typename GridAtlas<>::IndexConverter & indexConverter ,
	std::vector< Data >& texelModulatedStiffness ,
	bool verbose=false
)
{
	//Update Boundary Texels
	Miscellany::Timer timer;
	ThreadPool::ParallelFor( 0 , Channels , [&]( unsigned int , size_t c ){ boundaryCellBasedStiffnessRHSMatrix[c].Multiply( &cellSharpenningMask[0] , &boundaryTexelValues[c][0] ); } );
	ThreadPool::ParallelFor( 0 , indexConverter.numBoundary() , [&]( unsigned int , size_t i ){ for( int c=0 ; c<Channels ; c++ ) texelModulatedStiffness[ indexConverter.boundaryToSupported(i) ][c] = boundaryTexelValues[c][i]; } );
	if( verbose ) printf( "\tBoundary updated: %.3f(s) \n" , timer.elapsed() );
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

	if( verbose ) timer.reset();
	unsigned int threads = ThreadPool::NumThreads();
	ThreadPool::ParallelFor
		(
			0 , threads ,
			[&]( unsigned int , size_t t )
			{
				int firstLine = (int)(interiorTexelToCellLines.size() * t) / threads;
				int lastLine = (int)(interiorTexelToCellLines.size() * (t+1)) / threads;
				for( int r=firstLine ; r<lastLine ; r++ ) UpdateRow(r);
			}
		);

	if( verbose ) printf( "\t Interior update %.4f\n" , timer.elapsed() );
}
template< class Real >
void CellStiffnessToTexelStiffness
(
	const std::vector< Real >& cellSharpenningMask ,
	const std::vector< InteriorTexelToCellLine >& interiorTexelToCellLines ,
	const std::vector< Real >& interiorTexelToCellCoeffs ,
	SparseMatrix< Real , int > boundaryCellBasedStiffnessRHSMatrix ,
	std::vector< Real > boundaryTexelValues ,
	const typename GridAtlas<>::IndexConverter & indexConverter ,
	std::vector< Real >& texelModulatedStiffness ,
	bool verbose=false
)
{
	//Update Boundary Texels
	Miscellany::Timer timer;
	boundaryCellBasedStiffnessRHSMatrix.Multiply( &cellSharpenningMask[0] , &boundaryTexelValues[0] );
	ThreadPool::ParallelFor( 0 , boundaryToSupported.size() , [&]( unsigned int , size_t i ){ texelModulatedStiffness[ boundaryToSupported[i] ] = boundaryTexelValues[i]; } );
	if( verbose ) printf( "\tBoundary updated: %.3f(s) \n" , timer.elapsed() );
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

	if( verbose ) timer.reset();
	unsigned int threads = ThreadPool::NumThreads();
	ThreadPool::ParallelFor
		(
			0 , threads ,
			[&]( unsigned int , size_t t )
			{
				int firstLine = (int)(interiorTexelToCellLines.size() * t) / threads;
				int lastLine = (int)(interiorTexelToCellLines.size() * (t+1) ) / threads;
				for( int r=firstLine ; r<lastLine ; r++ ) UpdateRow(r);
			}
		);
	if( verbose ) printf( "\t Interior update %.4f\n" , timer.elapsed() );
}

template< class Real , class DataType , class DirectSolver >
void VCycle( std::vector< MultigridLevelVariables< DataType > > &variables , const std::vector< SystemCoefficients< Real > > &coefficients , const std::vector< MultigridLevelIndices< Real > > &indices , VCycleSolvers< DirectSolver > &vCycleSolvers , bool verbose , bool detailVerbose )
{
	int levels = (int)variables.size();

	//Set x0[i = 1 : levels -1] = 0
	Miscellany::Timer timer;
	for (int i = 1; i < levels; i++){
		memset(&variables[i].x[0], 0, variables[i].x.size() * sizeof(DataType));
	}
	if( verbose ) printf("Zero arrays %.4f \n" , timer.elapsed() );

	//Reduction phase
	for( int i=0 ; i<levels-1 ; i++ )
	{
		const SystemCoefficients< Real > & _coefficients = coefficients[i];
		const MultigridLevelIndices<Real> & _indices = indices[i];
		MultigridLevelVariables<DataType> & _variables = variables[i];

		const MultigridLevelIndices<Real> & nextLevelIndices = indices[i + 1];
		MultigridLevelVariables<DataType> & nextLevelVariables = variables[i + 1];

		if( verbose ) printf( "Level %d\n" , i );

		Miscellany::Timer tmr;
		RelaxationAndResidual( _coefficients.deepCoefficients , _coefficients.boundaryDeepMatrix , vCycleSolvers.boundary[i] , _indices.boundaryToSupported , _indices.threadTasks , _variables.rhs , _variables.x , _variables.boundary_rhs , _variables.boundary_value , _variables.variable_boundary_value , _coefficients.boundaryBoundaryMatrix , _variables.residual , 2 , detailVerbose );
		if( verbose ) printf("Relaxation  + Residual %.4f\n" , tmr.elapsed() );

		if( verbose ) tmr.reset();
		MultiplyByRestriction( _indices.boundaryRestriction , nextLevelIndices.boundaryToSupported , nextLevelVariables.boundary_value , nextLevelIndices.restrictionLines, _variables.residual, nextLevelVariables.rhs , detailVerbose );
		if( verbose ) printf( "Restriction %.4f\n" , tmr.elapsed() );
	}

	//Prolongation phase
	for( int i=levels-1 ; i>=0 ; i-- )
	{
		if( verbose ) printf( "Level %d\n" , i );
		if( i<levels-1 )
		{
			Miscellany::Timer tmr;

			const SystemCoefficients< Real > & _coefficients = coefficients[i];
			const MultigridLevelIndices<Real> & _indices = indices[i];
			MultigridLevelVariables<DataType> & _variables = variables[i];

			if( verbose ) tmr.reset();
			Relaxation( _coefficients.deepCoefficients , _coefficients.boundaryDeepMatrix , vCycleSolvers.boundary[i] , _indices.boundaryToSupported , _indices.threadTasks , _variables.rhs , _variables.x , _variables.boundary_rhs , _variables.boundary_value , _variables.variable_boundary_value , 2 , true , detailVerbose );
			if( verbose ) printf( "Gauss Seidel %.4f\n" , tmr.elapsed() );
		}
		else if( i==levels-1 )
		{
			MultigridLevelVariables<DataType> & _variables = variables[i];
			Miscellany::Timer tmr;
			if( verbose ) tmr.reset();
			solve( vCycleSolvers.coarse , _variables.x , _variables.rhs );
			if( verbose ) printf( "Direct solver %.4f\n" , tmr.elapsed() );
		}

		if( i>0 )
		{
			MultigridLevelVariables<DataType> & _variables = variables[i];
			const MultigridLevelIndices<Real> & previousLevelIndices = indices[i - 1];
			MultigridLevelVariables<DataType> & previousLevelVariables = variables[i - 1];

			Miscellany::Timer tmr;
			AddProlongation<Real>(previousLevelIndices.prolongationLines, _variables.x, previousLevelVariables.x, detailVerbose);
			if( verbose ) printf( "Prolongation %.4f\n" , tmr.elapsed() );
		}
	}
}
