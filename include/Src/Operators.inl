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

///////////////////////////////
// MassAndStiffnessOperators //
///////////////////////////////

template< typename MatrixReal >
template< typename Data >
void MassAndStiffnessOperators< MatrixReal >::operator()( double mWeight , double sWeight , const std::vector< Data > &in , std::vector< Data > &out , bool add ) const
{
	return operator()( mWeight , sWeight , GetPointer( in ) , GetPointer( out ) , add );
}

template< typename MatrixReal >
template< typename Data >
std::vector< Data > MassAndStiffnessOperators< MatrixReal >::operator()( double mWeight , double sWeight , const std::vector< Data > &in ) const
{
	std::vector< Data > out( rows() );
	operator()( mWeight , sWeight , in , out );
	return out;
}

template< typename MatrixReal >
template< typename Data >
void MassAndStiffnessOperators< MatrixReal >::operator()( double mWeight , double sWeight , ConstPointer( Data ) in , Pointer( Data ) out , bool add ) const
{
	if( add )
	{
		if     ( mWeight==1.0 && sWeight==0. ) return _evaluate< true , true , false >( mWeight , sWeight , in , out );
		else if( mWeight==0.0 && sWeight==1. ) return _evaluate< true ,false , true >( mWeight , sWeight , in , out );
		else                                   return _evaluate< true ,false , false >( mWeight , sWeight , in , out );
	}
	else
	{
		if     ( mWeight==1.0 && sWeight==0. ) return _evaluate< false , true , false >( mWeight , sWeight , in , out );
		else if( mWeight==0.0 && sWeight==1. ) return _evaluate< false , false , true >( mWeight , sWeight , in , out );
		else                                   return _evaluate< false , false , false >( mWeight , sWeight , in , out );
	}
}

template< typename MatrixReal >
template< typename Data >
void MassAndStiffnessOperators< MatrixReal >::mass( const std::vector< Data > &in , std::vector< Data > &out , bool add ) const
{
	if( add ) return _evaluate< true  , true , false >( 1 , 0 , GetPointer( in ) , GetPointer( out ) );
	else      return _evaluate< false , true , false >( 1 , 0 , GetPointer( in ) , GetPointer( out ) );
}

template< typename MatrixReal >
template< typename Data >
void MassAndStiffnessOperators< MatrixReal >::stiffness( const std::vector< Data > &in , std::vector< Data > &out , bool add ) const
{
	if( add ) return _evaluate< true  , false , true >( 0 , 1 , GetPointer( in ) , GetPointer( out ) );
	else      return _evaluate< false , false , true >( 0 , 1 , GetPointer( in ) , GetPointer( out ) );
}

template< typename MatrixReal >
template< typename Data >
std::vector< Data > MassAndStiffnessOperators< MatrixReal >::mass( const std::vector< Data > &in ) const
{
	std::vector< Data > out( rows() );
	_evaluate< false , true , false >( 1 , 0 , in , out );
	return out;
}

template< typename MatrixReal >
template< typename Data >
std::vector< Data > MassAndStiffnessOperators< MatrixReal >::stiffness( const std::vector< Data > &in ) const
{
	std::vector< Data > out( rows() );
	_evaluate< false , false , true >( 0 , 1 , in , out );
	return out;
}

template< typename MatrixReal >
template< typename Data >
void MassAndStiffnessOperators< MatrixReal >::mass( ConstPointer( Data ) in , Pointer( Data ) out , bool add ) const
{
	if( add ) return _evaluate< true  , true , false >( 1 , 0 , in , out );
	else      return _evaluate< false , true , false >( 1 , 0 , in , out );
}

template< typename MatrixReal >
template< typename Data >
void MassAndStiffnessOperators< MatrixReal >::stiffness( ConstPointer( Data ) in , Pointer( Data ) out , bool add ) const
{
	if( add ) return _evaluate< true  , false , true >( 0 , 1 , in , out );
	else      return _evaluate< false , false , true >( 0 , 1 , in , out );
}

template< typename MatrixReal >
template< typename OutReal >
Eigen::SparseMatrix< OutReal > MassAndStiffnessOperators< MatrixReal >::operator()( double mWeight , double sWeight ) const
{
	if     ( mWeight==1. && sWeight==0. ) return _matrix< true , false , OutReal >( 1 , 0 );
	else if( mWeight==0. && sWeight==1. ) return _matrix< false , true , OutReal >( 0 , 1 );
	else                                  return _matrix< false , false , OutReal >( mWeight , sWeight );
}

template< typename MatrixReal >
template< typename OutReal >
Eigen::SparseMatrix< OutReal > MassAndStiffnessOperators< MatrixReal >::mass( void ) const{ return _matrix< true , false , OutReal >( 1 , 0 ); }

template< typename MatrixReal >
template< typename OutReal >
Eigen::SparseMatrix< OutReal > MassAndStiffnessOperators< MatrixReal >::stiffness( void ) const{ return _matrix< false , true , OutReal >( 0 , 1 ); }


template< typename MatrixReal >
template< bool Add , bool Mass , bool Stiffness , typename Data >
void MassAndStiffnessOperators< MatrixReal >::_evaluate( double mWeight , double sWeight , ConstPointer( Data ) in , Pointer( Data ) out ) const
{
	static_assert( !( Mass && Stiffness ) , "[ERROR] Mass and Stiffness can't both be enabled" );

	// Process boundary
	{
		unsigned int numBoundaryVariables = static_cast< unsigned int >( indexConverter.numBoundary() );

		// Pull out the boundary values from the input array
		std::vector< Data > outBoundaryValues( numBoundaryVariables );
		std::vector< Data >  inBoundaryValues( numBoundaryVariables );
		for( unsigned int i=0 ; i<numBoundaryVariables ; i++ ) inBoundaryValues[i] = in[ static_cast< unsigned int >( indexConverter.boundaryToCombined( AtlasBoundaryTexelIndex(i) ) ) ];

		if constexpr( Mass )
		{
			// Perform the boundary -> boundary multiplication
			massCoefficients.boundaryBoundaryMatrix.template Multiply< Data , MatrixReal >( &inBoundaryValues[0] , &outBoundaryValues[0] , Add ? MULTIPLY_ADD : 0 );
			// Perform the interior -> boundary multiplication
			massCoefficients.boundaryDeepMatrix.template Multiply< Data , MatrixReal >( &in[0] , &outBoundaryValues[0] , MULTIPLY_ADD );
		}
		else if constexpr( Stiffness )
		{
			// Perform the boundary -> boundary multiplication
			stiffnessCoefficients.boundaryBoundaryMatrix.template Multiply< Data , MatrixReal >( &inBoundaryValues[0] , &outBoundaryValues[0] , Add ? MULTIPLY_ADD : 0 );
			// Perform the interior -> boundary multiplication
			stiffnessCoefficients.boundaryDeepMatrix.template Multiply< Data , MatrixReal >( &in[0] , &outBoundaryValues[0] , MULTIPLY_ADD );
		}
		else
		{
			std::vector< Data > _outBoundaryValues( numBoundaryVariables );

			// Perform the boundary -> boundary multiplication
			massCoefficients.boundaryBoundaryMatrix.template Multiply< Data , MatrixReal >( &inBoundaryValues[0] , &_outBoundaryValues[0] );
			// Perform the interior -> boundary multiplication
			massCoefficients.boundaryDeepMatrix.template Multiply< Data , MatrixReal >( &in[0] , &_outBoundaryValues[0] , MULTIPLY_ADD );
			if( Add ) ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( size_t i ){ outBoundaryValues[i] += _outBoundaryValues[i] * mWeight; } );
			else      ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( size_t i ){ outBoundaryValues[i]  = _outBoundaryValues[i] * mWeight; } );

			//  Perform the boundary -> boundary multiplication
			stiffnessCoefficients.boundaryBoundaryMatrix.template Multiply< Data , MatrixReal >( &inBoundaryValues[0] , &_outBoundaryValues[0] );
			// Perform the interior -> boundary multiplication
			stiffnessCoefficients.boundaryDeepMatrix.template Multiply< Data , MatrixReal >( &in[0] , &_outBoundaryValues[0] , MULTIPLY_ADD );
			ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( size_t i ){ outBoundaryValues[i] += _outBoundaryValues[i] * sWeight; } );
		}

		// Write the boundary values back into the output array
		if( Add ) ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( unsigned int , size_t i ){ out[ static_cast< unsigned int >(indexConverter.boundaryToCombined( AtlasBoundaryTexelIndex((unsigned int)i) ) ) ] += outBoundaryValues[i]; } );
		else      ThreadPool::ParallelFor( 0 , numBoundaryVariables , [&]( unsigned int , size_t i ){ out[ static_cast< unsigned int >(indexConverter.boundaryToCombined( AtlasBoundaryTexelIndex((unsigned int)i) ) ) ] = outBoundaryValues[i]; } );
	}

	// Process interior
	{
		// Perform the interior -> interior multiplication
		auto UpdateRow = [&]( int r )
			{

				Data* _out = out + rasterLines[r].lineStartIndex;
				const Data* _inCurr = in + rasterLines[r].lineStartIndex;
				const Data* _inPrev = in + rasterLines[r].prevLineIndex;
				const Data* _inNext = in + rasterLines[r].nextLineIndex;

				int lineLength = ( rasterLines[r].lineEndIndex - rasterLines[r].lineStartIndex + 1 );
				size_t lineDeepStart = static_cast< size_t >( rasterLines[r].coeffStartIndex );
				const MatrixReal * _deepM = &massCoefficients.deepCoefficients[ 10 * lineDeepStart ];
				const MatrixReal * _deepS = &stiffnessCoefficients.deepCoefficients[ 10 * lineDeepStart ];

				for( int i=0 ; i<lineLength ; _deepM+=10 , _deepS+=10 , i++ )
				{
					if constexpr( Mass )
					{
						if constexpr( Add )
							_out[i] +=
							(
								_inPrev[i-1] * _deepM[0] +
								_inPrev[i+0] * _deepM[1] +
								_inPrev[i+1] * _deepM[2] +
								_inCurr[i-1] * _deepM[3] +
								_inCurr[i+0] * _deepM[4] +
								_inCurr[i+1] * _deepM[5] +
								_inNext[i-1] * _deepM[6] +
								_inNext[i+0] * _deepM[7] +
								_inNext[i+1] * _deepM[8]
								);
						else
							_out[i] =
							(
								_inPrev[i-1] * _deepM[0] +
								_inPrev[i+0] * _deepM[1] +
								_inPrev[i+1] * _deepM[2] +
								_inCurr[i-1] * _deepM[3] +
								_inCurr[i+0] * _deepM[4] +
								_inCurr[i+1] * _deepM[5] +
								_inNext[i-1] * _deepM[6] +
								_inNext[i+0] * _deepM[7] +
								_inNext[i+1] * _deepM[8]
								);
					}
					else if constexpr( Stiffness )
					{
						if constexpr( Add )
							_out[i] +=
							(
								_inPrev[i-1] * _deepS[0] +
								_inPrev[i+0] * _deepS[1] +
								_inPrev[i+1] * _deepS[2] +
								_inCurr[i-1] * _deepS[3] +
								_inCurr[i+0] * _deepS[4] +
								_inCurr[i+1] * _deepS[5] +
								_inNext[i-1] * _deepS[6] +
								_inNext[i+0] * _deepS[7] +
								_inNext[i+1] * _deepS[8]
								);
						else
							_out[i] =
							(
								_inPrev[i-1] * _deepS[0] +
								_inPrev[i+0] * _deepS[1] +
								_inPrev[i+1] * _deepS[2] +
								_inCurr[i-1] * _deepS[3] +
								_inCurr[i+0] * _deepS[4] +
								_inCurr[i+1] * _deepS[5] +
								_inNext[i-1] * _deepS[6] +
								_inNext[i+0] * _deepS[7] +
								_inNext[i+1] * _deepS[8]
								);

					}
					else
					{
						if constexpr( Add )
							_out[i] =
							(
								_inPrev[i-1] * ( _deepM[0]*mWeight + _deepS[0]*sWeight ) +
								_inPrev[i+0] * ( _deepM[1]*mWeight + _deepS[1]*sWeight ) +
								_inPrev[i+1] * ( _deepM[2]*mWeight + _deepS[2]*sWeight ) +
								_inCurr[i-1] * ( _deepM[3]*mWeight + _deepS[3]*sWeight ) +
								_inCurr[i+0] * ( _deepM[4]*mWeight + _deepS[4]*sWeight ) +
								_inCurr[i+1] * ( _deepM[5]*mWeight + _deepS[5]*sWeight ) +
								_inNext[i-1] * ( _deepM[6]*mWeight + _deepS[6]*sWeight ) +
								_inNext[i+0] * ( _deepM[7]*mWeight + _deepS[7]*sWeight ) +
								_inNext[i+1] * ( _deepM[8]*mWeight + _deepS[8]*sWeight )
								);
						else
							_out[i] +=
							(
								_inPrev[i-1] * ( _deepM[0]*mWeight + _deepS[0]*sWeight ) +
								_inPrev[i+0] * ( _deepM[1]*mWeight + _deepS[1]*sWeight ) +
								_inPrev[i+1] * ( _deepM[2]*mWeight + _deepS[2]*sWeight ) +
								_inCurr[i-1] * ( _deepM[3]*mWeight + _deepS[3]*sWeight ) +
								_inCurr[i+0] * ( _deepM[4]*mWeight + _deepS[4]*sWeight ) +
								_inCurr[i+1] * ( _deepM[5]*mWeight + _deepS[5]*sWeight ) +
								_inNext[i-1] * ( _deepM[6]*mWeight + _deepS[6]*sWeight ) +
								_inNext[i+0] * ( _deepM[7]*mWeight + _deepS[7]*sWeight ) +
								_inNext[i+1] * ( _deepM[8]*mWeight + _deepS[8]*sWeight )
								);
					}
				}
			};

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
	}
}

template< typename MatrixReal >
template< bool Mass , bool Stiffness , typename OutReal >
Eigen::SparseMatrix< OutReal > MassAndStiffnessOperators< MatrixReal >::_matrix( double mWeight , double sWeight ) const
{
	static_assert( !( Mass && Stiffness ) , "[ERROR] Mass and Stiffness can't both be enabled" );

	std::vector< Eigen::Triplet< MatrixReal > > triplets;

	size_t entries = 0;
	{
		for( unsigned int r=0 ; r<rasterLines.size() ; r++ )
		{
			const RasterLine & line = rasterLines[r];
			entries += ( line.lineEndIndex - line.lineStartIndex + 1 ) * 9;
		}

		entries += massCoefficients.boundaryDeepMatrix.Entries();
		entries += massCoefficients.boundaryBoundaryMatrix.Entries();
	}

	triplets.reserve( entries );

	for( unsigned int r=0 ; r<rasterLines.size() ; r++ )
	{
		const RasterLine & line = rasterLines[r];
		size_t deepOffset = static_cast< size_t >( line.coeffStartIndex );
		int lineLength = line.lineEndIndex - line.lineStartIndex + 1;
		size_t prevCol = static_cast< size_t >( line.prevLineIndex );
		size_t current = static_cast< size_t >( line.lineStartIndex );
		size_t nextCol = static_cast< size_t >( line.nextLineIndex );
		const MatrixReal *mCoefficients = massCoefficients.deepCoefficients.data() + deepOffset * 10;
		const MatrixReal *sCoefficients = stiffnessCoefficients.deepCoefficients.data() + deepOffset * 10;


		for( int i=0 ; i<lineLength ; i++ , prevCol++ , current++ , nextCol++ , mCoefficients+=10 , sCoefficients+=10 )
		{
			if constexpr( Mass )
			{
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( prevCol-1 ) , mCoefficients[0] );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( prevCol   ) , mCoefficients[1] );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( prevCol+1 ) , mCoefficients[2] );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( current-1 ) , mCoefficients[3] );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( current   ) , mCoefficients[4] );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( current+1 ) , mCoefficients[5] );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( nextCol-1 ) , mCoefficients[6] );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( nextCol   ) , mCoefficients[7] );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( nextCol+1 ) , mCoefficients[8] );
			}
			else if constexpr( Stiffness )
			{
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( prevCol-1 ) , sCoefficients[0] );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( prevCol   ) , sCoefficients[1] );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( prevCol+1 ) , sCoefficients[2] );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( current-1 ) , sCoefficients[3] );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( current   ) , sCoefficients[4] );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( current+1 ) , sCoefficients[5] );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( nextCol-1 ) , sCoefficients[6] );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( nextCol   ) , sCoefficients[7] );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( nextCol+1 ) , sCoefficients[8] );
			}
			else
			{
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( prevCol-1 ) , mCoefficients[0]*mWeight + sCoefficients[0]*sWeight );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( prevCol   ) , mCoefficients[1]*mWeight + sCoefficients[1]*sWeight );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( prevCol+1 ) , mCoefficients[2]*mWeight + sCoefficients[2]*sWeight );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( current-1 ) , mCoefficients[3]*mWeight + sCoefficients[3]*sWeight );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( current   ) , mCoefficients[4]*mWeight + sCoefficients[4]*sWeight );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( current+1 ) , mCoefficients[5]*mWeight + sCoefficients[5]*sWeight );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( nextCol-1 ) , mCoefficients[6]*mWeight + sCoefficients[6]*sWeight );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( nextCol   ) , mCoefficients[7]*mWeight + sCoefficients[7]*sWeight );
				triplets.emplace_back( static_cast< int >( current ) , static_cast< int >( nextCol+1 ) , mCoefficients[8]*mWeight + sCoefficients[8]*sWeight );
			}
		}
	}

	for( unsigned int i=0 ; i<indexConverter.numBoundary() ; i++ )
	{
		size_t globalIndex = static_cast< size_t >( indexConverter.boundaryToCombined( AtlasBoundaryTexelIndex(i) ) );

		if constexpr( Mass )
		{
			for( unsigned int j=0 ; j<massCoefficients.boundaryDeepMatrix.RowSize(i) ; j++ )
				triplets.emplace_back( static_cast< int >( globalIndex ) , massCoefficients.boundaryDeepMatrix[i][j].N , massCoefficients.boundaryDeepMatrix[i][j].Value );

			for( unsigned int j=0 ; j<massCoefficients.boundaryBoundaryMatrix.RowSize(i) ; j++ )
			{
				size_t neighbourGlobalIndex = static_cast< size_t >( indexConverter.boundaryToCombined( AtlasBoundaryTexelIndex( massCoefficients.boundaryBoundaryMatrix[i][j].N ) ) );
				triplets.emplace_back( static_cast< int >( globalIndex ) , static_cast< int >( neighbourGlobalIndex ) , massCoefficients.boundaryBoundaryMatrix[i][j].Value );
			}
		}
		else if constexpr( Stiffness )
		{
			for( unsigned int j=0 ; j<stiffnessCoefficients.boundaryDeepMatrix.RowSize(i) ; j++ )
				triplets.emplace_back( static_cast< int >( globalIndex ) , stiffnessCoefficients.boundaryDeepMatrix[i][j].N , stiffnessCoefficients.boundaryDeepMatrix[i][j].Value );

			for( unsigned int j=0 ; j<stiffnessCoefficients.boundaryBoundaryMatrix.RowSize(i) ; j++ )
			{
				size_t neighbourGlobalIndex = static_cast< size_t >( indexConverter.boundaryToCombined( AtlasBoundaryTexelIndex( stiffnessCoefficients.boundaryBoundaryMatrix[i][j].N ) ) );
				triplets.emplace_back( static_cast< int >( globalIndex ) , static_cast< int >( neighbourGlobalIndex ) , stiffnessCoefficients.boundaryBoundaryMatrix[i][j].Value );
			}
		}
		else
		{
			for( unsigned int j=0 ; j<massCoefficients.boundaryDeepMatrix.RowSize(i) ; j++ )
				triplets.emplace_back( static_cast< int >( globalIndex ) , massCoefficients.boundaryDeepMatrix[i][j].N , massCoefficients.boundaryDeepMatrix[i][j].Value*mWeight + stiffnessCoefficients.boundaryDeepMatrix[i][j].Value*sWeight );

			for( unsigned int j=0 ; j<massCoefficients.boundaryBoundaryMatrix.RowSize(i) ; j++ )
			{
				size_t neighbourGlobalIndex = static_cast< size_t >( indexConverter.boundaryToCombined( AtlasBoundaryTexelIndex( massCoefficients.boundaryBoundaryMatrix[i][j].N ) ) );
				triplets.emplace_back( static_cast< int >( globalIndex ) , static_cast< int >( neighbourGlobalIndex ) , massCoefficients.boundaryBoundaryMatrix[i][j].Value*mWeight + stiffnessCoefficients.boundaryBoundaryMatrix[i][j].Value*sWeight );
			}
		}
	}

	Eigen::SparseMatrix< OutReal > M( indexConverter.numCombined() , indexConverter.numCombined() );
	M.setFromTriplets( triplets.begin() , triplets.end() );
	return M;
}

////////////////////////
// DivergenceOperator //
////////////////////////
template< typename MatrixReal >
template< bool SanityCheck >
std::vector< typename DivergenceOperator< MatrixReal >::DivergenceRasterLine >
DivergenceOperator< MatrixReal >::DivergenceRasterLine::GetRasterLines
(
	const Map< SimplexIndex< 1 , AtlasTexelIndex > , unsigned int > &edgeToIndex ,
	const std::vector< RasterLine > &rasterLines
)
{
	std::vector< typename DivergenceOperator< MatrixReal >::DivergenceRasterLine > divergenceRasterLines( rasterLines.size() );
	for( unsigned int i=0 ; i<rasterLines.size() ; i++ )
	{
		const RasterLine & line = rasterLines[i];
		DivergenceRasterLine & divLine = divergenceRasterLines[i];
		divLine.texelStart = line.lineStartIndex;
		divLine.texelEnd = line.lineEndIndex;
		divLine.deepCoefficientsStart = line.coeffStartIndex;


		{
			SimplexIndex< 1 , AtlasTexelIndex > prevEdgeKey( line.prevLineIndex-1 , line.prevLineIndex );
			auto iter = edgeToIndex.find( prevEdgeKey );
			if constexpr( SanityCheck ) if( iter==edgeToIndex.end() ) MK_THROW( "Edge not found" );
			divLine.prevEdgeRowStart = iter->second;
		}

		{
			SimplexIndex< 1 , AtlasTexelIndex > currEdgeKey( line.lineStartIndex-1 , line.lineStartIndex );
			auto iter = edgeToIndex.find( currEdgeKey );
			if constexpr( SanityCheck ) if( iter==edgeToIndex.end() ) MK_THROW( "Edge not found" );
			divLine.currEdgeRowStart = iter->second;
		}

		{
			SimplexIndex< 1 , AtlasTexelIndex > nextEdgeKey( line.nextLineIndex-1 , line.nextLineIndex );
			auto iter = edgeToIndex.find( nextEdgeKey );
			if constexpr( SanityCheck ) if( iter==edgeToIndex.end() ) MK_THROW( "Edge not found" );
			divLine.nextEdgeRowStart = iter->second;
		}
	}
	return divergenceRasterLines;
}

template< typename Real >
template< typename Data >
void DivergenceOperator< Real >::operator()( ConstPointer( Data ) edgeValues , Pointer( Data ) texelDivergence , bool add ) const
{
	// Update Boundary Texels 
	if( add ) boundaryMatrix.Multiply( &edgeValues[0] , &texelDivergence[0] , MULTIPLY_ADD );
	else      boundaryMatrix.Multiply( &edgeValues[0] , &texelDivergence[0] );

	// Update Deep Texels

	auto UpdateRow = [&]( int r )
		{
			Data* out = texelDivergence + rasterLines[r].texelStart;
			const Data* previousRowEdges = edgeValues + rasterLines[r].prevEdgeRowStart;
			const Data* currentRowEdges = edgeValues + rasterLines[r].currEdgeRowStart;
			const Data* nextRowEdges = edgeValues + rasterLines[r].nextEdgeRowStart;

			const Real * coeff = deepCoefficients.data() + static_cast< unsigned int >( rasterLines[r].deepCoefficientsStart ) * 12;
			int lineLength = rasterLines[r].texelEnd - rasterLines[r].texelStart + 1;

			if( add )
			{
				for( int i=0 ; i<lineLength ; coeff+=12 , previousRowEdges+=2 , currentRowEdges+=2 , nextRowEdges+=2 , i++ )
				{
					out[i] +=
						previousRowEdges[0] * coeff[0] +
						previousRowEdges[1] * coeff[1] +
						previousRowEdges[2] * coeff[2] +
						previousRowEdges[3] * coeff[3] +
						previousRowEdges[5] * coeff[4] +

						currentRowEdges[0] * coeff[5] +
						currentRowEdges[1] * coeff[6] +
						currentRowEdges[2] * coeff[7] +
						currentRowEdges[3] * coeff[8] +
						currentRowEdges[5] * coeff[9] +

						nextRowEdges[0] * coeff[10] +
						nextRowEdges[2] * coeff[11];
				}
			}
			else
			{
				for( int i=0 ; i<lineLength ; coeff+=12 , previousRowEdges+=2 , currentRowEdges+=2 , nextRowEdges+=2 , i++ )
				{
					out[i] =
						previousRowEdges[0] * coeff[0] +
						previousRowEdges[1] * coeff[1] +
						previousRowEdges[2] * coeff[2] +
						previousRowEdges[3] * coeff[3] +
						previousRowEdges[5] * coeff[4] +

						currentRowEdges[0] * coeff[5] +
						currentRowEdges[1] * coeff[6] +
						currentRowEdges[2] * coeff[7] +
						currentRowEdges[3] * coeff[8] +
						currentRowEdges[5] * coeff[9] +

						nextRowEdges[0] * coeff[10] +
						nextRowEdges[2] * coeff[11];
				}
			}

			// Edge indexing
			//		 ---0---
			//		|  
			//		1  
			//		|  


			// Reduced interior texel neighbour edge indexing
			//		--0-----2-- 
			//		|    |    |
			//		1    3    4
			//		|    |    |
			//		--5----7--
			//		|    |    |
			//		6    8    9
			//		|    |    |
			//      --10---11--


		};

	ThreadPool::ParallelFor( 0 , rasterLines.size() , [&]( unsigned int , size_t r ){ UpdateRow((unsigned int)r); } );
}

template< typename Real >
template< typename Data >
void DivergenceOperator< Real >::operator()( const std::vector< Data > & edgeValues , std::vector< Data > & texelDivergence , bool add ) const
{
	return operator()( GetPointer( edgeValues ) , GetPointer( texelDivergence ) , add );
}

template< typename Real >
template< typename Data >
std::vector< Data > DivergenceOperator< Real >::operator()( const std::vector< Data > &edgeValues ) const
{
	std::vector< Data > texelDivergence( boundaryMatrix.rows );
	operator()( edgeValues , texelDivergence , false );
	return texelDivergence;
}

template< typename Real >
template< typename OutReal >
Eigen::SparseMatrix< OutReal > DivergenceOperator< Real >::operator()( void ) const
{
	std::vector< Eigen::Triplet< OutReal > > triplets;
	Eigen::SparseMatrix< OutReal > M( boundaryMatrix.rows , edges.size() );

	size_t deepEntries = 0;
	for( unsigned int r=0 ; r<rasterLines.size() ; r++ ) deepEntries += ( rasterLines[r].texelEnd - rasterLines[r].texelStart + 1 ) + 12;
	triplets.reserve( deepEntries + boundaryMatrix.Entries() );

	for( unsigned int i=0 ; i<boundaryMatrix.rows ; i++ ) for( unsigned int j=0 ; j<boundaryMatrix.rowSizes[i] ; j++ ) triplets.emplace_back( i , boundaryMatrix[i][j].N , static_cast< OutReal >( boundaryMatrix[i][j].Value ) );

	for( unsigned int r=0 ; r<rasterLines.size() ; r++ )
	{
		size_t row = rasterLines[r].texelStart;
		size_t prevCol = static_cast< unsigned int >( rasterLines[r].prevEdgeRowStart );
		size_t currCol = static_cast< unsigned int >( rasterLines[r].currEdgeRowStart );
		size_t nextCol = static_cast< unsigned int >( rasterLines[r].nextEdgeRowStart );
		const Real * coeff = deepCoefficients.data() + static_cast< unsigned int >( rasterLines[r].deepCoefficientsStart ) * 12;
		int lineLength = rasterLines[r].texelEnd - rasterLines[r].texelStart + 1;

		for( int i=0 ; i<lineLength ; coeff+=12 , prevCol+=2 , currCol+=2 , nextCol+=2 , i++ )
		{
			triplets.emplace_back( static_cast< int >( row+i ) , static_cast< int >( prevCol+0 ) , static_cast< OutReal >( coeff[ 0] ) );
			triplets.emplace_back( static_cast< int >( row+i ) , static_cast< int >( prevCol+1 ) , static_cast< OutReal >( coeff[ 1] ) );
			triplets.emplace_back( static_cast< int >( row+i ) , static_cast< int >( prevCol+2 ) , static_cast< OutReal >( coeff[ 2] ) );
			triplets.emplace_back( static_cast< int >( row+i ) , static_cast< int >( prevCol+3 ) , static_cast< OutReal >( coeff[ 3] ) );
			triplets.emplace_back( static_cast< int >( row+i ) , static_cast< int >( prevCol+5 ) , static_cast< OutReal >( coeff[ 4] ) );
			triplets.emplace_back( static_cast< int >( row+i ) , static_cast< int >( currCol+0 ) , static_cast< OutReal >( coeff[ 5] ) );
			triplets.emplace_back( static_cast< int >( row+i ) , static_cast< int >( currCol+1 ) , static_cast< OutReal >( coeff[ 6] ) );
			triplets.emplace_back( static_cast< int >( row+i ) , static_cast< int >( currCol+2 ) , static_cast< OutReal >( coeff[ 7] ) );
			triplets.emplace_back( static_cast< int >( row+i ) , static_cast< int >( currCol+3 ) , static_cast< OutReal >( coeff[ 8] ) );
			triplets.emplace_back( static_cast< int >( row+i ) , static_cast< int >( currCol+5 ) , static_cast< OutReal >( coeff[ 9] ) );
			triplets.emplace_back( static_cast< int >( row+i ) , static_cast< int >( nextCol+0 ) , static_cast< OutReal >( coeff[10] ) );
			triplets.emplace_back( static_cast< int >( row+i ) , static_cast< int >( nextCol+2 ) , static_cast< OutReal >( coeff[11] ) );
		}
	}
	M.setFromTriplets( triplets.begin() , triplets.end() );
	return M;
}

/////////////////////////
// OperatorInitializer //
/////////////////////////

template< unsigned int Samples , bool SanityCheck , typename GeometryReal , typename MatrixReal >
void OperatorInitializer::_InitializeChart
(
	const ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > &texture_metrics ,
	const AtlasChart< GeometryReal > &atlasChart ,
	const GridChart< GeometryReal > &gridChart ,
	const typename GridAtlas<>::IndexConverter & indexConverter ,
	const std::vector< AtlasInteriorOrBoundaryNodeIndex >& fineBoundaryIndex ,
	std::vector< MatrixReal >& deepMassCoefficients ,
	std::vector< MatrixReal >& deepStiffnessCoefficients ,
	std::vector< MassAndStiffnessCoefficient< MatrixReal > > & boundaryBoundaryMassAndStiffness ,
	std::vector< MassAndStiffnessCoefficient< MatrixReal > > & boundaryDeepMassAndStiffness ,
	bool computeDivergence ,
	Map< SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > , AtlasRefinedBoundaryEdgeIndex > & fineBoundaryEdgeIndex,
	Map< SimplexIndex< 1 , AtlasTexelIndex > , unsigned int > & coarseEdgeIndex,
	std::vector< Eigen::Triplet< MatrixReal > > & boundaryDeepDivergenceTriplets,
	std::vector< Eigen::Triplet< MatrixReal > > & boundaryBoundaryDivergenceTriplets,
	std::vector< MatrixReal > & deepDivergenceCoefficients
)
{
	ExplicitIndexVector< ChartInteriorCellIndex , SquareMatrix< GeometryReal , 4 > > cellStiffness( gridChart.numInteriorCells() );
	ExplicitIndexVector< ChartInteriorCellIndex , SquareMatrix< GeometryReal , 4 > > cellMass( gridChart.numInteriorCells() );

	ExplicitIndexVector< ChartBoundaryTriangleIndex , SquareMatrix< GeometryReal , 6 > > triangleElementStiffness( static_cast< unsigned int >(gridChart.endBoundaryTriangleIndex) );
	ExplicitIndexVector< ChartBoundaryTriangleIndex , SquareMatrix< GeometryReal , 6 > > triangleElementMass( static_cast< unsigned int >(gridChart.endBoundaryTriangleIndex) );

	ExplicitIndexVector< ChartInteriorCellIndex , SquareMatrix< GeometryReal , 4 > > cellDivergence;
	if( computeDivergence ) cellDivergence.resize( gridChart.numInteriorCells() );
	ExplicitIndexVector< ChartBoundaryTriangleIndex , Matrix< GeometryReal , 6 , 15 > > triangleElementDivergence;
	if( computeDivergence ) triangleElementDivergence.resize( static_cast< unsigned int >(gridChart.endBoundaryTriangleIndex) );

	//// Rasterize
	int zeroAreaElementCount = 0;
	GeometryReal precision_error = (GeometryReal)PRECISION_EPSILON;

	// Node indexing
	//		0 ---- 1
	//		|	   |
	//		|	   |
	//		3------2

	int cellCornerPairs[12] = { 0,1,0,3,0,2,1,3,3,2,1,2 };

	// Cell edge indexing
	//		+---0--+
	//		|\    /|
	//		| 2  / |
	//		|  \/  5
	//		1  /\  |
	//		| 3  \ |
	//		|/    \|
	//      +---4--+

	int reducedCellCornerPairs[8] = { 0,1,0,3,3,2,1,2 };

	// Reduced cell edge indexing
	//		+--0---+
	//		|      |
	//		1      3
	//		|      |
	//		+--2---+


	auto InUnitSquare =   [&]( Point2D< GeometryReal > p ){ return !( p[0]<0-precision_error || p[1]<0-precision_error || p[0]>1+precision_error || p[1]>1+precision_error ); };
	auto InUnitTriangle = [&]( Point2D< GeometryReal > p ){ return !( p[0]<0-precision_error || p[1]<0-precision_error || ( p[0]+p[1] )>1+precision_error ); };
	auto CellInTriangle = [&]( int i , int j , const std::vector< Point2D< GeometryReal > >& vertices )
		{
			Point2D< GeometryReal > points[] = { gridChart.nodePosition(i,j) , gridChart.nodePosition(i+1,j) , gridChart.nodePosition(i+1,j+1) , gridChart.nodePosition(i,j+1) };
			unsigned int count = 0;
			for( unsigned int i=0 ; i<3 ; i++ )
			{
				SimplexIndex< 1 > eIndex = OutgoingEdgeIndex(i);
				EdgeEquation< GeometryReal > eq = EdgeEquation< GeometryReal >( vertices[ eIndex[0] ] , vertices[ eIndex[1] ] );
				for( unsigned int k=0 ; k<4 ; k++ ) if( eq( points[k] )<0 ) count++;
			}
			return count==0 || count==12;
		};

	// Compute the square-root of the weights to make taking the weighted dot-product faster
	GeometryReal _integrator_sampleWeight[Samples];
	for( int s=0 ; s<Samples ; s++ ) _integrator_sampleWeight[s] = (GeometryReal)sqrt( TriangleIntegrator< Samples >::Weights[s] );
	SquareMatrix< GeometryReal , 2 > cell_to_texture_differential;
	cell_to_texture_differential(0,0) = (GeometryReal)gridChart.cellSizeW;
	cell_to_texture_differential(1,1) = (GeometryReal)gridChart.cellSizeH;
	cell_to_texture_differential(0,1) = cell_to_texture_differential(1,0) = 0;

	SquareMatrix< GeometryReal , 4 > interior_cell_mass;
	SquareMatrix< GeometryReal , 2 > interior_cell_stiffnesses[4][4];
	SquareMatrix< GeometryReal , 2 > grad_edge_products[4][4];

	{
		std::vector< Point2D< GeometryReal > > polygon = { Point2D< GeometryReal >(0,0) , Point2D< GeometryReal >(1,0) , Point2D< GeometryReal >(1,1) , Point2D< GeometryReal >(0,1) };
		for( unsigned int p=2 ; p<polygon.size() ; p++ )
		{
			Point2D< GeometryReal > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };
			Point2D< GeometryReal > fragment_samples[Samples];
			for( unsigned int s=0 ; s<Samples ; s++ ) fragment_samples[s] = polygon[0] + dm[0] * (GeometryReal)TriangleIntegrator<Samples>::Positions[s][0] + dm[1] * (GeometryReal)TriangleIntegrator<Samples>::Positions[s][1];

			// Integrate scalar product and gradient field
			GeometryReal sampleValues[Samples][4];
			Point2D< GeometryReal > sampleGradients[Samples][4];
			Point2D< GeometryReal > __sampleGradients[Samples][4];

			for( unsigned int s=0 ; s<Samples ; s++ )
			{
				BilinearElementValuesAndGradients( fragment_samples[s] , sampleValues[s] , sampleGradients[s] );
				for( int k=0 ; k<4 ; k++ )
				{
					sampleValues   [s][k] *= _integrator_sampleWeight[s];
					__sampleGradients[s][k] = sampleGradients[s][k];
					sampleGradients[s][k] *= _integrator_sampleWeight[s];
				}
			}

			for( unsigned int k=0 ; k<4 ; k++ ) for( unsigned int l=0 ; l<4 ; l++ ) for( unsigned int s=0 ; s<Samples ; s++ )
			{
				interior_cell_mass(l,k) += sampleValues[s][k] * sampleValues[s][l] / 2;
				for( int m=0 ; m<2 ; m++ ) for( int n=0 ; n<2 ; n++ ) interior_cell_stiffnesses[l][k](m,n) += sampleGradients[s][l][m] * sampleGradients[s][k][n] / 2;
			}

			if( computeDivergence )
			{
				Point2D< GeometryReal > sampleVectorFields[Samples][4];
				for( unsigned int s=0 ; s<Samples ; s++ )
				{
					ReducedVectorFieldBasis( fragment_samples[s] , sampleVectorFields[s] );
					for( unsigned int k=0 ; k<4 ; k++ ) sampleVectorFields[s][k] *= _integrator_sampleWeight[s];
				}
				for( int k=0 ; k<4 ; k++ ) for( int l=0 ; l<4 ; l++ ) for( unsigned int s=0 ; s<Samples ; s++ )
					for( unsigned int m=0 ; m<2 ; m++ ) for( int n=0 ; n<2 ; n++ ) grad_edge_products[l][k](m,n) += sampleVectorFields[s][k][m] * sampleGradients[s][l][n] / 2;
			}
		}
	}

	using Index = RegularGrid< 2 >::Index;
	using Range = RegularGrid< 2 >::Range;
	Range cellRange;
	cellRange.second[0] = gridChart.width-1;
	cellRange.second[1] = gridChart.height-1;

	auto GetSimplex = [&]( unsigned int t )
		{
			Simplex< double , 2 , 2 > simplex;
			for( unsigned int i=0 ; i<=2 ; i++ )
			{
				simplex[i] = atlasChart.vertex( ChartMeshVertexIndex( atlasChart.triangleIndex( ChartMeshTriangleIndex(t) )[i] ) ) - gridChart.corner;
				simplex[i][0] /= gridChart.cellSizeW , simplex[i][1] /= gridChart.cellSizeH;
			}
			return simplex;
		};

	// Add the bounday cell contribution
	for( size_t i=0 ; i<gridChart.cellIndices.size() ; i++ )
	{
		if( gridChart.cellIndices[i].boundary!=ChartBoundaryCellIndex(-1) )
		{
			Index I( static_cast< unsigned int >( i % gridChart.cellIndices.res(0) ) , static_cast< unsigned int >( i / gridChart.cellIndices.res(0) ) );
			auto TextureToCell = [&]( Point2D< GeometryReal > p ){ return Point2D< GeometryReal >( (GeometryReal)( p[0] / gridChart.cellSizeW ) - I[0] , (GeometryReal)( p[1] / gridChart.cellSizeH ) - I[1] ); };

			ChartBoundaryCellIndex boundaryIndex = gridChart.cellIndices[i].boundary;
			const std::vector< BoundaryIndexedTriangle< GeometryReal > > & boundaryTriangles = gridChart.boundaryTriangles[boundaryIndex];

			// Iterate over all elements associated with the cell
			for( unsigned int bt=0 ; bt<boundaryTriangles.size() ; bt++ )
			{
				BoundaryIndexedTriangle< GeometryReal > element = boundaryTriangles[bt];
				ChartMeshTriangleIndex tIdx = element.sourceIdx;

				// Grab the vertices in the cell frame
				std::vector< Point2D< GeometryReal > > element_vertices(3);
				for( unsigned int ii=0 ; ii<3 ; ii++ ) element_vertices[ii] = TextureToCell( element[ii] );
				ChartBoundaryTriangleIndex boundaryTriangleIndex = element.idx;

				IndexedPolygon< GeometryReal > polygon;

				SetIndexedPolygonFromBoundaryTriangle( element , polygon );

				// Convert the polygon vertices from the texture frame to the cell frame
				for( int ii=0 ; ii<polygon.size() ; ii++ ) polygon[ii] = TextureToCell( polygon[ii] );

				SquareMatrix< GeometryReal , 2 > element_to_cell_differential , cell_to_element_differential;
				Point2D< GeometryReal > dm[] = { element_vertices[1]-element_vertices[0] , element_vertices[2]-element_vertices[0] };
				for( unsigned int x=0 ; x<2 ; x++ ) for( unsigned int y=0 ; y<2 ; y++ ) element_to_cell_differential(x,y) = dm[x][y];
				cell_to_element_differential = element_to_cell_differential.inverse();
				auto CellToElement = [&]( Point2D< GeometryReal > p ){ return cell_to_element_differential * ( p - element_vertices[0] ); };
				// Convert the polygon vertices from the cell frame to the element frame
				for( int ii=0 ; ii<polygon.size() ; ii++ ) polygon[ii] = CellToElement( polygon[ii] );

				SquareMatrix< GeometryReal , 2 > cell_metric = cell_to_texture_differential.transpose() * texture_metrics[ tIdx ] * cell_to_texture_differential;
				SquareMatrix< GeometryReal , 2 > element_metric = element_to_cell_differential.transpose() * cell_metric * element_to_cell_differential;
				SquareMatrix< GeometryReal , 2 > element_metric_inverse = element_metric.inverse();
				GeometryReal element_area_scale_factor = sqrt( element_metric.determinant() );
				SquareMatrix< GeometryReal , 6 > polygonStiffness , polygonMass;
				Matrix< GeometryReal , 6 , 15 > polygonDivergence;
				GeometryReal polygonArea = 0;

				for( unsigned int p=2 ; p<polygon.size() ; p++ )
				{
					Point2D< GeometryReal > d[] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };

					SquareMatrix< GeometryReal , 2 > fragment_to_element_differential;
					for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = d[x][y];
					GeometryReal fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );
					GeometryReal fragment_area = fragment_to_element_area_scale_factor * element_area_scale_factor / 2;
					if( fragment_area>0 )
					{
						polygonArea += fragment_area;

						Point2D< GeometryReal > fragment_samples[Samples];
						for( int s=0 ; s<Samples ; s++ )
						{
							fragment_samples[s] = polygon[0] + d[0] * (GeometryReal)TriangleIntegrator<Samples>::Positions[s][0] + d[1] * (GeometryReal)TriangleIntegrator<Samples>::Positions[s][1];
							if constexpr( SanityCheck ) if( !InUnitTriangle( fragment_samples[s] ) ) MK_THROW( "Boundary sample out of unit right triangle! (" , fragment_samples[s][0] , " " , fragment_samples[s][1] , ")" );
							{
								fragment_samples[s][0] = std::max< GeometryReal >( fragment_samples[s][0] , 0 );
								fragment_samples[s][1] = std::max< GeometryReal >( fragment_samples[s][1] , 0 );
								GeometryReal excess = ( fragment_samples[s][0] + fragment_samples[s][1] ) - 1;
								if( excess>0 ) fragment_samples[s][0] -= excess/2 , fragment_samples[s][1] -= excess/2;
							}
						}
						GeometryReal sampleValues[Samples][6];
						Point2D< GeometryReal > sampleGradients[Samples][6] , _sampleGradients[Samples][6], __sampleGradients[Samples][6];
						for( int s=0 ; s<Samples ; s++ )
						{
							QuadraticElement::ValuesAndDifferentials( fragment_samples[s] , sampleValues[s] , sampleGradients[s] );
							for( int k=0 ; k<6 ; k++ )
							{
								sampleValues[s][k] *= _integrator_sampleWeight[s];
								__sampleGradients[s][k] = sampleGradients[s][k];
								sampleGradients[s][k] *= _integrator_sampleWeight[s];
								_sampleGradients[s][k] = element_metric_inverse * sampleGradients[s][k];										
							}
						}

						for( unsigned int k=0 ; k<6 ; k++ ) for( unsigned int l=0 ; l<=k ; l++ )
						{
							GeometryReal vIntegral=0 , gIntegral=0;
							for( int s=0 ; s<Samples ; s++ )
							{
								vIntegral += sampleValues[s][k] * sampleValues[s][l];
								gIntegral += Point2D< GeometryReal >::Dot( sampleGradients[s][k] , _sampleGradients[s][l] );
							}
							polygonMass(l,k) += vIntegral * fragment_area;
							polygonStiffness(l,k) += gIntegral * fragment_area;
						}

						if( computeDivergence )
						{
							Point2D< GeometryReal > sampleVectorFields[Samples][15];

							for( int k=1 ; k<6 ; k++ ) for( int l=0 ; l<k ; l++ )
							{
								int edgeId = (k*(k - 1)) / 2 + l;
								for( int s=0 ; s<Samples ; s++ )
									sampleVectorFields[s][edgeId] = sampleValues[s][k] * __sampleGradients[s][l] - sampleValues[s][l] * __sampleGradients[s][k];
							}

							for( int k=0 ; k<15 ; k++ ) for( int l=0 ; l<6 ; l++ )
							{
								GeometryReal dIntegral = 0;
								for( int s=0 ; s<Samples ; s++ ) dIntegral += Point2D< GeometryReal >::Dot( sampleVectorFields[s][k] , _sampleGradients[s][l] );
								polygonDivergence(l,k) += dIntegral * fragment_area;
							}
						}
					}
					else
					{
						zeroAreaElementCount++;
						MK_WARN( "Zero area polygon at cell " , gridChart.cornerCoords[0] + I[0] , " " , gridChart.cornerCoords[1] + I[1] );
					}
				}
				for( unsigned int k=0 ; k<6 ; k++ ) for( unsigned int l=0 ; l<k ; l++ )
				{
					polygonMass(k,l) = polygonMass(l,k);
					polygonStiffness(k,l) = polygonStiffness(l,k);
				}

				GeometryReal integratedPolygonMass = 0;
				for( int k=0 ; k<6 ; k++ ) for( int l=0 ; l<6 ; l++ ) integratedPolygonMass += polygonMass(k,l);
				if constexpr( SanityCheck ) if( fabs( integratedPolygonMass - polygonArea )>precision_error ) MK_WARN( "Out of precision" );
				{
					triangleElementMass[ boundaryTriangleIndex ] += polygonMass;
					triangleElementStiffness[ boundaryTriangleIndex ] += polygonStiffness;
					if( computeDivergence ) triangleElementDivergence[ boundaryTriangleIndex ] += polygonDivergence;
				}
			}
		}
	}

	for( size_t t=0 ; t<atlasChart.numTriangles() ; t++ )
	{
		Point2D< GeometryReal > tPos[3];
		SimplexIndex< 2 , ChartMeshVertexIndex > tri = atlasChart.triangleIndex( ChartMeshTriangleIndex(t) );
		for( unsigned int i=0 ; i<3 ; i++ ) tPos[i] = atlasChart.vertex( tri[i] ) - gridChart.corner;

		SquareMatrix< GeometryReal , 2 > texture_metric = texture_metrics[ ChartMeshTriangleIndex(t) ];
		SquareMatrix< GeometryReal , 2 > cell_metric = cell_to_texture_differential.transpose() * texture_metric * cell_to_texture_differential;
		SquareMatrix< GeometryReal , 2 > cell_metric_inverse = cell_metric.inverse();
		GeometryReal cell_area_scale_factor = (GeometryReal)sqrt( cell_metric.determinant() );

		std::vector< Point2D< GeometryReal > > parametricVertices(3);
		parametricVertices[0] = tPos[0] , parametricVertices[1] = tPos[1] , parametricVertices[2] = tPos[2];

		IndexedTriangle< GeometryReal > atlasTriangle;
		for( unsigned int k=0 ; k<3 ; k++ )
		{
			atlasTriangle.vertices[k] = tPos[k];
			atlasTriangle.atlasEdgeIndices[k] = atlasChart.atlasEdge( GetChartMeshHalfEdgeIndex( ChartMeshTriangleIndex(t) , k ) );
			atlasTriangle.vertexIndices[k] = atlasChart.triangleIndex( ChartMeshTriangleIndex(t) )[k];
			atlasTriangle.atlasVertexParentEdge[k] = AtlasMeshEdgeIndex( -1 );
		}

		auto Kernel = [&]( Index I )
			{
				int i = I[0] , j = I[1];
				auto TextureToCell = [&]( Point2D< GeometryReal > p ){ return Point2D< GeometryReal >( (GeometryReal)( p[0] / gridChart.cellSizeW ) - i , (GeometryReal)( p[1] / gridChart.cellSizeH ) - j ); };

				ChartInteriorCellIndex interiorIndex = gridChart.cellIndices(i,j).interior;
				ChartBoundaryCellIndex boundaryIndex = gridChart.cellIndices(i,j).boundary;
				if constexpr( SanityCheck ) if( interiorIndex!=ChartInteriorCellIndex(-1) && boundaryIndex!=ChartBoundaryCellIndex(-1) ) MK_THROW( "Cell simultaneously interior and boundary" );

				// If the cell is entirely within the triangle...
				if( CellInTriangle( i , j , parametricVertices ) && interiorIndex!=ChartInteriorCellIndex(-1) )
				{
					SquareMatrix< GeometryReal , 4 > polygonMass , polygonStiffness;

					polygonMass = interior_cell_mass * cell_area_scale_factor;
					for( unsigned int k=0 ; k<4 ; k++ ) for( unsigned int l=0 ; l<=k ; l++ )
						polygonStiffness(l,k) = polygonStiffness(k,l) = SquareMatrix< GeometryReal , 2 >::Dot( cell_metric_inverse , interior_cell_stiffnesses[k][l] ) * cell_area_scale_factor;
					cellMass[interiorIndex] += polygonMass;
					cellStiffness[interiorIndex] += polygonStiffness;

					if( computeDivergence )
					{
						SquareMatrix< GeometryReal , 4 > polygonDivergence;
						for( unsigned int k=0 ; k<4 ; k++ ) for( unsigned int l=0 ; l<4 ; l++ )
							polygonDivergence(l,k) = SquareMatrix< GeometryReal , 2 >::Dot( cell_metric_inverse , grad_edge_products[l][k] ) * cell_area_scale_factor;
						cellDivergence[interiorIndex] += polygonDivergence;
					}
				}
				else if( interiorIndex!=ChartInteriorCellIndex(-1) )
				{
					// For interior cells, the cell and the element are the same thing
					auto TextureToElement = TextureToCell;

					SquareMatrix< GeometryReal , 2 > element_metric = cell_metric , element_metric_inverse = cell_metric_inverse;
					GeometryReal element_area_scale_factor = cell_area_scale_factor;
					CellClippedTriangle< GeometryReal > polygon = parametricVertices;

					// Clip the triangle to the cell
					if( ClipTriangleToPrimalCell( polygon , i , j , gridChart.cellSizeW , gridChart.cellSizeH ) )
					{
						// Transform the polygon vertices into the coordinate frame of the cell
						for( unsigned int ii=0 ; ii<polygon.size() ; ii++ ) polygon[ii] = TextureToElement( polygon[ii] );
						SquareMatrix< GeometryReal , 4 > polygonStiffness , polygonMass;
						SquareMatrix< GeometryReal , 4 > polygonDivergence;

						for( unsigned int p=2 ; p<polygon.size() ; p++ )
						{
							Point2D< GeometryReal > dm[2] = { polygon[p-1]-polygon[0] , polygon[p]-polygon[0] };

							SquareMatrix< GeometryReal , 2 > fragment_to_element_differential;
							for( int x=0 ; x<2 ; x++ ) for( int y=0 ; y<2 ; y++ ) fragment_to_element_differential(x,y) = dm[x][y];
							GeometryReal fragment_to_element_area_scale_factor = fabs( fragment_to_element_differential.determinant() );

							Point2D< GeometryReal > fragment_samples[Samples];
							for( int s=0 ; s<Samples ; s++ )
							{
								fragment_samples[s] = polygon[0] + dm[0] * (GeometryReal)TriangleIntegrator<Samples>::Positions[s][0] + dm[1] * (GeometryReal)TriangleIntegrator<Samples>::Positions[s][1];
								if constexpr( SanityCheck ) if( !InUnitSquare( fragment_samples[s] ) ) MK_THROW( "Interior sample out of unit box! (" , fragment_samples[s][0] , " " , fragment_samples[s][1] , ")" );
							}

							// Integrate scalar product and gradient field
							// Make the code more efficient by:
							// -- pre-multiplying the values and gradients by the square-root of the quadrature weight
							// -- pre-multiplying the gradients by the inverse of the metric
							// so that the computation within the inner loop is faster.
							GeometryReal fragment_area = element_area_scale_factor * fragment_to_element_area_scale_factor / 2;
							GeometryReal sampleValues[Samples][4];
							Point2D< GeometryReal > sampleGradients[Samples][4] , _sampleGradients[Samples][4], __sampleGradients[Samples][4];
							for( unsigned int s=0 ; s<Samples ; s++ )
							{
								BilinearElementValuesAndGradients( fragment_samples[s] , sampleValues[s] , sampleGradients[s] );
								for( unsigned int k=0 ; k<4 ; k++ )
								{
									sampleValues[s][k] *= _integrator_sampleWeight[s];
									__sampleGradients[s][k] = sampleGradients[s][k];
									sampleGradients[s][k] *= _integrator_sampleWeight[s];
									_sampleGradients[s][k] = element_metric_inverse * sampleGradients[s][k];
								}
							}
							for( unsigned int k=0 ; k<4 ; k++ ) for( unsigned int l=0 ; l<=k ; l++ )
							{
								GeometryReal vIntegral=0 , gIntegral=0;
								for( int s=0 ; s<Samples ; s++ )
								{
									vIntegral += sampleValues[s][k] * sampleValues[s][l];
									gIntegral += Point2D< GeometryReal >::Dot( sampleGradients[s][l] , _sampleGradients[s][k] );
								}
								polygonMass(l,k) += vIntegral * fragment_area;
								polygonStiffness(l,k) += gIntegral * fragment_area;
							}
							if( computeDivergence )
							{
								Point2D< GeometryReal > sampleVectorFields[Samples][4];
								for( unsigned int s=0 ; s<Samples ; s++ )
								{
									ReducedVectorFieldBasis(fragment_samples[s], sampleVectorFields[s]);
									for( unsigned int k=0 ; k<4 ; k++ ) sampleVectorFields[s][k] *= _integrator_sampleWeight[s];
								}
								for( unsigned int k=0 ; k<4 ; k++ ) for( unsigned int l=0 ; l<4 ; l++ )
								{
									GeometryReal dIntegral = 0;
									for( unsigned int s=0 ; s<Samples ; s++ ) dIntegral += Point2D< GeometryReal >::Dot( sampleVectorFields[s][k] , _sampleGradients[s][l] );
									polygonDivergence(l,k) += dIntegral * fragment_area;
								}
							}
						}
						for( unsigned int k=0 ; k<4 ; k++ ) for( unsigned int l=0 ; l<k ; l++ )
						{
							polygonMass(k,l) = polygonMass(l,k);
							polygonStiffness(k,l) = polygonStiffness(l,k);
						}
						cellMass[interiorIndex] += polygonMass;
						cellStiffness[interiorIndex] += polygonStiffness;
						if( computeDivergence ) cellDivergence[interiorIndex] += polygonDivergence;
					}
				}
			};
		Rasterizer2D::RasterizeSupports< true , true >( GetSimplex( static_cast< unsigned int >(t) ) , Kernel , cellRange );
	}

	if( zeroAreaElementCount ) MK_WARN( "Element with zero area = " , zeroAreaElementCount );

	int offset_i[4] = { 0, 1 , 1 ,0 };
	int offset_j[4] = { 0, 0 , 1 ,1 };
	int cellOffset[4] = { 3, 2, 0, 1 };



	// Node indexing
	//		0 ---- 1
	//		|	   |
	//		|	   |
	//		3------2

	int cellEdgeCornerOffset[24] =
	{
		13,  9,  0,  4,
		14, 10,  1,  5,
		15, 11,  2 , 6,
		16, 12,  3,  7,
		19, 18,  9, 13,
		17, 14,  5,  8
	};

	// Interior texel neighbour edge indexing
	//		+--0---+---4--+ 
	//		|\    /|\    /|  
	//		| 2  / | 6  / |
	//		|  \/  |  \/  |
	//		1  /\  5  /\  8
	//		| 3  \ | 7  \ |
	//		|/    \|/    \|
	//		+--9---+--13--+
	//		|\    /|\    /|
	//		| 11 / | 15 / |
	//		|  \/  |  \/  |
	//		10 /\ 14  /\ 17
	//		| 12 \ | 16 \ |
	//		|/    \|/	 \|
	//		+--18--+--19--+


	int reducedCellEdgeCornerOffset[24] =
	{
		7,  5,  0,  2,
		8,  6,  1,  3,
		11, 10,  5,  7,
		9,  8,  3,  4
	};


	// Reduced interior texel neighbour edge indexing
	//		+-0--+--2-+ 
	//		|    |    |
	//		1    3    4
	//		|    |    |
	//		+-5--+-7--+
	//		|    |    |
	//		6    8    9
	//		|    |    |
	//      +-10-+-11-+

	auto NeighbourOffset = [&]( int k , int l ){ return ( offset_j[l] - offset_j[k] + 1 ) * 3 + ( offset_i[l] - offset_i[k] + 1 ); };
	for( unsigned int i=0 ; i<gridChart.interiorCellCoveredTexelBilinearElementIndices.size() ; i++ )
	{
		const BilinearElementIndex< AtlasTexelIndex > & indicesCombined = gridChart.interiorCellCombinedTexelBilinearElementIndices[ ChartInteriorCellIndex(i) ];
		const BilinearElementIndex< AtlasCoveredTexelIndex > & indicesInterior = gridChart.interiorCellCoveredTexelBilinearElementIndices[ ChartInteriorCellIndex(i) ];

		AtlasCellIndex cellIndex = gridChart.chartToAtlasCombinedCellIndex( gridChart.interiorCellIndexToCombinedCellIndex[i] );

		unsigned int cellCoarseEdgeIndex[4];
		bool coarseEdgeIndexInitialized = false;

		for( unsigned int k=0 ; k<4 ; k++ )
		{
			AtlasTexelIndex currentNode = indicesCombined[k];
			AtlasBoundaryTexelIndex _currentBoundaryIndex = indexConverter.combinedToBoundary( currentNode );
			AtlasInteriorTexelIndex _currentInteriorIndex = indexConverter.combinedToInterior( currentNode );
			if( _currentInteriorIndex!=AtlasInteriorTexelIndex(-1) ) //Interior
			{
				for( unsigned int l=0 ; l<4 ; l++ )
				{
					deepMassCoefficients[ 10*static_cast< unsigned int >(_currentInteriorIndex) + NeighbourOffset(k,l) ] = (MatrixReal)( cellMass[ ChartInteriorCellIndex(i) ](k,l) + deepMassCoefficients[ 10*static_cast< unsigned int >(_currentInteriorIndex) + NeighbourOffset(k,l) ] );
					deepStiffnessCoefficients[ 10*static_cast< unsigned int >(_currentInteriorIndex) + NeighbourOffset(k,l) ] = (MatrixReal)( cellStiffness[ ChartInteriorCellIndex(i) ](k,l) + deepStiffnessCoefficients[ 10*static_cast< unsigned int >(_currentInteriorIndex) + NeighbourOffset(k,l) ] );
				}
				if( computeDivergence ) for( int l=0 ; l<4 ; l++ ) deepDivergenceCoefficients[ 12*static_cast< unsigned int >(_currentInteriorIndex) + reducedCellEdgeCornerOffset[4*l+k] ] = (MatrixReal)( cellDivergence[ ChartInteriorCellIndex(i) ](k,l) + deepDivergenceCoefficients[ 12*static_cast< unsigned int >(_currentInteriorIndex) + reducedCellEdgeCornerOffset[4*l+k] ] );
			}
			else if( _currentBoundaryIndex!=AtlasBoundaryTexelIndex(-1) )
			{
				for( unsigned int l=0 ; l<4 ; l++ )
				{
					AtlasTexelIndex neighborNode = indicesCombined[l];
					AtlasBoundaryTexelIndex neighborBoundaryIndex = indexConverter.combinedToBoundary( neighborNode );
					AtlasInteriorTexelIndex neighborInteriorIndex = indexConverter.combinedToInterior( neighborNode );
					if( neighborInteriorIndex!=AtlasInteriorTexelIndex(-1) )
					{
						boundaryDeepMassAndStiffness.emplace_back( static_cast< unsigned int >(_currentBoundaryIndex) , static_cast< unsigned int >(neighborNode) , (MatrixReal)cellMass[ ChartInteriorCellIndex(i) ](k,l) , (MatrixReal)cellStiffness[ ChartInteriorCellIndex(i) ](k,l) );
					}
					else if( neighborBoundaryIndex!=AtlasBoundaryTexelIndex(-1) )
					{
						int r = static_cast< unsigned int >( fineBoundaryIndex[ static_cast< unsigned int >(indicesInterior[k]) ] );
						int c = static_cast< unsigned int >( fineBoundaryIndex[ static_cast< unsigned int >(indicesInterior[l]) ] );
						if( r<=c ) boundaryBoundaryMassAndStiffness.emplace_back( r , c , (MatrixReal)cellMass[ ChartInteriorCellIndex(i) ](k,l) , (MatrixReal)cellStiffness[ ChartInteriorCellIndex(i) ](k,l) );
					}
					else if constexpr( SanityCheck ) MK_THROW( "Expected supported index" );
				}
				if( computeDivergence )
				{
					if( !coarseEdgeIndexInitialized )
					{
						for( unsigned int l=0 ; l<4 ; l++ )
						{
							AtlasTexelIndex edgeSourceCoarseIndex = indicesCombined[reducedCellCornerPairs[ 2*l+0 ] ];
							AtlasTexelIndex edgeTargetCoarseIndex = indicesCombined[reducedCellCornerPairs[ 2*l+1 ] ];
							SimplexIndex< 1 , AtlasTexelIndex > coarseEdgeKey( edgeSourceCoarseIndex , edgeTargetCoarseIndex );
							if constexpr( SanityCheck ) if( coarseEdgeIndex.find(coarseEdgeKey)==coarseEdgeIndex.end() ) MK_THROW( "Fine edge not found" );
							cellCoarseEdgeIndex[l] = coarseEdgeIndex[coarseEdgeKey];
						}
						coarseEdgeIndexInitialized = true;
					}
					for( unsigned int l=0 ; l<4 ; l++ ) boundaryDeepDivergenceTriplets.push_back( Eigen::Triplet< MatrixReal >( static_cast< unsigned int >(currentNode) , cellCoarseEdgeIndex[l] , (MatrixReal)cellDivergence[ ChartInteriorCellIndex(i) ](k,l) ) );
				}
			}
			else if constexpr( SanityCheck ) MK_THROW( "Expected supported index" );
		}
	}


	for( unsigned int c=0 ; c<gridChart.boundaryTriangles.size() ; c++ )
	{
		AtlasCellIndex cellIndex = gridChart.chartToAtlasCombinedCellIndex( gridChart.boundaryCellIndexToCombinedCellIndex[c] );

		const std::vector< BoundaryIndexedTriangle< GeometryReal > > & boundaryTriangles = gridChart.boundaryTriangles[ ChartBoundaryCellIndex(c) ];
		for( unsigned int b=0 ; b<boundaryTriangles.size() ; b++ )
		{
			ChartBoundaryTriangleIndex boundaryTriangleIndex = boundaryTriangles[b].idx;
			const QuadraticElement::Index & indices = boundaryTriangles[b].indices;
			QuadraticElement::Index fineTriangleElementIndices;
			for( unsigned int k=0 ; k<6 ; k++ ) fineTriangleElementIndices[k] = fineBoundaryIndex[ static_cast< unsigned int >(indices[k]) ];

			for( unsigned int k=0 ; k<6 ; k++ ) for( unsigned int l=0 ; l<6 ; l++ )
			{
				int r = static_cast< unsigned int >( fineTriangleElementIndices[k] );
				int c = static_cast< unsigned int >( fineTriangleElementIndices[l] );
				if( r<=c ) boundaryBoundaryMassAndStiffness.emplace_back( r , c , (MatrixReal)triangleElementMass[boundaryTriangleIndex](l,k) , (MatrixReal)triangleElementStiffness[boundaryTriangleIndex](l,k) );
			}

			if( computeDivergence ) for( unsigned int k=1 ; k<6 ; k++ ) for( unsigned int l=0 ; l<k ; l++ )
			{
				unsigned int edgeId = ( k*(k-1 ) ) / 2 + l;
				AtlasInteriorOrBoundaryNodeIndex edgeSourceFineIndex = fineTriangleElementIndices[k];
				AtlasInteriorOrBoundaryNodeIndex edgeTargetFineIndex = fineTriangleElementIndices[l];
				MatrixReal edgeSign = (MatrixReal)1.;
				if( edgeTargetFineIndex<edgeSourceFineIndex )
				{
					std::swap( edgeSourceFineIndex , edgeTargetFineIndex );
					edgeSign = (MatrixReal)-1.;
				}

				auto iter = fineBoundaryEdgeIndex.find( SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex >( edgeSourceFineIndex , edgeTargetFineIndex ) );
				if constexpr( SanityCheck ) if( iter==fineBoundaryEdgeIndex.end() ) MK_THROW( "Fine edge not found" );
				for( unsigned int n=0 ; n<6 ; n++ )
				{
					AtlasInteriorOrBoundaryNodeIndex fineNodeIndex = fineTriangleElementIndices[n];
					boundaryBoundaryDivergenceTriplets.push_back( Eigen::Triplet< MatrixReal >( static_cast< unsigned int >(fineNodeIndex) , static_cast< unsigned int >(iter->second) , (MatrixReal)triangleElementDivergence[boundaryTriangleIndex](n,edgeId) * edgeSign ) );
				}
			}
		}
	}
}

template< unsigned int Samples , bool SanityCheck , typename GeometryReal , typename MatrixReal >
void OperatorInitializer::_Initialize
(
	const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > &parameterMetric ,
	const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > &atlasCharts ,
	const GridAtlas< GeometryReal , MatrixReal > &gridAtlas ,
	const std::vector< AtlasInteriorOrBoundaryNodeIndex > &fineBoundaryIndex ,
	int numFineBoundaryNodes ,
	std::vector< MatrixReal >& deepMassCoefficients ,
	std::vector< MatrixReal >& deepStiffnessCoefficients ,
	SparseMatrix< MatrixReal , int >& boundaryBoundaryMassMatrix ,
	SparseMatrix< MatrixReal , int >& boundaryBoundaryStiffnessMatrix ,
	SparseMatrix< MatrixReal , int >& boundaryDeepMassMatrix ,
	SparseMatrix< MatrixReal , int >& boundaryDeepStiffnessMatrix ,
	bool computeDivergence ,
	Map< SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > , AtlasRefinedBoundaryEdgeIndex > & fineBoundaryEdgeIndex,
	Map< SimplexIndex< 1 , AtlasTexelIndex > , unsigned int > & coarseEdgeIndex,
	std::vector< Eigen::Triplet< MatrixReal > > & boundaryDeepDivergenceTriplets,
	std::vector< Eigen::Triplet< MatrixReal > > & boundaryBoundaryDivergenceTriplets,
	std::vector< MatrixReal >& deepDivergenceCoefficients
)
{
	auto MergeTriplets = [] ( const std::vector< std::vector< Eigen::Triplet< MatrixReal > > > &inTriplets , std::vector< Eigen::Triplet< MatrixReal > > &outTriplets )
		{
			for( unsigned int i=0 ; i<inTriplets.size() ; i++ ) outTriplets.insert( outTriplets.end() , inTriplets[i].begin() , inTriplets[i].end() );
		};
	auto MergeMassAndStiffnessCoefficients = [] ( const std::vector< std::vector< MassAndStiffnessCoefficient< MatrixReal > > > &in , std::vector< MassAndStiffnessCoefficient< MatrixReal > > &out )
		{
			for( unsigned int i=0 ; i<in.size() ; i++ ) out.insert( out.end() , in[i].begin() , in[i].end() );
		};

	const ExplicitIndexVector< ChartIndex , GridChart< GeometryReal > > &gridCharts = gridAtlas.gridCharts;
	const typename GridAtlas<>::IndexConverter & indexConverter = gridAtlas.indexConverter;

	std::vector< MassAndStiffnessCoefficient< MatrixReal > > boundaryBoundaryMassAndStiffness;
	std::vector< MassAndStiffnessCoefficient< MatrixReal > > boundaryDeepMassAndStiffness;

	std::vector< std::vector< MassAndStiffnessCoefficient< MatrixReal > > > _boundaryBoundaryMassAndStiffness( ThreadPool::NumThreads() );
	std::vector< std::vector< MassAndStiffnessCoefficient< MatrixReal > > > _boundaryDeepMassAndStiffness    ( ThreadPool::NumThreads() );
	std::vector< std::vector< Eigen::Triplet< MatrixReal > > > _boundaryDeepDivergenceTriplets            ( ThreadPool::NumThreads() );
	std::vector< std::vector< Eigen::Triplet< MatrixReal > > > _boundaryBoundaryDivergenceTriplets        ( ThreadPool::NumThreads() );

	ThreadPool::ParallelFor
	(
		0 , gridCharts.size() ,
		[&]( unsigned int thread , size_t i )
		{
		_InitializeChart< Samples , SanityCheck >
				(
					parameterMetric[ ChartIndex(i) ] ,
					atlasCharts[ ChartIndex(i) ] ,
					gridCharts[ ChartIndex(i) ] ,
					indexConverter ,
					fineBoundaryIndex ,
					deepMassCoefficients ,
					deepStiffnessCoefficients , 
					_boundaryBoundaryMassAndStiffness[thread] ,
					_boundaryDeepMassAndStiffness[thread] ,
					computeDivergence ,
					fineBoundaryEdgeIndex ,
					coarseEdgeIndex ,
					_boundaryDeepDivergenceTriplets[thread] ,
					_boundaryBoundaryDivergenceTriplets[thread] ,
					deepDivergenceCoefficients
				);
		}
	);

	MergeMassAndStiffnessCoefficients( _boundaryBoundaryMassAndStiffness , boundaryBoundaryMassAndStiffness );
	MergeMassAndStiffnessCoefficients( _boundaryDeepMassAndStiffness , boundaryDeepMassAndStiffness );
	MergeTriplets( _boundaryDeepDivergenceTriplets , boundaryDeepDivergenceTriplets );
	MergeTriplets( _boundaryBoundaryDivergenceTriplets , boundaryBoundaryDivergenceTriplets );

	SetSparseSymmetricMatrices( boundaryBoundaryMassAndStiffness , numFineBoundaryNodes , boundaryBoundaryMassMatrix , boundaryBoundaryStiffnessMatrix );
	SetSparseMatrices( boundaryDeepMassAndStiffness , static_cast< unsigned int >(gridAtlas.endBoundaryTexelIndex) , static_cast< unsigned int >(gridAtlas.endCombinedTexelIndex) , false , boundaryDeepMassMatrix , boundaryDeepStiffnessMatrix );
}

template< unsigned int Samples , bool SanityCheck , typename GeometryReal , typename MatrixReal >
void OperatorInitializer::_Initialize
(
	MassAndStiffnessOperators< MatrixReal > & massAndStiffnessOperators ,
	const GridAtlas< GeometryReal , MatrixReal > & gridAtlas ,
	const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > &parameterMetric ,
	const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > &atlasCharts ,
	const BoundaryProlongationData< MatrixReal > &boundaryProlongation ,
	bool computeDivergence ,
	DivergenceOperator< MatrixReal > & divergenceOperator
)
{
	massAndStiffnessOperators.indexConverter = gridAtlas.indexConverter;
	massAndStiffnessOperators.rasterLines = gridAtlas.rasterLines;

	//(2) Initialize mass and stiffness
	massAndStiffnessOperators.massCoefficients.deepCoefficients.resize( 10 * static_cast< unsigned int >( gridAtlas.endInteriorTexelIndex ) , 0 );
	massAndStiffnessOperators.stiffnessCoefficients.deepCoefficients.resize( 10 * static_cast< unsigned int >( gridAtlas.endInteriorTexelIndex ) , 0 );
	Map< SimplexIndex< 1 , AtlasTexelIndex > , unsigned int > edgeToIndex;
	if( computeDivergence )
	{
		InitializeIntraChartEdgeIndexing< SanityCheck >( gridAtlas.gridCharts , edgeToIndex );
		divergenceOperator.deepCoefficients.resize( 20 * static_cast< unsigned int >( gridAtlas.endInteriorTexelIndex ) , 0 );
	}

	SparseMatrix< MatrixReal , int > fineBoundaryBoundaryMassMatrix;
	SparseMatrix< MatrixReal , int > fineBoundaryBoundaryStiffnessMatrix;

	std::vector< Eigen::Triplet< MatrixReal > > boundaryDivergenceTriplets;
	std::vector< Eigen::Triplet< MatrixReal > > boundaryBoundaryDivergenceTriplets;
	Map< SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > , AtlasRefinedBoundaryEdgeIndex > fineBoundaryEdgeIndex;

	if( computeDivergence ) InitializeFineBoundaryEdgeIndexing( boundaryProlongation.fineBoundaryIndex , fineBoundaryEdgeIndex , gridAtlas.gridCharts );

	SparseMatrix< MatrixReal , int > fineBoundaryCellStiffnessRHSMatrix[3];
	std::vector< Point3D< MatrixReal > > fineBoundarySignal;

	_Initialize< Samples , SanityCheck >( parameterMetric , atlasCharts , gridAtlas , boundaryProlongation.fineBoundaryIndex , boundaryProlongation.numFineBoundaryNodes , massAndStiffnessOperators.massCoefficients.deepCoefficients , massAndStiffnessOperators.stiffnessCoefficients.deepCoefficients , fineBoundaryBoundaryMassMatrix , fineBoundaryBoundaryStiffnessMatrix , massAndStiffnessOperators.massCoefficients.boundaryDeepMatrix , massAndStiffnessOperators.stiffnessCoefficients.boundaryDeepMatrix , computeDivergence , fineBoundaryEdgeIndex , edgeToIndex , boundaryDivergenceTriplets , boundaryBoundaryDivergenceTriplets , divergenceOperator.deepCoefficients );

	{
		SparseMatrix< MatrixReal , int > temp = fineBoundaryBoundaryMassMatrix * boundaryProlongation.coarseBoundaryFineBoundaryProlongation;
		massAndStiffnessOperators.massCoefficients.boundaryBoundaryMatrix = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction * temp;
	}
	{
		SparseMatrix< MatrixReal , int > temp = fineBoundaryBoundaryStiffnessMatrix * boundaryProlongation.coarseBoundaryFineBoundaryProlongation;
		massAndStiffnessOperators.stiffnessCoefficients.boundaryBoundaryMatrix = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction * temp;
	}

	{
		std::vector< MatrixReal > in ( massAndStiffnessOperators.massCoefficients.boundaryBoundaryMatrix.Rows() , (MatrixReal)1. );
		std::vector< MatrixReal > out( massAndStiffnessOperators.massCoefficients.boundaryBoundaryMatrix.Rows() , (MatrixReal)0. );
		massAndStiffnessOperators.massCoefficients.boundaryBoundaryMatrix.Multiply( GetPointer(in) , GetPointer(out) );
		for( int i=0 ; i<out.size() ; i++ ) if( out[i]==0 )
			if( massAndStiffnessOperators.massCoefficients.boundaryBoundaryMatrix.RowSize(i)==0 ) MK_WARN( "Emptry row at boundary index " , i , ". Try running with jittering." );
			else                                                                                  MK_WARN( "Zero row at boundary index " , i , ". Try running with jittering." );
	}

	if( computeDivergence )
	{
		SparseMatrix< MatrixReal , int > fineBoundaryBoundaryDivergenceMatrix = SetSparseMatrix( boundaryBoundaryDivergenceTriplets , boundaryProlongation.numFineBoundaryNodes , (int)fineBoundaryEdgeIndex.size() , false );

		Map< SimplexIndex< 1 , AtlasInteriorOrBoundaryNodeIndex > , AtlasRefinedBoundaryEdgeIndex > boundaryCoarseEdgeIndex;
		std::vector< unsigned int > boundaryCoarseEdgeToGlobalEdge;

		InitializeBoundaryEdgeIndexing( massAndStiffnessOperators.massCoefficients.boundaryBoundaryMatrix , gridAtlas.indexConverter , edgeToIndex , boundaryCoarseEdgeToGlobalEdge , boundaryCoarseEdgeIndex );

		SparseMatrix< MatrixReal , int > boundaryCoarseToFineBoundaryOneFormProlongation;
		InitializeBoundaryCoarseToFineBoundaryOneFormProlongation< SanityCheck >( boundaryProlongation.coarseBoundaryFineBoundaryProlongation , boundaryCoarseEdgeIndex , fineBoundaryEdgeIndex , boundaryCoarseToFineBoundaryOneFormProlongation );

		SparseMatrix< MatrixReal , int > temp = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction * fineBoundaryBoundaryDivergenceMatrix;
		SparseMatrix< MatrixReal , int > boundaryBoundaryDivergenceMatrix = temp *  boundaryCoarseToFineBoundaryOneFormProlongation;
		const typename GridAtlas<>::IndexConverter & indexConverter = gridAtlas.indexConverter;
		for( int i=0 ; i<boundaryBoundaryDivergenceMatrix.Rows() ; i++ )
		{
			AtlasTexelIndex supportedIndex = indexConverter.boundaryToCombined( AtlasBoundaryTexelIndex(i) );
			for( int j=0 ; j<boundaryBoundaryDivergenceMatrix.RowSize(i) ; j++ )
				boundaryDivergenceTriplets.emplace_back( static_cast< unsigned int >(supportedIndex) , boundaryCoarseEdgeToGlobalEdge[ boundaryBoundaryDivergenceMatrix[i][j].N ] , boundaryBoundaryDivergenceMatrix[i][j].Value );
		}
		divergenceOperator.boundaryMatrix = SetSparseMatrix( boundaryDivergenceTriplets , static_cast< unsigned int >( gridAtlas.endCombinedTexelIndex ) , (int)edgeToIndex.size() , false );
		divergenceOperator.rasterLines = DivergenceOperator< MatrixReal >::DivergenceRasterLine::template GetRasterLines< SanityCheck >( edgeToIndex , gridAtlas.rasterLines );
		divergenceOperator.edges.resize( edgeToIndex.size() );
		for( auto edgeIter=edgeToIndex.begin() ; edgeIter!=edgeToIndex.end() ; edgeIter++ ) divergenceOperator.edges[ (*edgeIter).second ] = (*edgeIter).first;
	}
}

template< bool SanityCheck , typename GeometryReal , typename MatrixReal >
void OperatorInitializer::_Initialize
(
	unsigned int samples ,
	MassAndStiffnessOperators< MatrixReal > & massAndStiffnessOperators ,
	const GridAtlas< GeometryReal , MatrixReal > & gridAtlas ,
	const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > &parameterMetric ,
	const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > &atlasCharts ,
	const BoundaryProlongationData< MatrixReal > &boundaryProlongation ,
	bool computeDivergence ,
	DivergenceOperator< MatrixReal > & divergenceOperator
)
{
	switch( samples )
	{
	case  1: return _Initialize<  1 , SanityCheck >( massAndStiffnessOperators , gridAtlas , parameterMetric , atlasCharts , boundaryProlongation , computeDivergence , divergenceOperator );
	case  3: return _Initialize<  3 , SanityCheck >( massAndStiffnessOperators , gridAtlas , parameterMetric , atlasCharts , boundaryProlongation , computeDivergence , divergenceOperator );
	case  6: return _Initialize<  6 , SanityCheck >( massAndStiffnessOperators , gridAtlas , parameterMetric , atlasCharts , boundaryProlongation , computeDivergence , divergenceOperator );
	case 12: return _Initialize< 12 , SanityCheck >( massAndStiffnessOperators , gridAtlas , parameterMetric , atlasCharts , boundaryProlongation , computeDivergence , divergenceOperator );
	case 24: return _Initialize< 24 , SanityCheck >( massAndStiffnessOperators , gridAtlas , parameterMetric , atlasCharts , boundaryProlongation , computeDivergence , divergenceOperator );
	case 32: return _Initialize< 32 , SanityCheck >( massAndStiffnessOperators , gridAtlas , parameterMetric , atlasCharts , boundaryProlongation , computeDivergence , divergenceOperator );
	default: MK_THROW( "Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles: " , samples );
	}
}

template< typename GeometryReal , typename MatrixReal >
void OperatorInitializer::Initialize
(
	unsigned int samples ,
	MassAndStiffnessOperators< MatrixReal > & massAndStiffnessOperators ,
	const GridAtlas< GeometryReal , MatrixReal > & gridAtlas ,
	const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > &parameterMetric ,
	const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > &atlasCharts ,
	DivergenceOperator< MatrixReal > & divergenceOperator ,
	bool sanityCheck
)
{
	BoundaryProlongationData< MatrixReal > boundaryProlongation;
	InitializeBoundaryProlongationData( gridAtlas , boundaryProlongation , sanityCheck );
	if( sanityCheck ) _Initialize< true  >( samples , massAndStiffnessOperators , gridAtlas , parameterMetric , atlasCharts , boundaryProlongation , true , divergenceOperator );
	else              _Initialize< false >( samples , massAndStiffnessOperators , gridAtlas , parameterMetric , atlasCharts , boundaryProlongation , true , divergenceOperator );
}

template< typename GeometryReal , typename MatrixReal >
void OperatorInitializer::Initialize
(
	unsigned int samples ,
	MassAndStiffnessOperators< MatrixReal > & massAndStiffnessOperators ,
	const GridAtlas< GeometryReal , MatrixReal > & gridAtlas ,
	const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > &parameterMetric ,
	const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > &atlasCharts ,
	bool sanityCheck
)
{
	DivergenceOperator< MatrixReal > divergenceOperator;
	BoundaryProlongationData< MatrixReal > boundaryProlongation;
	InitializeBoundaryProlongationData( gridAtlas , boundaryProlongation , sanityCheck );
	if( sanityCheck ) _Initialize< true  >( samples , massAndStiffnessOperators , gridAtlas , parameterMetric , atlasCharts , boundaryProlongation , false , divergenceOperator );
	else              _Initialize< false >( samples , massAndStiffnessOperators , gridAtlas , parameterMetric , atlasCharts , boundaryProlongation , false , divergenceOperator );
}

template< typename GeometryReal , typename MatrixReal >
void OperatorInitializer::Initialize
(
	unsigned int samples ,
	MassAndStiffnessOperators< MatrixReal > & massAndStiffnessOperators ,
	const GridAtlas< GeometryReal , MatrixReal > & gridAtlas ,
	const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > &parameterMetric ,
	const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > &atlasCharts ,
	const BoundaryProlongationData< MatrixReal > &boundaryProlongation ,
	bool sanityCheck
)
{
	DivergenceOperator< MatrixReal > divergenceOperator;
	if( sanityCheck ) _Initialize< true  >( samples , massAndStiffnessOperators , gridAtlas , parameterMetric , atlasCharts , boundaryProlongation , false , divergenceOperator );
	else              _Initialize< false >( samples , massAndStiffnessOperators , gridAtlas , parameterMetric , atlasCharts , boundaryProlongation , false , divergenceOperator );
}

template< typename GeometryReal , typename MatrixReal >
void OperatorInitializer::Initialize
(
	unsigned int samples ,
	MassAndStiffnessOperators< MatrixReal > & massAndStiffnessOperators ,
	const GridAtlas< GeometryReal , MatrixReal > & gridAtlas ,
	const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > &parameterMetric ,
	const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > &atlasCharts ,
	const BoundaryProlongationData< MatrixReal > &boundaryProlongation ,
	DivergenceOperator< MatrixReal > & divergenceOperator ,
	bool sanityCheck
)
{
	if( sanityCheck ) _Initialize< true  >( samples , massAndStiffnessOperators , gridAtlas , parameterMetric , atlasCharts , boundaryProlongation , true , divergenceOperator );
	else              _Initialize< false >( samples , massAndStiffnessOperators , gridAtlas , parameterMetric , atlasCharts , boundaryProlongation , true , divergenceOperator );
}


template< typename GeometryReal , typename MatrixReal , typename SampleType >
void OperatorInitializer::Initialize
(
	unsigned int samples ,
	MassAndStiffnessOperators< MatrixReal > & massAndStiffnessOperators ,
	const GridAtlas< GeometryReal , MatrixReal > & gridAtlas ,
	const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > &parameterMetric ,
	const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > &atlasCharts ,
	Integrator< MatrixReal , SampleType > & integrator ,
	unsigned int numSamples ,
	bool approximate ,
	bool sanityCheck
)
{
	BoundaryProlongationData< MatrixReal > boundaryProlongation;
	DivergenceOperator< MatrixReal > divergenceOperator;
	InitializeBoundaryProlongationData( gridAtlas , boundaryProlongation , sanityCheck );
	if( sanityCheck ) _Initialize< true  >( samples , massAndStiffnessOperators , gridAtlas , parameterMetric , atlasCharts , boundaryProlongation , false , divergenceOperator );
	else              _Initialize< false >( samples , massAndStiffnessOperators , gridAtlas , parameterMetric , atlasCharts , boundaryProlongation , false , divergenceOperator );

	// Set integrator
	ExplicitIndexVector< AtlasInteriorCellIndex , std::pair< unsigned int , unsigned int > > interiorCellLineIndex;

	if( sanityCheck ) InitializeGridAtlasInteriorCellLines< true  >( gridAtlas.gridCharts , integrator.interiorCellLines , interiorCellLineIndex );
	else              InitializeGridAtlasInteriorCellLines< false >( gridAtlas.gridCharts , integrator.interiorCellLines , interiorCellLineIndex );
	if( interiorCellLineIndex.size()!=static_cast< unsigned int >( gridAtlas.endInteriorCellIndex) ) MK_THROW( "Inconsistent number of interior cells! Expected " , gridAtlas.endInteriorCellIndex , " . Result " , interiorCellLineIndex.size() , "." );

	integrator.indexConverter = gridAtlas.indexConverter;

	integrator.coarseBoundaryFineBoundaryProlongation = boundaryProlongation.coarseBoundaryFineBoundaryProlongation;
	integrator.fineBoundaryCoarseBoundaryRestriction = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction;
	integrator.samples.resize( integrator.interiorCellLines.size() );

	if( sanityCheck )
	{
		switch( numSamples )
		{
		case  1: InitializeIntegration<  1 , true >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case  3: InitializeIntegration<  3 , true >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case  6: InitializeIntegration<  6 , true >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case 12: InitializeIntegration< 12 , true >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case 24: InitializeIntegration< 24 , true >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case 32: InitializeIntegration< 32 , true >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		default: MK_THROW( "Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles: " , numSamples );
		}
	}
	else
	{
		switch( numSamples )
		{
		case  1: InitializeIntegration<  1 , false >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case  3: InitializeIntegration<  3 , false >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case  6: InitializeIntegration<  6 , false >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case 12: InitializeIntegration< 12 , false >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case 24: InitializeIntegration< 24 , false >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case 32: InitializeIntegration< 32 , false >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		default: MK_THROW( "Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles: " , numSamples );
		}
	}

	integrator.samples.sort();
}

template< typename GeometryReal , typename MatrixReal , typename SampleType >
void OperatorInitializer::Initialize
(
	unsigned int samples ,
	MassAndStiffnessOperators< MatrixReal > & massAndStiffnessOperators ,
	const GridAtlas< GeometryReal , MatrixReal > & gridAtlas ,
	const ExplicitIndexVector< ChartIndex , ExplicitIndexVector< ChartMeshTriangleIndex , SquareMatrix< GeometryReal , 2 > > > &parameterMetric ,
	const ExplicitIndexVector< ChartIndex , AtlasChart< GeometryReal > > &atlasCharts ,
	DivergenceOperator< MatrixReal > & divergenceOperator ,
	Integrator< MatrixReal , SampleType > & integrator ,
	unsigned int numSamples ,
	bool approximate ,
	bool sanityCheck
)
{
	BoundaryProlongationData< MatrixReal > boundaryProlongation;
	InitializeBoundaryProlongationData( gridAtlas , boundaryProlongation );
	_Initialize( samples , massAndStiffnessOperators , gridAtlas , parameterMetric , atlasCharts , boundaryProlongation , true , divergenceOperator );

	// Set integrator
	ExplicitIndexVector< AtlasInteriorCellIndex , std::pair< unsigned int , unsigned int > > interiorCellLineIndex;

	InitializeGridAtlasInteriorCellLines( gridAtlas.gridCharts , integrator.interiorCellLines , interiorCellLineIndex );
	if( interiorCellLineIndex.size()!=static_cast< unsigned int >( gridAtlas.endInteriorCellIndex) ) MK_THROW( "Inconsistent number of interior cells! Expected " , gridAtlas.endInteriorCellIndex , " . Result " , interiorCellLineIndex.size() , "." );

	integrator.indexConverter = gridAtlas.indexConverter;

	integrator.coarseBoundaryFineBoundaryProlongation = boundaryProlongation.coarseBoundaryFineBoundaryProlongation;
	integrator.fineBoundaryCoarseBoundaryRestriction = boundaryProlongation.fineBoundaryCoarseBoundaryRestriction;
	integrator.samples.resize( integrator.interiorCellLines.size() );

	if( sanityCheck )
	{
		switch( numSamples )
		{
		case  1: InitializeIntegration<  1 , true >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case  3: InitializeIntegration<  3 , true >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case  6: InitializeIntegration<  6 , true >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case 12: InitializeIntegration< 12 , true >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case 24: InitializeIntegration< 24 , true >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case 32: InitializeIntegration< 32 , true >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		default: MK_THROW( "Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles: " , numSamples );
		}
	}
	else
	{
		switch( numSamples )
		{
		case  1: InitializeIntegration<  1 , false >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case  3: InitializeIntegration<  3 , false >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case  6: InitializeIntegration<  6 , false >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case 12: InitializeIntegration< 12 , false >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case 24: InitializeIntegration< 24 , false >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		case 32: InitializeIntegration< 32 , false >( parameterMetric , atlasCharts , gridAtlas.gridCharts , interiorCellLineIndex , boundaryProlongation.fineBoundaryIndex , integrator.samples , approximate ) ; break;
		default: MK_THROW( "Only 1-, 3-, 6-, 12-, 24-, and 32-point quadrature supported for triangles: " , numSamples );
		}
	}

	integrator.samples.sort();
}

////////////////
// Integrator //
////////////////

template< typename Real , typename SampleType >
template< typename OutData , typename InData >
typename Integrator< Real , SampleType >::template Scratch< OutData , InData > Integrator< Real , SampleType >::getScratch( void ) const
{
	size_t numCoarseBoundaryNodes = fineBoundaryCoarseBoundaryRestriction.rows;
	size_t numFineBoundaryNodes = coarseBoundaryFineBoundaryProlongation.rows;

	Scratch< OutData , InData > scratch;
	scratch._coarseBoundaryPrimal.resize( numCoarseBoundaryNodes );
	scratch._coarseBoundaryDual  .resize( numCoarseBoundaryNodes );
	scratch._fineBoundaryPrimal  .resize( numFineBoundaryNodes );
	scratch._fineBoundaryDual    .resize( numFineBoundaryNodes );

	return scratch;
}

template< typename Real , typename SampleType >
template< typename OutData , typename InData , typename SampleFunction /* = std::function< OutData ( InData , SquareMatrix< Real , 2 > ) */ >
void Integrator< Real , SampleType >::operator()( const std::vector< InData > & primal , const SampleFunction & SF , std::vector< OutData > & dual , bool add ) const
{
	Scratch< OutData , InData > scratch = getScratch< OutData , InData >();
	return operator()( primal , SF , scratch , dual , add );
}

template< typename Real , typename SampleType >
template< typename OutData , typename InData , typename SampleFunction /* = std::function< OutData ( InData , SquareMatrix< Real , 2 > ) */ >
std::vector< OutData > Integrator< Real , SampleType >::operator()( const std::vector< InData > &primal , const SampleFunction & SF ) const
{
	Scratch< OutData , InData > scratch = getScratch< OutData , InData >();
	return operator()( primal , SF , scratch );
}

template< typename Real , typename SampleType >
template< typename OutData , typename InData , typename SampleFunction /* = std::function< OutData ( InData , SquareMatrix< Real , 2 > ) */ >
void Integrator< Real , SampleType >::operator()( ConstPointer( InData ) primal , const SampleFunction & SF , Pointer( OutData ) dual , bool add ) const
{
	Scratch< OutData , InData > scratch = getScratch< OutData , InData >();
	return operator()( primal , SF , scratch , dual , add );
}


template< typename Real , typename SampleType >
template< typename OutData , typename InData , typename SampleFunction /* = std::function< OutData ( InData , SquareMatrix< Real , 2 > ) */ >
void Integrator< Real , SampleType >::operator()( const std::vector< InData > & primal , const SampleFunction & SF , Scratch< OutData , InData > & scratch , std::vector< OutData > & dual , bool add ) const
{
	return operator()( GetPointer( primal ) , SF , scratch , GetPointer( dual ) , add );
}

template< typename Real , typename SampleType >
template< typename OutData , typename InData , typename SampleFunction /* = std::function< OutData ( InData , SquareMatrix< Real , 2 > ) */ >
std::vector< OutData > Integrator< Real , SampleType >::operator()( const std::vector< InData > &primal , const SampleFunction & SF , Scratch< OutData , InData > & scratch ) const
{
	std::vector< OutData > dual( indexConverter.numCombined() );
	operator()( primal , SF , scratch , dual );
	return dual;
}

template< typename Real , typename SampleType >
template< typename OutData , typename InData , typename SampleFunction /* = std::function< OutData ( InData , SquareMatrix< Real , 2 > ) */ >
void Integrator< Real , SampleType >::operator()( ConstPointer( InData ) primal , const SampleFunction & SF , Scratch< OutData , InData > & scratch , Pointer( OutData ) dual , bool add ) const
{
	// Extract the coarse (texel) boundary values
	ThreadPool::ParallelFor( 0 , indexConverter.numBoundary() , [&]( size_t i ){ scratch._coarseBoundaryPrimal[i] = primal[ static_cast< unsigned int >(indexConverter.boundaryToCombined( AtlasBoundaryTexelIndex(i) ) ) ]; } );

	// Prolong the coarse (texel) boundary values to the fine (node) boundary values
	coarseBoundaryFineBoundaryProlongation.Multiply( GetPointer( scratch._coarseBoundaryPrimal ) , GetPointer( scratch._fineBoundaryPrimal ) );

	// Clear the accumulators
	if( !add ) ThreadPool::ParallelFor( 0 , indexConverter.numCombined() , [&]( size_t i ){ dual[i] = OutData{}; } );
	ThreadPool::ParallelFor( 0 , scratch._fineBoundaryDual.size() , [&]( size_t i ){ scratch._fineBoundaryDual[i] = OutData{}; } );

	Integrate< Real >( interiorCellLines , samples , primal , GetPointer( scratch._fineBoundaryPrimal ) , SF , dual , GetPointer( scratch._fineBoundaryDual ) );

	// Restrict the fine (node) boundary values to the coarse (texel) boundary values
	fineBoundaryCoarseBoundaryRestriction.Multiply( GetPointer( scratch._fineBoundaryDual ) , GetPointer( scratch._coarseBoundaryDual ) );

	ThreadPool::ParallelFor( 0 , indexConverter.numBoundary() , [&]( size_t i ){ dual[ static_cast< unsigned int >( indexConverter.boundaryToCombined( AtlasBoundaryTexelIndex(i) ) ) ] += scratch._coarseBoundaryDual[i]; } );
}
