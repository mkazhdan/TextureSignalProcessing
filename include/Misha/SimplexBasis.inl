/*
Copyright (c) 2022, Michael Kazhdan
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

////////////////
// MultiIndex //
////////////////
template< unsigned int Size , typename Index , bool SmallestFirst >
template< typename ... UInts >
MultiIndex< Size , Index , SmallestFirst >::MultiIndex( UInts ... indices )
{
	static_assert( sizeof ... ( UInts )==Size , "[ERROR] Wrong number of indices" );
	const Index _indices[] = { (Index)indices... };
	_init( _indices );
}
template< unsigned int Size , typename Index , bool SmallestFirst >
bool MultiIndex< Size , Index , SmallestFirst >::operator < ( const MultiIndex &idx ) const
{
	for( unsigned int i=0 ; i<Size ; i++ )
	{
		if( _indices[i]<idx._indices[i] ) return true;
		else if( _indices[i]>idx._indices[i] ) return false;
	}
	return false;
};

template< unsigned int Size , typename Index , bool SmallestFirst >
bool MultiIndex< Size , Index , SmallestFirst >::operator == ( const MultiIndex &idx ) const
{
	return !( (*this)<idx ) && !( idx<(*this) );
}

template< unsigned int Size , typename Index , bool SmallestFirst >
void MultiIndex< Size , Index , SmallestFirst >::_init( const Index indices[] )
{
	memcpy( _indices , indices, sizeof(unsigned int) * Size );
	if( SmallestFirst ) std::sort( _indices , _indices + Size , []( unsigned int v1 , unsigned int v2 ){ return v1<v2; } );
	else                std::sort( _indices , _indices + Size , []( unsigned int v1 , unsigned int v2 ){ return v1>v2; } );
}

template< unsigned int Size , typename Index , bool SmallestFirst >
std::ostream &operator << ( std::ostream &os , const MultiIndex< Size , Index , SmallestFirst > &idx )
{
	os << "{ ";
	for( unsigned int d=0 ; d<Size ; d++ )
	{
		os << idx[d];
		if( d!=Size-1 ) os << " , ";
	}
	os << " }";
	return os;
}

////////////////////////
// RightSimplex::Data //
////////////////////////
template< unsigned int Dim >
RightSimplex< Dim >::Data::Data( void )
{
	for( unsigned int d=0 ; d<Dim ; d++ ) vertices[d+1][d] = 1;

	SimplexIndex< Dim , unsigned int > s;
	for( int k=0 ; k<=Dim ; k++ ) s[k] = k;
	for( int i=0 ; i<=Dim ; i++ )
	{
		faces[i] = s.face( i );
		if( (i&1)==0 ) std::swap( faces[i][0] , faces[i][1] );
	}
}

//////////////////
// RightSimplex //
//////////////////

template< unsigned int Dim >
const Point< double , Dim > &RightSimplex< Dim >::Vertex( unsigned int v )
{
	static Data data;
	return data.vertices[v];
}

template< unsigned int Dim >
const SimplexIndex< Dim-1 , unsigned int > &RightSimplex< Dim >::Face( unsigned int f )
{
	static Data data;
	return data.faces[f];
}

template< unsigned int Dim >
Matrix< double , Dim+1 , Dim > RightSimplex< Dim >::AffineTransform( Permutation< Dim+1 > p )
{
	Permutation< Dim+1 > pInv = p.inverse();
	// Define A to be the matrix taking barycentric coordinates to homogenous coordinates:
	// 	       | -1 -1 -1 ... -1 1 |
	// 	       |  1  0  0 ...  0 0 |
	// 	       |  0  1  0 ...  0 0 |
	// 	       |  0  0  1 ...  0 0 |
	// 	   A = |          .        |
	// 	       |           .       |
	// 	       |            .      |
	// 	       |  0  0  0      1 0 |
	// Then the matrix realizing the permutation is obtained by:
	// -- mapping from barycentric coordinates to homogenous coordinates,
	// -- permuting the homogenous coordinates,
	// -- mapping from the homogenous coordinates back to the barycentric coordinates.
	SquareMatrix< double , Dim+1 > A , A_inv;
	for( unsigned int i=0 ; i<Dim ; i++ ) A(i,0) = -1 , A(i,1+i) = 1;
	A(Dim,0) = 1;
	A = A.inverse() * p.template toMatrix< double >() * A;

	Matrix< double , Dim+1 , Dim > _A;
	for( unsigned int i=0 ; i<=Dim ; i++ ) for( unsigned int j=0 ; j<Dim ; j++ ) _A(i,j) = A(i,j);
	return _A;
}

template< unsigned int Dim >
template< unsigned int EmbeddingDimension >
SquareMatrix< double , Dim > RightSimplex< Dim >::Metric( const Simplex< double , EmbeddingDimension , Dim > &s )
{
	return Metric< EmbeddingDimension >( [&]( unsigned int idx ){ return s[idx]; } );
}

template< unsigned int Dim >
template< unsigned int EmbeddingDimension , typename PositionFunctor >
SquareMatrix< double , Dim > RightSimplex< Dim >::Metric( PositionFunctor positionFunctor )
{
	SquareMatrix< double , Dim > g;
	for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ )
		g(i,j) = Point< double , EmbeddingDimension >::Dot( positionFunctor(i+1)-positionFunctor(0) , positionFunctor(j+1)-positionFunctor(0) );
	return g;
}

template< unsigned int Dim >
template< unsigned int SubDim >
SquareMatrix< double , SubDim > RightSimplex< Dim >::RestrictedMetric( const SquareMatrix< double , Dim > &g , SimplexIndex< SubDim , unsigned int > &subSimplex )
{
	SquareMatrix< double , SubDim > _g;
	for( unsigned int i=0 ; i<SubDim ; i++ )
		for( unsigned int j=0 ; j<SubDim ; j++ )
		{
			/*
			_g(i,j) = < v[ subSimplex[i+1] ] - v[ subSimplex[0] ] , v[ subSimplex[j+1] ] - v[ subSimplex[0] ] >

			= < v[ subSimplex[i+1] ]-v[0] - v[ subSimplex[0] ]+v[0] , v[ subSimplex[j+1] ]-v[0] - v[ subSimplex[0] ]+v[0] >

			= < v[ subSimplex[i+1] ]-v[0] , v[ subSimplex[j+1] ]-v[0] >
			- < v[ subSimplex[i+1] ]-v[0] , v[ subSimplex[0]   ]-v[0] >
			- < v[ subSimplex[0]   ]-v[0] , v[ subSimplex[j+1] ]-v[0] >
			+ < v[ subSimplex[0]   ]-v[0] , v[ subSimplex[0]   ]-v[0] >

			= g( subSimplex[i+1]-1 , subSimplex[j+1]-1 )
			- g( subSimplex[i+1]-1 , subSimplex[0  ]-1 )
			- g( subSimplex[0  ]-1 , subSimplex[j+1]-1 )
			+ g( subSimplex[0  ]-1 , subSimplex[0  ]-1 )
			*/
			if( subSimplex[i+1] && subSimplex[j+1] ) _g(i,j) += g( subSimplex[i+1]-1 , subSimplex[j+1]-1 );
			if( subSimplex[i+1] && subSimplex[0]   ) _g(i,j) -= g( subSimplex[i+1]-1 , subSimplex[0]  -1 );
			if( subSimplex[j+1] && subSimplex[0]   ) _g(i,j) -= g( subSimplex[0]  -1 , subSimplex[j+1]-1 );
			if( subSimplex[0]                      ) _g(i,j) += g( subSimplex[0]  -1 , subSimplex[0]  -1 );
		}
	return _g;
}

template< unsigned int Dim >
Point< double , Dim > RightSimplex< Dim >::FaceNormal( unsigned int f , SquareMatrix< double , Dim > g )
{
	Point< double , Dim > v[ Dim>1 ? Dim-1 : 1 ];
	SimplexIndex< Dim-1 , unsigned int > s = Face(f);
	for( unsigned int k=0 ; k<Dim-1 ; k++ ) v[k] = Vertex( s[k+1] ) - Vertex( s[0] );
	Point< double , Dim > n = Point< double , Dim >::CrossProduct( v );
	n = g.inverse() * n;

	return n / sqrt( Point< double , Dim >::Dot( n , g * n ) );
}

template< unsigned int Dim >
template< unsigned int Degree , typename Real >
double RightSimplex< Dim >::Integral( const PFunction< Degree , Real > &P , SquareMatrix< double , Dim > g )
{
	return P.integrateUnitRightSimplex() * (Real)sqrt( g.determinant() );
}

template< unsigned int Dim >
template< unsigned int Degree , typename Real >
typename RightSimplex< Dim >::template PVectorField< (Degree>1) ? Degree-1 : 0 , Real > RightSimplex< Dim >::Gradient( const PFunction< Degree , Real > &P , SquareMatrix< double , Dim > g )
{
	PVectorField< (Degree>1) ? Degree-1 : 0 , Real > G , D;
	SquareMatrix< double , Dim > gInv = g.inverse();
	for( int k=0 ; k<Dim ; k++ ) D[k] = P.d( k );
	for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) G[i] += D[j] * gInv(j,i);
	return G;
}

template< unsigned int Dim >
template< unsigned int Degree , typename Real >
typename RightSimplex< Dim >::template PMatrixField< (Degree>2) ? Degree-2 : 0 , Real > RightSimplex< Dim >::Hessian( const PFunction< Degree , Real > &P , SquareMatrix< double , Dim > g )
{
	PMatrixField< (Degree>2) ? Degree-2 : 0 , Real > H , D2;
	SquareMatrix< double , Dim > gInv = g.inverse();
	for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) D2(i,j) = P.d(i).d(j);
	// H = gInv * D2 * gInv
	// H(i,j) = \sum_k gInv(k,j) * ( D2 * gInv )(i,k)
	//        = \sum_k gInv(k,j) * \sum_l D2(k,l) * gInv(i,l)
	//        = \sum_{k,l} gInv(k,j) * D2(k,l) * gInv(i,l)
	for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ )
		for( int k=0 ; k<Dim ; k++ ) for( int l=0 ; l<Dim ; l++ )
			H(i,j) += gInv(k,j) * D2(k,l) * gInv(i,l);
	return H;
}

template< unsigned int Dim >
template< unsigned int Degree , typename Real >
typename RightSimplex< Dim >::template PFunction< (Degree>1) ? Degree-1 : 0 , Real > RightSimplex< Dim >::Divergence( const PVectorField< Degree , Real > &V , SquareMatrix< double , Dim > g )
{
	RightSimplex< Dim >::PFunction< (Degree>1) ? Degree-1 : 0 , Real > D;
	for( int k=0 ; k<Dim ; k++ ) D += V[k].d(k);
	return D;
}

template< unsigned int Dim >
template< unsigned int Degree , typename Real >
typename RightSimplex< Dim >::template PFunction< Degree , Real > RightSimplex< Dim >::VectorFieldComponent( const PVectorField< Degree , Real > &V , Point< double , Dim > v , SquareMatrix< double , Dim > g )
{
	RightSimplex< Dim >::PFunction< Degree , Real > C;
	Point< double , Dim > g_v = g * v;
	g_v /= Point< double , Dim >::Dot( v , g * v );
	for( int k=0 ; k<Dim ; k++ ) C += V[k] * g_v[k];

	return C;
}

template< unsigned int Dim >
template< unsigned int Degree , typename Real >
typename RightSimplex< Dim >::template PFunction< Degree , Real > RightSimplex< Dim >::PushForward( const PFunction< Degree , Real > &P , const Simplex< double , Dim , Dim > &s )
{
	// (x_1,x_2,...,x_n) -> s[0] + ( s[1]-s[0] ) * x_1 + ... + ( s[n]-s[0] ) * x_n
	SquareMatrix< double , Dim+1 > _A;
	_A(Dim,Dim) = 1;
	for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) _A(i,j) = s[i+1][j] - s[0][j];
	for( int i=0 ; i<Dim ; i++ ) _A(Dim,i) = s[0][i];
	_A = _A.inverse();

	Matrix< double , Dim , Dim > A;
	Point< double , Dim > c;
	for( int k=0 ; k<Dim ; k++ ) c[k] = _A(Dim,k);
	for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) A(i,j) = _A(i,j);
	return P( A , c );
}

template< unsigned int Dim >
template< unsigned int Degree , typename Real >
typename RightSimplex< Dim >::template PVectorField< Degree , Real > RightSimplex< Dim >::PushForward( const PVectorField< Degree , Real > &V , const Simplex< double , Dim , Dim > &s )
{
	PVectorField< Degree , Real > W;
	// (x_1,x_2,...,x_n) -> s[0] + ( s[1]-s[0] ) * x_1 + ... + ( s[n]-s[0] ) * x_n
	SquareMatrix< double , Dim+1 > H;
	Matrix< double , Dim , Dim > A , A_inv;
	Point< double , Dim > c_inv;

	H(Dim,Dim) = 1;
	for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) H(i,j) = s[i+1][j] - s[0][j] , A(i,j) = s[i+1][j] - s[0][j];
	for( int i=0 ; i<Dim ; i++ ) H(Dim,i) = s[0][i];
	H = H.inverse();

	for( int k=0 ; k<Dim ; k++ ) c_inv[k] = H(Dim,k);
	for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) A_inv(i,j) = H(i,j);
	PVectorField< Degree , Real > _W;
	for( int k=0 ; k<Dim ; k++ ) _W[k] = V[k]( A_inv , c_inv );

	for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) W[i] += _W[j] * A(j,i);

	return W;
}

template< unsigned int Dim >
Point< double , Dim > RightSimplex< Dim >::PushForward( Point< double , Dim > p , bool tangentVector , const Simplex< double , Dim , Dim > &s )
{
	Point< double , Dim > q;
	for( int k=0 ; k<Dim ; k++ ) q += ( s[k+1] - s[0] ) * p[k];
	if( !tangentVector ) q += s[0];
	return q;
}

template< unsigned int Dim >
void RightSimplex< Dim >::GSOrthogonalize( Point< double , Dim > frame[Dim] , SquareMatrix< double , Dim > g )
{
	for( int i=0 ; i<Dim ; i++ )
	{
		for( int j=0 ; j<i ; j++ ) frame[i] -= frame[j] * Point< double , Dim >::Dot( frame[i] , g * frame[j] );
		frame[i] /= sqrt( Point< double , Dim >::Dot( frame[i] , g * frame[i] ) );
	}
}


///////////////////////////
// SimplexElements::Data //
///////////////////////////
template< unsigned int Dim , unsigned int Degree >
SimplexElements< Dim , Degree >::Data::Data( void )
{
	unsigned int indices[Degree];
	_initialize( indices , Dim );
	// Get the elements
	{
		Point< double , Dim > positions[ NodeNum ];
		for( unsigned int i=0 ; i<NodeNum ; i++ ) positions[i] = nodePosition( i );

		Point< double , NodeNum > values;
		for( unsigned int i=0 ; i<NodeNum ; i++ )
		{
			values[i] = 1;
			Point< double , Polynomial::Polynomial< Dim , Degree , double >::NumCoefficients > coefficients( Polynomial::Polynomial< Dim , Degree , double >::EvaluationMatrix( positions ).inverse() * values );
			elements[i] = Polynomial::Polynomial< Dim , Degree , double >( coefficients );
			values[i] = 0;
		}
	}

	// Get the differentials
	for( unsigned int i=0 ; i<NodeNum ; i++ ) for( int k=0 ; k<Dim ; k++ ) differentials[i][k] = elements[i].d( k );

	// Get the restriction of the differentials to the faces
	{
		// For each face, compute the affine map taking the (Dim-1)-dimensional face into the simplex
		Matrix< double , Dim , Dim > A[Dim+1];
		for( unsigned int f=0 ; f<=Dim ; f++ )
		{
			// To parameterize the face, use the indexing as returned by RightSimplex< Dim >::Face
			SimplexIndex< Dim-1 , unsigned int > faceIndex = RightSimplex< Dim >::Face(f);

			// Create a simplex whose first Dim vertices are the indices of the face and whose last index is the index opposite the face
			Simplex< double , Dim , Dim > face;
			for( unsigned int d=0 ; d<Dim ; d++ ) face[d] = RightSimplex< Dim >::Vertex( faceIndex[d] );
			face[Dim] = RightSimplex< Dim >::Vertex( f );

			// Create the affine transformation mapping the right triangle corresponding to the face into the right triangle corresponding to the simplex
			for( int j=0 ; j<Dim ; j++ )
			{
				for( int i=0 ; i<Dim-1 ; i++ ) A[f](i,j) = face[i+1][j] - face[0][j];
				A[f]( Dim-1 , j ) = face[0][j];
			}
		}
		for( unsigned int i=0 ; i<NodeNum ; i++ ) for( int k=0 ; k<Dim ; k++ ) for( int f=0 ; f<=Dim ; f++ )
			faceDifferentials[i][f][k] = differentials[i][k].template operator()< Dim >( A[f] );
	}
}

template< unsigned int Dim , unsigned int Degree >
template< unsigned int D >
void SimplexElements< Dim , Degree >::Data::_initialize( unsigned int indices[Degree] , unsigned int max )
{
	for( unsigned int i=0 ; i<=max ; i++ ) 
	{
		indices[ Degree-1-D ] = i;
		if constexpr( D==0 )
		{
			unsigned int idx = NodeIndex( indices );
			for( unsigned int d=0 ; d<Degree ; d++ ) nodeEndPoints[idx][d] = indices[d];
		}
		else _initialize< D-1 >( indices , i );
	}
}

template< unsigned int Dim , unsigned int Degree >
Point< double , Dim > SimplexElements< Dim , Degree >::Data::nodePosition( unsigned int idx ) const
{
	Point< double , Dim > p;
	for( int d=0 ; d<Degree ; d++ ) p += RightSimplex< Dim >::Vertex( nodeEndPoints[idx][d] );
	return p / Degree;
}

/////////////////////
// SimplexElements //
/////////////////////

template< unsigned int Dim , unsigned int Degree >
void SimplexElements< Dim , Degree >::SetElements( Polynomial::Polynomial< Dim , Degree , double > elements[ NodeNum ] )
{
	static Data data;
	for( unsigned int n=0 ; n<NodeNum ; n++ ) elements[n] = data.elements[n];
}

template< unsigned int Dim , unsigned int Degree >
double SimplexElements< Dim , Degree >::Volume( SquareMatrix< double , Dim > g )
{
	double v = sqrt( fabs( g.determinant() ) );
	for( unsigned int d=2 ; d<=Dim ; d++ ) v /= d;
	return v;
}


template< unsigned int Dim , unsigned int Degree >
typename SimplexElements< Dim , Degree >::SystemMatrix SimplexElements< Dim , Degree >::MassMatrix( SquareMatrix< double , Dim > g )
{
	static SquareMatrix< double , NodeNum > _M;
	static bool firstTime = true;
	if( firstTime )
	{
		Polynomial::Polynomial< Dim , Degree , double > elements[ NodeNum ];
		SetElements( elements );
		for( unsigned int i=0 ; i<NodeNum ; i++ ) for( unsigned int j=0 ; j<NodeNum ; j++ ) _M(i,j) = RightSimplex< Dim >::template Integral< 2*Degree , double >( elements[i] * elements[j] );
		firstTime = false;
	}
	return _M * sqrt( g.determinant() );
}

template< unsigned int Dim , unsigned int Degree >
typename SimplexElements< Dim , Degree >::SystemMatrix SimplexElements< Dim , Degree >::GradientSquareNormMatrix( SquareMatrix< double , Dim > g )
{
	static SquareMatrix< double , Dim > _S[NodeNum][NodeNum];
	static const int _Degree = Degree>1 ? Degree : 1;
	static bool firstTime = true;
	if( firstTime )
	{
		Point< Polynomial::Polynomial< Dim , _Degree-1 , double > , Dim > G[ NodeNum ];
		{
			Polynomial::Polynomial< Dim , Degree , double > elements[ NodeNum ];
			SetElements( elements );
			for( unsigned int i=0 ; i<NodeNum ; i++ ) G[i] = RightSimplex< Dim >::Gradient( elements[i] );
		}
		for( unsigned int i=0 ; i<NodeNum ; i++ ) for( unsigned int j=0 ; j<NodeNum ; j++ )
		{
			SquareMatrix< double , Dim > s;
			for( int k=0 ; k<Dim ; k++ ) for( int l=0 ; l<Dim ; l++ ) _S[i][j](k,l) = RightSimplex< Dim >::Integral( G[i][k] * G[j][l] );
		}

		firstTime = false;
	}

	SquareMatrix< double , Dim > gInv = g.inverse();
	SquareMatrix< double , NodeNum > S;
	for( unsigned int i=0 ; i<NodeNum ; i++ ) for( unsigned int j=0 ; j<NodeNum ; j++ ) S(i,j) = ( _S[i][j] * gInv ).trace();
	return S * sqrt( g.determinant() );
}

template< unsigned int Dim , unsigned int Degree >
typename SimplexElements< Dim , Degree >::SystemMatrix SimplexElements< Dim , Degree >::HessianSquareNormMatrix( SquareMatrix< double , Dim > g )
{
	SquareMatrix< double , NodeNum > H;
	Polynomial::Polynomial< Dim , Degree , double > elements[ NodeNum ];
	typename RightSimplex< Dim >::template PMatrixField< (Degree>2) ? Degree-2 : 0 , double > elementHessians[ NodeNum ];

	SetElements( elements );
	for( unsigned int i=0 ; i<NodeNum ; i++ ) elementHessians[i] = RightSimplex< Dim >::Hessian( elements[i] , g );
	for( unsigned int i=0 ; i<NodeNum ; i++ ) for( unsigned int j=0 ; j<NodeNum ; j++ )
		for( unsigned int k=0 ; k<Dim ; k++ ) for( unsigned int l=0 ; l<Dim ; l++ )
			for( unsigned int m=0 ; m<Dim ; m++ ) for( unsigned int n=0 ; n<Dim ; n++ ) 
				H(i,j) += RightSimplex< Dim >::Integral( elementHessians[i](k,l) * elementHessians[j](m,n) , g ) * g(k,m) * g(l,n);
	return H;
}

template< unsigned int Dim , unsigned int Degree >
typename SimplexElements< Dim , Degree >::SystemMatrix SimplexElements< Dim , Degree >::LaplacianSquareNormMatrix( SquareMatrix< double , Dim > g )
{
	SquareMatrix< double , NodeNum > B;
	static const int _Degree = Degree>2 ? Degree : 2;
	Polynomial::Polynomial< Dim , Degree , double > elements[ NodeNum ];
	Polynomial::Polynomial< Dim , (Degree>2) ? Degree-2 : 0 , double > elementLaplacians[ NodeNum ];

	SetElements( elements );
	for( unsigned int i=0 ; i<NodeNum ; i++ ) elementLaplacians[i] = RightSimplex< Dim >::Divergence( RightSimplex< Dim >::Gradient( elements[i] , g ) , g );

	for( unsigned int i=0 ; i<NodeNum ; i++ ) for( unsigned int j=0 ; j<NodeNum ; j++ ) B(i,j) = RightSimplex< Dim >::Integral( elementLaplacians[i] * elementLaplacians[j] , g );

	return B;
}

template< unsigned int Dim , unsigned int Degree > 
unsigned int SimplexElements< Dim , Degree >::__Choose( unsigned int D , unsigned int K )
{
	if( !K ) return 1;
	else if( !D ) return 0;
	else return ( __Choose( D-1 , K-1 ) * D ) / K;
}

template< unsigned int Dim , unsigned int Degree > 
unsigned int SimplexElements< Dim , Degree >::_Choose( unsigned int D , unsigned int K )
{
	if( K>D/2 ) return __Choose( D , D-K );
	else        return __Choose( D , K );
}

template< unsigned int Dim , unsigned int Degree > 
unsigned int SimplexElements< Dim , Degree >::NodeIndex( const unsigned int v[Degree] )
{
	unsigned int _v[Degree];
	memcpy( _v , v , sizeof(_v) );
	std::sort( _v , _v + Degree , []( unsigned int v1 , unsigned int v2 ){ return v1>v2; } );
	unsigned int idx = 0;
	for( int d=0 ; d<Degree ; d++ ) idx += _Choose( _v[d]+(Degree-1-d) , Degree-d );
	return idx;

}

template< unsigned int Dim , unsigned int Degree >
void SimplexElements< Dim , Degree >::FactorNodeIndex( unsigned int nodeIndex , unsigned int v[Degree] )
{
	static Data data;
	for( int d=0 ; d<Degree ; d++ ) v[d] = data.nodeEndPoints[ nodeIndex ][d];
}

template< unsigned int Dim , unsigned int Degree >
Point< double , Dim > SimplexElements< Dim , Degree >::NodePosition( unsigned int idx )
{
	static Data data;
	return data.nodePosition( idx );
}

template< unsigned int Dim , unsigned int Degree >
Matrix< Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > , SimplexElements< Dim , Degree >::NodeNum , Dim+1 > SimplexElements< Dim , Degree >::FaceGradientOrthogonalComponents( SquareMatrix< double , Dim > g )
{
	unsigned int indices[Dim+1];
	SimplexIndex< Dim , unsigned int > si;
	for( unsigned int d=0 ; d<=Dim ; d++ ) indices[d] = d;
	return FaceGradientOrthogonalComponents( g , indices );
}
template< unsigned int Dim , unsigned int Degree >
Matrix< Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > , SimplexElements< Dim , Degree >::NodeNum , Dim+1 > SimplexElements< Dim , Degree >::FaceGradientOrthogonalComponents( SquareMatrix< double , Dim > g , const unsigned int *indices )
{
	typedef Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > FaceD;
	Matrix< FaceD , NodeNum , Dim+1 > faceGradientOrthogonalComponents;

	SquareMatrix< double , Dim > gInv = g.inverse();
	for( unsigned int n=0 ; n<NodeNum ; n++ )
	{
		for( unsigned int f=0 ; f<=Dim ; f++ )
		{
			// Compute the gradient of the n-th basis function, restricted to face f
			FaceD faceGrad;
			{
				FaceD faceDifferential = FaceDifferential(n,f);
				for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) faceGrad[i] += faceDifferential[j] * gInv(j,i);
			}

			// Create the face
			SimplexIndex< Dim-1 , unsigned int > faceIndex = RightSimplex< Dim >::Face(f);

			// Compute the permutation sorting indices from smallest to largest
			Permutation< Dim > p( [&]( unsigned int i , unsigned int j ){ return indices[ faceIndex[i] ]<indices[ faceIndex[j] ]; } );
			Permutation< Dim > pInverse = p.inverse();

			Simplex< double , Dim , Dim > face;
			for( unsigned int d=0 ; d<Dim ; d++ ) face[d] = RightSimplex< Dim >::Vertex( faceIndex[ pInverse[d] ] );
			face[Dim] = RightSimplex< Dim >::Vertex( f );

			// Get the orientation of the face within the simplex
			bool evenParity = ( p.parity() & 1 )==0;

			// Compute an orthonormal frame
			Point< double , Dim > gFrame[Dim];
			{
				Point< double , Dim > frame[Dim];
				for( unsigned int d=0 ; d<Dim ; d++ ) frame[d] = face[d+1] - face[0];
				RightSimplex< Dim >::GSOrthogonalize( frame , g );
				if( evenParity ) frame[Dim-1] = -frame[Dim-1];
				for( unsigned int d=0 ; d<Dim ; d++ ) gFrame[d] = g * frame[d];
			}

			// Obtain the coefficients by projecting the gradient field onto the orthonormal frame
			for( unsigned int i=0 ; i<Dim ; i++ )
			{
				Polynomial::Polynomial< Dim-1 , Degree-1 , double > C;
				for( unsigned int j=0 ; j<Dim ; j++ ) C += ( faceGrad[j] * gFrame[i][j] );
				faceGradientOrthogonalComponents(n,f)[i] = C;
			}
		}
	}
	return faceGradientOrthogonalComponents;
}

template< unsigned int Dim , unsigned int Degree >
Matrix< Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > , SimplexElements< Dim , Degree >::NodeNum , Dim+1 > SimplexElements< Dim , Degree >::FaceGradients( SquareMatrix< double , Dim > g )
{
	unsigned int indices[Dim+1];
	SimplexIndex< Dim , unsigned int > si;
	for( unsigned int d=0 ; d<=Dim ; d++ ) indices[d] = d;
	return FaceGradients( g , indices );
}
template< unsigned int Dim , unsigned int Degree >
Matrix< Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > , SimplexElements< Dim , Degree >::NodeNum , Dim+1 > SimplexElements< Dim , Degree >::FaceGradients( SquareMatrix< double , Dim > g , const unsigned int *indices )
{
	typedef Point< Polynomial::Polynomial< Dim-1 , Degree-1 , double > , Dim > FaceD;
	Matrix< FaceD , NodeNum , Dim+1 > faceGradients;

	SquareMatrix< double , Dim > gInv = g.inverse();
	for( unsigned int n=0 ; n<NodeNum ; n++ ) for( unsigned int f=0 ; f<=Dim ; f++ )
	{
		// Compute the gradient of the n-th basis function, restricted to face f
		FaceD faceDifferential = FaceDifferential(n,f);
		for( int i=0 ; i<Dim ; i++ ) for( int j=0 ; j<Dim ; j++ ) faceGradients(n,f)[i] += faceDifferential[j] * gInv(j,i);
	}
	return faceGradients;
}

template< unsigned int Dim , unsigned int Degree >
Point< SquareMatrix< double , Dim-1 > , Dim+1 > SimplexElements< Dim , Degree >::FaceMetrics( SquareMatrix< double , Dim > g )
{
	Point< SquareMatrix< double , Dim-1 > , Dim+1 > faceMetrics;

	for( unsigned int f=0 ; f<=Dim ; f++ )
	{
		// Get the face
		SimplexIndex< Dim-1 , unsigned int > faceIndex = RightSimplex< Dim >::Face(f);

		// Get the restriction of the metric to the face
		faceMetrics[f] = RightSimplex< Dim >::RestrictedMetric( g , faceIndex );
	}
	return faceMetrics;
}

template< unsigned int Dim , unsigned int Degree >
Point< SquareMatrix< double , Dim-1 > , Dim+1 > SimplexElements< Dim , Degree >::FaceMetrics( SquareMatrix< double , Dim > g , const unsigned int *indices )
{
	Point< SquareMatrix< double , Dim-1 > , Dim+1 > faceMetrics = FaceMetrics( g );

	for( unsigned int f=0 ; f<=Dim ; f++ )
	{
		// Get the face
		SimplexIndex< Dim-1 , unsigned int > faceIndex = RightSimplex< Dim >::Face(f);
		Permutation< Dim > p( [&]( unsigned int i , unsigned int j ){ return indices[ faceIndex[i] ]<indices[ faceIndex[j] ]; } );
		Matrix< double , Dim , Dim-1 > A = RightSimplex< Dim-1 >::AffineTransform( p );
		SquareMatrix< double , Dim-1 > L;
		for( unsigned int i=0 ; i<Dim-1 ; i++ ) for( unsigned int j=0 ; j<Dim-1 ; j++ ) L(i,j) = A(i,j);
		faceMetrics[f] = L.transpose() * faceMetrics[f] * L;
	}
	return faceMetrics;
}
