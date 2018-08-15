/*
Copyright (c) 2015, Michael Kazhdan
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
#ifndef FEM_INCLUDED
#define FEM_INCLUDED

#include <string.h>
#include "SparseMatrix.h"
#include "Geometry.h"
#include "Array.h"
#include <algorithm>
namespace FEM
{
	//////////////
	// Elements //
	//////////////
	enum
	{
		ELEMENT_VERTEX ,
		ELEMENT_EDGE ,
		ELEMENT_TRIANGLE ,
		ELEMENT_COUNT
	};
	static const char* ElementNames[] = { "vertex" , "edge" , "triangle" };
	template< unsigned int Type > struct ElementInfo
	{
		static const unsigned int ElementsPerTriangle;
	};
	static void TestElementType( unsigned int ElementType , const char* header , bool forceFailure )
	{
		if( ElementType>=ELEMENT_COUNT ) fprintf( stderr , "[ERROR] %s: Unrecognized element type: 0 <= %d < %d\n" , header , ElementType , ELEMENT_COUNT ) , exit( 0 );
		else if( forceFailure ) fprintf( stderr , "[ERROR] %s: Element type unsupported: %d (%s)\n" , header , ElementType , ElementNames[ElementType] ) ,  exit( 0 );
	}

	//////////////////////////////
	// Types of Basis Functions //
	//////////////////////////////
	enum
	{
		BASIS_0_WHITNEY ,				// 0-form represented as the linear combination of hat functions (1 value per vertices)
		BASIS_1_CONFORMING ,			// 1-form represented as the linear combination of gradients and 90-degree rotated gradients of hat functions (2 values per vertex)
		BASIS_1_WHITNEY ,				// 1-form represented as the linear combination of the Whitney 1-form basis functions (1 value per edge)
		BASIS_1_TRIANGLE_CONSTANT ,		// 1-form represented as piecewise (per triangle) constant cotangent vectors (2 values per triangle)
		BASIS_2_WHITNEY ,				// 2-form represented as the linear combination of the Whitney  2-form basis functions (1 value per triangle)
		BASIS_2_VERTEX_CONSTANT ,		// 2-form represented as piecewise (around vertex) constant density fields (1 value per vertex)
		BASIS_COUNT
	};
	static const char* BasisNames[] = { "scalar whitney" , "vector conforming" , "vector whitney" , "vector triangle constant" , "density whitney" , "density vertex constant" };
	template< unsigned int Type > struct BasisInfo
	{
		static const unsigned int CoefficientsPerElement;
		static const unsigned int Dimension;
		static const bool Lumpable;
		static const unsigned int ElementType;
		static const unsigned int Coefficients;
		static const bool Singular;
	};
	static void TestBasisType( unsigned int BasisType , const char* header , bool forceFailure )
	{
		if( BasisType>=BASIS_COUNT ) fprintf( stderr , "[ERROR] %s: Unrecognized basis type: 0 <= %d < %d\n" , header , BasisType , BASIS_COUNT ) , exit( 0 );
		else if( forceFailure ) fprintf( stderr , "[ERROR] %s: Basis type unsupported: %d (%s)\n" , header , BasisType , BasisNames[BasisType] ) ,  exit( 0 );
	}

	///////////////////////////////
	// Tangent/Cotangent Vectors //
	///////////////////////////////
	enum
	{
		VECTOR_PRIMAL ,
		VECTOR_DUAL
	};
	// Defining Tangent/Cotangent vector classes
	template< class Real , int ID >
	struct Point2DWrapper : public Point2D< Real >
	{
		Point2DWrapper( void ) : Point2D< Real >() {;}
		Point2DWrapper( Real v1 , Real v2 ) : Point2D< Real >( v1 , v2 ) {;}
		Point2DWrapper( const Point2D< Real >& v ) : Point2D< Real >( v ) {;}
		Point2DWrapper  operator +  ( const Point2DWrapper& v ) const { return Point2D< Real >::operator + ( v ); }
		Point2DWrapper  operator -  ( const Point2DWrapper& v ) const { return Point2D< Real >::operator - ( v ); }
		Point2DWrapper& operator += ( const Point2DWrapper& v )       { Point2D< Real >::operator += ( v ) ; return *this; }
		Point2DWrapper& operator -= ( const Point2DWrapper& v )       { Point2D< Real >::operator -= ( v ) ; return *this; }
		Point2DWrapper  operator *  ( Real s ) const { return Point2D< Real >::operator * ( s ); }
		Point2DWrapper  operator /  ( Real s ) const { return Point2D< Real >::operator / ( s ); }
		Point2DWrapper& operator *= ( Real s )       { Point2D< Real >::operator *= ( s ) ; return *this; }
		Point2DWrapper& operator /= ( Real s )       { Point2D< Real >::operator /= ( s ) ; return *this; }
		Point2DWrapper  operator -  ( void ) const { return Point2D< Real >::operator - ( ); }
	};
	template< class Real > using   TangentVector = Point2DWrapper< Real , VECTOR_PRIMAL >;
	template< class Real > using CotangentVector = Point2DWrapper< Real , VECTOR_DUAL   >;

	// Operators for moving between the primal and dual spaces
	template< class Real >   TangentVector< Real > Dual( const SquareMatrix< Real , 2 >& tensor , CotangentVector< Real > v ){ return tensor.inverse() * v; }
	template< class Real > CotangentVector< Real > Dual( const SquareMatrix< Real , 2 >& tensor ,   TangentVector< Real > v ){ return tensor * v; }

	// Operators for measuring stuff
	template< class Real > Real Dot( const SquareMatrix< Real , 2 >& tensor ,   TangentVector< Real > v1 ,   TangentVector< Real > v2 ){ return Point2D< Real >::Dot( v1 , tensor*v2 ); }
	template< class Real > Real Dot( const SquareMatrix< Real , 2 >& tensor , CotangentVector< Real > v1 , CotangentVector< Real > v2 ){ return Point2D< Real >::Dot( v1 , tensor.inverse()*v2 ); }
	template< class Real > Real SquareLength( const SquareMatrix< Real , 2 >& tensor ,   TangentVector< Real > v ){ return Dot( tensor , v , v ); }
	template< class Real > Real SquareLength( const SquareMatrix< Real , 2 >& tensor , CotangentVector< Real > v ){ return Dot( tensor , v , v ); }
	template< class Real > Real Angle( const SquareMatrix< Real , 2 >& tensor ,   TangentVector< Real > v1 ,   TangentVector< Real > v2 );
	template< class Real > Real Angle( const SquareMatrix< Real , 2 >& tensor , CotangentVector< Real > v1 , CotangentVector< Real > v2 );
	template< class Real >   TangentVector< Real > Rotate90( const SquareMatrix< Real , 2 >& tensor ,   TangentVector< Real > v );
	template< class Real > CotangentVector< Real > Rotate90( const SquareMatrix< Real , 2 >& tensor , CotangentVector< Real > v );
	template< class Real > SquareMatrix< Real , 2 >   J( const SquareMatrix< Real , 2 >& tensor );
	template< class Real > SquareMatrix< Real , 2 > CoJ( const SquareMatrix< Real , 2 >& tensor );
	template< class Real > Real Area( const SquareMatrix< Real , 2 >& tensor , const Point2D< Real > triangle[3] );
	template< class Real > Real Area( const SquareMatrix< Real , 2 >& tensor , Point2D< Real > v1 , Point2D< Real > v2 , Point2D< Real > v3 );

	// Return the conformal/authalic metric that is closest to the given metric
	template< class Real > SquareMatrix< Real , 2 > MakeConformal( const SquareMatrix< Real , 2 >& sourceTensor , const SquareMatrix< Real , 2 >& targetTensor );
	template< class Real > SquareMatrix< Real , 2 > MakeAuthalic ( const SquareMatrix< Real , 2 >& sourceTensor , const SquareMatrix< Real , 2 >& targetGensor );

	// Compute a symmetric square root of the metric tensor
	template< class Real > SquareMatrix< Real , 2 > TensorRoot( const SquareMatrix< Real , 2 >& tensor );

	///////////////////////////////////////////
	// Right Triangle (with a metric tensor) //
	///////////////////////////////////////////
	template< class Real >
	struct RightTriangle
	{
		// (Some of) the different centers one can define for a triangle
		enum
		{
			CENTER_BARYCENTRIC ,
			CENTER_CIRCUMCENTRIC ,
			CENTER_INCENTRIC ,
			CENTER_ISOGONIC ,
			CENTER_COUNT
		};
		static const char* CenterNames[];

		// Coordinates of the corners/edge-centers/edge-directions
		// The vertices are ordered { (0,0) , (1,0) , (0,1) }.
		// The edges are ordered so that the i-th edge is opposite the i-th vertex (w/ ccw orientation)
		static const Point2D< Real > Corners[];
		static const Point2D< Real > EdgeMidpoints[];
		static const CotangentVector< Real > CornerGradients[];
		static const   TangentVector< Real > EdgeDirections[];

		// Evaluate the scalar/vector/density fields represented by the coefficients.
		// Note that the tensor is required for the conforming basis to define the 90-degree rotation and for the 2-forms to remove the volume
		template< unsigned int BasisType , class V > static V                    EvaluateScalarField        (                                          ConstPointer( V ) coefficients , const Point2D< Real >& position );
		template< unsigned int BasisType , class V > static CotangentVector< V > EvaluateScalarFieldGradient(                                          ConstPointer( V ) coefficients , const Point2D< Real >& position );
		template< unsigned int BasisType , class V > static CotangentVector< V > EvaluateCovectorField      ( const SquareMatrix< Real , 2 >& tensor , ConstPointer( V ) coefficients , const Point2D< Real >& position );
		template< unsigned int BasisType , class V > static V                    EvaluateDensityField       ( const SquareMatrix< Real , 2 >& tensor , ConstPointer( V ) coefficients , const Point2D< Real >& position );

		// Compute the (possible lumped/weighted) mass matrix
		template< unsigned int BasisType > static SquareMatrix< Real , BasisInfo< BasisType >::Coefficients >  GetMassMatrix( const SquareMatrix< Real , 2 >& tensor );
		template< unsigned int BasisType > static SquareMatrix< Real , BasisInfo< BasisType >::Coefficients >  GetMassMatrix( const SquareMatrix< Real , 2 >& tensor , const SquareMatrix< Real , 2 >& newTensor );
		template< unsigned int BasisType > static Point< Real , BasisInfo< BasisType >::Coefficients > GetDiagonalMassMatrix( const SquareMatrix< Real , 2 >& tensor );

		// Compute the differential operator
		template< unsigned int InBasisType , unsigned int OutBasisType > static Matrix< Real , BasisInfo< InBasisType >::Coefficients , BasisInfo< OutBasisType >::Coefficients > GetDMatrix( const SquareMatrix< Real , 2 >& tensor );

		// Compute the masks that indicate which entries could be non-zero
		template< unsigned int BasisType > static SquareMatrix< unsigned char , BasisInfo< BasisType >::Coefficients > GetMassMask( bool useTensor );
		template< unsigned int InBasisType , unsigned int OutBasisType > static Matrix< unsigned char , BasisInfo< InBasisType >::Coefficients , BasisInfo< OutBasisType >::Coefficients > GetDMask( bool& redundant );

		// Compute the integrals
		template< unsigned int BasisType > static Real Integrate( const SquareMatrix< Real , 2 >& tensor , ConstPointer( Real ) linear );
		template< unsigned int BasisType > static Point< Real , BasisInfo< BasisType >::Coefficients > IntegrationDual( const SquareMatrix< Real , 2 >& tensor , ConstPointer( Real ) linear );
		template< unsigned int BasisType > static Point< Real , BasisInfo< BasisType >::Coefficients > IntegrationDual( const SquareMatrix< Real , 2 >& tensor , ConstPointer( CotangentVector< Real > ) linear );
		template< unsigned int BasisType > static Point< Real , BasisInfo< BasisType >::Coefficients > IntegrationDual( const SquareMatrix< Real , 2 >& tensor , Point2D< Real > p , Point2D< Real > q );

		static Point2D< Real > Center( const SquareMatrix< Real , 2 >& tensor , int centerType );
		static Point3D< Real > CenterAreas( const SquareMatrix< Real , 2 >& tensor , int centerType );

		static Point3D< Real > SubTriangleAreas( const SquareMatrix< Real , 2 >& tensor , Point2D< Real > center );
	};

	// Position on a mesh
	template< class Real >
	struct SamplePoint
	{
		int tIdx;
		Point2D< Real > p;
		SamplePoint( void ){ ; }
		SamplePoint( int tIdx , const Point2D< Real > p ){ this->tIdx = tIdx , this->p = p; }
	};
	// Position and tangent direction on a mesh
	template< class Real >
	struct HermiteSamplePoint : public SamplePoint< Real >
	{
		using SamplePoint< Real >::tIdx;
		using SamplePoint< Real >::p;
		TangentVector< Real > v;
		HermiteSamplePoint( void ){ ; }
		HermiteSamplePoint( int tIdx , const Point2D< Real >& p , const TangentVector< Real > v=TangentVector< Real >() ){ this->tIdx = tIdx , this->p = p , this->v = v; }
		HermiteSamplePoint( const SamplePoint< Real >& p , const Point2D< Real > v=TangentVector< Real >() ){ tIdx = p.tIdx , this->p = p.p , this->v = v; }
	};
	// Transformations taking points/directions in one coordinate frame to points/directions in the other
	template< class Real >
	struct CoordinateXForm : public Group< CoordinateXForm< Real > >
	{
		SquareMatrix< Real , 2 > linear;
		Point2D< Real > constant;
		CoordinateXForm( void ){ constant = Point2D< Real >() , linear = SquareMatrix< Real , 2 >::Identity(); }
		// (A,s) * p = A*p + s
		Point2D< Real > operator() ( const Point2D< Real >& p ) const { return linear*p + constant; }
		Point2D< Real > operator * ( const Point2D< Real >& p ) const { return linear*p + constant; }

		// (A,s) * (B,t) * p = (A,s) * (B*p + t) = A*B*p + (A*t + s)
		void SetIdentity( void ){ constant = Point2D< Real >() , linear = SquareMatrix< Real , 2 >::Identity(); }
		void Multiply( const CoordinateXForm& xForm ){ constant += linear * xForm.constant , linear *= xForm.linear; }
		void Invert( void ){ linear = linear.inverse() , constant = - linear * constant; }
	};

	// An abstract class that can be sampled at a point on the mesh and returns a tangent vector
	template< class Real >
	struct TangentVectorField{ virtual TangentVector< Real > operator() ( const SamplePoint< Real >& p ) const = 0; };
	// A realization of the abstract class by a vector field that is constant per triangle.
	template< class Real >
	struct SampledTangentVectorField : public TangentVectorField< Real >
	{
		SampledTangentVectorField( ConstPointer( TangentVector< Real > ) samples ) : _samples( samples ){ ; }
		TangentVector< Real > operator()( const SamplePoint< Real >& p ) const { return _samples[p.tIdx]; }
	protected:
		ConstPointer( TangentVector< Real > ) _samples;
	};

	// A structure that maps from half-edges to edges and vice-versa
	// -- Given an edge index, e, the value returned by operator[e] is an array of integers corresponding to the two half-edges incident on the edge.
	//    (The orientation is given by the first and the second is -1 if the edge is on the boundary).
	// -- Given a half edge index, he, the value of edge( he ) is the associated edge.
	//    (To find if the half-edge is aligned with the orientation of the edge, a reference to a boolean can be passed in, or the isAligned method can be called.)
	// -- Given a half edge index, he, the value of opposite( he ) is the opposite half edge index.
	//    (Or -1 if the he is on the boundary.)
	struct EdgeMap
	{
		EdgeMap( ConstPointer( TriangleIndex ) triangles , size_t tCount );
		~EdgeMap( void ){ FreePointer( _e2he ) ; FreePointer( _he2e ); }
		unsigned int size( void ) const { return _eCount; }
		// Given an edge index, this returns a pointer to the two half-edges adjacent on that edge. (The first is positively aligned with the edge.)
		const int* operator[]( unsigned int eIdx ) const { return _e2he + eIdx*2; }
		// Given a half-edge, this returns the index of the associated edge and indicates whether or not the alignment is consistent
		int edge( unsigned int he , bool& aligned ) const { if( _he2e[he]<0 ){ aligned = false ; return -_he2e[he]-1; } else { aligned = true ; return _he2e[he]-1; } }
		int edge( unsigned int he ) const { return std::max< int >( _he2e[he] , -_he2e[he] )-1; }
		bool isAligned( unsigned int he ) const { return _he2e[he]>0; }
		int opposite( unsigned int he ) const { int e = edge(he) ; return he==_e2he[ 2*e+0 ] ? _e2he[ 2*e+1 ] : _e2he[ 2*e+0 ]; }
	protected:
		unsigned int _eCount;
		Pointer( int ) _e2he;
		Pointer( int ) _he2e;
		template< typename Real > friend struct RiemannianMesh;
	};
	// This structure represents a Riemmanian mesh, with the triangles giving the connectivity and the square (symmetric) matrices giving the metric
	template< class Real >
	struct RiemannianMesh
	{
	protected:
		Pointer( TriangleIndex ) _triangles;
		Pointer( SquareMatrix< Real , 2 > ) _g;
		size_t _tCount , _vCount;
		EdgeMap _edgeMap;
	public:
		ConstPointer( TriangleIndex ) triangles( void     ) const { return _triangles   ; }
		const        TriangleIndex&   triangles( size_t t ) const { return _triangles[t]; }
		size_t tCount( void ) const { return _tCount; }
		size_t vCount( void ) const { return _vCount; }
		size_t eCount( void ) const { return _edgeMap.size(); }
		const SquareMatrix< Real , 2 >& g( size_t t ) const { return _g[t]; }
		int opposite( unsigned int he ) const { return _edgeMap.opposite( he ); }
		RiemannianMesh( Pointer( TriangleIndex ) t , size_t tC );
		~RiemannianMesh( void );

		template< class Vertex > void setMetricFromEmbedding( ConstPointer( Vertex ) vertices );
		void setMetricFromEdgeLengths( ConstPointer( Real ) edgeLengths );
		void setMetricFromSquareEdgeLengths( ConstPointer( Real ) squareEdgeLengths );
		void makeUnitArea( void );
		Real area( void ) const;
		Real area( int tIdx ) const;
		Real squareEdgeLength( int heIdx ) const;


		// Get the associated index, given the triangle and the index of the element within the triangle
		template< unsigned int BasisType > size_t dimension( void ) const;
		template< unsigned int BasisType > unsigned int index( unsigned int t , unsigned int idx ) const;
		template< unsigned int BasisType > unsigned int index( unsigned int t , unsigned int idx , bool& isAligned ) const;

		template< unsigned int BasisType , class V >                  V   evaluateScalarField        ( ConstPointer( V ) coefficients , const SamplePoint< Real >& p ) const;
		template< unsigned int BasisType , class V > CotangentVector< V > evaluateScalarFieldGradient( ConstPointer( V ) coefficients , const SamplePoint< Real >& p ) const;
		template< unsigned int BasisType , class V > CotangentVector< V > evaluateCovectorField      ( ConstPointer( V ) coefficients , const SamplePoint< Real >& p ) const;

		/////////////////////////
		// Geometric Operators //
		/////////////////////////
		template< unsigned int BasisType > SparseMatrix< Real , int > massMatrix( bool lump=false , ConstPointer( SquareMatrix< Real , 2 > ) newTensors = NullPointer< SquareMatrix< Real , 2 > >() ) const;

		// Integrate the piecewise linear function over the mesh
		Real getIntegral( ConstPointer( Real ) coefficients ) const;

		CoordinateXForm< Real >  exp( ConstPointer( CoordinateXForm< Real > ) xForms , HermiteSamplePoint< Real >& p , Real eps=(Real)0 ) const;
		CoordinateXForm< Real > flow( ConstPointer( CoordinateXForm< Real > ) xForms , const TangentVectorField< Real >& vf , Real flowTime , SamplePoint< Real >& p ,                                  Real minStepSize , Real eps=(Real)0 , std::vector< SamplePoint< Real > >* path=NULL ) const;
		Real                    flow( ConstPointer( CoordinateXForm< Real > ) xForms , const TangentVectorField< Real >& vf , Real flowTime , SamplePoint< Real >& p , CoordinateXForm< Real >& xForm , Real minStepSize , Real eps=(Real)0 ) const;

		Pointer( CoordinateXForm< Real > ) getCoordinateXForms( void ) const;
		void                               setCoordinateXForms( Pointer( CoordinateXForm< Real > ) xForms ) const;
		CoordinateXForm< Real > getVertexCoordinateXForm( ConstPointer( CoordinateXForm< Real > ) xForms , int t , int v ) const;
		CoordinateXForm< Real > xForm( unsigned int he ) const;
		std::vector< unsigned int > getVertexCorners( int t , int v ) const;
		Real getVertexConeAngle( int t , int v ) const;

		bool edgeFlip( int e , Real eps=0 );
		bool isVoronoiEdge( unsigned int e , Real eps=0 ) const;
	};

	template< class Real , unsigned int BasisType >
	struct TangentVectorFieldWrapper : public TangentVectorField< Real >
	{
		TangentVectorFieldWrapper( const RiemannianMesh< Real >* mesh , ConstPointer( Real ) coefficients , bool precomputeInverses=false );
		~TangentVectorFieldWrapper( void );
		TangentVector< Real > operator() ( const SamplePoint< Real >& p ) const;
	protected:
		const RiemannianMesh< Real >* _mesh;
		Pointer( SquareMatrix< Real , 2 > ) _gInverse;
		ConstPointer( Real ) _coefficients;
	};
}
template< class Real > const char* FEM::RightTriangle< Real >::CenterNames[] = { "barycentric" , "circumcentric" , "incentric" , "isogonic" };

template< class Real > const Point2D< Real > FEM::RightTriangle< Real >::Corners        [] = { Point2D< Real >(0,0) , Point2D< Real >(1,0) , Point2D< Real >(0,1) };
template< class Real > const Point2D< Real > FEM::RightTriangle< Real >::EdgeMidpoints  [] = { Point2D< Real >((Real)0.5,(Real)0.5) , Point2D< Real >(0,(Real)0.5) , Point2D< Real >((Real)0.5,0) };
template< class Real > const FEM::CotangentVector< Real > FEM::RightTriangle< Real >::CornerGradients[] = { FEM::CotangentVector< Real >( Point2D< Real >(-1,-1) ) , FEM::CotangentVector< Real >( Point2D< Real >( 1, 0) ) , FEM::CotangentVector< Real >( Point2D< Real >( 0,1) ) };
template< class Real > const FEM::  TangentVector< Real > FEM::RightTriangle< Real >::EdgeDirections [] = { FEM::  TangentVector< Real >( Point2D< Real >(-1, 1) ) , FEM::  TangentVector< Real >( Point2D< Real >( 0,-1) ) , FEM::  TangentVector< Real >( Point2D< Real >( 1,0) ) };

template<> const unsigned int FEM::ElementInfo< FEM::ELEMENT_VERTEX   >::ElementsPerTriangle = 3;
template<> const unsigned int FEM::ElementInfo< FEM::ELEMENT_EDGE     >::ElementsPerTriangle = 3;
template<> const unsigned int FEM::ElementInfo< FEM::ELEMENT_TRIANGLE >::ElementsPerTriangle = 1;

template<> const unsigned int FEM::BasisInfo< FEM::BASIS_0_WHITNEY           >::CoefficientsPerElement = 1;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_1_CONFORMING        >::CoefficientsPerElement = 2;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_1_WHITNEY           >::CoefficientsPerElement = 1;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_1_TRIANGLE_CONSTANT >::CoefficientsPerElement = 2;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_2_WHITNEY           >::CoefficientsPerElement = 1;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_2_VERTEX_CONSTANT   >::CoefficientsPerElement = 1;

template<> const unsigned int FEM::BasisInfo< FEM::BASIS_0_WHITNEY           >::Coefficients = 3;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_1_CONFORMING        >::Coefficients = 6;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_1_WHITNEY           >::Coefficients = 3;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_1_TRIANGLE_CONSTANT >::Coefficients = 2;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_2_WHITNEY           >::Coefficients = 1;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_2_VERTEX_CONSTANT   >::Coefficients = 3;

template<> const unsigned int FEM::BasisInfo< FEM::BASIS_0_WHITNEY           >::Dimension = 0;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_1_CONFORMING        >::Dimension = 1;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_1_WHITNEY           >::Dimension = 1;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_1_TRIANGLE_CONSTANT >::Dimension = 1;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_2_WHITNEY           >::Dimension = 2;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_2_VERTEX_CONSTANT   >::Dimension = 2;

template<> const bool FEM::BasisInfo< FEM::BASIS_0_WHITNEY           >::Lumpable = true;
template<> const bool FEM::BasisInfo< FEM::BASIS_1_CONFORMING        >::Lumpable = false;
template<> const bool FEM::BasisInfo< FEM::BASIS_1_WHITNEY           >::Lumpable = true;
template<> const bool FEM::BasisInfo< FEM::BASIS_1_TRIANGLE_CONSTANT >::Lumpable = false;
template<> const bool FEM::BasisInfo< FEM::BASIS_2_WHITNEY           >::Lumpable = true;
template<> const bool FEM::BasisInfo< FEM::BASIS_2_VERTEX_CONSTANT   >::Lumpable = true;

template<> const unsigned int FEM::BasisInfo< FEM::BASIS_0_WHITNEY           >::ElementType = FEM::ELEMENT_VERTEX;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_1_CONFORMING        >::ElementType = FEM::ELEMENT_VERTEX;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_1_WHITNEY           >::ElementType = FEM::ELEMENT_EDGE;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_1_TRIANGLE_CONSTANT >::ElementType = FEM::ELEMENT_TRIANGLE;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_2_WHITNEY           >::ElementType = FEM::ELEMENT_TRIANGLE;
template<> const unsigned int FEM::BasisInfo< FEM::BASIS_2_VERTEX_CONSTANT   >::ElementType = FEM::ELEMENT_VERTEX;

template<> const bool FEM::BasisInfo< FEM::BASIS_0_WHITNEY           >::Singular = false;
template<> const bool FEM::BasisInfo< FEM::BASIS_1_CONFORMING        >::Singular = true;
template<> const bool FEM::BasisInfo< FEM::BASIS_1_WHITNEY           >::Singular = false;
template<> const bool FEM::BasisInfo< FEM::BASIS_1_TRIANGLE_CONSTANT >::Singular = false;
template<> const bool FEM::BasisInfo< FEM::BASIS_2_WHITNEY           >::Singular = false;
template<> const bool FEM::BasisInfo< FEM::BASIS_2_VERTEX_CONSTANT   >::Singular = false;

#include "FEM.inl"
#endif // FEM_INCLUDED