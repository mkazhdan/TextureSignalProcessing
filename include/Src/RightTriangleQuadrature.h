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
#ifndef RIGHT_TRIANGLE_CUADRATURE
#define RIGHT_TRIANGLE_CUADRATURE

#include <Misha/Geometry.h>

template< unsigned int Samples >
struct SegmentIntegrator
{
	static const double Weights[];
	static const double Positions[];
};
template<> const double SegmentIntegrator<1>::Weights[] = { 1. };
template<> const double SegmentIntegrator<2>::Weights[] = { 0.5 , 0.5 };
template<> const double SegmentIntegrator<3>::Weights[] = { 5./18 , 8./18 , 5./18 };
template<> const double SegmentIntegrator<4>::Weights[] = { ( 18.-sqrt(30) ) / 72 , ( 18.+sqrt(30) ) / 72 , ( 18.+sqrt(30) ) / 72 , ( 18.-sqrt(30) ) / 72 };
template<> const double SegmentIntegrator<5>::Weights[] = { ( 322. - 13.*sqrt(70.) )/900. , ( 322. + 13.*sqrt(70.) )/900. , 128./225 , ( 322. + 13.*sqrt(70.) )/900. , ( 322. - 13.*sqrt(70.) )/900. };
template<> const double SegmentIntegrator<1>::Positions[] = {                      1./2                      };
template<> const double SegmentIntegrator<2>::Positions[] = { 1./2 - sqrt(1./12) ,        1./2 + sqrt(1./12) };
template<> const double SegmentIntegrator<3>::Positions[] = { 1./2 - sqrt(3./20) , 1./2 , 1./2 + sqrt(3./20) };
template<> const double SegmentIntegrator<4>::Positions[] = { 1./2 - sqrt( 3./28 + 2./28 * sqrt(6./5) ) , 1./2 - sqrt( 3./28 - 2./28 * sqrt(6./5) ) , 1./2 + sqrt( 3./28 - 2./28 * sqrt(6./5) ) , 1./2 + sqrt( 3./28 + 2./28 * sqrt(6./5) ) };
template<> const double SegmentIntegrator<5>::Positions[] = { 1./2 - 1./6 * sqrt( 5. + 2. * sqrt(10./7) ) , 1./2 - 1./6 * sqrt( 5. - 2. * sqrt(10./7) ) , 1./2 , 1./2 + 1./6 * sqrt( 5. - 2. * sqrt(10./7) ) , 1./2 + 1./6 * sqrt( 5. + 2. * sqrt(10./7) ) };

template< unsigned int Samples >
struct QuadIntegrator
{
	QuadIntegrator( const Point2D< double > corners[2][2] )
	{
		// P(x,y) = c(0,0)*(1-x)*(1-y) + c(0,1)*(1-x)*y + c(1,1)*x*y + c(1,0)*x*(1-y)
		// dP(x,y) = [ -c(0,0)*(1-y) - c(0,1)*y + c(1,1)*y + c(1,0)*(1-y) , -c(0,0)*(1-x) + c(0,1)*(1-x) + c(1,1)*x  - c(1,0)*x ]
		for( int i=0 ; i<Samples ; i++ ) for( int j=0 ; j<Samples ; j++ )
		{
			const double x = SegmentIntegrator< Samples >::Positions[i] , y = SegmentIntegrator< Samples >::Positions[j];
			const double d[2][2] =
			{
				{ ( corners[1][0][0] - corners[0][0][0] )*(1.-y) + ( corners[1][1][0] - corners[0][1][0] )*y , ( corners[1][0][1] - corners[0][0][1] )*(1.-y) + ( corners[1][1][1] - corners[0][1][1] )*y } ,
				{ ( corners[0][1][0] - corners[0][0][0] )*(1.-x) + ( corners[1][1][0] - corners[1][0][0] )*x , ( corners[0][1][1] - corners[0][0][1] )*(1.-x) + ( corners[1][1][1] - corners[1][0][1] )*x }
			};

			positions[i*Samples+j] = corners[0][0] * (1.-x) * (1.-y) + corners[0][1] * (1.-x) * (   y) + corners[1][1] * (   x) * (   y) + corners[1][0] * (   x) * (1.-y) ;
			weights[i*Samples+j] = SegmentIntegrator< Samples >::Weights[i] * SegmentIntegrator< Samples >::Weights[j] * ( d[0][0] * d[1][1] - d[0][1] * d[1][0] );
		}
	};
	double weights[Samples*Samples];
	Point2D< double > positions[Samples*Samples];
};

template< unsigned int Samples >
struct TriangleIntegrator
{
	static const double Weights[];
	static const Point2D< double > Positions[];
};
template<> const double TriangleIntegrator<1>::Weights[] =
{
	1.
};
template<> const Point2D< double > TriangleIntegrator<1>::Positions[] =
{
	Point2D< double >( 1./3 , 1./3 )
};
template<> const double TriangleIntegrator<3>::Weights[] =
{
	1./3 ,
	1./3 ,
	1./3
};
template<> const Point2D< double > TriangleIntegrator<3>::Positions[] =
{
	Point2D< double >( 1./6 , 1./6 ) ,
	Point2D< double >( 1./6 , 2./3 ) ,
	Point2D< double >( 2./3 , 1./6 )
};
template<> const double TriangleIntegrator<6>::Weights[] =
{
	0.109951743655322 ,
	0.109951743655322 ,
	0.109951743655322 ,
	0.223381589678011 ,
	0.223381589678011 ,
	0.223381589678011
};
template<> const Point2D< double > TriangleIntegrator<6>::Positions[] =
{
	Point2D< double >( 0.816847572980459  , 0.0915762135097705 ),
	Point2D< double >( 0.0915762135097705 , 0.816847572980459  ),
	Point2D< double >( 0.0915762135097705 , 0.0915762135097705 ),
	Point2D< double >( 0.108103018168070  , 0.445948490915965  ),
	Point2D< double >( 0.445948490915965  , 0.108103018168070  ),
	Point2D< double >( 0.445948490915965  , 0.445948490915965  )
};
template<> const double TriangleIntegrator<32>::Weights[] =
{
	1.1887566790227083e-01 ,
	1.5044412520664885e-01 ,
	1.2632909284531338e-01 ,
	1.0192984975357525e-01 ,
	9.4999150650614317e-02 ,
	4.4981492398316447e-02 ,
	7.9147211585943858e-02 ,
	1.1997941465421234e-01 ,
	1.0670416609764186e-01 ,
	6.1058344824144795e-02 ,
	8.2563774790925248e-02 ,
	9.6297610073814668e-02 ,
	9.1875684331583440e-02 ,
	6.1150555208077911e-02 ,
	4.3370170834023010e-02 ,
	1.0829374522633514e-01 ,
	5.5887468639759713e-02 ,
	1.3351800054734712e-02 ,
	1.5428984747249670e-02 ,
	5.0198346855370224e-02 ,
	5.6291117210426664e-02 ,
	4.1240008239364231e-02 ,
	1.4239502872161450e-02 ,
	1.3691069308687381e-02 ,
	1.9309417484872689e-02 ,
	3.7090960843213061e-02 ,
	3.6967371622461546e-02 ,
	2.1018653471205032e-02 ,
	9.7760996293200769e-03 ,
	5.6339308919459923e-02 ,
	4.9808146403015403e-02 ,
	2.1361687315256585e-02
};
template<> const Point2D< double > TriangleIntegrator<32>::Positions[] =
{
	Point2D< double >( 3.7986021093401956e-01 , 2.1078525939140391e-01 ) ,
	Point2D< double >( 3.0141709320909305e-01 , 4.0978657777002531e-01 ) ,
	Point2D< double >( 5.5802528953120256e-01 , 2.1377743253005960e-01 ) ,
	Point2D< double >( 1.2512299505810387e-01 , 6.1938125736255578e-01 ) ,
	Point2D< double >( 2.1117939909804934e-01 , 2.4498296509349016e-01 ) ,
	Point2D< double >( 8.5431474947580432e-01 , 7.1871496101589105e-02 ) ,
	Point2D< double >( 7.1788185898052326e-01 , 2.0376848107772977e-01 ) ,
	Point2D< double >( 4.6631787462323071e-01 , 4.0896380449124475e-01 ) ,
	Point2D< double >( 2.5015500335339214e-01 , 6.2768261568031403e-01 ) ,
	Point2D< double >( 7.9955384841381316e-02 , 8.2600331401756000e-01 ) ,
	Point2D< double >( 7.1008125956836521e-01 , 6.4413220382260550e-02 ) ,
	Point2D< double >( 4.9732063377796598e-01 , 7.0566724344036824e-02 ) ,
	Point2D< double >( 2.6077068256562896e-01 , 9.5428585810584610e-02 ) ,
	Point2D< double >( 8.9602705800587434e-02 , 1.1638649906727733e-01 ) ,
	Point2D< double >( 2.3088148766115757e-02 , 7.4918973979067949e-01 ) ,
	Point2D< double >( 1.2953296900433620e-01 , 4.2260565743346001e-01 ) ,
	Point2D< double >( 9.3448087604440955e-02 , 2.4345813394879973e-01 ) ,
	Point2D< double >( 9.5526919357006035e-01 , 2.3551733249578710e-02 ) ,
	Point2D< double >( 8.4593539837314391e-01 , 1.5406460162685609e-01 ) ,
	Point2D< double >( 6.1600929617267497e-01 , 3.6118159118967208e-01 ) ,
	Point2D< double >( 3.9316510319604808e-01 , 5.8168921474014745e-01 ) ,
	Point2D< double >( 1.8920633061715936e-01 , 7.8860171922313160e-01 ) ,
	Point2D< double >( 4.3010560106405471e-02 , 9.4547507322097091e-01 ) ,
	Point2D< double >( 8.5815888421533082e-01 , 0.0000000000000000e+00 ) ,
	Point2D< double >( 6.2731531923241179e-01 , 0.0000000000000000e+00 ) ,
	Point2D< double >( 3.6384660446077510e-01 , 1.4566514788346974e-02 ) ,
	Point2D< double >( 1.5557066896897953e-01 , 2.1152223383121949e-02 ) ,
	Point2D< double >( 2.9754117496841759e-02 , 2.7110971356255786e-02 ) ,
	Point2D< double >( 0.0000000000000000e+00 , 9.2734897448394982e-01 ) ,
	Point2D< double >( 2.5716283623693881e-02 , 5.4444667627192522e-01 ) ,
	Point2D< double >( 2.4506286636990005e-02 , 3.3212908394764507e-01 ) ,
	Point2D< double >( 9.2296909059649268e-03 , 1.4604496167217568e-01 )
};

#endif// RIGHT_TRIANGLE_CUADRATURE
