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
#pragma once
#include <Src/RightTriangleQuadrature.h>

class TextureNodeInfo {
public:
	TextureNodeInfo() {
		tId = -1;
		baricentricCoords = Point2D<double>(0, 0);
		ci = -1;
		cj = -1;
		chartId = -1;
		isInterior = false;
	}
	TextureNodeInfo(int _tId, Point2D<double> _baricentricCoords, int _ci, int _cj, int _chartId, bool _isInterior) {
		tId = _tId;
		baricentricCoords = _baricentricCoords;
		ci = _ci;
		cj = _cj;
		chartId = _chartId;
		isInterior = _isInterior;
	}
	int tId;
	Point2D<double> baricentricCoords;
	int ci;
	int cj;
	int chartId;
	bool isInterior;
};


class CellIndex
{
protected:
	unsigned int v[4];
public:
	CellIndex(void) { v[0] = v[1] = v[2] = v[3] = 0; }
	CellIndex(unsigned int v0, unsigned int v1, unsigned int v2, unsigned int v3) { v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3; }
	unsigned int &operator[](unsigned int idx) { return v[idx]; }
	unsigned int  operator[](unsigned int idx) const { return v[idx]; }
};


class TriangleElementIndex {
protected:
	unsigned int v[6];
public:
	TriangleElementIndex(void) { v[0] = v[1] = v[2] = v[3] = v[4] = v[5] = 0; }
	TriangleElementIndex(unsigned int v0, unsigned int v1, unsigned int v2, unsigned int v3, unsigned int v4, unsigned int v5) { v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3; v[4] = v4; v[5] = v5; }
	unsigned int &operator[](unsigned int idx) { return v[idx]; }
	unsigned int  operator[](unsigned int idx) const { return v[idx]; }
};


double QuadraticElementValue(int elementIndex, Point2D<double> pos) {
	switch (elementIndex) {
	case 0:
		return 2 * pos[0] * pos[0] + 4 * pos[0] * pos[1] + 2 * pos[1] * pos[1] - 3 * pos[0] - 3 * pos[1] + 1;
		break;
	case 1:
		return 2 * pos[0] * pos[0] - 1 * pos[0];
		break;
	case 2:
		return 2 * pos[1] * pos[1] - 1 * pos[1];
		break;
	case 3:
		return 4 * pos[0] * pos[1];
		break;
	case 4:
		return -4 * pos[0] * pos[1] - 4 * pos[1] * pos[1] + 4 * pos[1];
		break;
	case 5:
		return -4 * pos[0] * pos[0] - 4 * pos[0] * pos[1] + 4 * pos[0];
		break;
	default:
		printf("Element out of bounds! \n");
		return 0;
		break;
	}
}

Point2D<double> QuadraticElementGradient(int elementIndex, Point2D<double> pos) {
	switch (elementIndex) {
	case 0:
		return  Point2D<double>(4 * pos[0] + 4 * pos[1] - 3.0, 4 * pos[0] + 4 * pos[1] - 3.0);
		break;
	case 1:
		return Point2D<double>(4 * pos[0] - 1.0, 0.0);
		break;
	case 2:
		return Point2D<double>(0, 4 * pos[1] - 1.0);
		break;
	case 3:
		return Point2D<double>(4 * pos[1], 4 * pos[0]);
		break;
	case 4:
		return  Point2D<double>(-4 * pos[1], -4 * pos[0] - 8 * pos[1] + 4.0);
		break;
	case 5:
		return Point2D<double>(-8 * pos[0] - 4 * pos[1] + 4, -4 * pos[0]);
		break;
	default:
		printf("Element out of bounds! \n");
		return Point2D<double>(0, 0);
		break;
	}
}
void QuadraticElementValuesAndGradients( Point2D< double > pos , double values[] , Point2D< double > gradients[] )
{
	double xx = pos[0]*pos[0] , xy = pos[0]*pos[1] , yy = pos[1]*pos[1] , x=pos[0] , y=pos[1];
	values[0] =  2*xx + 4*xy + 2*yy - 3*x - 3*y + 1 , gradients[0] = Point2D< double >(  4*x + 4*y - 3 ,  4*x + 4*y - 3 );
	values[1] =  2*xx               - 1*x           , gradients[1] = Point2D< double >(  4*x       - 1 ,              0 );
	values[2] =                2*yy       - 1*y     , gradients[2] = Point2D< double >(              0 ,  4*y       - 1 );
	values[3] =         4*xy                        , gradients[3] = Point2D< double >(        4*y     ,  4*x           );
	values[4] =        -4*xy - 4*yy       + 4*y     , gradients[4] = Point2D< double >(      - 4*y     , -4*x - 8*y + 4 );
	values[5] = -4*xx - 4*xy        + 4*x           , gradients[5] = Point2D< double >( -8*x - 4*y + 4 , -4*x           );
}

Point2D<double> QuadraticGradient(double values[6], Point2D<double> pos) {
	Point2D<double> cumValue(0, 0);
	for (int k = 0; k < 6; k++) cumValue += (QuadraticElementGradient(k, pos)*values[k]);
	return cumValue;
}

double QuadraticValue(double values[6], Point2D<double> pos) {
	double cumValue = 0;
	for (int k = 0; k < 6; k++) cumValue += (QuadraticElementValue(k, pos)*values[k]);
	return cumValue;
}

double LinearElementValue(int elementIndex, Point2D<double> pos) {
	switch (elementIndex) {
	case 0:
		return 1.0 - pos[0] - pos[1];
		break;
	case 1:
		return pos[0];
		break;
	case 2:
		return pos[1];
		break;
	default:
		printf("Element out of bounds! \n");
		return 0;
		break;
	}
}

Point2D<double> LinearElementGradient(int elementIndex) {
	switch (elementIndex) {
	case 0:
		return Point2D<double>(-1, -1);
		break;
	case 1:
		return Point2D<double>(1, 0);
		break;
	case 2:
		return Point2D<double>(0, 1);
		break;
	default:
		printf("Element out of bounds! \n");
		return Point2D<double>(0, 0);
		break;
	}
}

double BilinearValue(double f[4], Point2D<double> pos) {
	return f[0] * (1.0 - pos[0]) * (1.0 - pos[1]) + f[1] * pos[0] * (1.0 - pos[1]) + f[2] * pos[0] * pos[1] + f[3] * (1.0 - pos[0]) * pos[1];
}


Point2D<double> BilinearGradient(double f[4], Point2D<double> pos) {
	return Point2D<double>(pos[1] - 1.0, pos[0] - 1.0)*f[0] + Point2D<double>(1.0 - pos[1], -pos[0])*f[1] + Point2D<double>(pos[1], pos[0])*f[2] + Point2D<double>(-pos[1], 1.0 - pos[0])*f[3];
}

double BilinearElementValue(int elementIndex, Point2D<double> pos) {
	switch (elementIndex) {
	case 0:
		return (1.0 - pos[0]) * (1.0 - pos[1]);
		break;
	case 1:
		return pos[0] * (1.0 - pos[1]);
		break;
	case 2:
		return pos[0] * pos[1];
		break;
	case 3:
		return (1.0 - pos[0]) * pos[1];
		break;
	default:
		printf("Element out of bounds! \n");
		return 0;
		break;
	}
}


Point2D<double> BilinearElementGradient(int elementIndex, Point2D<double> pos) {
	switch (elementIndex) {
	case 0:
		return Point2D<double>(pos[1] - 1.0, pos[0] - 1.0);
		break;
	case 1:
		return Point2D<double>(1.0 - pos[1], -pos[0]);
		break;
	case 2:
		return Point2D<double>(pos[1], pos[0]);
		break;
	case 3:
		return  Point2D<double>(-pos[1], 1.0 - pos[0]);
		break;
	default:
		printf("Element out of bounds! \n");
		return Point2D<double>(0, 0);
		break;
	}
}
void BilinearElementValuesAndGradients( Point2D< double > pos , double values[] , Point2D< double > gradients[] )
{
	double x1 = 1.-pos[0] , y1 = 1.-pos[1] , x2 = pos[0] , y2 = pos[1];
	values[0] = x1*y1 , gradients[0] = Point2D< double >(-y1,-x1);
	values[1] = x2*y1 , gradients[1] = Point2D< double >( y1,-x2);
	values[2] = x2*y2 , gradients[2] = Point2D< double >( y2, x2);
	values[3] = x1*y2 , gradients[3] = Point2D< double >(-y2, x1);
}

