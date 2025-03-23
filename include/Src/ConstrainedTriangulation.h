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
#ifndef CONSTRAINED_TRIANGULATION_INCLUDED
#define CONSTRAINED_TRIANGULATION_INCLUDED

#include <Misha/Geometry.h>

#include "Triangle.h"

namespace MishaK
{
	template< typename GeometryReal >
	void TriangulatePolygon( const std::vector< Point2D< GeometryReal > > &vertices , std::vector< SimplexIndex< 2 > > &outputTriangles )
	{
		struct JonathanShewchuk::triangulateio< GeometryReal > in, out;

		/* Define input points. */

		in.numberofpoints = (int)vertices.size();

		in.pointlist = (REAL *)malloc(in.numberofpoints * 2 * sizeof(REAL));
		in.pointmarkerlist = (int *)malloc(in.numberofpoints * sizeof(int));
		for (int i = 0; i < in.numberofpoints; i++) {
			in.pointlist[2 * i] = vertices[i][0];
			in.pointlist[2 * i + 1] = vertices[i][1];
			in.pointmarkerlist[i] = 1; //Check boundary markers documentation
		}

		in.numberofsegments = (int)vertices.size();
		in.segmentlist = (int *)malloc(in.numberofsegments * 2 * sizeof(int));
		in.segmentmarkerlist = (int *)malloc(in.numberofsegments * sizeof(int));
		for (int i = 0; i < in.numberofsegments; i++) {
			in.segmentlist[2 * i] = i + 1;
			in.segmentlist[2 * i + 1] = i < (in.numberofsegments - 1) ? (i + 2) : 1;
			in.segmentmarkerlist[i] = 1;
		}
		in.numberofholes = 0;
		in.numberofregions = 0;
		in.numberofpointattributes = 0;

		out.pointlist = (REAL *)NULL;
		out.trianglelist = (int *)NULL;
		out.segmentlist = (int *)NULL;
		out.pointmarkerlist = (int *)NULL;
		out.triangleattributelist = (REAL *)NULL;
		out.segmentmarkerlist = (int *)NULL;
		/* Refine the triangulation according to the attached */
		/*   triangle area constraints.                       */

		char triSwitches[] = "pQ";
		triangulate( triSwitches , &in , &out , (struct JonathanShewchuk::triangulateio< GeometryReal > *) NULL );

		outputTriangles.resize(out.numberoftriangles);
		for (int i = 0; i < out.numberoftriangles; i++) outputTriangles[i] = SimplexIndex< 2 >(out.trianglelist[3 * i] - 1, out.trianglelist[3 * i + 1] - 1, out.trianglelist[3 * i + 2] - 1);

		free(in.pointlist);
		free(in.pointmarkerlist);
		free(in.segmentlist);
		free(in.segmentmarkerlist);


		free(out.pointlist);
		free(out.trianglelist);
		free(out.segmentlist);
		free(out.pointmarkerlist);
		free(out.triangleattributelist);
		free(out.segmentmarkerlist);
	}
}
#endif //CONSTRAINED_TRIANGULATION_INCLUDED
