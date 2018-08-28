<center><h2>Gradient Domain Texture Processing (Version 2.00)</h2></center>
<center>
<a href="#LINKS">links</a>
<a href="#EXECUTABLES">executables</a>
<a href="#USAGE">usage</a>
<a href="#COMPILATION">compilation</a>
<a href="#CHANGES">changes</a>
<!--
<a href="#SUPPORT">support</a>
-->
</center>
<hr>
This software supports gradient-domain signal processing within a texture atlas. Supported applications include:
<UL>
<LI>(localized) texture smoothing and sharpening,</LI>
<LI>vector-field visualization akin to line-integral convolution,</LI>
<LI>computation of single-source geodesics, and</LI>
<LI>simulation of reaction-diffusion following the Gray-Scott model.</LI>
</UL>
<hr>
<a name="LINKS"><b>LINKS</b></a><br>
<ul>
<b>Papers:</b>
<a href="http://www.cs.jhu.edu/~misha/MyPapers/SIG18.pdf">[Prada, Kazhdan, Chuang, and Hoppe, 2018]</a>,
<a href="https://en.wikipedia.org/wiki/Line_integral_convolution">[Cabral and Leedom, 1993]</a>,
<a href="https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/">[Crane, Weischedel, and Wardetzky, 2013]</a>
<br>
<b>Executables: </b>
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version1.00/TSP.x64.zip">Win64</a><br>
<b>Source Code:</b>
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version1.00/TSP.Source.zip">ZIP</a> <a href="https://github.com/mkazhdan/TextureSignalProcessing">GitHub</a><br>
<B>Data:</B>
<A HREF="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/TSP.Data.zip">ZIP</A><br>
<b>Older Versions:</b>
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version1.00/">V1</a>
</ul>
<hr>
<a name="EXECUTABLES"><b>EXECUTABLES</b></a><br>
<ul>
<dl>
<details>
<summary>
<font size="+1"><b>TextureFiltering</b></font>:
Supports the (localized) smoothing and sharpening of a texture by solving a screened Poisson equation which gives the signal whose values match the input and whose gradients match the modulated gradients of the input. If no output texture is specified, the executable will launch an interactive viewer that supports local "painting" of gradient modulation values and prescription of a global interpolation weight.<BR>
In the interactive viewer the modulation can be set globally by dragging the slider on the top left.<BR>
The modulation can be set locally by holding the [SHIFT] key down and either dragging with the left mouse button (to sharpen) or the right mouse button (to smooth).
</summary>
<dt><b>--in</b> &lt;<i>input mesh and texture names</i>&gt;</dt>
<dd> These two strings specify the the names of the mesh and the texture image.<br>
The input mesh is assumed to be in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, giving the set of vertices with the x-, y-, and z-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i>, and <i>z</i> the set of polygons encoded by two lists. The first gives the indices of the vertices in the polygon (integers). The second gives the texture coordinates at each polygon corner (pairs of floats).<br>
The input texture is assumed to be an image if the file extension is <I>png</I>, <I>jpg</I>, or <I>jpeg</I>, and a normal map if the extension is <I>normap</I>.
</dd>

<dt>[<b>--out</b> &lt;<i>output texture</i>&gt;]</dt>
<dd> This string is the name of the file to which the processed texture will be written.</B>
</dd>

<dt>[<b>--outVCycles</b> &lt;<i>output v-cycles</i>&gt;]</dt>
<dd> This integer specifies the number of v-cycles to use if the processed texture is output to a file and a direct solver is not used.</B>
The default value for this parameter is 6.
</dd>

<dt>[<b>--interpolation</b> &lt;<i>interpolation weight</i>&gt;]</dt>
<dd> This floating point values gives the interpolation weight.<BR>
The default value for this parameter is 1000.
</dd>

<dt>[<b>--modulation</b> &lt;<i>gradient modulation</i>&gt;]</dt>
<dd> This floating point values gives the (uniform) gradient modulation.<BR>
The default value for this parameter is 1.
</dd>

</dd><dt>[<b>--useDirectSolver</B>]</dt>
<dd> If enabled, this flag specifies that a direct solver should be used (instead of the default multigrid solver).
</dd>

</dd><dt>[<b>--jitter</B>]</dt>
<dd> If enabled, this flag specifies that the texture coordinates should be jittered slightly. (This is used to avoid singular situations when mesh vertices fall directly on edges in the texture grid. In such a situation, the executable will issue a warning <B>"Zero row at index ..."</B>.)
</dd>

</details>
</dl>
</ul>


<ul>
<dl>
<details>
<summary>
<font size="+1"><b>LineIntegralConvolution</b></font>:
Supports the <a href="https://en.wikipedia.org/wiki/Line_integral_convolution">line-integral-convolution</A> visualization of a vector-field by defining a new metric on the surface that stretches distances along the vector-field values, diffusing a random color texture with respect to the new anisotropic metric, and then sharpening the resulting signal.
If no output texture is specified, the executable will launch an interactive viewer that supports iteratively stepping through the diffusion.<BR>
Hit [SPACE] to start the iterative solver or hit "+" to advance one iteration at a time.
</summary>
<dt><b>--in</b> &lt;<i>input mesh name</i>&gt;</dt>
<dd> This string specifies the name of the mesh.<br>
The input mesh is assumed to be in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, giving the set of vertices with the x-, y-, and z-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i>, and <i>z</i> the set of polygons encoded by two lists. The first gives the indices of the vertices in the polygon (integers). The second gives the texture coordinates at each polygon corner (pairs of floats).<br>
</dd>

<dt>[<b>--vf</b> &lt;<i>vector-field file</i>&gt;]</dt>
<DD>This string specifies the file containing the vector-field for visualization. (If this parameter is not specified, the principal curvature direction is used.)<BR>
This file is assumed to be in binary, with the first four bytes storing an integer representing the number of vectors (this should be equal to the number of triangles in the mesh) followed by the list of vectors.
The latter are encoded using double-precision floating point values and should be <I>8</I>*<I>num_triangles</I>*<I>dim</I> bytes, with <I>num_triangles</I> the number of triangles/vectors and <I>dim</I> the dimension of vector-field. (The value of <I>dim</I> is equal to two if the <B>--intrinsicVF</B> is specified an three otherwise.)
</DD>

</dd><dt>[<b>--intrinsicVF</B>]</dt>
<dd> If enabled and a vector-field is specified, this flag indicates that the vector values are represented with two values per vector, using an intrinsic frame. Specifically, for triangle ( <I>v</I><SUB>0</SUB> , <I>v</I><SUB>1</SUB> , <I>v</I><SUB>2</SUB> ), the two-dimensional coefficients ( <I>x</I> , <I>y</I> ) correspond to the three-dimensional tangent vector ( <I>x</I>&middot;(<I>v</I><SUB>1</SUB>-<I>v</I><SUB>0</SUB>) , <I>y</I>&middot;(<I>v</I><SUB>2</SUB>-<I>v</I><SUB>0</SUB>) ).
</dd>

<dt>[<b>--out</b> &lt;<i>output texture</i>&gt;]</dt>
<dd> This string is the name of the file to which the line-integral-convolution texture will be written.</B>
</dd>

<dt>[<b>--outVCycles</b> &lt;<i>output v-cycles</i>&gt;]</dt>
<dd> This integer specifies the number of v-cycles to use if the processed texture is output to a file and a direct solver is not used.</B>
The default value for this parameter is 10.
</dd>

<dt>[<b>--licInterpolation</b> &lt;<i>line-integral-convolution interpolation weight</i>&gt;]</dt>
<dd> This floating point values gives the interpolation weight used for the line-integral-convolution.<BR>
The default value for this parameter is 10000.
</dd>

<dt>[<b>--sharpInterpolation</b> &lt;<i>sharpening interpolation weight</i>&gt;]</dt>
<dd> This floating point values gives the interpolation weight used for sharpening the line-integral-convolution results.<BR>
The default value for this parameter is 10000.
</dd>

<dt>[<b>--modulation</b> &lt;<i>sharpening gradient modulation</i>&gt;]</dt>
<dd> This floating point values gives the gradient modulation used for sharpening the line-integral-convolution results.<BR>
The default value for this parameter is 100.
</dd>

<dt>[<b>--width</b> &lt;<i>output texture width</i>&gt;]</dt>
<dd> This integers specifies the width of the output texture.</B>
The default value for this parameter is 2048.
</dd>

<dt>[<b>--height</b> &lt;<i>output texture height</i>&gt;]</dt>
<dd> This integers specifies the height of the output texture.</B>
The default value for this parameter is 2048.
</dd>

</dd><dt>[<b>--minor</B>]</dt>
<dd> If enabled, this flag specifies that the directions of minimal principal curvature should be used to define the vector-field (instead of the default maximal principal curvature directions).
</dd>

</dd><dt>[<b>--useDirectSolver</B>]</dt>
<dd> If enabled, this flag specifies that a direct solver should be used (instead of the default multigrid solver).
</dd>

</dd><dt>[<b>--jitter</B>]</dt>
<dd> If enabled, this flag specifies that the texture coordinates should be jittered slightly. (This is used to avoid singular situations when mesh vertices fall directly on edges in the texture grid. In such a situation, the executable will issue a warning <B>"Zero row at index ..."</B>.)
</dd>

</details>
</dl>
</ul>


<ul>
<dl>
<details>
<summary>
<font size="+1"><b>Geodesics</b></font>:
An interactive tool for visualization of single-source geodesics using the <A HREF="https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/">heat method</A>.<BR>
In the interactive viewer the source can be set by holding the [SHIFT] key down and clicking/dragging with either mouse button.
</summary>
<dt><b>--in</b> &lt;<i>input mesh name</i>&gt;</dt>
<dd> This string specifies the the name of the mesh.<br>
The input mesh is assumed to be in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, giving the set of vertices with the x-, y-, and z-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i>, and <i>z</i> the set of polygons encoded by two lists. The first gives the indices of the vertices in the polygon (integers). The second gives the texture coordinates at each polygon corner (pairs of floats).<br>
</dd>

<dt>[<b>--interpolation</b> &lt;<i>diffusion interpolation weight</i>&gt;]</dt>
<dd> This floating point values gives the interpolation weight used for diffusing the initial delta function.<BR>
The default value for this parameter is 10000.
</dd>

<dt>[<b>--width</b> &lt;<i>output texture width</i>&gt;]</dt>
<dd> This integers specifies the width of the output texture.</B>
The default value for this parameter is 1024.
</dd>

<dt>[<b>--height</b> &lt;<i>output texture height</i>&gt;]</dt>
<dd> This integers specifies the height of the output texture.</B>
The default value for this parameter is 1024.
</dd>

</dd><dt>[<b>--useDirectSolver</B>]</dt>
<dd> If enabled, this flag specifies that a direct solver should be used (instead of the default multigrid solver).
</dd>

</dd><dt>[<b>--jitter</B>]</dt>
<dd> If enabled, this flag specifies that the texture coordinates should be jittered slightly. (This is used to avoid singular situations when mesh vertices fall directly on edges in the texture grid. In such a situation, the executable will issue a warning <B>"Zero row at index ..."</B>.)
</dd>

</details>
</dl>
</ul>

<ul>
<dl>
<details>
<summary>
<font size="+1"><b>ReactionDiffusion</b></font>:
Performs simulation of reaction-diffusion based on the <A HREF="http://www.karlsims.com/rd.html">Gray-Scott model</A>.
If no output texture is specified, the executable will launch an interactive viewer that supports iteratively stepping through the reaction-diffusion process.<BR>
Hit [SPACE] to start the reaction-diffusion process or hit "+" to advance one step at a time.
</summary>
<dt><b>--in</b> &lt;<i>input mesh name</i>&gt;</dt>
<dd> This string specifies the the name of the mesh.<br>
The input mesh is assumed to be in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, giving the set of vertices with the x-, y-, and z-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i>, and <i>z</i> the set of polygons encoded by two lists. The first gives the indices of the vertices in the polygon (integers). The second gives the texture coordinates at each polygon corner (pairs of floats).<br>
</dd>

<dt>[<b>--out</b> &lt;<i>output texture</i>&gt;]</dt>
<dd> This string is the name of the file to which the reaction-diffusion texture will be written.</B>
</dd>

<dt>[<b>--outSteps</b> &lt;<i>output reaction-diffusion steps</i>&gt;]</dt>
<dd> This integer specifies the number of reaction-diffusion steps to be taken.</B>
The default value for this parameter is 1000.
</dd>

<dt>[<b>--width</b> &lt;<i>output texture width</i>&gt;]</dt>
<dd> This integers specifies the width of the output texture.</B>
The default value for this parameter is 512.
</dd>

<dt>[<b>--height</b> &lt;<i>output texture height</i>&gt;]</dt>
<dd> This integers specifies the height of the output texture.</B>
The default value for this parameter is 512.
</dd>

</dd><dt>[<b>--useDirectSolver</B>]</dt>
<dd> If enabled, this flag specifies that a direct solver should be used (instead of the default multigrid solver).
</dd>

</dd><dt>[<b>--jitter</B>]</dt>
<dd> If enabled, this flag specifies that the texture coordinates should be jittered slightly. (This is used to avoid singular situations when mesh vertices fall directly on edges in the texture grid. In such a situation, the executable will issue a warning <B>"Zero row at index ..."</B>.)
</dd>

</dd><dt>[<b>--dots</B>]</dt>
<dd> If enabled, this flag specifies that the feed/kill parameters for dot-formation should be used. Otherwise, the feed/kill parameters for stripes are used.
</dd>

</details>
</dl>
</ul>

<hr>
<a name="USAGE"><b>USAGE EXAMPLES (WITH SAMPLE DATA)</b></a><br>
For testing purposes, a number of <A HREF="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/TSP.Data.zip">textured mapped models</A> are provided (using the <U>.ply</U> extension).
Of these, <I>David</I> and <I>Julius</I> include normal maps (using the <U>.normap</U> extension) and <I>Fertility</I> includes the eight harmonic vector-fields (using the <U>.vf</U> extension).

<ul>
<dl>
<details>
<summary>
<font size="+1"><b>TextureFiltering</b></font>
</summary>
To run this executable you must specify the input mesh as well as the texture itself:
<blockquote><code>% Bin/*/TextureFiltering --in ../TSP.Data/David/david.ply ../TSP.Data/David/david.normap</code></blockquote>
This opens a viewer allowing the user to prescribe both global gradient modulation weights (through the slider) and local modulation weights (through a paint-brush interface, by depressing the [SHIFT] key and dragging with the left mouse button to sharpen and the right mouse button to smooth).<BR>
You can also bypass the viewer and output a globally sharpened/smoothed texture to a file:
<blockquote><code>% Bin/*/TextureFiltering --in ../TSP.Data/Julius/julius.ply ../TSP.Data/Julius/julius.normap --out julius.smooth.normap --modulation 0 --interpolation 100</code></blockquote>
Here a modulation weight less than 1 indicates that gradients should be dampened (resulting in smoothing) and a small interpolation weight reduces the interpolation penalty, exaggerating the smoothing.
</details>
</dl>

<dl>
<details>
<summary>
<font size="+1"><b>LineIntegralConvolution</b></font>
</summary>
To run this executable you must specify the input mesh:
<blockquote><code>% Bin/*/LineIntegralConvolution --in ../TSP.Data/Fertility/fertility.ply</code></blockquote>
This opens a viewer visualizing a vector-field by performing anisotropic diffusion to simulate line-integral-convolution. (To start the iterative solver, press the [SPACE] key.) By default, the vector-field used is defined by the (maximal) principal curvature directions.<BR>
You can also explicitly prescribe the vector-field:
<blockquote><code>% Bin/*/LineIntegralConvolution --in ../TSP.Data/Fertility/fertility.ply --vf ../TSP.Data/Fertility/harmonic-001.vf --intrinsicVF</code></blockquote>
(The <b>--intrinsicVF</b> flag is required because the vector-field in the file is represented using two intrinsic coordinates per triangle instead of three extrinsic ones.)<BR>
You can also bypass the viewer and output the line-integral-convolution texture to a file:
<blockquote><code>% Bin/*/LineIntegralConvolution --in ../TSP.Data/Hand/hand.ply --minimal --out hand.minimal.jpg</code></blockquote>
Here a visualization of the minimal principal curvature directions is written out as a texture image.
</details>
</dl>

<dl>
<details>
<summary>
<font size="+1"><b>Geodesics</b></font>
</summary>
To run this executable you must specify the input mesh:
<blockquote><code>% Bin/*/Geodesics --in ../TSP.Data/Bunny/bunny.ply</code></blockquote>
This opens a viewer allowing the user to prescribe the source of the geodesic by holding the [SHIFT] button and clicking on the source location with either mouse button.
</details>
</dl>


<dl>
<details>
<summary>
<font size="+1"><b>ReactionDiffusion</b></font>
</summary>
To run this executable you must specify the input mesh:
<blockquote><code>% Bin/*/ReactionDiffusion --in ../TSP.Data/Camel/camel.ply</code></blockquote>
This opens a viewer visualizing the "stripes" reaction-diffusion process. (To start the process, press the [SPACE] key.)<BR>
You can also bypass the viewer and output the reaction-diffusion texture to a file:
<blockquote><code>% Bin/*/ReactionDiffusion --in ../TSP.Data/David/david.ply --out david.dots.jpg --dots --outSteps 2000</code></blockquote>
Here a "dots" pattern is written out to an image. (Empirically, we have found that this reaction-diffusion process takes more steps to converge, hence the larger number of steps.)
</details>
</dl>

</ul>

<hr>
<details>
<summary>
<a name="COMPILATION"><b>COMPILATION AND EXECUTION</b></a><br>
</summary>
<UL>
<LI>The Windows executables require both the <B>glew</B> and <B>glut</B> dynamically linked libraries to run. These can be found <A HREF="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/TSP.DLLs.zip">here</A> and should be included either in the directory with the executables, or in the directory from which the executables are run.</LI>
<LI>Compiling under Windows requires both the <B>glew</B> and <B>glut</B> libraries. These can be found <A HREF="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/TSP.LIBs.zip">here</A> and should be placed in the output directory for linkage.</LI></LI>
</UL>
</details>

<hr>
<details>
<summary>
<a name="CHANGES"><b>HISTORY OF CHANGES</b></a><br>
</summary>
<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version2.00/">Version 2</a>:
<ol>
<li> Added support for reaction-diffusion based on the Gray-Scott model.</li></ol>
</details>


<hr>
<a name="SUPPORT"><b>SUPPORT</b></a><br>
This work genersouly supported by NSF grant #1422325.

<hr>
<a href="http://www.cs.jhu.edu/~misha">HOME</a>
