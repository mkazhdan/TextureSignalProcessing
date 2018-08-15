<center><h2>Gradient Domain Texture Processing (Version 1.00)</h2></center>
<center>
<a href="#LINKS">links</a>
<a href="#EXECUTABLES">executables</a>
<!--
<a href="#USAGE">usage</a>
-->
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
<LI>vector-field visualization akin to line-integral convolution, and</LI>
<LI>computation of single-source geodesics.</LI>
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
<A HREF="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/TSP.Data.zip">ZIP</A>
<!--
<b>Older Versions:</b>
-->
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
The input mesh is assumed to be in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, giving the set of vertices with the x-, y-, and z-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i>, and the set of polygons encoded by two lists. The first gives the indices of the vertices in the polygon (integers). The second gives the texture coordinates at each polygon corner (pairs of floats).<br>
The input texture is assumed to be an image if the file extension is <I>png</I>, <I>jpg</I>, or <I>jpeg</I>, and a normal map if the extension is <I>normap</I>.
</dd>

<dt>[<b>--out</b> &lt;<i>output texture</i>&gt;]</dt>
<dd> This string is the name of the file to which the processed texture will be written.</B>
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
Creates a <a href="https://en.wikipedia.org/wiki/Line_integral_convolution">line integral convolution</A> visualization of a vector-field by defining a new metric on the surface that stretches distances along the vector-field values, diffuses a random color texture with respect to the new anisotropic metric, and then sharpens the resulting signal.
</summary>
<dt><b>--in</b> &lt;<i>input mesh name</i>&gt;</dt>
<dd> This string specifies the name of the mesh.<br>
The input mesh is assumed to be in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, giving the set of vertices with the x-, y-, and z-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i>, and the set of polygons encoded by two lists. The first gives the indices of the vertices in the polygon (integers). The second gives the texture coordinates at each polygon corner (pairs of floats).<br>
</dd>

<dt>[<b>--vf</b> &lt;<i>vector field file</i>&gt;]</dt>
<DD>This string specifies the file containing the vector field for visualization. (If this parameter is not specified, the principal curvature direction is used.)<BR>
This file is assumed to be in binary, with the first four bytes storing an integer representing the number of vectors (this should be equal to the number of triangles in the mesh) followed by the list of vectors.
The latter are encoded using double-precision floating point values and should be <I>8</I>*<I>num_triangles</I>*<I>dim</I> bytes, with <I>num_triangles</I> the number of triangles/vectors and <I>dim</I> the dimension of vector field. (The value of <I>dim</I> is equal to two if the <B>--intrinsicVF</B> is specified an three otherwise.)
</DD>

<dt>[<b>--out</b> &lt;<i>output texture</i>&gt;]</dt>
<dd> This string is the name of the file to which the processed texture will be written.</B>
</dd>

<dt>[<b>--licInterpolation</b> &lt;<i>line integral convoluation interpolation weight</i>&gt;]</dt>
<dd> This floating point values gives the interpolation weight used for the line integral convolution.<BR>
The default value for this parameter is 10000.
</dd>

<dt>[<b>--sharpInterpolation</b> &lt;<i>sharpening interpolation weight</i>&gt;]</dt>
<dd> This floating point values gives the interpolation weight used for sharpening the line integral convolution results.<BR>
The default value for this parameter is 10000.
</dd>

<dt>[<b>--modulation</b> &lt;<i>sharpening gradient modulation</i>&gt;]</dt>
<dd> This floating point values gives the gradient modulation used for sharpening the line integral convolution results.<BR>
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

</dd><dt>[<b>--intrinsicVF</B>]</dt>
<dd> If enabled and a vector field is specified, this flag indicates that the vector values are represented with two values per vector field, using an intrinsic frame. Specifically, for triangle ( <I>v</I><SUB>0</SUB> , <I>v</I><SUB>1</SUB> , <I>v</I><SUB>2</SUB> ), the two-dimensional coefficients ( <I>x</I> , <I>y</I> ) correspond to the three-dimensional tangent vector ( <I>x</I>&middot;(<I>v</I><SUB>1</SUB>-<I>v</I><SUB>0</SUB>) , <I>y</I>&middot;(<I>v</I><SUB>2</SUB>-<I>v</I><SUB>0</SUB>) ).
</dd>

</dd><dt>[<b>--minor</B>]</dt>
<dd> If enabled, this flag specifies that the directions of minimal principal curvature should be used to define the vector field (instead of the default maximal principal curvature directions).
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
An interactive tool for visualizatin single-source geodesics using the <A HREF="https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/">heat method</A>.<BR>
In the interactive viewer the source can be set by holding the [SHIFT] key down and clicking/dragging with either mouse button.
</summary>
<dt><b>--in</b> &lt;<i>input mesh name</i>&gt;</dt>
<dd> This string specifies the the name of the mesh.<br>
The input mesh is assumed to be in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format, giving the set of vertices with the x-, y-, and z-coordinates of the positions encoded by the properties <i>x</i>, <i>y</i>, and the set of polygons encoded by two lists. The first gives the indices of the vertices in the polygon (integers). The second gives the texture coordinates at each polygon corner (pairs of floats).<br>
</dd>

<dt>[<b>--out</b> &lt;<i>output texture</i>&gt;]</dt>
<dd> This string is the name of the file to which the processed texture will be written.</B>
</dd>

<dt>[<b>--licInterpolation</b> &lt;<i>line integral convoluation interpolation weight</i>&gt;]</dt>
<dd> This floating point values gives the interpolation weight used for the line integral convolution.<BR>
The default value for this parameter is 10000.
</dd>

<dt>[<b>--sharpInterpolation</b> &lt;<i>sharpening interpolation weight</i>&gt;]</dt>
<dd> This floating point values gives the interpolation weight used for sharpening the line integral convolution results.<BR>
The default value for this parameter is 10000.
</dd>

<dt>[<b>--modulation</b> &lt;<i>sharpening gradient modulation</i>&gt;]</dt>
<dd> This floating point values gives the gradient modulation used for sharpening the line integral convolution results.<BR>
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
<dd> If enabled, this flag specifies that the directions of minimal principal curvature should be used to define the vector field (instead of the default maximal principal curvature directions).
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

<!--
<hr>
<a name="USAGE"><b>USAGE EXAMPLES (WITH SAMPLE DATA)</b></a><br>

<ul>
<dl>
<details>
<summary>
<font size="+1"><b>TextureFiltering</b></font>
</summary>
For testing purposes, a number of <A HREF="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Data.zip">datasets</A> are provided.
<ol>
</details>
</dl>
</ul>
-->

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
<!--
<a href="http://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version3/">Version 3</a>:
<ol>
<li> The implementation of the <b>--samplesPerNode</b> parameter has been modified so that a value of "1" more closely corresponds to a distribution with one sample per leaf node.
</li><li> The code has been modified to support compilation under MSVC 2010 and the associated solution and project files are now provided. (Due to a bug in the Visual Studios compiler, this required modifying the implementation of some of the bit-shifting operators.)
</li></ol>
-->
</details>


<hr>
<a name="SUPPORT"><b>SUPPORT</b></a><br>
This work genersouly supported by NSF grant #1422325.

<hr>
<a href="http://www.cs.jhu.edu/~misha">HOME</a>
