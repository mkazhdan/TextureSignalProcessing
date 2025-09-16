<center><h2>Gradient Domain Texture Processing (Version 9.55)</h2></center>
<center>
<a href="#LINKS">links</a>
<a href="#EXECUTABLES">executables</a>
<a href="#USAGE">usage</a>
<a href="#LIBRARY">library</a>
<a href="#COMPILATION">compilation</a>
<a href="#CHANGES">changes</a>
<a href="#SUPPORT">support</a>
</center>
<hr>
This software supports gradient-domain signal processing within a texture atlas. Supported applications include:
<UL>
<LI>(localized) texture smoothing and sharpening
<LI>vector-field visualization akin to line-integral convolution
<LI>computation of single-source geodesics
<LI>simulation of reaction-diffusion following the Gray-Scott model
<LI>dilation of texture into the gutter
<LI>masking of gutter/interior/boundary texels
<LI>solving for the smoothest interpolant within a prescribed subset of texels
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
<a href="TSP.x64.zip">Win64</a><br>
<b>Source Code:</b>
<a href="TSP.Source.zip">ZIP</a> <a href="https://github.com/mkazhdan/TextureSignalProcessing">GitHub</a><br>
<B>Data:</B>
<A HREF="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/TSP.Data.zip">ZIP</A><br>
<b>Older Versions:</b>
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version9.50/">V9.50</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version9.10/">V9.10</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version9.05/">V9.05</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version9.00/">V9.00</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version8.00/">V8.00</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version7.00/">V7.00</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version6.06/">V6.06</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version6.00/">V6.00</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version5.01/">V5.01</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version5.00/">V5.00</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.75/">V4.75</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.60/">V4.60</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.50/">V4.50</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.08/">V4.08</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.07/">V4.07</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.06/">V4.06</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.05/">V4.05</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.03/">V4.03</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.02/">V4.02</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.01/">V4.01</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.00/">V4.00</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version3.00/">V3.00</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version2.00/">V2.00</a>,
<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version1.00/">V1.00</a>
</ul>
<hr>
<a name="EXECUTABLES"><b>EXECUTABLES</b></a><br>
<ul>
These applications support reading in textures meshes in one of two formats:
<UL>
<LI><A href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</A> (to support multiple charts, texture is assumed to be encoded at wedges/corners rather than vertices.)
<LI><A HREF="https://www.fileformat.info/format/wavefrontobj/egff.htm">Wavefront OBJ</A>
</UL>
Input textures are assumed to be images in <I>png</I>, <I>jpg</I>, or <I>jpeg</I> format.


<dl>
<details>
<summary>
<font size="+1"><b>TextureFiltering</b></font>:
Supports the (localized) smoothing and sharpening of a texture by solving a screened Poisson equation which gives the signal whose values match the input and whose gradients match the modulated gradients of the input. If no output texture is specified, the executable will launch an interactive viewer that supports local "painting" of gradient modulation values and prescription of a global interpolation weight.<BR>
In the interactive viewer the modulation can be set globally by dragging the slider on the top left.<BR>
The modulation can be set locally by holding the [SHIFT] key down and either dragging with the left mouse button (to smooth) or the right mouse button (to sharpen).
</summary>
<dt><b>--in</b> &lt;<i>input mesh and texture names</i>&gt;</dt>
<dd> These two strings specify the the names of the mesh and the texture image.
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

</dd><dt>[<b>--jitter</B> &lt;<i>random seed</i>&gt;]</dt>
<dd> If specified, this integer value is used to seed the random number generation for jittering. (This is used to avoid singular situations when mesh vertices fall directly on edges in the texture grid. In such a situation, the executable will issue a warning <B>"Zero row at index ..."</B>.)
</dd>

</dd><dt>[<b>--useDirectSolver</B>]</dt>
<dd> If enabled, this flag specifies that a direct solver should be used (instead of the default multigrid solver).
</dd>


</details>
</dl>
</ul>


<ul>
<dl>
<details>
<summary>
<font size="+1"><b>TextureStitching</b></font>:
Supports the stitching together of multiple (possibly overlapping) textures by solving a screened Poisson equation with value constraints defined by the input texture.
The interactive viewer runs in two modes:
<OL>
<LI> A user specifies a single composite texture and (optionally) a mask file indicating when texels in the composite come from the same source.
In this case gradient constraints are obtained by copying gradients from the composite whenever the two texels defining an edge come from the same source, and setting the gradient constraint to zero along edges coming from different sources. If no mask file is provided, a default mask is created by assigning texels the same color if and only if they are covered by the same chart.<BR>
The viewer shows the stitched texture on the left and the composite texture on the right.
<LI> A user specifies multiple partial texture files and corresponding confidence masks.
In this case gradient constraints are obtained by blending gradients from the different inputs, weighted by confidence, and setting gradients to zero in regions where there are no textures with non-zero confidence.
The viewer shows the stitched texture on the left and a partial texture on the right. The user can selectively replace blended value/gradient constraints with the values/gradients from the partial texture by holding the [SHIFT] key down and dragging over the region to be in-painted.
</OL>
</summary>
<dt><b>--in</b> &lt;<i>input mesh and composite texture</i>&gt;</dt>
<dd> These two strings specify the names of the mesh and the texture image.
</dd>

<dt>[<b>--mask</b> &lt;<i>input mask</i>&gt;]</dt>
<dd> This string specifies the name of the mask image.<br>
Black pixels in the mask file should be used to denote regions where the texel value is unkown. (Results may be unpredictable if it is encoded using lossy compression.)
</dd>

<dt>[<b>--out</b> &lt;<i>output texture</i>&gt;]</dt>
<dd> This string is the name of the file to which the stitched texture will be written.</B>
</dd>

<dt>[<b>--outVCycles</b> &lt;<i>output v-cycles</i>&gt;]</dt>
<dd> This integer specifies the number of v-cycles to use if the stitched texture is output to a file and a direct solver is not used.<BR>
The default value for this parameter is 6.
</dd>

<dt>[<b>--interpolation</b> &lt;<i>interpolation weight</i>&gt;]</dt>
<dd> This floating point values gives the interpolation weight.<BR>
The default value for this parameter is 100.
</dd>

<!--
<dt>[<b>--dilateBounaries</b> &lt;<i>dilation radius</i>&gt;]</dt>
<dd> This integer values gives the radius by which the boundaries of the segments should be dilated before stithing is performed.<BR>
The default value for this parameter is -1, indicating no dilation.
</dd>
-->

</dd><dt>[<b>--jitter</B> &lt;<i>random seed</i>&gt;]</dt>
<dd> If specified, this integer value is used to seed the random number generation for jittering. (This is used to avoid singular situations when mesh vertices fall directly on edges in the texture grid. In such a situation, the executable will issue a warning <B>"Zero row at index ..."</B>.)
</dd>

</dd><dt>[<b>--useDirectSolver</B>]</dt>
<dd> If enabled, this flag specifies that a direct solver should be used (instead of the default multigrid solver).
</dd>

</dd><dt>[<b>--multi</B>]</dt>
<dd> If enabled, this flag specifies that the second and third arguments to the <b>--in</b> parameter are to be interpreted as format specifiers for the textures confidence map files.<BR>
<B>Note:</B> If this flat is enabled, the input masks must be specified using the <b>--mask</b> parameter.
</dd>


</details>
</dl>
</ul>


<ul>
<dl>
<details>
<summary>
<font size="+1"><b>LineIntegralConvolution</b></font>:
Supports the <a href="https://en.wikipedia.org/wiki/Line_integral_convolution">line-integral-convolution</A> visualization of a vector-field through <A HREF="https://dl.acm.org/citation.cfm?id=614456">anisotropic diffusion</A> by defining a new metric on the surface that stretches distances along the vector-field values, diffusing a random color texture with respect to the new anisotropic metric, and then sharpening the resulting signal.
If no output texture is specified, the executable will launch an interactive viewer that supports iteratively stepping through the diffusion.<BR>
Hit [SPACE] to start the iterative solver or hit "+" to advance one iteration at a time.
</summary>
<dt><b>--in</b> &lt;<i>input mesh name</i>&gt;</dt>
<dd> This string specifies the name of the mesh.
</dd>

<dt>[<b>--inVF</b> &lt;<i>vector-field file</i>&gt;]</dt>
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

</dd><dt>[<b>--jitter</B> &lt;<i>random seed</i>&gt;]</dt>
<dd> If specified, this integer value is used to seed the random number generation for jittering. (This is used to avoid singular situations when mesh vertices fall directly on edges in the texture grid. In such a situation, the executable will issue a warning <B>"Zero row at index ..."</B>.)
</dd>

</dd><dt>[<b>--minor</B>]</dt>
<dd> If enabled, this flag specifies that the directions of minimal principal curvature should be used to define the vector-field (instead of the default maximal principal curvature directions).
</dd>

</dd><dt>[<b>--useDirectSolver</B>]</dt>
<dd> If enabled, this flag specifies that a direct solver should be used (instead of the default multigrid solver).
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
<dd> This string specifies the the name of the mesh.
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

</dd><dt>[<b>--jitter</B> &lt;<i>random seed</i>&gt;]</dt>
<dd> If specified, this integer value is used to seed the random number generation for jittering. (This is used to avoid singular situations when mesh vertices fall directly on edges in the texture grid. In such a situation, the executable will issue a warning <B>"Zero row at index ..."</B>.)
</dd>

</dd><dt>[<b>--useDirectSolver</B>]</dt>
<dd> If enabled, this flag specifies that a direct solver should be used (instead of the default multigrid solver).
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
<dd> This string specifies the the name of the mesh.
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

</dd><dt>[<b>--jitter</B> &lt;<i>random seed</i>&gt;]</dt>
<dd> If specified, this integer value is used to seed the random number generation for jittering. (This is used to avoid singular situations when mesh vertices fall directly on edges in the texture grid. In such a situation, the executable will issue a warning <B>"Zero row at index ..."</B>.)
</dd>

</dd><dt>[<b>--useDirectSolver</B>]</dt>
<dd> If enabled, this flag specifies that a direct solver should be used (instead of the default multigrid solver).
</dd>

</dd><dt>[<b>--dots</B>]</dt>
<dd> If enabled, this flag specifies that the feed/kill parameters for dot-formation should be used. Otherwise, the feed/kill parameters for stripes are used.
</dd>

</details>
</dl>
</ul>


<ul>
<dl>
<details>
<summary>
<font size="+1"><b>TextureDilation</b></font>:
Dilates the texture by sampling across the chart seam and pulling interpolated texture values into the gutter texels.
</summary>
<dt><b>--in</b> &lt;<i>input mesh and texture names</i>&gt;</dt>
<dd> These two strings specify the the names of the mesh and the texture image.
</dd>

<dt>[<b>--out</b> &lt;<i>output texture mask</i>&gt;]</dt>
<dd> This string is the name of the file to which the dilated texture will be written.</B>
</dd>

<dt>[<b>--radius</b> &lt;<i>dilation radius</i>&gt;]</dt>
<dd> This integer values gives the radius by which the texture should be dilated.<BR>
The default value for this parameter is 0, indicating no dilation.
</dd>

<dt>[<b>--verbose</b>]</dt>
<dd> If this flag is enabled, performance information is printed to <CODE>stdout</CODE>.
</dd>
</dl>
</ul>

</details>
</dl>
</ul>


<ul>
<dl>
<details>
<summary>
<font size="+1"><b>TextureMasking</b></font>:
Identifies the active texels within a texture mask of prescribed resolution, with gutter texels rendered in red and texels supporting a triangle rendered in blue.
</summary>
<dt><b>--in</b> &lt;<i>input mesh name</i>&gt;</dt>
<dd> This string specifies the the name of the mesh.
</dd>

<dt><b>--res</b> &lt;<i>texture width and texture height</i>&gt;</dt>
<dd> These two integers specify the width and height of the mask.
</dd>

<dt>[<b>--out</b> &lt;<i>output texture mask</i>&gt;]</dt>
<dd> This string is the name of the file to which the texture mask will be written.</B>
</dd>

<dt>[<b>--rasterizer</b>]</dt>
<dd> This integer specifies the type of information to be rasterized. Valid values are:
<UL>
<LI><b>0</b>: active -- all texels whose support overlaps the texture atlas
<LI><b>1</b>: boundary -- all texels whose support overlaps the texture chart boundaries
<LI><b>2</b>: id -- colors texels by the triangle covering them
<LI><b>3</b>: node incidence count (unsigned) -- colors texels by the number of triangles sitting over them
<LI><b>4</b>: node incidence count (signed) -- colors texels by the number of positively oriented triangles over them, minus the number of negatively oriented triangles over them
</UL>
The default value for this parameter is <b>0</B> (active)
</dd>

</details>
</dl>
</ul>


<ul>
<dl>
<details>
<summary>
<font size="+1"><b>SeamStitcher</b></font>:
Ensures an as-smooth-as-possible transition, either across chart boundaries or in regions expressly specified by an input mask, by replacing the values of texels with values minimizing the Dirichlet energy (subject to the constraint that the values of other texels remaine fixed).
</summary>
<dt><b>--in</b> &lt;<i>input mesh and texture</i>&gt;</dt>
<dd> These two strings specify the names of the mesh and the texture image.
</dd>

<dt>[<b>--mask</b> &lt;<i>input mask</i>&gt;]</dt>
<dd> This string specifies the name of the mask image used to define texels whose value is to be replaced.<br>
Black pixels in the mask file indicate that the associated texels values should be updated to produce an as-smooth-as-possible image. (Results may be unpredictable if it is encoded using lossy compression.)<BR>
If specified, the mask resolution should match that of the input texture.<br>
If not specified, values of texels at chart boundaries will be updated.
</dd>

<dt>[<b>--out</b> &lt;<i>output texture</i>&gt;]</dt>
<dd> This string is the name of the file to which the smoothed texture will be written.</B>
</dd>

<dt>[<b>--exterior</b>]</dt>
<dd> If this flag is enabled (and a <b>--mask</b> parmeter is not specified), only the values of boundary texels whose centers are outside the charts are updated.
</dd>

<dt>[<b>--verbose</b>]</dt>
<dd> If this flag is enabled, performance information is printed to <CODE>stdout</CODE>.
</dd>


</details>
</dl>
</ul>

<hr>
<a name="USAGE"><b>USAGE EXAMPLES (WITH SAMPLE DATA)</b></a><br>
For testing purposes, a number of <A HREF="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/TSP.Data.zip">textured mapped models</A> are provided (using the <U>.ply</U> extension).
Of these, <I>David</I> and <I>Julius</I> include normal maps (using the <U>.normap</U> extension), <I>Fertility</I> includes the eight harmonic vector-fields (using the <U>.vf</U> extension), and <I>Rooster</I> uses (partial) texture maps as well as a mask image and confidence maps.

<ul>

<dl>
<details>
<summary>
<font size="+1"><b>TextureFiltering</b></font>
</summary>
To run this executable you must specify the input mesh as well as the texture itself:
<blockquote><code>% Bin/*/TextureFiltering --in ../TSP.Data/David/david.ply ../TSP.Data/David/david.normap</code></blockquote>
This opens a viewer allowing the user to prescribe both global gradient modulation weights (through the slider) and local modulation weights (through a paint-brush interface, by depressing the [SHIFT] key and dragging with the left mouse button to smooth and the right mouse button to sharpen).<BR>
You can also bypass the viewer and output a globally sharpened/smoothed texture to a file:
<blockquote><code>% Bin/*/TextureFiltering --in ../TSP.Data/Julius/julius.ply ../TSP.Data/Julius/julius.normap --out julius.smooth.normap --modulation 0 --interpolation 100</code></blockquote>
Here a modulation weight less than 1 indicates that gradients should be dampened (resulting in smoothing) and a small interpolation weight reduces the interpolation penalty, exaggerating the smoothing.
</details>
</dl>

<dl>
<details>
<summary>
<font size="+1"><b>TextureStitching</b></font>
</summary>
This viewer can be run in one of two modes:
<OL>
<LI>
In addition to the input mesh, specify a (single) composite texture and mask.
If adjacent texels share the same mask color, they are assumed to come from the same source, and the gradient between them is preserved.
Otherwise, the gradient is set to zero. Additionally, a mask color of black is reserved to indicate that the texel value is unknown.<BR>
For example, running
<blockquote><code>% Bin/*/TextureFiltering --in Rooster/rooster.ply ../TSP.Data/Rooster/texels.png --mask ../TSP.Data/Rooster/mask.png</code></blockquote>
opens a viewer showing the stitched texture on the left and the composite texture on the right.
<LI>
In addition to the input mesh, specify (multiple) partial textures and associated confidence maps.
The code blends the gradients in regions of overlap, with weights determined by the mask.
Texel and confidence file names are specified using integer format specifiers, with zero-indexing.
Colors are transformed to scalar confidence values by computing the gray-scale value and normalizing to the range [0,1].<br>
For example, running
<blockquote><code>% Bin/*/TextureFiltering --in Rooster/rooster.ply ../TSP.Data/Rooster/texels-%02d.png --mask ../TSP.Data/Rooster/mask-%02d.png --multi</code></blockquote>
opens a viewer showing the stitched texture on the left and the first partial textures on the right.<BR>
Pressing the 't' key toggles forward through the partial textures and pressing 'T' toggles backwards.<BR>
Holding [SHIFT] and clicking on the stitched model replaces the blended gradients under the paint-brush with the gradients from the currently visualized partial-texture.<BR>
</OL>
You can also bypass the viewer and output the stitched texture to a file:
<blockquote><code>% Bin/*/TextureStitching --in Rooster/rooster.ply ../TSP.Data/Rooster/texels-%02d.png --mask ../TSP.Data/Rooster/mask-%02d.png --multi --out stitched.png</code></blockquote>
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
<blockquote><code>% Bin/*/LineIntegralConvolution --in ../TSP.Data/Fertility/fertility.ply --inVF ../TSP.Data/Fertility/harmonic-001.vf --intrinsicVF</code></blockquote>
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
<a name="library"><b>HEADER-ONLY LIBRARY</b></a><br>
<UL>
<DL>
<font size="+1"><b>include/Src/GradientDomain.h</b></font>
In addition to executables, the gradient-domain processing code can be interfaced through either the <CODE>GradientDomain</CODE> class or the <CODE>HierarchicalGradientDomain</CODE> class declared in <CODE>include/Src/GradientDomain.h</CODE>.
<UL>

<DETAILS>
<SUMMARY>
<font size="+1"><CODE>GradientDomain</CODE></font>:
Once constructed, this object can be used to query the active texels/edges, as well as create the standard mass/stiffness/divergence matrices and apply the mass/stiffness/operators to the texel data. In the descriptions below, the template parameter <CODE>Real</CODE> is the floating point type used to represent data (typically <code>double</code>) and <CODE>Solver</CODE> is the class used to factor and solve the sparse system of linear equations (<CODE>Eigen::SimplicialLDLT</CODE> by default).
</SUMMARY>
<B>Code Description</B>:<br>
<UL>
The code performs basic gradient domain processing applications including texture smoothing/sharpening and stitching.
This is done by solving for the output texture values which simultaneouly fit value and derivative constraints.
<UL>
<LI> Value constraints are described by specifying the desired values at the texels, defined to be the input texel values.
<LI> Gradient constraints are described by specifying the desired differences across edges between texels, defined to be the dampened/amplified differences between input texture values (and zerod out if the texels come from different patches, in the case of stitching).
</UL>
</UL>

<B>Code walk-through</B>:<br>
<UL>
  The details of the implementation can be found in the <code>GradientDomain.example.cpp</code> code.
  <UL>
    <LI><U>Lines 113-125</U>: The texture-mapped geometry and the texture image (as well as a mask image describing when texels belong to the same patch, for stitching) are read in.
    <LI><U>Lines 143-157</U>: The <CODE>GradientDomain</CODE> object is constructed, passing in the resolution of the mesh as well as functor giving the indices of embedding/texture-vertices for each corner, functors giving the positions of embedding/texture-vertices, the texture image resolution, and the number of quadrature points per triangle used for integration (valid values are 1, 3, 6, 12, 24, and 32).
    <LI><U>Lines 161-166</U>: The input texture values are read from the image into a <CODE>std::vector</CODE>, using the member functions <CODE>GradientDomain::numNodes</CODE> to get the number of (active) texels in the texture map and <CODE>GradientDomain::node</CODE> to get the coordinates of the texel within the image.
    <LI><U>Lines 168-195</U>: The constraints to the linear system are constructed, specifying the target values and gradients:
    <UL>
      <LI><U>Lines 172-173</U>: The target value constraints are constructed by applying the mass matrix to the input texel values.
      <LI><U>Lines 175-190</U>: The target gradient constraints are obtained by computing the target per-edge differences and then computing the divergence:
      <UL>
        <LI><U>Lines 177-186</U>: The target edge differences are obtained by iterating over the edges, computing the difference between the input texel values at the end-points, and scaling by the gradient modulation value (and zeroing out the difference in the case the end-points are assigned different IDs, in the case of stiching). To this end, the member function <CODE>GradientDomain::numEdges</CODE> gives the number of edges, and the member function <CODE>GradientDomain::edge</CODE> returns the indices of the edge's two end-points.
        <LI><U>Lines 188-189</U>: The target gradient constraints are obtained applying the divergence operator to the computed edge differences.
      </UL>
      <LI><U>Lines 192-193</U>: The target value and gradient weights are combined using the weights specified by the user.
    </UL>
    <LI><U>Lines 197-199</U>: The system matrix is constructed by taking the weighted combination of the mass and stiffness matrices (using the same weights for combining the value and gradient constraints).
    <LI><U>Lines 201-211</U>: The system matrix is factored.
    <LI><U>Lines 213-221</U>: The values for the individual image channels are computed by solving the linear system.
    <LI><U>Lines 223-228</U>: The output texel values are written from the <CODE>std::vector</CODE> back into the texture image.
  </UL>
</UL>

<B>Assumptions</B>:<BR>
<UL>
The code make a number of assumptions about the input geometry:
<UL>
<LI>The code <I>should</I> support non-injective texture mappings.
<LI>For numerical purposes neither surface nor texture triangles should be degenerate.
<LI>The indexing of surface vertices is such that the topology implied by the vertex indexing matches that of the surface.
<LI>The indexing of texture vertices is such that the topology implied by the vertex indexing matches the toplogy of the texture atlas. (i.e. A single surface vertex can be associated with different texture vertices if the associated corners are in different charts.)
</UL>
</UL>
</DETAILS>

<DETAILS>
<SUMMARY>
<font size="+1"><CODE>HierarchicalGradientDomain</CODE></font>:
This class derives from <CODE>GradientDomain</CODE> and additionally supports an iterative multigrid solver for computing the solution to the gradient domain problem without requiring a potentially expensive matrix factorization.
</SUMMARY>
The interface is simlar to that of <CODE>GradientDomain</CODE> and the details of the implementation can be found in <CODE>HierarchicalGradientDomain.example.cpp</CODE>. What follows is a list of the minor changes needed in order to access and use the hierarchical solver:
<UL>
<LI><U>Line 152-167</U>: The interface uses the class <CODE>HierarchicalGradientDomain</CODE>, whose constructor takes an additional parameter describing the number of levels in the multigrid hierarchy.<BR>
<B>Note</B>: Because the class is also responsible for relaxing the gradient-domain system, it is templated off the type of the direct solver used to solve the system at the coarsest resolution (second argument) and the type of data being solved for (third argument). <BR>
[Compare to <U>Lines 143-157</U> of <CODE>GradientDomainExample.cpp</CODE>]
<LI><U>Line 170</U>: The <CODE>HierarchicalGradientDomain</CODE> class maintains its own constraint and solution vectors. So rather than allocating the constraint and solution vectors separately, the interface uses pointers to the data maintained in the <CODE>hgd</code> object.<BR>
[Compare to <U>Line 159</U> of <CODE>GradientDomainExample.cpp</CODE>]
<LI><U>Lines 208-210</U>: The <CODE>HierarchicalGradientDomain</CODE> class maintains its own representation of the system matrix. This is initialized by specifying the mass and stiffness weights for the system.<BR>
[Compare to <U>Lines 197-199</U> of <CODE>GradientDomainExample.cpp</CODE>]
<LI><U>Line 212-214</U>: The solution is relaxed by running the prescribed number of v-cycles with the prescribed number of Gauss-Seidel relaxations per hierarchy level.<BR>
[Compare to <U>Lines 201-221</U> of <CODE>GradientDomainExample.cpp</CODE>]
</UL>
</DETAILS>
</UL>


</DL>
</UL>

<hr>
<details>
<summary>
<a name="COMPILATION"><b>COMPILATION AND EXECUTION</b></a><br>
</summary>
<UL>
<LI>The Windows executables require both the <B>glew</B> and <B>glut</B> dynamically linked libraries to run. These can be found <A HREF="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/TSP.DLLs.zip">here</A> and should be included either in the directory with the executables, or in the directory from which the executables are run.</LI>
<LI>Compiling under Windows requires both the <B>glew</B> and <B>glut</B> libraries. These can be found <A HREF="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/TSP.LIBs.zip">here</A> and should be placed in the output directory for linkage.</LI></LI>
<LI>Compilation requires a linear solver. By default, we use the <CODE>LDLt</CODE> implementation provided by  <A HREF="https://eigen.tuxfamily.org/">Eigen</A>. If you have <A HREF="https://www.intel.com/content/www/us/en/docs/oneapi/programming-guide/2024-1/intel-oneapi-math-kernel-library-onemkl.html">Intel's oneMKL</A>, we encourage you to use Eigen's <CODE>Pardiso</CODE> implementation. To to this you will need to enable the <CODE>USE_EIGEN_PARDISO</CODE> flag in <CODE>include/Src/PreProcessing.h</CODE>
</UL>
</details>

<hr>
<details>
<summary>
<a name="CHANGES"><b>HISTORY OF CHANGES</b></a><br>
</summary>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version2.00/">Version 2.00</a>:
<ul>
<li> Added support for reaction-diffusion based on the Gray-Scott model.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version3.00/">Version 3.00</a>:
<ul>
<li> Added support for texture stitching.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.00/">Version 4.00</a>:
<ul>
<li> Added <CODE>Makefile.no_visual</CODE> to allow building texture filtering/stitching applications without visualizations.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.01/">Version 4.01</a>:
<ul>
<li> Added support for reading <code>.obj</code> files.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.02/">Version 4.02</a>:
<ul>
<li> Added support for mask visualization.
<li> Switched exceptions to warnings.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.03/">Version 4.03</a>:
<ul>
<li> Added support for segment boundary dilation in the <CODE>TextureStitching</CODE> code.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.05/">Version 4.05</a>:
<ul>
<li> Modified the <code>--jitter</code> flag to take a random seed.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.06/">Version 4.06</a>:
<ul>
<li> Added support for visualizing weights when using multi-stitching.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.07/">Version 4.07</a>:
<ul>
<li> Added support for providing a separate low-frequency signal for texture processing.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.08/">Version 4.08</a>:
<ul>
<li> Removing numerical issues in loop construction.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.50/">Version 4.50</a>:
<ul>
<li> Code clean-up
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.60/">Version 4.60</a>:
<ul>
<li> Added <B>--mask</B> for specifying mask(s) to support default cross-chart smoothing.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version4.75/">Version 4.75</a>:
<ul>
<li> Added support for 16-bit png files.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version5.00/">Version 5.00</a>:
<ul>
<li> Separated out OpenMP depeendency.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version5.01/">Version 5.01</a>:
<ul>
<li> Add support for processing polygonal (i.e. not necessarily triangular) faces.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version6.00/">Version 6.00</a>:
<ul>
<li> Add funcionality for dilating the texture map and masking out active texels.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version6.05/">Version 6.05</a>:
<ul>
<li> Fixed off-by-half-pixel issue with texture rasterization.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version6.06/">Version 6.06</a>:
<ul>
<LI> Cleaned up texture dilation
<LI> Made non-template, header-only functions <code>inline</code>.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version7.00/">Version 7.00</a>:
<ul>
<LI> Modified to support non-bijectve texture maps.
<LI> Removed dependence on <CODE>Triangle</CODE> code.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version8.00/">Version 8.00</a>:
<ul>
<LI> Added header-only library for standard geometry-prcoessing interfaces, wrapped in <CODE>include/Src/GradientDomain.h</CODE>
<LI> Added example code showing how to use the libary in <CODE>GradientDomain.example.cpp</CODE>.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version9.00/">Version 9.00</a>:
<ul>
<LI> Added <CODE>SeamStitcher</CODE> executable.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version9.05/">Version 9.05</a>:
<ul>
<LI> Added the option to provide a mask to the <CODE>SeamStitcher</CODE> executable to indicate which texels are to be locked.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version9.10/">Version 9.10</a>:
<ul>
<LI> Added options for <CODE>TextureMasking</CODE>.
<LI> Added more numerical stability.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version9.10/">Version 9.10</a>:
<ul>
<LI> Added header-only library for a multigrid solver supporting standard geometry-prcoessing interfaces, wrapped in <CODE>include/Src/GradientDomain.h</CODE>
<LI> Added example code showing how to use the libary in <CODE>HierarchicalGradientDomain.example.cpp</CODE>.
</ul>

<a href="http://www.cs.jhu.edu/~misha/Code/TextureSignalProcessing/Version9.55/">Version 9.55</a>:
<ul>
<LI> Added support for reading/writing <CODE>.exr</CODE> images.
</ul>

</details>


<hr>
<a name="SUPPORT"><b>SUPPORT</b></a><br>
This work genersouly supported by NSF grant #1422325.

<hr>
<a href="http://www.cs.jhu.edu/~misha">HOME</a>