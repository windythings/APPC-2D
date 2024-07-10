## `MESH2D: Delaunay-based mesh generation in MATLAB`

`MESH2D` is a `MATLAB` / `OCTAVE`-based unstructured mesh-generator for two-dimensional polygonal geometries, providing a range of relatively simple, yet effective two-dimensional meshing algorithms. `MESH2D` includes variations on the "classical" Delaunay refinement technique, a new "Frontal"-Delaunay refinement scheme, a non-linear mesh optimisation method, and auxiliary mesh and geometry pre- and post-processing facilities. 

Algorithms implemented in `MESH2D` are "provably-good" - ensuring convergence, geometrical and topological correctness, and providing guarantees on algorithm termination and worst-case element quality bounds. Support for user-defined "mesh-spacing" functions and "multi-part" geometry definitions is also provided, allowing `MESH2D` to handle a wide range of complex domain types and user-defined constraints. `MESH2D` typically generates very high-quality output, appropriate for a variety of finite-volume/element type applications.

`MESH2D` is a simplified version of my <a href="https://github.com/dengwirda/jigsaw-matlab/">`JIGSAW`</a> mesh-generation algorithm (a `C++` code). `MESH2D` aims to provide a straightforward `MATLAB` / `OCTAVE` implementation of these Delaunay-based triangulation and mesh optimisation techniques. 

### `Code Structure`

`MESH2D` is a pure `MATLAB` / `OCATVE` package, consisting of a core library + associated utilities:

    MESH2D::
    ├── MAIN-DIR. -- core MESH2D library functions. See REFINE2, SMOOTH2, TRIDEMO, etc.
    ├── aabb-tree -- support for fast spatial indexing, via tree-based data-structures.
    ├── geom-util -- geometry processing, repair, etc.
    ├── hfun-util -- mesh-spacing definitions, limiters, etc.
    ├── hjac-util -- solver for Hamilton-Jacobi eqn's.
    ├── mesh-ball -- circumscribing balls, orthogonal balls etc.
    ├── mesh-cost -- mesh cost/quality functions, etc.
    ├── mesh-file -- mesh i/o via ASCII serialisation.
    ├── mesh-util -- meshing/triangulation utility functions.
    └── poly-test -- fast inclusion test for polygons.

### `License`

This program may be freely redistributed under the condition that the copyright notices (including this entire header) are not removed, and no compensation is received through use of the software.  Private, research, and institutional use is free.  You may distribute modified versions of this code `UNDER THE CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR NOTICE IS GIVEN OF THE MODIFICATIONS`. Distribution of this code as part of a commercial system is permissible `ONLY BY DIRECT ARRANGEMENT WITH THE AUTHOR`. (If you are not directly supplying this code to a customer, and you are instead telling them how they can obtain it for free, then you are not required to make any arrangement with me.) 

`DISCLAIMER`:  Neither I nor the University of Sydney warrant this code in any way whatsoever.  This code is provided "as-is" to be used at your own risk.

### `References`

If you make use of `MESH2D` please include a reference to the following! `MESH2D` is designed to provide a simple and easy-to-understand implementation of Delaunay-based mesh-generation techniques. For a much more advanced, and fully three-dimensional mesh-generation library, see the <a href="https://github.com/dengwirda/jigsaw-matlab/">`JIGSAW`</a> package. `MESH2D` makes use of the <a href="https://github.com/dengwirda/aabb-tree">`AABBTREE`</a> and <a href="https://github.com/dengwirda/find-tria">`FINDTRIA`</a> packages to compute efficient spatial queries and intersection tests. 

`[1]` - Darren Engwirda, <a href="http://hdl.handle.net/2123/13148">Locally-optimal Delaunay-refinement and optimisation-based mesh generation</a>, Ph.D. Thesis, School of Mathematics and Statistics, The University of Sydney, September 2014.

`[2]` - Darren Engwirda, Unstructured mesh methods for the Navier-Stokes equations, Honours Thesis, School of Aerospace, Mechanical and Mechatronic Engineering, The University of Sydney, November 2005.


