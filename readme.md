[![Double bubbles sans toil and trouble: discrete circulation-preserving vortex sheets for soap films and foams](http://www.cs.columbia.edu/cg/doublebubbles/title.jpg)](http://www.cs.columbia.edu/cg/doublebubbles/)

SoapFilm3D is an open source project for the physical simulation of soap films and bubbles. It is licensed under Clear BSD License for academic and non-commercial use (other licenses may be obtained by contacting the faculty of the Columbia Computer Graphics Group or a Columbia University licensing officer).

It is the test program accompanying Da, Fang, et al. "Double bubbles sans toil and trouble: discrete circulation-preserving vortex sheets for soap films and foams." ACM Transactions on Graphics (TOG) 34.4 (2015): 149. (http://www.cs.columbia.edu/cg/doublebubbles/), including the multimaterial mesh-based surface tracking library LosTopo. 

This program is built by standard procedures using CMAKE (http://www.cmake.org).

Dependencies
--------------------
The following external libraries are required:

Eigen (http://eigen.tuxfamily.org)
OpenGL and GLUT (http://www.opengl.org/resources/libraries/glut/)
LAPACK (http://www.netlib.org/lapack/)
BLAS (http://www.netlib.org/blas/)
libPNG (https://libpng.sourceforge.io/)

On Mac OS X or Linux-based systems, most of the dependencies are either included, or can be easily installed with Homebrew (https://brew.sh) or the APT package handling utility. 

On Windows you may need manually download and install some of them.

Compilation
-----------------
SoapFilm3D has been tested with Clang (under Mac OS X), and GCC 4.8+ (under Linux). Currently the MSVC compiler is not supported.

To compile SoapFilm3D, you'll need CMake on Mac OS X or Linux.

1. make a directory, say, *build*, with *mkdir build*, enter the *build* directory, type *cmake ..*
2. Optionally you can adjust the options with *ccmake .*
3. type *make* to compile the code. For speeding up the compilation process you may use *make -j*.

Run the Demo
--------------------
Running the executables without command line arguments will display usage
information. All the data files for testing are located in the assets folder.



