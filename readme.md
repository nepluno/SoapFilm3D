[![Double bubbles sans toil and trouble: discrete circulation-preserving vortex sheets for soap films and foams](http://www.cs.columbia.edu/cg/doublebubbles/title.jpg)](http://www.cs.columbia.edu/cg/doublebubbles/)

SoapFilm3D is an open source project for the physical simulation of soap films and bubbles. It is cross-platform (Mac OS X, Linux, Windows, and more), and licensed under Clear BSD License for academic and non-commercial use (other licenses may be obtained by contacting the faculty of the Columbia Computer Graphics Group or a Columbia University licensing officer).

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
SoapFilm3D has been tested with Clang (under Mac OS X), GCC 4.8+ (under Linux), and Microsoft Visual Studio (under Windows 10).

To compile SoapFilm3D, you'll need CMake on Mac OS X or Linux, or CMake-GUI (https://cmake.org) on Windows.

On Mac OS X or Linux:
1. make a directory, say, *build*, with *mkdir build*, enter the *build* directory, type *cmake ..*
2. Optionally you can adjust the options with *ccmake .*
3. type *make* to compile the code. For speeding up the compilation process you may use *make -j*.

On Windows:
1. open CMake-GUI, enter the correct directory for source code and build. Then click *Configure*, choose your installed version of the Microsoft Visual Studio.
2. after configuration you may find several libraries not found, check the *Advanced* box and locate those missing libraries manually. Please make sure you have picked the libraries corresponding to the architecture you have selected (say, 32-bit libraries for x86, and 64-bit libraries for x64).
3. click generate after fixing all missing variables to generate your Visual Studio solution.
4. open the solution and compile the code.

Run the Demo
--------------------
Running the executables without command line arguments will display usage
information. All the data files for testing are located in the assets folder.



