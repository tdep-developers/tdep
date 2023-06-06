title: Installation
author: Olle Hellman

This package is written mainly in fortran, with parts in c, c++, and python. The minimum requirements are

*	fortran compiler. Tested with gcc 6+, ifort 16+. Several f2003/2008 features are extensively used, so the compiler needs to be rather new.
*  	MPI
*	fftw
* 	hdf5
*	python

In addition, to make use of all features you might want to install

*	gnuplot, tested with 4+ on OSX with [aquaterm](http://sourceforge.net/projects/aquaterm)
*	this documentation is generated with [FORD](https://github.com/cmacmackin/ford).

@note I do all the development in OSX, and install all libraries/dependencies with [homebrew](http://brew.sh). For reference, these are the things I installed, and referred to in the sample important_settings:
`brew install gcc fftw --with-fortran openmpi gnuplot --with-aquaterm hdf5 --with-fortran --with-fortran2003` If you do the same, installation should be a breeze.

To make matters a little complicated, I recently started using submodules. Support for submodules from compilers is slightly hit or miss, but the following versions seem to work properly as far as I can tell:

* `gfortran 8.2` (should work from version 6, haven't bothered checking older)
* `ifort 18.0.3` (compiler bug in `18.0.1` and `18.0.2`)
* `ifort 17.0.1`,`ifort 17.0.2`
* `ifort 16.0.3` (compiler bug in `16.0.0`)

This is by no means a complete list, but these are the compiler versions I have had access to so far.

### Download

@todo Will fix once I make the repository public.

### Dependencies

#### MPI

Nothing special is required by the MPI installation, except that it needs to support `MPI_IN_PLACE`. I have never had any issues with specific MPI implementations.

#### hdf5

It is likely that the cluster you use already has hdf5 installed. If not, compiling is a breeze (except on Cray systems) usually, I use

```bash
./configure FC=XX CC=XX --with-fortran --with-fortran2003 --prefix=XX
make
make install
```

Note that the fortran compiler needs to be exactly the same as the one used to compile TDEP later on.

### Compiling

The build system is a little unorthodox, but straightforward. You need to create a file called `important_settings`, and put it in the root directory of tdep. Several examples are provided, start with the one that seems closest to what you have. An example follows here, for OSX:

```bash
{!important_settings.osx_gfortran_accelerate!}
```
once that is set up, run

```
./build_things.sh
```

And things should start compiling. It will stop at any error, you probably have to install some library or adjust some path in `important_settings`.

### Setting paths

once `build_things.sh` has finished, it should have printed a `bashrc_tdep` file. The rest of the usage guides assumes you put this in your `.bashrc` or equivalent. It sets the path, manpath and a few minor things, contains the following lines:

```bash
{!bashrc_tdep!}
```

If this is set up, you are good to go.

### Complicated dependencies

#### CGAL

You might have noticed that there is an optional dependency on [CGAL](http://http://www.cgal.org/) in the provided `important_settings`. I wrote a Fortran interface to (a very small subset of) CGAL that can be used to generate meshes, triangulations and other things. The reason it is not enabled by default is that making CGAL and Fortran talk to each other is not straightforward. I can only describe how I got it to work. I don't know if all of these steps are necessary, or if there is some other thing I did that got it to work. This is all for OSX 10.11.4

* get the latest CGAL source, I used 4.11, and read the manual installation.
* when configuring cgal with cmake, specify C and c++ compilers manually to `gcc-8` and `g++-8`, installed via homebrew. (I never got it working with clang, my guess is that linking becomes easier if the fortran, C and c++ compilers are all gcc.)
* configure cgal to produce static libraries (untick the box that says `BUILD_SHARED_LIBS`)
* add `-DBOOST_PARAMETER_MAX_ARITY=12` to the compile line

If this made no sense whatsoever, you probably do not want to try and enable this.
