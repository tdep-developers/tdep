Installation
===

`TDEP` is written in modern Fortran, and requires a working Fortran environment. All system-specific settings are saved in a file called `important_settings`. There are some example `important_settings.target` template files saved, we advise to pick the one closest to your target system and adjust the respective fields.j

## Requirements

- Fortran compiler. Several f2003/2008 features are extensively used, so the compiler should not be super ancient.
- BLAS and LAPACK need to be installed
- FFTW needs to be installed
- MPI needs to be installed
- HDF5 needs to be installed
- python should be installed
- gnuplot should be installed

**Note: The `build_things.sh` script assumes that TDEP was cloned and lives in a git repository. If you wish obtain TDEP in another you have to adjust the script.**

## Example

For example on Macbook with `gfortran-12`, you would pick `important_settings.apple_silicon_gfortran12`,

```bash
cp important_settings.apple_silicon_gfortran12 important_settings
```

and adjust in particular the `PATH_TO_HDF5_LIB` and `PATH_TO_HDF5_INC`.

Then run the build script via 

```bash
bash build_things.sh --nthreads_make 4
```

## Check your installation

We advise to run the examples in `tdep/examples` to test your installation after compilation.

## HDF5

If you need to install HDF5 manually, it can be [downloaded here](https://www.hdfgroup.org/downloads/hdf5/). Installing it should go like

```bash
./configure FC=FC CC=CC --with-fortran --with-fortran2003 --prefix=XX
make
make install
```

where `FC` and `CC` should point to the same Fortran/C compilers you are using to install TDEP.

# Instructions for specific platforms

- [macOS](#macOS)
- [Supercomputers](#Supercomputers)
- [Platform-agnostic installation using Anaconda](#Anaconda)

## macOS

If you are using [`Homebrew`](https://brew.sh/), you can install all dependencies via `brew install`, e.g., something like

```bash
brew install gcc openmpi gnuplot hdf5
```

Then proceed as [in the example above](#Example). Check out the [`important_settings.osx_gfortran_accelerate`](./important_settings.osx_gfortran_accelerate) file as well for reference.

## Supercomputers

This will depend on the supercomputer you work with, but there should be no big problem per se since the requirements are standard software.

There are two template settings files you can look into:

- [`important_settings.sigma`](./important_settings.sigma) is a template file for a supercomputer with traditional Intel architecture and Intel compilers + MKL math library.
- [`important_settings.dardel`](./important_settings.dardel) is for a Cray supercomputer with AMD CPUs where `gfortran` is used to compile.

## Anaconda

One convenient, (mostly) platform-agnostic way to install TDEP is to use [Anaconda](https://anaconda.org/).

### Prepare environment


Create a `conda` environment with 

```bash
conda create --prefix /path/to/where/you/want/to/create/the/conda_env python=3.10
```

Activate:

```bash
conda activate /path/to/where/you/want/to/create/the/conda_env/
```

Install dependencies

```bash
conda install -c conda-forge gfortran openmpi-mpifort scalapack fftw hdf5 
```

### Install

Copy `important_settings.conda` to `important_settings` and adjust the `PREFIX`, i.e.,

```
...
PREFIX=/path/to/where/you/want/to/create/the/conda_env
...
```

**Note: On some apple-silicon devices you need to point to the `libm` library in `/Library/Developer/CommandLineTools`, see [troubleshooting section below](#libary-not-found-for--lm).**

Run

```bash
./build_things.sh --nthreads_make 4
```

This should be it.

# Troubleshooting

## "libary not found for -lm"

If you see the error

```
...
ld: library not found for -lm
collect2: error: ld returned 1 exit status
```

on an Apple device with silicon chip (M1, M2, etc.), try adding the following to your `important_settings`:

```
FCFLAGS_EXTRA="-L/Library/Developer/CommandLineTools//SDKs/MacOSX13.3.sdk/usr/lib/"
```

where the given path should point to you MacOS xcode installation `lib` path.

## "not a git repository"

If you see the error

```
fatal: not a git repository (or any of the parent directories): .git
```

it means you did not clone TDEP from github. In that case, either clone it or adjust the `build_things.sh` script.

