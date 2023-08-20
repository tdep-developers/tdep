Installation
===

`TDEP` is written in modern Fortran, and requires a working Fortran environment. All system-specific settings are saved in a file called `important_settings`. There are some example `important_settings.target` template files saved, we advise to pick the one closest to your target system and adjust the respective fields.j

## Requirements

- Fortran compiler. Several f2003/2008 features are extensively used, so the compiler should not be super ancient.
- BLAS and LAPACK need to be installed
- FFTW needs to be installed
- MPI needs to be installed
- HDF5 needs to be installed **with Fortran support**.
- python should be installed
- gnuplot should be installed

**Note: The `build_things.sh` script assumes that TDEP was cloned and lives in a git repository. If you wish obtain TDEP in another you have to adjust the script.**

**If you have a package manager, `homebrew`, `apt-get`, `yay`, `pacman`, you name it, getting these dependencies should be straightforward.**

For example, you would run

```bash
sudo apt-get install lapack blas hdf5
```

and confirm that libraries and include files add up at a place where you find them, e.g., `/usr/lib` or `/usr/local/lib`. **This will be the paths you specify in the important settings.**

## Example

- Pick e.g. the `important_settings.gfortran` and copy it to `important_settings`:

```bash
cp important_settings.gfortran important_settings
```

- adjust:

  - `FORTRAN_COMPILER`: your fortran compiler

  - `FCFLAGS`: compiler settings to make `gfortran` (or whatever other compiler you are using) work

  - `PATH_TO_BLASLAPACK_LIB`: Path to your BLAS and LAPACK installation. Probably `/usr/local` or similar – when in doubt find it via

    ```bash
    find /usr | grep lapack
    ```

    and see where your LAPACK library is actually located on your computer

  - `PATH_TO_BLASLAPACK_INC`: Same but for the INCLUDE files.

  - `PATH_TO_FFTW_LIB`, `PATH_TO_MPI_LIB`, `PATH_TO_HDF5_LIB` same for FFTW, MPI, HDF5

Then run the build script via 

```bash
bash build_things.sh --nthreads_make 4
```

If everything works find, the last lines of the output should be something like:

```bash
...

Printing bashrc_tdep, append these lines to your .bashrc for stuff to work nicely
MANPATH=$MANPATH:/path/to/tdep/man
PATH=$PATH:/path/to/tdep/bin
export MANPATH
export PATH
alias gnuplot='gnuplot -persist'

Everything should be set up and good to go!
```

i.e., add the respective lines to your `.bashrc` and you are all set up!

**If problems occur, please look at the [Troubleshooting section below](#Troubleshooting). If you cannot fix the error, please reach out, e.g., via the [issue tracker](https://github.com/tdep-developers/tdep/issues).**

## Check your installation

We advise to run the tests in [`./tests`](./tests)  to check your installation.

## HDF5

If you need to install HDF5 manually, it can be [downloaded here](https://www.hdfgroup.org/downloads/hdf5/). Installing it should go like

```bash
./configure FC=FC CC=CC --enable-fortran --prefix=XX
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

## "error: unrecognized command line option '-fallow-argument-mismatch'; did you mean '-Wno-argument-mismatch'?"

If you see the error

```
... error: unrecognized command line option '-fallow-argument-mismatch'; did you mean '-Wno-argument-mismatch'?
```

it likely means that you are using an older Fortran version, and should replace `-fallow-argument-mismatch` with `-Wno-argument-mismatch` in your `important_settings` as suggested.
