Installation
===

`TDEP` is written in modern Fortran, and requires a working Fortran environment. All system-specific settings are saved in a file called `important_settings`. There are some example `important_settings.target` template files saved, we advise to pick the one closest to your target system and adjust the respective fields.

### Prerequisitess

- BLAS and LAPACK need to be installed
- FFTW needs to be installed
- MPI needs to be installed
- HDF5 needs to be installed

### Example

For example on Macbook with `gfortran-12`, you would pick `important_settings.apple_silicon_gfortran12`,

```bash
cp important_settings.apple_silicon_gfortran12 important_settings
```

and adjust in particular the `PATH_TO_HDF5_LIB` and `PATH_TO_HDF5_INC`.

Then run the build script via 

```bash
bash build_things.sh --nthreads_make 4
```

### Check your installation

We advise to run the examples in `tdep/examples` to test your installation after compilation.

## Instructions for specific platforms

- [Platform-agnostic installation using `Anaconda`](#Anaconda)

### Anaconda

One convenient, (mostly) platform-agnostic way to install TDEP is to use [Anaconda](https://anaconda.org/).

#### Prepare environment


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

#### Install

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

## Troubleshooting

### "libary not found for -lm"

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
