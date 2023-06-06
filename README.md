### Temperature Dependent Effective Potentials (TDEP)


Briefly summarized, the package provides all the tools you need to build accurate model Hamiltonians for finite temperature lattice dynamics from first principles. In addition, there are codes to extract numerous physical properties from these model Hamiltonians, including but not limited to:

* Temperature dependent phonon frequencies
* Free energy including anharmonic contributions
* Phonon spectral functions, linewidths and shifts
* Lattice thermal conductivity


It is implementing the methods published here:

*	[Hellman, O., Abrikosov, I. A., & Simak, S. I. (2011). Lattice dynamics of anharmonic solids from first principles. Physical Review B, 84(18), 180301.](http://doi.org/10.1103/PhysRevB.84.180301)

*	[Hellman, O. & Abrikosov, I. A. (2013). Temperature-dependent effective third-order interatomic force constants from first principles. Physical Review B, 88(14), 144301.](http://doi.org/10.1103/PhysRevB.88.144301)

*	[Hellman, O., Steneteg, P., Abrikosov, I. A., & Simak, S. I. (2013). Temperature dependent effective potential method for accurate free energy calculations of solids. Physical Review B, 87(10), 104111.](http://doi.org/10.1103/PhysRevB.87.104111)

Installation instructions and manual has its own [page](http://ollehellman.github.io), including example workflows and theoretical background. This software is distributed under the MIT license. If you use it, please cite the respective publications.

## Other things to look at

* Tutorials can be found here: https://github.com/tdep-developers/tdep-tutorials
* `tdeptools` a package to facilitate working with TDEP and perform additional postprocessing can be found here: https://github.com/flokno/tools.tdep

## Installation

`TDEP` is written in modern Fortran, and requires a working Fortran environment. All system-specific settings are saved in a file called `important_settings`. There are some example `important_settings.target` template files saved, we advise to pick the one closest to your target system and adjust the respective fields.

### Prerequisites

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

That's it, happy computing!
