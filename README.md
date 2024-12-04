Temperature Dependent Effective Potentials (TDEP)
===

[![DOI](https://joss.theoj.org/papers/10.21105/joss.06150/status.svg)](https://doi.org/10.21105/joss.06150)
![GitHub Release](https://img.shields.io/github/v/release/tdep-developers/tdep)
![GitHub License](https://img.shields.io/github/license/tdep-developers/tdep)

Briefly summarized, the package provides all the tools you need to build accurate model Hamiltonians for finite temperature lattice dynamics from first principles. TDEP includes several programs for different tasks:

- `generate_structure`: Generate supercells of target size, with options to make them as cubic as possible to maximize the real-space cutoff for the force constants.

- `canonical_configuration`: Create supercells with thermal displacements from an initial guess or existing force constants, using Monte Carlo sampling from a classical or quantum canonical distribution.

- `extract_forceconstants`: Obtain (effective) harmonic force constants from a set of supercell snapshots with displaced positions and forces. Optionally fit higher-order force constants.

- `phonon_dispersion_relations`: Calculate phonon dispersion relations and related harmonic thermodynamic properties from the second-order force constants.

- `thermal_conductivity`: Compute thermal transport in the mode-coupling formalism including third- and fourth-order anharmonicity.

- `lineshape`: Compute phonon spectral functions including lifetime broadening and shifts for single q-points, q-point meshes, or q-point paths in the Brillouin zone. The grid mode computes _spectral_ thermal transport properties as well.

- `thermal_conductivity_2023`: Compute thermal transport by solving the phonon Boltzmann transport equation with perturbative treatment of third-order anharmonicity. Legacy implementation, the significantly improved program thermal_conductivity should be used!

More details, examples, and theoretical background can be found in the [online documentation](https://tdep-developers.github.io/tdep/program). See [below](#how-to-cite) which references should be cited for which program.

## Tutorials

You can find a range of tutorials for realistic research workflow using TDEP [in a dedicated repository](https://github.com/tdep-developers/tdep-tutorials).

## Installation

[Please find installation instructions in the TDEP repository.](https://github.com/tdep-developers/tdep/blob/main/INSTALL.md)

## Report bugs and issues

[Please use our github issue tracker to report any problems.](https://github.com/tdep-developers/tdep/issues) **Please make sure to include input/output and log files so that we can reproduce and investigate the question.**

## Contribute

[Please find instructions in the repository.](https://github.com/tdep-developers/tdep/blob/main/CONTRIBUTING.md)

## Other things to look at

* `tdeptools` a package to facilitate working with TDEP and perform additional postprocessing can be found here: [https://github.com/flokno/tools.tdep](https://github.com/flokno/tools.tdep)

## How to cite

This software is distributed under the MIT license. If you use it, please consider citing

[F. Knoop *et al.*, J. Open Source Softw **9**(94), 6150 (2024)](https://joss.theoj.org/papers/10.21105/joss.06150)

and the respective publications for the algorithms that were used:

`canonical_configuration`

- Classical statistics: [D. West and S. K. Estreicher, Phys Rev Lett **96**, 115504 (2006)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.96.115504)
- Quantum statistics: [N. Shulumba, O. Hellman, and A. J. Minnich, Phys. Rev. Lett. **119**, 185901 (2017)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.185901)

`extract_forceconstants`

- Second order: [O. Hellman *et al.*, Phys Rev B **87**, 104111 (2013)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.87.104111)

- Third order: [O. Hellman and I. A. Abrikosov, Phys Rev B **88**, 144301 (2013)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.88.144301)
- Fourth order: [A. H. Romero *et al.*, Phys Rev B **91**, 214310 (2015)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.91.214310)

`thermal_conductivity`

- Method: [D. A. Broido *et al.*, Appl Phys Lett **91**, 231922 (2007)](https://doi.org/10.1063/1.2822891)

- Implementation: [O. Hellman and D.A. Broido, Phys. Rev. B **90**, 134309 (2014)](https://dx.doi.org/10.1103/physrevb.90.134309)
- Fourth-order contributions: [J. Klarbring *et al.*, Phys Rev Lett **125**, 045701 (2020)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.045701)
- Off-diagonal contributions: 
    - [Simoncelli, M. & Marzari, N. & Mauri, F. (2019), Nat. Phys. **15** 803-819](https://doi.org/10.1038/s41567-019-0520-x)
    - [Isaeva, L & Barbalinardo, G. & Donadio, D. & Baroni, S., Nat. Comm. **10** 3853 (2019)](https://doi.org/10.1038/s41467-019-11572-4)
    - [Fiorentino, A. & Baroni, S, Phys. Rev. B, **107**, 054311 (2023)](https://doi.org/10.1103/PhysRevB.107.054311)


`lineshape`

- [A. H. Romero *et al.*, Phys Rev B **91**, 214310 (2015)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.91.214310)
- [N. Shulumba, O. Hellman, and A. J. Minnich, Phys. Rev. Lett. **119**, 185901 (2017)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.185901)
- Grid mode for spectral transport: [Đ. Dangić *et al.*, Npj Comput Mater **7**, 57 (2021)](https://www.nature.com/articles/s41524-021-00523-7)

## Troubleshooting

Some common issues:

### Symmetry errors

TDEP is very strict about crystal symmetries. In phonopy world, the symmetry precision is about `1e-10`. If you see an error like

```
ERROR
exit code 4: symmetry error
```

chances are high that your structure input files are not perfectly symmetric and consistent. **Precise input structures are a prerequisite for using TDEP successfully.**
