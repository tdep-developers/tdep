Temperature Dependent Effective Potentials (TDEP)
===


Briefly summarized, the package provides all the tools you need to build accurate model Hamiltonians for finite temperature lattice dynamics from first principles. TDEP includes several programs which we briefly introduce here. More details can be found in the [online documentation](https://tdep-developers.github.io/tdep/).

- `generate_structure`: Generate supercells of target size, with options to make them as cubic as possible to maximize the real-space cutoff for the force constants.

- `canonical_configuration`: Create supercells with thermal displacements from an initial guess or existing force constants, using Monte Carlo sampling from a classical or quantum canonical distribution.

- `extract_forceconstants`: Obtain (effective) harmonic force constants from a set of supercell snapshots with displaced positions and forces. Optionally fit higher-order force constants.

- `phonon_dispersion_relations`: Calculate phonon dispersion relations and related harmonic thermodynamic properties from the second-order force constants.

- `thermal_conductivity`: Compute thermal transport by solving the phonon Boltzmann transport equation with perturbative treatment of third-order anharmonicity.

- `lineshape`: Compute phonon spectral functions including lifetime broadening and shifts for single q-points, q-point meshes, or q-point paths in the Brillouin zone. The grid mode computes _spectral_ thermal transport properties as well.

## Manual

[The manual has its own page, including example workflows and theoretical background.](https://tdep-developers.github.io/tdep/)

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

This software is distributed under the MIT license. If you use it, please cite the respective publications:

- `generate_structure`:

  - [O. Hellman, I. A. Abrikosov, and S. I. Simak, Phys Rev B **84**, 180301 (2011)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.84.180301)

- `canonical_configuration`: 

  - Classical statistics: [D. West and S. K. Estreicher, Phys Rev Lett **96**, 115504 (2006)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.96.115504)
  - Quantum statistics: [N. Shulumba, O. Hellman, and A. J. Minnich, Phys. Rev. Lett. **119**, 185901 (2017)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.185901)

- `extract_forceconstants`: 

  - Second order: [O. Hellman *et al.*, Phys Rev B **87**, 104111 (2013)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.87.104111)

  - Third order: [O. Hellman and I. A. Abrikosov, Phys Rev B **88**, 144301 (2013)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.88.144301)

  - Fourth order: [A. H. Romero *et al.*, Phys Rev B **91**, 214310 (2015)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.91.214310)

- `phonon_dispersion_relations`: 

  -  [O. Hellman *et al.*, Phys Rev B **87**, 104111 (2013)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.87.104111)

- `thermal_conductivity`: 

  - Method: [D. A. Broido *et al.*, Appl Phys Lett **91**, 231922 (2007)](https://doi.org/10.1063/1.2822891)
  - Implementation: [A. H. Romero *et al.*, Phys Rev B **91**, 214310 (2015)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.91.214310)

  - Fourth order contributions: [J. Klarbring *et al.*, Phys Rev Lett **125**, 045701 (2020)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.045701)

- `lineshape`: 

  - [A. H. Romero *et al.*, Phys Rev B **91**, 214310 (2015)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.91.214310)
  - [N. Shulumba, O. Hellman, and A. J. Minnich, Phys. Rev. Lett. **119**, 185901 (2017)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.185901)
  - Grid mode for spectral transport: [Đ. Dangić *et al.*, Npj Comput Mater **7**, 57 (2021)](https://www.nature.com/articles/s41524-021-00523-7)
