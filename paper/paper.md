---
title: "TDEP: Temperature Dependent Effective Potenials"
tags:
  - Fortran
  - Physics
  - Phonons
  - Temperature
  - Anharmonicity
  - Neutron spectroscopy
authors:
  - name: Florian Knoop
    orcid: 0000-0002-7132-039X
    affiliation: 1
  - name: Igor Abrikosov
    orcid: 0000-0001-7551-4717
    affiliation: 1
  - name: Sergei Simak 
    orcid: 0000-0002-1320-389X
    affiliation: 1
  - name: Olle Hellman
    orcid: 0000-0002-3453-2975
    affiliation: 2
affiliations:
  - name: Linköping University, Linköping, Sweden
    index: 1
  - name:  Weizmann Institute of Science, Rehovot, Israel
    index: 2
date: August 2023
bibliography: paper.bib
---

# Introduction

Many properties of materials are influenced by temperature, i.e., the vibrational motion of nuclei.

expand


# Statement of need
- TDEP contains the reference implementation for the method of the same name
- Clean and rigorous implementation of physical symmetries and invariances
- Fast implementation in Fortran
- code-agnostic input/output files that are very simple to provide
- spectral response code

# Features

- `crystal_structure_info`: Symmetry analysis (seekpath?)

- `generate_structure`: generate supercells for materials simulations (algorithm to find cubic-as-possible cells of target size)

- `extract_forceconstants`: force constants fit obeying physical invariances for finite differences or effective force constants

- `canonical_configuration`: harmonic Monte Carlo samples generator for classical and quantum distributions of nuclei

- `phonon_dispersion_relations`: Calculate phonon dispersion relations and related harmonic thermodynamic properties

- `thermal_conductivity`: Compute thermal transport by solving the phonon Boltzmann transport equation with perturbative treatment of anharmonicity

- `lineshape`: Compute phonon spectral functions including lifetime broadening and shifts

  

# Summary

list

cite examples where TDEP has been used

# Acknowledgements
- SeRC
- VR grants

# References
