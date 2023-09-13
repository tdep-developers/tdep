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

Properties of materials change with temperature, i.e., the vibrational motion of electrons and nuclei. In a static thermal equilibrium, temperature determines the structural phase, the density, and many mechanical properties. Out of thermal equilibrium, when applying a thermal gradient or an external spectroscopic probe such as a light or neutron beam, temperature influences the response of the material to the perturbation, for example its ability to conduct heat through the sample, or the lineshape of the spectroscopic signal. Temperature is therefore at the core of materials properties both in applied and fundamental sciences.

In _ab initio_ materials modeling, the electronic temperature contribution is straightforward to include through appropriate occupation of the electronic states, whereas the nuclear contribution needs to be accounted for explicitly: This can be done by performing _molecular dynamics_ (MD) simulations which aims at _numerically_ reproducing the thermal nuclear motion in an atomistic simulation, and obtain temperature-dependent observables in equilibrium through static thermal expectation values, or out of equilibrium from time-dependent correlation functions or by directly observing the thermal relaxation dynamics [cite Tuckerman].

The alternative way is to construct approximate model Hamiltonians for the nuclear subsystem which enable an _analytic_ treatment of the nuclear motion, typically by leveraging perturbation theory starting from an exactly solvable solution. This approach is called _lattice dynamics_, and the exact reference solution is given in terms of _phonons_ which are eigensolutions to a Hamiltonian of _harmonic_ form, i.e., _quadratic_ in the nuclear displacements. Based on this starting point, _anharmonic_ contributions are included by extending the harmonic reference Hamiltonian with terms that are cubic, quartic, and possibly higher-order in the displacements. These terms can be included via established perturbative techniques. However, as the complexity of the perturbative expansion grows very quickly with system size and perturbation order, this approach is limited to low-order corrections including the third, sometimes fourth order in practice. It follows that the lattice dynamics approach is _not_ formally exact, and a perfect quantitative agreement with experimental observables is not possible in principle. However, in practice the phonon picture is a good starting point for a wide range of materials, and lattice dynamics can give _very good_ quantitative results, while providing _excellent_ qualitative microscopic insight into the physical mechanisms that drive certain phenomena as detailed below. Compared to MD simulations, one can say that the _accuracy_ is potentially limited because more approximations are made when describing the nuclear dynamics, but the _precision_ can be excellent because the statistical error can be tightly controlled in the analytic setting, whereas MD approaches, in particular purely _ab initio_ MD, will typically suffer from much larger statistical noise.

The Temperature Dependent Effective Potentials (TDEP) method is a framework to construct and solve these lattice dynamics model Hamiltonians for a variety of properties, and the TDEP code described here is the respective reference implementation. 


# Statement of need
The lattice dynamics model Hamiltonians can be constructed in a purely perturbative, temperature-independent way by constructing a Taylor expansion around the energetic minimum positions, or in a temperature-dependent way by constructing _effective_ model Hamiltonians, which furthermore can be done _self-consistently_. A detailed conceptual comparison is given in Ref. [Castellano2023]. The motivation for TDEP derives from being able to describe systems at conditions where the Taylor expansion approach is invalid because the average atomic positions do _not_ coincide with minima of the potential energy.

- 
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


## Overview of results

- 

# Summary

list

cite examples where TDEP has been used

# Acknowledgements
- SeRC
- VR grants

# References
