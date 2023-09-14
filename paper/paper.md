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

The alternative way is to construct approximate model Hamiltonians for the nuclear subsystem which enable an _analytic_ treatment of the nuclear motion, typically by leveraging perturbation theory starting from an exactly solvable solution. This approach is called _lattice dynamics_, and the exact reference solution is given in terms of _phonons_ which are eigensolutions to a Hamiltonian of _harmonic_ form, i.e., _quadratic_ in the nuclear displacements. Based on this starting point, _anharmonic_ contributions are included by extending the harmonic reference Hamiltonian with terms that are cubic, quartic, and possibly higher-order in the displacements. These terms can be included via established perturbative techniques. However, as the complexity of the perturbative expansion grows very quickly with system size and perturbation order, this approach is limited to low-order corrections including the third, sometimes fourth order in practice. It follows that the lattice dynamics approach is _not_ formally exact, and a perfect quantitative agreement with experimental observables is not possible in principle. However, in practice the phonon picture is a good starting point for a wide range of materials, and lattice dynamics can provide _excellent_ qualitative microscopic insight into the physical mechanisms that drive certain phenomena, and often even very good quantitative results, as detailed below. Compared to MD simulations, one can say that the _accuracy_ is potentially limited because more approximations are made when describing the nuclear dynamics, but the _precision_ can be excellent because the statistical error can be tightly controlled in the analytic setting, whereas MD approaches, in particular purely _ab initio_ MD, will typically suffer from much larger statistical noise.

The parameters in the lattice dynamics Hamiltonian, typically called _force constants_, can be obtained in a purely perturbative, temperature-independent way by constructing a Taylor expansion around the energetic minimum position. This idea is more than a century old and traces back to Born and von Karman [CITE]. Alternatively, temperature-dependent _effective_, renormalized model Hamiltonians are used in situations where the quadratic term in a bare Taylor expansion is not positive definite, i.e., the average atomic position does _not_ coincide with a minimum of the potential. The classical example ist the ⁴He problem in which a Taylor expansion does not yield well-defined phonons in dense solid He [CITE deBoer1948], whereas well-defined phonons that describe ⁴He satisfactorily can be obtained from an effective, positive-definite second order Hamiltonian, as was shown by Born and coworkers in the 1950's [CITE Born, Hooton]. At that time, it was most practical to obtain these effective phonons in a self-consistent way, hence denoted _self-consistent phonon theory_. An excellent review is given in Ref. [CITE Klein1972].


# Statement of need
The Temperature Dependent Effective Potentials (TDEP) method is a framework to construct and solve temperature-dependent, effective lattice dynamics model Hamiltonians for a variety of properties, and the TDEP code described here is the respective reference implementation.

As discussed in the introduction, the theoretical foundation of (self-consistent) phonon theory is well-established since decades. More recent developments evolve around finding practical ways of implementing this theory in computer simulations, typically based on _density functional theory_ (DFT). This has led to a variety of approaches that tackle the self-consistent phonon problem for anharmonic and dynamically stabilized systems [CITE SCAILD/qSCAILD, SCHA/SSCHA, SCP/Alamode]. Another development was the _non_ self-consistent construction of effective Hamiltonians by optimizing the force constants to the fully anharmonic dynamics observed during MD simulations [cite Levy1984, Dove1986]. The idea was extended later in the TDEP method to describe phonons in dynamically stabilized systems like Zirconium in the high-temperature bcc phase based on _ab initio_ MD simulations [cite Hellman2011, 2013]. A self-consistent extension to TDEP was later proposed in the form of _stochastic_ TDEP (sTDEP), where instead of MD simulations, thermal samples are created from the model Hamiltonian itself, and the force constants are optimized iteratively until self-consistency [Shulumba2017, Benshalom2022].

While effective phonons capture the effect of anharmonic frequency renormalization, they are still non-interacting quasiparticles with inifite lifetime, or infinitesimal linewidth, respectively. The effect of linewidth broadening due to phonon-phonon interactions can be included by using higher-order force constants up to fourth order [CITE Cowley, Hellman2013, Feng2016 ...]. These can be used to get better approximations to the free energy [CITE Wallace], describe thermal transport [CITE Broido, Romero], and linewidth broadening in spectroscopic experiments [CITE neutrons + Raman]. A formal justification in terms of mode-coupling theory, as well as a detailed comparison between bare perturbation theory with force constants from a Taylor expansion, self-consistent effective, and non-self-consistent effective approaches was recently given by some of the authors in Ref. [CITE Castellano].

To extract force constants from thermal snapshots efficiently, TDEP employs the spacegroup symmetry of a given system to rigorously reduce the free parameters in the model to an irreducible set _before_ fitting the model parameters. For example, this reduces the number of harmonic force constants in a bcc lattice with 4x4x4 supercell (128 atoms) from 147456 to only 11 unknowns after employing symmetry arguments. This can speed up the convergence by several orders of magnitude when comparing to a _post hoc_ symmetrization of the force constants [cite Hellman2013b]. Furthermore, additional lattice dynamics sum rules such as acoustic (translational) and rotational invariances, as well the Huang invariances which ensure the correct number and symmetry of elastic constants in the long-wavelength limit [CITE BornHuang]. While TDEP was probably the first approach to exploit all these constraints in a general way for arbitrary systems, other codes have adopted this practice by now [CITE hiphive, Ponce].

TDEP delivers a clean and fast Fortran implementation with message passing interface (MPI) parallelism for both extracting force constants, as well as analyzing (effective) harmonic phonons and explicit anharmonic properties such as thermal transport and the full phonon spectral function across the Brillouin zone, from which spectroscopic properties can be obtained. Explicitly incorporating dielectric response properties for light scattering experiments such as infrared and Raman was recently proposed [CITE Benshalom]. We highlight some applications and results below.

Additionally, TDEP provides tools to prepare and organize _ab initio_ supercell simulations, e.g., analyzing the crystal symmetry, finding good supercell sizes, visualizing the pair distribution functions from MD simulations, and creating thermal snapshots for accelerated and self-consistent sampling. Each program is fully documented with background information, and an extensive set of realistic research workflow tutorials is available as well. A list of the most important available features and respective programs is given below.

Another distinctive feature of TDEP is the use of plain input and output files which are code agnostic and easy to create or parse. These are either plain text formats, established human-readbale formats like csv, or self-documented HDF5 files for larger datasets. Thanks to exploiting the symmetry of force constants, the respective output files are very compact, even for higher-order force constants.

A separate python library for interfacing with different DFT and force field codes through the atomic simulation environment (ASE) [CITE Ask], as well as processing and further analysis of TDEP output files is available as well [CITE tdeptools].

## Features

Here we list the most important codes that are shipped with the TDEP code and explain their purpose. Are more detailed explanation of all features can be found in the online documentation.

- `extract_forceconstants`: Obtain force constants up to fourths order from a set of snapshots with positions and forces.

- `phonon_dispersion_relations`: Calculate phonon dispersion relations and related harmonic thermodynamic properties from the second-order force constants.

- `thermal_conductivity`: Compute thermal transport by solving the phonon Boltzmann transport equation with perturbative treatment of third-order anharmonicity.
- `lineshape`: Compute phonon spectral functions including lifetime broadening and shifts for single q-points, q-point meshes, or q-point paths in the Brillouin zone.
- `canonical_configuration`: Create supercells with thermal displacements from the force constants via Monte Carlo sampling from a classical and quantum canonical distribution.
- `generate_structure`: Generate supercells of target size that are as cubic as possible to maximize the largest possible real-space cutoff for the force constants.



## Overview of results

- some results?

# Summary

The TDEP method is a versatile and efficient approach to perform temperature-dependent materials simulations from first principles. This comprises thermodynamic properties in classical and quantum ensembles, and several response properties ranging from thermal transport to Neutron and Raman spectroscopy. A stable and fast reference implementation is given in the software package of the same name, which we have described here. The underlying theoretical framework and foundation was briefly sketched with an emphasis on discerning the conceptual difference between bare and effective phonon theory in self-consistent and non-self-consistent formulations. References to a more in-depth discussion of the theory is given in the introduction.

# Acknowledgements
- SeRC
- VR grants

# References
