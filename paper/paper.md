---
title: "TDEP: Temperature Dependent Effective Potenials"
tags:
  - Fortran
  - Physics
  - Phonons
  - Temperature
  - Anharmonicity
  - Thermal transport
  - Neutron spectroscopy
  - Raman spectroscopy
authors:
  - name: Florian Knoop
    orcid: 0000-0002-7132-039X
    affiliation: 1
  - name: Nina Shulumba
    orcid: 0000-0002-2374-7487
    affiliation: 1
  - name: Aloïs Castellano
    orcid: 0000-0002-8783-490X
    affiliation: 3
  - name: J. P. Alvarinhas Batista
    orcid: 0000-0002-3314-249X
    affiliation: 3
  - name: Roberta Farris
    orcid: 0000-0001-6710-0100
    affiliation: 7
  - name: Matthieu J. Verstraete
    orcid: 0000-0001-6921-5163
    affiliation: 3, 8
  - name: Matthew Heine
    orcid: 0000-0002-4882-6712
    affiliation: 5
  - name: David Broido
    orcid: 0000-0003-0182-4450
    affiliation: 5
  - name: Dennis S. Kim
    orcid: 0000-0002-5707-2609
    affiliation: 6
  - name: Johan Klarbring
    orcid: 0000-0002-6223-5812
    affiliation: 1, 4
  - name: Igor Abrikosov
    orcid: 0000-0001-7551-4717
    affiliation: 1
  - name: Sergei I. Simak
    orcid: 0000-0002-1320-389X
    affiliation: 1, 9
  - name: Olle Hellman
    orcid: 0000-0002-3453-2975
    affiliation: 2
affiliations:
  - name: Theoretical Physics Division, Department of Physics, Chemistry and Biology (IFM), Linköping University, SE-581 83 Linköping, Sweden
    index: 1
  - name:  Weizmann Institute of Science, Rehovot, Israel
    index: 2
  - name:  Nanomat group, QMAT center, CESAM research unit and European Theoretical Spectroscopy Facility, Université de Liège, allée du 6 août, 19, B-4000 Liège, Belgium
    index: 3
  - name: Department of Materials, Imperial College London, South Kensington Campus, London SW7 2AZ, UK
    index: 4
  - name:  Department of Physics, Boston College, Chestnut Hill, MA 02467, USA
    index: 5
  - name: College of Letters and Science, Department of Chemistry and Biochemistry, University of California, Los Angeles (UCLA), California 90025, USA
    index: 6
  - name: Catalan Institute of Nanoscience and Nanotechnology - ICN2 (BIST and CSIC), Campus UAB, 08193 Bellaterra (Barcelona), Spain
    index: 7
  - name: ITP, Physics Department, University of Utrecht, 3584 CC Utrecht, the Netherlands
    index: 8
  - name: Department of Physics and Astronomy, Uppsala University, SE-75120 Uppsala, Sweden
    index: 9
date: 2023
bibliography: literature.bib
---

# Summary

The Temperature Dependent Effective Potential (TDEP) method is a versatile and efficient approach to include temperature in _ab initio_ materials simulations based on phonon theory. TDEP can be used to describe thermodynamic properties in classical and quantum ensembles, and several response properties ranging from thermal transport to Neutron and Raman spectroscopy. A stable and fast reference implementation is given in the software package of the same name described here. The underlying theoretical framework and foundation is briefly sketched with an emphasis on discerning the conceptual difference between bare and effective phonon theory, in both self-consistent and non-self-consistent formulations. References to numerous applications and more in-depth discussions of the theory are given.

# Introduction

The properties of materials change both qualitatively and quantitatively with temperature, i.e., the macroscopic manifestation of the microscopic vibrational motion of electrons and nuclei. In thermal equilibrium, temperature influences the structural phase, the density, and many mechanical properties. Out of thermal equilibrium, for instance when applying a thermal gradient or an external spectroscopic probe such as a light or neutron beam, temperature influences the response of the material to the perturbation, for example its ability to conduct heat, or the lineshape of the spectroscopic signal. Temperature is therefore at the core of both applied and fundamental materials science.

In _ab initio_ materials modeling, the contribution of _electronic_ temperature is straightforward to include through appropriate occupation of the electronic states, whereas the nuclear contribution needs to be accounted for explicitly. This can be done by performing _molecular dynamics_ (MD) simulations which aim at _numerically_ reproducing the thermal nuclear motion in an atomistic simulation, and obtaining temperature-dependent observables either in equilibrium, through averaging, or out of equilibrium from time-dependent correlation functions or by directly observing the relaxation dynamics.

An alternative strategy is to construct approximate model Hamiltonians for the nuclear subsystem, which enables an _analytic_ description of the nuclear motion by leveraging perturbation theory starting from an exactly solvable lowest order model. In this approach, the starting point is given in terms of _phonons_ which are eigensolutions of a Hamiltonian of _harmonic_ form, i.e., _quadratic_ in the nuclear displacements. _Anharmonic_ contributions can be included via established perturbative techniques, in practice up to quartic terms. Higher-order contributions are elusive because the complexity and number of terms in the perturbative expansion grows very quickly with system size and perturbation order. The lattice dynamics approach is therefore not formally exact. However, the phonon picture is useful for describing a wide range of materials properties in practice, and often reaches excellent accuracy in comparison to experiment while providing precise microscopic insight into the underlying physical phenomena.

The chemical bonding in the lattice dynamics Hamiltonian is represented through _force constants_. These can be obtained in a purely perturbative, temperature-independent way by constructing a Taylor expansion of the interatomic potential energy about the periodically arranged atom positions in the crystal. This idea is more than a century old and traces back to Born and von Karman [@Born.1912]. Alternatively, temperature-dependent, _effective_ model Hamiltonians are used in situations where the quadratic term in a bare Taylor expansion is not positive definite, i.e., the average atomic position does _not_ coincide with a minimum of the potential. The classic example is the $^4$He problem in which a Taylor expansion led to imaginary phonon frequencies in the dense solid phase [@Boer.1948]. Born and coworkers solved this problem in the 1950's by developing a _self-consistent phonon theory_ in which an effective, positive-definite Hamiltonian yielding well-defined phonons is obtained self-consistently using a variational principle [@Born.1951; @Hooton.2010mfn; @Hooton.2010]. An excellent historical review of this development is given in Ref. [@Klein.1972].

While the theoretical foundation of (self-consistent) phonon theory has been well-established for decades, more recent developments are concerned with implementing this theory in computer simulations, typically based on _density functional theory_ (DFT) [@Hohenberg.1964; @Kohn.1965]. This has led to a variety of approaches that tackle the self-consistent phonon problem for anharmonic and dynamically stabilized systems [@Souvatzis.2008; @Errea.2013; @Tadano.2015; @Roekeghem.2021; @Monacelli.2021]. Another development was the _non_ self-consistent construction of effective Hamiltonians by optimizing the force constants to the fully anharmonic dynamics observed during MD simulations [@Levy.1984; @Dove.1986]. The idea was extended in the TDEP method to describe phonons in dynamically stabilized systems like Zirconium in the high-temperature bcc phase based on _ab initio_ MD simulations [@Hellman.2011; @Hellman.2013]. A self-consistent extension to TDEP was later proposed in the form of _stochastic_ TDEP (sTDEP), where thermal samples are created from the model Hamiltonian itself instead of using MD, and the force constants are optimized iteratively until self-consistency [@Shulumba.2017; @Benshalom.2022]. sTDEP furthermore allows to include nuclear quantum effects in materials with light elements in a straightforward way [@Shulumba.20179s8e; @Laniel.2022].

Effective phonons capture anharmonic frequency renormalization, but they are still non-interacting quasiparticles with infinite lifetime, or equivalently infinitesimal linewidth. The effect of linewidth broadening due to anharmonic phonon-phonon interactions can be included by using higher-order force constants up to third or fourth order [@Cowley.1963; @Hellman.2013oi5; @Feng.2016]. These can be used to get better approximations to the free energy [@Wallace.1972], describe thermal transport [@Broido.2007; @Romero.2015; @Dangić.2021], and linewidth broadening in spectroscopic experiments [@Romero.2015; @Kim.2018; @Benshalom.2022]. In practice, it was noted that the renormalized phonon quasiparticles interact more weakly than bare phonons. This means that the effective approach remains applicable in systems with strong anharmonicity where the bare phonon quasiparticle picture becomes invalid [@Ravichandran.2018]. A formal justification in terms of mode-coupling theory, as well as a detailed comparison between bare perturbation theory with force constants from a Taylor expansion, self-consistent effective, and non-self-consistent effective approaches was recently given by some of us in Ref. [@Castellano.2023]. Explicitly incorporating dielectric response properties for light scattering experiments such as infrared and Raman was recently proposed [@Benshalom.2022].


# Statement of need
The TDEP open source code is the reference implementation for the TDEP method introduced above. It delivers a clean and fast Fortran implementation with message passing interface (MPI) parallelism both for constructing and solving effective lattice-dynamics Hamiltonians. This allows for materials simulations of simple elemental solids up to complex compounds with reduced symmetry under realistic conditions.

To extract force constants from thermal snapshots efficiently, TDEP employs the permutation and spacegroup symmetries of a given system to reduce the free parameters in the model to an irreducible set, before fitting them [@Esfarjani.2008]. For example, this reduces the number of harmonic force constants of a 4x4x4 supercell of a bcc lattice (128 atoms) from 147456 to only 11 unknowns. This can speed up the convergence by several orders of magnitude when comparing to a _post hoc_ symmetrization of the force constants [@Hellman.2013oi5]. Further lattice dynamics sum rules are enforced after fitting, i.e., acoustic (translational) and rotational invariances, as well as the Huang invariances, which ensure the correct number and symmetry of elastic constants in the long-wavelength limit [@Born.1954]. While TDEP was one of the first numerical approaches exploiting all these constraints in a general way for arbitrary systems, other codes have adopted this practice by now [@Eriksson.2019; @Lin.2022].

Another distinctive feature of TDEP is the use of plain input and output files which are code agnostic and easy to create and parse. These are either plain text formats, established human-readable formats like csv, or self-documented HDF5 files for larger datasets. Thanks to the exploitation of the force constant symmetries, the respective output files are very compact, even for anharmonic force constants.

Additionally, TDEP provides tools to prepare and organize _ab initio_ supercell simulations, e.g., analyzing the crystal symmetry, finding good simulation (super)cells, visualizing the pair distribution functions from MD simulations, and creating thermal snapshots for accelerated and self-consistent sampling. Each program is fully documented with background information, and an extensive set of realistic research workflow tutorials is available as well. A list of the most important available features and respective programs is given below.

## Features

Here we list the most important codes that are shipped with the TDEP code, explain their purpose, and list the respective references in the literature. A more detailed explanation of all features can be found in the online documentation.

- `generate_structure`: Generate supercells of target size, with options to make them as cubic as possible to maximize the real-space cutoff for the force constants [@Hellman.2011].

- `canonical_configuration`: Create supercells with thermal displacements from an initial guess or existing force constants, using Monte Carlo sampling from a classical or quantum canonical distribution [@West.2006; @Shulumba.2017]. Self-consistent sampling with sTDEP is explained in detail in [@Benshalom.2022].

- `extract_forceconstants`: Obtain (effective) harmonic force constants from a set of supercell snapshots with displaced positions and forces [@Hellman.2013]. Optionally fit higher-order force constants [@Hellman.2013oi5], or dielectric tensor properties [@Benshalom.2022].

- `phonon_dispersion_relations`: Calculate phonon dispersion relations and related harmonic thermodynamic properties from the second-order force constants [@Hellman.2013], including Grüneisen parameters from third-order force constants [@Hellman.2013oi5].

- `thermal_conductivity`: Compute thermal transport by solving the phonon Boltzmann transport equation with perturbative treatment of third-order anharmonicity [@Broido.2007; @Romero.2015].

- `lineshape`: Compute phonon spectral functions including lifetime broadening and shifts for single q-points, q-point meshes, or q-point paths in the Brillouin zone [@Romero.2015; @Shulumba.2017]. The grid mode computes _spectral_ thermal transport properties [@Dangić.2021].

A separate python library for interfacing with different DFT and force field codes through the atomic simulation environment (ASE) [@Larsen.2017], as well as processing and further analysis of TDEP output files is available as well [@tdeptools].

We note that parts of the TDEP method have been implemented in other code packages as well [@Bottin.2020rn5].

# Acknowledgements

F.K. acknowledges support from the Swedish Research Council (VR) program 2020-04630, and the Swedish e-Science Research Centre (SeRC). Work at Boston College was supported by the U.S. Department of Energy (DOE), Office of Science, Basic Energy Sciences (BES) under Award #DE-SC0021071. MJV, AC, and JPB acknowledge funding by ARC project DREAMS (G.A. 21/25-11) funded by Federation Wallonie Bruxelles and ULiege, and by Excellence of Science project CONNECT number 40007563 funded by FWO and FNRS. R.F. acknowledges financial support by the Spanish State Research Agency under grant number PID2022-139776NB-C62 funded by MCIN/AEI/ 10.13039/501100011033 and by ERDF A way of making Europe. J. K. acknowledges support from the Swedish Research Council (VR) program 2021-00486.

# References
