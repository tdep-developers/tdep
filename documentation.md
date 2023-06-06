project: TDEP
src_dir: ./src/extract_forceconstants
src_dir: ./src/figure_out_symmetry_things
src_dir: ./src/crystal_structure_info
src_dir: ./src/phonon_dispersion_relations
src_dir: ./src/atomic_distribution
src_dir: ./src/thermal_conductivity
src_dir: ./src/samples_from_md
src_dir: ./src/lineshape
src_dir: ./src/generate_structure
src_dir: ./src/canonical_configuration
src_dir: ./src/convert_phonopy_to_forceconstant
src_dir: ./src/convert_abinit_ddb_to_forceconstant
src_dir: ./src/pack_simulation
src_dir: ./src/refine_structure
src_dir: ./src/generate_kpoints
page_dir: ./src/manual
output_dir: ./doc
project_github: https://github.com/ollehellman/tdep-devel
summary: Temperature dependent effective potential 1.2
author: Olle Hellman
github: https://github.com/ollehellman
predocmark: >
display: public
display: private
source: false
graph: false
media_dir: src/media
favicon: src/media/matematico_favicon.png
author_pic: ../src/media/snappingturtle.jpg
version: 1.2
exclude: options.f90
exclude: io.f90
exclude: celldimensions.f90
exclude: fftwrapper.f90
exclude: interpolatesqe.f90
exclude: type_acf.f90
exclude: type_fcmd.f90
exclude: relax_structure.f90
exclude: test_constraints.f90
exclude: velocitydos.f90
exclude: correlationfunction.f90
exclude: mean_square_displacement.f90
exclude: pair_distribution.f90
exclude: pairmapping.f90
exclude: timedistance_correlation.f90
exclude: vectordist.f90
exclude: mfp.f90
exclude: mfp_cumulative_integrations.f90
exclude: pbe.f90
exclude: pbe_kappa.f90
exclude: pbe_qs.f90
exclude: pbe_scf.f90
exclude: phononevents.f90
exclude: scatteringstrengths.f90
exclude: type_lompi.f90
exclude: phonondamping.f90
exclude: phonondamping_gaussian.f90
exclude: phonondamping_tetrahedron.f90
exclude: scatteringrates.f90
exclude: type_gridsim.f90
exclude: gridenergy.f90
exclude: type_equation_of_state.f90
exclude: epot.f90
exclude: energy_dos.f90
exclude: energy_aux.f90
exclude: energy.f90
exclude: dipole.f90
exclude: dipole_extra_1.f90
exclude: dipole_extra_2.f90
exclude: solvers.f90
exclude: solvers_aux.f90
exclude: solvers_lasso.f90
exclude: forceconstant_solvers.F90

### What is this?

The TDEP package is a collection of tools for finite temperature lattice dynamics. Features include, but are not limited to temperature dependent phonon frequencies, anharmonic free energy and lattice thermal conductivity. The package is released under the MIT license, available on [github](http://github.com/ollehellman/tdep-devel).

### What can you do with this code?

The capabilites are growing constantly. A good place to start would be to read through the papers that have used TDEP in the past, the list can be found [here](page/9_publications.html).

### How to?

If you are new, there is a [getting started guide](page/index.html). There are example workflows, for simple and more complicated things.

### Interfaces

This is not a stand-alone package. TDEP interfaces with atomistic codes, such as DFT codes or classical forcefields. The basic requirement is that must be able to calculate forces and energies from a given atomic configuration. So far, TDEP has successfully been used in combination with

* VASP
* Abinit
* FHI-AIMS
* Quantum Espresso
* LAMMPS

Adding support for additional atomistic packages is straightforward, as long as the package can output energies, positions, and forces.

### Developers

TDEP is developed by Olle Hellman.

<!--
 Here are [all](lists/notelist.html) TODOs that are left.
-->
