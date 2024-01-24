# Programs

[`generate_structure`: Generate supercells of target size, with options to make them as cubic as possible to maximize the real-space cutoff for the force constants.](generate_structure.md)

[`canonical_configuration`: Create supercells with thermal displacements from an initial guess or existing force constants, using Monte Carlo sampling from a classical or quantum canonical distribution.](canonical_configuration.md)

[`extract_forceconstants`: Obtain (effective) harmonic force constants from a set of supercell snapshots with displaced positions and forces. Optionally fit higher-order force constants.](extract_forceconstants.md)

[`phonon_dispersion_relations`: Calculate phonon dispersion relations and related harmonic thermodynamic properties from the second-order force constants.](phonon_dispersion_relations.md)

[`thermal_conductivity`: Compute thermal transport by solving the phonon Boltzmann transport equation with perturbative treatment of third-order anharmonicity.](thermal_conductivity.md)

[`lineshape`: Compute phonon spectral functions including lifetime broadening and shifts for single q-points, q-point meshes, or q-point paths in the Brillouin zone. The grid mode computes _spectral_ thermal transport properties as well.](lineshape.md)

[`anharmonic_free_energy`: Compute the free energy including anharmonic contributions.](anharmonic_free_energy.md)

[`atomic_distribution`: Compute pair distribution functions.](atomic_distribution.md)

[`pack_simulation`: Clean, compress and store  simulation data.](pack_simulation.md)

[`samples_from_md`: Pick samples from an MD simulation in a clever way.](samples_from_md.md)

[`crystal_structure_info`: Report which crystal structure and spacegroup TDEP sees.](crystal_structure_info.md)

[`refine_structure`: Clean up input structures with imprecise symmetry.](refine_structure.md)