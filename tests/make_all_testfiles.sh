folders="anharmonic_free_energy/
atomic_distribution/
canonical_configuration/
crystal_structure_info/
dump_dynamical_matrices/
extract_forceconstants/
generate_structure/
lineshape/
pack_simulation/
phonon_dispersion_relations/
thermal_conductivity/
phasespace_surface/
"

for folder in ${folders}
do
        export PATH="../../bin/:$PATH"
	pushd $folder
	make testfiles
	popd
done
