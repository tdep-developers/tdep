folders="anharmonic_free_energy/
atomic_distribution/
canonical_configuration/
crystal_structure_info/
extract_forceconstants/
generate_structure/
lineshape/
pack_simulation/
phonon_dispersion_relations/
refine_structure/
samples_from_md/
thermal_conductivity/
"

for folder in ${folders}
do
	pushd $folder
	make testfiles
	popd
done
