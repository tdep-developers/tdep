folders="anharmonic_free_energy/
atomic_distribution/
canonical_configuration/
crystal_structure_info/
dump_dynamical_matrices/
convert_XX_to_forceconstant/
extract_forceconstants/
generate_structure/
lineshape/
pack_simulation/
phonon_dispersion_relations/
thermal_conductivity_2023/
thermal_conductivity/
"

echo "RUN TESTS"

source 00-set_path.sh

for folder in ${folders}
do
	echo "RUN ${folder}"
	pushd $folder
	make testfiles
	popd
	echo
done
