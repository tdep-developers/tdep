test:
	extract_forceconstants -rc2 3 -rc3 2
	ln -sf outfile.forceconstant infile.forceconstant
	ln -sf outfile.forceconstant_thirdorder infile.forceconstant_thirdorder
	mpirun phonon_dispersion_relations --dos -p

