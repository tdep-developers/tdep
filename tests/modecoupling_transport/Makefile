testfiles:
	extract_forceconstants -rc2 0 -rc3 3 --polar
	ln -sf outfile.forceconstant infile.forceconstant
	ln -sf outfile.forceconstant_thirdorder infile.forceconstant_thirdorder
	modecoupling_transport -qg 5 5 5

clean:
	rm -f outfile.*
	rm -f infile.forceconstant*
