testfiles:
	extract_forceconstants -rc2 0 -rc3 3 -rc4 0 --polar
	ln -sf outfile.forceconstant infile.forceconstant
	ln -sf outfile.forceconstant_thirdorder infile.forceconstant_thirdorder
	ln -sf outfile.forceconstant_fourthorder infile.forceconstant_fourthorder
	thermal_conductivity -qg 5 5 5 --fourthorder -qg4ph 3 3 3 --seed 42

clean:
	rm -f outfile.*
	rm -f infile.forceconstant*
