cases = tetragonal_ab_average

testfiles:
	for case in $(cases); do \
		(cd $$case && average_structure > outfile.log 2>&1; echo $$? > outfile.exitcode); \
	done
