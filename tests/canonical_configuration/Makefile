testfiles:
	canonical_configuration --modes | tail -n 12 | tee canonical_configuration.dat
	rm contcar*

test: testfiles
	pytest
