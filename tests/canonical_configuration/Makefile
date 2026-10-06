testfiles:
	canonical_configuration --modes | tail -n 12 | tee canonical_configuration.dat
	rm contcar*
	# --seed: seed 1 twice and seed 0, 2 configurations each
	for run in 1_a 1_b 0; do \
		canonical_configuration -n 2 --seed $${run%_*} > /dev/null && \
		cat contcar_conf0001 contcar_conf0002 > seed_$${run}.dat; \
		rm -f contcar*; \
	done
	# no --seed: two runs, 1 s apart
	for run in a b; do \
		canonical_configuration -n 2 > /dev/null && \
		cat contcar_conf0001 contcar_conf0002 > noseed_$${run}.dat; \
		rm -f contcar*; \
		sleep 1; \
	done

test: testfiles
	pytest
