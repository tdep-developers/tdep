source 00-ln_files.sh

phasespace_surface \
        --verbose \
	--qpoint_grid 11 11 11 \
        --qpoint 0.0 0.0 1E-4 \
	--intensities \
        | tee phasespace_surface.log
        # --povray \
