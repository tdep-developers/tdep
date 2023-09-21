CMDbib = biber --tool --output_align --output_indent=2 --output_fieldcase=lower --output-legacy-dates --output-field-replace=journaltitle:journal

bib:
	${CMDbib} paper_tdep_joss.bib
	mv paper_tdep_joss_bibertool.bib literature.bib
