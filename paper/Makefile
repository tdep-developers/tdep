CMDbib = biber --tool --output_align --output_indent=2 --output_fieldcase=lower --output-legacy-dates --output-field-replace=journaltitle:journal

final:
	${CMDbib} literature_final.bib
	mv literature_final_bibertool.bib literature_final.bib

bib:
	${CMDbib} paper_tdep_joss.bib
	mv paper_tdep_joss_bibertool.bib literature.bib
