.PHONY: all setup archive \
	nextflow \
	tecan_validations \
	networks \
	clean

scripts/nextflow: 
	curl -s https://get.nextflow.io | bash
	mv nextflow scripts/nextflow

nextflow: scripts/nextflow \
		scripts/run_pipeline.nf \
		scripts/run_pipeline.nfconfig
	unset SINGULARITY_CACHEDIR && \
	$< run scripts/run_pipeline.nf \
		-c scripts/run_pipeline.nfconfig \
		-resume -ansi-log false \
		-with-dag reports/dag.html

archive: data.zip
data.zip: data
	zip -9 -r $@ $<
	
unfold_archive: 
	@echo "ARE YOU SURE ABOUT THIS? THIS RE-WRITES THE data FOLDER FROM data.zip"
	read something
	@echo "REALLY? YOU better have data.zip here..."
	read something
	unzip data.zip 

all: 
	@echo "none"

lightclean:
	@echo "do a clean? if not cancel cntrl D"
	@read DO_IT
	rm -rf scripts/.ipynb_checkpoints
	rm -rf scripts/*.png
	rm -rf scripts/*.pdf
	rm -rf scripts/*.jpeg
	rm -rf scripts/*_cache
	rm -rf scripts/*_files

clean: lightclean
	@echo "do a clean? if not cancel cntrl D"
	@read DO_IT
	rm -rf tmp
	rm -rf .nextflow
	rm -rf work
	rm -rf reports
