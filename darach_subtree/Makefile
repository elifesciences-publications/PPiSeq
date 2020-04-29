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

archive: data.zip-piece-00
data.zip: data
	zip -9 -r $@ $<
data.zip-piece-00: data.zip
	split --bytes=200M -d $< data.zip-piece-
	
data-archive.zip: data.zip-piece-00
	@echo "I'm assuming that if you have the first piece, you have the last"
	@echo "piece. Also, I'm going to rename it as 'data-archive.zip' to not"
	@echo "over write it. So 'data.zip' is what I made, then I split it,"
	@echo "then you're putting those back together into 'data-archive.zip',"
	@echo "and then you unfold it in the below step."
	cat data.zip-piece-* > data-archive.zip

unfold_archive: data-archive.zip
	@echo "ARE YOU SURE ABOUT THIS? THIS unzips `data-archive.zip` to re-write "
	@echo "the `data` folder."
	read something
	@echo "REALLY? YOU better have data-archive.zip here..."
	read something
	unzip data-archive.zip -d data


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
