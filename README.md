# scylla
single sample processing

NB about initial files and pipelines:
 - the kraken_server and kraken_client modules are placeholders, but they allow files to be produced locally for testing on
 - this pipeline is supposed to run on single samples but as a result I have not been careful about e.g. whether channels will do funny things if given a channel with multiple fastq
 - I want to be able to run kraken with different databases and to collect the results into a single report html
 - Need to think further about the directory structure for results
 - no files/modules/pipelines are complete
 - need to improve on the files modified from the wf-metagenomics pipeline
 - nextflow run ./subworkflows/kraken_pipeline.nf --fastq test/test_data/barcode02/barcode02.fq works
