// workflow to run kraken, check for human, run qc checks and generate html report for a single sample fastq
include { get_params_and_versions } from '../modules/get_params_and_versions'

workflow ingest {
    take:
        fastq
        metadata
        database_set
    main:
        get_params_and_versions()
        //clean_metadata(metadata)
        //kraken2_pipeline(fastq, database_set)
        //dehumanize(fastq,kraken_assignments)
        //qc_checks(dehumanize.output)
        //denovo_assemble(fastq, kraken_report)
        //generate_report(fastq,metadata,kraken_report)
    emit:
        clean_fastq
        clean_metadata
        html_report

}