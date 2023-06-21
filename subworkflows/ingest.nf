// workflow to run kraken, check for human, run qc checks and generate html report for a single sample fastq
include { get_params_and_versions } from '../modules/get_params_and_versions'
include { kraken_pipeline } from '../subworkflows/kraken_pipeline'
include { extract_taxa } from '../modules/extract_taxa'

workflow ingest {
    take:
        unique_id
        fastq
        //metadata
    main:
        // get_params_and_versions()
        //clean_metadata(metadata)
        kraken_pipeline(unique_id, fastq)

        if (params.paired) {
                fastq_1 = file(params.fastq1, type: "file", checkIfExists:true)

                fastq_2 = file(params.fastq2, type: "file", checkIfExists:true)
        } else { 
                fastq_1 = file(params.fastq, type: "file", checkIfExists:true)
                fastq_2 = None
        }

        extract_taxa(unique_id, fastq_1, fastq_2, kraken_pipeline.out.kraken_assignments, kraken_pipeline.out.kraken_report, kraken_pipeline.out.bracken_report)
    emit:
        html_report = kraken_pipeline.out.report

}
