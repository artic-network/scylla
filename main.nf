include { ingest } from './subworkflows/ingest'
include { fastp_single; fastp_paired; paired_concatenate } from './modules/preprocess'

workflow {
    unique_id = "${params.unique_id}"

    // check input fastq exists and run fastp

    ingest(params.unique_id)
}