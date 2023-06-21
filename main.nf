include { ingest } from './subworkflows/ingest'

workflow {
    unique_id = "${params.unique_id}"

    // check input fastq exists and run fastp

    ingest(params.unique_id)
}