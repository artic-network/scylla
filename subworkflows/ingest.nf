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
        get_params_and_versions()
        //clean_metadata(metadata)
        kraken_pipeline(unique_id, fastq)
        extract_taxa(unique_id, fastq, kraken_pipeline.out.kraken_assignments, kraken_pipeline.out.kraken_report, kraken_pipeline.out.bracken_report)
    emit:
        html_report = kraken_pipeline.out.report

}

workflow {
    // check input fastq exists
    input_fastq = file("${params.fastq}", type: "file", checkIfExists:true)

    unique_id = "${params.unique_id}"
    if (unique_id == "null") {
        unique_id = "${input_fastq.simpleName}"
    }

    ingest(unique_id, input_fastq)
}