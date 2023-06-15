// workflow to run kraken, check for human, run qc checks and generate html report for a single sample fastq
include { get_params_and_versions } from '../modules/get_params_and_versions'
include { kraken_pipeline } from '../subworkflows/kraken_pipeline'
include { extract_taxa } from '../modules/extract_taxa'
include { fastp_single; fastp_paired} from '../modules/fastp'

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
    // check input fastq exists and run fastp
    
    if (params.paired) {
        input_fastq_1 = file(params.fastq1, type: "file", checkIfExists:true)
        input_fastq_2 = file(params.fastq2, type: "file", checkIfExists:true)
        fastp_paired(params.unique_id, input_fastq_1, input_fastq_2)
        fastp_paired.out.processed_fastq
            .set {processed_fastq}
    } else {
        input_fastq = file(params.fastq, type: "file", checkIfExists:true)
        fastp_single(params.unique_id, input_fastq)
        fastp_single.out.processed_fastq
            .set {processed_fastq}
    }

    ingest(unique_id, processed_fastq)
}