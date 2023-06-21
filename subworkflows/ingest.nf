// workflow to run kraken, check for human, run qc checks and generate html report for a single sample fastq
include { get_params_and_versions } from '../modules/get_params_and_versions'
include { kraken_pipeline } from '../subworkflows/kraken_pipeline'
include { extract_reads } from '../modules/extract_taxa'

workflow ingest {
    take:
        unique_id
        //metadata
    main:

        if (params.paired) {
            Channel.of(file(params.fastq1, type: "file", checkIfExists:true))
                .set {input_fastq_1_ch}

            Channel.of(file(params.fastq2, type: "file", checkIfExists:true))
                .set {input_fastq_2_ch}

            fastp_paired(unique_id, input_fastq_1_ch, input_fastq_2_ch)

            paired_concatenate(unique_id, fastp_paired.out.processed_fastq_1, fastp_paired.out.processed_fastq_2)

            paired_concatenate.out.concatenated_fastq
                .set {processed_fastq}
        } else {
            Channel.of(file(params.fastq, type: "file", checkIfExists:true))
                .set {input_fastq_ch}

            fastp_single(unique_id, input_fastq_ch)
            fastp_single.out.processed_fastq
            .set {processed_fastq}
        }

        // get_params_and_versions()
        //clean_metadata(metadata)
        kraken_pipeline(unique_id, processed_fastq)

        if (params.paired) {
                fastq_1 = fastp_paired.out.processed_fastq_1

                fastq_2 = fastp_paired.out.processed_fastq_2
        } else { 
                fastq_1 = proccessed_fastq
                fastq_2 = None
        }

        extract_reads(unique_id, fastq_1, fastq_2, kraken_pipeline.out.kraken_assignments, kraken_pipeline.out.kraken_report, kraken_pipeline.out.bracken_report)
    emit:
        html_report = kraken_pipeline.out.report

}
