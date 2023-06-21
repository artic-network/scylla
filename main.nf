include { ingest } from './subworkflows/ingest'
include { fastp_single; fastp_paired; paired_concatenate } from './modules/preprocess'

workflow {
    unique_id = "${params.unique_id}"

    // check input fastq exists and run fastp
    
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

    ingest(params.unique_id, processed_fastq)
}