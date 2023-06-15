include { ingest } from './subworkflows/ingest'
include { fastp_single; fastp_paired} from '../modules/fastp'

workflow {
    // check input fastq exists and run fastp
    
    if (params.paired) {
        Channel.of(file(params.fastq1, type: "file", checkIfExists:true))
            .set {input_fastq_1_ch}

        Channel.of(file(params.fastq2, type: "file", checkIfExists:true))
            .set {input_fastq_2_ch}

        fastp_paired(params.unique_id, input_fastq_1_ch, input_fastq_2_ch)
        fastp_paired.out.processed_fastq
            .set {processed_fastq}
    } else {
        Channel.of(file(params.fastq, type: "file", checkIfExists:true))
            .set {input_fastq_ch}

        fastp_single(params.unique_id, input_fastq_ch)
        fastp_single.out.processed_fastq
            .set {processed_fastq}
    }

    ingest(params.unique_id, processed_fastq)
}