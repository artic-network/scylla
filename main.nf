include { ingest } from './subworkflows/ingest'

workflow {
    unique_id = "${params.unique_id}"

    if (unique_id == "null") {
        if (params.fastq) {
            unique_id = "${fastq.simpleName}"
        } else if (params.fastq_dir) {
            unique_id = "${fastq_dir.simpleName}"
        } else if (params.paired & params.fastq1 & params.fastq2) {
            unique_id = "${fastq1.simpleName}"
        } else {
            exit 1, "One of fastq, fastq_dir or fastq1 and fastq2 need to be provided -- aborting"
        }
    }

    // check input fastq exists and run fastp
    ingest(unique_id)
}