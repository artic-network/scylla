process fastp_paired {

    label "process_medium"

    publishDir "${params.outdir}/${unique_id}/preprocess/", mode: params.publish_dir_mode

    container "${params.wf.container}:${params.wf.container_version}"

    errorStrategy { task.exitStatus in [255, 10] ? "ignore" : "terminate" }

    input:
    tuple val(unique_id), path(fastq_1), path(fastq_2)

    output:
    tuple val(unique_id), path("${unique_id}_1.fastp.fastq.gz"), path("${unique_id}_2.fastp.fastq.gz"), emit: fastq
    path "${unique_id}.fastp.json"

    script:
    """
    mkfifo fastp_out_1 fastp_out_2
    crabz -p ${task.cpus} -f gzip -o ${unique_id}_1.fastp.fastq.gz fastp_out_1 &
    BG1=\$!
    crabz -p ${task.cpus} -f gzip -o ${unique_id}_2.fastp.fastq.gz fastp_out_2 &
    BG2=\$!

    # `|| FASTP_STATUS=\$?` keeps the failure in-band: process.shell sets
    # `bash -euo pipefail`, so a bare `fastp` failure would abort the script
    # right here, leaving the crabz readers above orphaned and blocked on their
    # FIFOs. A command on the left of `||` is exempt from `set -e`.
    FASTP_STATUS=0
    fastp \\
        --in1 ${fastq_1} \\
        --in2 ${fastq_2} \\
        --out1 fastp_out_1 \\
        --out2 fastp_out_2 \\
        --json ${unique_id}.fastp.json \\
        --low_complexity_filter \\
        --thread ${task.cpus} \\
        2> ${unique_id}.fastp.log || FASTP_STATUS=\$?

    if [ "\$FASTP_STATUS" -ne 0 ]; then
        # fastp never opened - or stopped writing to - its output FIFOs, so the
        # crabz readers would wait for a writer forever.
        kill \$BG1 \$BG2 2>/dev/null || true
        wait \$BG1 2>/dev/null || true
        wait \$BG2 2>/dev/null || true
        exit \$FASTP_STATUS
    fi

    # fastp succeeded, so both crabz readers definitely had a writer: any
    # non-zero status here is a real compression failure, and swallowing it would
    # publish a truncated .fastq.gz as a successful task, since the read count
    # check below only inspects fastp's own JSON.
    CRABZ_STATUS=0
    wait \$BG1 || CRABZ_STATUS=\$?
    wait \$BG2 || CRABZ_STATUS=\$?
    if [ "\$CRABZ_STATUS" -ne 0 ]; then
        echo "ERROR: compression of fastp output failed (crabz exit \$CRABZ_STATUS)" >&2
        exit 11
    fi

    if [ -s ${unique_id}.fastp.json ]; then
        READS=\$(jq '.filtering_result.passed_filter_reads' ${unique_id}.fastp.json)
    else
        READS=0
    fi
    if [ "\$READS" -eq 0 ]; then
        exit 10
    fi
    """
}

process fastp_single {

    label "process_medium"

    publishDir "${params.outdir}/${unique_id}/preprocess/", mode: params.publish_dir_mode

    container "${params.wf.container}:${params.wf.container_version}"

    errorStrategy { task.exitStatus in [255, 10] ? "ignore" : "terminate" }

    input:
    tuple val(unique_id), path(fastq)

    output:
    tuple val(unique_id), path("${unique_id}.fastp.fastq.gz"), emit: fastq
    path "${unique_id}.fastp.json"

    script:

    """
    set -o pipefail
    fastp \\
        --in1 ${fastq} \\
        --stdout \\
        --json ${unique_id}.fastp.json \\
        --thread ${task.cpus} \\
        --disable_adapter_trimming \\
        --low_complexity_filter \\
        --qualified_quality_phred 10 \\
        2> ${unique_id}.fastp.log \\
        | crabz -p ${task.cpus} -f gzip -o ${unique_id}.fastp.fastq.gz

    READS=\$(jq '.filtering_result.passed_filter_reads' ${unique_id}.fastp.json)
    if [ "\$READS" -eq 0 ]; then
        exit 10
    fi
    """
}

process paired_concatenate {

    label "process_low"

    errorStrategy { task.exitStatus in [5, 8] ? 'ignore' : 'terminate' }

    publishDir "${params.outdir}/${unique_id}/preprocess/", mode: params.publish_dir_mode

    container "${params.wf.container}:${params.wf.container_version}"

    input:
    tuple val(unique_id), path(processed_fastq_1), path(processed_fastq_2)

    output:
    tuple val(unique_id), path("${unique_id}.concatenated.fastq.gz"), emit: concatenated_fastq

    script:
    """
    set -o pipefail
    concatenate_reads.py --no-interleave \\
        ${processed_fastq_1} ${processed_fastq_2} \\
        --strict \\
        | crabz -p ${task.cpus} -f gzip -o ${unique_id}.concatenated.fastq.gz
    """
}

workflow preprocess {
    take:
    unique_id

    main:
    if (params.paired) {
        fastq1 = file(params.fastq1, type: "file", checkIfExists: true)
        fastq2 = file(params.fastq2, type: "file", checkIfExists: true)
        input_ch = Channel.from([[unique_id, fastq1, fastq2]])

        fastp_paired(input_ch)
        fastp_paired.out.fastq.set { processed_fastq_ch }

        paired_concatenate(fastp_paired.out.fastq)
        paired_concatenate.out.concatenated_fastq.set { combined_fastq_ch }
    }
    else if (params.fastq) {
        fastq = file(params.fastq, type: "file", checkIfExists: true)
        input_ch = Channel.from([[unique_id, fastq]])

        fastp_single(input_ch)

        fastp_single.out.fastq
            .tap { processed_fastq_ch }
            .set { combined_fastq_ch }
    }
    else if (params.fastq_dir) {
        fastqdir = file("${params.fastq_dir}", type: "dir", checkIfExists: true)
        Channel.fromPath(fastqdir / "*.f*q*", type: "file")
            .set { fastq_ch }
        input_ch = fastq_ch.map { fastq -> [unique_id, fastq] }

        fastp_single(input_ch)
        fastp_single.out.fastq
            .map { unique_id, fastq -> [unique_id + ".fq.gz", fastq] }
            .collectFile()
            .map { it -> [it.simpleName, it] }
            .tap { processed_fastq_ch }
            .set { combined_fastq_ch }
    }

    emit:
    processed_fastq = processed_fastq_ch
    combined_fastq  = combined_fastq_ch
}
