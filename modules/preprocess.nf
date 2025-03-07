process check_single_fastq {

    label "process_low"

    publishDir "${params.outdir}/${unique_id}/preprocess/*.fixed.*", mode: "copy"
    publishDir "${params.outdir}/${unique_id}/preprocess/*.R*.fq", mode: "copy"

    conda "bioconda::pyfastx=2.01"
    container "biocontainers/pyfastx:2.0.1--py39h3d4b85c_0"

    input:
        tuple val(unique_id), path(fastq)

    output:
        tuple val(unique_id), path("*.fixed.fq"), optional: true, emit: single_fastq
        tuple val(unique_id), path("*.R1.fq"), path("*.R2.fq"), optional: true, emit: paired_fastq

    script:
    """
    check_reads.py --fastq ${fastq}
    """
}

process fastp_paired {
    
    label "process_medium"

    publishDir "${params.outdir}/${unique_id}/preprocess/", mode: "copy"

    container "${params.wf.container}:${params.wf.container_version}"

    errorStrategy {task.exitStatus in [255, 10] ? "ignore" : "terminate"}

    input:
        tuple val(unique_id), path(fastq_1), path(fastq_2)

    output:
        tuple val(unique_id), path("${unique_id}_1.fastp.fastq.gz"), path("${unique_id}_2.fastp.fastq.gz"), emit: fastq
        path "${unique_id}.fastp.json"

    script:
    """
    fastp \\
        --in1 ${fastq_1} \\
        --in2 ${fastq_2} \\
        --out1 ${unique_id}_1.fastp.fastq \\
        --out2 ${unique_id}_2.fastp.fastq \\
        --json ${unique_id}.fastp.json \\
        --low_complexity_filter \\
        --thread $task.cpus \\
        2> ${unique_id}.fastp.log

    if [ -s ${unique_id}_1.fastp.fastq ]; then
        bgzip --threads $task.cpus -c ${unique_id}_1.fastp.fastq > ${unique_id}_1.fastp.fastq.gz
    else
        exit 10
    fi

    if [ -s ${unique_id}_2.fastp.fastq ]; then
        bgzip --threads $task.cpus -c ${unique_id}_2.fastp.fastq > ${unique_id}_2.fastp.fastq.gz
    else
        exit 10
    fi    
    """

}

process fastp_single {
    
    label "process_medium"

    publishDir "${params.outdir}/${unique_id}/preprocess/", mode: "copy"

    container "${params.wf.container}:${params.wf.container_version}"

    errorStrategy {task.exitStatus in [255, 10] ? "ignore" : "terminate"}

    input:
        tuple val(unique_id), path(fastq)

    output:
        tuple val(unique_id), path("${unique_id}.fastp.fastq.gz"), emit: fastq
        path "${unique_id}.fastp.json"

    script:

    """
    fastp \\
        --in1 ${fastq} \\
        --out1 ${unique_id}.fastp.fastq \\
        --json ${unique_id}.fastp.json \\
        --thread $task.cpus \\
        --disable_adapter_trimming \\
        --low_complexity_filter \\
        --qualified_quality_phred 10 \\
        2> ${unique_id}.fastp.log

    if [ -s ${unique_id}.fastp.fastq ]; then
        bgzip --threads $task.cpus -c ${unique_id}.fastp.fastq > ${unique_id}.fastp.fastq.gz
    else
        exit 10
    fi
    """
}

process paired_concatenate {

    label "process_low"

    errorStrategy {task.exitStatus in [5, 8] ? 'ignore' : 'terminate'}

    publishDir "${params.outdir}/${unique_id}/preprocess/", mode: "copy"

    container "${params.wf.container}:${params.wf.container_version}"

    input:
        tuple val(unique_id), path(processed_fastq_1), path(processed_fastq_2)

    output:
        tuple val(unique_id), path("${unique_id}.concatenated.fastq.gz"), emit: concatenated_fastq

    script:
    """
    concatenate_reads.py --no-interleave \\
        ${processed_fastq_1} ${processed_fastq_2} \\
        --strict \\
        > ${unique_id}.concatenated.fastq

    
    bgzip --threads $task.cpus -c ${unique_id}.concatenated.fastq > ${unique_id}.concatenated.fastq.gz
    """
}

workflow preprocess {
    take:
        unique_id
    main:
        if (params.paired) {
            fastq1 = file(params.fastq1, type: "file", checkIfExists:true)
            fastq2 = file(params.fastq2, type: "file", checkIfExists:true)
            input_ch = Channel.from([[unique_id, fastq1, fastq2]])

            fastp_paired(input_ch)
            fastp_paired.out.fastq
                .set { processed_fastq_ch }

            paired_concatenate(fastp_paired.out.fastq)
            paired_concatenate.out.concatenated_fastq
                .set {combined_fastq_ch}

        } else if (params.fastq) {
            fastq = file(params.fastq, type: "file", checkIfExists:true)
            input_ch = Channel.from([[unique_id, fastq]])

            check_single_fastq(input_ch)
            fastp_single(input_ch)

            fastp_single.out.fastq
                .tap {processed_fastq_ch}
                .set {combined_fastq_ch}

        } else if (params.fastq_dir) {
            fastqdir = file("${params.fastq_dir}", type: "dir", checkIfExists:true)
            Channel.fromPath( fastqdir / "*.f*q*", type: "file")
                .set {fastq_ch}
            input_ch = fastq_ch.map{fastq -> [unique_id, fastq]}

            check_single_fastq(input_ch)
            fastp_single(input_ch)
            fastp_single.out.fastq.map{ unique_id, fastq -> [unique_id + ".fq.gz", fastq]}
                .collectFile()
                .map{ it -> [it.simpleName, it] }
                .tap {processed_fastq_ch}
                .set {combined_fastq_ch}
        }
    emit:
        processed_fastq = processed_fastq_ch
        combined_fastq = combined_fastq_ch
}
