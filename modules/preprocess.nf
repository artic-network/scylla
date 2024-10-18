process fastp_paired {
    
    label "process_medium"

    publishDir "${params.outdir}/${unique_id}/preprocess/", mode: "copy"

    container "${params.wf.container}:${params.wf.container_version}"

    errorStrategy {task.exitStatus in [255, 10] ? "ignore" : "terminate"}

    input:
        val unique_id
        path fastq_1 
        path fastq_2

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
        val unique_id
        path fastq

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
            Channel.of(file(params.fastq1, type: "file", checkIfExists:true))
                .set {input_fastq_1_ch}

            Channel.of(file(params.fastq2, type: "file", checkIfExists:true))
                .set {input_fastq_2_ch}

            fastp_paired(unique_id, input_fastq_1_ch, input_fastq_2_ch)
            fastp_paired.out.fastq
                .set { processed_fastq_ch }

            paired_concatenate(fastp_paired.out.fastq)

            paired_concatenate.out.concatenated_fastq
                .set {combined_fastq_ch}
        } else if (params.fastq) {
            Channel.of(file(params.fastq, type: "file", checkIfExists:true))
                .set {input_fastq_ch}

            fastp_single(unique_id, input_fastq_ch)
            fastp_single.out.fastq
                .tap {processed_fastq_ch}
                .set {combined_fastq_ch}
        } else if (params.fastq_dir) {
            fastqdir = file("${params.fastq_dir}", type: "dir", checkIfExists:true)
            Channel.fromPath( fastqdir / "*.f*q*", type: "file")
                .set {input_fastq_ch}

            fastp_single(unique_id, input_fastq_ch)
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

workflow {
    unique_id = "${params.unique_id}"

    preprocess(unique_id)
}