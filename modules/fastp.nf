process fastp_paired {
    
    label 'process_medium'

    publishDir "${params.out_dir}/${unique_id}/preprocess/", mode: 'copy'

    container "biowilko/scylla@${params.wf.container_sha}"

    input:
        val unique_id
        val prefix
        path fastq_1
        path fastq_2

    output:
        val prefix
        path "${prefix}_1.fastp.fastq.gz", optional: true
        path "${prefix}_2.fastp.fastq.gz", optional: true
        path "${prefix}.fastp.fastq.gz", emit: processed_fastq
        path "${prefix}.fastp.json"

    script:
    """
    fastp \\
        --in1 ${fastq_1} \\
        --in2 ${fastq_2} \\
        --out1 ${prefix}_1.fastp.fastq \\
        --out2 ${prefix}_2.fastp.fastq \\
        --merged_out ${prefix}.fastp.fastq \\
        -m \\
        --json ${prefix}.fastp.json \\
        --thread $task.cpus \\
        2> ${prefix}.fastp.log
    
    if [ -s ${prefix}_1.fastp.fastq ]; then
        bgzip --threads $task.cpus -c ${prefix}_1.fastp.fastq > ${prefix}_1.fastp.fastq.gz
    fi

    if [ -s ${prefix}_2.fastp.fastq ]; then
        bgzip --threads $task.cpus -c ${prefix}_2.fastp.fastq > ${prefix}_2.fastp.fastq.gz
    fi

    if [ -s ${prefix}.fastp.fastq ]; then
        bgzip --threads $task.cpus -c ${prefix}.fastp.fastq > ${prefix}.fastp.fastq.gz
    fi

    """

}

process fastp_single {
    
    label 'process_medium'

    publishDir "${params.out_dir}/${unique_id}/preprocess/", mode: 'copy'

    container "biowilko/scylla@${params.wf.container_sha}"

    input:
        val unique_id
        val prefix
        path fastq

    output:
        val prefix
        path "${prefix}.fastp.fastq.gz", emit: processed_fastq
        path "${prefix}.fastp.json"

    script:

    """
    fastp \\
        --in1 ${fastq} \\
        --out1 ${prefix}.fastp.fastq \\
        --json ${prefix}.fastp.json \\
        --thread $task.cpus \\
        2> ${prefix}.fastp.log

    if [ -s ${prefix}.fastp.fastq ]; then
        bgzip --threads $task.cpus -c ${prefix}.fastp.fastq > ${prefix}.fastp.fastq.gz
    fi
    """
}