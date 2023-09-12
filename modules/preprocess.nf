process fastp_paired {
    
    label 'process_medium'

    publishDir "${params.outdir}/${unique_id}/preprocess/", mode: 'copy'

    container "biowilko/scylla@${params.wf.container_sha}"

    input:
        val unique_id
        path fastq_1 
        path fastq_2

    output:
        val unique_id
        path "${unique_id}_1.fastp.fastq.gz", emit: processed_fastq_1
        path "${unique_id}_2.fastp.fastq.gz", emit: processed_fastq_2
        path "${unique_id}.fastp.json"

    script:
    """
    fastp \\
        --in1 ${fastq_1} \\
        --in2 ${fastq_2} \\
        --out1 ${unique_id}_1.fastp.fastq \\
        --out2 ${unique_id}_2.fastp.fastq \\
        --json ${unique_id}.fastp.json \\
        --thread $task.cpus \\
        2> ${unique_id}.fastp.log

    if [ -s ${unique_id}_1.fastp.fastq ]; then
        bgzip --threads $task.cpus -c ${unique_id}_1.fastp.fastq > ${unique_id}_1.fastp.fastq.gz
    fi

    if [ -s ${unique_id}_2.fastp.fastq ]; then
        bgzip --threads $task.cpus -c ${unique_id}_2.fastp.fastq > ${unique_id}_2.fastp.fastq.gz
    fi    
    """

}

process fastp_single {
    
    label 'process_medium'

    publishDir "${params.outdir}/${unique_id}/preprocess/", mode: 'copy'

    container "biowilko/scylla@${params.wf.container_sha}"

    input:
        val unique_id
        path fastq

    output:
        val unique_id
        path "${unique_id}.fastp.fastq.gz", emit: processed_fastq
        path "${unique_id}.fastp.json"

    script:

    """
    fastp \\
        --in1 ${fastq} \\
        --out1 ${unique_id}.fastp.fastq \\
        --json ${unique_id}.fastp.json \\
        --thread $task.cpus \\
        2> ${unique_id}.fastp.log

    if [ -s ${unique_id}.fastp.fastq ]; then
        bgzip --threads $task.cpus -c ${unique_id}.fastp.fastq > ${unique_id}.fastp.fastq.gz
    fi
    """
}

process paired_concatenate {

    label 'process_low'

    publishDir "${params.outdir}/${unique_id}/preprocess/", mode: 'copy'

    container "biowilko/scylla@${params.wf.container_sha}"

    input:
        val unique_id
        path processed_fastq_1
        path processed_fastq_2

    output:
        val unique_id
        path "${unique_id}.concatenated.fastq.gz", emit: concatenated_fastq

    script:
    """
    concatenate_reads.py --no-interleave \\
        ${processed_fastq_1} ${processed_fastq_2} \\
        > ${unique_id}.concatenated.fastq

    
    bgzip --threads $task.cpus -c ${unique_id}.concatenated.fastq > ${unique_id}.concatenated.fastq.gz
    """
}