process fastp_paired {
    
    label 'process_medium'

    publishDir "${params.out_dir}/preprocess/", mode: 'copy'

    conda "bioconda::fastp=0.23.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0' :
        'biocontainers/fastp:0.23.4--h5f740d0_0' }"

    input:
        val(prefix)
        path(fastq_1)
        path(fastq_2)

    output:
        val(prefix)
        path("${prefix}_1.fastp.fastq.gz"), optional: true
        path("${prefix}_2.fastp.fastq.gz"), optional: true
        path("${prefix}.fastp.fastq.gz"), emit: processed_fastq
        path("${prefix}.fastp.json")

    script:
    """
    fastp \\
        --in1 ${fastq_1} \\
        --in2 ${fastq_2} \\
        --out1 ${prefix}_1.fastp.fastq \\
        --out2 ${prefix}_2.fastp.fastq \\
        --merged-out ${prefix}.fastp.fastq \\
        --json ${prefix}.fastp.json \\
        --thread $task.cpus \\
        --detect_adapter_for_pe \\
        2> ${prefix}.fastp.log
    
    if [ -s ${prefix}_1.fastp.fastq ]; then
        bgzip --threads -c $task.cpus > ${prefix}_1.fastp.fastq.gz
    fi

    if [ -s ${prefix}_2.fastp.fastq ]; then
        bgzip --threads -c $task.cpus > ${prefix}_2.fastp.fastq.gz
    fi

    if [ -s ${prefix}.fastp.fastq ]; then
        bgzip --threads -c $task.cpus > ${prefix}.fastp.fastq.gz
    fi

    """

}

process fastp_single {
    
    label 'process_medium'

    publishDir "${params.out_dir}/preprocess/", mode: 'copy'

    conda "bioconda::fastp=0.23.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0' :
        'biocontainers/fastp:0.23.4--h5f740d0_0' }"

    input:
        val(prefix)
        path(fastq)

    output:
        val(prefix)
        path("${prefix}.fastp.fastq.gz"), emit: processed_fastq
        path("${prefix}.fastp.json")

    script:

    """
    fastp \\
        --in1 ${fastq} \\
        --out1 ${prefix}.fastp.fastq \\
        --json ${prefix}.fastp.json \\
        --thread $task.cpus \\
        2> ${prefix}.fastp.log

    if [ -s ${prefix}.fastp.fastq ]; then
        bgzip --threads -c $task.cpus > ${prefix}.fastp.fastq.gz
    fi
    """
}