process fastp_paired {
    tag "$meta.id"
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
        path("${prefix}_1.fastp.fastq.gz")
        path("${prefix}_2.fastp.fastq.gz")
        path("${prefix}.fastp.fastq.gz"), emit: processed_fastq
        path("${prefix}.fastp.json")

    script:
    """
    fastp \\
        --in1 ${fastq_1} \\
        --in2 ${fastq_2} \\
        --out1 ${prefix}_1.fastp.fastq.gz \\
        --out2 ${prefix}_2.fastp.fastq.gz \\
        --merged-out ${prefix}.fastp.fastq.gz \\
        --json ${prefix}.fastp.json \\
        --thread $task.cpus \\
        --detect_adapter_for_pe \\
        2> ${prefix}.fastp.log
    """

}

process fastp_single {
    tag "$meta.id"
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
        --out1 ${prefix}.fastp.fastq.gz \\
        --json ${prefix}.fastp.json \\
        --thread $task.cpus \\
        2> ${prefix}.fastp.log
    """
}