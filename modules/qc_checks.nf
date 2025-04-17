// module for other ingest qc checks
process check_single_fastq {

    label "process_low"

    publishDir "${params.outdir}/${unique_id}/preprocess/", mode: "copy", pattern: "*.fixed.*"
    publishDir "${params.outdir}/${unique_id}/preprocess/", mode: "copy", pattern: "*.R*.fastq"

    conda "bioconda::pyfastx=2.01"
    container "biocontainers/pyfastx:2.0.1--py39h3d4b85c_0"

    input:
        tuple val(unique_id), path(fastq)

    output:
        tuple val(unique_id), path("*.fixed.fastq"), optional: true, emit: single_fastq
        tuple val(unique_id), path("*.R1.fastq"), path("*.R2.fastq"), optional: true, emit: paired_fastq

    script:
    """
    check_reads.py --fastq ${fastq}
    """
}


process read_stats {
    
    label "process_medium"

    conda "nanoporetech::fastcat=0.15.1"
    container "${params.wf.container}:${params.wf.container_version}"

    input:
        tuple val(unique_id), path(fastq)
    output:
        tuple val(unique_id), path("${fastq.simpleName}.stats"), emit: stats
    script:
    """
    fastcat -s "${fastq.simpleName}" \
            -r "${fastq.simpleName}.stats" \
            "${fastq}" > filtered.fastq
    """
}

process publish_stats {
    
    label "process_low"

    publishDir "${params.outdir}/${unique_id}/qc", mode: 'copy'

    container "${params.wf.container}:${params.wf.container_version}"

    input:
        tuple val(unique_id), path(stats)
    output:
        tuple val(unique_id), path("raw_reads.stats")
    script:
    """
    mv ${stats} "raw_reads.stats"
    """
}

workflow qc_checks {
    take:
        input_ch
    main:
        if (params.paired) {
            input_ch.map{ unique_id, fastq1, fastq2 -> [unique_id, fastq1, unique_id, fastq2] }
                    .flatten()
                    .collate( 2 )
                    .set{ fastq_ch }
        } else {
            check_single_fastq(input_ch)
            fastq_ch = input_ch
        }
        read_stats(fastq_ch)
        read_stats.out.collectFile(keepHeader:true, skip: 1)
                  .map{ it -> [it.simpleName, it] }
                  .set{ stats_ch }
        publish_stats(stats_ch)
    emit:
        publish_stats.out
}
