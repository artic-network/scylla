// module for other ingest qc checks

process read_stats {
    
    label "process_medium"

    conda "epi2melabs::kraken2-server=0.1.3"
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

    conda "conda-forge::pandas=1.2.4 conda-forge::numpy=1.20.3"
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

workflow {
    unique_id = "${params.unique_id}"
    fastq = file(params.fastq, type: "file", checkIfExists:true)
    if (unique_id == "null") {
       unique_id = "${fastq.simpleName}"
    }

    qc_ch = [unique_id, fastq]
    qc_checks(qc_ch)
}