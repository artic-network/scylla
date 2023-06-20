// module for other ingest qc checks

process read_stats {
    
    label "process_medium"

    conda "epi2melabs::kraken2-server=0.1.3"
    container "biowilko/scylla@${params.wf.container_sha}"

    input:
        path fastq
    output:
        path "${fastq.simpleName}.stats", emit: stats
    script:
    """
    fastcat -s "${fastq.simpleName}" \
            -r "${fastq.simpleName}.stats" \
            "${fastq}" > filtered.fastq
    """
}

process combine_stats {
    
    label "process_low"

    publishDir "${params.out_dir}/${unique_id}/qc", mode: 'copy'

    conda "conda-forge::pandas=1.2.4 conda-forge::numpy=1.20.3"
    container "biowilko/scylla@${params.wf.container_sha}"
    

    input:
        val unique_id
        path stats
    output:
        path "raw_reads.stats", emit: stats
        path "raw_reads.json", emit: json
    script:
    """
    cp ${stats} "raw_reads.stats"
    fastcat_histogram.py \
            --sample_id "raw_reads" \
            "raw_reads.stats" "raw_reads.json"
    """
}

workflow qc_checks {
    take:
        unique_id
        fastq
    main:
        read_stats(fastq)
        all_stats = read_stats.out.collectFile(name:"all_stats.txt", keepHeader:true, skip: 1)
        combine_stats(unique_id, all_stats)
    emit:
        combine_stats.out.json
}

workflow {
    qc_checks(params.unique_id, params.fastq)
}