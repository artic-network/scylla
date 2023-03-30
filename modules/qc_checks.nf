// module for other ingest qc checks

process read_stats {
    publishDir "${params.out_dir}/${unique_id}/qc", mode: 'copy'
    input:
        val unique_id
        path fastq
    output:
        path "raw_reads.stats", emit: stats
        path "raw_reads.json", emit: json
    script:
    """
    fastcat -s "raw_reads" \
            -r "raw_reads.stats" \
            "${fastq}" > filtered.fastq
    $projectDir/../bin/fastcat_histogram.py \
            --sample_id "raw_reads" \
            "raw_reads.stats" "raw_reads.json"
    """
}

workflow qc_checks {
    take:
        unique_id
        fastq
    main:
        read_stats(unique_id, fastq)

    emit:
        read_stats.out.json
}

workflow {
    qc_checks(params.unique_id, params.fastq)
}