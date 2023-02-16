// module for other ingest qc checks

process read_stats {
    publishDir "${params.out_dir}/${sample_id}/stats", mode: 'copy'
    input:
        val sample_id
        path fastq
    output:
        path "${sample_id}.stats", emit: stats
        path "${sample_id}.json", emit: json
    script:
    """
    fastcat -s "${sample_id}" \
            -r "${sample_id}.stats" \
            "${fastq}" > filtered.fastq
    $projectDir/../bin/fastcat_histogram.py \
            --sample_id "${sample_id}" \
            "${sample_id}.stats" "${sample_id}.json"
    """
}

workflow qc_checks {
    take:
        sample_id
        fastq
    main:
        read_stats(sample_id, fastq)

    emit:
        read_stats.out.json
}

workflow {
    qc_checks()
}