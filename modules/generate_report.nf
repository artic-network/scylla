// module to make a single sample report during ingest containing kraken summary and qc stats
process make_report {
    publishDir path: "${params.out_dir}/${unique_id}", mode: 'copy'
    maxForks 1
    cpus 1
    input:
        val unique_id
        path stats
        path lineages
        path template
    output:
        path "scylla-report.html", emit: report_html
    script:
        report_name = "scylla-report.html"
    """
    $projectDir/../bin/single_sample_report.py \
        "${report_name}" \
        --stats ${stats} \
        --lineages ${lineages} \
        --report_template "${template}"
    """
}

workflow generate_report {
    take:
        unique_id
        stats
        bracken_jsons
    main:
        // Acquire report template
        template = file("$projectDir/../bin/report-visualisation.html")

        make_report(unique_id, stats, bracken_jsons, template)
    emit:
        make_report.out
}
