// module to make a single sample report during ingest containing kraken summary and qc stats
process make_report {
    
    label "process_low"

    publishDir path: "${params.out_dir}/${unique_id}/", mode: 'copy'
    
    conda "bioconda::biopython=1.78 anaconda::Mako=1.2.3"
    container "biowilko/scylla@${params.wf.container_sha}"
    
    input:
        val unique_id
        path stats
        path lineages
        path template
    output:
        path "${unique_id}_report.html", emit: report_html
    script:
        report_name = "${unique_id}"
    """
    make_report.py \
        --prefix "${report_name}" \
        --read_counts ${stats} \
        --assignments ${lineages} \
        --version "${workflow.manifest.version}" \
        --template "${template}"
    """
}

workflow generate_report {
    take:
        unique_id
        stats
        bracken_jsons
    main:
        // Acquire report template
        template = file("$baseDir/bin/scylla.mako.html")

        make_report(unique_id, stats, bracken_jsons, template)
    emit:
        make_report.out
}
