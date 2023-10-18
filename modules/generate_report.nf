// module to make a single sample report during ingest containing kraken summary and qc stats
process make_report {
    
    label "process_low"

    publishDir path: "${params.outdir}/${unique_id}/", mode: 'copy'
    
    conda "bioconda::biopython=1.78 anaconda::Mako=1.2.3"
    container "${params.wf.container}@${params.wf.container_sha}"
    
    input:
        tuple val(unique_id), path(stats), path(lineages)
        path template
    output:
        tuple val(unique_id), path("${unique_id}_report.html")
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
        report_ch
    main:
        // Acquire report template
        template = file("$baseDir/bin/scylla.mako.html")

        make_report(report_ch, template)
    emit:
        make_report.out
}
