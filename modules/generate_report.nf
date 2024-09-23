// module to make a single sample report during ingest containing kraken summary and qc stats
process make_report {
    
    label "process_single"
    label "process_more_memory"

    errorStrategy "retry"
    maxRetries 3

    publishDir "${params.outdir}/${unique_id}/", mode: 'copy'
    
    conda "anaconda::Mako=1.2.3"
    container "${params.wf.container}:${params.wf.container_version}"
    
    input:
        tuple val(unique_id), path(stats), path(lineages), path(warnings)
        path template
    output:
        tuple val(unique_id), path("${unique_id}_report.html")
    script:
        report_name = "${unique_id}"
        if ( params.run_sourmash ){
            classifier = "Sourmash"
            classification_database = "${params.sourmash_db_name}"
        } else {
            classifier = "Kraken"
            classification_database = "${params.database_set}"
        }
    """
    make_report.py \
        --prefix "${report_name}" \
        --read_counts ${stats} \
        --assignments ${lineages} \
        --warnings ${warnings} \
        --version "${workflow.manifest.version}" \
        --classifier "${classifier}" \
        --classification_database "${classification_database}" \
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
