include { start_server } from '../modules/kraken_server'
include { stop_server } from '../modules/kraken_server'
include { qc_checks } from '../modules/qc_checks'
include { generate_report } from '../modules/generate_report'


kraken_compute = params.threads == 1 ? 1 : params.threads - 1

process kraken2_client {
    label "scylla"
    publishDir path: "${params.out_dir}/${sample_id}/kraken", mode: 'copy'

    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    // retry if server responds out of resource
    errorStrategy = { task.exitStatus in [8] ? 'retry' : 'finish' }
    maxForks kraken_compute
    input:
        val sample_id
        path fastq
    output:
        path "${sample_id}.kraken_report.txt", emit: report
        path "${sample_id}.kraken_assignments.tsv", emit: assignments
    script:
    """
    kraken2_client \
        --port $params.port --report "${sample_id}.kraken_report.txt" \
        --sequence ${fastq} > "${sample_id}.kraken_assignments.tsv"
    """
}

process determine_bracken_length {
    label "scylla"
    input:
        path database
    output:
        env BRACKEN_LENGTH
    """
    if [[ -f "${database}"/database${params.bracken_length}mers.kmer_distrib ]]; then
        BRACKEN_LENGTH="${params.bracken_length}"
    else
        cd "${database}"
        BRACKEN_LENGTH=\$(ls -v1 *.kmer_distrib | tail -1 | sed -e "s@^database@@" -e "s@mers.kmer_distrib@@")
        cd ..
    fi
    """
}

process bracken {
    publishDir path: "${params.out_dir}/${sample_id}/bracken", mode: 'copy'
    input:
        val sample_id
        path kraken_report
        path database
        val bracken_length
    output:
        path "${sample_id}.bracken_summary.txt", emit: summary
        path "${sample_id}.bracken_report.txt", emit: report
    """
    bracken \
      -d "${database}" \
      -i "${kraken_report}" \
      -r "${bracken_length}" \
      -l "${params.bracken_level}" \
      -o "${sample_id}.bracken_summary.txt" \
      -w "${sample_id}.bracken_report.txt"
    """
}

process bracken_to_json {
    publishDir path: "${params.out_dir}/${sample_id}/bracken", mode: 'copy'
    input:
        val sample_id
        path kraken_report
        path taxonomy_dir
        val database_name
        path bracken_summary
    output:
        path "${sample_id}.${database_name}.json"
    """
    cat "${sample_id}.bracken_summary.txt" | cut -f2,6 | tail -n+2 > taxacounts.txt
    cat "${sample_id}.bracken_summary.txt" | cut -f2 | tail -n+2 > taxa.txt
    taxonkit lineage --data-dir ${taxonomy_dir}  -R taxa.txt  > lineages.txt
    $projectDir/../bin/aggregate_lineages_bracken.py \\
            -i "lineages.txt" -b "taxacounts.txt" \\
            -u "${kraken_report}" \\
            -p "${sample_id}"
    file1=`cat *.json`
    echo "{"'"${database_name}"'": "\$file1"}" >> "${sample_id}.${database_name}.json"
    """
}


workflow run_kraken_and_bracken {
    take:
        sample_id
        fastq
        database_name
        database
        taxonomy
    main:
        kraken2_client(sample_id, fastq)
        bracken_length = determine_bracken_length(database)
        bracken(sample_id, kraken2_client.output.report, database, bracken_length)
        bracken_to_json(sample_id, kraken2_client.output.report, taxonomy, database_name, bracken.out.summary)
    emit:
        report = bracken.output.report
        assignments = kraken2_client.output.assignments
        json = bracken_to_json.out
}


