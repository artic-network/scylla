include { start_server } from '../modules/kraken_server'
include { stop_server } from '../modules/kraken_server'
include { qc_checks } from '../modules/qc_checks'
include { generate_report } from '../modules/generate_report'


kraken_compute = params.threads == 1 ? 1 : params.threads - 1

process kraken2_client {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 5
    label "scylla"
    publishDir path: "${params.out_dir}/${unique_id}/classifications", mode: 'copy'

    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    // retry if server responds out of resource
    errorStrategy = { task.exitStatus in [8] ? 'retry' : 'finish' }
    maxForks kraken_compute
    input:
        val unique_id
        path fastq
    output:
        path "${params.database_set}.kraken_report.txt", emit: report
        path "${params.database_set}.kraken_assignments.tsv", emit: assignments
    script:
    """
    kraken2_client \
        --port $params.port --report "${params.database_set}.kraken_report.txt" \
        --sequence ${fastq} > "${params.database_set}.kraken_assignments.tsv"
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
    publishDir path: "${params.out_dir}/${unique_id}/classifications", mode: 'copy'
    input:
        val unique_id
        path kraken_report
        path database
        val bracken_length
    output:
        path "${params.database_set}.bracken_summary.txt", emit: summary
        path "${params.database_set}.bracken_report.txt", emit: report
    """
    bracken \
      -d "${database}" \
      -i "${kraken_report}" \
      -r "${bracken_length}" \
      -l "${params.bracken_level}" \
      -o "${params.database_set}.bracken_summary.txt" \
      -w "${params.database_set}.bracken_report.txt"
    """
}

process bracken_to_json {
    publishDir path: "${params.out_dir}/${unique_id}/classifications", mode: 'copy'
    input:
        val unique_id
        path kraken_report
        path taxonomy_dir
        path bracken_summary
    output:
        path "${params.database_set}.bracken.json"
    """
    cat "${bracken_summary}" | cut -f2,6 | tail -n+2 > taxacounts.txt
    cat "${bracken_summary}" | cut -f2 | tail -n+2 > taxa.txt
    taxonkit lineage --data-dir ${taxonomy_dir}  -R taxa.txt  > lineages.txt
    $projectDir/../bin/aggregate_lineages_bracken.py \\
            -i "lineages.txt" -b "taxacounts.txt" \\
            -u "${kraken_report}" \\
            -p "temp_bracken"
    file1=`cat *.json`
    echo "{"'"${params.database_set}"'": "\$file1"}" >> "${params.database_set}.bracken.json"
    """
}


workflow run_kraken_and_bracken {
    take:
        unique_id
        fastq
        database
        taxonomy
    main:
        kraken2_client(unique_id, fastq)
        bracken_length = determine_bracken_length(database)
        bracken(unique_id, kraken2_client.output.report, database, bracken_length)
        bracken_to_json(unique_id, kraken2_client.output.report, taxonomy, bracken.out.summary)
    emit:
        bracken_report = bracken.output.report
        kraken_report = kraken2_client.output.report
        kraken_assignments = kraken2_client.output.assignments
        json = bracken_to_json.out
}


