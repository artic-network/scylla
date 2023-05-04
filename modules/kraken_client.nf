include { start_server } from '../modules/kraken_server'
include { stop_server } from '../modules/kraken_server'
include { qc_checks } from '../modules/qc_checks'
include { generate_report } from '../modules/generate_report'


kraken_compute = params.threads == 1 ? 1 : params.threads - 1

process kraken2_client {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 5
    label "scylla"

    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    // retry if server responds out of resource
    errorStrategy = { task.exitStatus in [8] ? 'retry' : 'finish' }
    maxForks kraken_compute
    input:
        path fastq
    output:
        path "${fastq.baseName}.kraken_report.txt", emit: report
        path "${fastq.baseName}.kraken_assignments.tsv", emit: assignments
    script:
    """
    kraken2_client \
        --port ${params.k2_port} --host-ip ${params.k2_host} \
        --report "${fastq.baseName}.kraken_report.txt" \
        --sequence ${fastq} > "${fastq.baseName}.kraken_assignments.tsv"
    """
}

process combine_kraken_outputs {
    publishDir path: "${params.out_dir}/${unique_id}/classifications", mode: 'copy'
    input:
        val unique_id
        path kraken_reports
        path kraken_assignments
    output:
        path "${params.database_set}.kraken_report.txt", emit: report
        path "${params.database_set}.kraken_assignments.tsv", emit: assignments
    script:
    if ( kraken_reports.size() == 1) {
        """
        mv ${kraken_reports[0]} "${params.database_set}.kraken_report.txt"
        mv ${kraken_assignments[0]} "${params.database_set}.kraken_assignments.tsv"
        """
    } else {
        """
        $projectDir/../bin/combine_kreports.py \\
                -r ${kraken_reports} \\
                -o "${params.database_set}.kraken_report.txt"

        cat ${kraken_assignments} > "${params.database_set}.kraken_assignments.tsv"
        """
    }
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

// this fails if the kraken file input is empty - currently have no check that it is populated
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
        //fastq = Channel.fromPath(input_fastq)
        kraken2_client(fastq)
        combine_kraken_outputs(unique_id, kraken2_client.out.report.collect(), kraken2_client.out.assignments.collect())
        bracken_length = determine_bracken_length(database)
        bracken(unique_id, combine_kraken_outputs.output.report, database, bracken_length)
        bracken_to_json(unique_id, combine_kraken_outputs.output.report, taxonomy, bracken.out.summary)
    emit:
        bracken_report = bracken.output.report
        kraken_report = combine_kraken_outputs.output.report
        kraken_assignments = combine_kraken_outputs.output.assignments
        json = bracken_to_json.out
}


