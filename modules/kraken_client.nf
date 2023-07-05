include { qc_checks } from '../modules/qc_checks'
include { generate_report } from '../modules/generate_report'

process kraken2_client {
    
    label 'process_low'
    label 'error_retry'

    conda "epi2melabs::kraken2-server=0.1.3"
    container "biowilko/scylla@${params.wf.container_sha}"

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

    label 'process_single'

    container "biowilko/scylla@${params.wf.container_sha}"

    publishDir path: "${params.outdir}/${unique_id}/classifications", mode: 'copy'
    
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
        combine_kreports.py \\
                -r ${kraken_reports} \\
                -o "${params.database_set}.kraken_report.txt"

        cat ${kraken_assignments} > "${params.database_set}.kraken_assignments.tsv"
        """
    }
}

process determine_bracken_length {
    label "process_low"

    conda "anaconda::sed=4.8"
    container "biowilko/scylla@${params.wf.container_sha}"

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
    
    label 'process_low'

    publishDir path: "${params.outdir}/${unique_id}/classifications", mode: 'copy'

    conda "bioconda::bracken=2.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bracken:2.7--py39hc16433a_0':
        'biocontainers/bracken:2.7--py39hc16433a_0' }"

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
    
    label "process_low"

    publishDir path: "${params.outdir}/${unique_id}/classifications", mode: 'copy'
    
    conda "bioconda::biopython=1.78 anaconda::Mako=1.2.3"
    container "biowilko/scylla@${params.wf.container_sha}"


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
    aggregate_lineages_bracken.py \\
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


