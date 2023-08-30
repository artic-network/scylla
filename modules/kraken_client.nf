include { qc_checks } from '../modules/qc_checks'
include { generate_report } from '../modules/generate_report'

kraken_compute = params.kraken_clients == 1 ? 1 : params.kraken_clients - 1

process kraken2_client {
    
    label 'process_low'
    label 'error_retry'
    maxForks kraken_compute

    conda "epi2melabs::kraken2-server=0.1.3"
    container "${params.wf.container}@${params.wf.container_sha}"
    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}

    publishDir path: "${params.outdir}/${unique_id}/classifications", mode: 'copy'

    input:
        val unique_id
        path fastq
    output:
        path "${params.database_set}.kraken_report.txt", emit: report
        path "${params.database_set}.kraken_assignments.tsv", emit: assignments
    script:
    """
    kraken2_client \
        --port ${params.k2_port} --host-ip ${params.k2_host} \
        --report "${params.database_set}.kraken_report.txt" \
        --sequence ${fastq} > "${params.database_set}.kraken_assignments.tsv"
    """
}


process determine_bracken_length {
    label "process_low"

    conda "anaconda::sed=4.8"
    container "${params.wf.container}@${params.wf.container_sha}"

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
    container "${params.wf.container}@${params.wf.container_sha}"

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
    container "${params.wf.container}@${params.wf.container_sha}"


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

process kraken_to_json {

    label "process_low"

    publishDir path: "${params.outdir}/${unique_id}/classifications", mode: 'copy'

    conda "bioconda::biopython=1.78 anaconda::Mako=1.2.3"
    container "${params.wf.container}@${params.wf.container_sha}"


    input:
        val unique_id
        path kraken_report
        path taxonomy_dir
    output:
        path "${params.database_set}.kraken.json"

    """
    awk '{ print \$5 "\t" \$3 }' "${kraken_report}" | tail -n+3 > taxacounts.txt
    cat "${kraken_report}" | cut -f5 | tail -n+3 > taxa.txt
    taxonkit lineage --data-dir ${taxonomy_dir}  -R taxa.txt  > lineages.txt
    aggregate_lineages_bracken.py \\
            -i "lineages.txt" -b "taxacounts.txt" \\
            -u "${kraken_report}" \\
            -p "temp_kraken"
    file1=`cat *.json`
    echo "{"'"${params.database_set}"'": "\$file1"}" >> "${params.database_set}.kraken.json"
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
        kraken2_client(unique_id, fastq)
        if (params.run_bracken) {
            bracken_length = determine_bracken_length(database)
            bracken(unique_id, kraken2_client.out.report, database, bracken_length)
            bracken_to_json(unique_id, kraken2_client.out.report, taxonomy, bracken.out.summary)
            out_json = bracken_to_json.out
            out_report = bracken.output.report
        } else {
            kraken_to_json(unique_id, kraken2_client.out.report, taxonomy)
            out_json = kraken_to_json.out
            out_report = kraken2_client.out.report
        }
    emit:
        bracken_report = out_report
        kraken_report = kraken2_client.out.report
        kraken_assignments = kraken2_client.out.assignments
        json = out_json
}


