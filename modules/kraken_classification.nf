include { kraken2_client } from '../modules/kraken_client'
include { get_server } from '../modules/kraken_server'
include { stop_server } from '../modules/kraken_server'


process determine_bracken_length {
    label "process_low"

    conda "anaconda::sed=4.8"
    container "${params.wf.container}:${params.wf.container_version}"

    input:
        tuple val(database_name), path(database)
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

    label "process_low"

    errorStrategy {"ignore"}

    publishDir "${params.outdir}/${unique_id}/classifications", mode: "copy"

    conda "bioconda::bracken=2.7"
    container "biocontainers/bracken:2.9--py39h1f90b4d_0"

    input:
        tuple val(unique_id), val(db_name), path(kraken_report)
        tuple val(database_name), path(database)
        val bracken_length
    output:
        tuple val(unique_id), path("${database_name}.bracken_report.txt"), path("${database_name}.bracken_summary.txt"), emit: summary
        tuple val(unique_id), path("${database_name}.bracken_report.txt"), emit: report
    """
    bracken \
          -d "${database}" \
          -i "${kraken_report}" \
          -r "${bracken_length}" \
          -l "${params.bracken_level}" \
          -o "${database_name}.bracken_summary.txt" \
          -w "${database_name}.bracken_report.txt"
    """
}

process kraken_to_json {

    label "process_low"

    publishDir "${params.outdir}/${unique_id}/classifications", mode: "copy"

    conda "bioconda::taxonkit=0.15.1 python=3.10"
    container "${params.wf.container}:${params.wf.container_version}"


    input:
        tuple val(unique_id), val(database_name), path(kraken_report)
        path taxonomy_dir
    output:
       tuple val(unique_id), val(database_name), path("${database_name}.kraken.json")

    """
    cat "${kraken_report}" | cut -f5,3 | tail -n+3 > taxacounts.txt
    cat "${kraken_report}" | cut -f5 | tail -n+3 > taxa.txt
    taxonkit lineage --data-dir ${taxonomy_dir}  -R taxa.txt  > lineages.txt
    aggregate_lineages_bracken.py \\
            -i "lineages.txt" -b "taxacounts.txt" \\
            -u "${kraken_report}" \\
            -p "temp_kraken"
    file1=`cat *.json`
    echo "{"'"${database_name}"'": "\$file1"}" >> "${database_name}.kraken.json"
    """
}

workflow run_kraken_and_bracken {
    take:
        fastq_ch
        server_ch
        database_ch
    main:
        kraken2_client(server_ch, fastq_ch)
        if (params.run_bracken) {
            bracken_length = determine_bracken_length(database_ch)
            bracken(kraken2_client.out.report, database_ch, bracken_length)
        }

    emit:
        assignments = kraken2_client.out.assignments
        kreport = kraken2_client.out.report
}

workflow kraken_classify {
    take:
        fastq_ch
        database_key
        raise_server
    main:
        get_server(database_key, raise_server)

        run_kraken_and_bracken(fastq_ch, get_server.out.server, get_server.out.database)

        if (raise_server)
            stop_server(get_server.out.server, run_kraken_and_bracken.out.kreport.collect())
    emit:
        assignments = run_kraken_and_bracken.out.assignments
        kreport = run_kraken_and_bracken.out.kreport
}
