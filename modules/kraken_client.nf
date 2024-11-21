
kraken_compute = params.kraken_clients == 1 ? 1 : params.kraken_clients - 1

process kraken2_client {
    
    label "process_single"
    label "error_retry"
    maxForks kraken_compute

    conda "nanoporetech::kraken2-server=0.1.7"
    container "${params.wf.container}:${params.wf.container_version}"
    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}

    publishDir "${params.outdir}/${unique_id}/classifications", mode: "copy"

    input:
        tuple val(unique_id), path(fastq)

    output:
        tuple val(unique_id), path("${params.database_set}.kraken_assignments.tsv"), emit: assignments
        tuple val(unique_id), path("${params.database_set}.kraken_report.txt"), emit: report

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
    container "${params.wf.container}:${params.wf.container_version}"

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
    
    label "process_low"

    errorStrategy {"ignore"}

    publishDir "${params.outdir}/${unique_id}/classifications", mode: "copy"

    conda "bioconda::bracken=2.7"
    container "biocontainers/bracken:2.9--py39h1f90b4d_0"

    input:
        tuple val(unique_id), path(kraken_report)
        path database
        val bracken_length
    output:
        tuple val(unique_id), path("${params.database_set}.bracken_report.txt"), path("${params.database_set}.bracken_summary.txt"), emit: summary
        tuple val(unique_id), path("${params.database_set}.bracken_report.txt"), emit: report
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

    publishDir "${params.outdir}/${unique_id}/classifications", mode: "copy"
    
    conda "bioconda::taxonkit=0.15.1 python=3.10"
    container "${params.wf.container}:${params.wf.container_version}"

    input:
        tuple val(unique_id), path(bracken_report), path(bracken_summary)
        path taxonomy_dir
    output:
        tuple val(unique_id), path("${params.database_set}.bracken.json")

    """
    cat "${bracken_summary}" | cut -f2,6 | tail -n+2 > taxacounts.txt
    cat "${bracken_summary}" | cut -f2 | tail -n+2 > taxa.txt
    taxonkit lineage --data-dir ${taxonomy_dir}  -R taxa.txt  > lineages.txt
    aggregate_lineages_bracken.py \\
            -i "lineages.txt" -b "taxacounts.txt" \\
            -u "${bracken_report}" \\
            -p "temp_bracken"
    file1=`cat *.json`
    echo "{"'"${params.database_set}"'": "\$file1"}" >> "${params.database_set}.bracken.json"
    """
}

process kraken_to_json {

    label "process_low"

    publishDir "${params.outdir}/${unique_id}/classifications", mode: "copy"

    conda "bioconda::taxonkit=0.15.1 python=3.10"
    container "${params.wf.container}:${params.wf.container_version}"


    input:
        tuple val(unique_id), path(kraken_report)
        path taxonomy_dir
    output:
       tuple val(unique_id), path("${params.database_set}.kraken.json")

    """
    cat "${kraken_report}" | cut -f5,3 | tail -n+3 > taxacounts.txt
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
        fastq_ch
        database
        taxonomy
    main:
        unique_id = "test"
        kraken2_client(fastq_ch)
        if (params.run_bracken) {
            bracken_length = determine_bracken_length(database)
            bracken(kraken2_client.out.report, database, bracken_length)

        }
        kraken_to_json(kraken2_client.out.report, taxonomy)
        out_json = kraken_to_json.out
        out_report = kraken2_client.out.report
    emit:
        assignments = kraken2_client.out.assignments
        kreport = out_report
        json = out_json
}

