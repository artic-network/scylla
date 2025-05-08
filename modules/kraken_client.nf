kraken_compute = params.kraken_clients == 1 ? 1 : params.kraken_clients - 1

process kraken2_client {

    label "process_single"
    label "error_retry"
    maxForks kraken_compute

    conda "nanoporetech::kraken2-server=0.1.7"
    container "${params.wf.container}:${params.wf.container_version}"
    containerOptions { workflow.profile != "singularity" ? "--network host" : "" }

    publishDir "${params.outdir}/${unique_id}/classifications", mode: "copy"

    input:
    tuple val(database_name), val(host), val(port)
    tuple val(unique_id), path(fastq)

    output:
    tuple val(unique_id), val(database_name), path("${database_name}.kraken_assignments.tsv"), emit: assignments
    tuple val(unique_id), val(database_name), path("${database_name}.kraken_report.txt"), emit: report

    script:
    """
    kraken2_client \
        --port ${port} --host-ip ${host} \
        --report "${database_name}.kraken_report.txt" \
        --sequence ${fastq} > "${database_name}.kraken_assignments.tsv"
    """
}
