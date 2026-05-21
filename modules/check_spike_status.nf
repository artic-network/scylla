// module to calculate and filter spike ins

// NB these could be a named set of taxa, or could correspond to a reference file

process check_spike_ins {
    label "process_single"

    conda "bioconda::mappy=2.26"
    container "biocontainers/mappy:2.26--py310h83093d7_1"

    publishDir "${params.outdir}/${unique_id}/qc/", mode: 'copy', pattern: "spike*.json"
    publishDir "${params.outdir}/${unique_id}/classifications", mode: "copy", overwrite: true, pattern: "*.json"

    input:
    tuple val(unique_id), val(database_name), path(kreport), path(reads), path(spike_mapping_stats)
    val spike_ins
    path spike_in_dict
    path spike_in_ref_dir

    output:
    tuple val(unique_id), path("spike_summary.json"), emit: status, optional: true
    tuple val(unique_id), path("spike_count_summary.json"), emit: counts, optional: true
    tuple val(unique_id), path("${kreport.baseName}*.json"), emit: kreport

    script:
    """
        handle_spike_ins.py \
            -r ${kreport} \
            --spike_ins ${spike_ins} \
            --spike_in_dict ${spike_in_dict} \
            --spike_in_ref_dir ${spike_in_ref_dir} \
            --idxstats ${spike_mapping_stats} \
            --save_json
        """
}

workflow check_spike_status {
    take:
    kreport_ch
    fastq_ch
    spike_mapping_stats_ch
    
    main:
    spike_ins = "${params.spike_ins}"
    println(spike_ins)
    spike_in_dict = file("${params.spike_in_dict}", type: "file", checkIfExists: true)
    spike_in_ref_dir = file("${params.spike_in_ref_dir}", type: "dir", checkIfExists: true)

    kreport_ch.join(fastq_ch).join(spike_mapping_stats_ch) .set { input_ch }

    check_spike_ins(input_ch, spike_ins, spike_in_dict, spike_in_ref_dir)

    empty_file = file("${baseDir}/resources/empty_file")
    kreport_ch
        .map { unique_id, database_name, kreport -> [unique_id, empty_file] }
        .concat(check_spike_ins.out.status)
        .collectFile()
        .map { f -> [f.simpleName, f] }
        .set { status_ch }

    emit:
    status_ch
}
