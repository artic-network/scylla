// module to check for HCID

process rammap_hcid {

    label "process_medium"

    container "quay.io/biowilko/rammap:1.1.2"

    input:
    tuple val(unique_id), val(database_name), path(kreport), path(reads)
    path hcid_refs

    output:
    tuple val(unique_id), val(database_name), path(kreport), path("hcid.mmp.sam")

    script:
    preset = ""
    if (params.read_type == "illumina") {
        preset = "sr"
    }
    else if (params.paired) {
        preset = "sr"
    }
    else {
        preset = "map-ont"
    }
    """
    rammap -ax ${preset} --secondary=no -N 1 -t ${task.cpus} --sam-hit-only ${hcid_refs} ${reads} > hcid.mmp.sam
    """
}

process check_hcid {

    label "process_single"
    label "process_more_memory"

    conda "bioconda::simplesam=0.1.4.1 bioconda::pyfastx=2.2.0"
    container "community.wave.seqera.io/library/pyfastx_simplesam:9161c822eef64e5a"

    publishDir "${params.outdir}/${unique_id}/qc/", mode: params.publish_dir_mode

    input:
    tuple val(unique_id), val(database_name), path(kreport), path(ref_sam)
    path taxonomy
    path hcid_defs
    path hcid_refs

    output:
    tuple val(unique_id), path("*warning.json"), emit: warnings, optional: true
    tuple val(unique_id), path("*reads.fq"), emit: reads, optional: true
    tuple val(unique_id), path("hcid.counts.csv"), emit: counts

    script:
    """
    check_hcid.py \
        -k ${kreport} \
        -t ${taxonomy} \
        -i ${hcid_defs} \
        -s ${ref_sam} \
        -d ${hcid_refs} \
        -p "hcid"
    """
}

workflow check_hcid_status {
    take:
    kreport_ch
    fastq_ch
    taxonomy

    main:

    hcid_defs = file("${projectDir}/resources/hcid/hcid.json")
    hcid_refs = file("${projectDir}/resources/hcid/hcid_refs.fa.gz")

    kreport_ch.join(fastq_ch).set { input_ch }
    rammap_hcid(input_ch, hcid_refs)
    check_hcid(rammap_hcid.out, taxonomy, hcid_defs, hcid_refs)
    check_hcid.out.warnings.set { warning_ch }

    emit:
    warning_ch
}
