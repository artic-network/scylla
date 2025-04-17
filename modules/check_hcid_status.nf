// module to check for HCID

process minimap2_hcid {

    label "process_low"

    conda "bioconda::minimap2=2.28"
    container "community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797"

    input:
        tuple val(unique_id), val(database_name), path(kreport), path(reads)
        path hcid_refs
    output:
        tuple val(unique_id), val(database_name), path(kreport), path(reads), path("hcid.mmp.sam")
    script:
        preset = ""
        if ( params.read_type == "illumina") {
            preset = "sr"
        } else if ( params.paired ) {
            preset = "sr"
        } else {
            preset = "map-ont"
        }
        """
        minimap2 -ax ${preset} ${hcid_refs} ${reads} --secondary=no -N 1 -t ${task.cpus} --sam-hit-only > hcid.mmp.sam
        """

}

process check_hcid {

    label "process_single"

    conda "bioconda::simplesam=0.1.4.1 bioconda::pyfastx=2.2.0"
    container "community.wave.seqera.io/library/pyfastx_simplesam:9161c822eef64e5a"

    publishDir "${params.outdir}/${unique_id}/qc/", mode: 'copy'

    input:
    tuple val(unique_id), val(database_name), path(kreport), path(reads), path(ref_sam)
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
            -r ${reads} \
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

    hcid_defs = file("${projectDir}/resources/hcid.json")
    hcid_refs = file("${projectDir}/resources/hcid_refs.fa.gz")

    kreport_ch.join(fastq_ch).set { input_ch }
    minimap2_hcid(input_ch,hcid_refs)
    check_hcid(minimap2_hcid.out, taxonomy, hcid_defs, hcid_refs)
    check_hcid.out.warnings.set { warning_ch }

    emit:
    warning_ch
}
