// module to check for HCID

process check_hcid {

    label "process_single"

    conda "bioconda::mappy=2.26"
    container "biocontainers/mappy:2.26--py310h83093d7_1"

    publishDir "${params.outdir}/${unique_id}/qc/", mode: 'copy'

    input:
        tuple val(unique_id), path(kreport), path(reads)
        path taxonomy
        path hcid_defs
        path hcid_refs
    output:
        tuple val(unique_id), path("*.warning.json"), emit: warnings, optional: true
        tuple val(unique_id), path("hcid.counts.csv"), emit: counts
    script:
        preset = ""
        if ( params.read_type == "illumina") {
            preset = "--illumina"
        } else if ( params.paired ) {
            preset = "--illumina"
        }
        """
        check_hcid.py \
            -k ${kreport} \
            -r ${reads} \
            -t ${taxonomy} \
            -i ${hcid_defs} \
            -d ${hcid_refs} \
            -p "hcid" ${preset}
        """
}

workflow check_hcid_status {
    take:
        kreport_ch
        fastq_ch
        taxonomy
    main:
        hcid_defs = file("$projectDir/resources/hcid.json")
        hcid_refs = file("$projectDir/resources/hcid_refs.fa.gz")

        kreport_ch.join(fastq_ch).set{input_ch}
        check_hcid(input_ch, taxonomy, hcid_defs, hcid_refs)

        empty_file = file("$projectDir/resources/empty_file")
        kreport_ch.map{unique_id, kreport -> [unique_id, empty_file]}
            .concat(check_hcid.out.warnings)
            .collectFile()
            .map{f -> [f.simpleName, f]}
            .set{warning_ch}
    emit:
        warning_ch
}

workflow {
    unique_id = "${params.unique_id}"
    fastq = file(params.fastq, type: "file", checkIfExists:true)
    kreport = file(params.kraken_report, type: "file", checkIfExists:true)
    if (unique_id == "null") {
       unique_id = "${fastq.simpleName}"
    }
    kreport_ch = Channel.of([unique_id, kreport])
    fastq_ch = Channel.of([unique_id, fastq])

    taxonomy_dir = file(params.taxonomy, type: "dir", checkIfExists:true)

    check_hcid_status(kreport_ch, fastq_ch, taxonomy_dir)
}

