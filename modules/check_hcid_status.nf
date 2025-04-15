// module to check for HCID

process check_hcid {

    label "process_single"

    conda "bioconda::mappy=2.28 bioconda::pyfastx=2.1.0"
    container "community.wave.seqera.io/library/mappy_pyfastx:b4cc4b80f5e5decf"

    publishDir "${params.outdir}/${unique_id}/qc/", mode: 'copy', pattern: '*.warning.json', optional:true
    publishDir "${params.outdir}/${unique_id}/qc/", mode: 'copy', pattern: '*.reads.fq', optional:true
    publishDir "${params.outdir}/${unique_id}/qc/", mode: 'copy', pattern: 'hcid.counts.csv'


    input:
        tuple val(unique_id), val(database_name), path(kreport), path(reads)
        path taxonomy
        path hcid_defs
        path hcid_refs
    output:
        tuple val(unique_id), path("*warning.json"), emit: warnings
        tuple val(unique_id), path("*reads.fq"), emit: reads
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

        PATTERN=(*.warning.json)
        if [ ! -f \${PATTERN[0]} ]; then
            echo "Found no warning files"
            touch no_warning.json
            touch no_reads.fq
        fi
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
        check_hcid.out.warnings.set{ warning_ch }
    emit:
        warning_ch
}
