// module to extract reads and de novo assemble top taxa

process extract_reads {
    publishDir path: "${params.out_dir}/${unique_id}/reads_by_taxa", mode: 'copy', pattern: "*.fa"
    input:
        val unique_id
        path fastq
        path kraken_assignments
        path kraken_report
        path bracken_report
    output:
        path("reads.${tax_id}.f*")
    script:
    """
    $projectDir/../bin/extract_kraken_reads.py \
        -s ${fastq} \
        -k ${kraken_assignments} \
        -r ${kraken_report} \
        -b ${bracken_report} \
        -p reads \
        --include_children \
        --max_human ${params.max_human_reads_before_rejection} \
        --min_count_descendants ${params.assembly_min_reads} \
        --rank ${params.assembly_rank} \
        --min_percent ${params.assembly_min_percent}
    """
}

process denovo_assemble {
    errorStrategy 'ignore'
    publishDir path: "${params.out_dir}/${unique_id}/assemblies_by_taxa/", mode: 'copy'
    input:
        val unique_id
        path fastq
    output:
        path "${fastq.baseName}"
    script:
    mode = "${params.read_type}"
    if( mode == 'ont' )
        """
            flye \
                --nano-raw ${fasta} \
                -o "${fasta.baseName}"
        """
    else if( mode == 'illumina' )
        """
            megahit \
                --r ${fasta} \
                -o "${fasta.baseName}"
        """
    else
        error "Invalid alignment read_type: ${mode} - must be one of [illumina, ont]"
}

workflow assemble_taxa {
    take:
        unique_id
        fastq
        kraken_assignments
        kraken_report
        bracken_report
    main:
        extract_reads(unique_id, fastq, kraken_assignments, kraken_report, bracken_report)
        denovo_assemble(unique_id, extract_reads_for_id.out)
}

workflow {
    fastq = file("${params.fastq}")
    if (!fastq.exists()) {
            throw new Exception("--fastq: File doesn't exist, check path.")
        }

    unique_id = "${params.unique_id}"
    if ("${params.unique_id}" == "null") {
        unique_id = "${fastq.simpleName}"
    }

    kraken_assignments = file("${params.kraken_assignments}")
    bracken_report = file("${params.bracken_report}")
    kraken_report = file("${params.kraken_report}")

    assemble_taxa(unique_id, fastq, kraken_assignments, kraken_report, bracken_report)
}