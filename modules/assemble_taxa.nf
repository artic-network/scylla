// module to extract reads and de novo assemble top taxa

process extract_ids {
    input:
        path bracken_report
    output:
        path tax_ids
    script:
    """
    $projectDir/../bin/extract_ids.py \
        --bracken_report ${bracken_report} \
        --target_rank ${params.assembly_rank} \
        --min_num_reads ${params.assembly_min_reads} \
        --max_n ${params.assembly_max_n} \
        --min_percent ${params.assembly_min_percent} \
        > tax_ids
    """
}

process extract_reads_for_id {
    publishDir path: "${params.out_dir}/${sample_id}/reads_by_taxa", mode: 'copy', pattern: "*.fa"
    input:
        val tax_id
        val sample_id
        path fastq
        path kraken_assignments
        path bracken_report
    output:
        tuple val(sample_id), val(tax_id), path("${fastq.baseName}.${tax_id}.fa")
    script:
    """
    $projectDir/../bin/extract_kraken_reads.py \
        -s ${fastq} \
        -k ${kraken_assignments} \
        -t ${tax_id} \
        -o "${fastq.baseName}.${tax_id}.fa" \
        --include-children \
        --report ${bracken_report}
    """
}

process denovo_assemble {
    errorStrategy 'ignore'
    publishDir path: "${params.out_dir}/${sample_id}/assemblies/", mode: 'copy'
    input:
        tuple val(sample_id), val(tax_id), path(fasta)
    output:
        path "${fasta.baseName}"
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
        sample_id
        fastq
        kraken_assignments
        bracken_report
    main:
        extract_ids(bracken_report)
        tax_ids = extract_ids.out.splitText().map{it -> it.trim()}
        extract_reads_for_id(tax_ids, sample_id, fastq, kraken_assignments,bracken_report)
        denovo_assemble(extract_reads_for_id.out)
}

workflow {
    fastq = file("${params.fastq}")
    if (!fastq.exists()) {
            throw new Exception("--fastq: File doesn't exist, check path.")
        }

    sample_id = "${params.sample_id}"
    if (sample_id == "null") {
        sample_id = "${fastq.simpleName}"
    }

    kraken_assignments = file("${params.kraken_assignments}")
    bracken_report = file("${params.bracken_report}")
    assemble_taxa(sample_id, fastq, kraken_assignments, bracken_report)
}