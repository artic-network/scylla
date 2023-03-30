// module to extract reads and de novo assemble top taxa


// it is possible that no files would be extracted if there were no subsets of reads which matched the criteria
// also note that the reads extracted don't match up with bracken abundance reestimates, although we do use those
// as more accurate numbers when deciding what to pull out (bracken doesn't provide read break down)
// probably want to count up how many have been found here for run log
// ALSO this step will currently "fail" with exitcode 2 if the number of human reads found exceeds the number specified
// in config so could be good dehuman sanity check
process extract_reads {
    publishDir path: "${params.out_dir}/${unique_id}/reads_by_taxa", mode: 'copy'
    input:
        val unique_id
        path fastq
        path kraken_assignments
        path kraken_report
        path bracken_report
    output:
        path "reads.*.f*.gz"
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

    for f in \$(ls reads.*.f*)
      do
        gzip \$f
      done
    """
}

// Very lossy step - we do not expect to get successful de novo assemblies, so we ignore any failures. These could
// either be as a result of the pipeline failing for one of many reasons, or because the file generated at the end
// of a successful pipeline run is empty, in which case it would be meaningless to copy
// I've not had a successful test case with the very small numbers of reads/database on my laptop
process denovo_assemble {
    errorStrategy 'ignore'
    publishDir path: "${params.out_dir}/${unique_id}/assemblies_by_taxa/", mode: 'copy'
    input:
        val unique_id
        path fastq
    output:
        path "assembly.*.f*.gz"
    script:
    mode = "${params.read_type}"
    outfile = "${fastq.name}".replaceAll("reads","assembly")
    if( mode == 'ont' )
        """
            flye \
                --nano-raw ${fastq} \
                -o "assembly"
            if [ -s assembly/final.contigs.fa ]
            then
                mv assembly/final.contigs.fa ${outfile}
                gzip ${outfile}
            fi
        """
    else if( mode == 'illumina' )
        """
            megahit \
                --r ${fastq} \
                -o "assembly"
            if [ -s assembly/final.contigs.fa ]
            then
                mv assembly/final.contigs.fa ${outfile}
                gzip ${outfile}
            fi
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
        denovo_assemble(unique_id, extract_reads.out.flatten())
}

workflow {
    fastq = file("${params.fastq}", type: "file", checkIfExists:true)

    unique_id = "${params.unique_id}"
    if ("${params.unique_id}" == "null") {
        unique_id = "${fastq.simpleName}"
    }

    kraken_assignments = file("${params.kraken_assignments}")
    bracken_report = file("${params.bracken_report}")
    kraken_report = file("${params.kraken_report}")

    assemble_taxa(unique_id, fastq, kraken_assignments, kraken_report, bracken_report)
}