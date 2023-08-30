// module to extract reads and de novo assemble top taxa


// it is possible that no files would be extracted if there were no subsets of reads which matched the criteria
// also note that the reads extracted don't match up with bracken abundance reestimates, although we do use those
// as more accurate numbers when deciding what to pull out (bracken doesn't provide read break down)
// probably want to count up how many have been found here for run log
// ALSO this step will currently "fail" with exitcode 2 if the number of human reads found exceeds the number specified
// in config so could be good dehuman sanity check
process extract_paired_reads {
    
    label 'process_high'

    errorStrategy {task.exitStatus == 2 ? 'ignore' : 'terminate'}

    publishDir path: "${params.outdir}/${unique_id}/reads_by_taxa", mode: 'copy'

    conda 'bioconda::biopython=1.78 bioconda::tabix=1.11'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        val unique_id
        path fastq1
        path fastq2
        path kraken_assignments
        path bracken_report
        path taxonomy_dir
    output:
        path "reads.*.f*.gz", emit: reads
        path "reads_summary.json", emit: summary
    script:
        """
        extract_kraken_reads.py \
            -s1 ${fastq1} \
            -s2 ${fastq2} \
            -k ${kraken_assignments} \
            -r ${bracken_report} \
            -t ${taxonomy_dir} \
            -p reads \
            --include_children \
            --max_human ${params.max_human_reads_before_rejection} \
            --min_count_descendants ${params.assembly_min_reads} \
            --rank ${params.assembly_rank} \
            --min_percent ${params.assembly_min_percent}

        for f in \$(ls reads.*.f*)
          do
            bgzip --threads $task.cpus \$f
          done        
        """
}

process extract_reads {

    label 'process_high'

    errorStrategy {task.exitStatus == 2 ? 'ignore' : 'terminate'}

    publishDir path: "${params.outdir}/${unique_id}/reads_by_taxa", mode: 'copy'

    conda 'bioconda::biopython=1.78 bioconda::tabix=1.11'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        val unique_id
        path fastq
        path kraken_assignments
        path bracken_report
        path taxonomy_dir
    output:
        path "reads.*.f*.gz", emit: reads
        path "reads_summary.json", emit: summary
    script:
        """
        extract_kraken_reads.py \
            -s ${fastq} \
            -k ${kraken_assignments} \
            -r ${bracken_report} \
            -t ${taxonomy_dir} \
            -p reads \
            --include_children \
            --max_human ${params.max_human_reads_before_rejection} \
            --min_count_descendants ${params.assembly_min_reads} \
            --rank ${params.assembly_rank} \
            --min_percent ${params.assembly_min_percent}

        PATTERN=(reads.*.f*)
        if [ -f \${PATTERN[0]} ]; then
            for f in \$(ls reads.*.f*)
              do
                bgzip --threads $task.cpus \$f
              done
        else
            exit 2
        fi
        """
}




// workflow extract_taxa {
//     take:
//         unique_id
//         fastq_1
//         fastq_2
//         kraken_assignments
//         kraken_report
//         bracken_report
//     main:
//         extract_reads(unique_id, fastq_1, fastq_2, kraken_assignments, kraken_report, bracken_report)
// }

// workflow {
//     unique_id = "${params.unique_id}"
//     if ("${params.unique_id}" == "null") {
//         unique_id = "${fastq.simpleName}"
//     }

//     kraken_assignments = file("${params.kraken_assignments}")
//     bracken_report = file("${params.bracken_report}")
//     kraken_report = file("${params.kraken_report}")

//     if (params.paired) {
//             fastq_1 = file(params.fastq1, type: "file", checkIfExists:true)
//             fastq_2 = file(params.fastq2, type: "file", checkIfExists:true)
//             extract_paired_reads(unique_id, fastq_1, fastq_2, kraken_assignments, kraken_report, bracken_report)
//     } else {
//             fastq = file(params.fastq, type: "file", checkIfExists:true)
//             extract_reads(unique_id, fastq, kraken_assignments, kraken_report, bracken_report)
//     }
// }