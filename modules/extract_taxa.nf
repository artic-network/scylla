// module to extract reads and de novo assemble top taxa


// it is possible that no files would be extracted if there were no subsets of reads which matched the criteria
// also note that the reads extracted don't match up with bracken abundance reestimates, although we do use those
// as more accurate numbers when deciding what to pull out (bracken doesn't provide read break down)
// probably want to count up how many have been found here for run log
// ALSO this step will currently "fail" with exitcode 2 if the number of human reads found exceeds the number specified
// in config so could be good dehuman sanity check

process split_kreport {

    label 'process_single'

    conda 'bioconda::biopython=1.78'
    container "biocontainers/pyfastx:2.0.1--py39h3d4b85c_0"

    input:
        tuple val(unique_id), path(kreport)
    output:
        tuple val(unique_id), path("*.kreport_split.txt")
    script:
        """
        split_kraken_report.py \
            -r ${kreport} \
            --splits ${params.kreport_splits}
        """
}

process extract_paired_reads {
    
    label 'process_single'
    label 'process_high_memory'

    errorStrategy {task.exitStatus in 2..3 ? 'ignore' : 'terminate'}

    publishDir path: "${params.outdir}/${unique_id}/reads_by_taxa", pattern: "reads_summary.json", mode: 'copy'

    conda 'bioconda::biopython=1.78 bioconda::tabix=1.11'
    container "biocontainers/pyfastx:2.0.1--py39h3d4b85c_0"

    input:
        tuple val(unique_id), path(fastq1), path(fastq2), path(kraken_assignments), path(kreport), val(min_reads), val(min_percent)
        path taxonomy_dir
    output:
        path("reads.*.f*"), emit: reads
        path "reads_summary.json", emit: summary
    script:
        """
        extract_kraken_reads.py \
            -s1 ${fastq1} \
            -s2 ${fastq2} \
            -k ${kraken_assignments} \
            -r ${kreport} \
            -t ${taxonomy_dir} \
            -p reads \
            --include_children \
            --max_human ${params.max_human_reads_before_rejection} \
            --min_count_descendants ${min_reads} \
            --rank ${params.extract_rank} \
            --min_percent ${min_percent}

        PATTERN=(reads.*.f*)
        if [ ! -f \${PATTERN[0]} ]; then
            echo "Found no output files - maybe there weren't any for this sample"
            exit 3
        fi
        """
}

process extract_reads {

    label 'process_single'
    label 'process_high_memory'
    
    errorStrategy {task.exitStatus in 2..3 ? 'ignore' : 'terminate'}

    publishDir path: "${params.outdir}/${unique_id}/reads_by_taxa", pattern: "reads_summary.json", mode: 'copy'

    conda 'bioconda::biopython=1.78 bioconda::tabix=1.11'
    container "biocontainers/pyfastx:2.0.1--py39h3d4b85c_0"

    input:
        tuple val(unique_id), path(fastq), path(kraken_assignments), path(kreport), val(min_reads), val(min_percent)
        path taxonomy_dir
    output:
        path("reads.*.f*"), emit: reads
        path "reads_summary.json", emit: summary
    script:
        """
        extract_kraken_reads.py \
            -s ${fastq} \
            -k ${kraken_assignments} \
            -r ${kreport} \
            -t ${taxonomy_dir} \
            -p reads \
            --include_children \
            --max_human ${params.max_human_reads_before_rejection} \
            --min_count_descendants ${min_reads} \
            --rank ${params.extract_rank} \
            --min_percent ${min_percent}

        PATTERN=(reads.*.f*)
        if [ ! -f \${PATTERN[0]} ]; then
            echo "Found no output files - maybe there weren't any for this sample"
            exit 3
        fi
        """
}

process bgzip_extracted_taxa {
      
      label 'process_medium'
  
      publishDir path: "${params.outdir}/${params.unique_id}/reads_by_taxa", mode: 'copy'
  
      conda 'bioconda::biopython=1.78 bioconda::tabix=1.11'
      container "${params.wf.container}:${params.wf.container_version}"
  
      input:
          path(read_files)
      output:
          path "reads.*.f*.gz"
      script:
          """
          for f in \$(ls reads.*.f*)
            do
            bgzip --threads $task.cpus \$f
            done
          """
}

// process check_reads {

//     label 'process_low'

//     errorStrategy {task.exitStatus == 2 ? 'ignore' : 'terminate'}

//     publishDir path: "${params.outdir}/${unique_id}/reads_by_taxa", mode: 'copy'

//     conda 'bioconda::biopython=1.78 bioconda::tabix=1.11'
//     container "${params.wf.container}:${params.wf.container_version}"

//     input:
//         tuple val(unique_id), path(reads)
//     output:
//         tuple val(unique_id), path("reads*.f*.gz")
//     script:
//         """
//         PATTERN=(reads.*.f*)
//         if [ -f \${PATTERN[0]} ]; then
//             for f in \$(ls \${PATTERN[0]})
//               do
//                 if [ ! -s \$f ]
//                   then
//                     rm \$f
//                   else
//                     bgzip --threads $task.cpus \$f
//                 fi
//               done
//         fi

//         PATTERN=(reads.*.f*.gz)
//         if [ ! -f \${PATTERN[0]} ]; then
//             echo "Found no output files - maybe there weren't any for this sample"
//             exit 2
//         fi
//         """
// }


workflow extract_taxa {
    take:
        unique_id
        fastq_ch
        assignments_ch
        kreport_ch
        taxonomy_dir
    main:
        thresholds = params.extract_thresholds
        split_kreport(kreport_ch)
        split_kreport.out.transpose()
            .map { unique_id, kreport -> [unique_id, kreport, kreport.simpleName, thresholds.containsKey(kreport.simpleName)] }
            .branch { unique_id, kreport, key, status ->
                valid: status
                    return tuple( unique_id, kreport, key )
                invalid: !status
                    return tuple( unique_id, kreport, "default" )
                }
            .set { result }
        result.valid.concat(result.invalid)
                    .map { unique_id, kreport, key -> [unique_id, kreport, thresholds.get(key,"false").get("min_reads","false"), thresholds.get(key,"false").get("min_percent","false")] }
                    .set{ kreport_params_ch }

        fastq_ch.combine(assignments_ch, by: 0)
                .combine(kreport_params_ch, by: 0)
                .set{ extract_ch }

        if ( params.paired ){
            extract_paired_reads(extract_ch, taxonomy_dir)
            extract_paired_reads.out.reads
                .set {extracted_taxa}
        } else {
            extract_reads(extract_ch, taxonomy_dir)
            extract_reads.out.reads
                .set {extracted_taxa}            
        }
        bgzip_extracted_taxa(extracted_taxa)

}


workflow {
    unique_id = "${params.unique_id}"
    fastq = file(params.fastq, type: "file", checkIfExists:true)
    assignments = file(params.kraken_assignments, type: "file", checkIfExists:true)
    kreport = file(params.kraken_report, type: "file", checkIfExists:true)
    if (unique_id == "null") {
       unique_id = "${fastq.simpleName}"
    }

    fastq_ch = Channel.of([unique_id, fastq])
    assignments_ch = Channel.of([unique_id, assignments])
    kreport_ch = Channel.of([unique_id, kreport])
    taxonomy_dir = file(params.taxonomy, type: "dir", checkIfExists:true)

    extract_taxa(unique_id, fastq_ch, assignments_ch, kreport_ch, taxonomy_dir)
}

