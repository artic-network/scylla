// module to extract reads and de novo assemble top taxa


// it is possible that no files would be extracted if there were no subsets of reads which matched the criteria
// also note that the reads extracted don't match up with bracken abundance reestimates, although we do use those
// as more accurate numbers when deciding what to pull out (bracken doesn't provide read break down)
// probably want to count up how many have been found here for run log
// ALSO this step will currently "fail" with exitcode 2 if the number of human reads found exceeds the number specified
// in config so could be good dehuman sanity check

process split_kreport {

    label 'process_single'

    conda 'bioconda::biopython=1.78
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        tuple val(unique_id), path(kreport)
    output:
        tuple val(unique_id), path("*.kreport_splits.txt")
    script:
        """
        split_kraken_report.py \
            -r ${kreport} \
            --splits ${params.kreport_splits}
        """
}

process extract_paired_reads {
    
    label 'process_high'

    errorStrategy {task.exitStatus == 2 ? 'ignore' : 'terminate'}

    publishDir path: "${params.outdir}/${unique_id}/reads_by_taxa", pattern: "reads_summary.json", mode: 'copy'

    conda 'bioconda::biopython=1.78 bioconda::tabix=1.11'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        tuple val(unique_id), path(fastq1), path(fastq2), path(kraken_assignments), path(kreport), val(min_reads), val(min_percent)
        path taxonomy_dir
    output:
        tuple val(unique_id), path("reads.*.f*"), emit: reads
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
        """
}

process extract_reads {

    label 'process_high'

    errorStrategy {task.exitStatus == 2 ? 'ignore' : 'terminate'}

    publishDir path: "${params.outdir}/${unique_id}/reads_by_taxa", pattern: "reads_summary.json", mode: 'copy'

    conda 'bioconda::biopython=1.78 bioconda::tabix=1.11'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        tuple val(unique_id), path(fastq), path(kraken_assignments), path(kreport), val(min_reads), val(min_percent)
        path taxonomy_dir
    output:
        tuple val(unique_id), path("reads.*.f*"), emit: reads
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
        """
}


process check_reads {

    label 'process_low'

    errorStrategy {task.exitStatus == 2 ? 'ignore' : 'terminate'}

    publishDir path: "${params.outdir}/${unique_id}/reads_by_taxa", mode: 'copy'

    conda 'bioconda::biopython=1.78 bioconda::tabix=1.11'
    container "${params.wf.container}@${params.wf.container_sha}"

    input:
        tuple val(unique_id), path(reads)
    output:
        tuple val(unique_id), path("reads_*.f*.gz")
    script:
        """
        PATTERN=(reads_*.f*)
        if [ -f \${PATTERN[0]} ]; then
            for f in \$(ls reads_*.f*)
              do
                if [ ! -s \$f ]
                  then
                    rm \$f
                  else
                    bgzip --threads $task.cpus \$f
                fi
              done
        fi

        PATTERN=(reads_*.f*.gz)
        if [ ! -f \${PATTERN[0]} ]; then
            echo "Found no output files - maybe there weren't any for this sample"
            exit 2
        fi
        """
}




workflow extract_taxa {
    take:
        fastq_ch
        assignments_ch
        kreport_ch
        taxonomy_dir
    main:
        split_kreport(kreport_ch)
        split_kreport.out.transpose()
                        .branch {
                            x: params.extract_thresholds.containsKey(it[1].simpleName)
                                return [it[0], it[1], params.extract_thresholds.get(it.simpleName, false).get("min_reads"),params.extract_thresholds.get(it.simpleName, false).get("min_percent")]
                            y:
                                return [it[0], it[1], params.extract_thresholds.get("default", false).get("min_reads"),params.extract_thresholds.get("default", false).get("min_percent")]
                            }
                            .set(split_kreport_ch)
        fastq_ch.join(assignments_ch)
                .join(split_kreport_ch)
                .set{ extract_ch }
        if (params.paired) {
            reads_ch = extract_paired_reads(ch_extract, taxonomy_dir)
        } else {
            reads_ch = extract_reads(ch_extract, taxonomy_dir)
        }
        check_reads(reads_ch.out.reads)
    emit:
        reads = check_reads.out
        summary = reads_ch.out.summary
