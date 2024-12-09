// module to extract reads and de novo assemble top taxa


// it is possible that no files would be extracted if there were no subsets of reads which matched the criteria
// also note that the reads extracted don't match up with bracken abundance reestimates, although we do use those
// as more accurate numbers when deciding what to pull out (bracken doesn"t provide read break down)
// probably want to count up how many have been found here for run log
// ALSO this step will currently "fail" with exitcode 2 if the number of human reads found exceeds the number specified
// in config so could be good dehuman sanity check

process split_kreport {

    label "process_single"

    conda "python=3.10"
    container "biocontainers/python:3.10"

    publishDir "${params.outdir}/${unique_id}/classifications", mode: "copy", overwrite: false, pattern: "*.json"

    input:
        tuple val(unique_id), val(database_name), path(kreport)
    output:
        tuple val(unique_id), path("*.kreport_split.txt"), emit: reports
        tuple val(unique_id), path("*.json"), emit: json
    script:
        """
        split_kraken_report.py \
            -r ${kreport} \
            --splits ${params.kreport_splits} \
            --save_json
        """
}

process extract_taxa_paired_reads {
    
    label "process_single"
    label "process_more_memory"

    errorStrategy {task.exitStatus in 2..3 ? "ignore" : "retry"}
    maxRetries 3

    conda "bioconda::pyfastx=2.01"
    container "biocontainers/pyfastx:2.0.1--py39h3d4b85c_0"

    input:
        tuple val(unique_id), path(fastq1), path(fastq2), val(database_name), path(kraken_assignments), path(kreport), val(taxon_rank), val(min_reads), val(min_percent)
        path taxonomy_dir
    output:
        tuple val(unique_id), path("*.fastq"), emit: reads
        tuple val(unique_id), path("${kreport}_summary.json"), emit: summary
    script:
        extra = ""
        if ( params.reject_human )
            extra += " --max_human ${params.max_human_reads_before_rejection}"
        """
        extract_taxa_from_reads.py \
            -s1 ${fastq1} \
            -s2 ${fastq2} \
            -k ${kraken_assignments} \
            -r ${kreport} \
            -t ${taxonomy_dir} \
            -p ${kreport} \
            --include_children \
            --min_count_descendants ${min_reads} \
            --rank ${taxon_rank} \
            --min_percent ${min_percent} \
            ${extra}

        PATTERN=(*.f*q)
        if [ ! -f \${PATTERN[0]} ]; then
            echo "Found no output files - maybe there weren't any for this sample"
            exit 3
        fi
        """
}

process extract_taxa_reads {

    label "process_single"
    label "process_more_memory"
    
    errorStrategy {task.exitStatus in 2..3 ? "ignore" : "retry"}
    maxRetries 3

    conda "bioconda::pyfastx=2.01"
    container "biocontainers/pyfastx:2.0.1--py39h3d4b85c_0"

    input:
        tuple val(unique_id), path(fastq), val(database_name), path(kraken_assignments), path(kreport), val(taxon_rank), val(min_reads), val(min_percent)
        path taxonomy_dir
    output:
        tuple val(unique_id), path("*.f*q"), emit: reads
        tuple val(unique_id), path("${kreport}_summary.json"), emit: summary
    script:
        extra = ""
        if ( params.reject_human )
            extra += " --max_human ${params.max_human_reads_before_rejection}"
        """
        extract_taxa_from_reads.py \
            -s ${fastq} \
            -k ${kraken_assignments} \
            -r ${kreport} \
            -t ${taxonomy_dir} \
            -p ${kreport} \
            --include_children \
            --min_count_descendants ${min_reads} \
            --rank ${taxon_rank} \
            --min_percent ${min_percent} \
            ${extra}

        PATTERN=(*.f*q)
        if [ ! -f \${PATTERN[0]} ]; then
            echo "Found no output files - maybe there weren't any for this sample"
            exit 3
        fi
        """
}

process extract_paired_virus_and_unclassified {

    label "process_single"
    label "process_more_memory"

    errorStrategy {task.exitStatus in 2..3 ? "ignore" : "retry"}
    maxRetries 3

    conda "bioconda::pyfastx=2.01"
    container "biocontainers/pyfastx:2.0.1--py39h3d4b85c_0"

    input:
        tuple val(unique_id), path(fastq1), path(fastq2), val(database_name), path(kraken_assignments), path(kreport)
        path taxonomy_dir
    output:
        tuple val(unique_id), path("*.fastq"), emit: reads
        tuple val(unique_id), path("*_summary.json"), emit: summary
    script:
        """
        extract_fraction_from_reads.py \
            -s1 ${fastq1} \
            -s2 ${fastq2} \
            -k ${kraken_assignments} \
            -t ${taxonomy_dir} \
            -p "virus_and_unclassified" \
            --taxid 10239 0 \
            --include_unclassified
        """
}

process extract_virus_and_unclassified {

    label "process_single"
    label "process_more_memory"

    errorStrategy {task.exitStatus in 2..3 ? "ignore" : "retry"}
    maxRetries 3

    conda "bioconda::pyfastx=2.01"
    container "biocontainers/pyfastx:2.0.1--py39h3d4b85c_0"

    input:
        tuple val(unique_id), path(fastq), val(database_name), path(kraken_assignments), path(kreport)
        path taxonomy_dir
    output:
        tuple val(unique_id), path("*.fastq"), emit: reads
        tuple val(unique_id), path("*_summary.json"), emit: summary
    script:
        """
        extract_fraction_from_reads.py \
            -s ${fastq} \
            -k ${kraken_assignments} \
            -t ${taxonomy_dir} \
            -p "virus_and_unclassified" \
            --taxid 10239 0 \
            --include_unclassified
        """
}


process extract_paired_virus {

    label 'process_single'
    label 'process_more_memory'

    errorStrategy {task.exitStatus in 2..3 ? 'ignore' : 'retry'}
    maxRetries 3

    conda 'bioconda::biopython=1.78 bioconda::tabix=1.11'
    container "biocontainers/pyfastx:2.0.1--py39h3d4b85c_0"

    input:
        tuple val(unique_id), path(fastq1), path(fastq2), val(database_name), path(kraken_assignments), path(kreport)
        path taxonomy_dir
    output:
        tuple val(unique_id), path("*.fastq"), emit: reads
        tuple val(unique_id), path("*_summary.json"), emit: summary
    script:
        """
        extract_fraction_from_reads.py \
            -s1 ${fastq1} \
            -s2 ${fastq2} \
            -k ${kraken_assignments} \
            -t ${taxonomy_dir} \
            -p "virus" \
            --taxid 10239
        """
}

process extract_virus {

    label 'process_single'
    label 'process_more_memory'

    errorStrategy {task.exitStatus in 2..3 ? 'ignore' : 'retry'}
    maxRetries 3

    conda 'bioconda::biopython=1.78 bioconda::tabix=1.11'
    container "biocontainers/pyfastx:2.0.1--py39h3d4b85c_0"

    input:
        tuple val(unique_id), path(fastq), val(database_name), path(kraken_assignments), path(kreport)
        path taxonomy_dir
    output:
        tuple val(unique_id), path("*.fastq"), emit: reads
        tuple val(unique_id), path("*_summary.json"), emit: summary
    script:
        """
        extract_fraction_from_reads.py \
            -s ${fastq} \
            -k ${kraken_assignments} \
            -t ${taxonomy_dir} \
            -p "virus" \
            --taxid 10239
        """
}


process extract_paired_dehumanised {

    label "process_single"
    label "process_more_memory"

    errorStrategy {task.exitStatus in 2..3 ? "ignore" : "retry"}
    maxRetries 3

    conda "bioconda::pyfastx=2.01"
    container "biocontainers/pyfastx:2.0.1--py39h3d4b85c_0"

    input:
        tuple val(unique_id), path(fastq1), path(fastq2), val(database_name), path(kraken_assignments), path(kreport)
        path taxonomy_dir
    output:
        tuple val(unique_id), path("*.fastq"), emit: reads
        tuple val(unique_id), path("*_summary.json"), emit: summary
    script:
        """
        extract_fraction_from_reads.py \
            -s1 ${fastq1} \
            -s2 ${fastq2} \
            -k ${kraken_assignments} \
            -t ${taxonomy_dir} \
            -p "human_filtered" \
            --exclude \
            --taxid ${params.taxid_human}
        """
}

process extract_dehumanised {

    label "process_single"
    label "process_more_memory"

    errorStrategy {task.exitStatus in 2..3 ? "ignore" : "retry"}
    maxRetries 3

    conda "bioconda::pyfastx=2.01"
    container "biocontainers/pyfastx:2.0.1--py39h3d4b85c_0"

    input:
        tuple val(unique_id), path(fastq), val(database_name), path(kraken_assignments), path(kreport)
        path taxonomy_dir
    output:
        tuple val(unique_id), path("*.fastq"), emit: reads
        tuple val(unique_id), path("*_summary.json"), emit: summary
    script:
        """
        extract_fraction_from_reads.py \
            -s ${fastq} \
            -k ${kraken_assignments} \
            -t ${taxonomy_dir} \
            -p "human_filtered" \
            --exclude \
            --taxid ${params.taxid_human}
        """
}

process bgzip_extracted_taxa {
      
      label "process_medium"
  
      publishDir "${params.outdir}/${unique_id}/${prefix}", mode: "copy"
  
      conda "bioconda::tabix=1.11"
      container "${params.wf.container}:${params.wf.container_version}"
  
      input:
          tuple val(unique_id), path(read_files)
          val(prefix)
      output:
          tuple val(unique_id), path("*.f*q.gz")
          tuple val(unique_id), path("viral*_and_unclassified*.f*q.gz"), emit: virus, optional:true
      script:
          """
          for f in \$(ls *.f*q)
            do
            bgzip --threads $task.cpus \$f
            done
          """
}

process merge_read_summary {

    label "process_single"

    publishDir "${params.outdir}/${unique_id}/${prefix}", pattern: "reads_summary_combined.json", mode: "copy"

    container "${params.wf.container}:${params.wf.container_version}"

    input:
        tuple val(unique_id), path(reads_summary)
        val(prefix)
    output:
        path "reads_summary_combined.json"
    
    """
    jq -s '[.[][]]' *_summary.json > reads_summary_combined.json
    """
}


workflow extract_taxa {
    take:
        fastq_ch
        assignments_ch
        kreport_ch
        taxonomy_dir
    main:
        thresholds = params.extract_thresholds
        split_kreport(kreport_ch)
        split_kreport.out.reports.transpose()
            .map { unique_id, kreport -> [unique_id, kreport, kreport.simpleName, thresholds.containsKey(kreport.simpleName)] }
            .branch { unique_id, kreport, key, status ->
                valid: status
                    return tuple( unique_id, kreport, key )
                invalid: !status
                    return tuple( unique_id, kreport, "default" )
                }
            .set { result }
        result.valid.concat(result.invalid)
                    .map { unique_id, kreport, key -> [unique_id, kreport, thresholds.get(key,"false").get("taxon_rank","false"), thresholds.get(key,"false").get("min_reads","false"), thresholds.get(key,"false").get("min_percent","false")] }
                    .set{ kreport_params_ch }

        assignments_ch.combine(kreport_params_ch, by:0).set{ classify_ch }
        fastq_ch.combine(classify_ch, by: 0)
                .set{ extract_ch }

        if ( params.paired ){
            extract_taxa_paired_reads(extract_ch, taxonomy_dir)
            extract_taxa_paired_reads.out.reads
                .set {extracted_taxa}
            extract_taxa_paired_reads.out.summary
                .groupTuple()
                .set {reads_summary_ch}
        } else {
            extract_taxa_reads(extract_ch, taxonomy_dir)
            extract_taxa_reads.out.reads
                .set {extracted_taxa}
            extract_taxa_reads.out.summary
                .groupTuple()
                .set {reads_summary_ch}       
        }
        bgzip_extracted_taxa(extracted_taxa, "reads_by_taxa")
        merge_read_summary(reads_summary_ch, "reads_by_taxa")
    emit:
        kraken_json = split_kreport.out.json
}


workflow extract_fractions {
    take:
        fastq_ch
        assignments_ch
        kreport_ch
        taxonomy_dir
    main:
         assignments_ch.combine(kreport_ch, by:[0,1]).set{ classify_ch }
         fastq_ch.combine(classify_ch, by: 0)
                 .set{ full_extract_ch }

        if ( params.paired ){
            extract_paired_dehumanised(full_extract_ch, taxonomy_dir)
            extract_paired_virus_and_unclassified(full_extract_ch, taxonomy_dir)
            extract_paired_virus(full_extract_ch, taxonomy_dir)
            extract_paired_dehumanised.out.reads
                .concat(extract_paired_virus_and_unclassified.out.reads, extract_paired_virus.out.reads)
                .set {extracted_fractions}
            extract_paired_dehumanised.out.summary
                 .concat(extract_paired_virus_and_unclassified.out.summary, extract_paired_virus.out.summary)
                 .groupTuple()
                 .set {fractions_summary_ch}
        } else {
            extract_dehumanised(full_extract_ch, taxonomy_dir)
            extract_virus_and_unclassified(full_extract_ch, taxonomy_dir)
            extract_virus(full_extract_ch, taxonomy_dir)
            extract_dehumanised.out.reads
                .concat(extract_virus_and_unclassified.out.reads, extract_virus.out.reads)
                .set {extracted_fractions}
            extract_dehumanised.out.summary
                .concat(extract_virus_and_unclassified.out.summary, extract_virus.out.summary)
                .groupTuple()
                .set {fractions_summary_ch}
        }
        bgzip_extracted_taxa(extracted_fractions, "read_fractions")
        merge_read_summary(fractions_summary_ch, "read_fractions")
    emit:
        virus = bgzip_extracted_taxa.out.virus
}

workflow extract_virus_fraction {
    take:
        fastq_ch
        assignments_ch
        kreport_ch
        taxonomy_dir
    main:
         assignments_ch.combine(kreport_ch, by:[0,1]).set{ classify_ch }
         fastq_ch.combine(classify_ch, by: 0)
                 .set{ full_extract_ch }

        if ( params.paired ){
            extract_paired_virus_and_unclassified(full_extract_ch, taxonomy_dir)
            extract_paired_virus_and_unclassified.out.reads.set{extracted_fractions}
        } else {
            extract_virus_and_unclassified(full_extract_ch, taxonomy_dir)
            extract_virus_and_unclassified.out.reads.set{extracted_fractions}
        }
        bgzip_extracted_taxa(extracted_fractions, "read_fractions")
    emit:
        virus = bgzip_extracted_taxa.out.virus
}


workflow extract_all {
    take:
        fastq_ch
        assignments_ch
        kreport_ch
        taxonomy_dir
    main:
        extract_taxa(fastq_ch, assignments_ch, kreport_ch, taxonomy_dir)
        extract_fractions(fastq_ch, assignments_ch, kreport_ch, taxonomy_dir)
    emit:
        kraken_json = extract_taxa.out.kraken_json
        virus = extract_fractions.out.virus
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
    assignments_ch = Channel.of([unique_id, "Viral", assignments])
    kreport_ch = Channel.of([unique_id, "Viral", kreport])
    taxonomy_dir = file(params.taxonomy, type: "dir", checkIfExists:true)

    extract_all(fastq_ch, assignments_ch, kreport_ch, taxonomy_dir)
}

