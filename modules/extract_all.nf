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

    publishDir "${params.outdir}/${unique_id}/classifications", mode: params.publish_dir_mode, overwrite: false, pattern: "*.json"

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

    label "process_low"
    label "process_more_memory"

    errorStrategy { task.exitStatus in 2..3 ? "ignore" : "retry" }
    maxRetries 3

    publishDir "${params.outdir}/${unique_id}/reads_by_taxa", mode: params.publish_dir_mode

    conda "bioconda::pyfastx=2.3.1 conda-forge::numpy=2.5.2 bioconda::htslib=1.24 conda-forge::crabz"
    container "community.wave.seqera.io/library/htslib_pyfastx_numpy_which_pruned:a5bd0b38d76bbba6"

    input:
    tuple val(unique_id), path(fastq1), path(fastq2), val(database_name), path(kraken_assignments), path(kreports), val(report_config_json)
    path taxonomy_dir

    output:
    tuple val(unique_id), path("*.f*q.gz"), emit: reads
    tuple val(unique_id), path("*_summary.json"), emit: summary

    script:
    extra = ""
    if (params.reject_human) {
        extra += " --max_human ${params.max_human_reads_before_rejection}"
    }
    """
        cat <<'REPORT_CONFIG_EOF' > report_config.json
${report_config_json}
REPORT_CONFIG_EOF

        extract_taxa_from_reads.py \
            -s1 ${fastq1} \
            -s2 ${fastq2} \
            -k ${kraken_assignments} \
            -t ${taxonomy_dir} \
            -p ${unique_id} \
            --report_config report_config.json \
            --include_children \
            ${extra}

        PATTERN=(*.f*q)
        if [ ! -f \${PATTERN[0]} ]; then
            echo "Found no output files - maybe there weren't any for this sample"
            exit 3
        fi

        for f in \$(ls *.f*q); do
            crabz -p ${task.cpus} -f gzip -I \$f
        done
        """
}

process extract_taxa_reads {

    label "process_medium"
    label "process_more_memory"

    errorStrategy { task.exitStatus in 2..3 ? "ignore" : "retry" }
    maxRetries 3

    publishDir "${params.outdir}/${unique_id}/reads_by_taxa", mode: params.publish_dir_mode

    conda "bioconda::pyfastx=2.3.1 conda-forge::numpy=2.5.2 bioconda::htslib=1.24 conda-forge::crabz"
    container "community.wave.seqera.io/library/htslib_pyfastx_numpy_which_pruned:a5bd0b38d76bbba6"

    input:
    tuple val(unique_id), path(fastq), val(database_name), path(kraken_assignments), path(kreports), val(report_config_json)
    path taxonomy_dir

    output:
    tuple val(unique_id), path("*.f*q.gz"), emit: reads
    tuple val(unique_id), path("*_summary.json"), emit: summary

    script:
    extra = ""
    if (params.reject_human) {
        extra += " --max_human ${params.max_human_reads_before_rejection}"
    }
    """
        cat <<'REPORT_CONFIG_EOF' > report_config.json
${report_config_json}
REPORT_CONFIG_EOF

        extract_taxa_from_reads.py \
            -s ${fastq} \
            -k ${kraken_assignments} \
            -t ${taxonomy_dir} \
            -p ${unique_id} \
            --report_config report_config.json \
            --include_children \
            ${extra}

        PATTERN=(*.f*q)
        if [ ! -f \${PATTERN[0]} ]; then
            echo "Found no output files - maybe there weren't any for this sample"
            exit 3
        fi

        for f in \$(ls *.f*q); do
            crabz -p ${task.cpus} -f gzip -I \$f
        done
        """
}

process extract_fractions_paired_reads {

    label "process_low"
    label "process_more_memory"

    errorStrategy { task.exitStatus in 2..3 ? "ignore" : "retry" }
    maxRetries 3

    publishDir "${params.outdir}/${unique_id}/read_fractions", mode: params.publish_dir_mode

    conda "bioconda::pyfastx=2.3.1 conda-forge::numpy=2.5.2 bioconda::htslib=1.24 conda-forge::crabz"
    container "community.wave.seqera.io/library/htslib_pyfastx_numpy_which_pruned:a5bd0b38d76bbba6"

    input:
    tuple val(unique_id), path(fastq1), path(fastq2), val(database_name), path(kraken_assignments), path(kreport), val(fraction_config_json)
    path taxonomy_dir

    output:
    tuple val(unique_id), path("*.f*q.gz"), emit: reads
    tuple val(unique_id), path("*_summary.json"), emit: summary
    tuple val(unique_id), path("viral*_and_unclassified*.f*q.gz"), emit: virus, optional: true

    script:
    """
        cat <<'FRACTION_CONFIG_EOF' > fraction_config.json
${fraction_config_json}
FRACTION_CONFIG_EOF

        extract_fraction_from_reads.py \
            -s1 ${fastq1} \
            -s2 ${fastq2} \
            -k ${kraken_assignments} \
            -t ${taxonomy_dir} \
            --fraction_config fraction_config.json

        for f in \$(ls *.f*q); do
            crabz -p ${task.cpus} -f gzip -I \$f
        done
        """
}

process extract_fractions_reads {

    label "process_medium"
    label "process_more_memory"

    errorStrategy { task.exitStatus in 2..3 ? "ignore" : "retry" }
    maxRetries 3

    publishDir "${params.outdir}/${unique_id}/read_fractions", mode: params.publish_dir_mode

    conda "bioconda::pyfastx=2.3.1 conda-forge::numpy=2.5.2 bioconda::htslib=1.24 conda-forge::crabz"
    container "community.wave.seqera.io/library/htslib_pyfastx_numpy_which_pruned:a5bd0b38d76bbba6"

    input:
    tuple val(unique_id), path(fastq), val(database_name), path(kraken_assignments), path(kreport), val(fraction_config_json)
    path taxonomy_dir

    output:
    tuple val(unique_id), path("*.f*q.gz"), emit: reads
    tuple val(unique_id), path("*_summary.json"), emit: summary
    tuple val(unique_id), path("viral*_and_unclassified*.f*q.gz"), emit: virus, optional: true

    script:
    """
        cat <<'FRACTION_CONFIG_EOF' > fraction_config.json
${fraction_config_json}
FRACTION_CONFIG_EOF

        extract_fraction_from_reads.py \
            -s ${fastq} \
            -k ${kraken_assignments} \
            -t ${taxonomy_dir} \
            --fraction_config fraction_config.json

        for f in \$(ls *.f*q); do
            crabz -p ${task.cpus} -f gzip -I \$f
        done
        """
}

process merge_read_summary {

    label "process_single"

    publishDir "${params.outdir}/${unique_id}/${prefix}", pattern: "reads_summary_combined.json", mode: params.publish_dir_mode

    container "${params.wf.container}:${params.wf.container_version}"

    input:
    tuple val(unique_id), path(reads_summary)
    val prefix

    output:
    path "reads_summary_combined.json"

    script:
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

    // Build one report_config JSON per sample describing every kreport split
    // for that sample, so extract_taxa_reads/extract_taxa_paired_reads can
    // extract from all of them in a single pass over the FASTQ, instead of
    // one process (and one full FASTQ read) per split.
    split_kreport.out.reports
        .map { unique_id, kreports ->
            def kreport_list = kreports instanceof List ? kreports : [kreports]
            def config = kreport_list.collect { kreport ->
                def key = thresholds.containsKey(kreport.simpleName) ? kreport.simpleName : "default"
                def t = thresholds.get(key)
                [
                    report: kreport.name,
                    rank: t.taxon_rank,
                    min_count_descendants: t.min_reads,
                    min_percent: t.min_percent,
                ]
            }
            [unique_id, kreport_list, groovy.json.JsonOutput.toJson(config)]
        }
        .set { report_config_ch }

    assignments_ch.combine(report_config_ch, by: 0).set { classify_ch }
    fastq_ch
        .combine(classify_ch, by: 0)
        .set { extract_ch }

    if (params.paired) {
        extract_taxa_paired_reads(extract_ch, taxonomy_dir)
        extract_taxa_paired_reads.out.summary
            .groupTuple()
            .set { reads_summary_ch }
    }
    else {
        extract_taxa_reads(extract_ch, taxonomy_dir)
        extract_taxa_reads.out.summary
            .groupTuple()
            .set { reads_summary_ch }
    }
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
    // Every fraction here is defined the same way for every sample (unlike the per-kreport-split
    // report_config in extract_taxa), so the config can be built once rather than per-sample.
    fraction_config = groovy.json.JsonOutput.toJson(
        [
            [
                prefix: "virus_and_unclassified",
                taxid: ["10239", "0"],
                exclude: false,
                include_unclassified: true,
            ],
            [prefix: "virus", taxid: ["10239"], exclude: false, include_unclassified: false],
            [
                prefix: "human_filtered",
                taxid: [params.taxid_human.toString()],
                exclude: true,
                include_unclassified: false,
            ],
        ]
    )

    assignments_ch.combine(kreport_ch, by: [0, 1]).set { classify_ch }
    fastq_ch
        .combine(classify_ch, by: 0)
        .map { it + [fraction_config] }
        .set { full_extract_ch }

    if (params.paired) {
        extract_fractions_paired_reads(full_extract_ch, taxonomy_dir)
        extract_fractions_paired_reads.out.virus.set { virus_ch }
        // Unlike extract_taxa (which may group several kreport-split runs together), each
        // sample here already gets exactly one extract_fractions_*_reads call, whose single
        // "*_summary.json" output already resolves to all 3 fractions' summary files at once -
        // no groupTuple() needed (and grouping a channel with one row per key here would just
        // wrap the already-flat file list in another list).
        extract_fractions_paired_reads.out.summary.set { fractions_summary_ch }
    }
    else {
        extract_fractions_reads(full_extract_ch, taxonomy_dir)
        extract_fractions_reads.out.virus.set { virus_ch }
        extract_fractions_reads.out.summary.set { fractions_summary_ch }
    }
    merge_read_summary(fractions_summary_ch, "read_fractions")

    emit:
    virus = virus_ch
}

workflow extract_virus_fraction {
    take:
    fastq_ch
    assignments_ch
    kreport_ch
    taxonomy_dir

    main:
    fraction_config = groovy.json.JsonOutput.toJson(
        [
            [
                prefix: "virus_and_unclassified",
                taxid: ["10239", "0"],
                exclude: false,
                include_unclassified: true,
            ]
        ]
    )

    assignments_ch.combine(kreport_ch, by: [0, 1]).set { classify_ch }
    fastq_ch
        .combine(classify_ch, by: 0)
        .map { it + [fraction_config] }
        .set { full_extract_ch }

    if (params.paired) {
        extract_fractions_paired_reads(full_extract_ch, taxonomy_dir)
        extract_fractions_paired_reads.out.virus.set { virus_ch }
    }
    else {
        extract_fractions_reads(full_extract_ch, taxonomy_dir)
        extract_fractions_reads.out.virus.set { virus_ch }
    }

    emit:
    virus = virus_ch
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
    virus       = extract_fractions.out.virus
}
