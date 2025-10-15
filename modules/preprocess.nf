process prepare_spike_references {
    label "process_single"

    conda "python=3.10"
    container "biocontainers/python:3.10"

    input:
    val spike_ins
    path spike_in_dict
    path spike_in_ref_dir

    output:
    path "combined_spikes.fa", emit: combined_ref, optional: true

    script:
    """
    prep_spike_refs.py --spike_ins ${spike_ins} --spike_in_dict ${spike_in_dict} --spike_in_ref_dir ${spike_in_ref_dir} -o combined_spikes.fa
    """
}

process remove_spike_paired {
    label "process_medium"

    conda "bioconda::minimap2 bioconda::samtools"
    container "community.wave.seqera.io/library/minimap2_samtools:c2863f226d833ac8"
    publishDir "${params.outdir}/${unique_id}/classifications", mode: "copy", overwrite: true, pattern: "*spikes*.fastq"
    publishDir "${params.outdir}/${unique_id}/qc", mode: "copy", overwrite: true, pattern: "*.txt"

    input:
    tuple val(unique_id), path(fastq1), path(fastq2)
    path combined_reference

    output:
      tuple val(unique_id), path("${unique_id}_removed_1.fastq"), path("${unique_id}_removed_2.fastq"), emit: removed_paired
      tuple val(unique_id), path("${unique_id}_spikes_1.fastq"), path("${unique_id}_spikes_2.fastq"), emit: spike_reads
      tuple val(unique_id), path("read_count_removed.txt"), path("read_count_spikes.txt"), emit: stats

    script:
    """
        minimap2 -ax sr $combined_reference ${fastq1} ${fastq2} | samtools view -bS - | samtools sort -o sorted.bam

        samtools view -b -f 4 sorted.bam | samtools fastq -1 ${unique_id}_removed_1.fastq -2 ${unique_id}_removed_2.fastq -
        samtools view -b -F 4 sorted.bam | samtools fastq -1 ${unique_id}_spikes_1.fastq -2 ${unique_id}_spikes_2.fastq -

        samtools view -c -f 4 sorted.bam > read_count_removed.txt
        samtools view -c -F 4 sorted.bam > read_count_spikes.txt
    """
}

process remove_spike_single {
    label "process_medium"

    conda "bioconda::minimap2 bioconda::samtools"
    container "community.wave.seqera.io/library/minimap2_samtools:c2863f226d833ac8"

    publishDir "${params.outdir}/${unique_id}/classifications", mode: "copy", overwrite: true, pattern: "*spikes*.fastq"
    publishDir "${params.outdir}/${unique_id}/qc", mode: "copy", overwrite: true, pattern: "*.txt"

    input:
    tuple val(unique_id), path(fastq)
    path combined_reference

    output:
      tuple val(unique_id), path("${unique_id}_removed.fastq"), emit: removed_single
      tuple val(unique_id), path("${unique_id}_spikes.fastq"), emit: spike_reads
      tuple val(unique_id), path("read_count_unmapped.txt"), path("read_count_mapped.txt"), emit: stats

    script:
    """
        minimap2 -ax map-ont $combined_reference ${fastq} | samtools view -bS - | samtools sort -o sorted.bam

        samtools view -b -f 4 sorted.bam | samtools fastq - > ${unique_id}_removed.fastq
        samtools view -b -F 4 sorted.bam | samtools fastq - > ${unique_id}_spikes.fastq

        samtools view -c -f 4 sorted.bam > read_count_unmapped.txt
        samtools view -c -F 4 sorted.bam > read_count_mapped.txt
    """
}

process fastp_paired {

    label "process_medium"

    publishDir "${params.outdir}/${unique_id}/preprocess/", mode: "copy"

    container "${params.wf.container}:${params.wf.container_version}"

    errorStrategy { task.exitStatus in [255, 10] ? "ignore" : "terminate" }

    input:
    tuple val(unique_id), path(fastq_1), path(fastq_2)

    output:
    tuple val(unique_id), path("${unique_id}_1.fastp.fastq.gz"), path("${unique_id}_2.fastp.fastq.gz"), emit: fastq
    path "${unique_id}.fastp.json"

    script:
    """
    fastp \\
        --in1 ${fastq_1} \\
        --in2 ${fastq_2} \\
        --out1 ${unique_id}_1.fastp.fastq \\
        --out2 ${unique_id}_2.fastp.fastq \\
        --json ${unique_id}.fastp.json \\
        --low_complexity_filter \\
        --thread ${task.cpus} \\
        2> ${unique_id}.fastp.log

    if [ -s ${unique_id}_1.fastp.fastq ]; then
        bgzip --threads ${task.cpus} -c ${unique_id}_1.fastp.fastq > ${unique_id}_1.fastp.fastq.gz
    else
        exit 10
    fi

    if [ -s ${unique_id}_2.fastp.fastq ]; then
        bgzip --threads ${task.cpus} -c ${unique_id}_2.fastp.fastq > ${unique_id}_2.fastp.fastq.gz
    else
        exit 10
    fi    
    """
}

process fastp_single {

    label "process_medium"

    publishDir "${params.outdir}/${unique_id}/preprocess/", mode: "copy"

    container "${params.wf.container}:${params.wf.container_version}"

    errorStrategy { task.exitStatus in [255, 10] ? "ignore" : "terminate" }

    input:
    tuple val(unique_id), path(fastq)

    output:
    tuple val(unique_id), path("${unique_id}.fastp.fastq.gz"), emit: fastq
    path "${unique_id}.fastp.json"

    script:

    """
    fastp \\
        --in1 ${fastq} \\
        --out1 ${unique_id}.fastp.fastq \\
        --json ${unique_id}.fastp.json \\
        --thread ${task.cpus} \\
        --disable_adapter_trimming \\
        --low_complexity_filter \\
        --qualified_quality_phred 10 \\
        2> ${unique_id}.fastp.log

    if [ -s ${unique_id}.fastp.fastq ]; then
        bgzip --threads ${task.cpus} -c ${unique_id}.fastp.fastq > ${unique_id}.fastp.fastq.gz
    else
        exit 10
    fi
    """
}

process paired_concatenate {

    label "process_low"

    errorStrategy { task.exitStatus in [5, 8] ? 'ignore' : 'terminate' }

    publishDir "${params.outdir}/${unique_id}/preprocess/", mode: "copy"

    container "${params.wf.container}:${params.wf.container_version}"

    input:
    tuple val(unique_id), path(processed_fastq_1), path(processed_fastq_2)

    output:
    tuple val(unique_id), path("${unique_id}.concatenated.fastq.gz"), emit: concatenated_fastq

    script:
    """
    concatenate_reads.py --no-interleave \\
        ${processed_fastq_1} ${processed_fastq_2} \\
        --strict \\
        > ${unique_id}.concatenated.fastq

    
    bgzip --threads ${task.cpus} -c ${unique_id}.concatenated.fastq > ${unique_id}.concatenated.fastq.gz
    """
}

workflow preprocess {
    take:
    unique_id

    main:
    if (params.paired) {
        fastq1 = file(params.fastq1, type: "file", checkIfExists: true)
        fastq2 = file(params.fastq2, type: "file", checkIfExists: true)
        input_ch = Channel.from([[unique_id, fastq1, fastq2]])

        if (params.spike_ins) {
              spike_in_dict = file(params.spike_in_dict, type: "file", checkIfExists: true)
              spike_in_ref_dir = file(params.spike_in_ref_dir, type: "dir", checkIfExists: true)

              prepare_spike_references(params.spike_ins, spike_in_dict, spike_in_ref_dir)
              remove_spike_paired(input_ch, prepare_spike_references.out.combined_ref)
              spike_removed_ch = remove_spike_paired.out.removed_paired
          } else {
              spike_removed_ch = input_ch
          }

        fastp_paired(spike_removed_ch)
        fastp_paired.out.fastq.set { processed_fastq_ch }

        paired_concatenate(fastp_paired.out.fastq)
        paired_concatenate.out.concatenated_fastq.set { combined_fastq_ch }
    }
    else if (params.fastq) {
        fastq = file(params.fastq, type: "file", checkIfExists: true)
        input_ch = Channel.from([[unique_id, fastq]])

        if (params.spike_ins) {
              spike_in_dict = file(params.spike_in_dict, type: "file", checkIfExists: true)
              spike_in_ref_dir = file(params.spike_in_ref_dir, type: "dir", checkIfExists: true)

              prepare_spike_references(params.spike_ins, spike_in_dict, spike_in_ref_dir)
              remove_spike_single(input_ch, prepare_spike_references.out.combined_ref)
              spike_removed_ch = remove_spike_single.out.removed_single
          } else {
              spike_removed_ch = input_ch
          }

        fastp_single(spike_removed_ch)

        fastp_single.out.fastq
            .tap { processed_fastq_ch }
            .set { combined_fastq_ch }
    }
    else if (params.fastq_dir) {
        fastqdir = file("${params.fastq_dir}", type: "dir", checkIfExists: true)
        Channel.fromPath(fastqdir / "*.f*q*", type: "file")
            .set { fastq_ch }
        input_ch = fastq_ch.map { fastq -> [unique_id, fastq] }

        if (params.spike_ins) {
              spike_in_dict = file(params.spike_in_dict, type: "file", checkIfExists: true)
              spike_in_ref_dir = file(params.spike_in_ref_dir, type: "dir", checkIfExists: true)

              prepare_spike_references(params.spike_ins, spike_in_dict, spike_in_ref_dir)
              remove_spike_single(input_ch, prepare_spike_references.out.combined_ref)
              spike_removed_ch = remove_spike_single.out.removed_single
          } else {
              spike_removed_ch = input_ch
          }

        fastp_single(spike_removed_ch)
        fastp_single.out.fastq
            .map { unique_id, fastq -> [unique_id + ".fq.gz", fastq] }
            .collectFile()
            .map { it -> [it.simpleName, it] }
            .tap { processed_fastq_ch }
            .set { combined_fastq_ch }
    }

    emit:
    processed_fastq = processed_fastq_ch
    combined_fastq  = combined_fastq_ch
}
