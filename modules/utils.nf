// Module for utility functions

import groovy.json.JsonBuilder

include { paired_concatenate    } from '../modules/preprocess'

process get_versions {

    conda 'environment.yml'
    container "${params.wf.container}:${params.wf.container_version}"
    publishDir "${params.tracedir}", mode: 'copy'
    cpus 1

    input:
    val unique_id

    output:
    path "*.txt"

    script:
    """
    conda list > "versions_${unique_id}.txt"
    echo "${workflow.manifest.version}" > workflow_version_${unique_id}.txt
    """
}


process get_params {
    container "${params.wf.container}:${params.wf.container_version}"

    publishDir "${params.tracedir}", mode: 'copy'
    cpus 1

    input:
    val unique_id

    output:
    path "params_${unique_id}.log"

    script:
    def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '${paramsJSON}' > "params_${unique_id}.log"
    """
}

workflow get_params_and_versions {
    take:
    unique_id

    main:
    software_versions = get_versions(unique_id)
    workflow_params = get_params(unique_id)

    emit:
    software_versions
    workflow_params
}


workflow get_fastq_ch {
    take:
    unique_id

    main:
    if (params.run_dir) {
        run_dir = file("${params.run_dir}", type: "dir", checkIfExists: true)
        if (params.paired) {
            fastq_ch = Channel.fromFilePairs("${run_dir}/*_R{1,2}*.f*q*", type: "file", checkIfExists: true)
        }
        else {
            fastq_ch = Channel.fromPath("${run_dir}/*", type: "dir", checkIfExists: true, maxDepth: 1).map { [it.baseName, get_fq_files_in_dir(it)] }
        }
    }
    else if (params.paired) {
        fastq1 = file(params.fastq1, type: "file", checkIfExists: true)
        fastq2 = file(params.fastq2, type: "file", checkIfExists: true)
        fastq_ch = Channel.from([[unique_id, fastq1, fastq2]])
    }
    else if (params.fastq) {
        fastq = file(params.fastq, type: "file", checkIfExists: true)
        fastq_ch = Channel.from([[unique_id, fastq]])
    }
    else if (params.fastq_dir) {
        fastqdir = file("${params.fastq_dir}", type: "dir", checkIfExists: true)
        Channel.fromPath(fastqdir / "*.f*q*", type: "file")
            .set { fastq_ch }
        fastq_ch = fastq_ch.map { fastq -> [unique_id, fastq] }
    }

    emit:
    fastq_ch
}

workflow get_fastq_channels {
    take:
    unique_id

    main:
    if (params.run_dir) {
        run_dir = file("${params.run_dir}", type: "dir", checkIfExists: true)
        if (params.paired) {
            fastq_ch = Channel.fromFilePairs("${run_dir}/*_R{1,2}*.f*q*", type: "file", checkIfExists: true)
            
            paired_concatenate(fastq_ch)
            paired_concatenate.out.concatenated_fastq.set { combined_fastq_ch }
        } 
        else {
            fastq_ch = Channel.fromPath("${run_dir}/*", type: "dir", checkIfExists: true, maxDepth: 1).map { [it.baseName, get_fq_files_in_dir(it)] }
            fastq_ch.tap { combined_fastq_ch }
        }
    }
    else if (params.paired) {
        fastq1 = file(params.fastq1, type: "file", checkIfExists: true)
        fastq2 = file(params.fastq2, type: "file", checkIfExists: true)
        fastq_ch = Channel.from([[unique_id, fastq1, fastq2]])

        paired_concatenate(fastq_ch)
        paired_concatenate.out.concatenated_fastq.set { combined_fastq_ch }
    }
    else if (params.fastq) {
        fastq = file(params.fastq, type: "file", checkIfExists: true)
        fastq_ch = Channel.from([[unique_id, fastq]])

        fastq_ch.tap { combined_fastq_ch }
    }
    else if (params.fastq_dir) {
        fastqdir = file("${params.fastq_dir}", type: "dir", checkIfExists: true)
        Channel.fromPath(fastqdir / "*.f*q*", type: "file")
            .set { input_ch }
        fastq_ch = input_ch.map { fastq -> [unique_id, fastq] }

        fastq_ch
            .map { unique_id, fastq -> [unique_id + ".fq.gz", fastq] }
            .collectFile()
            .map { it -> [it.simpleName, it] }
            .set { combined_fastq_ch }
    }

    emit:
    processed_fastq = fastq_ch
    combined_fastq  = combined_fastq_ch
}