// Module to print params and versions of software

import groovy.json.JsonBuilder

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
    echo '$paramsJSON' > "params_${unique_id}.log"
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

workflow {
    get_params_and_versions("${params.unique_id}")
}