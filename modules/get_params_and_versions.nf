// Module to print params and versions of software

import groovy.json.JsonBuilder

process get_versions {

    
    publishDir "${params.out_dir}/${unique_id}/execution", mode: 'copy'
    cpus 1
    input:
        val unique_id
    output:
        path "versions.txt"
    script:
    """
    kraken2 --version | head -n 1 | sed 's/ version /,/' >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    taxonkit version | sed 's/ /,/' >> versions.txt
    """
}


process get_params {
    publishDir "${params.out_dir}/${unique_id}/execution", mode: 'copy'
    cpus 1
    input:
        val unique_id
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}

workflow get_params_and_versions {
    take:
        unique_ud
    main:
        software_versions = get_versions(unique_id)
        workflow_params = get_params(unique_id)
    emit:
        software_versions
        workflow_params
}

workflow {
    get_params_and_versions()
}