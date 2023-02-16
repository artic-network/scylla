// Module to print params and versions of software

import groovy.json.JsonBuilder

process get_versions {
    label "scylla"
    publishDir "${params.out_dir}/execution", mode: 'copy'
    cpus 1
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
    publishDir "${params.out_dir}/execution", mode: 'copy'
    cpus 1
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
    main:
        software_versions = get_versions()
        workflow_params = get_params()
    emit:
        software_versions
        workflow_params
}

workflow {
    get_params_and_versions()
}