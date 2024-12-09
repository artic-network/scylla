process unpack_taxonomy {
    label "process_single"
    storeDir "${params.store_dir}"
    input:
        path taxonomy
    output:
        path "taxonomy_dir"
    """
    if [[ "${taxonomy}" == *.tar.gz ]]
    then
        mkdir taxonomy_dir
        tar xf "${taxonomy}" -C taxonomy_dir
    elif [ -d "${taxonomy}" ]
    then
        mv "${taxonomy}" taxonomy_dir
    else
        echo "Error: taxonomy is neither .tar.gz nor a dir"
        echo "Exiting".
        exit 1
    fi
    """
}


workflow setup_taxonomy {
    main:
        if (params.taxonomy) {
            taxonomy = file(params.taxonomy, type: "dir", checkIfExists:true)
        } else {
            input_taxonomy = file("${params.store_dir}/taxonomy_dir")
            if (input_taxonomy.isEmpty()) {
                taxonomy = unpack_taxonomy(params.default_taxonomy)
            } else {
                taxonomy = input_taxonomy
            }
        }
    emit:
        taxonomy = taxonomy
}


workflow {
    setup_taxonomy()
}

