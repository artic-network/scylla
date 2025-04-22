process unpack_database {
    label "process_single"
    storeDir "${params.store_dir}/kraken/${database_name}"
    input:
        val(database_name)
        path(database_path)
    output:
        tuple val(database_name), path("database_dir")

    script:
    """
    if [[ "${database_path}" == *.tar.gz ]]
    then
        mkdir database_dir
        tar xf "${database_path}" -C database_dir
    elif [ -d "${database_path}" ]
    then
        mv "${database_path}" database_dir
    else
        echo "Error: database is neither .tar.gz nor a dir"
        echo "Exiting".
        exit 1
    fi
    """
}


workflow setup_database {
    take:
        database_key
    main:
        databases = params.kraken_database
        database_data = databases.get(database_key)
        if (!databases.containsKey(database_key) || !database_data) {
            keys = databases.keySet()
            throw new Exception("Source $database_key is invalid, must be one of $keys")
        }
        database_name = database_data.get("name")
        database_path = database_data.get("path")

        if (!database_path) {
            sources = params.database_sets
            source_data = sources.get(database_name, false)
            if (!sources.containsKey(database_name) || !source_data) {
                keys = sources.keySet()
                throw new Exception("Source $database_name is invalid, must be one of $keys")
            }
        }

        // Grab database files
        if (database_path) {
            found_database = file(database_path, type: "dir", checkIfExists:true)
            database = Channel.of([database_name, found_database])
        } else {
            default_database = source_data.get("database", false)
            if (!default_database) {
                throw new Exception(
                    "Error: Source $database_name does not include a database for "
                    + "use with kraken, please choose another source or "
                    + "provide a custom database.")
            }

            input_database = file("${params.store_dir}/kraken/$database_name/database_dir")
            if (input_database.isEmpty()) {
                database = unpack_database(database_name, default_database)
            } else {
                database = Channel.of([database_name, input_database])
            }
        }
    emit:
        database = database
}


workflow {
    setup_database("default")
}

