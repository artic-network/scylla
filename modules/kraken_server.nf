// This file contains workflow to startup kraken_server for a given database,
// downloading and unpacking the database if required.
// Also workflow to shutdown kraken_server for a given database
//
// Notes on CPU resource of kraken server and client:
//
// - server will use as much resource as number of clients running,
//   plus one extra thread => set maxForks of clients to threads - 1.
//   (not doing so gives gRPC server: "Server Threadpool Exhausted")
// - we need potentially one extra request thread to allow stop request
//   to be handled
// - cannot start so many clients (or other processes) such that
//   server never starts from Nextflow executor limit
// - we'd like to leave some resource for downstream processes such that we
//   get reporting updated frequently
// - this might all be considered a bit inefficient - but we are set
//   up for real-time dynamism not speed
//



process unpack_database {
    label "process_single"
    storeDir "${params.store_dir}/${params.database_set}"
    input:
        path database
    output:
        path "database_dir"
    """
    if [[ "${database}" == *.tar.gz ]]
    then
        mkdir database_dir
        tar xf "${database}" -C database_dir
    elif [ -d "${database}" ]
    then
        mv "${database}" database_dir
    else
        echo "Error: database is neither .tar.gz nor a dir"
        echo "Exiting".
        exit 1
    fi
    """
}

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

kraken_compute = params.kraken_clients == 1 ? 1 : params.kraken_clients - 1

process kraken_server {
    label "process_long"
    memory { 8.GB * task.attempt }
    cpus params.threads
    container "${params.wf.container}@${params.wf.container_version}"
    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    input:
        path database
    output:
        val true
    script:
    """
    # we add one to requests to allow for stop signal
    kraken2_server \
        --max-requests ${kraken_compute + 1} --thread-pool ${params.server_threads}\
        --port ${params.k2_port} \
        --host-ip ${params.k2_host} \
        --db ./${database}/
    """
}


process stop_kraken_server {
    label "process_single"
    container "${params.wf.container}@${params.wf.container_version}"
    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    // this shouldn't happen, but we'll keep retrying
    // errorStrategy = { task.exitStatus in [8, 14] && task.attempt < 3 ? 'retry' : 'ignore' }
    errorStrategy 'ignore'
    input:
        val stop
    """
    kraken2_client --port ${params.k2_port} --host-ip ${params.k2_host} --shutdown
    """
}

workflow start_server {
    take:
        database
    main:
        kraken_server(database)
    emit:
        server = kraken_server.out
}

workflow stop_server {
    take:
        server
        stop
    main:
        stop_kraken_server(stop)
}

workflow {
    // Grab database files
    if (params.database) {
        database = file(params.database, type: "dir", checkIfExists:true)
    } else {
            sources = params.database_sets
            source_name = params.database_set
            source_data = sources.get(source_name, false)
            if (!sources.containsKey(source_name) || !source_data) {
                keys = sources.keySet()
                throw new Exception("Source $params.source is invalid, must be one of $keys")
            }

            default_database = source_data.get("database", false)
            if (!default_database) {
                throw new Exception(
                    "Error: Source $source_name does not include a database for "
                    + "use with kraken, please choose another source or "
                    + "provide a custom database.")
            }

            input_database = file("${params.store_dir}/${params.database_set}/database_dir")
            if (input_database.isEmpty()) {
                database = unpack_database(default_database)
            } else {
                database = input_database
            }
    }

    start_server(params.database)
}
