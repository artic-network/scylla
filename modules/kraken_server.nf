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



process unpackDatabase {
    label "wfmetagenomics"
    cpus 1
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

process unpackTaxonomy {
    cpus 1
    storeDir "${params.store_dir}/${taxonomy.simpleName}"
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

kraken_compute = params.threads == 1 ? 1 : params.threads - 1

process kraken_server {
    errorStrategy 'ignore'
    label "scylla"
    cpus params.threads
    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    input:
        path database
    output:
        val true
    script:
    """
    # we add one to requests to allow for stop signal
    kraken2_server \
        --max-requests ${kraken_compute + 1} --port ${params.port} \
        --db ./${database}/
    """
}


process stop_kraken_server {
    label "scylla"
    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    // this shouldn't happen, but we'll keep retrying
    // errorStrategy = { task.exitStatus in [8, 14] && task.attempt < 3 ? 'retry' : 'ignore' }
    errorStrategy 'ignore'
    input:
        val stop
    """
    kraken2_client --port $params.port --shutdown
    """
}

workflow start_server {
    take:
        database
        taxonomy
    main:
        input_database = file("${params.store_dir}/${params.database_set}/database_dir")
        if (input_database.isEmpty()) {
            database = unpackDatabase(database)
        } else {
            database = input_database
        }

        input_taxonomy = file("${params.store_dir}/${taxonomy.simpleName}/taxonomy_dir")
        if (input_database.isEmpty()) {
            taxonomy = unpackTaxonomy(taxonomy)
        } else {
            taxonomy = input_taxonomy
        }
        kraken_server(database)
    emit:
        database = database
        taxonomy = taxonomy
        server = kraken_server.out
}

workflow stop_server {
    take:
        stop
    main:
        stop_kraken_server(stop)
}