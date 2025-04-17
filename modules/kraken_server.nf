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
include { setup_database } from '../modules/setup_database'

kraken_compute = params.kraken_clients == 1 ? 1 : params.kraken_clients - 1

process start_kraken_server {
    label "process_long"
    label "process_high_memory"
    cpus {
        Integer max_local_threads = workflow.session.config?.executor?.$local?.cpus ?: \
                    Runtime.getRuntime().availableProcessors()
        if (max_local_threads == 1) {
            throw new Exception("Cannot run kraken_server and kraken_client at the same time as the local executor appears to be configured with only one CPU.")
        } else if (max_local_threads == 2) {
            // run the server single threaded and expect one client
            log.info("Automatically set kraken2 classification server threads to 1 to ensure a classification client can be run.")
            1
        } else {
            // remove one thread for at least one client, and another for other business
            log.info("Set kraken2 classification server threads to ${Math.min(params.threads, max_local_threads - 2)} to ensure a classification client can be run.")
            Math.min(params.threads, max_local_threads - 2)
        }
    }

    conda "nanoporetech::kraken2-server=0.1.7"
    container "${params.wf.container}:${params.wf.container_version}"
    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    input:
        tuple val(database_name), path(database), val(k2_host), val(k2_port)
    output:
        tuple val(database_name), val(k2_host), val(k2_port)
    script:
    """
    # we add one to requests to allow for stop signal
    kraken2_server \
        --max-requests ${kraken_compute + 1} --thread-pool ${task.cpus}\
        --port ${k2_port} \
        --host-ip ${k2_host} \
        --db ./${database}/ \
        --confidence ${params.kraken_confidence}
    """
}


process stop_kraken_server {
    label "process_single"
    conda "nanoporetech::kraken2-server=0.1.7"
    container "${params.wf.container}:${params.wf.container_version}"
    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    // this shouldn't happen, but we'll keep retrying
    // errorStrategy = { task.exitStatus in [8, 14] && task.attempt < 3 ? 'retry' : 'ignore' }
    errorStrategy 'ignore'
    input:
        tuple val (database_name), val(k2_host), val(k2_port)
        val (stop)
    """
    kraken2_client --port ${k2_port} --host-ip ${k2_host} --shutdown
    """
}

workflow get_server_address {
    take:
        server_key
    main:
        kraken_servers = params.kraken_database
        server_data = kraken_servers.get(server_key)
        if (!kraken_servers.containsKey(server_key) || !server_data) {
            keys = kraken_servers.keySet()
            throw new Exception("Source $server_key is invalid, must be one of $keys")
        }
        name = server_data.get("name")
        host = server_data.get("host")
        port = server_data.get("port")
        server_address = [name, host, port]
        //println(server_address)
    emit:
        server_address = server_address
}


workflow get_server {
    take:
        database_key
        raise_server
    main:
        setup_database(database_key)
        get_server_address(database_key)
        database_ch = setup_database.out.database
        server_ch = get_server_address.out.server_address
        if (raise_server){
            println("Raising server for ${database_key}")
            combined_ch = database_ch.combine(server_ch, by: 0)
            //combined_ch.view()
            start_kraken_server(combined_ch)
        }
    emit:
        database = database_ch
        server = server_ch
}

workflow get_default_server {
    take:
        raise_server
    main:
        get_server("default", raise_server)
    emit:
        database = get_server.out.database
        server = get_server.out.server
}

workflow get_viral_server {
    take:
        raise_server
    main:
        get_server("viral", raise_server)
    emit:
        database = get_server.out.database
        server = get_server.out.server
}

workflow stop_server {
    take:
       server
       stop
    main:
        stop_kraken_server(server, stop)
}

workflow stop_default_server {
    take:
       server
       stop
    main:
        println("Stopping server for ${params.kraken_database.default.name}")
        stop_kraken_server(server, stop)
}

workflow stop_viral_server {
    take:
       server
       stop
    main:
        println("Stopping server for ${params.kraken_database.viral.name}")
        stop_kraken_server(server, stop)
}

workflow {
    get_server("default", true)
}
