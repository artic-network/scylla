include { run_kraken_and_bracken } from '../modules/kraken_client'
include { unpack_database } from '../modules/kraken_server'
include { unpack_taxonomy } from '../modules/kraken_server'
include { start_server } from '../modules/kraken_server'
include { stop_server } from '../modules/kraken_server'

workflow kraken_setup {
    take:
        raise_server
    main:
        // Check source param is valid if applicable
        if (!params.database or !params.taxonomy) {
            sources = params.database_sets
            source_name = params.database_set
            source_data = sources.get(source_name, false)
            if (!sources.containsKey(source_name) || !source_data) {
                keys = sources.keySet()
                throw new Exception("Source $params.source is invalid, must be one of $keys")
            }
        }

        // Grab taxonomy files
        if (params.taxonomy) {
            taxonomy = file(params.taxonomy, type: "dir", checkIfExists:true)
        } else {
            default_taxonomy = source_data.get("taxonomy", false)
            if (!default_taxonomy) {
                throw new Exception(
                    "Error: Source $source_name does not include a taxonomy for "
                    + "use with kraken, please choose another source or"
                    + "provide a custom taxonomy.")
            }
            input_taxonomy = file("${params.store_dir}/${params.database_set}/taxonomy_dir")
            if (input_taxonomy.isEmpty()) {
                taxonomy = unpack_taxonomy(default_taxonomy)
            } else {
                taxonomy = input_taxonomy
            }
        }

        // Grab database files
        if (params.database) {
            database = file(params.database, type: "dir", checkIfExists:true)
        } else {
            default_database = source_data.get("database", false)
            if (!default_database) {
                throw new Exception(
                    "Error: Source $source_name does not include a database for "
                    + "use with kraken, please choose another source or "
                    + "provide a custom database.")
            }

            input_database = file("${params.store_dir}/kraken/${params.database_set}/database_dir")
            if (input_database.isEmpty()) {
                database = unpack_database(default_database)
            } else {
                database = input_database
            }
        }

        if (raise_server) {
            start_server(database)
            server = start_server.out.server
        } else {
            server = null
        }
    emit:
        database = database
        taxonomy = taxonomy
        server = server
}

workflow kraken_end {
    take:
        server
        stop
    main:
        stop_server(server, stop)
}

workflow kraken_classify {
    take:
        fastq_ch
        raise_server
    main:
        kraken_setup(raise_server)

        run_kraken_and_bracken(fastq_ch, kraken_setup.out.database, kraken_setup.out.taxonomy)

        if (raise_server)
            kraken_end(kraken_setup.out.server, run_kraken_and_bracken.out.kreport.collect())

    emit:
        assignments = run_kraken_and_bracken.out.assignments
        kreport = run_kraken_and_bracken.out.kreport
        json = run_kraken_and_bracken.out.json
        taxonomy = kraken_setup.out.taxonomy
}

workflow {
    unique_id = "${params.unique_id}"

    // check input fastq exists
    if (params.fastq) {
        fastq = file("${params.fastq}", type: "file", checkIfExists:true)
        if (unique_id == "null") {
            unique_id = "${fastq.simpleName}"
        }
        input_fastq = Channel.fromPath(fastq)
    } else if (params.fastq_dir) {
        fastqdir = file("${params.fastq_dir}", type: "dir", checkIfExists:true)
        if (unique_id == "null") {
            unique_id = "${fastqdir.simpleName}"
        }
        input_fastq = Channel.fromPath( fastqdir / "*.f*q*", type: "file")
    } else {
        exit 1, "One of fastq or fastq_dir need to be provided -- aborting"
    }

    input_fastq.map { it -> [unique_id, it] }.set { fastq_ch }
    kraken_classify(fastq_ch, "${params.raise_server}")
}

