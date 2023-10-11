include { run_kraken_and_bracken } from '../modules/kraken_client'
include { unpack_database } from '../modules/kraken_server'
include { unpack_taxonomy } from '../modules/kraken_server'
include { start_server } from '../modules/kraken_server'
include { stop_server } from '../modules/kraken_server'
include { qc_checks } from '../modules/qc_checks'
include { generate_report } from '../modules/generate_report'

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

            input_database = file("${params.store_dir}/${params.database_set}/database_dir")
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

workflow kraken_pipeline {
    take:
        unique_id
        fastq
        raise_server
    main:
        kraken_setup(raise_server)

        run_kraken_and_bracken(unique_id, fastq, kraken_setup.out.database, kraken_setup.out.taxonomy)

        if (raise_server)
            kraken_end(kraken_setup.out.server, run_kraken_and_bracken.out.bracken_report.collect())

        qc_checks(unique_id, fastq)
        if (params.additional_bracken_jsons) {
            Channel.of(file(params.additional_bracken_jsons, type: "file", checkIfExists:true))
                .concat(run_kraken_and_bracken.out.json)
                .unique {it.getName()}
                .flatten()
                .set { bracken_jsons }
        } else {
            run_kraken_and_bracken.out.json
                .flatten()
                .set { bracken_jsons }
        }

        generate_report(unique_id, qc_checks.out, bracken_jsons )
    emit:
        bracken_report = run_kraken_and_bracken.out.bracken_report
        kraken_report = run_kraken_and_bracken.out.kraken_report
        kraken_assignments = run_kraken_and_bracken.out.kraken_assignments
        report = generate_report.out
        taxonomy = kraken_setup.out.taxonomy
        unique_id = unique_id
        fastq = fastq
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

    kraken_pipeline(unique_id, input_fastq)
}


