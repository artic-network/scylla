include { start_server } from '../modules/kraken_server'
include { stop_server } from '../modules/kraken_server'
include { run_kraken_and_bracken } from '../modules/kraken_client'
include { qc_checks } from '../modules/qc_checks'
include { generate_report } from '../modules/generate_report'

workflow kraken_pipeline {
    take:
        unique_id
        fastq
    main:

        // Check source param is valid
        sources = params.database_sets
        source_name = params.database_set
        source_data = sources.get(source_name, false)
        if (!sources.containsKey(source_name) || !source_data) {
            keys = sources.keySet()
            throw new Exception("Source $params.source is invalid, must be one of $keys")
        }

        // Grab taxonomy files
        taxonomy = file(sources[source_name]["taxonomy"], type: "file")
        if (params.taxonomy) {
            log.info("Checking custom taxonomy mapping exists")
            taxonomy = file(params.taxonomy, type: "dir", checkIfExists:true)
        }

        source_database = source_data.get("database", false)
        if (!source_database) {
            throw new Exception(
                "Error: Source $source_name does not include a database for "
                + "use with kraken, please choose another source, "
                + "provide a custom database or disable kraken.")
        }
        database = file(source_database, type: "file")

        start_server(database, taxonomy)
        run_kraken_and_bracken(unique_id, fastq, start_server.out.database, start_server.out.taxonomy)
        stop_server(run_kraken_and_bracken.out.bracken_report.collect())
        qc_checks(unique_id, fastq)
        all_bracken_jsons = Channel.fromPath("${params.out_dir}/${unique_id}/classifications/*.bracken.json")
                    .concat(run_kraken_and_bracken.out.json)
                    .unique {it.getName()}
                    .collect()
        generate_report(unique_id, qc_checks.out, all_bracken_jsons )
    emit:
        bracken_report = run_kraken_and_bracken.out.bracken_report
        kraken_report = run_kraken_and_bracken.out.kraken_report
        kraken_assignments = run_kraken_and_bracken.out.kraken_assignments
        report = generate_report.out
}

workflow {
    // check input fastq exists
    input_fastq = file("${params.fastq}", type: "file", checkIfExists:true)

    unique_id = "${params.unique_id}"
    if (unique_id == "null") {
        unique_id = "${input_fastq.simpleName}"
    }

    kraken_pipeline(unique_id, input_fastq)
}


