include { start_server } from '../modules/kraken_server'
include { stop_server } from '../modules/kraken_server'
include { run_kraken_and_bracken } from '../modules/kraken_client'
include { qc_checks } from '../modules/qc_checks'
include { generate_report } from '../modules/generate_report'

workflow {
    // check input fastq exists
    input_fastq = file("${params.fastq}")
    if (!input_fastq.exists()) {
            throw new Exception("--fastq: File doesn't exist, check path.")
        }
    sample_id = "${params.sample_id}"
    if (sample_id == "null") {
        sample_id = "${input_fastq.simpleName}"
    }

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
    run_kraken_and_bracken(sample_id, input_fastq, source_name, start_server.out.database, start_server.out.taxonomy)
    stop_server(run_kraken_and_bracken.out.report.collect())
    qc_checks(sample_id, input_fastq)
    generate_report(sample_id, qc_checks.out, run_kraken_and_bracken.out.json)
}


