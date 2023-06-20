include { run_kraken_and_bracken } from '../modules/kraken_client'
include { qc_checks } from '../modules/qc_checks'
include { generate_report } from '../modules/generate_report'

workflow kraken_pipeline {
    take:
        unique_id
        fastq
    main:

        // Check source param is valid
        // sources = params.database_sets
        // source_name = params.database_set
        // source_data = sources.get(source_name, false)
        // if (!sources.containsKey(source_name) || !source_data) {
        //     keys = sources.keySet()
        //     throw new Exception("Source $params.source is invalid, must be one of $keys")
        // }

        // Grab taxonomy files
        // taxonomy = file(sources[source_name]["taxonomy"], type: "file")
        // if (params.taxonomy) {
        //     log.info("Checking custom taxonomy mapping exists")
            
        // }

        taxonomy = file(params.taxonomy_dir, type: "dir", checkIfExists:true)

        // if (database) {
        //     source_database = database
        // } else {
        //     source_database = source_data.get("database", false)
        //     if (!source_database) {
        //         throw new Exception(
        //             "Error: Source $source_name does not include a database for "
        //             + "use with kraken, please choose another source, "
        //             + "provide a custom database or disable kraken.")
        //     }
        database = file(params.db, type: "dir", checkIfExists:true)

        // start_server(database, taxonomy)
        run_kraken_and_bracken(unique_id, fastq, database, taxonomy)

        // stop_server(run_kraken_and_bracken.out.bracken_report.collect())
        qc_checks(unique_id, fastq)
        if (params.additional_bracken_jsons) {
            Channel.of(file(params.additional_bracken_jsons, type: "file", checkIfExists:true))
                .concat(run_kraken_and_bracken.out.json)
                .unique {it.getName()}
                .collect()
                .set { bracken_jsons }
        } else {
            run_kraken_and_bracken.out.json
                .collect()
                .set { bracken_jsons }
        }

        generate_report(unique_id, qc_checks.out, bracken_jsons )
    emit:
        bracken_report = run_kraken_and_bracken.out.bracken_report
        kraken_report = run_kraken_and_bracken.out.kraken_report
        kraken_assignments = run_kraken_and_bracken.out.kraken_assignments
        report = generate_report.out
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
        input_fastq = Channel.fromPath( fastqdir / "*.f*q", type: "file")
    } else {
        exit 1, "One of fastq or fastq_dir need to be provided -- aborting"
    }

    kraken_pipeline(unique_id, input_fastq)
}


