include { classify           } from '../subworkflows/classify'
include { kraken_to_json     } from '../modules/kraken_classification'
include { qc_checks          } from '../modules/qc_checks'
include { check_hcid_status  } from '../modules/check_hcid_status'
include { check_spike_status } from '../modules/check_spike_status'
include { setup_taxonomy     } from '../modules/setup_taxonomy'
include { generate_report    } from '../modules/generate_report'


workflow classify_and_report {
    take:
    fastq_ch
    concat_fastq_ch
    raise_server

    main:
    qc_checks(fastq_ch)

    classify(fastq_ch, concat_fastq_ch, raise_server)

    setup_taxonomy()
    check_hcid_status(classify.out.kreport, concat_fastq_ch, setup_taxonomy.out.taxonomy)

    if (params.spike_ins) {
        check_spike_status(classify.out.kreport, concat_fastq_ch)
    }

    kraken_to_json(classify.out.kreport, setup_taxonomy.out.taxonomy)
    qc_checks.out
        .join(kraken_to_json.out, failOnMismatch: true, failOnDuplicate: true)
        .join(check_hcid_status.out, failOnMismatch: false, failOnDuplicate: true, remainder: true)
        .map { uniqueId, stats, database_name, lineages, warnings ->
            [
                uniqueId,
                stats,
                database_name,
                lineages,
                warnings ? warnings : [],
            ]
        }
        .set { report_ch }
    generate_report(report_ch)

    emit:
    assignments = classify.out.assignments
    kreport     = classify.out.kreport
    report      = generate_report.out
    taxonomy    = setup_taxonomy.out.taxonomy
}

workflow {
    unique_id = "${params.unique_id}"

    // check input fastq exists
    if (params.fastq) {
        fastq = file("${params.fastq}", type: "file", checkIfExists: true)
        if (unique_id == "null") {
            unique_id = "${fastq.simpleName}"
        }
        input_fastq = Channel.fromPath(fastq)
    }
    else if (params.fastq_dir) {
        fastqdir = file("${params.fastq_dir}", type: "dir", checkIfExists: true)
        if (unique_id == "null") {
            unique_id = "${fastqdir.simpleName}"
        }
        input_fastq = Channel.fromPath(fastqdir / "*.f*q*", type: "file")
    }
    else {
        exit(1, "One of fastq or fastq_dir need to be provided -- aborting")
    }

    input_fastq.map { it -> [unique_id, it] }.set { fastq_ch }
    input_fastq.map { it -> [unique_id, it] }.set { concat_fastq_ch }
    classify_and_report(fastq_ch, concat_fastq_ch, "default", params.raise_server)
}
