include { kraken_classify } from '../modules/kraken_classification'
include { qc_checks } from '../modules/qc_checks'
include { generate_report } from '../modules/generate_report'


workflow classify_and_report {
    take:
        fastq_ch
        raise_server
    main:
        qc_checks(fastq_ch)
        kraken_classify(fastq_ch, raise_server)

        if (params.additional_bracken_jsons) {
            jsons = Channel.of(file(params.additional_bracken_jsons, type: "file", checkIfExists:true))
            fastq_ch.map{ it -> it[0] }
                .combine(jsons)
                .concat(kraken_classify.out.json)
                .groupTuple()
                .set { classified_jsons }
        } else {
            kraken_classify.out.json
                .set { classified_jsons }
        }

        qc_checks.out.join(classified_jsons).set { report_ch }
        generate_report( report_ch )
    emit:
        assignments = kraken_classify.out.assignments
        kreport = kraken_classify.out.kreport
        report = generate_report.out
        taxonomy = kraken_classify.out.taxonomy

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

    fastq_ch = [unique_id, input_fastq]
    classify_and_report(fastq_ch)
}


