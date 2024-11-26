include { kraken_classify; kraken_to_json } from '../modules/kraken_classification'
include { sourmash_classify } from '../modules/sourmash_classification'
include { centrifuge_classify } from '../modules/centrifuge_classification'
include { qc_checks } from '../modules/qc_checks'
include { check_hcid_status } from '../modules/check_hcid_status'
include { check_spike_status } from '../modules/check_spike_status'
include { setup_taxonomy } from '../modules/setup_taxonomy'
include { generate_report } from '../modules/generate_report'


workflow classify_and_report {
    take:
        fastq_ch
        concat_fastq_ch
        database_key
        raise_server
    main:
        qc_checks(fastq_ch)
        kraken_classify(concat_fastq_ch, database_key, raise_server)

        if (params.run_sourmash){
            sourmash_classify(concat_fastq_ch)
        }
        if (params.run_centrifuge){
            centrifuge_classify(concat_fastq_ch)
        }

        setup_taxonomy()
        check_hcid_status(kraken_classify.out.kreport, concat_fastq_ch, setup_taxonomy.out.taxonomy)
        
        if (params.spike_ins) {
            check_spike_status(kraken_classify.out.kreport, concat_fastq_ch)
        }

        kraken_to_json(kraken_classify.out.kreport, setup_taxonomy.out.taxonomy)
        qc_checks.out.combine(kraken_to_json.out, by: 0)
            .join(check_hcid_status.out).set { report_ch }
        generate_report( report_ch )
    emit:
        assignments = kraken_classify.out.assignments
        kreport = kraken_classify.out.kreport
        report = generate_report.out
        taxonomy = setup_taxonomy.out.taxonomy

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
    input_fastq.map { it -> [unique_id, it] }.set { concat_fastq_ch }
    classify_and_report(fastq_ch, concat_fastq_ch, "default", params.raise_server)
}


