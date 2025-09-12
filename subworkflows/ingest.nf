// workflow to run kraken, check for human, run qc checks and generate html report for a single sample fastq
include { get_params_and_versions } from '../modules/utils'
include { preprocess              } from '../modules/preprocess'
include { classify_and_report     } from '../subworkflows/classify_and_report'
include { extract_fractions; extract_all } from '../modules/extract_all'

workflow ingest {
    take:
    unique_id

    main:
    get_params_and_versions(unique_id)

    preprocess(unique_id)

    classify_and_report(preprocess.out.processed_fastq, preprocess.out.combined_fastq, params.raise_server)
    if (params.skip_extract_taxa)
        extract_fractions(preprocess.out.processed_fastq, classify_and_report.out.assignments, classify_and_report.out.kreport, classify_and_report.out.taxonomy)
    else
        extract_all(preprocess.out.processed_fastq, classify_and_report.out.assignments, classify_and_report.out.kreport, classify_and_report.out.taxonomy)


    emit:
    html_report = classify_and_report.out.report
}
