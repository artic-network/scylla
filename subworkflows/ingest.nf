// workflow to run kraken, check for human, run qc checks and generate html report for a single sample fastq
include { get_params_and_versions } from '../modules/get_params_and_versions'
include { preprocess } from '../modules/preprocess'
include { classify_and_report } from '../subworkflows/classify_and_report'
include { extract_taxa } from '../modules/extract_taxa'
include { classify_novel_taxa; classify_novel_taxa_paired } from '../modules/classify_novel_taxa'

workflow ingest {
    take:
        unique_id
    main:
        get_params_and_versions(unique_id)

        preprocess(unique_id)

        classify_and_report(preprocess.out.combined_fastq, params.raise_server)
        extract_taxa(preprocess.out.processed_fastq, classify_and_report.out.assignments, classify_and_report.out.kreport, classify_and_report.out.taxonomy)

        if (params.classify_novel_viruses) {
            if (params.paired) {
                classify_novel_taxa_paired(unique_id, fastq_1, fastq_2)
            } else {
                classify_novel_taxa(unique_id, fastq)
            }
        }

    emit:
        html_report = classify_and_report.out.report

}
