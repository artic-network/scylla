include { kraken_classify; kraken_to_json } from '../modules/kraken_classification'
include { sourmash_classify } from '../modules/sourmash_classification'
include { centrifuge_classify } from '../modules/centrifuge_classification'
include { extract_virus_fraction } from '../modules/extract_all'
include { preprocess } from '../modules/preprocess'
include { setup_taxonomy } from '../modules/setup_taxonomy'

workflow default_classify {
    take:
        concat_fastq_ch
        raise_server
    main:
        kraken_classify(concat_fastq_ch, "default", raise_server)

        if (params.run_sourmash){
            sourmash_classify(concat_fastq_ch)
        }
        if (params.run_centrifuge){
            centrifuge_classify(concat_fastq_ch)
        }
    emit:
        assignments = kraken_classify.out.assignments
        kreport = kraken_classify.out.kreport
}

workflow viral_classify {
    take:
        concat_fastq_ch
        raise_server
    main:
        kraken_classify(concat_fastq_ch, "viral", raise_server)
    emit:
        assignments = kraken_classify.out.assignments
        kreport = kraken_classify.out.kreport
}

process merge_classifications {
    label "process_single"

    conda "bioconda::mappy=2.28 bioconda::pyfastx=2.1.0"
    container "community.wave.seqera.io/library/mappy_pyfastx:b4cc4b80f5e5decf"

    publishDir "${params.outdir}/${unique_id}/classifications/", mode: 'copy'

    input:
        tuple val(unique_id), path(default_assignments), path(default_kreport), path(viral_assignments), path(viral_kreport)
    output:
        tuple val(unique_id), val("merged"), path("merged.kraken_assignments.tsv"), emit: assignments
        tuple val(unique_id), val("merged"), path("merged.kraken_report.txt"), emit: kreport
    script:
        """
        ../../../bin/merge.py \
          -a ${default_assignments} ${viral_assignments} \
          -r ${default_kreport} ${viral_kreport}
        """
}

workflow classify {
    take:
        fastq_ch
        concat_fastq_ch
        raise_server
    main:
        default_classify(concat_fastq_ch, raise_server)

        if (params.run_viral_reclassification) {
            setup_taxonomy()
            extract_virus_fraction(fastq_ch, default_classify.out.assignments, default_classify.out.kreport, setup_taxonomy.out)
            viral_classify(extract_virus_fraction.out.virus, raise_server)
            default_classify.out.assignments
                .join(default_classify.out.kreport, by:[0,1])
                .map{ unique_id, database_name, assignments, kreport -> [unique_id, assignments, kreport]}
                .set{ default_ch }
            viral_classify.out.assignments
                .join(viral_classify.out.kreport, by:[0,1])
                .map{ unique_id, database_name, assignments, kreport -> [unique_id, assignments, kreport]}
                .set{ viral_ch }
            default_ch.join( viral_ch , by:0).set{ merge_ch }
            merge_classifications(merge_ch)
            assignments = merge_classifications.out.assignments
            kreport = merge_classifications.out.kreport
        } else {
            assignments = default_classify.out.assignments
            kreport = default_classify.out.kreport
        }

    emit:
        assignments = assignments
        kreport = kreport
}

workflow {
    unique_id = "${params.unique_id}"

    // check input fastq exists
    if (unique_id == "null") {
        if (params.fastq) {
            fastq = file(params.fastq, type: "file", checkIfExists:true)
            unique_id = "${fastq.simpleName}"
        } else if (params.fastq_dir) {
            fastq_dir = file(params.fastq_dir, type: "dir", checkIfExists:true)
            unique_id = "${fastq_dir.simpleName}"
        } else if (params.paired && params.fastq1 && params.fastq2) {
            fastq1 = file(params.fastq1, type: "file", checkIfExists:true)
            unique_id = "${fastq1.simpleName}"
        } else if (params.run_dir) {
            run_dir = file(params.run_dir, type: "dir", checkIfExists:true)
            unique_id = "${run_dir.simpleName}"
        } else {
            exit 1, "One of fastq, fastq_dir or fastq1 and fastq2 need to be provided -- aborting"
        }
    }

    preprocess(unique_id)
    classify(preprocess.out.processed_fastq, preprocess.out.combined_fastq, params.raise_server)
}


