// workflow to run kraken, check for human, run qc checks and generate html report for a single sample fastq
include { get_params_and_versions ; get_fastq_channels } from '../modules/utils'

include { preprocess              } from '../modules/preprocess'
include { qc_checks               } from '../modules/qc_checks'
include { centrifuge_classify     } from '../modules/centrifuge_classification'
include { sourmash_classify       } from '../modules/sourmash_classification'
include { check_hcid_status       } from '../modules/check_hcid_status'
include { check_spike_status      } from '../modules/check_spike_status'
include { extract_all             } from '../modules/extract_all'
include { classify_virus_fastq    } from '../modules/classify_novel_viruses'
include { classify ; reclassify   } from '../subworkflows/classify'

workflow run_module {
    take:
    unique_id

    main:
    get_params_and_versions(unique_id)
    get_fastq_channels(unique_id)
    fastq_ch = get_fastq_channels.out.processed_fastq

    if (params.module == "preprocess") {
        preprocess(unique_id)
    }
    else if (params.module == "qc_checks") {
        qc_checks(fastq_ch)
    }
    else if (params.module == "centrifuge_classification") {
        centrifuge_classify(fastq_ch)
    }
    else if (params.module == "kraken_classification") {
        classify(fastq_ch, get_fastq_channels.out.combined_fastq, params.raise_server)
    }
    else if (params.module == "sourmash_classification") {
        sourmash_classify(fastq_ch)
    }
    else if (params.module == "check_hcid_status") {
        kreport = file(params.kraken_report, type: "file", checkIfExists: true)
        kreport_ch = Channel.of([unique_id, "default", kreport])
        taxonomy_dir = file(params.taxonomy, type: "dir", checkIfExists: true)
        check_hcid_status(kreport_ch, fastq_ch, taxonomy_dir)
    }
    else if (params.module == "check_spike_status") {
        kreport = file(params.kraken_report, type: "file", checkIfExists: true)
        kreport_ch = Channel.of([unique_id, "default", kreport])
        check_spike_status(kreport_ch, fastq_ch)
    }
    else if (params.module == "extract_all") {
        assignments = file(params.kraken_assignments, type: "file", checkIfExists: true)
        assignments_ch = Channel.of([unique_id, "default", assignments])
        kreport = file(params.kraken_report, type: "file", checkIfExists: true)
        kreport_ch = Channel.of([unique_id, "default", kreport])
        taxonomy_dir = file(params.taxonomy, type: "dir", checkIfExists: true)
        extract_all(fastq_ch, assignments_ch, kreport_ch, taxonomy_dir)
    }
     else if (params.module == "kraken_reclassification") {
        assignments = file(params.kraken_assignments, type: "file", checkIfExists: true)
        assignments_ch = Channel.of([unique_id, "default", assignments])
        kreport = file(params.kraken_report, type: "file", checkIfExists: true)
        kreport_ch = Channel.of([unique_id, "default", kreport])
        reclassify(fastq_ch, assignments_ch, kreport_ch, params.raise_server)
    }
    else if (params.module == "classify_novel_viruses") {
        classify_virus_fastq(fastq_ch)
    }
    else {
        exit(1, "Supported modules are [preprocess, qc_checks, centrifuge_classification, kraken_classification, sourmash_classification, check_hcid_status, check_spike_status, extract_all, classify_novel_viruses]")
    }
}
