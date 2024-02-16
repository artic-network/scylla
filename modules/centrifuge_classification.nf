// This file contains workflow to classify with centrifuge
include { unpack_taxonomy } from '../modules/kraken_server'

process unpack_database {
    label "process_single"
    storeDir "${params.store_dir}/centrifuge"
    input:
        val database
    output:
        path "database_dir", emit: database
    """
    mkdir database_dir
    cd database_dir
    wget "${database}"
    tar xzf *tar.gz
    cd ..
    """
}

process centrifuge {

    label "process_higher_memory"

    container "docker.io/ontresearch/centrifuge:latest"

    publishDir "${params.outdir}/${unique_id}/classifications", mode: 'copy'

    input:
        tuple val(unique_id), path(fastq)
        path database
    output:
        tuple val(unique_id), path("centrifuge_assignments.tsv"), emit: assignments
        tuple val(unique_id), path("centrifuge_summary.tsv"), emit: summary
    script:
    """
    centrifuge -x "${database}/${params.centrifuge_db_name}" -U ${fastq} \
        -S centrifuge_assignments.tsv \
        --report-file centrifuge_summary.tsv
    """
}

process centrifuge_report {

    label 'process_low'

    container "docker.io/ontresearch/centrifuge:latest"

    publishDir "${params.outdir}/${unique_id}/classifications", mode: 'copy'

    input:
        tuple val(unique_id), path(assignments)
        path database
    output:
        tuple val(unique_id), path("centrifuge.kreport.txt"), emit: kreport
    script:
    """
    centrifuge-kreport -x "${database}/${params.centrifuge_db_name}" ${assignments} > centrifuge.kreport.txt
    """
}

workflow centrifuge_classify {
    take:
        fastq_ch
    main:
        if (params.centrifuge_database) {
                database = file("${params.centrifuge_database}/database_dir", type: "dir", checkIfExists:true)
        } else {
            stored_database = file("${params.store_dir}/centrifuge/database_dir", type: "dir")
            if (stored_database.isEmpty()) {
                unpack_database("${params.centrifuge_remote}")
                database = unpack_database.out.database
            } else {
                database = file("${params.store_dir}/centrifuge/database_dir", type: "dir", checkIfExists:true)
            }
        }

        centrifuge(fastq_ch, database)
        centrifuge_report(centrifuge.out.assignments, database)
    emit:
        kreport = centrifuge_report.out.kreport
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

    Channel.of(unique_id).combine(input_fastq).set{ fastq_ch }
    centrifuge_classify(fastq_ch)
}
