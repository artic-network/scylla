// This file contains workflow to classify with centrifuge
include { unpack_taxonomy } from '../modules/setup_taxonomy'

process unpack_database {
    label "process_single"
    storeDir "${params.store_dir}/centrifuge"
    input:
        path database
    output:
        path "database_dir", emit: database
    """
    if [[ "${database}" == *.tar.gz ]]
    then
        mkdir database_dir
        tar xf "${database}" -C database_dir
    elif [ -d "${database}" ]
    then
        mv "${database}" database_dir
    else
        echo "Error: database is neither .tar.gz nor a dir"
        echo "Exiting".
        exit 1
    fi
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
                remote = file("${params.centrifuge_remote}", checkIfExists:true)
                unpack_database(remote)
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
