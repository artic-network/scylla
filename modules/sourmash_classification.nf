// This file contains workflow to classify with sourmash

process unpack_database {
    label "process_single"
    storeDir "${params.store_dir}/sourmash"
    input:
        val database
    output:
        path "database_dir", emit: database
        path "lineages_dir", emit: lineages
    """
    mkdir database_dir
    cd database_dir
    for db in viral archaea bacteria protozoa fungi
    do
        curl -JLO "${database}-\$db-k${params.sourmash_k}.zip"
    done
    cd ..

    mkdir lineages_dir
    cd lineages_dir
    for db in viral protozoa
    do
        curl -JLO "${database}-\$db.lineages.csv.gz"
    done
    cd ..

    """
}


process sourmash_sketch_dna {

    label 'process_low'
    label 'error_retry'

    conda "bioconda::sourmash=4.8.4"
    container "quay.io/biocontainers/sourmash:4.8.4--hdfd78af_0"

    input:
        val unique_id
        path fastq
    output:
        path "${unique_id}.dna.sig.zip"
    script:
    """
    sourmash sketch dna ${fastq} -p k=${params.sourmash_k},dna,scaled=1000,abund \
                                        --name ${unique_id} -o "${unique_id}.dna.sig.zip"
    """
}

process sourmash_gather {

    label 'process_medium'
    label 'error_retry'

    conda "bioconda::sourmash=4.8.4"
    container "quay.io/biocontainers/sourmash:4.8.4--hdfd78af_0"

    input:
        val unique_id
        path sketch
        path database
    output:
        path "${unique_id}.k${params.sourmash_k}.gather.csv", emit: gather_csv
        path "${unique_id}.k${params.sourmash_k}.gather.txt", emit: gather_txt
    script:
    """
    sourmash gather ${sketch} database_dir/${params.sourmash_db_name}-{viral,archaea,bacteria,protozoa,fungi}-k${params.sourmash_k}.zip --dna --ksize ${params.sourmash_k} \
                     --threshold-bp ${params.sourmash_threshold_bp} \
                     -o "${unique_id}.k${params.sourmash_k}.gather.csv" > "${unique_id}.k${params.sourmash_k}.gather.txt"
    """
}

process sourmash_tax_metagenome {

    label 'process_low'
    label 'error_retry'

    conda "bioconda::sourmash=4.8.4"
    container "quay.io/biocontainers/sourmash:4.8.4--hdfd78af_0"
    publishDir path: "${params.outdir}/${unique_id}/classifications", mode: 'copy'

    input:
        val unique_id
        path gather_csv
        path lineages
    output:
        path "sourmash.k${params.sourmash_k}.kreport.txt", emit: kreport
        path "sourmash.k${params.sourmash_k}.summarized.csv", emit: summary
    script:
    """
    sourmash tax metagenome -g ${gather_csv} \
        -t lineages_dir/${params.sourmash_db_name}-{viral,archaea,bacteria,protozoa,fungi}.lineages.csv.gz \
        -o "sourmash.k${params.sourmash_k}" \
        --output-format krona csv_summary kreport \
        --rank species
    """
}

workflow sourmash_classify {
    take:
        unique_id
        fastq
    main:
        stored_database = file("${params.store_dir}/sourmash/database_dir", type: "dir")
        if (stored_database.isEmpty()) {
            unpack_database("${params.sourmash_database}")
            database = unpack_database.out.database
            lineages = unpack_database.out.database
        } else {
            database = file("${params.store_dir}/sourmash/database_dir", type: "dir", checkIfExists:true)
            lineages = file("${params.store_dir}/sourmash/lineages_dir", type: "dir", checkIfExists:true)
        }

        sourmash_sketch_dna(unique_id, fastq)
        sourmash_gather(unique_id, sourmash_sketch_dna.out, database)
        sourmash_tax_metagenome(unique_id, sourmash_gather.out.gather_csv, lineages)
    emit:
        genbank = sourmash_tax_metagenome.out.kreport
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

    sourmash_classify(unique_id, input_fastq)
}