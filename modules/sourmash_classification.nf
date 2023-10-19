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
    for db in viral archaea bacteria protozoa fungi
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
        tuple val(unique_id), path(fastq)
    output:
        tuple val(unique_id), path("${unique_id}.dna.sig.zip")
    script:
    """
    sourmash sketch dna ${fastq} -p k=${params.sourmash_k},dna,scaled=1000,abund \
                                        --name ${unique_id} -o "${unique_id}.dna.sig.zip"
    """
}

process sourmash_gather {

    label 'error_retry'
    memory { 60.GB * task.attempt }

    conda "bioconda::sourmash=4.8.4"
    container "quay.io/biocontainers/sourmash:4.8.4--hdfd78af_0"

    input:
        tuple val(unique_id), path(sketch)
        path database
    output:
        tuple val(unique_id), path("${unique_id}.k${params.sourmash_k}.gather.csv"), emit: gather_csv
        tuple val(unique_id), path("${unique_id}.k${params.sourmash_k}.gather.txt"), emit: gather_txt
    script:
    """
    sourmash gather ${sketch} database_dir/${params.sourmash_db_name}-{viral,protozoa}-k${params.sourmash_k}.zip --dna --ksize ${params.sourmash_k} \
                     --threshold-bp ${params.sourmash_threshold_bp} \
                     -o "${unique_id}.k${params.sourmash_k}.gather.csv" > "${unique_id}.k${params.sourmash_k}.gather.txt"
    """
}

process sourmash_tax_metagenome {

    label 'process_low'
    label 'error_retry'

    conda "bioconda::sourmash=4.8.4"
    container "quay.io/biocontainers/sourmash:4.8.4--hdfd78af_0"
    publishDir path: "${params.outdir}/${unique_id}/classifications", mode: 'copy', pattern: '*.csv'

    input:
        tuple val(unique_id), path(gather_csv)
        path lineages
    output:
        tuple val(unique_id), path("sourmash.k${params.sourmash_k}.unformatted.kreport.txt"), emit: kreport
        path "sourmash.k${params.sourmash_k}.summarized.csv", emit: summary
    script:
    """
    sourmash tax metagenome -g ${gather_csv} \
        -t lineages_dir/${params.sourmash_db_name}-{viral,protozoa}.lineages.csv.gz \
        -o "sourmash.k${params.sourmash_k}" \
        --output-format krona csv_summary kreport \
        --rank species
    mv "sourmash.k${params.sourmash_k}.kreport.txt" "sourmash.k${params.sourmash_k}.unformatted.kreport.txt"
    """
}

process sourmash_to_json {

    label "process_low"

    publishDir path: "${params.outdir}/${unique_id}/classifications", mode: 'copy'

    conda "bioconda::biopython=1.78 anaconda::Mako=1.2.3"
    container "${params.wf.container}@${params.wf.container_sha}"


    input:
        tuple val(unique_id), path(kraken_report)
        path taxonomy_dir
    output:
       tuple val(unique_id), path("${params.database_set}.kraken.json")
       tuple val(unique_id), path("sourmash.k${params.sourmash_k}.kreport.txt"), emit: kreport

    """
    reformat_sourmash_kreport.py -r ${kraken_report} -t ${taxonomy_dir} -o "reformatted_kreport.txt"

    cat "reformatted_kreport.txt" | cut -f5,3 | tail n+2 > taxacounts.txt
    cat "reformatted_kreport.txt" | cut -f5 | tail n+2 > taxa.txt
    taxonkit lineage --data-dir ${taxonomy_dir}  -R taxa.txt  > lineages.txt
    aggregate_lineages_bracken.py \\
            -i "lineages.txt" -b "taxacounts.txt" \\
            -u "reformatted_kreport.txt" \\
            -p "temp_kraken"
    file1=`cat *.json`
    echo "{"'Sourmash'": "\$file1"}" >> "Sourmash.kraken.json"
    """
}

workflow sourmash_classify {
    take:
        fastq_ch
    main:
        if (params.sourmash_database) {
                database = file("${params.sourmash_database}/database_dir", type: "dir", checkIfExists:true)
                lineages = file("${params.sourmash_database}/lineages_dir", type: "dir", checkIfExists:true)
        } else {
            stored_database = file("${params.store_dir}/sourmash/database_dir", type: "dir")
            if (stored_database.isEmpty()) {
                unpack_database("${params.sourmash_remote}")
                database = unpack_database.out.database
                lineages = unpack_database.out.database
            } else {
                database = file("${params.store_dir}/sourmash/database_dir", type: "dir", checkIfExists:true)
                lineages = file("${params.store_dir}/sourmash/lineages_dir", type: "dir", checkIfExists:true)
            }
        }

        input_taxonomy = file("${params.store_dir}/${params.database_set}/taxonomy_dir")
        if (input_taxonomy.isEmpty()) {
            taxonomy = unpack_taxonomy(default_taxonomy)
        } else {
            taxonomy = input_taxonomy
        }

        sourmash_sketch_dna(fastq_ch)
        sourmash_gather(sourmash_sketch_dna.out, database)
        sourmash_tax_metagenome(sourmash_gather.out.gather_csv, lineages)
        sourmash_to_json(sourmash_tax_metagenome.out.kreport, taxonomy)

    emit:
        kreport = sourmash_tax_metagenome.out.kreport
        json = sourmash_to_json.out
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
    sourmash_classify(fastq_ch)
}
