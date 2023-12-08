// This workflow is for trying to find and classify possibly novel RNA viruses that might have been overlooked by kraken. 
// 1) runs a metagenomic assembly on short or long reads,
// 2) runs VirBot or GeNomad classification on them.             
// Output: file with contigs identified as RNA viral 'viral_contigs.fa',
//     contig id to taxonomy mapping file 'tax_assignments.tsv',
//     (optional) {assembler}_assembly_stats.txt file with the performed assembly statistics -
//     all in the 'discovery' directory

// Related options and params:
// --classifier) 'virbot'|'genomad', default: virbot
// Notes: genomad requires longer contigs for its neural network to run,
//        so the contigs are filtered by length (>=2kbp) beforehand.
//        Therefore, it is adviced to use genomad on long read data or relatively contiguous assembly
//    overall, at times virbot has greater sensitivity while genomad might be more precise (?),
//        but generally they produce similar output
//    virbot does classification to the lowest level possible
//    genomad here is run on preset 'relaxed' (no end filtering)
//    for genomad, the assignments other than 'Riboviria' (including those labeled as 'Unclassified')
//        are excluded from the final taxonomy and contigs files
// --outdir) output directory of the pipeline; REQUIRED for testing
// --genomad_db) path to the directory with genomad database

include { assemble_taxa } from '../modules/assemble_taxa'

process run_virbot {
    label "process_medium"
    errorStrategy 'ignore'

    container "biowilko/virbot:1.0.0"
    conda "bioconda:virbot:1.0.0"
    publishDir "${params.outdir}/${unique_id}/discovery/${taxon}", mode: 'copy', saveAs: { it == "virbot/output.vb.fasta" ? "discovered_contigs.fa" : "tax_assignments.tsv" }

    input:
        tuple val(unique_id), val(taxon), path(contigs)
    output:
        path "virbot/output.vb.fasta"
        path "pos_contig_score.tsv"
    script:
    """
    VirBot.py --input "${contigs}" --output virbot
    awk -F',' 'BEGIN {print "contig_id\ttaxonomy"} NR>1{ gsub(/; /, ";", \$4) ; print \$1"\t"\$4 }' virbot/pos_contig_score.csv > pos_contig_score.tsv
    """
}

process download_genomad_database {
    label "process_single"

    container "biocontainers/genomad:1.7.1--pyhdfd78af_0"

    storeDir "${params.store_dir}/"
    output:
        path "genomad_db"
    script:
    """
    if [ -d "${params.store_dir}/genomad_db" ]
    then
        echo "Database already downloaded"
    else
        genomad download-database .
    fi
    """
}

process run_genomad {
    label "process_medium"
    errorStrategy 'ignore'

    container "biocontainers/genomad:1.7.1--pyhdfd78af_0"

    input:
        tuple val(unique_id), val(taxon), path(contigs)
        path genomad_db
    output:
        tuple val(unique_id), val(taxon), path("genomad/filtered_contigs_summary/filtered_contigs_virus.fna"), path("tax_assignments.tsv")
    script:
    """
    genomad end-to-end -t ${task.cpus} --disable-find-proviruses --relaxed \
        ${contigs} genomad ${genomad_db}
    if [ -s 'genomad/filtered_contigs_summary/filtered_contigs_virus_summary.tsv' ];
    then
        grep 'Riboviria' genomad/filtered_contigs_summary/filtered_contigs_virus_summary.tsv | \
        awk -F'\t' 'BEGIN{print "contig_id\ttaxonomy"} NR>1{print \$1"\t"\$11}' > tax_assignments.tsv
    else
        touch tax_assignments.tsv
    fi
    """
}

process filter_short_contigs {
    // preprocess assembly before genomad
    label "process_low"
    errorStrategy 'ignore'

    conda "bioconda::bbmap"
    container "biocontainers/bbmap:39.01--h5c4e2a8_0"
 
    input:
        tuple val(unique_id), val(taxon), path(contigs), path(tax_assignments)
    output:
        tuple val(unique_id), val(taxon), path("filtered_contigs.fa"), path(tax_assignments)
    script:
    """
    reformat.sh in=${contigs} out=filtered_contigs.fa minlength=2000
    """
}

process select_Riboviria {
    // postprocess genomad output selecting RNA viral contigs
    label "process_single"
    errorStrategy 'ignore'

    conda "bioconda::bbmap"
    container "biocontainers/bbmap:39.01--h5c4e2a8_0"

    publishDir "${params.outdir}/${unique_id}/discovery/{taxon}", mode: 'copy', saveAs: { it == "RNA_viral_contigs.fa" ? "discovered_contigs.fa" : "tax_assignments.tsv" }

    input:
        tuple val(unique_id), val(taxon), path(contigs), path(tax_assignments)
    output:
        tuple val(unique_id), val(taxon), path("RNA_viral_contigs.fa"), path(tax_assignments)
    script:
    """
    awk '{print \$1}' ${tax_assignments} > names.txt
    filterbyname.sh in=${contigs} out=RNA_viral_contigs.fa names=names.txt
    """
}

workflow classify_novel_taxa {
    take:
        taxon_fastq_ch
    main:
        assemble_taxa(taxon_fastq_ch)

        if ( params.classifier == 'virbot' ) {
            run_virbot(assemble_taxa.out.contigs)
        } else if ( params.classifier == 'genomad' ) {
            if (! params.genomad_db) {
                genomad_db_ch = download_genomad_database()
            } else {
                Channel.of(file(${params.genomad_db}, type: "dir", checkIfExists:true))
                    .set {genomad_db_ch}
            }
            filter_short_contigs(assemble_taxa.out.contigs)
            run_genomad(filter_short_contigs.out, genomad_db_ch)
            select_Riboviria(run_genomad.out)
        } else {
            error "Invalid classifier: ${params.classifier} - must be either 'virbot' or 'genomad'"
        }
}

workflow {
    if (params.paired){
        fastq = file("${params.fastq1}", type: "file", checkIfExists:true)
        fastq2 = file("${params.fastq2}", type: "file", checkIfExists:true)
    } else {
        fastq = file("${params.fastq}", type: "file", checkIfExists:true)
    }
    unique_id = "${params.unique_id}"
    if ("${params.unique_id}" == "null") {
        unique_id = "${fastq.simpleName}"
    }

    taxon = "all"
    if (params.paired)
        taxon_fastq_ch = Channel.of([unique_id, taxon, fastq, fastq2])
    else
        taxon_fastq_ch = Channel.of([unique_id, taxon, fastq])
    taxon_fastq_ch.view()

    classify_novel_taxa(taxon_fastq_ch)
}

