// This workflow is for trying to find and classify possibly novel RNA viruses that might have been overlooked by kraken. 
// 1) runs a metagenomic assembly on short or long reads,
// 2) runs VirBot or GeNomad classification on them.             
// Output: file with contigs identified as RNA viral 'viral_contigs.fa',
//     contig id to taxonomy mapping file 'tax_assignments.tsv',
//     (optional) {assembler}_assembly_stats.txt file with the performed assembly statistics -
//     all in the 'discovery' directory

// Related options and params:
// --paired) default: false
// --read_type) 'illumina'|'ont'
// --assembler) for short (Illumina) reads: 'megahit' (default), 'rnaspades'
//          for long (Oxford Nanopore) reads: 'rnabloom' (default), 'flye'
// Notes: megahit is faster while rnaspades tends to produce more contiguous assembly 
//      megahit here allows short contigs (>=150 bp) whereas rnaspades does not,
//        so in the case of low-quality sequincing data the assembly with
//        rnaspades might lead to lower recall
//    rnaspades results in the identification of more sequences of retroviral origin
//    flye assembly requires high coverage with reads (at least ~x30) to complete,
//        so it was common to fail, and FAILURE IS NOT HANDLED.
//      rnabloom will fail if there are very few reads, again, FAILURE IS NOT HANDLED.
//        Also it might produce an empty fasta file if all the resulting contigs have
//        length below the threshold of 200 bp. This will fail with genomad but not with virbot 
//    all mentioned assemblers (rnabloom, flye, megahit, spades) should work on gzipped files
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
// --write_assembly_stats) default: true
//
// (for workflow test run:)
// --fastq) path to the file with reads; REQUIRED for testing
// --unique_id) sample or run unique id


// 1. Assemble reads into contigs

process assemble_flye {

    label "process_high"

    conda "bioconda::flye=2.9"
    container "biocontainers/flye:2.9--py39h6935b12_1"

    input:
        val unique_id
        path fastq

    output:
        path "flye/final.contigs.fa"

    script:
    """
    flye --nano-raw ${fastq} --meta -t ${task.cpus} -outdir "flye"
    """
}

process assemble_rnabloom {

    label "process_high"

    conda "bioconda::rnabloom"
    container "${params.wf.container}:${params.wf.container_version}"

    input:
        val unique_id
        path fastq, name: "reads.fastq.gz"

    output:
        path "rnabloom/rnabloom.transcripts.fa"

    script:
    """
    rnabloom -long reads.fastq.gz -t ${task.cpus} -outdir "rnabloom"
    """
}


process assemble_megahit {
    
    label "process_high"

	conda "bioconda::megahit"
	container "biocontainers/mulled-v2-0f92c152b180c7cd39d9b0e6822f8c89ccb59c99:8ec213d21e5d03f9db54898a2baeaf8ec729b447-0"

    input:
        val unique_id
        path fastq

    output:
        path "megahit/final.contigs.fa"

    script:
    """
    megahit -r ${fastq} -m 0.5 --min-contig-len 100 \
            --presets meta-sensitive -t ${task.cpus} -o megahit
    """
}

process assemble_rnaspades {

    label "process_high"
    container "biocontainers/spades:3.15.5--h95f258a_1"

    input:
        val unique_id
        path fastq

    output:
        path "rnaspades/transcripts.fasta"

    script:
    """
    spades.py --rna -s ${fastq} -t ${task.cpus} -o rnaspades
    """
}


process assemble_megahit_paired {
    
    label "process_high"

    conda "bioconda::megahit"
	container "biocontainers/mulled-v2-0f92c152b180c7cd39d9b0e6822f8c89ccb59c99:8ec213d21e5d03f9db54898a2baeaf8ec729b447-0"

    input:
        val unique_id
        path fastq_1
        path fastq_2

    output:
        path "megahit/final.contigs.fa"

    script:
    """
    megahit -1 ${fastq_1} -2 ${fastq_2} \
            -m 0.5 --min-contig-len 150 --presets meta-sensitive -t 8 -o megahit
    """
}

process assemble_rnaspades_paired {

    label "process_high"

    container "biocontainers/spades:3.15.5--h95f258a_1"

    input:
        val unique_id
        path fastq_1
        path fastq_2

    output:
        path "rnaspades/transcripts.fasta"

    script:
    """
    spades.py --rna -1 ${fastq_1} -2 ${fastq_2} -t ${task.cpus} -o rnaspades
    """
}



// 1a. Count assembly statistics

process gen_assembly_stats {

    label "process_single"

    conda "bioconda::bbmap"
    container "biocontainers/bbmap:39.01--h92535d8_1"

    publishDir "${params.outdir}/${unique_id}/discovery", mode: 'copy', saveAs: { filename -> "${params.assembler}_${filename}" }

    input:
        val unique_id
        path contigs

    output:
        path("assembly_stats.txt")

    script:
    """
    stats.sh in=${contigs} out="assembly_stats.txt"
    """
}

// 2. Run identification

process run_virbot {

    label "process_medium"

    container "biowilko/virbot:0.1.0"

    publishDir "${params.outdir}/${unique_id}/discovery", mode: 'copy', saveAs: { it == "virbot/output.vb.fasta" ? "viral_contigs.fa" : "tax_assignments.tsv" }

    input:
        val unique_id
        path contigs

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

    container "biocontainers/genomad:1.7.1--pyhdfd78af_0"

    input:
        val unique_id
        path contigs, name: 'filtered_contigs.fa'
        path genomad_db

    output:
        path "genomad/filtered_contigs_summary/filtered_contigs_virus.fna", emit: contigs
        path "tax_assignments.tsv", emit: taxonomy

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

    conda "bioconda::bbmap"

    container "biocontainers/bbmap:39.01--h5c4e2a8_0"
 
    input:
        val unique_id
        path contigs

    output:
        path 'filtered_contigs.fa'

    script:
    """
    reformat.sh in=${contigs} out=filtered_contigs.fa minlength=2000
    """
}

process select_Riboviria {
    // postprocess genomad output selecting RNA viral contigs

    label "process_single"

    conda "bioconda::bbmap"

    container "biocontainers/bbmap:39.01--h5c4e2a8_0"

    publishDir "${params.outdir}/${unique_id}/discovery", mode: 'copy', saveAs: { it == "RNA_viral_contigs.fa" ? "viral_contigs.fa" : "tax_assignments.tsv" }

    input:
        val unique_id
        path tax_ids
        path viral_contigs
    
    output:
        path tax_ids
        path "RNA_viral_contigs.fa"

    script:
    """
    awk '{print \$1}' ${tax_ids} > names.txt
    filterbyname.sh in=${viral_contigs} out=RNA_viral_contigs.fa names=names.txt
    """
}

workflow classify_novel_taxa {

    take:
        unique_id
        fastq

    main:
        // specify the assembler
        if ( params.read_type == 'ont' ) {
            params.assembler = 'rnabloom'
        } else if ( params.read_type == 'illumina' ) {
            params.assembler = 'megahit'
        } else {
            error "Invalid specification of read_type: ${params.read_type} - must be one of [ont, illumina]"
        }

        // 1. Assembly reads into contigs
        if ( params.read_type == 'ont' ) {
            if ( params.assembler == 'flye' ) {
                assemble_flye(unique_id,fastq)
                    .set{ contigs }
            } else if ( params.assembler == 'rnabloom' ) {
                assemble_rnabloom(unique_id,fastq)
                    .set{ contigs }
            } else {
                error "Invalid assembler specification for ont reads: ${params.assembler} - must be one of [flye, rnabloom]"
            }

        } else if ( params.read_type == 'illumina' ) {
            if ( params.assembler == 'megahit' ) {
                assemble_megahit(unique_id,fastq)
                    .set{ contigs }
            } else if ( params.assembler == 'rnaspades' ) {
                assemble_rnaspades(unique_id,fastq)
                    .set{ contigs }
            } else {
                error "Invalid assembler specification for Illumina reads: ${params.assembler} - must be one of [megahit, rnaspades]"
            }
        } else {
            error "Invalid read_type specification: ${params.read_type} - must be one of [ont, illumina]"
        }

        // 1a. Write assembly statistics
        if ( params.write_assembly_stats ) {
            gen_assembly_stats(unique_id,contigs)
        }

        // 2. Run viral identification
        if ( params.classifier == 'virbot' ) {
            run_virbot(unique_id,contigs)
        } else if ( params.classifier == 'genomad' ) {
            if (! params.genomad_db) {
                download_genomad_database()
                download_genomad_database.out.set {genomad_db_ch}
            } else {
                Channel.of(file(${params.genomad_db}, type: "dir", checkIfExists:true))
                    .set {genomad_db_ch}
            }
                filter_short_contigs(unique_id,contigs)
                run_genomad(unique_id,filter_short_contigs.out, genomad_db_ch)
                select_Riboviria(unique_id,run_genomad.out.taxonomy,run_genomad.out.contigs)
        } else {
            error "Invalid classifier: ${params.classifier} - must be either 'virbot' or 'genomad'"
        }
}

workflow classify_novel_taxa_paired {
    take:
        unique_id
        fastq_1
        fastq_2

    main:
        // 1. Assembly reads into contigs
        if (params.assembler == 'megahit') {
            assemble_megahit_paired(unique_id, fastq_1, fastq_2)
                .set{ contigs }
        } else if (params.assembler == 'rnaspades') {
            assemble_rnaspades_paired(unique_id, fastq_1, fastq_2)
                .set{ contigs }
        } else {
            error "Invalid assembler specification for Illumina reads: ${params.assembler} - must be one of [megahit, rnaspades]"
        }

        // 1a. Write assembly statistics
        if ( params.write_assembly_stats ) {
            gen_assembly_stats(unique_id,contigs)
        }

        // 2. Run viral identification
        if ( params.classifier == 'virbot' ) {
            run_virbot(unique_id,contigs)
        } else if ( params.classifier == 'genomad' ) {
            if (! params.genomad_db) {
                download_genomad_database()
                download_genomad_database.out.set {genomad_db_ch}
            } else {
                Channel.of(file(${params.genomad_db}, type: "dir", checkIfExists:true))
                    .set {genomad_db_ch}
            }
            filter_short_contigs(unique_id,contigs)
            run_genomad(unique_id,filter_short_contigs.out, genomad_db_ch)
            select_Riboviria(unique_id,run_genomad.out.taxonomy,run_genomad.out.contigs)
        } else {
            error "Invalid classifier: ${params.classifier} - must be either 'virbot' or 'genomad'"
        }
}


workflow {
    unique_id = "${params.unique_id}"

    // check input fastq exists
    if (params.fastq) {
        fastq = file("${params.fastq}", type: "file", checkIfExists:true)
        if (! unique_id) {
            unique_id = "${fastq.simpleName}"
        }
        input_fastq = Channel.fromPath(fastq)
    } else if (params.fastq_dir) {
        fastqdir = file("${params.fastq_dir}", type: "dir", checkIfExists:true)
        if (! unique_id) {
            unique_id = "${fastqdir.simpleName}"
        }
        input_fastq = Channel.fromPath( fastqdir / "*.f*q*", type: "file")
    } else {
        exit 1, "One of fastq or fastq_dir need to be provided -- aborting"
    }

    classify_novel_taxa(unique_id,fastq)
}

