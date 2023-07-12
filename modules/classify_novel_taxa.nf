// This workflow is for trying to find and classify possibly novel RNA viruses that might have been overlooked by kraken. 
// 1) runs a metagenomic assembly on short or long reads,
// 2) runs VirBot or GeNomad classification on them.
// The output currently is the "discovery" folder with a "virbot" and/or "genomad" output directory and 
//			     {assembler}_assembly_stats.txt file with the performed assembly statistics
// Output: file with contigs identified as RNA viral 'viral_contigs.fa'
//         contig id to taxonomy mapping file 'tax_ide', in the 'discovery' directory

// Related options and params:
// --paired) default: false
// --read_type) 'illumina'|'ont'
// --assembler) for short (Illumina) reads: 'megahit' (default), 'rnaspades'
//              for long (Oxford Nanopore) reads: 'rnabloom' (default), 'flye'
// Notes: megahit is faster while rnaspades tends to produce more contiguous assembly 
//	  megahit here allows short contigs (>=150 bp) whereas rnaspades does not,
//		so in the case of low-quality sequincing data the assembly with
//		rnaspades might lead to lower recall
//        rnaspades results in the identification of more sequences of retroviral origin
//        flye assembly requires high coverage with reads (at least ~x30) to complete, 
//		so it was common to fail, and FAILURE IS NOT HANDLED.
//	  rnabloom will fail if there are very few reads, again, FAILURE IS NOT HANDLED.
//		Also it might produce an empty fasta file if all the resulting contigs have
//		length below the threshold of 200 bp. This will fail with genomad but not with virbot 
// --classifier) 'virbot'|'genomad', default: virbot
// Notes: genomad requires longer contigs for its neural network to run,
//		so the contigs are filtered by length (>=2kbp) beforehand.
//		Therefore, it is adviced to use genomad on long read data or relatively contiguous assembly
//        overall, at times virbot has greater sensitivity while genomad might be more precise (?), 
//		but generally they produce similar output
//        virbot does classification to the lowest level possible
//        genomad here is run on preset 'relaxed' (no end filtering)
//        for genomad, the assignments other than 'Riboviria' (including those labeled as 'Unclassified') 
//		are excluded from the final taxonomy and contigs files
// --out_dir) output directory of the pipeline; REQUIRED for testing
// --genomad_db) path to the directory with genomad database
// --write_assembly_stats) default: false
//
// (for workflow test run:)
// --fastq) path to the file with reads; REQUIRED for testing
// --unique_id) sample or run unique id

nextflow.enable.dsl = 2


// set parameters

params.fastq = "/localdisk/home/s2420489/mscape_test/data/CAMB_28_6_23_B2_SISPArapid/CAMB_28_6_23_B2_SISPArapid.fq"
params.unique_id = "test1_ont"
params.out_dir = "/localdisk/home/s2420489/nextflow/mscape/test_pipeline"

//params.classifier = "virbot"
params.genomad_db = "/localdisk/home/s2420489/software/genomad_db"

params.paired = false
params.read_type = 'ont'
if ( params.read_type == 'ont' ) {
	params.assembler = 'rnabloom'
} else if ( read_type == 'illumina' ) {
	params.assembler = 'megahit'
} else {
	error "Invalid specification of read_type: ${params.read_type} - must be one of [ont, illumina]"
}


// main pipeline workflow

workflow classify_novel_taxa_single {

	take:
	unique_id
	fastq
		
	main:
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


	gen_assembly_stats(unique_id,contigs)

	if ( params.classifier == 'virbot' ) {
		run_virbot(unique_id,contigs)
	} else if ( params.classifier == 'genomad' ) {
                filter_short_contigs(unique_id,contigs)
                run_genomad(unique_id,filter_short_contigs.out)
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
	if (params.assembler == 'megahit') {
		assemble_megahit_paired(unique_id, fastq_1, fastq_2)
			.set{ contigs }
	} else if (params.assembler == 'rnaspades') {
		assemble_rnaspades_paired(unique_id, fastq_1, fastq_2)
			.set{ contigs }
	} else {
                error "Invalid assembler specification for Illumina reads: ${params.assembler} - must be one of [megahit, rnaspades]"
        }
	

	if ( params.write_assembly_stats ) {
        	gen_assembly_stats(unique_id,contigs)
	}

        if ( params.classifier == 'virbot' ) {
                run_virbot(unique_id,contigs)
        } else if ( params.classifier == 'genomad' ) {
		filter_short_contigs(unique_id,contigs)
                run_genomad(unique_id,filter_short_contigs.out)
		select_Riboviria(unique_id,run_genomad.out.taxonomy,run_genomad.out.contigs)
        } else {
                error "Invalid classifier: ${params.classifier} - must be either 'virbot' or 'genomad'"
        }
}

// test run workflow (on test data)

workflow {
	fastq = file("${params.fastq}", type: "file", checkIfExists:true)
	println(fastq)

	unique_id = "${params.unique_id}"
	if ("${params.unique_id}" == "null") {
		unique_id = "${fastq.simpleName}"
	}
	
	classify_novel_taxa_single(unique_id,fastq)
}



// processes

// 1. Assemble reads into contigs

process assemble_flye {

	label "process_high"

	conda "bioconda::flye=2.9"
	container "biowilko/scylla@${params.wf.container_sha}"

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
	container "biowilko/scylla@${params.wf.container_sha}"

        input:
        val unique_id
        path fastq, name: 'reads.fastq'

	output:
        path "rnabloom/rnabloom.transcripts.fa"

	script:
	"""
        rnabloom -long reads.fastq -t ${task.cpus} -outdir "rnabloom"
	"""
}


process assemble_megahit {
	
        label "process_high"

	conda "bioconda::megahit"
	container "biowilko/scylla@${params.wf.container_sha}"

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

	//----------ENV----------
	container "biowilko/scylla@${params.wf.container_sha}"

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

// paired input

process assemble_megahit_paired {
	
	label "process_high"

	conda "bioconda::megahit"
	container "biowilko/scylla@${params.wf.container_sha}"

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

	//----------ENV-----------
	container "biowilko/scylla@${params.wf.container_sha}"

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
	container "biowilko/scylla@${params.wf.container_sha}"

        publishDir "${params.out_dir}/${unique_id}/discovery", mode: 'copy', saveAs: { filename -> "${params.assembler}_${filename}" }

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

        // UPLOAD ENV & ADD CONTAINERS
	conda "/localdisk/home/s2420489/conda/virbot.yml"
	container "biowilko/scylla@${params.wf.container_sha}"

        publishDir "${params.out_dir}/${unique_id}/discovery", mode: 'copy', saveAs: { it == "output.vb.fasta" ? "viral_contigs.fa" : "tax_assignments.tsv" }

        input:
        val unique_id
	path contigs

        output:
        path "virbot/output.vb.fasta"
	path "virbot/pos_contig_score.tsv"

        script:
        """
        VirBot.py --input "${contigs}" --output virbot
	awk -F',' 'BEGIN {print "contig_id\ttaxonomy"} NR>1{ gsub(/; /, ";", $4) ; print $1"\t"$4 }' pos_contig_score.csv > pos_contig_score.tsv
        """
}

process run_genomad {

	label "process_medium"	

	//UPLOAD ENV
	conda "/localdisk/home/s2420489/conda/genomad.yml"
	container "biowilko/scylla@${params.wf.container_sha}"

	input:
        val unique_id
        path contigs

        output:
        path "genomad/*.transcripts_summary/*.transcripts_virus.fna", emit: contigs
	path "tax_assignments.tsv", emit: taxonomy

        script:
        """
        genomad end-to-end -t ${task.cpus} --disable-find-proviruses --relaxed \
		${contigs} genomad ${params.genomad_db}
	grep 'Riboviria' genomad/*.transcripts_summary/*.transcripts_virus_summary.tsv | \
		awk -F'\t' 'BEGIN{print "contig_id\ttaxonomy"} NR>1{print $1"\t"$11}' > tax_assignments.tsv
        """
}

// preprocess assembly before genomad
process filter_short_contigs {

	label "process_low"

        conda "bioconda::bbmap"
	container "biowilko/scylla@${params.wf.container_sha}"

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

// postprocess genomad output selecting RNA viral contigs
process select_Riboviria {

	label "process_single"

	conda "bioconda::bbmap"
	container "biowilko/scylla@${params.wf.container_sha}"

	publishDir "${params.out_dir}/${unique_id}/discovery", mode: 'copy', saveAs: { it == "RNA_viral_contigs.fa" ? "viral_contigs.fa" : "tax_assignments.tsv" }

	input:
	val unique_id
	path tax_ids
	path viral_contigs
	
	output:
	path tax_ids
	path "RNA_viral_contigs.fa"

	script:
	"""
	awk '{print $1}' ${tax_ids} > names.txt
	filterbyname.sh in=${viral_contigs} out=RNA_viral_contigs.fa names=names.txt
	"""
}
