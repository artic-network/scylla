// This workflow is for trying to find and classify possibly novel RNA viruses that might have been overlooked by kraken. 
// 1) runs a metagenomic assembly on short or long reads,
// 2) runs VirBot or GeNomad classification on them.
// The output currently is the "discovery" folder with a "virbot" and/or "genomad" output directory and 
//			     {assembler}_assembly_stats.txt file with the performed assembly statistics
// Desired output: file with viral contigs and simple tsv mapping, ideally the same for virbot and genomad

// Related options and params:
// --paired) default: false
// --read_type) 'illumina'|'ont'
// --assembler) for short (Illumina) reads: 'megahit' (default), 'rnaspades'
//              for long (Oxford Nanopore) reads: 'rnabloom' (default), 'flye'
// Notes: megahit is faster while rnaspades tends to produce more contiguous assembly 
//	  megahit here allows short contigs (>=150 bp) whereas rnaspades does not,
//		so in the case of low-quality sequincing data rnaspades might lead to lower recall
//        rnaspades results in the identification of more sequences of retroviral origin
//        flye assembly requires high coverage with reads (at least ~x30) to complete, 
//		so it is common to fail, and FAILURE IS NOT HANDLED.
// --classifier) 'virbot'|'genomad', default: both are run
// Notes: genomad requires longer contigs for its neural network to run,
//		so the contigs are filtered by length (>=2kbp) beforehand.
//		Therefore, it is adviced to use genomad on long read data or relatively contiguous assembly
//        overall, at times virbot has greater sensitivity while genomad is more precise, 
//		generally they produce similar output
//        virbot does classification to the lowest level possible
//        genomad here is run on preset 'relaxed' (no end filtering)
// --out_dir) output directory of the pipeline; REQUIRED for testing
// --genomad_db) path to the directory with genomad database
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
		run_genomad(unique_id,contigs)	
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
	
        gen_assembly_stats(unique_id,contigs)

        if ( params.classifier == 'virbot' ) {
                run_virbot(unique_id,contigs)
        } else if ( params.classifier == 'genomad' ) {
                run_genomad(unique_id,contigs)
        } else {
                error "Invalid classifier: ${params.classifier} - must be either 'virbot' or 'genomad'"
        }
}

// test run workflow (on test data)

workflow {
	fastq = file("${params.fastq}", type: "file", checkIfExists:true)
	println(fastq)

	test_conda()

	unique_id = "${params.unique_id}"
	if ("${params.unique_id}" == "null") {
		unique_id = "${fastq.simpleName}"
	}
	
	classify_novel_taxa_single(unique_id,fastq)
}



// processes

// 1. Assemble reads into contigs

process assemble_flye {

	conda "bioconda::flye=2.9"
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        	'https://depot.galaxyproject.org/singularity/flye:2.9--py39h6935b12_1' :
                'biocontainers/flye:2.9--py39h6935b12_1' }"

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

	conda "bioconda::rnabloom"

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
	
	conda "bioconda::megahit"

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

	//----------ENV----------

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
	
	conda "bioconda::megahit"

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

	//----------ENV-----------

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

        conda "bioconda::bbmap"

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
        // UPLOAD ENV & ADD CONTAINERS
	conda "/localdisk/home/s2420489/conda/virbot.yml"

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

	conda "/localdisk/home/s2420489/conda/genomad.yml"
	publishDir "${params.out_dir}/${unique_id}/discovery", mode: 'copy', saveAs: { it == "*.transcripts_virus.fna" ? "viral_contigs.fa" : "tax_assignments.tsv" }

	input:
        val unique_id
        path contigs

        output:
        path "genomad/*.transcripts_summary/*.transcripts_virus.fna"
	path "tax_assignments.tsv"

        script:
        """
        genomad end-to-end -t ${task.cpus} --disable-find-proviruses --relaxed \
		${contigs} genomad ${params.genomad_db}
	grep 'Riboviria' genomad/*.transcripts_summary/*.transcripts_virus_summary.tsv | \
		awk 'BEGIN{print "contig_id\ttaxonomy"} NR>1{print $1"\t"$11}' > tax_assignments.tsv
        """
}

// preprocess before genomad
process filter_short_contigs {
        conda "bioconda::bbmap"

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
