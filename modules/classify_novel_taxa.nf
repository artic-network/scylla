// This workflow is for trying to find and classify possibly novel RNA viruses that might have been overlooked by kraken. 
// 1) runs a metagenomic assembly on short or long reads,
// 2) runs VirBot or GeNomad classification on them.

// Related options:
// --paired) default: false
// --read_type) 'illumina'|'ont'
// --assembler) for short (Illumina) reads: 'megahit' (default), 'rnaspades'
//              for long (Oxford Nanopore) reads: 'rnabloom' (default), 'flye'
// Notes: megahit is faster while rnaspades tends to produce more contiguous assembly 
//	  megahit here allows short contigs (>=200 bp) whereas rnaspades does not,
//		so in the case of low-quality sequincing data rnaspades might lead to lower recall
//        rnaspades results in the identification of more sequences of retroviral origin
//        flye assembly requires high coverage with reads (at least ~x30) to complete, so it is common to fail
// --classifier) 'virbot'|'genomad', default: both are run
// Notes: genomad requires longer contigs for its neural network to run,
//		so the contigs are filtered by length (>=2kbp) beforehand.
//		Therefore, it is adviced to use genomad on long read data or relatively contiguous assembly
//        overall, at times virbot has greater sensitivity while genomad is more precise, 
//		generally they produce similar output
//        genomad here is run on preset 'relaxed' (no end filtering)
//        
//
// --out_dir) output directory of the pipeline
//
// (for workflow test run:)
// --fastq) path to the file with reads
// --unique_id) sample or run unique id

nextflow.enable.dsl = 2


// main pipeline workflow

workflow classify_novel_taxa {
	if ( params.paired ) {
		take:
		unique_id
		processed_fastq_1
		processed_fastq_2
	} 
	else {
		take:
		unique_id
		fastq
	}



	main:
	if ( params.paired ) {
		assemble_paired(unique_id, processed_fastq_1, processed_fastq_2)
		run_virbot(assemble_paired.out)
	}
        else {
		assemble(unique_id,fastq)
		run_virbot(assemble.out)
        }
	
}

// test run workflow (on test data)
workflow {
	fastq = file("${params.fastq}", type: "file", checkIfExists:true)

	unique_id = "${params.unique_id}"
	if ("${params.unique_id}" == "null") {
		unique_id = "${fastq.simpleName}"
	}
	
	classify_novel_taxa(unique_id,fastq)
}

process assemble {

	input:
	val unique_id
	path fastq

	output:
	path "final_contigs.fa"

	script:
	if ( params.read_type == 'ont' ) {
		if (params.assembler == 'flye') {
			conda "bioconda::flye=2.9"
    			container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        			'https://depot.galaxyproject.org/singularity/flye:2.9--py39h6935b12_1' :
        			'biocontainers/flye:2.9--py39h6935b12_1' }"

			"""

			flye \
        			--nano-raw ${fastq} \
        			-o "assembly"
    			if [ -s assembly/final.contigs.fa ]
    			then
        		mv assembly/final.contigs.fa ${outfile}
        		gzip ${outfile}
    			fi
		
			"""
		} else {
			conda "bioconda::rnabloom"

			script:
			"""
			
			"""
		}
	}
	else if ( params.read_type == 'illumina' ) {
		if ( params.assembler == 'megahit' ) {
			conda "bioconda::megahit"
                	"""
                	megahit -r ${fastq} -m 0.5 --min-contig-len 100 \
				--presets meta-sensitive -t ${task.cpus} -o megahit
                	"""
		} else if ( params.assembler == 'rnaspades' ) {
			"""
			"""
		}	
	}
}

process assemble_paired {

        input:
        val unique_id
        path fastq_1
	path fastq_2

        output:
        path "contigs.fa"

        script:
	if (params.assembler == 'megahit') {
       		conda "bioconda::megahit"
        	"""
		megahit -1 ${fastq_1} -2 ${fastq_2} \
			-m 0.5 --min-contig-len 200 --presets meta-sensitive -t 8 -o megahit
        	if [ -s megahit/final.contigs.fa ]
		then
		mv megahit/final.contigs.fa contigs.fa
		fi
		"""
	} else if (params.assembler == 'rnaspades') {
		// FAULTY CONDA RELEASE
		// USE PARAMS.SPADES_PATH ??

		"""
        	spades.py --rna -1 ${fastq_1} -2 ${fastq_2} -t ${task.cpus} -o rnaspades
        	if [ -s rnaspades/transcripts.fasta ]
                then
                mv rnaspades/transcripts.fasta contigs.fa
                fi
		"""
	}
// else - paired nanopore reads?? 
// else (wrong thing specified) - error
}

process gen_assembly_stats {
	
        conda "bioconda::bbmap"

        //publishDir "${params.out_dir}/${unique_id}/viral_rna_classifications", mode: 'copy', saveAs: { filename -> "${assembler}_${filename}" }

        input:
	val unique_id
	path contigs

        output:
        path("assembly_stats.txt")

        script:
        """
        stats.sh in=${contigs} out="assembly_stats.txt"
        """

process run_virbot {
        // UPLOAD ENV & ADD CONTAINERS
	conda "/localdisk/home/s2420489/conda/virbot.yml"

        //publishDir "${params.out_dir}/${unique_id}/viral_rna_classifications", mode: 'copy'

        input:
        val unique_id
	path contigs

        output:
        path "virbot"

        script:
        """
        VirBot.py --input "${contigs}" --output virbot &
        """
}

process run_genomad {
	
}


process some_example {

	// change dir name
	//publishDir path: "${params.out_dir}/${unique_id}/viral_rna_classifications", mode: 'copy'

	//conda ""
	//container??

	//input:
	
}
