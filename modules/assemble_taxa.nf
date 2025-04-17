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
// --write_assembly_stats) default: true

process assemble_flye {
    label "process_high"
    errorStrategy 'ignore'

    conda "bioconda::flye=2.9"
    container "biocontainers/flye:2.9--py39h6935b12_1"

    publishDir "${params.outdir}/${unique_id}/assemblies/${taxon}", mode: 'copy', pattern: "flye/assembly.fasta", saveAs: {filename -> "flye_assembled_contigs.fa"}
    input:
        tuple val(unique_id), val(taxon), path(fastq)
    output:
        tuple val(unique_id), val(taxon), path("flye/assembly.fasta")
    script:
    """
    flye --nano-raw ${fastq} --meta -t ${task.cpus} --out-dir "flye"
    """
}

process assemble_rnabloom {
    label "process_high"
    errorStrategy 'ignore'

    conda "bioconda::rnabloom"
    container "docker.io/jdelling7igfl/rnabloom:2.0.1"

    publishDir "${params.outdir}/${unique_id}/assemblies/${taxon}", mode: 'copy', pattern: "rnabloom/rnabloom.transcripts.fa", saveAs: {filename -> "rnabloom_assembled_contigs.fa"}
    input:
        tuple val(unique_id), val(taxon), path(fastq)
    output:
        tuple val(unique_id), val(taxon), path("rnabloom/rnabloom.transcripts.fa")
    script:
    """
    rnabloom -long ${fastq} -t ${task.cpus} -outdir "rnabloom"
    """
}

process assemble_megahit {
    label "process_high"
    errorStrategy 'ignore'

	conda "bioconda::megahit conda-forge::bzip2 conda-forge::libcxx=8.0"
	container "biocontainers/mulled-v2-0f92c152b180c7cd39d9b0e6822f8c89ccb59c99:8ec213d21e5d03f9db54898a2baeaf8ec729b447-0"

    publishDir "${params.outdir}/${unique_id}/assemblies/${taxon}", mode: 'copy', pattern: "megahit/final.contigs.fa", saveAs: {filename -> "megahit_assembled_contigs.fa"}
    input:
        tuple val(unique_id), val(taxon), path(fastq)
    output:
        tuple val(unique_id), val(taxon), path("megahit/final.contigs.fa")
    script:
    """
    megahit -r ${fastq} -m 0.5 --min-contig-len 100 \
            --presets meta-sensitive -t ${task.cpus} -o "megahit"
    """
}

process assemble_rnaspades {
    label "process_high"
    errorStrategy 'ignore'

    conda "bioconda::spades=3.15 conda-forge::pyyaml==3.12"
    container "biocontainers/spades:3.15.5--h95f258a_1"

    publishDir "${params.outdir}/${unique_id}/assemblies/${taxon}", mode: 'copy', pattern: "rnaspades/soft_filtered_transcripts.fasta", saveAs: {filename -> "rnaspades_assembled_contigs.fa"}

    input:
        tuple val(unique_id), val(taxon), path(fastq)
    output:
        tuple val(unique_id), val(taxon), path("rnaspades/soft_filtered_transcripts.fasta")
    script:
    """
    spades.py --rna -s ${fastq} -t ${task.cpus} -o "rnaspades"
    """
}

process assemble_megahit_paired {
    label "process_high"
    errorStrategy 'ignore'

    conda "bioconda::megahit conda-forge::bzip2 conda-forge::libcxx=8.0"
	container "biocontainers/mulled-v2-0f92c152b180c7cd39d9b0e6822f8c89ccb59c99:8ec213d21e5d03f9db54898a2baeaf8ec729b447-0"

    publishDir "${params.outdir}/${unique_id}/assemblies/${taxon}", mode: 'copy', pattern: "megahit/final.contigs.fa", saveAs: {filename -> "megahit_assembled_contigs.fa"}

    input:
        tuple val(unique_id), val(taxon), path(fastq_1), path(fastq_2)
    output:
        tuple val(unique_id), val(taxon), path("megahit/final.contigs.fa")
    script:
    """
    megahit -1 ${fastq_1} -2 ${fastq_2} \
            -m 0.5 --min-contig-len 150 --presets meta-sensitive -t ${task.cpus} -o megahit
    """
}

process assemble_rnaspades_paired {
    label "process_high"
    errorStrategy 'ignore'

    conda "bioconda::spades=3.15 conda-forge::pyyaml==3.12"
    container "biocontainers/spades:3.15.5--h95f258a_1"

    publishDir "${params.outdir}/${unique_id}/assemblies/${taxon}", mode: 'copy', pattern: "rnaspades/soft_filtered_transcripts.fasta", saveAs: {filename -> "rnaspades_assembled_contigs.fa"}

    input:
        tuple val(unique_id), val(taxon), path(fastq_1), path(fastq_2)
    output:
        tuple val(unique_id), val(taxon), path("rnaspades/soft_filtered_transcripts.fasta")
    script:
    """
    spades.py --rna -1 ${fastq_1} -2 ${fastq_2} -t ${task.cpus} -o rnaspades
    """
}

process generate_assembly_stats {
    label "process_single"
    errorStrategy 'ignore'

    conda "bioconda::bbmap"
    container "biocontainers/bbmap:39.01--h92535d8_1"

    publishDir "${params.outdir}/${unique_id}/assemblies/${taxon}", mode: 'copy', pattern: "assembly_stats.txt"
    input:
        tuple val(unique_id), val(taxon), path(contigs)
    output:
        tuple val(unique_id), val(taxon), path("assembly_stats.txt")
    script:
    """
    stats.sh in=${contigs} out="assembly_stats.txt"
    """
}

workflow assemble_taxa {
    take:
        taxon_fastq_ch
    main:
        // Set/check assembler choice
        assembler = "${params.assembler}"
        read_type = "${params.read_type}"

        if (!params.assembler){
            if ( params.read_type == 'ont' ) {
                assembler = 'rnabloom'
            } else if ( params.read_type == 'illumina' || params.paired) {
                assembler = 'megahit'
                read_type = 'illumina'
            } else {
                error "Invalid specification of read_type: ${params.read_type} - must be one of [ont, illumina]"
            }
        } else {
            if (read_type == 'ont' && assembler != "flye" && assembler != "rnabloom") {
                error "Invalid assembler specification for ont reads: ${assembler} - must be one of [flye, rnabloom]"
            }
            if (read_type == 'illumina' && assembler != "megahit" && assembler != "rnaspades") {
                error "Invalid assembler specification for illumina reads: ${assembler} - must be one of [megahit, rnaspades]"
            }
        }
        println "${read_type}, ${params.paired}, ${assembler}"

        // Assemble reads
        if ( assembler == 'flye' ) {
            assemble_flye(taxon_fastq_ch).set{ contigs }
        } else if ( assembler == 'rnabloom' ) {
            assemble_rnabloom(taxon_fastq_ch).set{ contigs }
        } else if ( assembler == 'megahit' && params.paired) {
            assemble_megahit_paired(taxon_fastq_ch).set{ contigs }
        } else if ( assembler == 'megahit' && !params.paired) {
            assemble_megahit(taxon_fastq_ch).set{ contigs }
        } else if ( assembler == 'rnaspades' && params.paired) {
            assemble_rnaspades_paired(taxon_fastq_ch).set{ contigs }
        } else if ( assembler == 'rnaspades' && !params.paired) {
            assemble_rnaspades(taxon_fastq_ch).set{ contigs }
        } else {
            contigs = Channel.empty()
        }

        // Compose stats
        if ( params.write_assembly_stats ) {
            generate_assembly_stats(contigs)
        }
    emit:
        contigs = contigs
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

    assemble_taxa(taxon_fastq_ch)
}