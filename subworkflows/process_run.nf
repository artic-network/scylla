// workflow to run kraken, check for human, run qc checks and generate html report for a single sample fastq
include { get_params_and_versions } from '../modules/get_params_and_versions'
include { kraken_pipeline } from '../subworkflows/kraken_pipeline'
include { extract_reads; extract_paired_reads } from '../modules/extract_taxa'
include { fastp_single; fastp_paired; paired_concatenate } from '../modules/preprocess'
include { classify_novel_taxa; classify_novel_taxa_paired } from '../modules/classify_novel_taxa'

EXTENSIONS = ["fastq", "fastq.gz", "fq", "fq.gz"]

ArrayList get_fq_files_in_dir(Path dir) {
    return EXTENSIONS.collect { file(dir.resolve("*.$it"), type: "file") } .flatten()
}

process move_or_compress {
    cpus 1
    input:
        // don't stage `input` with a literal because we check the file extension
        tuple val(barcode), path(input)
    output:
        val "${barcode}", emit: barcode_id
        path "${barcode}.all.fastq.gz", emit: barcode_fq
    script:
        """
        for file in $input
        do
            if [[ "\$file" == *.gz ]]
            then
                cat "\$file" >> "${barcode}.all.fastq.gz"
            else
                cat "\$file" | bgzip -@ $task.cpus >> "${barcode}.all.fastq.gz"
            fi
        done
        """
}

workflow process_barcode {
    take:
        barcode_fq
    main:
        //file(params.run_dir, type: "dir", checkIfExists:true)
        barcode_id = "${barcode_fq.simpleName}"
        //barcode_dir = barcode. map { get_fq_files_in_dir(it) }.flatten()
        //barcode_dir = file("${params.run_dir}/${barcode}", type: "dir", checkIfExists:true)
        println "Barcode ${barcode} with directory ${barcode_fq}"
        //Channel.fromPath( barcode_dir / "*.f*q*", type: "file")
        //                .set {input_fastq_ch}
        //barcode_fq = barcode_dir. map { get_fq_files_in_dir(it) }
        //barcode_fq = Channel.fromPath( "${barcode_dir}" / "*.f*q*", type: "file")
        //barcode_fq.view()
        //barcode_id = "${barcode_dir.getName()}"
        //barcode_fq = Channel.fromPath( "${params.run_dir}" / "${barcode}" / "*.f*q*", type: "file")
        //barcode_fq.view()
        //fastp_single(barcode_id, barcode_fq)
        //fastp_single.out.processed_fastq.collectFile(name: "${barcode_id}.fq.gz")
        //    .set {processed_fastq}
        //kraken_pipeline(barcode_id, processed_fastq)
        //extract_reads(barcode_id, processed_fastq, kraken_pipeline.out.kraken_assignments, kraken_pipeline.out.bracken_report, kraken_pipeline.out.taxonomy)
    //emit:

}


workflow process_run {
    take:
        unique_id
    main:
        run_dir = file("${params.run_dir}", type: "dir", checkIfExists:true)
        barcode_input = Channel.fromPath("${run_dir}/*", type: "dir", checkIfExists:true, maxDepth:1).map { [it.baseName, get_fq_files_in_dir(it)]}
        ch_input = move_or_compress(barcode_input)
        kraken_pipeline(ch_input.barcode_id, ch_input.barcode_fq)

}
