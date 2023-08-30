// workflow to run kraken, check for human, run qc checks and generate html report for a single sample fastq
include { get_params_and_versions } from '../modules/get_params_and_versions'
include { kraken_setup; kraken_pipeline; kraken_end } from '../subworkflows/kraken_pipeline'
include { extract_reads } from '../modules/extract_taxa'

EXTENSIONS = ["fastq", "fastq.gz", "fq", "fq.gz"]

ArrayList get_fq_files_in_dir(Path dir) {
    return EXTENSIONS.collect { file(dir.resolve("*.$it"), type: "file") } .flatten()
}

process move_or_compress {

    label "process_low"

    conda "bioconda::tabix==v1.11"
    container "biocontainers/tabix:v1.9-11-deb_cv1"
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
                cat "\$file" | gunzip | bgzip -@ $task.cpus >> "${barcode}.all.fastq.gz"
            else
                cat "\$file" | bgzip -@ $task.cpus >> "${barcode}.all.fastq.gz"
            fi
        done
        """
}

workflow process_barcode {
    take:
        barcode_id
        barcode_fq
    main:
        kraken_pipeline(barcode_id, barcode_fq, null)
        extract_reads(kraken_pipeline.out.unique_id, kraken_pipeline.out.fastq, kraken_pipeline.out.kraken_assignments, kraken_pipeline.out.bracken_report, kraken_pipeline.out.taxonomy)
    emit:
        barcode_report = kraken_pipeline.out.report
        barcode_id = barcode_id
}

workflow process_run {
    take:
        unique_id
    main:
        run_dir = file("${params.run_dir}", type: "dir", checkIfExists:true)
        barcode_input = Channel.fromPath("${run_dir}/*", type: "dir", checkIfExists:true, maxDepth:1).map { [it.baseName, get_fq_files_in_dir(it)]}
        ch_input = move_or_compress(barcode_input)

        if (params.raise_server)
            kraken_setup(params.raise_server)

        process_barcode(ch_input.barcode_id, ch_input.barcode_fq)
        process_barcode.out.barcode_id.collectFile(name: "${params.outdir}/${unique_id}/samples.csv", sort:true, newLine:true) { item -> "${item},${params.outdir}/${item}/${item}_report.html"}

        if (params.raise_server)
            kraken_end(kraken_setup.out.server, process_barcode.out.barcode_report.collect())
}
