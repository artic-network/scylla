// workflow to run kraken, check for human, run qc checks and generate html report for a single sample fastq
include { get_params_and_versions } from '../modules/get_params_and_versions'
include { kraken_setup; kraken_end } from '../modules/kraken_classification'
include { classify_and_report } from '../subworkflows/classify_and_report'
include { classify_virus_fastq } from '../modules/classify_novel_viruses'
include { extract_all } from '../modules/extract_all'

EXTENSIONS = ["fastq", "fastq.gz", "fq", "fq.gz"]

ArrayList get_fq_files_in_dir(Path dir) {
    return EXTENSIONS.collect { file(dir.resolve("*.$it"), type: "file") } .flatten()
}

process move_or_compress {

    label "process_low"

    errorStrategy {task.exitStatus == 2 ? 'ignore' : 'terminate'}

    conda "bioconda::tabix==v1.11"
    container "${params.wf.container}:${params.wf.container_version}"
    cpus 1

    input:
        // don't stage `input` with a literal because we check the file extension
        tuple val(barcode), path(input)
    output:
        tuple val(barcode), path("${barcode}.all.fastq.gz")
    script:
        // nb we rezip with bgzip because usually gzipped and this isn't suitable later
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
        if [ ! -f "${barcode}.all.fastq.gz" ]
        then
            echo "No fastq files"
            exit 2
        fi
        """
}

workflow process_barcode {
    take:
        barcode_ch
    main:
        classify_and_report(barcode_ch, barcode_ch, null)
        extract_all(barcode_ch, classify_and_report.out.assignments, classify_and_report.out.kreport, classify_and_report.out.taxonomy)
        if ( params.classify_novel_viruses ){
            classify_virus_fastq(extract_all.out.virus)
        }
    emit:
        report = classify_and_report.out.report
}

workflow process_run {
    take:
        unique_id
    main:
        get_params_and_versions(unique_id)

        run_dir = file("${params.run_dir}", type: "dir", checkIfExists:true)
        barcode_input = Channel.fromPath("${run_dir}/barcode*", type: "dir", checkIfExists:true, maxDepth:1).map { [it.baseName, get_fq_files_in_dir(it)]}
        move_or_compress(barcode_input)

        if (params.raise_server)
            kraken_setup(params.raise_server)

        process_barcode(move_or_compress.out)
        process_barcode.out.report
            .map{ barcode_id, barcode_report -> "barcode,filepath,sample_report\n${barcode_id},${barcode_id}/classifications/${params.database_set}.kraken_report.txt,${barcode_id}/${barcode_id}_report.html\n" }.collectFile(name: "${params.outdir}/${unique_id}/samples.csv", sort:true, keepHeader:true, skip:1)

        if (params.raise_server)
            kraken_end(kraken_setup.out.server, process_barcode.out.report.collect())
}
