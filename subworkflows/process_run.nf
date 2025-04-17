// workflow to run kraken, check for human, run qc checks and generate html report for a single sample fastq
include { get_params_and_versions } from '../modules/utils'
include { get_default_server; get_viral_server; stop_default_server; stop_viral_server } from '../modules/kraken_server'
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
        input_ch
    main:
        move_or_compress(input_ch).set{ barcode_ch }

        if (params.paired)
            input_ch.map{barcode_id, barcode_files -> [barcode_id, barcode_files[0], barcode_files[1]]}.set{ raw_ch }
        else
            raw_ch = barcode_ch

        classify_and_report(raw_ch, barcode_ch, false)
        extract_all(raw_ch, classify_and_report.out.assignments, classify_and_report.out.kreport, classify_and_report.out.taxonomy)
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
        if (params.paired)
            barcode_input = Channel.fromFilePairs("${run_dir}/*_R{1,2}*.f*q*", type: "file", checkIfExists:true)
        else
            barcode_input = Channel.fromPath("${run_dir}/*", type: "dir", checkIfExists:true, maxDepth:1).map { [it.baseName, get_fq_files_in_dir(it)]}

        if (params.raise_server){
            get_default_server(params.raise_server)
            default_ch = get_default_server.out.server
            if (params.run_viral_reclassification){
                get_viral_server(params.raise_server)
                viral_ch = get_viral_server.out.server
            }
            //default_ch.concat(viral_ch).view()
        }

        process_barcode(barcode_input)
        process_barcode.out.report
            .map{ barcode_id, barcode_report -> "barcode,filepath,sample_report\n${barcode_id},${barcode_id}/classifications/merged.kraken_report.txt,${barcode_id}/${barcode_id}_report.html\n" }.collectFile(name: "${params.outdir}/${unique_id}/samples.csv", sort:true, keepHeader:true, skip:1)

        if (params.raise_server){
            stop_default_server(default_ch, process_barcode.out.report.last())
            if (params.run_viral_reclassification)
                stop_viral_server(viral_ch, process_barcode.out.report.last())
        }
}
