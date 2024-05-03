// module to calculate and filter spike ins

// NB these could be a named set of taxa, or could correspond to a reference file

process map_to_spike_ins {
    label "process_single"

    conda "bowtie2=2.5.3"
    container "biocontainers/bowtie2:2.5.3--py39h6fed5c7_1"

    input:
        tuple val(unique_id), path(fastq)
        path(spike_in_ref_dir)
        val(spike_in_prefixes)
    output:
        tuple val(unique_id), path("spike_ins.sam"), emit: spike
        tuple val(unique_id), path("filtered.fq.gz"), emit: filtered_reads
        tuple val(unique_id), path("filter.log"), emit: log
    script:
        """
        mv ${fastq} start.fq.gz
        for prefix in ${spike_in_prefixes.join(' ')}
        do
            bowtie2 -x ${spike_in_ref_dir}/\$prefix start.fq.gz -S spike_ins_\$prefix.sam --no-unal --un-gz filtered.fq.gz >> filter.log
            mv filtered.fq.gz start.fq.gz
        done
        mv start.fq.gz filtered.fq.gz
        samtools merge spike_ins_*.sam -o spike_ins.sam
        """
}

process map_to_spike_ins_paired {
    label "process_single"

    conda "bowtie2=2.5.3 samtools=1.20"
    container "biocontainers/bowtie2:2.5.3--py39h6fed5c7_1"

    input:
        tuple val(unique_id), path(fastq1), path(fastq2)
        path(spike_in_ref_dir)
        val(spike_in_prefixes)
    output:
        tuple val(unique_id), path("spike_ins.sam"), emit: spike
        tuple val(unique_id), path("filtered_1.fq.gz"), path("filtered_2.fq.gz"), emit: filtered_reads
        tuple val(unique_id), path("filter.log"), emit: log
    script:
        """
        mv ${fastq1} start_1.fq.gz
        mv ${fastq2} start_2.fq.gz
        for prefix in ${spike_in_prefixes.join(' ')}
        do
            bowtie2 -x ${spike_in_ref_dir}/\$prefix -1 start_1.fq.gz -2 start_2.fq.gz -S "spike_ins_\$prefix.sam" --no-unal --un-conc-gz filtered_%.fq.gz >> filter.log
            mv filtered_1.fq.gz start_1.fq.gz
            mv filtered_2.fq.gz start_2.fq.gz
        done
        mv start_1.fq.gz filtered_1.fq.gz
        mv start_2.fq.gz filtered_2.fq.gz
        samtools merge spike_ins_*.sam -o spike_ins.sam
        """
}

workflow get_spike_in_channels {
    main:
        if (params.spike_ins) {
            spike_ins = "${params.spike_ins}".split(',') as List
            Channel.from( spike_ins )
                .branch { spike_in ->
                    valid: spike_in.matches("[0-9.]+")
                        return spike_in
                    ref: params.spike_in_refs.containsKey(spike_in)
                        return params.spike_in_refs.get(spike_in, false)
                    expand: params.spike_in_taxon_sets.containsKey(spike_in)
                        return params.spike_in_taxon_sets.get(spike_in, false)
                    other: true
                        return spike_in
                    }
                    .set { result }

            taxon_keys = params.spike_in_taxon_sets.keySet()
            ref_keys = params.spike_in_refs.keySet()
            result.other.map { spike_in -> throw new Exception("Named spike in $spike_in is invalid, must be one of $taxon_keys or $ref_keys")}

            result.valid.concat(result.expand).flatten().collect().set{ expanded_spike_ins }
            result.ref.flatten().collect().set{ ref_spike_ins }
        } else {
            expanded_spike_ins = Channel.empty()
            ref_spike_ins = Channel.empty()
        }
    emit:
        taxa = expanded_spike_ins
        refs = ref_spike_ins
}

workflow subtract_spike_ins_by_mapping {
    take:
        fastq_ch
        refs_ch
    main:
        ref_dir = file(params.spike_in_ref_dir, type: "dir", checkIfExists:true)
        if (params.paired) {
            println "paired and ref spike"
            map_to_spike_ins_paired(fastq_ch, ref_dir, refs_ch)
            filtered_ch = map_to_spike_ins_paired.out.filtered_reads
        } else {
            println "single and ref spike"
            map_to_spike_ins(fastq_ch, ref_dir, refs_ch)
            filtered_ch = map_to_spike_ins.out.filtered_reads
        }
        result_ch = filtered_ch.ifEmpty(fastq_ch)
    emit:
        fastq=result_ch
}

workflow {
    unique_id = "${params.unique_id}"
    if (params.fastq) {
        fastq = file(params.fastq, type: "file", checkIfExists:true)
        if (unique_id == "null") {
           unique_id = "${fastq.simpleName}"
        }
        fastq_ch = Channel.of([unique_id, fastq])
    } else if (params.fastq1) {
        fastq1 = file(params.fastq1, type: "file", checkIfExists:true)
        fastq2 = file(params.fastq2, type: "file", checkIfExists:true)
        if (unique_id == "null") {
           unique_id = "${fastq1.simpleName}"
        }
        fastq_ch = Channel.of([unique_id, fastq1, fastq2])
    }

    get_spike_in_channels()
    subtract_spike_ins_by_mapping(fastq_ch, get_spike_in_channels.out.refs)
    subtract_spike_ins_by_mapping.out.view()
}

