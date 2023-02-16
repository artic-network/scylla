// module to check for human contamination in reads and throw away

process dehumanizeWithMinimap2 {
    cpus 1
    input:
        path fastq
    output:
        path "${fastq.baseName}.clean.fastq"
    script:
        manifest=file({params.manifest})
        preset={params.preset}
    """
    $projectDir/../bin/simple_dehumanize.py \
      --manifest {manifest} \
      --preset {preset} \
      --clean "${fastq.baseName}.clean.fastq" \
      --dirty "${fastq.baseName}.dirty.fastq" \
      {fastq}
    """
}

process dehumanizeWithKraken {
    cpus 1
    input:
        path fastq
        path kraken_assignments
    output:
        path "${fastq.baseName}.clean.fastq"
    script:
    """
    $projectDir/../bin/simple_dehumanize.py \
      --kraken_assignments {kraken_assignments} \
      --clean "${fastq.baseName}.clean.fastq" \
      --dirty "${fastq.baseName}.dirty.fastq" \
      {fastq}
    """
}

workflow dehumanize {
    take:
        fastq
        kraken_assignments
    main:
        dehumanizeWithKraken(fastq, kraken_assignments)
    emit:
        clean = dehumanizeWithKraken.output

}

workflow {
    dehumanize()
}
