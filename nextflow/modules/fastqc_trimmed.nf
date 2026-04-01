#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process FASTQC_TRIMMED {

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path "*_fastqc.zip",  emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc --threads ${task.cpus} ${read1} ${read2}
    """
}
