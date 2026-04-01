#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process FASTQC {

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.zip",  emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """
}
