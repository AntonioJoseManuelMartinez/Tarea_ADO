#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process FEATURECOUNTS {
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path gtf

    output:
    path "${sample_id}.counts.txt",         emit: counts
    path "${sample_id}.counts.txt.summary", emit: summary

    script:
    """
    featureCounts \\
        -T ${task.cpus} \\
        -p \\
        -a ${gtf} \\
        -o ${sample_id}.counts.txt \\
        ${bam}
    """
}
