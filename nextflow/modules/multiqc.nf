#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path collected_files
    val  output_name

    output:
    path "${output_name}.html", emit: report
    path "${output_name}_data", emit: data

    script:
    """
    multiqc . --filename ${output_name}.html --force
    """
}
