#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process TRIM_GALORE {

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
    path "*_trimming_report.txt", emit: trimming_reports

    script:
    """
    trim_galore --fastqc --paired --cores ${task.cpus} ${reads}
    """
}
