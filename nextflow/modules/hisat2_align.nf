#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process HISAT2_ALIGN {
    input:
    tuple val(sample_id), path(read1), path(read2)
    path index_files   // lista de .ht2 copiados al workdir

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam
    path "${sample_id}.hisat2.log",                        emit: log

    script:
    // El prefijo es el nombre sin .1.ht2
    def idx = index_files[0].toString().replaceAll(/\.\d+\.ht2$/, '')
    """
    hisat2 -x ${idx} \\
           -1 ${read1} -2 ${read2} \\
           -p ${task.cpus} \\
           --new-summary --summary-file ${sample_id}.hisat2.log | \\
        samtools view -bS | \\
        samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam
    samtools index ${sample_id}.sorted.bam
    """
}
