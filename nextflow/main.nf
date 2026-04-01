#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FASTQC         } from './modules/fastqc.nf'
include { TRIM_GALORE    } from './modules/trim_galore.nf'
include { FASTQC_TRIMMED } from './modules/fastqc_trimmed.nf'
include { HISAT2_ALIGN   } from './modules/hisat2_align.nf'
include { FEATURECOUNTS  } from './modules/featurecounts.nf'
include { MULTIQC        } from './modules/multiqc.nf'

workflow {
    read_ch = Channel.fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.sample_id, [file(row.fastq_1), file(row.fastq_2)])
        }

    // Índice: recoger todos los .ht2 como lista
    hisat2_idx_ch = Channel.fromPath("${params.hisat2_index}*.ht2").collect()

    // 1. QC inicial
    FASTQC(read_ch)

    // 2. Recorte de adaptadores
    TRIM_GALORE(read_ch)

    // 3. QC post-recorte
    FASTQC_TRIMMED(TRIM_GALORE.out.trimmed_reads)

    // 4. Alineamiento + ordenar BAM (integrado)
    HISAT2_ALIGN(
        TRIM_GALORE.out.trimmed_reads,
        hisat2_idx_ch
    )

    // 5. Cuantificación
    FEATURECOUNTS(
        HISAT2_ALIGN.out.bam,
        file(params.gtf)
    )

    // 6. Informe MultiQC agregado
    multiqc_input = FASTQC.out.zip
        .mix(FASTQC.out.html)
        .mix(TRIM_GALORE.out.trimming_reports)
        .mix(FASTQC_TRIMMED.out.zip)
        .mix(FASTQC_TRIMMED.out.html)
        .mix(HISAT2_ALIGN.out.log)
        .mix(FEATURECOUNTS.out.summary)
        .collect()

    MULTIQC(multiqc_input, params.report_id)
}
