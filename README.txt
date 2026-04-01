# RNA-seq en sangre completa: análisis de expresión diferencial en Galaxy y Nextflow

Trabajo práctico de la asignatura **Análisis de Datos Ómicos** del Máster en
Bioinformática de la Universidad de Murcia (curso 2025/2026).

**Autores:** Antonio José Manuel Martínez · Álvaro Máximo González López

---

## Descripción

Análisis de expresión diferencial en RNA-seq de sangre completa comparando individuos
control frente a pacientes portadores de variantes en *ZCCHC8*, un componente del
complejo NEXT implicado en la maduración del ARN de la telomerasa. Los datos proceden
del repositorio GEO bajo el identificador
[GSE244915](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244915).

El trabajo se desarrolló en dos plataformas complementarias:

- **Galaxy Europe** — flujo completo: FastQC → STAR → featureCounts → DESeq2.
  El análisis diferencial y el enriquecimiento funcional se realizaron a partir
  de los ficheros de conteos y resultados exportados desde Galaxy.
- **Nextflow DSL2** — pipeline modular ejecutado en el clúster Dayhoff (SLURM),
  que reproduce las etapas de preprocesamiento y cuantificación incorporando
  además recorte de adaptadores y un informe agregado con MultiQC.

---

## Datos

| Condición | Accessions SRA |
|-----------|---------------|
| Control   | SRR26327155, SRR26327156, SRR26327157 |
| ZCCHC8    | SRR26327152, SRR26327153, SRR26327154 |

- Plataforma: Illumina NovaSeq 6000, *paired-end* 150 pb
- Material: sangre completa en tubos PAXgene, ARN total
- Librería: NEBNext Ultra II Directional RNA Library Prep Kit

---

## Estructura del repositorio
```
.
├── Codigo_Tarea_ADO.R               # Análisis diferencial, volcano plot y GO en R
├── ficheros_generados_Galaxy/       # Salidas exportadas desde Galaxy Europe
│   ├── Galaxy57-[SRR26327152].tabular
│   ├── Galaxy58-[SRR26327153].tabular
│   ├── Galaxy59-[SRR26327154].tabular
│   ├── Galaxy60-[SRR26327155].tabular
│   ├── Galaxy61-[SRR26327156].tabular
│   ├── Galaxy62-[SRR26327157].tabular
│   └── Galaxy66-[DESeq2 result file on dataset 57-62].tabular
├── imagenes_generadas/              # Figuras del análisis
│   ├── Galaxy67-[DESeq2 plots on dataset 57-62].pdf
│   ├── GO_BP_dotplot_up_control.png
│   ├── GO_BP_dotplot_up_ZCCHC8.png
│   └── volcano_plot_control_vs_ZCCHC8.png
└── nextflow/                        # Pipeline alternativo en Nextflow DSL2
    ├── main.nf
    ├── nextflow.config
    ├── modules/
    │   ├── fastqc.nf
    │   ├── fastqc_trimmed.nf
    │   ├── trim_galore.nf
    │   ├── hisat2_align.nf
    │   ├── featurecounts.nf
    │   └── multiqc.nf
    ├── scripts/                     # Jobs SLURM para preparación de datos
    │   ├── download_fastq.sh
    │   ├── download_genome.sh
    │   ├── build_index.sh
    │   └── run_nextflow.sh
    └── results/
        ├── counts/                  # Ficheros de conteos por muestra
        │   ├── control_1.counts.txt
        │   ├── control_2.counts.txt
        │   ├── control_3.counts.txt
        │   ├── zcchc8_1.counts.txt
        │   ├── zcchc8_2.counts.txt
        │   └── zcchc8_3.counts.txt
        └── multiqc/
            └── rnaseq_report.html   # Informe agregado de calidad
```

---

## Análisis en Galaxy

El flujo de trabajo se ejecutó en [Galaxy Europe](https://usegalaxy.eu/) e incluyó:

1. **FastQC** — inspección de calidad de las librerías crudas
2. **RNA STAR** — alineamiento al genoma humano hg38
3. **samtools flagstat** — evaluación básica del alineamiento
4. **featureCounts** — cuantificación de lecturas por gen
5. **DESeq2** — modelado estadístico de expresión diferencial (control vs. ZCCHC8)

Los ficheros exportados desde Galaxy se encuentran en `ficheros_generados_Galaxy/`
y constituyen la entrada del script de R.

---

## Análisis en R

El script `Codigo_Tarea_ADO.R` toma como entrada la tabla de resultados de DESeq2
exportada desde Galaxy y produce:

- Anotación de genes con **org.Hs.eg.db**
- Tablas de genes diferencialmente expresados (padj < 0.05, |log2FC| > 1)
- **Volcano plot** (control vs. ZCCHC8)
- Análisis de enriquecimiento GO Biological Process con **clusterProfiler**
- Dotplots de términos GO enriquecidos en cada grupo

Dependencias de R: `ggplot2`, `dplyr`, `clusterProfiler`, `org.Hs.eg.db`,
`enrichplot`, `AnnotationDbi`.

---

## Pipeline de Nextflow
```
FASTQC (raw) → TRIM_GALORE → FASTQC (trimmed) → HISAT2 → SAMTOOLS_SORT → FEATURECOUNTS → MULTIQC
```

### Herramientas y versiones

| Herramienta   | Versión | Función                     |
|---------------|---------|-----------------------------|
| FastQC        | 0.12.1  | Control de calidad          |
| Trim Galore   | 0.6.10  | Recorte de adaptadores      |
| HISAT2        | 2.2.2   | Alineamiento al genoma      |
| SAMtools      | 1.9     | Conversión y ordenación BAM |
| featureCounts | 2.0.6   | Cuantificación génica       |
| MultiQC       | 1.9     | Informe agregado de QC      |
| Nextflow      | 24.10.6 | Orquestación del pipeline   |

Referencia genómica: **GRCh38, Ensembl 109**

### Reproducción

**1. Crear el entorno Conda**
```bash
conda create -n rnaseq_env -c bioconda -c conda-forge \
    fastqc trim-galore hisat2 samtools subread multiqc nextflow -y
```

**2. Descargar los datos** (desde el directorio `nextflow/`)
```bash
mkdir -p logs
sbatch scripts/download_fastq.sh    # FASTQs desde ENA vía FTP
sbatch scripts/download_genome.sh   # Genoma GRCh38 y GTF desde Ensembl 109
```

**3. Construir el índice HISAT2**
```bash
sbatch scripts/build_index.sh       # ~2-4 h con 32 hilos, requiere ~180 GB RAM
```

**4. Lanzar el pipeline**
```bash
sbatch scripts/run_nextflow.sh
```

Los resultados se depositan en `results/`. La opción `-resume` permite reanudar
ejecuciones interrumpidas reutilizando los resultados en caché.

---

## Referencias principales

- Tummala *et al.* (2024). *Blood*, 143(10), 890–910.
- Kim *et al.* (2019). HISAT2. *Nature Biotechnology*, 37, 907–915.
- Liao *et al.* (2014). featureCounts. *Bioinformatics*, 30(7), 923–930.
- Love *et al.* (2014). DESeq2. *Genome Biology*, 15, 550.
- Di Tommaso *et al.* (2017). Nextflow. *Nature Biotechnology*, 35, 316–319.
- Ewels *et al.* (2016). MultiQC. *Bioinformatics*, 32(19), 3047–3048.
- Yu *et al.* (2012). clusterProfiler. *OMICS*, 16(5), 284–287.

> **Nota:** el directorio `data/` (FASTQs, genoma de referencia e índice HISAT2) no está incluido en el repositorio por su tamaño. 
> Puede generarse siguiendo los pasos de la sección anterior.
