#!/bin/bash
#SBATCH --job-name=hisat2_build
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=04:00:00
#SBATCH --output=logs/hisat2_build.log
#SBATCH --chdir=/home/alumno07/rnaseq_nextflow
#SBATCH --partition=eck-q

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq_env

GENOME="data/genome/genome.fa"
GTF="data/genome/annotation.gtf"
IDXDIR="data/index"
mkdir -p "$IDXDIR"

echo "Extrayendo splice sites y exones del GTF..."
hisat2_extract_splice_sites.py "$GTF" > "${IDXDIR}/splice_sites.txt"
hisat2_extract_exons.py "$GTF"        > "${IDXDIR}/exons.txt"

echo "Construyendo índice HISAT2 con 32 hilos..."
hisat2-build \
    -p 32 \
    --ss "${IDXDIR}/splice_sites.txt" \
    --exon "${IDXDIR}/exons.txt" \
    "$GENOME" \
    "${IDXDIR}/GRCh38_109_index"

echo "Comprimiendo índice para Nextflow..."
tar -czvf data/GRCh38_109_index.tar.gz -C "$IDXDIR" \
    $(ls ${IDXDIR}/GRCh38_109_index*.ht2 | xargs -n1 basename)

echo "Índice disponible en data/GRCh38_109_index.tar.gz"
ls -lh data/GRCh38_109_index.tar.gz
