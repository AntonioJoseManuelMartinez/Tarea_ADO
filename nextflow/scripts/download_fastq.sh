#!/bin/bash
#SBATCH --job-name=download_fastq
#SBATCH --partition=eck-q
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --chdir=/home/alumno07/rnaseq_nextflow
#SBATCH --output=logs/download_fastq.log

mkdir -p data/fastq_split

RUNS=(SRR26327152 SRR26327153 SRR26327154 SRR26327155 SRR26327156 SRR26327157)

for RUN in "${RUNS[@]}"; do
    echo "Obteniendo URLs para ${RUN}..."
    # Consultar ENA API para obtener las URLs exactas
    URLS=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RUN}&result=read_run&fields=fastq_ftp" \
           | grep -v "^run_accession" \
           | cut -f2)
    
    URL_1=$(echo "$URLS" | tr ';' '\n' | grep "_1.fastq.gz")
    URL_2=$(echo "$URLS" | tr ';' '\n' | grep "_2.fastq.gz")

    echo "  URL_1: $URL_1"
    echo "  URL_2: $URL_2"

    curl "ftp://${URL_1}" -o "data/fastq_split/${RUN}_1.fastq.gz"
    curl "ftp://${URL_2}" -o "data/fastq_split/${RUN}_2.fastq.gz"
    echo "  Listo: ${RUN}"
done

echo "Todas las descargas completadas."
