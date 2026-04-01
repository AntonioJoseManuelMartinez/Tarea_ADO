#!/bin/bash
#SBATCH --job-name=download_genome
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/download_genome.log
#SBATCH --chdir=/home/alumno07/rnaseq_nextflow
#SBATCH --partition=eck-q

mkdir -p data/genome logs
cd data/genome

echo "Descargando genoma GRCh38 primary assembly..."
curl -L "https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
     -o genome.fa.gz &

echo "Descargando anotación GTF Ensembl 109..."
curl -L "https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz" \
     -o annotation.gtf.gz &

wait
echo "Descomprimiendo..."
gunzip genome.fa.gz
gunzip annotation.gtf.gz

echo "Listo. Archivos disponibles en data/genome/:"
ls -lh
