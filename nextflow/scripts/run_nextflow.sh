#!/bin/bash
#SBATCH --job-name=rnaseq_nextflow
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/nextflow_main.log
#SBATCH --chdir=/home/alumno07/rnaseq_nextflow
#SBATCH --partition=eck-q

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq_env

echo "Lanzando pipeline Nextflow..."
nextflow run main.nf \
    -resume \
    -with-conda \
    -with-dag dag.png \
    > logs/nextflow_run.log 2>&1

echo "Pipeline finalizado."
