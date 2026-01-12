#!/bin/bash
#SBATCH --job-name=nf-PLNT                     # Job name
#SBATCH --mem=8G                               # Job memory request
#SBATCH --cpus-per-task=1                      # number of cpu per task
#SBATCH --time=96:00:00

export PYTHONPATH=$PYTHONPATH:$(pwd)/src
pytest tests

module load nextflow

nextflow run main.nf \
    -profile slurm \
    --genome data/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
    --annotation data/annotation.csv