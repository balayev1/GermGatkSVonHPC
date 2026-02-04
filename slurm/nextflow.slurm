#!/bin/bash
#SBATCH --job-name=GATK_SV_Master
#SBATCH --output=logs/nf_gatksv_%j.out
#SBATCH --error=logs/nf_gatksv_%j.err
#SBATCH --partition=asvnode1,msibigmem
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=balay011@umn.edu

# Make logs directory
mkdir -p logs

# Load required modules
module load nextflow
module load java/openjdk-17.0.2
module load singularity
module load R

# Singularity cachedir
export NXF_SINGULARITY_CACHEDIR="/scratch.global/balay011/singularity_cache"
mkdir -p $NXF_SINGULARITY_CACHEDIR

# Run the pipeline
nextflow run main.nf \
    -profile singularity \
    -resume