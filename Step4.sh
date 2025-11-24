#!/bin/bash
#SBATCH --time=14:00:00
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=4
#SBATCH --job-name=counts
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8

WORKDIR="/data/users/rsubramanian/rnaseq/ref_genome"
INDIR="$WORKDIR/mapped_reads"
INFILE="$WORKDIR/genome.gtf"
OUTFILE="$WORKDIR/counts.txt"


apptainer exec --bind /data/ /containers/apptainer/subread_2.0.1--hed695b0_0.sif featureCounts -p -T 12 -t exon -g gene_id -a $INFILE -o $OUTFILE $INDIR/*_mapped_sorted.bam
