#!/bin/bash
#SBATCH --array=1-12               # Array for 12 samples
#SBATCH --time=20:00:00           # Maximum runtime
#SBATCH --mem=16gb                # Memory allocation
#SBATCH --cpus-per-task=3         # Number of CPU cores per task
#SBATCH --job-name=sort_index_bam # Job name
#SBATCH --output=array_%J.out     # Standard output log
#SBATCH --error=array_%J.err      # Error log
#SBATCH --partition=pibu_el8      # Partition

WORKDIR="/data/users/rsubramanian/rnaseq/ref_genome/mapped_reads"
OUTDIR="$WORKDIR"

# Define sample names based on SLURM_ARRAY_TASK_ID
case $SLURM_ARRAY_TASK_ID in
1)
  SAMPLE="HER21"
  ;;
2)
  SAMPLE="HER22"
  ;;
3)
  SAMPLE="HER23"
  ;;
4)
  SAMPLE="NonTNBC1"
  ;;
5)
  SAMPLE="NonTNBC2"
  ;;
6)
  SAMPLE="NonTNBC3"
  ;;
7)
  SAMPLE="Normal1"
  ;;
8)
  SAMPLE="Normal2"
  ;;
9)
  SAMPLE="Normal3"
  ;;
10)
  SAMPLE="TNBC1"
  ;;
11)
  SAMPLE="TNBC2"
  ;;
12)
  SAMPLE="TNBC3"
  ;;
*)
  echo "Invalid SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
  exit 1
  ;;
esac

# Navigate to output directory
cd $OUTDIR

# Perform sorting and indexing for the specific sample
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif bash -c \
    "samtools sort -@ 3 $OUTDIR/${SAMPLE}_mapped.bam -o $OUTDIR/${SAMPLE}_mapped_sorted.bam; \
     samtools index $OUTDIR/${SAMPLE}_mapped_sorted.bam"

# Check for success and provide feedback
if [ $? -eq 0 ]; then
  echo "Sorting and indexing completed successfully for $SAMPLE"
else
  echo "Error: Sorting and indexing failed for $SAMPLE"
  exit 1
fi
