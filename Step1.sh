#fastqc
#!/bin/bash
#SBATCH --job-name=fastqc_job                # Job name
#SBATCH --output=fastqc_%j.out               # Output file
#SBATCH --error=fastqc_%j.err                # Error file
#SBATCH --partition=pshort_el8               # Specify partition
#SBATCH --time=01:00:00                      # Time limit: hh:mm:ss
#SBATCH --cpus-per-task=1                    # Number of CPU cores
#SBATCH --mem=1000MB                         # Memory per node

# Run FastQC on all .fastq.gz files
singularity exec /containers/apptainer/fastqc-0.12.1.sif fastqc \
  -o /data/users/rsubramanian/rnaseq/fastqc_results \
  -t 1 \
  /data/courses/rnaseq_course/breastcancer_de/reads/*.fastq.gz

##downloaded genome separately and uploaded it to my directory using WinSCP
#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=16
#SBATCH --job-name=index_files
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8

WORKDIR="/data/users/rsubramanian/rnaseq/ref_genome"
REFDIR="$WORKDIR"
INDEXDIR="$REFDIR/genome_index"

# Create required directories

mkdir -p $INDEXDIR


# Perform HISAT2 indexing
cd $REFDIR
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
hisat2-build -p 16 genome.fa $INDEXDIR/genome_index

echo "Indexing completed successfully!"