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
#######indexing using hisat2#####################
#downloaded genome separately and uploaded it to my directory using WinSCP
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

######generating mapped bam files ###################

#!/bin/bash
#SBATCH --array=1-12
#SBATCH --time=20:00:00
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=3
#SBATCH --job-name=mapped_bam_files
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8

WORKDIR="/data/users/rsubramanian/rnaseq/ref_genome"
OUTDIR="$WORKDIR/mapped_reads"
INDEXFILES="$WORKDIR/genome_index"

mkdir -p $OUTDIR

# Define sample and read files based on task ID
case $SLURM_ARRAY_TASK_ID in
1)
  SAMPLE="HER21"
  READ1="/data/courses/rnaseq_course/breastcancer_de/reads/HER21_R1.fastq.gz"
  READ2="/data/courses/rnaseq_course/breastcancer_de/reads/HER21_R2.fastq.gz"
  ;;
2)
  SAMPLE="HER22"
  READ1="/data/courses/rnaseq_course/breastcancer_de/reads/HER22_R1.fastq.gz"
  READ2="/data/courses/rnaseq_course/breastcancer_de/reads/HER22_R2.fastq.gz"
  ;;
3)
  SAMPLE="HER23"
  READ1="/data/courses/rnaseq_course/breastcancer_de/reads/HER23_R1.fastq.gz"
  READ2="/data/courses/rnaseq_course/breastcancer_de/reads/HER23_R2.fastq.gz"
  ;;
4)
  SAMPLE="NonTNBC1"
  READ1="/data/courses/rnaseq_course/breastcancer_de/reads/NonTNBC1_R1.fastq.gz"
  READ2="/data/courses/rnaseq_course/breastcancer_de/reads/NonTNBC1_R2.fastq.gz"
  ;;
5)
  SAMPLE="NonTNBC2"
  READ1="/data/courses/rnaseq_course/breastcancer_de/reads/NonTNBC2_R1.fastq.gz"
  READ2="/data/courses/rnaseq_course/breastcancer_de/reads/NonTNBC2_R2.fastq.gz"
  ;;
6)
  SAMPLE="NonTNBC3"
  READ1="/data/courses/rnaseq_course/breastcancer_de/reads/NonTNBC3_R1.fastq.gz"
  READ2="/data/courses/rnaseq_course/breastcancer_de/reads/NonTNBC3_R2.fastq.gz"
  ;;
7)
  SAMPLE="Normal1"
  READ1="/data/courses/rnaseq_course/breastcancer_de/reads/Normal1_R1.fastq.gz"
  READ2="/data/courses/rnaseq_course/breastcancer_de/reads/Normal1_R2.fastq.gz"
  ;;
8)
  SAMPLE="Normal2"
  READ1="/data/courses/rnaseq_course/breastcancer_de/reads/Normal2_R1.fastq.gz"
  READ2="/data/courses/rnaseq_course/breastcancer_de/reads/Normal2_R2.fastq.gz"
  ;;
9)
  SAMPLE="Normal3"
  READ1="/data/courses/rnaseq_course/breastcancer_de/reads/Normal3_R1.fastq.gz"
  READ2="/data/courses/rnaseq_course/breastcancer_de/reads/Normal3_R2.fastq.gz"
  ;;
10)
  SAMPLE="TNBC1"
  READ1="/data/courses/rnaseq_course/breastcancer_de/reads/TNBC1_R1.fastq.gz"
  READ2="/data/courses/rnaseq_course/breastcancer_de/reads/TNBC1_R2.fastq.gz"
  ;;
11)
  SAMPLE="TNBC2"
  READ1="/data/courses/rnaseq_course/breastcancer_de/reads/TNBC2_R1.fastq.gz"
  READ2="/data/courses/rnaseq_course/breastcancer_de/reads/TNBC2_R2.fastq.gz"
  ;;
12)
  SAMPLE="TNBC3"
  READ1="/data/courses/rnaseq_course/breastcancer_de/reads/TNBC3_R1.fastq.gz"
  READ2="/data/courses/rnaseq_course/breastcancer_de/reads/TNBC3_R2.fastq.gz"
  ;;
*)
  echo "Invalid SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
  exit 1
  ;;
esac

apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif bash -c \
    "hisat2 -p 3 -x $INDEXFILES -1 $READ1 -2 $READ2 -S $OUTDIR/$SAMPLE.sam 2> $OUTDIR/${SAMPLE}_hisat2_summary.log; \
     samtools view -S -b $OUTDIR/$SAMPLE.sam > $OUTDIR/${SAMPLE}_mapped.bam"



######## sorting mapped bam files###################

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


#########feature counts####################
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
