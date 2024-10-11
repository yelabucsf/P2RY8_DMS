#!/usr/bin/env bash
#$ -S /bin/bash  # run job as a Bash shell [IMPORTANT]
#$ -cwd          # run job in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l mem_free=4G     # job requires up to 4 GiB of RAM per slot
#$ -l scratch=100G      # job requires up to 100 GiB of local /scratch space
#$ -l h_rt=24:00:00   # job requires up to 24 hours of runtime
#$ -r y               # if job crashes, it should be restarted

# Load conda and activate myjupyter
module load CBI miniconda3
conda activate myjupyter

# Load samtools
module load CBI samtools

# Designed to be run from within a directory that contains the necessary reference fasta file and index
# Define the input directory, where the processed fastq files are located
input_dir="../Analysis/PostBBMerge"

# Define the output directory, where the merged bam files will be saved
output_dir="../Analysis/Mapped"

# Iterate over .fastq files in the input directory
for R1 in "$input_dir"/*_R1.fastq.gz; do

    # Extract the file name without extension
    file_wo_extension=$(basename "$R1" "_corrected_R1.fastq.gz")

    # Define the corresponding reverse read file
    R2="$input_dir/${file_wo_extension}_corrected_R2.fastq.gz"

    # Define the output file name for mapped reads
    output="$output_dir/${file_wo_extension}_mapped.bam"
    
# Run bbmap.sh command 
bbmap.sh in1="$R1" in2="$R2" outm="$output" sam=1.3

    echo "Mapped ${file_wo_extension}"
done

echo "Mapping complete"

# deactivate conda
conda deactivate

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                          # e.g. "did my job exceed its memory request?"
