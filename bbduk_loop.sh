#!/bin/bas#! /usr/bin/env bash
#$ -S /bin/bash  # run job as a Bash shell [IMPORTANT]
#$ -cwd          # run job in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l mem_free=2G     # job requires up to 2 GiB of RAM per slot
#$ -l scratch=100G      # job requires up to 100 GiB of local /scratch space
#$ -l h_rt=24:00:00   # job requires up to 24 hours of runtime
#$ -r y               # if job crashes, it should be restarted

# Load conda and activate myjupyter
module load CBI miniconda3
conda activate myjupyter

# Define the input directory where your .fastq files are located
input_dir="./TLJY04"

# Define the output directory where processed files will be saved
output_dir="./Analysis/PostBBDuk"

# Iterate over .fastq files in the input directory
for R1 in "$input_dir"/*_R1_001.fastq.gz; do
    # Extract the file name without extension
    file_wo_extension=$(basename "$R1" "_R1_001.fastq.gz")

    # Define the corresponding reverse read file
    R2="$input_dir/${file_wo_extension}_R2_001.fastq.gz"

    # Define the output file names for processed reads
    output_F="$output_dir/${file_wo_extension}_processed_R1.fastq.gz"
    output_R="$output_dir/${file_wo_extension}_processed_R2.fastq.gz"

    # Run bbduk.sh command
    bbduk.sh in1="$R1" in2="$R2" out1="$output_F" out2="$output_R" ref="/wynton/home/ye/laflamt/.conda/envs/myjupyter/opt/bbmap-38.18/resources/adapters.fa" ktrim=r k=23 mink=11 tpe tbo

    echo "Processed $file_wo_extension"
done

echo "Processing complete"

# deactivate conda
conda deactivate 

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                          # e.g. "did my job exceed its memory request?"
