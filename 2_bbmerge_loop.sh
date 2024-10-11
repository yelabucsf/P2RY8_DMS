#!/usr/bin/env bash
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

# Define the input directory, where the processed fastq files are located
input_dir="./Analysis/PostBBDuk"

# Define the output directory, where the merged bam files will be saved
output_dir="./Analysis/PostBBMerge"

# Iterate over .fastq files in the input directory
for R1 in "$input_dir"/*_R1.fastq.gz; do

    # Extract the file name without extension
    file_wo_extension=$(basename "$R1" "_processed_R1.fastq.gz")

    # Define the corresponding reverse read file
    R2="$input_dir/${file_wo_extension}_processed_R2.fastq.gz"

    # Define the output file name for corrected reads
    output_F="$output_dir/${file_wo_extension}_corrected_R1.fastq.gz"
    output_R="$output_dir/${file_wo_extension}_corrected_R2.fastq.gz"

# Run bbmerge.sh command. This will error-correct overlapping portion. 
bbmerge.sh in1="$R1" in2="$R2" out1="$output_F" out2="$output_R" ecco mix

    echo "Corrected ${file_wo_extension}"
done

echo "BBMerge correction complete"

# deactivate conda
conda deactivate

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                          # e.g. "did my job exceed its memory request?"

