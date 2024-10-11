#!/usr/bin/env bash
#$ -S /bin/bash  # run job as a Bash shell [IMPORTANT]
#$ -cwd          # run job in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l mem_free=8G     # job requires up to 8 GiB of RAM per slot
#$ -l scratch=100G      # job requires up to 100 GiB of local /scratch space
#$ -l h_rt=24:00:00   # job requires up to 24 hours of runtime
#$ -r y               # if job crashes, it should be restarted

# Load conda and GATK (need to use 4.3.0.0 not current 4.4.0.0 to be compatible with Wynton JRE)
module load CBI miniconda3 gatk/4.3.0.0

#Activate conda jupyter notebook in which to run GATK (need to have already created myjupyter)
conda activate myjupyter

# Define the input directory, where the processed fastq files are located
input_dir="./Analysis/Mapped"

# Define the output directory, where the merged bam files will be saved
output_dir="./Analysis/GATK_Output"

# Iterate over .bam files in the input directory
for Input_bam in "$input_dir"/*.bam; do

    # Extract the file name without extension
    file_wo_extension=$(basename "$Input_bam" ".bam")

    # Define the output file name prefix for GATK analyzed reads
    output="$output_dir/${file_wo_extension}"

# Run GATK AnalyzeSaturationMutagenesis command
gatk AnalyzeSaturationMutagenesis -I "$Input_bam" -R POP2_seq/POP2_seq.fa --orf 151-1224 -O "$output"

    echo "Analyzed ${file_wo_extension}"
done

echo "AnalyzeSaturationMutagenesis complete"

# deactivate conda
conda deactivate

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                          # e.g. "did my job exceed its memory request?"
