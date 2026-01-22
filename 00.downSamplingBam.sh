#!/bin/bash

# Function to display help message
show_help() {
cat << EOF
Usage: ${0##*/} [-h] <input_bam> <output_bam> <target_depth>

This script down-samples a BAM file to a target sequencing depth and prints logs to stdout.

    <input_bam>        Input BAM file.
    <output_bam>       Output BAM file.
    <target_depth>     Target sequencing depth (e.g., 5, 10, 20).

Options:
    -h, --help         Display this help and exit.

Example:
    ${0##*/} input.bam output.bam 10
    This will down-sample the input BAM to 10X depth and output the result as output.bam.
    Log information will be printed to stdout and can be redirected to a file if needed.
EOF
}

# Check if help is requested
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_help
    exit 0
fi

# Check if the number of input arguments is correct
if [ "$#" -ne 3 ]; then
  echo "Error: Incorrect number of arguments."
  show_help
  exit 1
fi

# Assign input parameters
input_bam=$1
output_bam=$2
target_depth=$3

# Check if the input BAM file exists
if [ ! -f "$input_bam" ]; then
  echo "Error: Input BAM file does not exist: $input_bam"
  exit 1
fi

# Calculate the current average depth of the input BAM file
current_depth=$(samtools depth "$input_bam" | awk '{sum+=$3} END {print sum/NR}')
echo "Current average depth: $current_depth"

# Exit if the current depth is zero
if (( $(echo "$current_depth == 0" | bc -l) )); then
  echo "Error: The current depth is 0, please check the input BAM file."
  exit 1
fi

# Calculate the sampling ratio based on target depth
sampling_ratio=$(echo "$target_depth / $current_depth" | bc -l)
echo "Sampling ratio: $sampling_ratio"

# Exit if the sampling ratio is greater than 1 (can't increase depth)
if (( $(echo "$sampling_ratio > 1" | bc -l) )); then
  echo "Error: Target depth is greater than current depth. Cannot increase depth."
  exit 1
fi

# Generate the down-sampled BAM file using the calculated sampling ratio
samtools view -s "$sampling_ratio" -b "$input_bam" > "$output_bam"

# Create an index for the new down-sampled BAM file
samtools index "$output_bam"

# Calculate the new average depth for the down-sampled BAM file
new_depth=$(samtools depth "$output_bam" | awk '{sum+=$3} END {print sum/NR}')
echo "New average depth for ${target_depth}x BAM: $new_depth"

# Completion message
echo "Downsampling complete. Output BAM: $output_bam"
