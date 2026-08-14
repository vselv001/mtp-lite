#!/bin/bash
# ------- PRE-REQUISITES
# Install the following packages in your conda environment before running this script
# hifiasm (channel: bioconda)
# command to install them:
# conda install -c bioconda hifiasm

READ_FILE="path/to/your/input_reads.fasta"  # Replace with the path to your input reads file
OUT_DIR="path/to/your/output_directory"  # Replace with the desired output directory
PREFIX="prefix_for_output_files"  # Replace with the desired prefix for output files
THREADS="number_of_threads"  # Replace with the number of threads you want to use

hifiasm -o $OUT_DIR/$PREFIX  -t $THREADS $READ_FILE
awk '/^S/{print ">"$2;print $3}' $OUT_DIR/$PREFIX.bp.p_ctg.gfa > $OUT_DIR/$PREFIX.asm.fasta
