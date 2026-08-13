#!/bin/bash
# ------- PRE-REQUISITES
# Install the following packages in your conda environment before running this script
# hifiasm (channel: bioconda)
# command to install them:
# conda install -c bioconda hifiasm

READ_FILE="/home/vselv001/MTPLite/output/mtpv1.1.fasta"
OUT_DIR="/home/vselv001/MTPLite/output/mtpv1.1"
PREFIX="mtpv1.1"
THREADS="40"

hifiasm -o $OUT_DIR/$PREFIX  -t $THREADS $READ_FILE
awk '/^S/{print ">"$2;print $3}' $OUT_DIR/$PREFIX.bp.p_ctg.gfa > $OUT_DIR/$PREFIX.asm.fasta
