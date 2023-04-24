#!/bin/bash

module load Conda/3
cd /mnt/home/kristybr/20230410_Amplicon_KRI13538_PE250

# Screen merged, filtered reads using HMMER
hmmsearch --cpu 32 --domtblout hmmOut1.out /mnt/home/kristybr/NifMAP/Resources/hmm_nuc_1160_nifH.hmm ./merged/merged_reads_filtered.fa &/dev/null