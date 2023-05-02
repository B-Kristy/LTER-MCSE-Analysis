#!/bin/bash

module load Conda/3

# Corrected AA sequences are in _corr_prot.fasta# Screen with hmm to identify all hitscd /mnt/home/kristybr/20230410_Amplicon_KRI13538_PE250hmmscan --cpu 16  --domtblout hmmOut2.out /mnt/home/kristybr/NifMAP/Resources/nifH_chlL_bchX.hmm ./nifH/corr_prot.fasta &> /dev/null