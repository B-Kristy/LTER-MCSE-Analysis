#!/bin/bash

cd /mnt/home/kristybr/20230410_Amplicon_KRI13538_PE250/merged
# Filter merged fastq sequences
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_filter merged_reads.fastq -fastq_trunclen 380 -fastq_maxee 1.0 -fastaout merged_reads_filtered.fastq
