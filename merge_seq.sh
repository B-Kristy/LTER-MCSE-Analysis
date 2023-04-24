#!/bin/bash

cd /mnt/home/kristybr/20230410_Amplicon_KRI13538_PE250
# Merge paired-end .fastq sequences using usearch
for i in *R1* .fastq; do /mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_mergepairs $i -fastqout merged_reads_${i} -sample $i
done

# Concatenate all merged reads into one file
cat *merged_reads* > /merged/merged_reads.fa

# Remove intermediate merge filescd /mnt/home/kristybr/20230410_Amplicon_KRI13538_PE250/rm *merged_reads*