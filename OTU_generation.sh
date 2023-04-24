#!/bin/bash
# First, dereplicate the sequences
cd /mnt/home/kristybr/20230410_Amplicon_KRI13538_PE250

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_uniques ./merged/merged_reads_filtered_hmm.fa -minuniquesize 2 -fastaout ./merged/unique.fa 

# Second, sort unique sequences by length
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -sortbylength ./merged/unique.fa -fastout ./merged/unique_sorted.fa

# Third, cluster sequences into OTUs at 97% similarity
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -cluster_fast ./merged/unique_sorted.fa -id 0.97 -centroids otus.fa -uc clusters.uc

