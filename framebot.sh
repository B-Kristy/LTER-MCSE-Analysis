#!/bin/bash
cd /mnt/home/kristybr/20230410_Amplicon_KRI13538_PE250

/mnt/home/kristybr/anaconda3/pkgs/rdptools-2.0.2-1/bin/FrameBot framebot -N -l 30 -i 0.4 -o ./nifH/ /mnt/home/kristybr/Framebot/refset/nifh_prot_ref.fasta otus.fa 
