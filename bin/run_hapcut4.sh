#!/bin/bash
#BSUB -J hapcut4
#BSUB -P acc_PBG
#BSUB -q alloc
#BSUB -n 4
#BSUB -R span[ptile=4]
#BSUB -R rusage[mem=12000]
#BSUB -W 08:00
#BSUB -m manda
#BSUB -o /sc/orga/scratch/webste01/mrsa_analysis/hapcut4.stdout
#BSUB -eo /sc/orga/scratch/webste01/mrsa_analysis/hapcut4.stderr
#BSUB -L /bin/bash


/sc/orga/work/webste01/gitrepos/hapcut/HAPCUT --fragments /sc/orga/scratch/webste01/mrsa_analysis/fragment_matrix_file --VCF /sc/orga/scratch/webste01/mrsa_analysis/scripts/tmp.vcf  --output /sc/orga/scratch/webste01/mrsa_analysis/output_haplotype_file4 --maxiter 100  --longreads 1  > /sc/orga/scratch/webste01/mrsa_analysis/hapcut4.log
