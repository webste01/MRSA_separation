#1. run blasr to align bams back to merged file
/sc/orga/projects/bashia02b/tools/bin/blasr_for_hapcut -nproc 12 all.fofn  /sc/orga/scratch/webste01/mrsa/mrsa_only/sequence/mrsa_only.fasta -sam -out all.sam -bestn 1 -clipping subread

#2. Convert to bam
all.bam

#3. Sort
all.sorted.bam 

#4. pass into fake hairs note the -r flag defines a region (this is a nice feature for subsetting the vcf/bams), vcf file must not have a header
#python fake_hairs_snps_ali.py -r unitig_61:1-150000 -i all_test.sorted.bam -v unitig_61_input_wtabs_noheader.vcf -o all_test.hairs

#5. Extract the hairs WORKS ON THE NEW BAMS, DONT NEED FAKE HAIRS

~/ThirdParty/HAPCUT-v0.7/extractHAIRS  --mbq 0  --maxIS 1000000 --noquality 7 --qvoffset 33 --VCF /sc/orga/scratch/webste01/mrsa_analysis/scripts/tmp.vcf --bam all.sorted.bam

/sc/orga/work/webste01/gitrepos/hapcut/extractHAIRS  --mbq 0  --maxIS 1000000 --noquality 7 --qvoffset 33 --VCF /sc/orga/scratch/webste01/mrsa_analysis/scripts/tmp.vcf --bam all.sorted.bam > fragment_matrix_file


#6. Run hapcut
/sc/orga/work/webste01/gitrepos/hapcut/HAPCUT --fragments fragment_matrix_file --VCF /sc/orga/scratch/webste01/mrsa_analysis/scripts/tmp.vcf  --output output_haplotype_file --maxiter maxiter  > hapcut.log

