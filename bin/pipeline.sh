#Run strain phasing pipeline from beginning () to end ()

#1. Align reads back to merged assembly with blasr
/sc/orga/projects/bashia02b/tools/bin/blasr_for_hapcut -nproc 12 all.fofn  /sc/orga/scratch/webste01/mrsa/mrsa_only/sequence/mrsa_only.fasta -sam -out all.sam -bestn 20 -clipping subread


#2. Sort amd convert to bam
samtools view -Sb all.sam  > all.bam
samtools sort all.bam -o sorted.bam


#3. Make vcf file
#VCF header info lines have to be structured this way:
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples with Fully Called Data">

python ../bin/syntheticStrainVcf.py 020180.variants.vcf 020042.variants.vcf --outfile mrsa_separation.vcf -s


#3. Extract the hairs from the bam file

/sc/orga/work/webste01/gitrepos/hapcut/extractHAIRS  --mbq 0  --maxIS 1000000 --noquality 7 --qvoffset 33 --VCF mrsa_separation.vcf --bam sorted.bam > fragment_matrix_file


#6. Run hapcut
/sc/orga/work/webste01/gitrepos/hapcut/HAPCUT --fragments fragment_matrix_file --VCF /sc/orga/scratch/webste01/mrsa_analysis/scripts/tmp.vcf  --output output_haplotype_file --maxiter maxiter  > hapcut.log

#7. Run read phaser
python /sc/orga/work/webste01/gitrepos/readphaser/pr.py -u unphased.fq -p phased.fq output_haplotype_file4 all.sorted.bam 

#8. Run greedy partitioner
python /sc/orga/work/webste01/gitrepos/ryan-neff-bashir-lab-rotation/scripts/greedy_partitioner.py  -h myhairs.hairs  -c output_haplotype_file4  -i all.sorted.bam -o output.ann.bam
