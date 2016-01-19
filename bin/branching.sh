#Get reads that map to the 1st and last 20kb of unitig 54
fastalength /sc/orga/scratch/webste01/mrsa/mrsa_only/sequence/mrsa_only.fasta | grep "unitig_54"
#20kb = 20000 base pairs
samtools view output.ann.bam unitig_54:343618-363618 > last_20kb_unitig54.bam
samtools view output.ann.bam unitig_54:0-20000 > first_20kb_unitig54.bam

fastalength /sc/orga/scratch/webste01/mrsa/mrsa_only/sequence/mrsa_only.fasta | grep "unitig_61"
#20kb = 20000 base pairs
samtools view -b output.ann.bam unitig_61:101841-121841  > last_20kb_unitig61.bam
samtools view -b output.ann.bam unitig_61:0-20000 > first_20kb_unitig61.bam

#Get read IDS
samtools view -F 4 first_20kb_unitig61.bam | cut -f1 | sort -u > first_utg61


fastalength /sc/orga/scratch/webste01/mrsa/mrsa_only/sequence/mrsa_only.fasta | grep "unitig_36"
#20kb = 20000 base pairs
samtools view output.ann.bam unitig_36:26518-46518 > last_20kb_unitig36.bam
samtools view output.ann.bam unitig_36:0-20000 > first_20kb_unitig36.bam


#Re-header the new sam files
#Save the header to a file:
samtools view -H output.ann.bam > all.sorted.header

#Add the header to the sam file
samtools reheader all.sam last_20kb_unitig61.bam -o last_20kb_unitig61_reheader.bam
samtools reheader all.sam last_20kb_unitig61.bam -o first_20kb_unitig61_reheader.bam
