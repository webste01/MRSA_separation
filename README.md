# MRSA_separation
Strain separation of two MRSA strains, using haplotype phasing and epigenetics

Immediate next steps in analysis:


#1. Identify reads that preferentially overlap one contig over another 

Input: Set of overlapping reads from the get_overlapping_reads.py script
	Data necessary for running this script is located in the /data dir (there is also a README in the data dir):
		/sc/orga/work/webste01/gitrepos/MRSA_separation/data/output.ann.bam
		/sc/orga/work/webste01/gitrepos/MRSA_separation/data/branching_tigs_fake.txt
		/sc/orga/work/webste01/gitrepos/MRSA_separation/data/reference.fasta
Output: File with read name and the contig name that is preferentially overlapped by that read

#2. Extract the haplotype for these reads (ZH tag in the data/output.ann.bam)
Input: the read names of the overlapping reads. 
Output: read names of the overlapping reads, the contigs they preferentially map to, and their haplotype

