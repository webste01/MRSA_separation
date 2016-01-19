import subprocess
import pysam
import sys
import os
import shelve

kb=20000
branch_file = "/sc/orga/scratch/webste01/mrsa_analysis/branching_tigs_fake.txt"
ann_bam = pysam.AlignmentFile("/sc/orga/scratch/webste01/mrsa_analysis/output.ann.bam", "rb")
ref = pysam.FastaFile("/sc/orga/scratch/webste01/mrsa/mrsa_only/sequence/mrsa_only.fasta")
ref_fas = "/sc/orga/scratch/webste01/mrsa/mrsa_only/sequence/mrsa_only.fasta"

#Get the branching unitigs from the branch file
def make_branching_dict(branch_file):
	'''
	make a dictionary keyed on branching unitigs, with values of neighboring unitigs from the celera assembler graph.
	'''
	branching_unitigs = {}
	with open(branch_file,"r") as branch:
    		a=branch.readline()
    		for l in branch:
			ll=l.split()[0]
			pred1,pred2,succ1,succ2 = l.split()[2:6]
			branching_unitigs[str(ll)]= pred1,pred2,succ1,succ2
	branch.close()
	return branching_unitigs

branching_tigs = make_branching_dict(branch_file)

#Split up the bam into two different bams according to haplotypes (ZH tag 1/2)
#os.system("samtools view output.ann.bam | grep 'ZH:Z:9278,1' | samtools view -bS -T /sc/orga/scratch/webste01/mrsa/mrsa_only/sequence/mrsa_only.fasta - > ann.hap1.bam")
#os.system("samtools view output.ann.bam | grep 'ZH:Z:9278,2' | samtools view -bS -T /sc/orga/scratch/webste01/mrsa/mrsa_only/sequence/mrsa_only.fasta - > ann.hap2.bam")

#index these bam files (necessary for pysam
#os.system("samtools index /sc/orga/scratch/webste01/mrsa_analysis/ann.hap1.bam")
#os.system("samtools index /sc/orga/scratch/webste01/mrsa_analysis/ann.hap2.bam")
hap1bam = pysam.AlignmentFile("/sc/orga/scratch/webste01/mrsa_analysis/ann.hap1.bam","rb")
hap2bam = pysam.AlignmentFile("/sc/orga/scratch/webste01/mrsa_analysis/ann.hap2.bam","rb")

haplo_dict={}
def make_hap_dict(annotated_bam,hap):
	'''
	make dictionary of read ids and the haplotypes they correspond to
	'''
	for bamread in annotated_bam.fetch():
		haplo_dict[bamread.query_name]=hap
	return haplo_dict

make_hap_dict(hap1bam,1)
make_hap_dict(hap2bam,2)

def get_overhang_reads(unitig, start_pos, end_pos):
	'''
	Gets the overhanging reads of a particular unitig, given a start and end position
	'''
	fasta = open("%s.fasta" %(unitig), "w")
	iter = ann_bam.fetch(unitig,start_pos, end_pos)
	ahs = [] #aligned hits
	for x in iter:
		if x.is_secondary:
			continue
		ahs.append(x)
		overhanging_reads_list.append(x.query_name)
		fasta.write(">%s\n%s\n" %(x.query_name, x.query_sequence))
	fasta.close()
	#Index the fasta
	os.system("samtools faidx %s.fasta" %(unitig))
	return ahs


#get overhanging reads for each unitig in the branching unitigs, create dictionary dict[unitig] = list of branching reads
ahs_start_reads={}
ahs_end_reads={}
overhanging_reads_list=[]
for unitig in branching_tigs:
	ahs_start_reads[unitig] = get_overhang_reads(unitig,0,20000)
	fastalength = subprocess.Popen(["fastalength %s | grep %s" %(ref_fas, unitig)],stdout=subprocess.PIPE,shell=True) 
	a = fastalength.communicate()
	end = int(a[0].split()[0])
	overlap = end - kb
	ahs_end_reads[unitig] = get_overhang_reads(unitig,overlap,end)


#Align branching unitigs back to the reference
#for unitig in branching_tigs.keys():
#	proc = subprocess.Popen(["blasr -bestn 20  %s.fasta %s -m 4 > %s_blasr_out" %(unitig,ref_fas,unitig)],stdout=subprocess.PIPE, shell=True)
#	proc.communicate()


#Make sure you load the Overlap.py from the mrsa_analysis directory
from Overlap import overlap
from Overlap import read_overlaps

c=0
blasr_file = "/sc/orga/scratch/webste01/mrsa_analysis/unitig_54_blasr_out"
with open(blasr_file,"r") as bf:
        for l in bf:
                name1, name2, score, pctiden, strand1, start1, end1, len1, strand2, start2, end2, len2 = l.split()[0:12]
		#check to see that the read is in an overhanging read
		if name1 in overhanging_reads_list:
			print "in list"
		o = overlap(name1, name2, score, pctiden, strand1, start1, end1, len1, strand2, start2, end2, len2)
		if o.hasFullOverlap():
			print name1, "overlapping"
		#else:
		#	print name1, "not full overlap"	
		#else:
	#		print "non overhanging"


#Check whether the overhanging reads are also overlapping reads





#Check whether the branching reads align to neighbors of the unitig in the graph

n=0
m=0
t=0
pos=0
neg=0
for unitig in branching_tigs:
	reads=[]
	blasr_f = open("%s_blasr_out" %(unitig),"r")
	for l in blasr_f:
		ll=l.split()
		t=t+1
		qual = float(ll[3])
		#format: qName tName score percentSimilarity qStrand qStart qEnd qLength tStrand tStart tEnd tLength mapQV
		if qual > 75: #Only use reads that are above a 75% similarity
			n = n+1
			if ll[8] == "0": #Template strand is positive
				pos=pos + 1
				if ll[10] < ll[11]:
					m=m+1
					reads.append(ll[0])
			elif ll[8] == "1": #Template strand is negative
				neg = neg + 1
		overhanging_reads[unitig]=reads


#See if these reads preferentially map to one contig or another, by looking at the blasr output (score from lowest to highest)
rank_file = open("rank_file_2.txt", "w")
for unitig in overhanging_reads:
		for rname in overhanging_reads[unitig]:
			proc = subprocess.Popen(["grep %s %s_blasr_out" %(rname,unitig)],stdout=subprocess.PIPE,shell=True)
			for line in proc.stdout:
       				maps_to = str(line.split()[1])
				score = int(line.split()[2])
				rank_file.write("%s %s %s %s\n" %(unitig, rname, maps_to, score))

rank_file.close()





#Go through the top ranked candidates, and check neighbor status in graph and whether the haplotypes are consistent.

#Simpligy the graph, merge contigs

