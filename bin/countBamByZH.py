import os 
import sys
import pysam
from collections import Counter
bam_fn = sys.argv[1]
contig_id = sys.argv[2]
bamfile = pysam.AlignmentFile(bam_fn)
hapCounter = Counter()
totalcount = 0
secondarycount = 0
for ah in bamfile.fetch(contig_id):
    totalcount += 1
    if ah.is_secondary:
        secondarycount += 1
    try:
        ah_hap = ah.get_tag('ZH').split(",")[1]
        hapCounter[ah_hap] += 1
    except:
        pass

print secondarycount, totalcount
for key, value in hapCounter.items(): print key, value


