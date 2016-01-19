import os, sys
import pysam

def getSubBam (bamfile, region, tag = "", edge=2000):
    chrom, start, end  = region
    tempreadsleft = bamfile.fetch(chrom, start-100, start-99)
    tempreadsright = bamfile.fetch(chrom, end+99, end+100)
    allreadsDict = {}
    reads = []
    readSet = set()
    for ar in tempreadsleft:
        allreadsDict[ar.qname] = ar.seq
        if ar.qend < len(ar.seq) - 500:
            reads.append((ar.qname, ar.seq))
            readSet.add(ar.qname)
        #if not ar.is_reverse:
        #    print region
         #   print ar.seq
         #   print ar.aligned_pairs
         #   print ar.qlen, ar.qend
        #break

    for ar in tempreadsright:
        allreadsDict[ar.qname] = ar.seq
        if ar.qstart > 500:
            if not ar.qname in readSet:
                reads.append((ar.qname, ar.seq))
    read_fa = "reads_%s_%s.fa" %("_".join(map(str, region)), tag)
    rout = open(read_fa, 'w')
    print >>rout, "\n".join(map(lambda x: ">%s\n%s" %(x[0], x[1]), reads))
    rout.close()
        
    reg = "%s:%i-%i" %(region[0], int(region[1])-edge, int(region[2])+edge) 
    ref_fa = "reg_%s.fa" %("_".join(map(str, region)))
    #faidxcmd = "samtools faidx /projects/HuPac/refs/ucsc.hg19.fasta %s > %s" %(reg, ref_fa)
    #os.system(faidxcmd)
    return allreadsDict.values(), reads

#def getReadsSpanning (bam, region):
#    outfn = 


#genomic_region = sys.argv[1]
#bamfns = sys.argv[2].split(",")
#regionfn = sys.argv[3]
#bam = pysam.Samfile(bamfn, "rb" )
regionfn = sys.argv[1]
bamfns = sys.argv[2:]
bamfiles = []
for bamfn in bamfns:
    bamfiles.append( pysam.Samfile( bamfn, "rb" ))

with open(regionfn) as f:
    for l in f:
        ll = l.strip().split()
        counter = 0
        rcount = []
        for bamfile in bamfiles:
            counter += 1
            readsfull, reads = getSubBam(bamfile, [ll[0], int(ll[1]), int(ll[2])], "bam%i" %(counter))
            rcount.append(len(readsfull))
            rcount.append(len(reads))
        print "%s\t%s" %("\t".join(map(str, ll[0:3])), "\t".join(map(str, rcount)))

        #getSpanningReads (outbam, region)
    
