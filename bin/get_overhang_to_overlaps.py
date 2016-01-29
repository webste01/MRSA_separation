import os, sys
import pysam
from Overlap import overlap
from collections import Counter

BLASR = "/sc/orga/projects/bashia02b/tools/bin/blasr"
NPROC = 1

def constructOverlap (bamfile, ah, edge, tolerance):
    qid = ah.query_name
    qs = ah.query_alignment_start
    qe = ah.query_alignment_end
    ql = ah.query_length
    tid = bamfile.getrname(ah.reference_id)
    ts = ah.reference_start
    te = ah.reference_end
    tl = ah.reference_length
    strand = 0
    score = ah.mapping_quality
    if ah.is_reverse:
        strand = 1
    ov = overlap(qid, tid, score, -1, strand, qs, qe, ql, 0, ts, te, tl)
    return ov

def checkReadOverhang(bamfile, ah, edge, tolerance):
    ov = constructOverlap(bamfile, ah, edge, tolerance)
    if ov.hasOverhang(edge, tolerance) and ov.hasFullOverlap(tolerance):
        return True
    return False

def findContigOverhangReads (bamfile, fastafile, contig, seqfn=False, fastafn=None, edge=2000, tolerance=200):
    
    end = fastafile.get_reference_length(contig)
    tempreadsleft = bamfile.fetch(contig, 0, edge)
    tempreadsright = bamfile.fetch(contig, end-edge, end)
    leftreadsDict = {}
    rightreadsDict = {}
    reads = []
    readSet = set()
    for ah in tempreadsleft:
        if checkReadOverhang(bamfile, ah, edge, tolerance):
            hap = 0
            try:
                hap = ah.get_tag('ZH').split(",")[1]
            except:
                pass
            leftreadsDict[ah.qname] = hap
            reads.append((ah.qname, ah.seq))
    for ah in tempreadsright:
        if checkReadOverhang(bamfile, ah, edge, tolerance):
            hap = 0
            try:
                hap = ah.get_tag('ZH').split(",")[1]
            except:
                pass
            rightreadsDict[ah.qname] = hap
            reads.append((ah.qname, ah.seq))

    if seqfn:
        # write out sequence to fasta
        seq_fa_fn = seqfn+".fa"
        rout = open(seq_fa_fn, 'w')
        print >>rout, "\n".join(map(lambda x: ">%s\n%s" %(x[0], x[1]), reads))
        rout.close()
        
        # run blasr
        blasrcmd = "%s -nproc %s -sam -out %s.sam -bestn 20 -noSplitSubreads -clipping soft %s %s" %(BLASR, NPROC, seqfn, seq_fa_fn, fastafn)
        # convert to bam
        samcmd = "samtools view -Sbh %s.sam > %s.bam" %(seqfn, seqfn)
        # sort
        sortcmd = "samtools sort %s.bam %s.sorted" %(seqfn, seqfn)
        # index
        indexcmd = "samtools index %s.sorted.bam" %(seqfn)
        for cmd in [blasrcmd, samcmd, sortcmd, indexcmd]:
            print >>sys.stderr, cmd
            os.system(cmd)

    return leftreadsDict, rightreadsDict

def AltAnchoredContigs (bamfile, bam_indexer,  fasta_file, alt_contigs, contig, read_dict, edge=2000, tolerance=200):
    '''
    Goal is to take in a list of alt contigs and potential
    reads and find out if the read is uniquely anchored to this 
    contig
    '''
    haplotypes = {}
    for read_name, hap in read_dict.items():
        read_indexer = bam_indexer.find(read_name)
        overlap_contigs = set()
        for ah in read_indexer:
            tid = bamfile.getrname(ah.reference_id)
            if tid != contig: # ignore if its the seed contig
                ov = constructOverlap(bamfile, ah, edge, tolerance)
                if ov.hasFullOverlap(): overlap_contigs.add(tid)
        if len(overlap_contigs) == 1:
            haplotypes.setdefault(hap, Counter())[tid] += 1
    return haplotypes

def printContigSummary(contig, leftjoins, rightjoins):
    print "#"*80
    print contig
    print "LEFT HAPLOTYPES:"
    for haplotype, contigCounts in leftjoins.items():
        print "-",haplotype
        for contig, contigCount in contigCounts.items():
            print "--", contig, contigCount
    print "RIGHT HAPLOTYPES:"
    for haplotype, contigCounts in rightjoins.items():
        print "-",haplotype
        for contig, contigCount in contigCounts.items():
            print "--", contig, contigCount
    print "#"*80
    
def read_branch_tigs (fn):
    """
    Ex format:
    unitig degree pred1 pred2 succ1 succ2
    unitig_54 4 unitig_55 unitig_65 unitig_13 unitig_74
    """
    bf_dict = {}
    with open(fn) as bf:
        bf.readline()
        for l in bf:
            ll = l.strip().split(" ")
            bf_dict[ll[0]]=ll[2:]
    return bf_dict
    
    
def run ():
    # index by read name
    bam_fn = sys.argv[1]
    fasta_fn = sys.argv[2]
    contig_pair_file = sys.argv[3]
    bamfile = pysam.AlignmentFile(bam_fn)
    fastafile = pysam.FastaFile(fasta_fn)
    contig_dict = read_branch_tigs (contig_pair_file)

    print >>sys.stderr, "Iterating over contigs"
    for contig, alt_contigs in contig_dict.items():
        print >>sys.stderr, contig, alt_contigs
        
        # create a temp sequence file 
        temp_seq_fn = "temp_blasr"
        temp_bam_fn = "temp_blasr.sorted.bam"
        
        # get overhang reads
        leftreads, rightreads = findContigOverhangReads(bamfile, fastafile, contig, seqfn=temp_seq_fn, fastafn=fasta_fn) 
        print >>sys.stderr, leftreads
        print >>sys.stderr, rightreads

        # index by read name
        print >>sys.stderr, "Indexing bamfile"
        smallbamfile = pysam.AlignmentFile(temp_bam_fn)
        indexer = pysam.IndexedReads(smallbamfile)
        indexer.build()

        # get alt anchors
        leftjoins = AltAnchoredContigs (smallbamfile, indexer, fastafile, alt_contigs, contig, leftreads)
        rightjoins = AltAnchoredContigs (smallbamfile, indexer, fastafile, alt_contigs, contig, rightreads)
        printContigSummary(contig, leftjoins, rightjoins)

run()        
