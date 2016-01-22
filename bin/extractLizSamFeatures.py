import os, sys
import pysam

samfn = sys.argv[1]

samfile = pysam.AlignmentFile(samfn, "r" )

for ah in samfile.fetch():
    qid = ah.query_name
    qs = ah.query_alignment_start
    qe = ah.query_alignment_end
    ql = ah.query_length
    tid = samfile.getrname(ah.reference_id)
    ts = ah.reference_start
    te = ah.reference_end
    tl = ah.reference_length
    strand = "+"
    if ah.is_reverse:
        strand = "-"
    print "%s %i %i %i %s %i %i %i %s" %(qid, qs, qe, ql, tid, ts, te, tl, strand)
