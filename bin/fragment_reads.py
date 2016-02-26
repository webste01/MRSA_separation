import os
import sys
from pbcore.io.FastaIO import FastaReader

inFN = sys.argv[1] # inputfofn
chunklen = int(sys.argv[2])

with open(inFN) as f:
    for l in f:
        fn = l.strip()
        print >>sys.stderr, fn
        for entry in FastaReader(fn):
            # fragment reads
            seq = entry.sequence
            sid = entry.name.split()[0]
            for i in range(0, len(seq), chunklen):
                newid = "%s/%i" %(sid, i)
                subseq = seq[i:i+chunklen]
                print ">%s\n%s" %(newid, subseq)


