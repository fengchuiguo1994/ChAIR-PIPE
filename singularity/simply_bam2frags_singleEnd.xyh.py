import sys
import pysam
from collections import OrderedDict
 
import gzip
def readFile(infile):
    if infile.endswith((".gz","gzip")):
        fin = gzip.open(infile,'rt')
    else:
        fin = open(infile,'r')
    return fin
        
def writeFile(outfile):
    if outfile.endswith((".gz","gzip")):
        fout = gzip.open(outfile,'wt')
    else:
        fout = open(outfile,'w')
    return fout

import pysam
def readSam(insamfile):
    if insamfile.endswith(".bam"):
        insam = pysam.AlignmentFile(insamfile,'rb')
    elif insamfile.endswith(".sam.gz"):
        insam = pysam.AlignmentFile(insamfile,'rb')
    elif insamfile.endswith(".sam"):
        insam = pysam.AlignmentFile(insamfile,'r')
    else:
        raise ValueError("the input sam/bam file is not end with sam or bam!")
    return insam
        
def writeSam(outsamfile,header):
    if outsamfile.endswith(".bam"):
        outsam = pysam.AlignmentFile(outsamfile,'wb',header=header)
    elif outsamfile.endswith(".sam"):
        outsam = pysam.AlignmentFile(outsamfile,'w',header=header)
    else:
        raise ValueError("the output sam/bam file is not end with sam or bam!")
    return outsam

fin = readSam(sys.argv[1])
fout = writeFile(sys.argv[2])
fdis = writeFile("{0}.len".format(sys.argv[2]))

read_file_counter = 0
frags = OrderedDict()

for read in fin:
    read_file_counter += 1
    chrom = read.reference_name
    if read.is_reverse:
        start = read.reference_start
        end = read.reference_end - 5
    else:
        start = read.reference_start + 4
        end = read.reference_end
    rlen = end - start
    fdis.write("{0}\n".format(rlen))
    barcode = read.get_tag("CB")
    frag_id = "{0}\t{1}\t{2}\t{3}".format(chrom,start,end,barcode)
    if frag_id not in frags:
        frags[frag_id] = 0
    frags[frag_id] += 1

fdis.close()
fin.close()

for k in frags.keys():
    fout.write("{0}\t{1}\n".format(k,frags[k]))
fout.close()
