import sys
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
fout = pysam.AlignmentFile("-", "w", template=fin)
for read1 in fin:
    read2 = next(fin)
    if read1.query_name != read2.query_name:
        print(read1)
        print(read2)
        sys.stderr.write("reads name are not match")
        sys.exit()
    if read1.reference_name != read2.reference_name:
        print(read1)
        print(read2)
        sys.stderr.write("reads chromosome are not match")
        sys.exit()
    start = read1.reference_start
    if start > read2.reference_start:
        start = read2.reference_start
    if start < 0:
        start = 0
    end = read1.reference_end
    if end < read2.reference_end:
        end = read2.reference_end

    a = pysam.AlignedSegment()
    a.query_name = read1.query_name
    a.query_sequence=""
    a.query_qualities = pysam.qualitystring_to_array("")
    a.flag = 0
    if read1.is_reverse:
        a.flag = 16
    a.reference_id = read1.reference_id
    a.reference_start = start
    a.mapping_quality = 60
    a.cigar = ((0,end-start),)
    a.tags = (("CB", read1.get_tag("CB")), ("ZB", read1.get_tag("ZB")))
    fout.write(a)
fout.close()