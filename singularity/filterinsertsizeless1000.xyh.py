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
fsingle = writeSam("{0}.singleton.bam".format(sys.argv[2]),header=fin.header)
fother = writeSam("{0}.other.bam".format(sys.argv[2]),header=fin.header)
fpet = writeSam("{0}.pet.bam".format(sys.argv[2]),header=fin.header)
for read1 in fin:
    read2 = next(fin)
    if read1.query_name != read2.query_name:
        print(read1)
        print(read2)
        sys.stderr.write("reads name are not match")
        sys.exit()
    if read1.reference_name == read2.reference_name:
        if abs(read1.isize) >= 1000:
            fpet.write(read1)
            fpet.write(read2)
        else:
            if read1.is_reverse == read2.is_reverse:
                fother.write(read1)
                fother.write(read2)
            else:
                fsingle.write(read1)
                fsingle.write(read2)
    else:
        fpet.write(read1)
        fpet.write(read2)

fsingle.close()
fpet.close()
fother.close()