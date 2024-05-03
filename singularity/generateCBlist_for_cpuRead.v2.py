import collections,time,os,sys,pysam
import re,gzip,pickle

#The readID in CPU output is sequential order of fastq entries 
#This code is to get CB based on this order
#CB list generated include all fastq reads that have noCB

try:
    file10x = sys.argv[1]
except IndexError as ie:
    raise SystemError("Error: Specify 10x bam file 1st arg\n")

try:
    fastq   = sys.argv[2]
except IndexError as ie:
    raise SystemError("Error: Specify fasq.gz file 2nd arg\n")

try:
    OUTNAME = sys.argv[3]
except IndexError as ie:
    raise SystemError("Error: Specify prefix output file at 3rd arg\n")

if not os.path.exists(file10x):
    raise SystemError("Error: %s does not exist\n"%(file10x))
if not os.path.exists(fastq):
    raise SystemError("Error: %s does not exist\n"%(fastq))



savedict = OUTNAME + '.CB_forCPU.pkl'

t0=time.time()

#Read 10x bam file
bam10x = pysam.Samfile(file10x, "rb")
iter10x = bam10x.fetch()

s = pysam.flagstat(file10x)
s = s.split('\n') #flagstat lines
nread = int(s[6].split()[0])

t1=time.time()
print('%s has %s reads; %.2f sec'%(file10x, nread, t1-t0))

CBdict = {}
for x in iter10x:
    TAG=dict(x.tags)
    r=x.query_name
    if ('CB' in TAG.keys()):
        CBdict[r]=TAG["CB"]

t2=time.time()
print('Dict of %d barcodes created: %.2f sec'%(len(CBdict), t2-t1))

#Read fastq
lcount = 0
lnocb  = 0
orderedCB = ['noCB'] * nread #the BC ordered by fastq record
k = -1

with gzip.open(fastq, 'rb') as fh:
    lines = []
    line = fh.readline()
    while line:
        lcount += 1 
        lines.append(line.decode('utf8').rstrip())
        if len(lines) == 4:
            k += 1
            r = lines[0].split(' ')[0][1:]
            if (r in CBdict.keys()):
                orderedCB[k] = CBdict[r]
            else:
                orderedCB[k] = 'noCB'
                lnocb += 1
            lines = []
        line = fh.readline()

t3=time.time()
print('%d out of %d FASTQ reads were noCB; time elapsed  %.2f secs'%(lnocb, lcount/4, t3-t2))


with open(savedict, 'wb') as f:
    pickle.dump(orderedCB, f)

t4=time.time()

print('%d CB fastq saved! Total elapse time : %.2f sec'%(len(orderedCB), t4-t0) )
