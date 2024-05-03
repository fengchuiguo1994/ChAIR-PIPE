import collections,time,os,sys
import pickle


#Given cpu output sam as stdout, add CB & ZB tags
#pickle list of CB must exist
#sam from cpu output use numerical order of fastq in the readID

try:
    CBpkl = sys.argv[1]
except IndexError as ie:
    raise SystemError("Error: Specify pickle filename that contains CB list at 1st arg!\n")
if not os.path.exists(CBpkl):
    raise SystemError("Error: %s does not exist\n"%(CBpkl))

#try:
#    cpubam = sys.argv[2]
#except IndexError as ie:
#    raise SystemError("Error: Specify CPU bam filename at arg2!\n")
#if not os.path.exists(cpubam):
#    raise SystemError("Error: %s does not exist\n"%(cpubam))
#
t0=time.time()

with open(CBpkl, 'rb') as f:
    CBlist = pickle.load(f)

t1=time.time()
nCB = len(CBlist)
print('List of %d barcodes loaded: %.2f sec'%(nCB, t1-t0), file=sys.stderr)

#fnocb = open('noCB.tmp.txt', 'w') 
#bam = pysam.Samfile(cpubam, "rb")
#iterbam = bam.fetch()
#
#t2=time.time()
#print('%s  loaded: %.2f sec'%(cpubam, t2-t1), file=sys.stderr)
z = {'1':'A', '2':'C', '3':'G', '4':'T', '5':'N'} #for ZB
oldNum = 0
for line in sys.stdin:
    bits = line.strip()
    if bits[0] != '@':
        body = bits.split('\t')
        n = len(body)-1
        rNum = int(body[0].split(':')[0])-1
        cb = CBlist[rNum]
        if (cb != 'noCB'):
            cbs = cb.split('-')
            zb = 'ZB:Z:' + cbs[0] + '-' + z[cbs[1]]
            body.insert(n,zb)
        cb = 'CB:Z:' + cb
        body.insert(n,cb) #ZB is for MarkDuplicates, if no CB then no need ZB
        samrec = '\t'.join(body)
        print(samrec)
    else:
        print(bits)


t2=time.time()

#fnocb.close()
