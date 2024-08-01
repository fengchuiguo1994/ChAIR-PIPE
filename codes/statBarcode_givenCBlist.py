import collections,time,os,sys
import pickle


#Given cpu output sam as stdout No HEADERS, count PET or tag per CB
#pickle list of CB must exist
#sam from cpu output use numerical order of fastq in the readID

try:
    CBpkl = sys.argv[1]
except IndexError as ie:
    raise SystemError("Error: Specify pickle filename that contains CB list at 1st arg!\n")
if not os.path.exists(CBpkl):
    raise SystemError("Error: %s does not exist\n"%(CBpkl))

t0=time.time()

with open(CBpkl, 'rb') as f:
    CBlist = pickle.load(f)

#print('List of %d barcodes loaded: %.2f sec'%(len(CBlist), t1-t0), file=sys.stderr)

#
#t2=time.time()

z = {'1':'A', '2':'C', '3':'G', '4':'T', '5':'N'} #for ZB
t1=time.time()
CBdict = {} # per barcode:[list of IDs]
for line in sys.stdin: #stream the no header samtools view# if bits[0] != '@': 
    bits = line.strip()
    body = bits.split('\t')
    readID = body[0].split(':')[0] #string numerical order of fastq readID
    rNum = int(readID)-1 #for 0-based python list
    cb = CBlist[rNum]
    if (cb not in CBdict.keys()): #This barcode has not been seen previously
        CBdict[cb] = [readID]
    else:
        CBdict[cb].append(readID)

for k,v in sorted(CBdict.items()):
    print("%s\t%d"%(k, len(set(v))))

#t2=time.time()
#print('%s dict length built in  %.2f sec'%(len(CBdict), t2-t1), file=sys.stderr)

