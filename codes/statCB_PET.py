import collections,time,os,sys
#import pickle
#Per barcode: total PETs, cisPETs, transPETs, cisPETs > 8kb, ratio (cisPETs/transPETs)

t0=time.time()

#INPUT: query sorted CPU bam with CB tag and deduped therefore clean from noCB
#verify the bam input: samtools view cb.nr.nsorted.bam |cut -f1 |awk -F":" '{print $1}' |uniq -c | awk '{print $1}' |sort |uniq -c

#list BC: [nPET, nTrans, nCis, >1kb, >8kb ] #use this order!
nstat = 5
d1k = 1000
d8k = 8000
CBdict = {}
oldID = 'new'
for line in sys.stdin: #stream the no header samtools view# if bits[0] != '@':
    bits = line.strip()
    body = bits.split('\t')
    readID = body[0].split(':')[0] #string numerical order of fastq readID
    str_match = list(filter(lambda x: 'CB:Z:' in x, body[10:]))
    cb = str_match[0].split(':Z:')[1]
    if (cb not in CBdict.keys()): #This barcode has not been seen previously
        CBdict[cb] = [0] * nstat
    if oldID != readID: #this is a PE mate
        chrom1 = body[2]
        chrom2 = body[6]
        cisdis = abs(int(body[8]))
        CBdict[cb][0] += 1 #PET count
        if chrom2 != '=':
            CBdict[cb][1] += 1 #Trans count
        else:
            CBdict[cb][2] += 1 #Cis count
            if (cisdis >= d1k):
                CBdict[cb][3] += 1 #Cis >1kb count
            if (cisdis >= d8k):
                CBdict[cb][4] += 1 #Far cis count

        oldID = readID  #set for next record

t1=time.time()

print("CorrectBarcode\tnPET\tnTrans\tnCis\t>1kb\t>8kb")
for k,v in sorted(CBdict.items()):
    print("%s\t%d\t%d\t%d\t%d\t%d"%(k,v[0], v[1], v[2], v[3], v[4] ))

t2=time.time()
print('Total time spent: %.2f sec'%(t2-t0), file=sys.stderr)

