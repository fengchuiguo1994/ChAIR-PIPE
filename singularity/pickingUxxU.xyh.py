import collections,time,os,sys
#From UxxU, pick highest AS (if exists) otherwise highest mapQ or random tied mapQ (col 5)
#Must be query sorted

t0=time.time()

#INPUT: query sorted CPU bam 
# minMAPQ = 30
minMAPQ = int(sys.argv[1])
lowcount = 0
oldNum = -1
newread = True
SAMrecs = {} # all records for one read num (readNum)
AS = {} # AS scoring per read
mapQs = {} # mapq storage if no AS per read

for line in sys.stdin: #stream with header samtools view# if bits[0] != '@':
    bits = line.strip()
    if bits[0] != '@':
        body = bits.split('\t')
        readNum = body[0].split(':')[0] #string numerical order of fastq readNum
        longName = body[0]
        str_match = list(filter(lambda x: 'AS:i:' in x, body[10:]))
        if len(SAMrecs) < 1: #no read recorded
            SAMrecs[longName] = bits
            oldNum = readNum
            if(len(str_match) > 0):
                AS[longName]  = int(str_match[0].split('i:')[1])
            else:
                mapQs[longName] = int(body[4])
        else:
            if readNum == oldNum:
                SAMrecs[longName] = bits #continue adding for the same read
                oldNum = readNum
                if(len(str_match) > 0):
                    AS[longName]  = int(str_match[0].split('i:')[1])
                else:
                    mapQs[longName] = int(body[4])
            else: 
                #print("Encountering new reads, but make decision here first: %d"%(len(SAMrecs)))
                if(len(AS) > 0):
                    x = max(AS, key=AS.get)
                    print(SAMrecs[x])
                else:
                    #print('Choose from mapQ for %s'%(longName))
                    x = max(mapQs, key=mapQs.get)
                    if mapQs[x] < minMAPQ:
                        lowcount += 1
                    else:
                        print(SAMrecs[x])
                SAMrecs = {} #reset
                AS = {} #reset
                mapQs = {} #reset
                oldNum = readNum
                SAMrecs[longName] = bits
                if(len(str_match) > 0):
                    AS[longName] = int(str_match[0].split('i:')[1])
                else:
                    mapQs[longName] = int(body[4])
    else:
        print(bits)

#The last bit is out of loop
if(len(AS) > 0):
    x = max(AS, key=AS.get)
    print(SAMrecs[x])
else:
    x = max(mapQs, key=mapQs.get)
    if mapQs[x] < minMAPQ:
        lowcount += 1
    else:
        print(SAMrecs[x])


t2=time.time()
print('#of reads discarded MAPQ<%d: %d'%(minMAPQ, lowcount), file=sys.stderr)
print('Total time spent: %.2f sec'%(t2-t0), file=sys.stderr)

