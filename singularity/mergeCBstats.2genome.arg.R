#Create a final table from 5 CB count tsv files
#Runs under bulkChiatac/  
args <- commandArgs(trailingOnly = TRUE);
RUN <- args[1]; 
genome <- args[2]; 

message('Combining stats for: ', RUN, ' ', genome)
#Collect all tsv files: UU, UxxU, UxxU.F2048.CB.nr
singlabel='.singlelinker.single'
pairlabel='.singlelinker.paired'
nonelabel='.none'


#CPU aligned sam.gz counts:
none = read.delim(paste0(RUN,nonelabel,'.statCB.totalPET.tsv'), check.names=FALSE, header=FALSE) #direct from fastq file
stag = read.delim(paste0(RUN,singlabel,'.statCB.totalPET.tsv'), check.names=FALSE, header=FALSE) #direct from fastq file
tpet = read.delim(paste0(RUN,pairlabel,'.statCB.totalPET.tsv'), check.names=FALSE, header=FALSE) #direct from fastq file

#message(paste0(RUN,nonelabel,'.statCB.totalPET.tsv'), ' rows: ', nrow(none))
#message(paste0(RUN,singlabel,'.statCB.totalPET.tsv'), ' rows: ', nrow(stag))
#message(paste0(RUN,pairlabel,'.statCB.totalPET.tsv'), ' rows: ', nrow(tpet))
colnames(none) = c('CorrectBarcode', 'NoLinker_tot')
colnames(stag) = c('CorrectBarcode', '1tag_tot')
colnames(tpet) = c('CorrectBarcode', '2tag_tot')

a=merge(none, stag, all=TRUE)
a=merge(a, tpet, all=TRUE)

#paired.UU.bam 

F2048pet=read.delim(paste0(genome, '/', RUN,'.',genome, pairlabel,'.UU.F2048.CB.PET.tsv'), check.names=FALSE, header=TRUE) #from -F 2048 paired.UU.bam
colnames(F2048pet)[2:6] = paste0(colnames(F2048pet)[2:6],'_',genome)
nrpet=read.delim(paste0(genome, '/', RUN,'.',genome, pairlabel,'.UU.F2048.CB.nrPET.tsv'), check.names=FALSE, header=TRUE) #non-redundant F2048 
colnames(nrpet)[2:6] = paste0(colnames(nrpet)[2:6],'_nr_',genome)
petbc = merge(F2048pet, nrpet, all=TRUE)

pairUxxU=read.delim(paste0(genome, '/', RUN,'.',genome, pairlabel,'.UxxU.F2048.CB.tsv'), check.names=FALSE, header=FALSE) #2-column BC stat -F 2048
pairUxxUnr=read.delim(paste0(genome, '/', RUN,'.',genome, pairlabel,'.UxxU.F2048.CB.nr.tsv'), check.names=FALSE, header=FALSE) #2-column BC stat -F 2048
singUxxU=read.delim(paste0(genome, '/', RUN,'.',genome, singlabel,'.UxxU.F2048.CB.tsv'), check.names=FALSE, header=FALSE) #2-column BC stat -F 2048
singUxxUnr=read.delim(paste0(genome, '/', RUN,'.',genome, singlabel,'.UxxU.F2048.CB.nr.tsv'), check.names=FALSE, header=FALSE) #2-column BC stat -F 2048
noneUxxU=read.delim(paste0(genome, '/', RUN,'.',genome, nonelabel,'.UxxU.F2048.CB.tsv'), check.names=FALSE, header=FALSE) #2-column BC stat -F 2048
noneUxxUnr=read.delim(paste0(genome, '/', RUN,'.',genome, nonelabel,'.UxxU.F2048.CB.nr.tsv'), check.names=FALSE, header=FALSE) #2-column BC stat -F 2048

F2048none=read.delim(paste0(genome, '/', RUN,'.',genome, nonelabel,'.UU.F2048.CB.tsv'), check.names=FALSE, header=FALSE) #from -F 2048 2-column BC stat
nrnone=read.delim(paste0(genome, '/', RUN,'.',genome, nonelabel,'.UU.F2048.CB.nr.tsv'), check.names=FALSE, header=FALSE) #from -F 2048 2-column BC stat

colnames(pairUxxU) = c('CorrectBarcode', paste0('2tag_Ux_',genome))
colnames(pairUxxUnr) = c('CorrectBarcode', paste0('2tag_Ux_nr_',genome))
colnames(singUxxU) = c('CorrectBarcode', paste0('1tag_Ux_',genome))
colnames(singUxxUnr) = c('CorrectBarcode', paste0('1tag_Ux_nr_',genome))
colnames(noneUxxU) = c('CorrectBarcode', paste0('NoLinker_Ux_',genome))
colnames(noneUxxUnr) = c('CorrectBarcode', paste0('NoLinker_Ux_nr_',genome))
colnames(F2048none) = c('CorrectBarcode', paste0('NoLinker_UU_',genome))
colnames(nrnone) = c('CorrectBarcode', paste0('NoLinker_UU_nr_',genome))


ux1  = merge(noneUxxU, singUxxU,  all=TRUE)
ux2  = merge(ux1, pairUxxU,  all=TRUE)
f2048.df = merge(F2048none, ux2,  all=TRUE)
ux1nr  = merge(noneUxxUnr, singUxxUnr,  all=TRUE)
ux2nr  = merge(ux1nr, pairUxxUnr,  all=TRUE)
f2048nr.df = merge(nrnone, ux2nr,  all=TRUE)
f2048all = merge(f2048.df, f2048nr.df, all=TRUE)

atab = merge(a, f2048all, all=TRUE)
atab = merge(atab, petbc, all=TRUE)
atab[is.na(atab)] = 0
write.table(atab, file=paste0(RUN, '.', genome, '.CBcounts.tsv'), quote=F, row.names=F, col.names=T, sep="\t")

#sums
tot = colSums(atab[,2:22])
totnocb = colSums(subset(atab, CorrectBarcode!='noCB')[,2:22])
write.table(rbind(tot, totnocb), file=paste0(RUN, '.', genome, '.CBsums.tsv'), quote=F, row.names=F, col.names=T, sep="\t")
