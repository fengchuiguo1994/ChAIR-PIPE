#!/bin/bash
# genome a string e.g. mm10 or GRCh38
# Execute this script under bulkChiatac folder
# Extract summary statistics from scChiatac output files
# 2022 tjongh
# The Jackson Laboratory for Genomic Medicine

die () {
 echo >&2 "$@"
 exit 1
}
[ $# -gt 1 ] || die "RUN & GENOME1 variables must be given!"

function waitFile {
    while [ ! -f $1 ]
    do
      echo "$1 is missing, wait for 5 minutes..."
      date
      sleep 300
    done
}

RUN=$1
GENOME1=$2
hybridName=$3
samdir='dump_intermediate'
samdir='.'

grun=${RUN}.${GENOME1}
pairlabel='singlelinker.paired'
singlabel='singlelinker.single'
nonelabel='none'
extbp=500

## Set the output file
out_file=${grun}.bulk_stats.tsv

samcmd='samtools flagstat'
#------------------------------------


echo Starting to sumarize ${grun}
echo "Output file (overwrite mode): $out_file" 

#`grep export ../slurm/${RUN}.map.${pairlabel}.slurm`
## Get library ID

cat ${RUN}.stat | grep "Filename" | \
    awk -F'[ \t.]' '{print "Library_ID\t"$2}' > ${out_file}

# Get reference genome
echo -e "Genome_origin\t"${GENOME1} >> ${out_file}


## PET count
# Get PET count
n_pet=$( cat ${RUN}.stat | grep "Total pairs" | awk -F'[ \t]' '{print $3}' )

# Write PET count
n_pet=$( printf "%'.f" ${n_pet} )
echo -e "Raw_reads\t"${n_pet} >> ${out_file}

## Get PET count with linker
read_link=$( cat ${RUN}.stat | grep "Linker detected" | \
    awk -F '[ \t]' '{print $3}' | xargs printf "%'.f")
echo -e "Linker_detected_PEreads\t"${read_link} >> ${out_file}

## Get PET count with linker
pet_link=$( cat ${RUN}.stat | grep "Single Linker 2 tags" | \
    awk -F '[ \t]' '{print $6}' | xargs printf "%'.f")
echo -e "Raw_PET_with_linker\t"${pet_link} >> ${out_file}

cbtsv="${grun}.CBcounts.tsv"
cbsumtsv="${grun}.CBsums.tsv"
waitFile $cbtsv #check if it exists, if not just wait
waitFile $cbsumtsv #check if it exists, if not just wait
noneRAW=$(awk 'NR>1 {sum +=$2} END {print sum}' $cbtsv | xargs printf "%'.f") #consistent with none.stat.xls
singRAW=$(awk 'NR>1 {sum +=$3} END {print sum}' $cbtsv | xargs printf "%'.f") #consistent with singlelinker.single.stat.xls
pairRAW=$(awk 'NR>1 {sum +=$4} END {print sum}' $cbtsv | xargs printf "%'.f") #consistent with singlelinker.paired.stat.xls

#echo -e "paired_category\t"${pairRAW} >> ${out_file} #the same as $pet_link
echo -e "Raw_1tag_with_linker\t"${singRAW} >> ${out_file}
echo -e "Raw_NoLinker\t"${noneRAW} >> ${out_file}

## Mapping
# Get mapped PET count 

# Get uniquely mapped PET count 
waitFile ${RUN}.$pairlabel.stat.xls #check if it exists, if not just wait
unique=$( cat ${RUN}.$pairlabel.stat.xls | grep "Uniquely Mapped" | \
    awk -F '[\t]' '{print $3}' | xargs printf "%'.f" )
punique=$( cat ${RUN}.$pairlabel.stat.xls | grep "Uniquely Mapped" | \
    awk -F '[\t]' '{print $4}' )
echo -e "PET_UU_all_mapped\t"${unique} >> ${out_file}
#echo -e "%_PET_UU_all_mapped\t"${punique} >> ${out_file}

uxxu=$( cat ${RUN}.$pairlabel.stat.xls | grep "Ux/xU" | \
    awk -F '[\t]' '{print $3}' | xargs printf "%'.f" )
puxxu=$( cat ${RUN}.$pairlabel.stat.xls | grep "Ux/xU" | \
    awk -F '[\t]' '{print $4}' )
echo -e "PET_UxxU_all_mapped\t"${uxxu} >> ${out_file}
#echo -e "%_PET_UxxU_mapped\t"${puxxu} >> ${out_file}

waitFile ${RUN}.$singlabel.stat.xls #check if it exists, if not just wait
u1tag=$( cat ${RUN}.$singlabel.stat.xls | grep "Ux/xU" | \
    awk -F '[\t]' '{print $3}' | xargs printf "%'.f" )
pu1tag=$( cat ${RUN}.$singlabel.stat.xls | grep "Ux/xU" | \
    awk -F '[\t]' '{print $4}' )
echo -e "1tag_UxxU_all_mapped\t"${u1tag} >> ${out_file}
#echo -e "%_1tag_UxxU_mapped\t"${pu1tag} >> ${out_file}


waitFile ${RUN}.$nonelabel.stat.xls #check if it exists, if not just wait
unone=$( cat ${RUN}.$nonelabel.stat.xls | grep "Uniquely Mapped" | \
    awk -F '[\t]' '{print $3}' | xargs printf "%'.f" )
punone=$( cat ${RUN}.$nonelabel.stat.xls | grep "Uniquely Mapped" | \
    awk -F '[\t]' '{print $4}' )
echo -e "NoLinker_UU_all_mapped\t"${unone} >> ${out_file}
#echo -e "%_NoLinker_UU_mapped\t"${punone} >> ${out_file}

uxxunone=$( cat ${RUN}.$nonelabel.stat.xls | grep "Ux/xU" | \
    awk -F '[\t]' '{print $3}' | xargs printf "%'.f" )
echo -e "NoLinker_UxxU_all_mapped\t"${uxxunone} >> ${out_file}

#Count hybrid
pairUUsam=$samdir/${RUN}.${pairlabel}.UU.${GENOME1}.sam.gz
pairUxxUsam=$samdir/${RUN}.${pairlabel}.UxxU.${GENOME1}.sam.gz
singUxxUsam=$samdir/${RUN}.${singlabel}.UxxU.${GENOME1}.sam.gz
noneUxxUsam=$samdir/${RUN}.${nonelabel}.UxxU.${GENOME1}.sam.gz
noneUUsam=$samdir/${RUN}.${nonelabel}.UU.${GENOME1}.sam.gz
pairUUmixsam=$samdir/${RUN}.${pairlabel}.UU.mix_${hybridName}.sam.gz
noneUUmixsam=$samdir/${RUN}.${nonelabel}.UU.mix_${hybridName}.sam.gz

waitFile $pairUUsam #check if it exists, if not just wait
waitFile $pairUxxUsam #check if it exists, if not just wait
waitFile $singUxxUsam #check if it exists, if not just wait
waitFile $noneUUsam  #check if it exists, if not just wait
waitFile $pairUUmixsam #
waitFile $noneUUmixsam #
waitFile $noneUxxUsam  #check if it exists, if not just wait

pairUUcount=$( $samcmd $pairUUsam | grep "read1" | \
    awk '{print $1}' | xargs printf "%'.f" ) 
pairUUmixcount=$( $samcmd $pairUUmixsam | grep "read1" | \
    awk '{print $1}' | xargs printf "%'.f" ) 
noneUUcount=$( $samcmd $noneUUsam | grep "read1" | \
    awk '{print $1}' | xargs printf "%'.f" ) 
noneUUmixcount=$( $samcmd $noneUUmixsam | grep "read1" | \
    awk '{print $1}' | xargs printf "%'.f" ) 
pairUxxUcount=$( $samcmd $pairUxxUsam | grep "primary mapped" | \
    awk '{print $1}' | xargs printf "%'.f" ) 
singUxxUcount=$( $samcmd $singUxxUsam | grep "primary mapped" | \
    awk '{print $1}' | xargs printf "%'.f" ) 
noneUxxUcount=$( $samcmd $noneUxxUsam | grep "primary mapped" | \
    awk '{print $1}' | xargs printf "%'.f" ) 


noneUxCB=$(awk 'NR==1 {for (i=1; i<=NF; ++i) if ($i ~ "NoLinker_Ux")  {n=i; break}} {print $n}'  $cbsumtsv |head -2 |tail -1 | xargs printf "%'.f") 
singUxCB=$(awk 'NR==1 {for (i=1; i<=NF; ++i) if ($i ~ "1tag_Ux")  {n=i; break}} {print $n}'  $cbsumtsv |head -2 |tail -1 | xargs printf "%'.f") 
pairUxCB=$(awk 'NR==1 {for (i=1; i<=NF; ++i) if ($i ~ "2tag_Ux")  {n=i; break}} {print $n}'  $cbsumtsv |head -2 |tail -1 | xargs printf "%'.f") 

echo -e "PET_UxxU_mapped_${GENOME1}\t"${pairUxxUcount} >> ${out_file}
echo -e "1tag_UxxU_mapped_${GENOME1}\t"${singUxxUcount} >> ${out_file}
echo -e "NoLinker_UxxU_mapped_${GENOME1}\t"${noneUxxUcount} >> ${out_file}
echo -e "PET_UU_q30_mix_${hybridName}\t"${pairUUmixcount} >> ${out_file}
echo -e "NoLinkerUU_q30_mix_${hybridName}\t"${noneUUmixcount} >> ${out_file}
echo -e "PET_UU_q30_${GENOME1}\t"${pairUUcount} >> ${out_file}
echo -e "NoLinker_UU_q30_${GENOME1}\t"${noneUUcount} >> ${out_file}
echo -e "PET_Ux_q30_${GENOME1}\t"${pairUxCB} >> ${out_file}
echo -e "NoLinker_Ux_q30_${GENOME1}\t"${noneUxCB} >> ${out_file}
echo -e "1tag_Ux_q30_${GENOME1}\t"${singUxCB} >> ${out_file}

pairMarkdupFile=$GENOME1/${grun}.${pairlabel}.UU.F2048.CB.nr.MarkDuplicates.txt #consistent with singlelinker.paired.UU.F2048.CB.PET.tsv
uxxuMarkdupFile=$GENOME1/${grun}.${pairlabel}.UxxU.F2048.CB.nr.MarkDuplicates.txt #
singMarkdupFile=$GENOME1/${grun}.${singlabel}.UxxU.F2048.CB.nr.MarkDuplicates.txt
noneMarkdupFile=$GENOME1/${grun}.${nonelabel}.UU.F2048.CB.nr.MarkDuplicates.txt
noneUxxUMarkdupFile=$GENOME1/${grun}.${nonelabel}.UxxU.F2048.CB.nr.MarkDuplicates.txt
#LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	SECONDARY_OR_SUPPLEMENTARY_RDS	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE
waitFile $pairMarkdupFile
waitFile $uxxuMarkdupFile
waitFile $singMarkdupFile
waitFile $noneMarkdupFile
pairExam=$(grep "Unknown Library"  $pairMarkdupFile | awk -F '[\t]' '{print $3}' |  xargs printf "%'.f" )
uxxuExam=$(grep "Unknown Library"  $uxxuMarkdupFile | awk -F '[\t]' '{print $2}' |  xargs printf "%'.f" )
singExam=$(grep "Unknown Library"  $singMarkdupFile | awk -F '[\t]' '{print $2}' |  xargs printf "%'.f" )
noneExam=$(grep "Unknown Library"  $noneMarkdupFile | awk -F '[\t]' '{print $3}' |  xargs printf "%'.f" )
noneUxxUExam=$(grep "Unknown Library"  $noneUxxUMarkdupFile | awk -F '[\t]' '{print $2}' |  xargs printf "%'.f" )

echo -e "PET_UU_q30_validCB_${GENOME1}\t"${pairExam} >> ${out_file} #consistent with CBcounts.tsv without noCB e.g.: $(awk 'NR>1 && $1 != "noCB" {sum +=$5} END {print sum}' $cbtsv | xargs printf "%'.f")
echo -e "PET_UxxU_q30_validCB_${GENOME1}\t"${uxxuExam} >> ${out_file}
echo -e "1tag_UxxU_q30_validCB_${GENOME1}\t"${singExam} >> ${out_file}
echo -e "NoLinker_UU_q30_validCB_${GENOME1}\t"${noneExam} >> ${out_file}
echo -e "NoLinker_UxxU_q30_validCB_${GENOME1}\t"${noneUxxUExam} >> ${out_file}

uuPairDupRate=$(grep "Unknown Library"  $pairMarkdupFile | awk -F '[\t]' '{print $9*100}' |  xargs printf "%.1f%%" )
uxxuPairDupRate=$(grep "Unknown Library"  $uxxuMarkdupFile | awk -F '[\t]' '{print $9*100}' |  xargs printf "%.1f%%" )
uxxuSingDupRate=$(grep "Unknown Library"  $singMarkdupFile | awk -F '[\t]' '{print $9*100}' |  xargs printf "%.1f%%" )
uuNoneDupRate=$(grep "Unknown Library"  $noneMarkdupFile | awk -F '[\t]' '{print $9*100}' |  xargs printf "%.1f%%" )
uxxuNoneDupRate=$(grep "Unknown Library"  $noneUxxUMarkdupFile | awk -F '[\t]' '{print $9*100}' |  xargs printf "%.1f%%" )
echo -e "Redundancy_PET_UU_${GENOME1}\t"${uuPairDupRate} >> ${out_file} 
echo -e "Redundancy_PET_UxxU_${GENOME1}\t"${uxxuPairDupRate} >> ${out_file} 
echo -e "Redundancy_1tag_UxxU_${GENOME1}\t"${uxxuSingDupRate} >> ${out_file} 
echo -e "Redundancy_NoLinker_UU_${GENOME1}\t"${uuNoneDupRate} >> ${out_file} 
echo -e "Redundancy_NoLinker_UxxU_${GENOME1}\t"${uxxuNoneDupRate} >> ${out_file} 

nrPETuu=$(awk 'NR>1 && $1 != "noCB" {sum +=$18} END {print sum}' $cbtsv | xargs printf "%'.f")
echo -e "nrPET_UU_${GENOME1}\t"${nrPETuu} >> ${out_file} 

nrPETtra=$(awk 'NR>1 && $1 != "noCB" {sum +=$19} END {print sum}' $cbtsv)
printf  "nrPET_trans_${GENOME1}\t%'.f\n" ${nrPETtra} >> ${out_file} 

nrPETcis=$(awk 'NR>1 && $1 != "noCB" {sum +=$20} END {print sum}' $cbtsv)
nrPET1kb=$(awk 'NR>1 && $1 != "noCB" {sum +=$21} END {print sum}' $cbtsv )
nrPET8kb=$(awk 'NR>1 && $1 != "noCB" {sum +=$22} END {print sum}' $cbtsv )

pet_ratio=$( echo "${nrPETcis} / ${nrPETtra}" | bc -l )
pet_ratio1=$( echo "${nrPET1kb} / ${nrPETtra}" | bc -l )
pet_ratio8=$( echo "${nrPET8kb} / ${nrPETtra}" | bc -l )

printf  "nrPET_cis\t%'.f\n" ${nrPETcis} >> ${out_file} 
printf  "nrPET_cis>1kb\t%'.f\n" ${nrPET1kb} >> ${out_file} 
printf  "nrPET_cis>8kb\t%'.f\n" ${nrPET8kb} >> ${out_file} 

printf  "Ratio_cis/trans_PET\t%.2f\n" ${pet_ratio} >> ${out_file} 
printf  "Ratio_cis>1kb/trans_PET\t%.2f\n" ${pet_ratio1} >> ${out_file} 
printf  "Ratio_cis>8kb/trans_PET\t%.2f\n" ${pet_ratio8} >> ${out_file} 



## Singleton
cis_file=$GENOME1/${grun}.*cis.chiasig.gz
tra_file=$GENOME1/${grun}.*trans.chiasig.gz
waitFile $cis_file
waitFile $tra_file
# Get singleton PET count (all)
singleton=$(zcat $GENOME1/${grun}.*chiasig.gz | awk '$7==1{print}' | wc -l)
singleton=$( printf "%'.f" ${singleton} )
echo -e "Singleton\t"$singleton >> ${out_file}

# Get singleton PET count (intra-chr)
intra_singleton=$(zcat ${cis_file} | awk '$7==1{print}' | wc -l)
intra_singleton=$( printf "%'.f" ${intra_singleton} )
echo -e "Intra-chr_singleton\t"$intra_singleton >> ${out_file}

# Get singleton PET count (inter-chr)
inter_singleton=$(zcat ${tra_file} | awk '$7==1{print}' | wc -l)
inter_singleton=$( printf "%'.f" ${inter_singleton} )
echo -e "Inter-chr_singleton\t"$inter_singleton >> ${out_file}

## Clusters (overall)
# Get cluster count
total_cluster_number=$(zcat $GENOME1/${grun}.*chiasig.gz | awk '$7 > 1 {print}' | wc -l)
total_cluster_number=$( printf "%'.f" ${total_cluster_number} )
echo -e "PET_cluster\t"${total_cluster_number} >> ${out_file}

intra_cluster=$(zcat ${cis_file} | awk '$7 > 1 {print}' | wc -l)
inter_cluster=$(zcat ${tra_file} | awk '$7 > 1 {print}' | wc -l)
cluster_ratio=$( echo "${intra_cluster} / ${inter_cluster}" | bc -l )

printf  "Cis_clusters\t%'.f\n" ${intra_cluster} >> ${out_file}
printf  "Trans_clusters\t%'.f\n" ${inter_cluster} >> ${out_file}
printf  "Ratio_cis/trans_clusters\t%.2f\n" ${cluster_ratio} >> ${out_file}

# Get intra-chr cluster count by number of PETs (1 - 10)
for i in $(seq 2 10)
do
        intra_pets_number=$(zcat ${cis_file} | \
            awk -v cutoff=${i} '$7 == cutoff {print}' | wc -l | \
            xargs printf "%'.f")

        echo -e "pets_number_"${i}"\t"${intra_pets_number} >> ${out_file}
done

# Get intra-chr cluster count with > 10 PETs
echo -e "pets_number>10\t"$(zcat ${cis_file} | awk '$7 >10 {print}' | \
    wc -l | xargs printf "%'.f") >> ${out_file}


## Clusters (inter-chr)
# Write inter-chr cluster count
inter_cluster=$( printf "%'.f" ${inter_cluster} )
echo -e "Inter-chr_PET_cluster\t"${inter_cluster} >> ${out_file}


# Get inter-chr cluster count by number of PETs (1 - 10)
for i in $(seq 2 10)
do
        inter_pets_number=$(zcat ${tra_file}  | \
            awk -v cutoff=${i} '$7 == cutoff {print}' | wc -l | \
            xargs printf "%'.f")

        echo -e "pets_number_"${i}"\t"${inter_pets_number} >> ${out_file}
done

# Get inter-chr cluster count with > 10 PETs
echo -e "pets_number>10\t"$(zcat ${tra_file}  | \
    awk '$7 >10 {print}' | wc -l | xargs printf "%'.f") >> ${out_file}

#significant interaction
cis_prefix="${grun}.e${extbp}.clusters.cis.chiasig"
chiasigoutput="$GENOME1/${cis_prefix}.BE2.sigf.interactions"
ns=$( grep -v : $chiasigoutput | wc -l )
ns=$( printf "%'.f" ${ns} )
echo -e "Significant_interaction\t"$ns >> ${out_file}

np=$( cat $GENOME1/${grun}_peaks.narrowPeak | wc -l )
n_peak=$( printf "%'.f" ${np} )
echo -e "Peak\t"$n_peak >> ${out_file}




