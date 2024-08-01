# ChIAIR-PIPE
Program for ChAIR data analysis

![](img/pipeline.png)

# REQUIREMENT
```
ChIAPET-pipeline (https://github.com/TheJacksonLaboratory/ChIA-PIPE)
BEDTools
SAMTools
DeepTools
Picard
MACS2
Cell-ranger
Cell-ranger-atac
Cell-ranger-arc

Downstream analysis:
Seurat
Signac
DoubletFinder
Harmony
Monocle
BandNorm
Dip-c
Nuc
CtG

Python
R
Perl

Slurm
```
Some of the software we use is the singularity version. The following four software can be downloaded from Baidu.com (https://pan.baidu.com/s/1KKC0DJ43wuuCdPwr09dgCg?pwd=bj2p   Password: bj2p)
```
cellranger_6.0.0.sif
cpu0.0.1a-r2.sumner.sif
GATK_latest.sif
juicer_1.22.01.sif
```
```
cp cellranger_6.0.0.sif cpu0.0.1a-r2.sumner.sif GATK_latest.sif juicer_1.22.01.sif singularity
```

# USAGE
```
1. Download data from CNCB (https://www.cncb.ac.cn/)
2. Changing the file name to match cellranger's requirements
3. bash submit_cpu10x_for_scChiatac.mm10.nobl.sh SCG0192_GT22-15872_SI-NA-D6.fastq.prefix > SCG0192_GT22-15872_SI-NA-D6.fastq.prefix.log 2>SCG0192_GT22-15872_SI-NA-D6.fastq.prefix.log.run.sh
```

# CONTACT
黄星宇 (Xingyu Huang, xingyu.huang@zju.edu.cn/huang182@live.cn)