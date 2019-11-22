export PATH=$PATH:/home/zhangfeng/disk/project/Cuffdiff/tools/cufflinks-2.2.1.Linux_x86_64
export PATH=$PATH:/home/zhangfeng/disk/project/Cuffdiff/tools/bowtie2
export PATH=$PATH:/home/zhangfeng/disk/project/Cuffdiff/tools/tophat-2.1.0.Linux_x86_64
export PATH=$PATH:/home/zhangfeng/disk/project/Cuffdiff/tools/samtools-1.2

cuffdiff='/home/zhangfeng/disk/project/Cuffdiff/tools/cufflinks-2.2.1.Linux_x86_64/cuffdiff'
cuffmerge='/home/zhangfeng/disk/project/Cuffdiff/tools/cufflinks-2.2.1.Linux_x86_64/cuffmerge'
cufflinks='/home/zhangfeng/disk/project/Cuffdiff/tools/cufflinks-2.2.1.Linux_x86_64/cufflinks'
pc_gtf='/home/zhangfeng/disk/project/data/annotation/mm9.gtf'
downloaddir='/home/zhangfeng/disk/project/data/download/'

human_gtf='/home/genomewide/annotation/hg19/Homo_sapiens.GRCh37.75.chr.gtf'
read1=/home/zhangfeng/disk/project/data/HXH_SEQ/HXH_lincRNA/AC_1.fastq
read2=/home/zhangfeng/disk/project/data/HXH_SEQ/HXH_lincRNA/AC_2.fastq


tophat2 -p 12 -o $read1\.tophatdir ../index/hg19 $read1 $read2
$cufflinks -p 8 -u -G $human_gtf -o $read1\.cufflinks $read1\.tophatdir/accepted_hits.bam

read1=/home/zhangfeng/disk/project/data/HXH_SEQ/HXH_lincRNA/AS_1.fastq
read2=/home/zhangfeng/disk/project/data/HXH_SEQ/HXH_lincRNA/AS_2.fastq
tophat2 -p 8 -o $read1\.tophatdir ../index/hg19 $read1 $read2
$cufflinks -p 8 -u -G $human_gtf -o $read1\.cufflinks $read1\.tophatdir/accepted_hits.bam

$cuffmerge -g $human_gtf -s ../index/hg19.fa -p 10  -o /home/zhangfeng/disk/project/data/HXH_SEQ/Cuffdiff/AC_AS/GTF  /home/zhangfeng/disk/project/data/HXH_SEQ/Cuffdiff/AC_AS/assemblies.txt

AC='/home/zhangfeng/disk/project/data/HXH_SEQ/HXH_lincRNA/AC_1.fastq.tophatdir/accepted_hits.bam'
AS='/home/zhangfeng/disk/project/data/HXH_SEQ/HXH_lincRNA/AS_1.fastq.tophatdir/accepted_hits.bam'
$cuffdiff -o /home/zhangfeng/disk/project/data/HXH_SEQ/Cuffdiff/AC_AS/Diff -p 10 -L AC,AS -u /home/zhangfeng/disk/project/data/HXH_SEQ/Cuffdiff/AC_AS/GTF/merged.gtf $AC $AS


###############*********************************************************************************************###############

mm10_index=/home/genomewide/refgenome/mm10/mm10
mm10_seq=/home/genomewide/refgenome/mm10/mm10.fa
mm10_gtf=/home/genomewide/annotation/mm10/Mus_musculus.GRCm38.87.chr.gtf
fastqc=/home/zhangfeng/tools/FastQC/fastqc
tophat2=/home/zhangfeng/disk/project/Cuffdiff/tools/tophat-2.1.0.Linux_x86_64/tophat2
cpu=8

wp=/home/yzj/ZF/0416/YN0301b/
`mkdir /home/yzj/ZF/0416/YN0301b/fastqc`
line=L2
read1=$wp'V300016995_L2_HK500HUMvtrTAAGRAAPEI-591_1.fq.gz'
read2=$wp'V300016995_L2_HK500HUMvtrTAAGRAAPEI-591_2.fq.gz'
$fastqc $read1 $read2 -o $wp'/fastqc'
$tophat2 -p $cpu -o $wp$line\.tophatdir $mm10_index $read1 $read2
 
`mkdir /home/yzj/ZF/0416/YN0301b/cufflinks`
cufflinks -p $cpu -u -G $mm10_gtf -o ${wp}cufflinks ${wp}$line\.tophatdir/accepted_hits.bam

mm10_index=/home/genomewide/refgenome/mm10/mm10
mm10_seq=/home/genomewide/refgenome/mm10/mm10.fa
mm10_gtf=/home/genomewide/annotation/mm10/Mus_musculus.GRCm38.87.chr.gtf

#realpath */cufflinks/transcripts.gtf > GTF/assemblies.txt

wp=/home/yzj/ZF/0416/
cuffmerge -g $mm10_gtf -s $mm10_seq -p 35 -o ${wp}GTF ${wp}GTF/assemblies.txt

N1=${wp}YN0117/L2.tophatdir/accepted_hits.bam
N2=${wp}YN0220/L2.tophatdir/accepted_hits.bam
N3=${wp}YN0226/L2.tophatdir/accepted_hits.bam
N4=${wp}YN0301a/L2.tophatdir/accepted_hits.bam
N5=${wp}YN0301b/L2.tophatdir/accepted_hits.bam

T1=${wp}T4S0222/L2.tophatdir/accepted_hits.bam
T2=${wp}T4S0222a/L2.tophatdir/accepted_hits.bam
T3=${wp}T4S0226b/L2.tophatdir/accepted_hits.bam
T4=${wp}TZ0114/L2.tophatdir/accepted_hits.bam
T5=${wp}TZ0117/L2.tophatdir/accepted_hits.bam

merged_gtf=${wp}GTF/merged.gtf
cuffdiff -o ${wp}DIFF_N_T -p 35 -L N,T -u $merged_gtf $N1,$N2,$N3,$N4,$N5 $T1,$T2,$T3,$T4,$T5
