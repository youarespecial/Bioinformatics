#************************************Bioinformatics/STAR_RSEM/0_tree.txt********************************#

.
├── data
│   ├── 1_download_srr.sh
│   ├── to_tar
│   └── U87MG.txt
├── res
│   ├── fastqc
│   ├── RSEM
│   └── STAR
└── scr
    └── qc_trim_map_count.sh

U87MG.txt：
SRR388226
SRR388227
SRR388228
SRR388229

strand.sh：
strand_file=$1
strand1=`sed -n '5p' $strand_file | sed 's/.*: //'`
strand2=`sed -n '6p' $strand_file | sed 's/.*: //'`
tag=$(echo "$strand1 - $strand2" | bc )
if [ $(echo "$tag < -.4"|bc) -eq 1 ]
then
    echo "0"
elif [ $(echo "$tag > .4"|bc) -eq 1 ]
then
    echo "1"
else
    echo "0.5"
fi

#************************************Bioinformatics/STAR_RSEM/1_download_srr.sh********************************#

parallel=/home/disk/pengying/tools/parallel/bin/parallel
fastqdump=/home/disk/RNAediting_Cancer/tools/sratoolkit.2.8.2-centos_linux64/bin/fastq-dump

cat U87MG.txt | $parallel $fastqdump --split-3 -O ./
list=`cat U87MG.txt`
for srr in ${list[@]}
do
    if [ -e /home/disk/pengying/ncbi/public/sra/${srr}* ]; then rm /home/disk/pengying/ncbi/public/sra/${srr}* ; fi
done



#************************************Bioinformatics/STAR_RSEM/2_mkidx.sh********************************#
perl get_iso.pl
‘’‘
open FO, ">knownIsoforms_to_sort.txt";

while (<FI>) {
    chomp;
    next if (/^#/);

    @line = split "\t";
    $anno = $line[8];

    ($gene_id)          = $_ =~ /gene_id \"([^\"]+)\"\;/;
    ($transcript_id)    = $_ =~ /transcript_id \"([^\"]+)\"\;/;

    if ($transcript_id) {
        print FO "$gene_id\t$transcript_id\n";
    }
}

system "sort -V knownIsoforms_to_sort.txt | uniq > knownIsoforms.txt";

close FI;
close FO;
’‘’



time STAR \
--runThreadN 35 \
--runMode genomeGenerate \
--genomeDir /local/pengying/data/anno/STAR \
--genomeFastaFiles /local/pengying/data/anno/hg38.fa \
--sjdbGTFfile /local/pengying/data/anno/Homo_sapiens.GRCh38.87.chr.gtf \
--sjdbOverhang 100


time /local/pengying/tools/RSEM-1.3.1/rsem-prepare-reference \
--gtf /local/pengying/data/anno/Homo_sapiens.GRCh38.87.chr.gtf \
--transcript-to-gene-map /local/pengying/data/anno/knownIsoforms.txt \
--star \
--star-path /local/pengying/tools/STAR-2.6.0a/bin/Linux_x86_64 \
-p 35 \
/local/pengying/data/anno/hg38.fa \
/local/pengying/data/anno/RSEM/hg10

#************************************Bioinformatics/STAR_RSEM/3_qc_trim_map_count.sh********************************#

## ******************** setting ********************
cpu=6
fastqc=/home/disk/pengying/tools/FastQC/fastqc
trimmomatic=/home/disk/pengying/tools/Trimmomatic-0.38/trimmomatic-0.38.jar
STAR=/home/disk/pengying/tools/STAR-2.6.0a/bin/Linux_x86_64_static/STAR
RSEM=/home/disk/pengying/tools/RSEM-1.3.1/rsem-calculate-expression
fastq_phred=/home/disk/pengying/bin/fastq_phred.pl
infer_experiment=/home/disk/pengying/tools/RSeQC-2.6.5/scripts/infer_experiment.py
strand_test=/home/disk/pengying/rna_editing/scr/strand.sh

STAR_index=/home/genomewide/RNA-seq_idx/hg38/STAR
RSEM_index=/home/genomewide/RNA-seq_idx/hg38/RSEM/hg38
RefSeq=/home/genomewide/RNA-seq_idx/hg38/hg38_RefSeq.bed
fastqc_res=/home/disk/pengying/project/rna-seq/star/res/fastqc
STAR_res=/home/disk/pengying/project/rna-seq/star/res/STAR
RSEM_res=/home/disk/pengying/project/rna-seq/star/res/RSEM
log_file=/home/disk/pengying/project/rna-seq/star/res/quantity_log.txt

cd /home/disk/pengying/project/rna-seq/star/data

## ******************** doing ********************
list=`ls *fastq`
for data in ${list[@]}
do
    if [[ $data =~ .*_1\.fastq ]]; then
        srr=${data%%_*}
        read1=${srr}_1.fastq
        read2=${srr}_2.fastq

        echo "$srr analysis start at "`date` >> $log_file
        $fastqc $read1 -o $fastqc_res -t $cpu
        $fastqc $read2 -o $fastqc_res -t $cpu
        unzip $fastqc_res/${data%%.fastq}_fastqc.zip -d $fastqc_res

        phred=`$fastq_phred $read1`
        echo -e "\tphred: $phred" >> $log_file
        headcrop=`grep "Per base sequence content" $fastqc_res/${srr}_1_fastqc/summary.txt | cut -f 1`
        echo -e "\theadcrop: $headcrop" >> $log_file
        if [ $headcrop = FAIL ] || [ $headcrop = WARN ]
        then
            java -jar $trimmomatic PE -phred$phred $read1 $read2 ${read1}.map ${read1}.unmap ${read2}.map ${read2}.unmap HEADCROP:12 SLIDINGWINDOW:5:20
        else
            java -jar $trimmomatic PE -phred$phred $read1 $read2 ${read1}.map ${read1}.unmap ${read2}.map ${read2}.unmap SLIDINGWINDOW:5:20
        fi
        
        mkdir $STAR_res/$srr/
        $STAR --runThreadN $cpu --twopassMode Basic --outSAMstrandField intronMotif --genomeDir $STAR_index --readFilesIn ${read1}.map ${read2}.map --outFileNamePrefix $STAR_res/$srr/ --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts TranscriptomeSAM

        $infer_experiment -i $STAR_res/$srr/Aligned.sortedByCoord.out.bam -r $RefSeq > $STAR_res/$srr/strand.txt
        strand=`sh $strand_test $STAR_res/$srr/strand.txt`
        echo -e "\tstrand: $strand" >> $log_file
        $RSEM -p $cpu --bam --paired-end --forward-prob $strand $STAR_res/$srr/Aligned.toTranscriptome.out.bam $RSEM_index $RSEM_res/$srr

        mv ${srr}*fastq to_tar
        rm -r ${srr}*map $fastqc_res/${srr}*fastqc $RSEM_res/${srr}.transcript.bam $RSEM_res/${srr}.stat $STAR_res/$srr
        echo "$srr done at "`date` >> $log_file
        echo >> $log_file
        
    elif [[ $data =~ [A-Za-z]+[0-9]+\.fastq ]]; then
        srr=${data%%.*}
        read1=${data}
        "$srr analysis start at "`date` >> $time_log
        $fastqc $data -o $fastqc_res -t $cpu
        unzip $fastqc_res/${data%%.fastq}_fastqc.zip -d $fastqc_res

        phred=`$fastq_phred.pl $read1`
        echo -e "\tphred: $phred" >> $log_file
        headcrop=`grep "Per base sequence content" $fastqc_res/${srr}_fastqc/summary.txt | cut -f 1`
        echo -e "\theadcrop: $headcrop" >> $log_file
        if [ $headcrop = FAIL ] || [ $headcrop = WARN ]
        then
            java -jar $trimmomatic SE -phred$phred $read1 ${read1}.map HEADCROP:12 SLIDINGWINDOW:5:20
        else
            java -jar $trimmomatic SE -phred$phred $read1 ${read1}.map SLIDINGWINDOW:5:20
        fi

        mkdir $STAR_res/$srr/
        $STAR --runThreadN $cpu --twopassMode Basic --outSAMstrandField intronMotif --genomeDir $STAR_index --readFilesIn ${read1}.map --outFileNamePrefix $STAR_res/$srr/ --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts TranscriptomeSAM

        $infer_experiment -i $STAR_res/$srr/Aligned.sortedByCoord.out.bam -r $RefSeq > $STAR_res/$srr/strand.txt
        strand=`sh $strand_test $STAR_res/$srr/strand.txt`
        echo -e "\tstrand: $strand" >> $log_file
        $RSEM -p $cpu --bam --forward-prob $strand $STAR_res/$srr/Aligned.toTranscriptome.out.bam $RSEM_index $RSEM_res/$srr

        mv ${srr}*fastq to_tar
        rm -r ${srr}*map $fastqc_res/${srr}*fastqc $RSEM_res/${srr}.transcript.bam $RSEM_res/${srr}.stat $STAR_res/$srr
        echo "$srr done at "`date` >> $log_file
        echo >> $log_file
    fi
done












#************************************Bioinformatics/STAR_RSEM/4_diffExp.R****************************************#

library(DESeq)

## get dataframe
myfile <- Sys.glob("*.genes.results")
mydata <- c()
mysrr <- c()

for (file in myfile) {
  curdata <- read.table(file, header = TRUE)$expected_count
  mydata <- cbind(mydata, curdata)
  srr <- substr(file, 1, 9)
  mysrr <- c(mysrr, srr)
}
colnames(mydata) <- mysrr
rownames(mydata) <- read.table(file, header = TRUE)$gene_id

## create CountDataSet
type <- factor(c(rep("wt",2), rep("kd",2)), levels = c("wt", "kd"))
database <- round(as.matrix(mydata))
cds <- newCountDataSet(database, type)

## calling differential expression
# ***** 1. Standard comparison between two experimental conditions *****
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
res <- nbinomTest(cds, "wt", "kd")

# ***** 2. Working partially without replicates *****
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
res <- nbinomTest(cds, "wt", "kd")

# ***** 3. Working without any replicates *****
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds, method="blind", fitType = "local", sharingMode="fit-only")
res <- nbinomTest(cds, "wt", "kd")

## summary and output
res <- res[order(res$padj),]
sum(res$padj <= 0.01, na.rm=TRUE)
write.csv(res, file="U87MG_wt_vs_kd_DESeq.csv")

## plot heatmap
myres <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 1 & (!is.na(res$padj)),]
myres <- myres[order(myres$padj),]
this_data <- cbind(myres$baseMeanA, myres$baseMeanB)
colnames(this_data) <- c('wt', 'kd')
rownames(this_data) <- myres$id
pheatmap(this_data, cluster_row = T, cluster_col = FALSE,
         cellwidth = 20, cellheight= 10, border_color = 'grey',
         main = "DiffExp Gene",
         cex.main=1.2)






#************************************RNA_Seq/DiffEXP/DESeq.R************************************#

library(DESeq)

## get dataframe
myfile <- Sys.glob("*.genes.results")
mydata <- c()
mysrr <- c()

for (file in myfile) {
  curdata <- read.table(file, header = TRUE)$expected_count
  mydata <- cbind(mydata, curdata)
  srr <- substr(file, 1, 9)
  mysrr <- c(mysrr, srr)
}
colnames(mydata) <- mysrr
rownames(mydata) <- read.table(file, header = TRUE)$gene_id

## create CountDataSet
type <- factor(c(rep("wt",2), rep("kd",2)))
database <- round(as.matrix(mydata))
cds <- newCountDataSet(database, type)

## calling differential expression
# ***** 1. Standard comparison between two experimental conditions *****
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
res <- nbinomTest(cds, "wt", "kd")

# ***** 2. Working partially without replicates *****
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
res <- nbinomTest(cds, "wt", "kd")

# ***** 3. Working without any replicates *****
cds <- estimateDispersions(cds, method="blind", sharingMode="fit-only")
res <- nbinomTest(cds, "wt", "kd")

## summary and output
res <- res[order(res$padj),]
sum(res$padj <= 0.01, na.rm=TRUE)
write.csv(res, file="U87MG_wt_vs_kd_DESeq.csv")

#************************************RNA_Seq/DiffEXP/DESeq2.R************************************#

library(DESeq2)
setwd("/home/fyh/Desktop/RNA_edit/RSEM_res/")
myfile<-Sys.glob("*.genes.results")
mydata<-c()
mysrr<-c()
for (file in myfile) {
    curdata<-read.table(file,header = TRUE)$expected_count
    #get count info
    mydata<-cbind(mydata,curdata)
    #get srr name info
    srr<-substr(file,1,9)
    mysrr<-c(mysrr,srr)
}
colnames(mydata)<-mysrr
rownames(mydata)<-read.table(file,header = TRUE)$gene_id
type<-factor(c(rep("wt",2),rep("kd",2)),levels = c("wt","kd"))
database<-round(as.matrix(mydata))
coldata<-data.frame(row.names = colnames(database),type)
dds<-DESeqDataSetFromMatrix(database,coldata,design = ~type)
dds<-DESeq(dds)
res = results(dds, contrast=c("type", "wt", "kd"))
res=res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj <= 0.01)
#head(diff_gene_deseq2)
write.csv(diff_gene_deseq2,file= "DEG_treat_vs_control.csv")
print("output to DEG_treat_vs_control.csv")



#************************************RNA_Seq/DiffEXP/cqn_edgeR/cqn.R************************************#
#import the package
library(cqn)
library(scales)

#import the data

#1.Count Matrix, colnames for normal&disease sample, rownames for gene
data("montgomery.subset")
#2.total mapped reads for every sample, make sure in the same order as colnames
data("sizeFactors.subset")
#3.gene length&GC content for every gene, make sure in the same order as rownames
data("uCovar")

#build the cqn object
cqn.subset=cqn(montgomery.subset,lengths = uCovar$length,x=uCovar$gccontent,sizeFactors = sizeFactors.subset,verbose = TRUE)

#get the Normalization value
#log2(RPM)=s(x)+s(log2(length)), x usully for GC content and length for gene length 
RPKM.cqn=cqn.subset$y+cqn.subset$offset

#import Normalization value into edgeR, this is the workflow for reference
library(edgeR)

#build the DGElList object
grp1=c("NA06985","NA06994","NA07037","NA10847","NA11920")
grp2=c("NA11918","NA11931","NA12003","NA12006","NA12287")
d.mont=DGEList(counts = montgomery.subset,lib.size = sizeFactors.subset,group = rep(c("grp1","grp2"),each=5),genes = uCovar)

#build the design matrix
design=model.matrix(~ d.mont$samples$group)
d.mont$offset=cqn.subset$glm.offset
d.mont.cqn=estimateGLMCommonDisp(d.mont,design = design)

#fit the glm model
efit.cqn=glmFit(d.mont.cqn,design = design)
elrt.cqn=glmLRT(efit.cqn,coef = 2)

#get the top DE gene
topTags(elrt.cqn)
summary(decideTestsDGE(elrt.cqn))


