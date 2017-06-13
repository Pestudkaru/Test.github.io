# RNA-seq on our data
blsd
library("GenomicFeatures")
library("GenomicAlignments")
library("DESeq2")
library("Rsamtools")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("genefilter")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("data.table")
library("WriteXLS")
library("reshape")

fin_info <- read.csv("/media/gunnar/Data/RNA-Seq_Experiments/Finni_data/Finni_samples.csv",header = T)
filenames <- file.path("/media/gunnar/Data/RNA-Seq_Experiments/Finni_data/05_Sorted_and_Indexed/", fin_info$bamfile)
file.exists(filenames)
bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])
gtffile <- ("/media/gunnar/Data/Genomes/Hooman/Homo_sapiens.GRCh38.85.gtf")
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))
(ebg <- exonsBy(txdb, by="gene"))

fin_se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )
#fin_se
#assayNames(fin_se)
#head(assay(fin_se), 3)
#colSums(assay(fin_se))
#rowRanges(fin_se)
#str(metadata(rowRanges(fin_se)))
#colData(fin_se)

(colData(fin_se) <- DataFrame(fin_info))
#sum_overlaps$Cell_line
#sum_overlaps$Treatment1
#sum_overlaps$Treatment2
#sum_overlaps$Sample_no
#How many cell/lines do i have with more than 1mil reads (also how many reads)
  Finni_dds <- DESeqDataSet(fin_se, design = ~ treatment2+)
nrow(Finni_dds)
Finni_dds <- Finni_dds[ rowSums(counts(Finni_dds)) > 16, ]
nrow(Finni_dds)
Finni_Deseq_dds <- DESeq(Finni_dds)

save(Finni_Deseq_dds,fin_se,file = "~/R/NGS R Scripts/Fionan/fin's_DESEQ_2.Rdata")
load("~/R/NGS R Scripts/Fionan/fin's_DESEQ_2.Rdata")

#Gene expression count and visualizing counts
countdata <- assay(fin_se)
nrow(countdata)
head(countdata)
norm_countdata <- counts(Finni_Deseq_dds,normalized=T)
nrow(norm_countdata)
### raw counts
#barplot(colSums(countdata),names.arg = paste(fin_info$name))
counts_data <- (colSums((countdata)/ 1e6, 1))
counts_data=as.data.frame(counts_data)
counts_data$names=(fin_info$name)
read_count_plot=ggplot(counts_data,aes(x=names,y=counts_data))
read_count_plot+geom_bar(stat="identity")+theme(axis.text.x=element_text(size=14, angle=45,vjust = 1,hjust = 1))

#getting the normalized count data
#barplot(colSums(norm_count),names.arg = paste(fin_info$name))
norm_counts_data=(colSums((norm_countdata)/ 1e6, 1))
norm_counts_data=as.data.frame(norm_counts_data)
norm_counts_data$names=(fin_info$name)
norm_counts_data_plot=ggplot(norm_counts_data,aes(x=names,y=norm_counts_data))
norm_counts_data_plot+geom_bar(stat="identity")+theme(axis.text.x=element_text(size=14, angle=45,vjust = 1,hjust = 1))
#A useful first step in an RNA-seq analysis is often to assess overall similarity between samples: Which samples are similar to each other, which are different? Does this fit to the expectation from the experiment???s design?
fin_rld <- rlog(Finni_Deseq_dds, blind=T)
sampleDists <- dist( t( assay(fin_rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
sampleDistMatrix1=sampleDistMatrix
rownames(sampleDistMatrix) <- paste(fin_rld$name)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
#CONTRAST can only be used for columns that were used to create the DDS DataFrame
#Log Fold Change is the effect size estimate. It tells us how much the gene xxxs expression seems to have changed due to treatment.
#This value is reported on a logarithmic scale to base 2: for example, a log2 fold change of 1.5 means that the gene???s expression is increased by a multiplicative factor of 2^1.5???2.82.
resultsNames(Finni_Deseq_dds)
Finni_Deseq_dds$treatment2
tulemused_U_L=results(Finni_Deseq_dds, contrast=c("treatment2","WT_LPS14h","KO_LPS14h")) 
mcols(tulemused_U_L, use.names=TRUE)
topGene <- rownames(tulemused_U_L)[which.min(tulemused_U_L$padj)]
#plotCounts(Finni_Deseq_dds, gene=topGene, intgroup=("treatment2"))
######## Single sig gene
data <- plotCounts(Finni_Deseq_dds, gene=topGene, intgroup=c("treatment2"), returnData=TRUE)
ggplot(data, aes(x=treatment2, y=count, fill=treatment2)) +
scale_y_log10() + geom_dotplot(binaxis="y", stackdir="center")+ ggtitle(topGene)+theme(axis.text.x=element_blank())+ labs(x = "Treatments")
# MA
plotMA(tulemused_U_L, ylim=c(-5,5))
with(tulemused_U_L[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")})
#How many sig genes at all
sum(tulemused_U_L$padj < 0.1, na.rm=TRUE)
#comparing 2 groups to each other.
par( mfrow = c( 1, 2 ) )
plot(log2(counts(Finni_Deseq_dds, normalized=TRUE)[,c(1,5)] + 1),
     pch=16, cex=0.3)
plot(assay(fin_rld)[,c(1,5)],
     pch=16, cex=0.3)
#PCA#
PCAdata <- plotPCA(fin_rld, intgroup = c( "treatment2"), returnData=TRUE)
percentVar <- round(100 * attr(PCAdata, "percentVar"))
PCAdata=cbind(PCAdata,fin_info$name,fin_info$treatment)
ggplot(PCAdata, aes(PC1, PC2, color=group,shape=fin_info$treatment)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

#####################################################################################
#                                                                                   #
#                                                                                   #
#                               UP- regulated                                       #
#                                                                                   #
#                                                                                   #
#####################################################################################
###DIfferentially upregulated genes, because of LPS, in 14H exp
Finni_Deseq_dds$treatment2
WT_14h=results(Finni_Deseq_dds,contrast = c("treatment2","WT_unst14h","WT_LPS14h"))
WT_14h_DT=as.data.table(WT_14h)
WT_14h_DT$ensemble_id=rownames(WT_14h)
#Since Deseq gives BH multiple testing correction as padj then lets consider a fraction of 10% false positives as acceptable. So anything with adjusted p value below 10% = 0.1 are significant. 
WT_14h_Sig=subset(WT_14h_DT, padj < 0.1)
#and again - everything with a log fold change of 1.5
WT_14h_Sig_Upreg=subset(WT_14h_Sig,log2FoldChange>0)
###################
Finni_Deseq_dds$treatment2
KO_14h=results(Finni_Deseq_dds,contrast = c("treatment2","KO_unst14h","KO_LPS14h"))
KO_14h_DT=as.data.table(KO_14h)
KO_14h_DT$ensemble_id=rownames(KO_14h)
#again 10%fdr
KO_14h_Sig=subset(KO_14h_DT, padj < 0.1)
#and again - everything with a log fold change of 1.5
KO_14h_Sig_Upreg=subset(KO_14h_DT,log2FoldChange>1.5)
#Genes common in both
Common_UP_Genes=WT_14h_Sig_Upreg[ensemble_id %in% KO_14h_Sig_Upreg$ensemble_id]
#Genes that are missing from KO
Missing_in_KO=subset(WT_14h_Sig_Upreg,!(ensemble_id%in%KO_14h_Sig_Upreg$ensemble_id))
#Genes that are missing from WT
Missing_in_WT=subset(KO_14h_Sig_Upreg,!(ensemble_id%in%WT_14h_Sig_Upreg$ensemble_id))
#TEST
#Test=rbind(WT_14h_Sig_Upreg,KO_14h_Sig_Upreg)
#Test=unique(Test$ensemble_id)
WriteXLS(Common_UP_Genes,"~/Desktop/up_common_14h.xlsx")
WriteXLS(Missing_in_KO,"~/Desktop/up_MIA_KO_14.xlsx")
WriteXLS(Missing_in_WT,"~/Desktop/up_MIA_WT_14.xlsx")
#####################################################################################
#                                                                                   #
#                                                                                   #
#                               Downregulated                                       #
#                                                                                   #
#                                                                                   #
#####################################################################################
###DIfferentially downregulated genes, because of LPS, in 14H exp
WT_14h_Sig=subset(WT_14h_DT, padj < 0.1)
#and again - everything with a log fold change of 1.5
WT_14h_Sig_Down=subset(WT_14h_Sig,log2FoldChange<(-1.5))
###################
#again 10%fdr
KO_14h_Sig=subset(KO_14h_DT, padj < 0.1)
#and again - everything with a log fold change of 1.5
KO_14h_Sig_Down=subset(KO_14h_DT,log2FoldChange<(-1.5))
#Genes common in both
Common_Down_Genes=WT_14h_Sig_Down[ensemble_id %in% KO_14h_Sig_Down$ensemble_id]
#Genes that are missing from KO
Missing_in_KO_down=subset(WT_14h_Sig_Down,!(ensemble_id%in%KO_14h_Sig_Down$ensemble_id))
#Genes that are missing from WT
Missing_in_WT_down=subset(KO_14h_Sig_Down,!(ensemble_id%in%WT_14h_Sig_Down$ensemble_id))
#TEST
#Test=rbind(WT_14h_Sig_Down,KO_14h_Sig_Down)
#Test=unique(Test$ensemble_id)
#length(Test)==nrow(rbind(Common_Down_Genes,Missing_in_KO_down,Missing_in_WT_down))
WriteXLS(Common_Down_Genes,"~/Desktop/down_common_14h.xlsx")
WriteXLS(Missing_in_KO_down,"~/Desktop/down_MIA_KO_14.xlsx")
WriteXLS(Missing_in_WT_down,"~/Desktop/down_MIA_WT_14.xlsx")

#####################################################################################
#                                                                                   #
#                                     2H                                            #
#                                                                                   #
#####################################################################################
#                                                                                   #
#                                                                                   #
#                               UP- regulated                                       #
#                                                                                   #
#                                                                                   #
#####################################################################################
###DIfferentially upregulated genes, because of LPS, in 14H exp
Finni_Deseq_dds$treatment2
WT_2h=results(Finni_Deseq_dds,contrast = c("treatment2","WT_unst2h","WT_LPS2h"))
WT_2h_DT=as.data.table(WT_2h)
WT_2h_DT$ensemble_id=rownames(WT_2h)
#Since Deseq gives BH multiple testing correction as padj then lets consider a fraction of 10% false positives as acceptable. So anything with adjusted p value below 10% = 0.1 are significant. 
WT_2h_Sig=subset(WT_2h_DT, padj < 0.1)
#and again - everything with a log fold change of 1.5
WT_2h_Sig_Upreg=subset(WT_2h_Sig,log2FoldChange>1)
###################
KO_2h=results(Finni_Deseq_dds,contrast = c("treatment2","KO_unst2h","KO_LPS2h"))
KO_2h_DT=as.data.table(KO_2h)
KO_2h_DT$ensemble_id=rownames(KO_2h)
#again 10%fdr
KO_2h_Sig=subset(KO_2h_DT, padj < 0.1)
#and again - everything with a log fold change of 1.5
KO_2h_Sig_Upreg=subset(KO_2h_DT,log2FoldChange>1)
#Genes common in both
Common_UP_Genes_2H=WT_2h_Sig_Upreg[ensemble_id %in% KO_2h_Sig_Upreg$ensemble_id]
nrow(Common_UP_Genes_2H)
#Genes that are missing from KO
Missing_in_KO_2H=subset(WT_2h_Sig_Upreg,!(ensemble_id%in%KO_2h_Sig_Upreg$ensemble_id))
nrow(Missing_in_KO_2H)
#Genes that are missing from WT
Missing_in_WT_2H=subset(KO_2h_Sig_Upreg,!(ensemble_id%in%WT_2h_Sig_Upreg$ensemble_id))
nrow(Missing_in_WT_2H)
#TEST
#Test=rbind(WT_2h_Sig_Upreg,KO_2h_Sig_Upreg)
#Test=unique(Test$ensemble_id)
#length(Test)==nrow(rbind(Common_UP_Genes_2H,Missing_in_KO_2H,Missing_in_WT_2H))
WriteXLS(Common_UP_Genes_2H,"~/Desktop/up_common_2h.xlsx")
WriteXLS(Missing_in_KO_2H,"~/Desktop/up_MIA_KO_2.xlsx")
WriteXLS(Missing_in_WT_2H,"~/Desktop/up_MIA_WT_2.xlsx")
#####################################################################################
#                                                                                   #
#                                     2H                                            #
#                               Downregulated                                       #
#                                                                                   #
#                                                                                   #
#####################################################################################
###DIfferentially downregulated genes, because of LPS, in 14H exp
WT_2h_Sig=subset(WT_2h_DT, padj < 0.1)
#and again - everything with a log fold change of 1.5
WT_2h_Sig_Down=subset(WT_2h_Sig,log2FoldChange<(-1.5))
###################
#again 10%fdr
KO_2h_Sig=subset(KO_2h_DT, padj < 0.1)
#and again - everything with a log fold change of 1.5
KO_2h_Sig_Down=subset(KO_2h_DT,log2FoldChange<(-1.5))
#Genes common in both
Common_Down_Genes_2H=WT_2h_Sig_Down[ensemble_id %in% KO_2h_Sig_Down$ensemble_id]
nrow(Common_Down_Genes_2H)
#Genes that are missing from KO
Missing_in_KO_down_2H=subset(WT_2h_Sig_Down,!(ensemble_id%in%KO_2h_Sig_Down$ensemble_id))
nrow(Missing_in_KO_down_2H)
#Genes that are missing from WT
Missing_in_WT_down_2H=subset(KO_2h_Sig_Down,!(ensemble_id%in%WT_2h_Sig_Down$ensemble_id))
nrow(Missing_in_WT_down_2H)
#TEST
#Test=rbind(WT_2h_Sig_Down,KO_2h_Sig_Down)
#Test=unique(Test$ensemble_id)
#length(Test)==nrow(rbind(Common_Down_Genes_2H,Missing_in_KO_down_2H,Missing_in_WT_down_2H))
WriteXLS(Common_Down_Genes_2H,"~/Desktop/down_common_2h.xlsx")
WriteXLS(Missing_in_KO_down_2H,"~/Desktop/down_MIA_KO_2h.xlsx")
WriteXLS(Missing_in_WT_down_2H,"~/Desktop/down_MIA_WT_2h.xlsx")

################### ONE MORE LIST TO SEE UPreg genes missing in KO and also missing in 2H



#    Main list- missing in KO 14H
Interesting_list=subset(Missing_in_KO,!(ensemble_id%in%Second_LPS_genes_UP$ensemble_id))
Interesting_list=subset(Interesting_list,!(ensemble_id%in%First_LPS_genes_UP$ensemble_id))
Interesting_list=subset(Interesting_list,!(ensemble_id%in%WT_2h_Sig_Upreg$ensemble_id))
Interesting_list=subset(Interesting_list,!(ensemble_id%in%KO_2h_Sig_Upreg$ensemble_id))
Interesting_list

Interesting_list$symbol <- mapIds(org.Hs.eg.db,
                       keys=Interesting_list$ensemble_id,
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
WriteXLS(Interesting_list,"~/Desktop/Interesting_list.xlsx")

#####################################################################################
#                                                                                   #
#                                                                                   #
#                               Veit's way of doing it                              #
#                                                                                   #
#                                                                                   #
#####################################################################################
Finni_Deseq_dds$treatment2
resultsNames(Finni_Deseq_dds)
data1=results(Finni_Deseq_dds,contrast = c ("treatment2","WT_LPS14h","KO_LPS14h"))
data1_DT=as.data.table(data1)
data1_DT$ensemble_id=rownames(data1)
data1_DT_sig=subset(data1_DT, padj < 0.1)
##    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ##
##    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ##
data2=results(Finni_Deseq_dds,contrast = c ("treatment2","WT_unst14h","KO_unst14h"))
data2_DT=as.data.table(data2)
data2_DT$ensemble_id=rownames(data2)
data2_DT_sig=subset(data2_DT, padj < 0.1)

controlled_for=subset(data1_DT_sig,!(ensemble_id %in% data2_DT_sig$ensemble_id))
Common_test=data1_DT_sig[ensemble_id %in% data2_DT_sig$ensemble_id]

Fin_countdata <- counts(Finni_Deseq_dds,normalized=T)
head(Fin_countdata)
Fin_countdata_dt=as.data.table(Fin_countdata)
head(class(Fin_countdata_dt$V1))
colnames(Fin_countdata_dt)=as.character(fin_info$name)
head(Fin_countdata_dt)
#Fin_countdata_dt=Fin_countdata_dt[,.(KO_1_LPS14h,KO_2_LPS14h,WT_1_LPS14h,WT_2_LPS14h)]
Fin_countdata_dt$ensembl_id=rownames(Fin_countdata)
#data1_DT_sig=data1_DT_sig[ order(data1_DT_sig$log2FoldChange,decreasing = T), ]
all=controlled_for$ensemble_id
new_list=c(all,"ENSG00000162711","ENSG00000136244","ENSG00000125538","ENSG00000232810","ENSG00000169245","ENSG00000185745")
# NLRP3, IL6, IL1b, TNF, CXCL10, IFIT1
interesting=Fin_countdata_dt[ensembl_id %in% new_list]
interesting=as.data.frame(interesting)
head(class(interesting$KO_2_unst2h))
class(interesting)
interesting3=t(interesting[1:14,1:16,])
head(class(interesting3[,9]))
interesting2=as.data.frame(interesting3)
interesting2$Sample_names=rownames(interesting3)
interesting4=interesting2
colnames(interesting4)=new_list
WriteXLS(interesting4,"~/Desktop/interesting.xlsx")
WriteXLS(table,"~/Desktop/table.xls")



katse=interesting2
colnames(katse)=c(new_list,"Sample")
katse=katse[1:16,1:15,]
katse= cbind(katse,table[,4:6])
gendata=aggregate(ENSG00000257335 ~ Time  + trt + Cell, FUN=mean, data=katse) 
ggplot(gendata, aes(x=Time, y=ENSG00000257335, colour=Cell)) + geom_line() + ggtitle("ENSG00000257335")

new_new=c("ENSG00000080546","ENSG00000115009","ENSG00000115414","ENSG00000137801","ENSG00000149294","ENSG00000170113","ENSG00000189058","ENSG00000257335","ENSG00000277734","ENSG00000162711","ENSG00000136244","ENSG00000125538","ENSG00000232810","ENSG00000169245","ENSG00000185745")
New_counts=Fin_countdata_dt[ensembl_id %in% new_new]
New_counts=as.data.frame(New_counts)
New_counts2=t(New_counts)
New_counts2=as.data.frame(New_counts2)
colnames(New_counts2)=new_new
table=read.csv("~/Scripts/Finnis_RNA_Seq/table.csv")
New_counts2=cbind(New_counts2,table[,2:4])
list=colnames(New_counts2)
"ENSG00000080546" "ENSG00000115009" "ENSG00000115414" "ENSG00000137801" "ENSG00000149294" "ENSG00000170113" "ENSG00000189058"
"ENSG00000257335" "ENSG00000277734" "ENSG00000162711" "ENSG00000136244" "ENSG00000125538" "ENSG00000232810" "ENSG00000169245"
"ENSG00000185745"      
jpeg("ENSG00000277734.jpg")
ggplot(New_counts2, aes(x=Time, y=ENSG00000277734, colour=Cell)) + geom_point() + 
  stat_smooth(se=FALSE,method="loess") +  scale_y_log10()+ggtitle("ENSG00000277734")
dev.off()


for (i in List){
  print(paste(i,"on äge!"))}

##    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ##
##    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ##
##    ##    ####    ##     2 Hours now  ####    ##    ####    ##    ####    ##    ####
Finni_Deseq_dds$treatment2
data3=results(Finni_Deseq_dds,contrast = c("treatment2","WT_LPS2h","KO_LPS2h"))
data3_DT=as.data.table(data3)
data3_DT$ensemble_id=rownames(data3)
data3_DT_sig=subset(data3_DT, padj < 0.1)
##    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ##
##    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ####    ##    ##
data4=results(Finni_Deseq_dds,contrast = c("treatment2","WT_unst2h","KO_unst2h"))
data4_DT=as.data.table(data4)
data4_DT$ensemble_id=rownames(data4)
data4_DT_sig=subset(data4_DT, padj < 0.1)
### ### ### ### ### 
plotMA(data1, ylim=c(-5,5))
abline(h=2)

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
fin_info <- read.csv("/media/gunnar/Data/RNA-Seq_Experiments/Finni_data/Finni_samples.csv",header = T)
filenames <- file.path("/media/gunnar/Data/RNA-Seq_Experiments/Finni_data/05_Sorted_and_Indexed/", fin_info$bamfile)
file.exists(filenames)
bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])
gtffile <- ("/media/gunnar/Data/Genomes/Hooman/Homo_sapiens.GRCh38.85.gtf")
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))
(ebg <- exonsBy(txdb, by="gene"))

fin_se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                            mode="Union",
                            singleEnd=FALSE,
                            ignore.strand=TRUE,
                            fragments=TRUE )
fin_se2=fin_se
(colData(fin_se2) <- DataFrame(fin_info))
Finni_dds2 <- DESeqDataSet(fin_se2, design = ~ treatment3)
nrow(Finni_dds2)
Finni_dds2 <- Finni_dds2[ rowSums(counts(Finni_dds2)) > 16, ]
nrow(Finni_dds2)
Finni_Deseq_dds2 <- DESeq(Finni_dds2)
save(Finni_Deseq_dds2,fin_se2,file = "~/R/NGS R Scripts/Fionan/fin's_DESEQ2_2.Rdata")
load("~/R/NGS R Scripts/Fionan/")
### all the DE genes
res=results(Finni_Deseq_dds2)
res
sig=subset(res,padj < 0.1)
sig
DE_genes=rownames(sig)  

load("~/R/NGS R Scripts/Fionan/fin's_DESEQ.Rdata")
plotDispEsts(Finni_Deseq_dds)
sizefactors=sizeFactors(Finni_Deseq_dds)
fin_se_dt=as.data.table(assay(fin_se2))
fin_se_dt2=as.data.table((counts(Finni_Deseq_dds,normalized=T)))

nrow(fin_se_dt)
head(fin_se_dt)
head(assay(fin_se2))
fin_se_dt$ensembl=rownames(fin_se2)
head(fin_se_dt)
fin_se_dt=fin_se_dt[ensembl %in% DE_genes]
head(fin_se_dt)
nrow(fin_se_dt)


countdata2=as.data.frame(fin_se_dt)
head(countdata2)
class(countdata2)
nrow(countdata2)
countdata=as.matrix(countdata2[,1:16])
head(countdata)
nrow(countdata)
rownames(countdata)=fin_se_dt$ensembl
head(countdata)
colnames(countdata)=fin_info$name
head(countdata)
coldata <- colData(fin_se2)
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = coldata,
                                  design = ~ treatment2)
nrow(ddsMat)
sizeFactors(ddsMat)=sizefactors
ddsMat_DESeq <- estimateDispersions(ddsMat_DESeq)
plotDispEsts(ddsMat_DESeq)
ddsMat_DESeq=DESeq(ddsMat)

norm_count=counts(ddsMat_DESeq,normalized=T)
norm_counts_data=(colSums((norm_count)/ 1e6, 1))
norm_counts_data=as.data.frame(norm_counts_data)
norm_counts_data$names=(fin_info$name)
norm_counts_data_plot=ggplot(norm_counts_data,aes(x=names,y=norm_counts_data))
norm_counts_data_plot+geom_bar(stat="identity")+theme(axis.text.x=element_text(size=14, angle=45,vjust = 1,hjust = 1))

# TIMEPOINT EXP
time_info=read.csv("/media/gunnar/Data/RNA-Seq_Experiments/Finni_data/timepoint.csv",stringsAsFactors = F)
time_info2=time_info[-c(3,7,11,15),]
time_info2$time=as.factor(time_info2$time)
time_info2$cell_type=as.factor(time_info2$cell_type)
time_info2$treatment=as.factor(time_info2$treatment)
class(time_info2$sample_no)
class(time_info2$cell_type)
class(time_info2$replicate)
class(time_info2$treatment)
class(time_info2$time)
class(time_info2$name)

filenames <- file.path("/media/gunnar/Data/RNA-Seq_Experiments/Finni_data/05_Sorted_and_Indexed/", time_info2$bamfile)
file.exists(filenames)
bamfiles <- BamFileList(filenames, yieldSize=2000000)
gtffile <- ("/media/gunnar/Data/Genomes/Hooman/Homo_sapiens.GRCh38.85.gtf")
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))
(ebg <- exonsBy(txdb, by="gene"))

time_se2 <- summarizeOverlaps(features=ebg, reads=bamfiles,
                            mode="Union",
                            singleEnd=FALSE,
                            ignore.strand=TRUE,
                            fragments=TRUE )


colData(time_se2) <- DataFrame(time_info2)
time_dds2 <- DESeqDataSet(time_se2, ~ cell_type + treatment)
nrow(time_dds2)
time_dds2 <- time_dds2[ rowSums(counts(time_dds2)) > 12, ]
nrow(time_dds2)
deseq_time2 <- DESeq(time_dds2)
resultsNames(deseq_time2)
save(time_se2,file="~/R/NGS R Scripts/Fionan/time_no14unstim.R")
load("~/R/NGS R Scripts/Fionan/time_no14unstim.R")
resTC <- results(deseq_time2)
head(resTC[order(resTC$padj),])

huvitavad=resTC[order(resTC$padj),]
huvit_dt=as.data.table(huvitavad)
huvit_dt$ensembl=rownames(huvitavad)
huvit_sig=huvit_dt[huvit_dt$padj<0.01]

list=huvit_sig$ensembl
list
first=rownames(resTC)=="ENSG00000102755"
second=rownames(resTC)=="ENSG00000120738"
third=rownames(resTC)=="ENSG00000132669"
fourth=rownames(resTC)=="ENSG00000170113"
sixth=ENSG00000155846
seventh=ENSG00000095951

resTC[first,]


data<- plotCounts(deseq_time2, which(rownames(resTC)=="ENSG00000170113"), 
                   intgroup=c("time","cell_type"), returnData=TRUE)
ggplot(data, aes(x=time, y=count, color=cell_type, group=cell_type)) + 
  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()+ggtitle("ENSG00000170113")


#####################################################################################
#####################################################################################
#                                                                                   #
#                                                                                   #
#                               2 variable model                                    #
#                                                                                   #
#                                                                                   #
#####################################################################################

fin_info <- read.csv("/media/gunnar/Data/RNA-Seq_Experiments/Finni_data/Finni_samples_2_variables.csv",header = T)
filenames <- file.path("/media/gunnar/Data/RNA-Seq_Experiments/Finni_data/05_Sorted_and_Indexed/", fin_info$bamfile)
file.exists(filenames)
bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])
gtffile <- ("/media/gunnar/Data/Genomes/Hooman/Homo_sapiens.GRCh38.85.gtf")
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character()))
(ebg <- exonsBy(txdb, by="gene"))

fin_se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                            mode="Union",
                            singleEnd=FALSE,
                            ignore.strand=TRUE,
                            fragments=TRUE )

(colData(fin_se) <- DataFrame(fin_info))
Finni_dds <- DESeqDataSet(fin_se, design = ~treatment3 + cell_type)
nrow(Finni_dds)
Finni_dds <- Finni_dds[ rowSums(counts(Finni_dds)) > 8, ]
nrow(Finni_dds)
Finni_Deseq_dds <- DESeq(Finni_dds)

Finni_rld <- rlog(Finni_Deseq_dds)
save(fin_se,Finni_Deseq_dds,Finni_rld,fin_info,file="/media/gunnar/Data/RNA-Seq_Experiments/Finni_data/2variable_data.Rdata")
load("/media/gunnar/Data/RNA-Seq_Experiments/Finni_data/2variable_data.Rdata")

sampleDists <- dist(t(assay(Finni_rld)))            
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( Finni_rld$treatment3, Finni_rld$cell_type, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
clustering_distance_rows=sampleDists,
clustering_distance_cols=sampleDists,
col=colors)

poisd <- PoissonDistance(t(counts(Finni_Deseq_dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste(Finni_rld$treatment3, Finni_rld$cell_type, sep="-" )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
clustering_distance_rows=poisd$dd,
clustering_distance_cols=poisd$dd,
col=colors)

plotPCA(Finni_rld, intgroup = c("treatment3", "cell_type"))

(Fres <- results(Finni_Deseq_dds))
mcols(Fres, use.names=TRUE)
summary(Fres)

Fres.dt=as.data.table(Fres)
Fres.dt$ensemble_id=rownames(Fres)
Fres.dt_sig=subset(Fres.dt, padj < 0.1)
# vector of old genes that came up in the old analysis
old_genes=c("ENSG00000080546","ENSG00000115009","ENSG00000115414","ENSG00000137801","ENSG00000149294","ENSG00000170113","ENSG00000189058","ENSG00000257335","ENSG00000277734")
# which ones are still there
intersext=res.dt_sig[ensemble_id%in%old_genes]
# so 6 are still in common: ENSG00000115009 ENSG00000115414 ENSG00000137801 ENSG00000149294 ENSG00000257335 ENSG00000277734
old_genes[!old_genes%in%res.dt_sig$ensemble_id]
# the three that are no longer considered significant: "ENSG00000080546" "ENSG00000170113" "ENSG00000189058"
missing=subset(res.dt_sig,!(ensemble_id%in%old_genes))
# and 26 new ones! 

res.dt_sig$symbol <- mapIds(org.Hs.eg.db,
                                  keys=res.dt_sig$ensemble_id,
                                  column="SYMBOL",
                                  keytype="ENSEMBL",
                                  multiVals="first")
WriteXLS(res.dt_sig,"~/Desktop/new.list.xls")


head(Fres)
Fresults=as.data.frame(Fres)
Fresults$Gene=mapIds(org.Hs.eg.db,keys=rownames(Fresults),column="SYMBOL",keytype="ENSEMBL",multiVals="first")
Fresults = mutate(Fresults, sig=ifelse(Fresults$padj<0.05, "FDR<0.1", "Not Sig"))
p = ggplot(Fresults, aes(log2FoldChange, -log10(pvalue)))+geom_point(aes(col=sig))+
  scale_color_manual(values=c("red", "black"))
p
p+geom_text_repel(data=subset(Fresults, padj<.05 & abs(log2FoldChange)>1), aes(label=Gene))
new_volcano_list=as.data.frame(subset(Fresults, padj<.1))
üüüüüüüüü=as.data.frame(subset(Fresults, padj<.1 & abs(log2FoldChange)>1))
WriteXLS(new_volcano_list,"~/Desktop/new_volcano_list.xls",row.names = T)
WriteXLS(üüüüüüüüü,"~/Desktop/new.list.xls",row.names = T)

# Make a basic volcano plot
#with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
#with(subset(res, padj<.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
#with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
#with(subset(res, padj<.1 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
# Label points with the textxy function from the calibrate plot
#library(calibrate)
#with(subset(res, padj<.1 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=rownames(res), cex=.8))
### INTERESTING GENES BASED ON VOLCANO PLOT
counts=counts(Finni_Deseq_dds,normalized=T)
head(counts)
lol=Fresults[Fresults$Gene=="MDM2",]
lol=counts[rownames(counts)=="ENSG00000135679",]
lol=as.data.frame(lol)
lol$names=as.character(fin_info$name)
lol$trt=fin_info$treatment
ggplot(lol,aes(x=names,y=lol,color=trt))+geom_point(alpha=0.9,shape=16)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
