library("DESeq2")
library(dplyr)
library(ggplot2)
library("ggbeeswarm")
library("RColorBrewer")
library("gplots")

directory <- "."
setwd(directory)

sampleFiles <- c("basic_1_S9_2/Res.counts","basic_2_S10_2/Res.counts","basic_3_S11_2/Res.counts","basic_aAm1_S12_2/Res.counts","basic_aAm2_S13_2/Res.counts","basic_aAm3_S14_2/Res.counts")
sampleFiles2 <- c("basic_1_S9_2C/Res.counts","basic_2_S10_2C/Res.counts","basic_3_S11_2C/Res.counts","basic_aAm1_S12_2C/Res.counts","basic_aAm2_S13_2C/Res.counts","basic_aAm3_S14_2C/Res.counts")
sampleFiles3 <- c("basic_1_S9_2/Res4.counts","basic_2_S10_2/Res4.counts","basic_3_S11_2/Res4.counts","basic_aAm1_S12_2/Res4.counts","basic_aAm2_S13_2/Res4.counts","basic_aAm3_S14_2/Res4.counts")

outputPrefix <- "_analysis"

sampleNames <- c("NT1","NT2","NT3","Am1","Am2","Am3")
sampleCondition <- c("NT","NT","NT","Am","Am","Am")
sampleReplicate <- c("one","two","three","one","two","three")
sampleReag <- c("one","two","three","one","two","three")

sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition, rep = sampleReplicate, reag = sampleReag)
sampleTable2 <- data.frame(sampleName = sampleNames, fileName = sampleFiles2, condition = sampleCondition, rep = sampleReplicate, reag = sampleReag)
sampleTable3 <- data.frame(sampleName = sampleNames, fileName = sampleFiles3, condition = sampleCondition, rep = sampleReplicate, reag = sampleReag)

ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory="/mnt/scratch/gurdon/cap76/Kei/Oryx3/", design=~condition)
ddsHTSeq2<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable2, directory="/mnt/scratch/gurdon/cap76/Kei/Oryx3/", design=~condition)
ddsHTSeq3<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable3, directory="/mnt/scratch/gurdon/cap76/Kei/Oryx3/", design=~condition)

dds <- DESeq(ddsHTSeq)
res <- results(dds)

dds2 <- DESeq(ddsHTSeq2)
res2 <- results(dds2)

dds3 <- DESeq(ddsHTSeq3)
res3 <- results(dds3)

#Extract out pairwise comparisons
res1<-results(dds, contrast=c("condition","NT","Am"))
res2<-results(dds2, contrast=c("condition","NT","Am"))
res3<-results(dds3, contrast=c("condition","NT","Am"))

resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata1)[1] <- 'gene'
write.csv(resdata1, file = "NTvsAm_mouse.csv")
write.table(as.data.frame(counts(dds),normalized=T), file = paste0("Normalized_counts.txt"), sep = '\t')

resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds2,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata2)[1] <- 'gene'
write.csv(resdata2, file = "NTvsAm_oryx.csv")
write.table(as.data.frame(counts(dds2),normalized=T), file = paste0("Normalized_counts2.txt"), sep = '\t')

resdata3 <- merge(as.data.frame(res3), as.data.frame(counts(dds3,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata3)[1] <- 'gene'
write.csv(resdata3, file = "NTvsAm_mouseoryx.csv")
write.table(as.data.frame(counts(dds3),normalized=T), file = paste0("Normalized_counts3.txt"), sep = '\t')

#Transform raw counts into normalized values. DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
rld <- rlogTransformation(dds, blind=F)
vsd <- varianceStabilizingTransformation(dds, blind=F)
write.table(as.data.frame(assay(rld),file = paste0("Rlog-transformed-counts.txt"), sep = '\t'))
write.table(as.data.frame(assay(vsd),file = paste0("vst-transformed-counts.txt"), sep = '\t'))

rld2 <- rlogTransformation(dds2, blind=F)
vsd2 <- varianceStabilizingTransformation(dds2, blind=F)
write.table(as.data.frame(assay(rld2),file = paste0("Rlog-transformed-counts2.txt"), sep = '\t'))
write.table(as.data.frame(assay(vsd2),file = paste0("vst-transformed-counts2.txt"), sep = '\t'))

rld3 <- rlogTransformation(dds3, blind=F)
vsd3 <- varianceStabilizingTransformation(dds3, blind=F)
write.table(as.data.frame(assay(rld3),file = paste0("Rlog-transformed-counts3.txt"), sep = '\t'))
write.table(as.data.frame(assay(vsd3),file = paste0("vst-transformed-counts3.txt"), sep = '\t'))

