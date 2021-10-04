library("DESeq2")
library(dplyr)
library(ggplot2)
library(stringr)
#library("ggbeeswarm")
#library("RColorBrewer")
#library("gplots")
#library(GenomicFeatures)

#txdb <- makeTxDbFromGFF("/mnt/scratch/surani/cap76/Genomes/UCSC/hg19/Annotation/Genes/genes.gtf",format="gtf")
#exons.list.per.gene <- exonsBy(txdb,by="gene")
#exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
#geneLength <- as.data.frame(unlist(exonic.gene.sizes))

directory <- "."
setwd(directory)

library("refGenome")

gtf1 <- ensemblGenome()
read.gtf(gtf1, filename="../Genomes/Mus_musculus_c3hhej.C3H_HeJ_v1.94.gtf")
genes_1 =  unique(gtf1@ev$gtf[ ,c("gene_id","gene_name")])

gtf2 <- ensemblGenome()
read.gtf(gtf2, filename="../Genomes/Mus_musculus_c57bl6nj.C57BL_6NJ_v1.94.gtf")
genes_2 =  unique(gtf2@ev$gtf[ ,c("gene_id","gene_name")])

gtf3 <- ensemblGenome()
read.gtf(gtf3, filename="../Genomes/Mus_musculus_dba2j.DBA_2J_v1.94.gtf")
genes_3 =  unique(gtf3@ev$gtf[ ,c("gene_id","gene_name")])

library("GenomicFeatures")

refgenes1 <- makeTxDbFromGFF("../Genomes/Mus_musculus_c3hhej.C3H_HeJ_v1.94.gtf",format="gtf")
refgenes2 <- makeTxDbFromGFF("../Genomes/Mus_musculus_c57bl6nj.C57BL_6NJ_v1.94.gtf",format="gtf")
refgenes3 <- makeTxDbFromGFF("../Genomes/Mus_musculus_dba2j.DBA_2J_v1.94.gtf",format="gtf")
refgenes4 <- makeTxDbFromGFF("/mnt/scratch/surani/cap76/Genomes/UCSC/mm10/Annotation/Genes/genes.gtf",format="gtf")

transcripts1 <- genes(refgenes1, columns=c("tx_id", "tx_name"))
tss1 <- as.data.frame(resize(transcripts1, width=1, fix='start'))
transcripts1 <- as.data.frame(transcripts1)
rownames(genes_1) <- genes_1[,1]
write.table(merge(tss1, genes_1, by ="row.names", all.x = TRUE), file = "TSS_C3H.bed",quote=FALSE,sep="\t",col.names=FALSE)
write.table(merge(transcripts1, genes_1, by ="row.names", all.x = TRUE), file = "gene_body_C3H.bed",quote=FALSE,sep="\t",col.names=FALSE)

transcripts2 <- genes(refgenes2, columns=c("tx_id", "tx_name"))
tss2 <- as.data.frame(resize(transcripts2, width=1, fix='start'))
transcripts2 <- as.data.frame(transcripts2)
rownames(genes_2) <- genes_2[,1]
write.table(merge(tss2, genes_2, by ="row.names", all.x = TRUE), file = "TSS_BL6.bed",quote=FALSE,sep="\t",col.names=FALSE)
write.table(merge(transcripts2, genes_2, by ="row.names", all.x = TRUE), file = "gene_body_BL6.bed",quote=FALSE,sep="\t",col.names=FALSE)

transcripts3 <- genes(refgenes3, columns=c("tx_id", "tx_name"))
tss3 <- as.data.frame(resize(transcripts3, width=1, fix='start'))
transcripts3 <- as.data.frame(transcripts3)
rownames(genes_3) <- genes_3[,1]
write.table(merge(tss3, genes_3, by ="row.names", all.x = TRUE), file = "TSS_DBA.bed",quote=FALSE,sep="\t",col.names=FALSE)
write.table(merge(transcripts3, genes_3, by ="row.names", all.x = TRUE), file = "gene_body_DBA.bed",quote=FALSE,sep="\t",col.names=FALSE)

transcripts4 <- genes(refgenes4, columns=c("tx_id", "tx_name"))
tss4 <- as.data.frame(resize(transcripts4, width=1, fix='start'))
transcripts4 <- as.data.frame(transcripts4)
write.table(as.data.frame(tss4), file = "TSS_mm10.bed",quote=FALSE,sep="\t",col.names=FALSE)
write.table(as.data.frame(transcripts4), file = "gene_body_mm10.bed",quote=FALSE,sep="\t",col.names=FALSE)


