library(stringr)
library(refGenome)
library('sleuth')
library("dplyr")  

setwd("/Users/christopherpenfold/Desktop/Kei/2Cell/")

base_dir <- "/Users/christopherpenfold/Desktop/Kei/2Cell/"
saveext <- "/Users/christopherpenfold/Desktop/Thorsten/FINAL/Mouse_ZA/"

#Load in TFs
SIGNAL<-read.table("/Users/christopherpenfold/Desktop/LigandReceptor.csv",sep=",",header = F)
SIGNAL1 <- SIGNAL$V1[SIGNAL$V2=="Receptor" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Receptor"]
SIGNAL2 <- SIGNAL$V1[SIGNAL$V2=="Ligand" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand"]
SIGNAL3 <- SIGNAL$V1[SIGNAL$V2=="ECM" | SIGNAL$V2=="ECM/Receptor/Ligand" | SIGNAL$V2=="ECM/Ligand" | SIGNAL$V2=="ECM/Receptor"]
TF<-read.table("/Users/christopherpenfold/Desktop/Toshiaki/UberA/TF.txt",header = F)
TF <- TF$V1

#Load in DE analysis
D0_2 <- read.csv("/Users/christopherpenfold/Desktop/Thorsten/FINAL/Mouse_ZA/ml2C.csv")
D0_1 <- read.csv("/Users/christopherpenfold/Desktop/Thorsten/FINAL/Mouse_ZA/m2C.csv")
D0_3 <- read.csv("/Users/christopherpenfold/Desktop/Thorsten/FINAL/Mouse_ZA/m4C.csv")

tclist <-  D0_1$X[which(D0_1$p_val_adj<0.01 & abs(D0_1$avg_logFC)>log(1.5) )]
ltclist <- D0_2$X[which(D0_2$p_val_adj<0.01 & abs(D0_2$avg_logFC)>log(1.5) )]
fclist <-  D0_3$X[which(D0_3$p_val_adj<0.01 & abs(D0_3$avg_logFC)>log(1.5) )]

tclist <- as.data.frame(tclist)
ltclist <- as.data.frame(ltclist)
fclist  <- as.data.frame(fclist )

tclist$var1 <- 1
ltclist$var1 <- 1
fclist$var1 <- 1

colnames(tclist)[1] <- "genes"
colnames(ltclist)[1] <-  "genes"
colnames(fclist)[1] <-  "genes"

#Load in genomes
gtf <- ensemblGenome()
read.gtf(gtf, filename="../Mus_musculus.GRCm38.96.gtf")
genes = unique(gtf@ev$gtf[ ,c("transcript_id","gene_name")])

read.gtf(gtf, filename="../Mus_musculus_dba2j.DBA_2J_v1.94.gtf")
genes2 = unique(gtf@ev$gtf[ ,c("transcript_id","projection_parent_transcript","gene_name")])

read.gtf(gtf, filename="../Mus_musculus_c3hhej.C3H_HeJ_v1.94.gtf")
genes3 = unique(gtf@ev$gtf[ ,c("transcript_id","projection_parent_transcript","gene_name")])

genes3$projection_parent_transcript <- str_remove(genes3$projection_parent_transcript,"[.]15")
genes3$projection_parent_transcript <- str_remove(genes3$projection_parent_transcript,"[.]14")
genes3$projection_parent_transcript <- str_remove(genes3$projection_parent_transcript,"[.]13")
genes3$projection_parent_transcript <- str_remove(genes3$projection_parent_transcript,"[.]12")
genes3$projection_parent_transcript <- str_remove(genes3$projection_parent_transcript,"[.]11")
genes3$projection_parent_transcript <- str_remove(genes3$projection_parent_transcript,"[.]10")
genes3$projection_parent_transcript <- str_remove(genes3$projection_parent_transcript,"[.]1")
genes3$projection_parent_transcript <- str_remove(genes3$projection_parent_transcript,"[.]2")
genes3$projection_parent_transcript <- str_remove(genes3$projection_parent_transcript,"[.]3")
genes3$projection_parent_transcript <- str_remove(genes3$projection_parent_transcript,"[.]4")
genes3$projection_parent_transcript <- str_remove(genes3$projection_parent_transcript,"[.]5")
genes3$projection_parent_transcript <- str_remove(genes3$projection_parent_transcript,"[.]6")
genes3$projection_parent_transcript <- str_remove(genes3$projection_parent_transcript,"[.]7")
genes3$projection_parent_transcript <- str_remove(genes3$projection_parent_transcript,"[.]8")
genes3$projection_parent_transcript <- str_remove(genes3$projection_parent_transcript,"[.]9")

genes2$projection_parent_transcript <- str_remove(genes2$projection_parent_transcript,"[.]15")
genes2$projection_parent_transcript <- str_remove(genes2$projection_parent_transcript,"[.]14")
genes2$projection_parent_transcript <- str_remove(genes2$projection_parent_transcript,"[.]13")
genes2$projection_parent_transcript <- str_remove(genes2$projection_parent_transcript,"[.]12")
genes2$projection_parent_transcript <- str_remove(genes2$projection_parent_transcript,"[.]11")
genes2$projection_parent_transcript <- str_remove(genes2$projection_parent_transcript,"[.]10")
genes2$projection_parent_transcript <- str_remove(genes2$projection_parent_transcript,"[.]1")
genes2$projection_parent_transcript <- str_remove(genes2$projection_parent_transcript,"[.]2")
genes2$projection_parent_transcript <- str_remove(genes2$projection_parent_transcript,"[.]3")
genes2$projection_parent_transcript <- str_remove(genes2$projection_parent_transcript,"[.]4")
genes2$projection_parent_transcript <- str_remove(genes2$projection_parent_transcript,"[.]5")
genes2$projection_parent_transcript <- str_remove(genes2$projection_parent_transcript,"[.]6")
genes2$projection_parent_transcript <- str_remove(genes2$projection_parent_transcript,"[.]7")
genes2$projection_parent_transcript <- str_remove(genes2$projection_parent_transcript,"[.]8")
genes2$projection_parent_transcript <- str_remove(genes2$projection_parent_transcript,"[.]9")

names(genes3)[2] <- 'transcID'
names(genes2)[2] <- 'transcID'
names(genes)[1]  <- 'transcID'

merge1 <- merge(genes2, genes, by="transcID")
merge2 <- merge(genes3, genes, by="transcID")

#Get mapping transcripts to genes
Genes <- data.frame(c(genes$transcID,merge1$transcript_id,merge2$transcript_id),c(genes$gene_name,merge1$gene_name.y,merge2$gene_name.y))

names(Genes)[1] <- "target_id"
names(Genes)[2] <- names(genes)[2]

#First compare early
samples <- c("R_r1_R1","R_r2_R1","R_r3_R1","RAm_RAm1_R1","RAm_RAm2_R1","RAm_RAm3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("NoAmn","NoAmn","NoAmn","Amn","Amn","Amn"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_allgene_RvsRAm.csv")
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_all_RvsRAm.csv")

# WT test (to obtain the beta values):
model <- sleuth_wt(so, 'timepointNoAmn','full')
results_wt <- sleuth_results(model, 'timepointNoAmn', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_all_b_RvsRAm.csv")


#Now late 2C
samples <- c("Am_Am1_R1","Am_Am2_R1","Am_Am3_R1","NT_NT1_R1","NT_NT2_R1","NT_NT3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("NoAmn","NoAmn","NoAmn","Amn","Amn","Amn"), stringsAsFactors=FALSE)
so2 <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so2 <- sleuth_fit(so2, ~timepoint, 'full')
#And reduced model
so2 <- sleuth_fit(so2, ~1, 'reduced')
so2 <- sleuth_lrt(so2, 'reduced', 'full')
sleuth_table_gene_2 <- sleuth_results(so2, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene_2, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_allgene_AmvsNT.csv")
sleuth_table <- sleuth_results(so2, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_all_AmvsNT.csv")

#WT test (to obtain the beta values):
model <- sleuth_wt(so2, 'timepointNoAmn','full')
results_wt <- sleuth_results(model, 'timepointNoAmn', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_b_AmvsNT.csv")


#Now NT 21 vs 30
samples <- c("RAm_RAm1_R1","RAm_RAm2_R1","RAm_RAm3_R1","NT_NT1_R1","NT_NT2_R1","NT_NT3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("21h","21h","21h","30h","30h","30h"), stringsAsFactors=FALSE)
so2 <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so2 <- sleuth_fit(so2, ~timepoint, 'full')
#And reduced model
so2 <- sleuth_fit(so2, ~1, 'reduced')
so2 <- sleuth_lrt(so2, 'reduced', 'full')
sleuth_table_gene_2 <- sleuth_results(so2, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene_2, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_allgene_NTvsRAm.csv")
sleuth_table <- sleuth_results(so2, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/leuth_table_gene_all_NTvsRAm.csv")

# WT test (to obtain the beta values):
model <- sleuth_wt(so2, 'timepoint30h','full')
results_wt <- sleuth_results(model, 'timepoint30h', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/leuth_table_gene_b_NTvsRAm.csv")


#21 vs 30 (Am)
samples <- c("R_r1_R1","R_r2_R1","R_r3_R1","Am_Am1_R1","Am_Am2_R1","Am_Am3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("21h","21h","21h","30h","30h","30h"), stringsAsFactors=FALSE)
so2 <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so2 <- sleuth_fit(so2, ~timepoint, 'full')
#And reduced model
so2 <- sleuth_fit(so2, ~1, 'reduced')
so2 <- sleuth_lrt(so2, 'reduced', 'full')
sleuth_table_gene_2 <- sleuth_results(so2, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene_2, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_allgene_RvsAm.csv")
sleuth_table_gene <- sleuth_results(so2, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_all_RvsAm.csv")
# WT test (to obtain the beta values):
model <- sleuth_wt(so2, 'timepoint30h','full')
results_wt <- sleuth_results(model, 'timepoint30h', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_b_RvsAm.csv")

##### vs C2 C12
#First compare early
samples <- c("R1_v3","R2_v3","R3_v3","R4_v3","R5_v3","R6_v3","RAm_RAm1_R1","RAm_RAm2_R1","RAm_RAm3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("C2C12","C2C12","C2C12","C2C12","C2C12","C2C12","Amn","Amn","Amn"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_all_C2C12vsRAm.csv")
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)


# WT test (to obtain the beta values):
model <- sleuth_wt(so, 'timepointC2C12','full')
results_wt <- sleuth_results(model, 'timepointC2C12', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_b__C2C12vsRAm.csv")

#First compare early
samples <- c("R1_v3","R2_v3","R3_v3","R4_v3","R5_v3","R6_v3","NT_NT1_R1","NT_NT2_R1","NT_NT3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("C2C12","C2C12","C2C12","C2C12","C2C12","C2C12","Amn","Amn","Amn"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_all_C2C12vsNT.csv")
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
model <- sleuth_wt(so, 'timepointC2C12','full')
results_wt <- sleuth_results(model, 'timepointC2C12', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_b__C2C12vsNT.csv")

#First compare early
samples <- c("R1_v3","R2_v3","R3_v3","R4_v3","R5_v3","R6_v3","Am_Am1_R1","Am_Am2_R1","Am_Am3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("C2C12","C2C12","C2C12","C2C12","C2C12","C2C12","NoAmn","NoAmn","NoAmn"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_all_C2C12vsAm.csv")
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

model <- sleuth_wt(so, 'timepointNoAmn','full')
results_wt <- sleuth_results(model, 'timepointNoAmn', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_b__C2C12vsAm.csv")

#First compare early
samples <- c("R1_v3","R2_v3","R3_v3","R4_v3","R5_v3","R6_v3","R_r1_R1","R_r2_R1","R_r3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("C2C12","C2C12","C2C12","C2C12","C2C12","C2C12","NoAmn","NoAmn","NoAmn"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_all_C2C12vsR.csv")
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)


model <- sleuth_wt(so, 'timepointNoAmn','full')
results_wt <- sleuth_results(model, 'timepointNoAmn', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_b__C2C12vsR.csv")

#################################
#First run Xen vs MEF as baseline
#################################

samples <- c("Am_Am1_R1","Am_Am2_R1","Am_Am3_R1","NT_NT1_R1","NT_NT2_R1","NT_NT3_R1")
samples <- c("RAm_RAm1_R1","RAm_RAm2_R1","RAm_RAm3_R1","NT_NT1_R1","NT_NT2_R1","NT_NT3_R1")

samples <- c("../Xenopus/NT_R1/","../Xenopus/NT_R2/","../Xenopus/NT_R3/","Am_Am1_R1","Am_Am2_R1","Am_Am3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("NT","NT","NT","Am","Am","Am"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE, num_cores = 1)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_XenvsAm.csv")

model <- sleuth_wt(so, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_b__XenvsAm.csv")

samples <- c("../Xenopus/NT_R1/","../Xenopus/NT_R2/","../Xenopus/NT_R3/","NT_NT1_R1","NT_NT2_R1","NT_NT3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("NT","NT","NT","Am","Am","Am"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE, num_cores = 1)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_XenvsNT.csv")

model <- sleuth_wt(so, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_b__XenvsNT.csv")

samples <- c("../Xenopus/NT_R1/","../Xenopus/NT_R2/","../Xenopus/NT_R3/","RAm_RAm1_R1","RAm_RAm2_R1","RAm_RAm3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("NT","NT","NT","Am","Am","Am"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE, num_cores = 1)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_XenvsRAm.csv")

model <- sleuth_wt(so, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_b__XenvsRAm.csv")

samples <- c("../Xenopus/NT_R1/","../Xenopus/NT_R2/","../Xenopus/NT_R3/","R_r1_R1","R_r2_R1","R_r3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("NT","NT","NT","Am","Am","Am"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE, num_cores = 1)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_XenvsR.csv")

model <- sleuth_wt(so, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_b__XenvsR.csv")


##### vs C2 C12
#First compare early
samples <- c("R1_v3","R2_v3","R3_v3","R4_v3","R5_v3","R6_v3","RAm_RAm1_R1","RAm_RAm2_R1","RAm_RAm3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("C2C12","C2C12","C2C12","C2C12","C2C12","C2C12","Amn","Amn","Amn"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_all_C2C12vsRAm.csv")
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)


# WT test (to obtain the beta values):
model <- sleuth_wt(so, 'timepointC2C12','full')
results_wt <- sleuth_results(model, 'timepointC2C12', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_b__C2C12vsRAm.csv")


#First compare early
samples <- c("R1_v3","R2_v3","R3_v3","R4_v3","R5_v3","R6_v3","NT_NT1_R1","NT_NT2_R1","NT_NT3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("C2C12","C2C12","C2C12","C2C12","C2C12","C2C12","Amn","Amn","Amn"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_all_C2C12vsNT.csv")
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

model <- sleuth_wt(so, 'timepointC2C12','full')
results_wt <- sleuth_results(model, 'timepointC2C12', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_b__C2C12vsNT.csv")


#C2C12 other test
samples <- c("C2C12_R1_v2","C2C12_R2_v2","C2C12_R3_v2","C2C12_R4_v2","C2C12_R5_v2","R6_v3","Am_Am1_R1","Am_Am2_R1","Am_Am3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("C2C12","C2C12","C2C12","C2C12","C2C12","C2C12","NoAmn","NoAmn","NoAmn"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_all_C2C12vsAm.csv")
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

model <- sleuth_wt(so, 'timepointNoAmn','full')
results_wt <- sleuth_results(model, 'timepointNoAmn', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_b__C2C12vsAm.csv")


#First compare early
samples <- c("R1_v3","R2_v3","R3_v3","R4_v3","R5_v3","R6_v3","R_r1_R1","R_r2_R1","R_r3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("C2C12","C2C12","C2C12","C2C12","C2C12","C2C12","NoAmn","NoAmn","NoAmn"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_all_C2C12vsR.csv")
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)


model <- sleuth_wt(so, 'timepointNoAmn','full')
results_wt <- sleuth_results(model, 'timepointNoAmn', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/2Cell/sleuth_table_gene_b__C2C12vsR.csv")
