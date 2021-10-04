library(stringr)
library(refGenome)
library('sleuth')


setwd("/Users/christopherpenfold/Desktop/Kei/2Cell/")

base_dir <- "/Users/christopherpenfold/Desktop/Kei/2Cell/"
saveext <- "/Users/christopherpenfold/Desktop/Thorsten/FINAL/Mouse_ZA/"


#Now run DEseq2 for various combinations
gtf <- ensemblGenome()
read.gtf(gtf, filename="../Mus_musculus.GRCm38.96.gtf")
genes = unique(gtf@ev$gtf[ ,c("transcript_id","gene_name")])

read.gtf(gtf, filename="../Mus_musculus_dba2j.DBA_2J_v1.94.gtf")
genes2 = unique(gtf@ev$gtf[ ,c("transcript_id","projection_parent_transcript","gene_name")])

read.gtf(gtf, filename="../Mus_musculus_c3hhej.C3H_HeJ_v1.94.gtf")
genes3 = unique(gtf@ev$gtf[ ,c("transcript_id","projection_parent_transcript","gene_name")])

genes3$projection_parent_transcript2 <- genes3$projection_parent_transcript
genes3$projection_parent_transcript2 <- str_remove(genes3$projection_parent_transcript2,"[.]10")
genes3$projection_parent_transcript2 <- str_remove(genes3$projection_parent_transcript2,"[.]1")
genes3$projection_parent_transcript2 <- str_remove(genes3$projection_parent_transcript2,"[.]2")
genes3$projection_parent_transcript2 <- str_remove(genes3$projection_parent_transcript2,"[.]3")
genes3$projection_parent_transcript2 <- str_remove(genes3$projection_parent_transcript2,"[.]4")
genes3$projection_parent_transcript2 <- str_remove(genes3$projection_parent_transcript2,"[.]5")
genes3$projection_parent_transcript2 <- str_remove(genes3$projection_parent_transcript2,"[.]6")
genes3$projection_parent_transcript2 <- str_remove(genes3$projection_parent_transcript2,"[.]7")
genes3$projection_parent_transcript2 <- str_remove(genes3$projection_parent_transcript2,"[.]8")
genes3$projection_parent_transcript2 <- str_remove(genes3$projection_parent_transcript2,"[.]9")

genes2$projection_parent_transcript2 <- genes2$projection_parent_transcript
genes2$projection_parent_transcript2 <- str_remove(genes2$projection_parent_transcript2,"[.]10")
genes2$projection_parent_transcript2 <- str_remove(genes2$projection_parent_transcript2,"[.]1")
genes2$projection_parent_transcript2 <- str_remove(genes2$projection_parent_transcript2,"[.]2")
genes2$projection_parent_transcript2 <- str_remove(genes2$projection_parent_transcript2,"[.]3")
genes2$projection_parent_transcript2 <- str_remove(genes2$projection_parent_transcript2,"[.]4")
genes2$projection_parent_transcript2 <- str_remove(genes2$projection_parent_transcript2,"[.]5")
genes2$projection_parent_transcript2 <- str_remove(genes2$projection_parent_transcript2,"[.]6")
genes2$projection_parent_transcript2 <- str_remove(genes2$projection_parent_transcript2,"[.]7")
genes2$projection_parent_transcript2 <- str_remove(genes2$projection_parent_transcript2,"[.]8")
genes2$projection_parent_transcript2 <- str_remove(genes2$projection_parent_transcript2,"[.]9")


names(genes3)[2] <- 'transcID'
names(genes2)[2] <- 'transcID'
names(genes)[1] <- 'transcID'

merge1 <- merge(genes2, genes, by="transcID")
merge2 <- merge(genes3, genes, by="transcID")


Genes <- data.frame(c(genes$transcID,merge1$transcript_id,merge2$transcript_id),c(genes$gene_name,merge1$gene_name.y,merge2$gene_name.y))


#Genes <- data.frame(c(genes$transcript_id,genes2$transcript_id,genes3$transcript_id),c(genes$gene_name,genes2$gene_name,genes3$gene_name))
names(Genes)[1] <- "target_id" #names(genes)[1]
names(Genes)[2] <- names(genes)[2]

#Genes$gene_id <- cbind(genes$gene_id,genes2$gene_name.x,genes3$gene_name.x)
#Genes$gene_name <- list(genes$gene_name,genes2$gene_name.y,genes3$gene_name.y)
#merge1 <- merge(genes2, genes, by="gene_id")

setwd("/Users/christopherpenfold/Desktop/Kei/Rerun/")

#ens <- ensemblGenome()
#read.gtf(ens, "Mus_musculus_c57bl6nj.C57BL_6NJ_v1.95.gtf")
#genes = ens@ev$genes[ ,c("gene_id","gene_name")]
#D1 <- read.csv(file="/Users/christopherpenfold/Desktop/Kei/transcripttogene.txt", header=0, sep="\t")
#names(D1)<-c('target_id','ens_gene','ext_gene')
#ens@ev$genes[ ,c(D1[:,1],D1[:,2])]

#names(Genes)[1] <- "target_id"
#names(Genes)[2] <- "gene_id"
#
base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"
samples <- c("DR_DR1","DR_DR2","DR_DR3","D_D1","D_D4","D_D5")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("DNA","DNA","DNA","noDNA","noDNA","noDNA"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_D_DR.csv")


model <- sleuth_wt(so, 'timepointnoDNA','full')
results_wt <- sleuth_results(model, 'timepointnoDNA', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_D_DR.csv")



samples <- c("1-8_S14_v2","1-16_S16_v2","1-18_S18_v2","1-15_S15_v2","1-17_S13_v2","1-17_S17_v2")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("A-ctrl","A-ctrl","A-ctrl","PBS","PBS","PBS"), stringsAsFactors=FALSE)
so2 <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so2 <- sleuth_fit(so2, ~timepoint, 'full')
#And reduced model
so2 <- sleuth_fit(so2, ~1, 'reduced')
so2 <- sleuth_lrt(so2, 'reduced', 'full')
sleuth_table_gene_2 <- sleuth_results(so2, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene_2, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_1-8_1-15.csv")
#models(so)

model <- sleuth_wt(so2, 'timepointPBS','full')
results_wt <- sleuth_results(model, 'timepointPBS', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_1-8_1-15.csv")




base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"
samples <- c("DR_DR1","DR_DR2","DR_DR3","D_D1","D_D4","D_D5")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("DNA","DNA","DNA","noDNA","noDNA","noDNA"), stringsAsFactors=FALSE)
so3 <- sleuth_prep(s2c, read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE, aggregation_column = 'gene_name')
#Fit full model
so3 <- sleuth_fit(so3, ~timepoint, 'full')
#And reduced model
so3 <- sleuth_fit(so3, ~1, 'reduced')
so3 <- sleuth_lrt(so3, 'reduced', 'full')
sleuth_table_gene_3 <- sleuth_results(so3, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene_3, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all3.csv")

samples <- c("1-8_S14_v2","1-16_S16_v2","1-18_S18_v2","1-15_S15_v2","1-17_S13_v2","1-17_S17_v2")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("A-ctrl","A-ctrl","A-ctrl","PBS","PBS","PBS"), stringsAsFactors=FALSE)
so4 <- sleuth_prep(s2c, read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE, aggregation_column = 'gene_name')
#Fit full model
so4 <- sleuth_fit(so4, ~timepoint, 'full')
#And reduced model
so4 <- sleuth_fit(so4, ~1, 'reduced')
so4 <- sleuth_lrt(so4, 'reduced', 'full')
sleuth_table_gene_4 <- sleuth_results(so4, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene_4, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all4.csv")
#models(so)

#############################
#Versus 2C
#############################
#samples <- c("R_r1_R1","R_r2_R1","R_r3_R1","RAm_RAm1_R1","RAm_RAm2_R1","RAm_RAm3_R1")
#kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
#s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("NoAmn","NoAmn","NoAmn","Amn","Amn","Amn"), stringsAsFactors=FALSE)
#samples <- c("Am_Am1_R1","Am_Am2_R1","Am_Am3_R1","NT_NT1_R1","NT_NT2_R1","NT_NT3_R1")
#kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
#s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("NoAmn","NoAmn","NoAmn","Amn","Amn","Amn"), stringsAsFactors=FALSE)




base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"
samples <- c("DR_DR1","DR_DR2","DR_DR3","../2Cell/R_r1_R1","../2Cell/R_r2_R1","../2Cell/R_r3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("DNA","DNA","DNA","noDNA","noDNA","noDNA"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_DR_R.csv")

model <- sleuth_wt(so, 'timepointnoDNA','full')
results_wt <- sleuth_results(model, 'timepointnoDNA', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_DR_R.csv")


base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"
samples <- c("D_D1","D_D4","D_D5","../2Cell/R_r1_R1","../2Cell/R_r2_R1","../2Cell/R_r3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("DNA","DNA","DNA","noDNA","noDNA","noDNA"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_D_R.csv")

model <- sleuth_wt(so, 'timepointnoDNA','full')
results_wt <- sleuth_results(model, 'timepointnoDNA', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_D_R.csv")







base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"
samples <- c("DR_DR1","DR_DR2","DR_DR3","../2Cell/RAm_RAm1_R1","../2Cell/RAm_RAm2_R1","../2Cell/RAm_RAm3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("DNA","DNA","DNA","noDNA","noDNA","noDNA"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_DR_RAm.csv")

model <- sleuth_wt(so, 'timepointnoDNA','full')
results_wt <- sleuth_results(model, 'timepointnoDNA', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_DR_RAm.csv")


base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"
samples <- c("D_D1","D_D4","D_D5","../2Cell/RAm_RAm1_R1","../2Cell/RAm_RAm2_R1","../2Cell/RAm_RAm3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("DNA","DNA","DNA","noDNA","noDNA","noDNA"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_D_RAm.csv")

model <- sleuth_wt(so, 'timepointnoDNA','full')
results_wt <- sleuth_results(model, 'timepointnoDNA', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_D_RAm.csv")





base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"
samples <- c("DR_DR1","DR_DR2","DR_DR3","../2Cell/Am_Am1_R1","../2Cell/Am_Am2_R1","../2Cell/Am_Am3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("DNA","DNA","DNA","noDNA","noDNA","noDNA"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_DR_Am.csv")

model <- sleuth_wt(so, 'timepointnoDNA','full')
results_wt <- sleuth_results(model, 'timepointnoDNA', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_DR_Am.csv")

base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"
samples <- c("D_D1","D_D4","D_D5","../2Cell/Am_Am1_R1","../2Cell/Am_Am2_R1","../2Cell/Am_Am3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("DNA","DNA","DNA","noDNA","noDNA","noDNA"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_D_Am.csv")

model <- sleuth_wt(so, 'timepointnoDNA','full')
results_wt <- sleuth_results(model, 'timepointnoDNA', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_D_Am.csv")








base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"
samples <- c("DR_DR1","DR_DR2","DR_DR3","../2Cell/NT_NT1_R1","../2Cell/NT_NT2_R1","../2Cell/NT_NT3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("DNA","DNA","DNA","noDNA","noDNA","noDNA"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_DR_NT.csv")

model <- sleuth_wt(so, 'timepointnoDNA','full')
results_wt <- sleuth_results(model, 'timepointnoDNA', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_DR_NT.csv")


base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"
samples <- c("D_D1","D_D4","D_D5","../2Cell/NT_NT1_R1","../2Cell/NT_NT2_R1","../2Cell/NT_NT3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("DNA","DNA","DNA","noDNA","noDNA","noDNA"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_D_NT.csv")

model <- sleuth_wt(so, 'timepointnoDNA','full')
results_wt <- sleuth_results(model, 'timepointnoDNA', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_D_NT.csv")










base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"
samples <- c("1-15_S15_v2","1-17_S13_v2","1-17_S17_v2","../2Cell/R_r1_R1","../2Cell/R_r2_R1","../2Cell/R_r3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("DNA","DNA","DNA","noDNA","noDNA","noDNA"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_1-15_R.csv")

model <- sleuth_wt(so, 'timepointnoDNA','full')
results_wt <- sleuth_results(model, 'timepointnoDNA', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_1-15_R.csv")


base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"
samples <- c("1-15_S15_v2","1-17_S13_v2","1-17_S17_v2","../2Cell/RAm_RAm1_R1","../2Cell/RAm_RAm2_R1","../2Cell/RAm_RAm3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("DNA","DNA","DNA","noDNA","noDNA","noDNA"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_1-15_RAm.csv")

model <- sleuth_wt(so, 'timepointnoDNA','full')
results_wt <- sleuth_results(model, 'timepointnoDNA', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_1-15_RAm.csv")











base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"
samples <- c("1-15_S15_v2","1-17_S13_v2","1-17_S17_v2","../2Cell/NT_NT1_R1","../2Cell/NT_NT2_R1","../2Cell/NT_NT3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("DNA","DNA","DNA","noDNA","noDNA","noDNA"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_1-15_NT.csv")

model <- sleuth_wt(so, 'timepointnoDNA','full')
results_wt <- sleuth_results(model, 'timepointnoDNA', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_1-15_NT.csv")


base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"
samples <- c("1-15_S15_v2","1-17_S13_v2","1-17_S17_v2","../2Cell/Am_Am1_R1","../2Cell/Am_Am2_R1","../2Cell/Am_Am3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("DNA","DNA","DNA","noDNA","noDNA","noDNA"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_1-15_Am.csv")

model <- sleuth_wt(so, 'timepointnoDNA','full')
results_wt <- sleuth_results(model, 'timepointnoDNA', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_1-15_Am.csv")


#models(so)

model <- sleuth_wt(so2, 'timepointPBS','full')
results_wt <- sleuth_results(model, 'timepointPBS', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_1-8_1-15.csv")



plot_bootstrap(so, "ENST00000263734", units = "est_counts", color_by = "condition")

#so <- sleuth_wt(so, which_beta="timepointt1") 

#

#so <- sleuth_prep(s2c, ~timepoint,read_bootstrap_tpm = TRUE, target_mapping = D1, extra_bootstrap_summary = TRUE, aggregation_column = 'ens_gene', gene_mode = TRUE)
#so <- sleuth_fit(so, ~timepoint, 'reduced')
#so <- sleuth_fit(so, ~1, 'reduced')
#so <- sleuth_lrt(so, 'reduced', 'full')

#write.table(tst, file = "/Users/christopherpenfold/Desktop/Kei/Kei/Rerun/GeneLevel.csv")
#tst<-sleuth_to_matrix(so,"obs_norm","tpm")
#sleuth_live(so)
#so <- sleuth_prep(s2c, ~condition, target_mapping = t2g,
#                  aggregation_column = 'ens_gene')
#so <- sleuth_prep(so, target_mapping = D1, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE)
#s2c <- data.frame(path=kal_dirs, sample=samples, target_mapping = D1), timepoint = c("ctrl","ctrl","ctrl","t1","t1","t1"), stringsAsFactors=FALSE)

#so <- sleuth_fit(so, ~timepoint, 'full')
#so <- sleuth_fit(so, ~timepoint, 'reduced')
#so <- sleuth_lrt(so, 'reduced', 'full')

#sleuth_table_tx <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = FALSE)
#sleuth_table_tx <- dplyr::filter(sleuth_table_tx, qval <= 0.05)
#head(sleuth_table_tx, 20)

library(ggplot2)
sleuth_table_tx <- dplyr::filter(sleuth_table_gene, pval <= 0.01)
#list2 <- head(sleuth_table_gene, 100)
for (i in 1:dim(sleuth_table_tx)[1]){
  plot_bootstrap(so, sleuth_table_tx$target_id[i], units = "est_counts", color_by = "timepoint")
  ggsave(filename=paste("/Users/christopherpenfold/Desktop/Kei/Rerun/","Plot_",i,"_DNA_vs_noDNA.pdf",sep=""))
}


library(ggplot2)
sleuth_table_tx <- dplyr::filter(sleuth_table_gene_2, pval <= 0.001)
#list2 <- head(sleuth_table_gene_2, 100)
for (i in 1:dim(sleuth_table_tx)[1]){
  plot_bootstrap(so2, sleuth_table_tx$target_id[i], units = "est_counts", color_by = "timepoint")
  ggsave(filename=paste("/Users/christopherpenfold/Desktop/Kei/Rerun/","Plot_",i,"_Amanitin_vs_PBS.pdf",sep=""))
}


