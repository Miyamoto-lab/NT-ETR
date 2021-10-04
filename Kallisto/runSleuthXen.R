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

names(Genes)[1] <- "target_id"
names(Genes)[2] <- names(genes)[2]

setwd("/Users/christopherpenfold/Desktop/Kei/Rerun/")
base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"

#################################
#First run Xen vs MEF as baseline
#################################
samples <- c("../Xenopus/NT_R1/","../Xenopus/NT_R2/","../Xenopus/NT_R3/","../Xenopus/MEF1/","../Xenopus/MEF2/","../Xenopus/MEF3/")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("NT","NT","NT","MEF","MEF","MEF"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE, num_cores = 1)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced modeldir
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_XenvsMEF.csv")

model <- sleuth_wt(so, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_XenvsMEF.csv")

samples <- c("../Xenopus/NT_R1/","../Xenopus/NT_R2/","../Xenopus/NT_R3/","../Xenopus/NT_T1/","../Xenopus/NT_T2/","../Xenopus/NT_T3/")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("NT","NT","NT","Ctrl","Ctrl","Ctrl"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE, num_cores = 1)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_XenvsCtrl.csv")

model <- sleuth_wt(so, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_XenvsCtrl.csv")


##############################################
#Now four cell
##############################################


base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"
samples <- c("../Xenopus/NT_R1/","../Xenopus/NT_R2/","../Xenopus/NT_R3/","D_D1","D_D4","D_D5")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("NT","NT","NT","noDNA","noDNA","noDNA"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_XNT_D.csv")

model <- sleuth_wt(so, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_XNT_D.csv")



###didnt run
samples <- c("../Xenopus/NT_R1/","../Xenopus/NT_R2/","../Xenopus/NT_R3/","1-15_S15_v2","1-17_S13_v2","1-17_S17_v2")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("NT","NT","NT","PBS","PBS","PBS"), stringsAsFactors=FALSE)
so2 <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so2 <- sleuth_fit(so2, ~timepoint, 'full')
#And reduced model
so2 <- sleuth_fit(so2, ~1, 'reduced')
so2 <- sleuth_lrt(so2, 'reduced', 'full')
sleuth_table_gene_2 <- sleuth_results(so2, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene_2, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_XNT_1-15.csv")
#models(so)

model <- sleuth_wt(so2, 'timepointPBS','full')
results_wt <- sleuth_results(model, 'timepointPBS', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_XNT_1-15.csv")



base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"
samples <- c("DR_DR1","DR_DR2","DR_DR3","../Xenopus/NT_R1/","../Xenopus/NT_R2/","../Xenopus/NT_R3/")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("DNA","DNA","DNA","NT","NT","NT"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_XNT_DR.csv")

model <- sleuth_wt(so, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_XNT_DR.csv")



#Didn't run???
#kal_files <- file.path("1-8_S14_v2", "abundance.h5")
#kal_files <- file.path("../Xenopus/NT_R2_v3", "abundance.h5")

#n_targets <- sapply(kal_files, function(file) {
#  ids <- rhdf5::h5read(file, "aux/ids")
#  n_targets <- length(ids)
#  n_targets
#})
#versions <- sapply(kal_files, function(file) {
#  file_name <- basename(file)
#  version <- rhdf5::h5read(file, "aux/kallisto_version")
#  version
#})



#"DR_DR1","DR_DR2","DR_DR3"

samples <- c("../2Cell/R1_v3","../2Cell/R2_v3","../2Cell/R3_v3","../2Cell/R4_v3","../2Cell/R5_v3","../2Cell/R6_v3","DR_DR1","DR_DR2","DR_DR3")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("C2C12","C2C12","C2C12","C2C12","C2C12","C2C12","NT","NT","NT"), stringsAsFactors=FALSE)
so2 <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so2 <- sleuth_fit(so2, ~timepoint, 'full')
#And reduced model
so2 <- sleuth_fit(so2, ~1, 'reduced')
so2 <- sleuth_lrt(so2, 'reduced', 'full')
sleuth_table_gene_2 <- sleuth_results(so2, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene_2, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_C2C12_vs_DR.csv")
#models(so)
model <- sleuth_wt(so2, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_C2C12_vs_DR.csv")





samples <- c("../2Cell/R1_v3","../2Cell/R2_v3","../2Cell/R3_v3","../2Cell/R4_v3","../2Cell/R5_v3","../2Cell/R6_v3","1-15_S15_v2","1-17_S13_v2","1-17_S17_v2")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("C2C12","C2C12","C2C12","C2C12","C2C12","C2C12","NT","NT","NT"), stringsAsFactors=FALSE)
so2 <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so2 <- sleuth_fit(so2, ~timepoint, 'full')
#And reduced model
so2 <- sleuth_fit(so2, ~1, 'reduced')
so2 <- sleuth_lrt(so2, 'reduced', 'full')
sleuth_table_gene_2 <- sleuth_results(so2, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene_2, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_C2C12_vs_1-15.csv")
#models(so)
model <- sleuth_wt(so2, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_C2C12_vs_1-15.csv")






samples <- c("../2Cell/R1_v3","../2Cell/R2_v3","../2Cell/R3_v3","../2Cell/R4_v3","../2Cell/R5_v3","../2Cell/R6_v3","DR_DR1","DR_DR2","DR_DR3")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("C2C12","C2C12","C2C12","C2C12","C2C12","C2C12","NT","NT","NT"), stringsAsFactors=FALSE)
so2 <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so2 <- sleuth_fit(so2, ~timepoint, 'full')
#And reduced model
so2 <- sleuth_fit(so2, ~1, 'reduced')
so2 <- sleuth_lrt(so2, 'reduced', 'full')
sleuth_table_gene_2 <- sleuth_results(so2, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene_2, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_C2C12_vs_DR.csv")
#models(so)
model <- sleuth_wt(so2, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_C2C12_vs_DR.csv")







samples <- c("../2Cell/R1_v3","../2Cell/R2_v3","../2Cell/R3_v3","../2Cell/R4_v3","../2Cell/R5_v3","../2Cell/R6_v3","1-8_S14_v2","1-16_S16_v2","1-18_S18_v2")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("C2C12","C2C12","C2C12","C2C12","C2C12","C2C12","NT","NT","NT"), stringsAsFactors=FALSE)
so2 <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so2 <- sleuth_fit(so2, ~timepoint, 'full')
#And reduced model
so2 <- sleuth_fit(so2, ~1, 'reduced')
so2 <- sleuth_lrt(so2, 'reduced', 'full')
sleuth_table_gene_2 <- sleuth_results(so2, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene_2, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_C2C12_vs_1-8.csv")
#models(so)
model <- sleuth_wt(so2, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_C2C12_vs_1-8.csv")








samples <- c("1-8_S14_v2","1-16_S16_v2","1-18_S18_v2","../Xenopus/NT_R1/","../Xenopus/NT_R2/","../Xenopus/NT_R3/")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("A-ctrl","A-ctrl","A-ctrl","NT","NT","NT"), stringsAsFactors=FALSE)
so2 <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so2 <- sleuth_fit(so2, ~timepoint, 'full')
#And reduced model
so2 <- sleuth_fit(so2, ~1, 'reduced')
so2 <- sleuth_lrt(so2, 'reduced', 'full')
sleuth_table_gene_2 <- sleuth_results(so2, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene_2, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_XNT_1-8.csv")
#models(so)
model <- sleuth_wt(so2, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_XNT_1-8.csv")



#############################################
#Now the same for MEFs
#############################################

#First 2C

samples <- c("../Xenopus/NT_R1/","../Xenopus/NT_R2/","../Xenopus/NT_R3/","../2Cell/R_r1_R1","../2Cell/R_r2_R1","../2Cell/R_r3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("NT","NT","NT","2CnoAm","2CnoAm","2CnoAm"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_XNT_2CR.csv")

model <- sleuth_wt(so, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_XNT_2CR.csv")




samples <- c("../Xenopus/NT_R1/","../Xenopus/NT_R2/","../Xenopus/NT_R3/","../2Cell/RAm_RAm1_R1","../2Cell/RAm_RAm2_R1","../2Cell/RAm_RAm3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("NT","NT","NT","2CAm","2CAm","2CAm"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_XNT_2CRAm.csv")

model <- sleuth_wt(so, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_XNT_2CRAm.csv")


samples <- c("../Xenopus/NT_R1/","../Xenopus/NT_R2/","../Xenopus/NT_R3/","../2Cell/Am_Am1_R1","../2Cell/Am_Am2_R1","../2Cell/Am_Am3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("NT","NT","NT","2CnoAm","2CnoAm","2CnoAm"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_XNT_2CAm.csv")

model <- sleuth_wt(so, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_XNT_2CAm.csv")


samples <- c("../Xenopus/NT_R1/","../Xenopus/NT_R2/","../Xenopus/NT_R3/","../2Cell/NT_NT1_R1","../2Cell/NT_NT2_R1","../2Cell/NT_NT3_R1")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("NT","NT","NT","2CAm","2CAm","2CAm"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_XNT_2CNT.csv")

model <- sleuth_wt(so, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_XNT_2CNT.csv")

#Now 4C


base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"
samples <- c("../Xenopus/NT_R1/","../Xenopus/NT_R2/","../Xenopus/NT_R3/","D_D1","D_D4","D_D5")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("NT","NT","NT","noDNA","noDNA","noDNA"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_XNT_D.csv")

model <- sleuth_wt(so, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_XNT_D.csv")



###didnt run
samples <- c("../Xenopus/NT_R1/","../Xenopus/NT_R2/","../Xenopus/NT_R3/","1-15_S15_v2","1-17_S13_v2","1-17_S17_v2")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("NT","NT","NT","PBS","PBS","PBS"), stringsAsFactors=FALSE)
so2 <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so2 <- sleuth_fit(so2, ~timepoint, 'full')
#And reduced model
so2 <- sleuth_fit(so2, ~1, 'reduced')
so2 <- sleuth_lrt(so2, 'reduced', 'full')
sleuth_table_gene_2 <- sleuth_results(so2, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene_2, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_XNT_1-15.csv")
#models(so)

model <- sleuth_wt(so, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_XNT_1-15.csv")



base_dir <- "/Users/christopherpenfold/Desktop/Kei/Rerun/"
samples <- c("DR_DR1","DR_DR2","DR_DR3","../Xenopus/NT_R1/","../Xenopus/NT_R2/","../Xenopus/NT_R3/")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("DNA","DNA","DNA","NT","NT","NT"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so <- sleuth_fit(so, ~timepoint, 'full')
#And reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_XNT_DR.csv")

model <- sleuth_wt(so, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_XNT_DR.csv")



#Didn't run???
#kal_files <- file.path("1-8_S14_v2", "abundance.h5")
#kal_files <- file.path("../Xenopus/NT_R2_v3", "abundance.h5")

#n_targets <- sapply(kal_files, function(file) {
#  ids <- rhdf5::h5read(file, "aux/ids")
#  n_targets <- length(ids)
#  n_targets
#})
#versions <- sapply(kal_files, function(file) {
#  file_name <- basename(file)
#  version <- rhdf5::h5read(file, "aux/kallisto_version")
#  version
#})


samples <- c("1-8_S14_v2","1-16_S16_v2","1-18_S18_v2","../Xenopus/NT_R1/","../Xenopus/NT_R2/","../Xenopus/NT_R3/")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("A-ctrl","A-ctrl","A-ctrl","NT","NT","NT"), stringsAsFactors=FALSE)
so2 <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, extra_bootstrap_summary = TRUE)
#Fit full model
so2 <- sleuth_fit(so2, ~timepoint, 'full')
#And reduced model
so2 <- sleuth_fit(so2, ~1, 'reduced')
so2 <- sleuth_lrt(so2, 'reduced', 'full')
sleuth_table_gene_2 <- sleuth_results(so2, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene_2, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_XNT_1-8.csv")
#models(so)
model <- sleuth_wt(so, 'timepointNT','full')
results_wt <- sleuth_results(model, 'timepointNT', test_type = 'wt')
write.table(results_wt, file = "/Users/christopherpenfold/Desktop/Kei/Rerun/sleuth_table_gene_all_b_XNT_1-8.csv")
