#Vamsi Chavali 
#Midsemester Project 


#libraries 
library(TCGAbiolinks)
library(DESeq2)
library(SummarizedExperiment)
library(maftools)
library(ggplot2)
library(survival)
library(survminer)

clinic <- data.table::fread("/Users/vbchavali/Documents/qbio_data_analysis_chavaliv/analysis_data/clinic.csv",
                            data.table = F)
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

gender_mask = clinic$gender == ""
print(gender_mask)
clinic_mod = clinic
clinic_mod$survival_time = ifelse(is.na(clinic$days_to_death),
                                  clinic$days_to_last_follow_up, 
                                  clinic$days_to_death )
clinic_mod$death_event = ifelse(clinic$vital_status == "Alive" , 0, 1 )
clinic_mod = clinic_mod[!gender_mask, ]

surv_object <- Surv(time = clinic_mod$days_to_death, 
                    event = clinic_mod$death_event)

dim(clinic_mod)

gender_fit <- surv_fit(surv_object ~ clinic_mod$gender,data = clinic_mod)

survplot = ggsurvplot(gender_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")

p_gender = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p_gender

boxplot(clinic_mod$survival_time ~ clinic$gender, xlab = "Gender", ylab = "Survival Time")

#Mutation Analysis with MafTools and Plots
maf_dataframe = data.table::fread("/Users/vbchavali/Documents/qbio_data_analysis_chavaliv/analysis_data/GDCdata/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv",
                                  data.table = F) 
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

maf_object <- read.maf(maf = maf_dataframe, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)

oncoplot(maf = maf_object, top = 30)

male_patients_ids = clinic$Tumor_Sample_Barcode[clinic$gender == "MALE"]
male_maf = subsetMaf(maf = maf_object,tsb = male_patients_ids)
female_patients_ids = clinic$Tumor_Sample_Barcode[clinic$gender == "FEMALE"]
female_maf = subsetMaf(maf = maf_object,tsb = female_patients_ids)

coOncoplot(m1 = male_maf, 
           m2 = female_maf, 
           m1Name = "Male Patients Mutation Data", 
           m2Name = "Female Patient Mutation Data")

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
sum_exp <- GDCprepare(query)

assays(sum_exp)$"HTSeq - Counts"
counts = assays(sum_exp)$"HTSeq - Counts"
bool_gender_na = is.na(colData(sum_exp)$gender)
patient_data = colData(sum_exp)[!bool_gender_na, ]
counts = counts[,!bool_gender_na]

patient_data$gender = factor(patient_data$gender, levels = 
                                     c("Male", "Female"))

#BRAF Plots 
mask_braf = rowData(sum_exp)$external_gene_name == "BRAF"
ensembl_braf_id = rowData(sum_exp)$ensembl_gene_id[mask_braf]
maf_braf = subsetMaf(maf = maf_object,genes = "BRAF")
maf_braf_patient_IDs = maf_braf@clinical.data$Tumor_Sample_Barcode
print(maf_braf_patient_IDs)
mask_braf = clinic_mod$Tumor_Sample_Barcode %in% c(maf_braf_patient_IDs)
clinic_mod_braf = clinic_mod[mask_braf,]
clinic_mod_braf$survival_time = ifelse(is.na(clinic_mod_braf$days_to_death),
                                  clinic_mod_braf$days_to_last_follow_up, 
                                  clinic_mod_braf$days_to_death )
clinic_mod_braf$death_event = ifelse(clinic_mod_braf$vital_status == "Alive" , 0, 1 )
#gender_mask_braf = clinic_mod_braf$gender == ""
#clinic_mod_braf = clinic_mod_braf[,!gender_mask_braf]

braf_object <- Surv(time = clinic_mod_braf$days_to_death, 
                    event = clinic_mod_braf$death_event)

braf_fit <- surv_fit(braf_object ~ clinic_mod_braf$gender, data = clinic_mod_braf)

brafplot = ggsurvplot(braf_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")
p_braf = brafplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p_braf

lollipopPlot2(m1 = male_maf, 
              m2 = female_maf, 
              m1_name = "Male Patients Mutation Data", 
              m2_name = "Female Patient Mutation Data",
              gene = "BRAF")

#KRAS Plots 
mask_kras = rowData(sum_exp)$external_gene_name == "KRAS"
ensembl_kras_id = rowData(sum_exp)$ensembl_gene_id[mask_kras]
maf_kras = subsetMaf(maf = maf_object,genes = "KRAS")
maf_kras_patient_IDs = maf_kras@clinical.data$Tumor_Sample_Barcode
mask_kras = clinic_mod$Tumor_Sample_Barcode %in% c(maf_kras_patient_IDs)
clinic_mod_kras = clinic_mod[mask_kras, ]
clinic_mod_kras$survival_time = ifelse(is.na(clinic_mod_kras$days_to_death),
                                       clinic_mod_kras$days_to_last_follow_up, 
                                       clinic_mod_kras$days_to_death )
clinic_mod_kras$death_event = ifelse(clinic_mod_kras$vital_status == "Alive" , 0, 1 )


kras_object <- Surv(time = clinic_mod_kras$days_to_death, 
                    event = clinic_mod_kras$death_event)

kras_fit <- surv_fit(kras_object ~ clinic_mod_kras$gender,data = clinic_mod_kras)

krasplot = ggsurvplot(kras_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")
p_kras = krasplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p_kras

lollipopPlot2(m1 = male_maf, 
              m2 = female_maf, 
              m1_name = "Male Patients Mutation Data", 
              m2_name = "Female Patient Mutation Data",
              gene = "KRAS")


#PIK3CA Plots 
mask_pik3ca = rowData(sum_exp)$external_gene_name == "PIK3CA"
ensembl_pik3ca_id = rowData(sum_exp)$ensembl_gene_id[mask_pik3ca]
maf_pik3ca = subsetMaf(maf = maf_object,genes = "PIK3CA")
maf_pik3ca_patient_IDs = maf_pik3ca@clinical.data$Tumor_Sample_Barcode
mask_pik3ca = clinic_mod$Tumor_Sample_Barcode %in% c(maf_pik3ca_patient_IDs)
clinic_mod_pik3ca = clinic_mod[mask_pik3ca, ]
clinic_mod_pik3ca$survival_time = ifelse(is.na(clinic_mod_pik3ca$days_to_death),
                                       clinic_mod_pik3ca$days_to_last_follow_up, 
                                       clinic_mod_pik3ca$days_to_death )
clinic_mod_pik3ca$death_event = ifelse(clinic_mod_pik3ca$vital_status == "Alive" , 0, 1 )

pik3ca_object <- Surv(time = clinic_mod_pik3ca$days_to_death, 
                    event = clinic_mod_pik3ca$death_event)

pik3ca_fit <- surv_fit(pik3ca_object ~ clinic_mod_pik3ca$gender,data = clinic_mod_pik3ca)

pik3caplot = ggsurvplot(pik3ca_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")
p_pik3ca = pik3caplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p_pik3ca

lollipopPlot2(m1 = male_maf, 
              m2 = female_maf, 
              m1_name = "Male Patients Mutation Data", 
              m2_name = "Female Patient Mutation Data",
              gene = "PIK3CA")


#SMAD2 Plots 
mask_smad2 = rowData(sum_exp)$external_gene_name == "SMAD2"
ensembl_smad2_id = rowData(sum_exp)$ensembl_gene_id[mask_smad2]
maf_smad2 = subsetMaf(maf = maf_object,genes = "SMAD2")
maf_smad2_patient_IDs = maf_smad2@clinical.data$Tumor_Sample_Barcode
mask_smad2 = clinic_mod$Tumor_Sample_Barcode %in% c(maf_smad2_patient_IDs)
clinic_mod_smad2 = clinic_mod[mask_smad2, ]
clinic_mod_smad2$survival_time = ifelse(is.na(clinic_mod_smad2$days_to_death),
                                         clinic_mod_smad2$days_to_last_follow_up, 
                                         clinic_mod_smad2$days_to_death )
clinic_mod_smad2$death_event = ifelse(clinic_mod_smad2$vital_status == "Alive" , 0, 1 )

smad2_object <- Surv(time = clinic_mod_smad2$days_to_death, 
                      event = clinic_mod_smad2$death_event)

smad2_fit <- surv_fit(smad2_object ~ clinic_mod_smad2$gender,data = clinic_mod_smad2)

smad2plot = ggsurvplot(smad2_fit, 
                        pval=TRUE, 
                        ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                        legend = "right")
p_smad2 = smad2plot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p_smad2

lollipopPlot2(m1 = male_maf, 
              m2 = female_maf, 
              m1_name = "Male Patients Mutation Data", 
              m2_name = "Female Patient Mutation Data",
              gene = "SMAD2")

#SMAD4 Plots 
mask_smad4 = rowData(sum_exp)$external_gene_name == "SMAD4"
ensembl_smad4_id = rowData(sum_exp)$ensembl_gene_id[mask_smad4]
maf_smad4 = subsetMaf(maf = maf_object,genes = "SMAD4")
maf_smad4_patient_IDs = maf_smad4@clinical.data$Tumor_Sample_Barcode
mask_smad4 = clinic_mod$Tumor_Sample_Barcode %in% c(maf_smad4_patient_IDs)
clinic_mod_smad4 = clinic_mod[mask_smad4, ]
clinic_mod_smad4$survival_time = ifelse(is.na(clinic_mod_smad4$days_to_death),
                                        clinic_mod_smad4$days_to_last_follow_up, 
                                        clinic_mod_smad4$days_to_death )
clinic_mod_smad4$death_event = ifelse(clinic_mod_smad4$vital_status == "Alive" , 0, 1 )


smad4_object <- Surv(time = clinic_mod_smad4$days_to_death, 
                     event = clinic_mod_smad4$death_event)

smad4_fit <- surv_fit(smad4_object ~ clinic_mod_smad4$gender,data = clinic_mod_smad4)

smad4plot = ggsurvplot(smad4_fit, 
                       pval=TRUE, 
                       ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                       legend = "right")
p_smad4 = smad4plot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p_smad4

lollipopPlot2(m1 = male_maf, 
              m2 = female_maf, 
              m1_name = "Male Patients Mutation Data", 
              m2_name = "Female Patient Mutation Data",
              gene = "SMAD4")

#Differential Expression Analysis
fc_threshold = 2 # set a threshold of at least a 2 fold increase 
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05

counts_mod = counts[c(ensembl_braf_id, ensembl_kras_id, ensembl_pik3ca_id, 
                  ensembl_smad2_id, ensembl_smad4_id),]

dds = DESeqDataSetFromMatrix(countData = counts_mod,
                             colData = patient_data,
                             design = ~gender)

dds_obj = DESeq(dds)
resultsNames(dds_obj)
results = results(dds_obj, format = "DataFrame", contrast = c("gender", "male", "female"))
print(results)

volcano_plot = ggplot(data = data.frame(results), 
                      aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = ifelse(log2FoldChange < -1 & padj < 0.05, "lower in males",
                                ifelse(log2FoldChange > 1 & padj < 0.05, "higher in males", "GENES"))),
             size = 0.5) + 
  theme_minimal() + # make things pretty +
  theme(legend.title = element_blank()) + 
  # next 2 lines draw lines at the thresholds
  geom_vline(xintercept=c(-log2(fc_threshold), log2(fc_threshold)), color="green") + 
  geom_hline(yintercept=-log10(p_threshold), color="blue") + 
  scale_color_discrete(type=c("red", "blue", "black")) +
  labs(x = "log2 Fold Change (Male/Female)",
       y = "-log10 Adjusted p-value")


volcano_plot


