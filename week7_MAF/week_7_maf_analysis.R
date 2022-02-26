BiocManager::install("maftools")
library(maftools)
library(TCGAbiolinks)
library(ggplot2)


#exercise 1.1 
# Done 

#exercise 1.2 
# replace the path!
# fread is a way to read in a data frame, but it's much faster than the built-in read.csv()
# by default, it reads it in as a slightly different data type than a data frame, so we set the data.table flag to false
clinic <- data.table::fread("/Users/vbchavali/Documents/qbio_data_analysis_chavaliv/analysis_data/clinic.csv",
                            data.table = F)

# rename the patient barcode to make it work with maftools
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

#exercise 1.3 
length(colnames(clinic))
# 76 data values in colnames(clinic)
str(colnames(clinic) == "Tumor_Sample_Barcode")
print(clinic$Tumor_Sample_Barcode)
# It contains the barcode ID of the data. 
print(colnames(clinic) == "Tumor_Sample_Barcode")
# There is only one true. 

#exercise 1.4
mutation_query <- GDCquery_Maf(tumor = "COAD", 
                               pipeline = "mutect2",
                               save.csv = TRUE)

maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)
#exercise 1.5 
# cleared environent 
# found appropriate csv file 
maf_dataframe = data.table::fread("/Users/vbchavali/Documents/qbio_data_analysis_chavaliv/analysis_data/GDCdata/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv",
                                  data.table = F) 
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
clinic <- data.table::fread("/Users/vbchavali/Documents/qbio_data_analysis_chavaliv/analysis_data/clinic.csv",
                            data.table = F)
maf_object <- read.maf(maf = maf_dataframe, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)

#exercise 2.1
maf_object
str(maf_object)
maf_object@clinical.data
maf_object@data
print(colnames(maf_object@clinical.data) == colnames(maf_object@data))
str(maf_object@data)
str(maf_object@clinical.data)
#Tumor_Sample_Barcode matches between these two variables. 

#exercise 3.1 
oncoplot(maf = maf_object,
         top = 20) 

ggsave("/Users/vbchavali/Documents/qbio_data_analysis_chavaliv/week7_maf/mutation.png")

#exercise 3.2 
# I chose TP53. This geneis a tumor suppressor gene, so it is not directly 
# related to cell growth or prolitferation. 

#exercise 3.3
clinic = maf_object@clinical.data
clinic$age_category = ifelse(clinic$age_at_initial_pathologic_diagnosis < 50, "Young", "Old")
young_patients_ids = clinic$Tumor_Sample_Barcode[clinic$age_category == "Young"]
young_maf = subsetMaf(maf = maf_object,tsb = young_patients_ids)
old_patients_ids = clinic$Tumor_Sample_Barcode[clinic$age_category == "Old"]
old_maf = subsetMaf(maf = maf_object,tsb = old_patients_ids)
rm(old_patient_ids)

#exercise 3.4 
coOncoplot(m1 = young_maf, 
           m2 = old_maf, 
           m1Name = "Young Patients Mutation Data", 
           m2Name = "Old Patient Mutation Data")

ggsave("/Users/vbchavali/Documents/qbio_data_analysis_chavaliv/week7_maf/youngvold.png")

#exercise 3.5 
# The genes do not follow the pattern I expected. 

#exercise 3.6 
lollipopPlot(maf_object, gene = "TP53")

ggsave("/Users/vbchavali/Documents/qbio_data_analysis_chavaliv/week7_maf/tp53.png")

#exercise 3.7 
lollipopPlot2(m1 = young_maf, 
              m2 = old_maf, 
              m1_name = "Young Patients Mutation Data", 
              m2_name = "Old Patient Mutation Data",
              gene = "TP53")
ggsave("/Users/vbchavali/Documents/qbio_data_analysis_chavaliv/week7_maf/tp53youngvoldlolli.png")

#exercise 3.8 
# 75% of samples without a mutation in A, have a mutation in B. 

#exercise 3.9 
# b = 7 c = 2 e = 37 f = 42 d = 35

#exercise 3.10 
# geneA is TP53
geneA_maf <- subsetMaf(maf = maf_object,
                       genes = "TP53")

# geneB is KRAS
geneB_maf <- subsetMaf(maf = maf_object, 
                       genes = "KRAS")
print(geneA_maf)
print(geneB_maf)

#exercise 3.11 
# subsetMAF() splits the data based on the selected paramgeters that are inputted. 
# FINISH THIS Q 

#exericise 3.12 
# 1. Access the barcodes of the patients with mutations in genes A and B
# bc stands for barcode
mut_bc_geneA = geneA_maf@clinical.data$Tumor_Sample_Barcode
mut_bc_geneB = geneB_maf@clinical.data$Tumor_Sample_Barcode

# 2. Get the lengths of these two vectors
num_mut_geneA = length(mut_bc_geneA)
num_mut_geneB = length(mut_bc_geneB)
print(num_mut_geneA) #213
print(num_mut_geneB) #163

# 3. Fill in the intersect here! Then get the nubmer of patients
mut_bc_geneAB = intersect(mut_bc_geneA, mut_bc_geneB)
num_mut_geneAB = length(mut_bc_geneAB)
print(num_mut_geneAB) #78

#exercise 3.13 
num_mut_geneA_only = num_mut_geneA - num_mut_geneAB
num_mut_geneB_only = num_mut_geneB - num_mut_geneAB
print(num_mut_geneA_only) #135
print(num_mut_geneB_only) #85

maf_object@clinical.data





#exericise 3.14
num_neither_mutation = 397 - num_mut_geneA - num_mut_geneB + num_mut_geneAB

contig_table <- matrix(c(num_mut_geneAB, 
                         num_mut_geneB_only,
                         num_mut_geneA_only,
                         num_neither_mutation), 
                       nrow=2)

# view the contingency table
contig_table

#exercise 3.15 
fe_results <- fisher.test(contig_table)
fe_results

#P-Value is 0.06543. We cannot reject the null hypothesis. 

#exercise 3.16 
# Let's write a function to make our code easier to understand
# Note that you can speed this up a lot more by optimizing this GetBarcodes function
if (!require(foreach)) install.packages("foreach")
if (!require(doParallel)) install.packages("doParallel")
if (!require(reshape2)) install.packages("reshape2")

library(foreach)
library(doParallel)

# how many computations you can do at once
registerDoParallel(detectCores())

# turns the data frame into the following structure:
# rows: patient IDs
# columns: genes
# cells: number of mutations in that gene
maf_object@data$Index = 1
maf_cast = reshape2::dcast(maf_object@data, Tumor_Sample_Barcode ~ Hugo_Symbol)

# Compute the signifiances
ComputePairwiseMutation_Parallel = function(maf_cast, geneA){
  
  # Tumor_Sample_Barcode is the first column name
  all_genes <- unique(colnames(maf_cast)[-1])
  # TRUE if has mutation, FALSE otherwise
  # the important speedup: this is much faster than identifying names
  # and calling intersect
  has_mutA = (maf_cast[[geneA]] != 0)
  n = length(all_genes)
  
  # with parallelization
  foreach(i = 1:n, .combine = "c") %dopar% {
    gene = all_genes[i]
    has_mutB = (maf_cast[[gene]] != 0)
    
    # we can use the built-in table() function
    contig_table = table(has_mutA, has_mutB)
    fe_results <- fisher.test(contig_table)
    pval = fe_results$p.value
    
    names(pval) = gene
    return(pval)
  }
}

p_vals = ComputePairwiseMutation_Parallel(maf_cast, "KRAS")

# 1. reformat the results by creating a data frame
results = data.frame(
  gene = names(p_vals),
  pvalues = p_vals
)


# 2. adjust p-values and create the padj column
# BH is the Benjamini-Hochberg method
results$padj = p.adjust(p = results$pvalues, method = "BH")

# 3. sort your results data frame
row_order = order(results$pvalues)
results = results[row_order,]

# 4. create your sig_results data frame
sig_results = results[results$padj < 0.05,]

data.table::fwrite(results,
                   file = "Users/vbchavali/Documents/qbio_data_analysis_chavaliv/week7_maf//KRAS_pairwise_results.csv")
                   
                   