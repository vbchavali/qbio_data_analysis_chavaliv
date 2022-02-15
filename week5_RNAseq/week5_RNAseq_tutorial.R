library(TCGAbiolinks)
library(SummarizedExperiment)

#exercise 2.1 
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
GDCdownload(query) # only need to download the data once! Comment this out once you have completed it once
sum_exp <- GDCprepare(query)

#exercise 2.2
assays(sum_exp)$"HTSeq - Counts"[1:5, 1:5]
str(sum_exp)

#exercise 2.3 
dim(colData(sum_exp))
dim(rowData(sum_exp))
dim(assays(sum_exp)$"HTSeq - Counts")

#exercise 2.4 
str(colData(sum_exp))
head(colData(sum_exp))
str(assays(sum_exp)$"HTSeq - Counts")


#exercise 2.5
sum_exp$age_at_diagnosis
#It seems to be age_at_diagnosis. 

#exercise 2.6 
sum_exp$age_at_diagnosis[1:10]
str(sum_exp$age_at_diagnosis)
#It looks like the data has units of days alive. 

#exercise 2.7
colData(sum_exp)$age_at_diagnosis = colData(sum_exp)$age_at_diagnosis/365.25
  
#exercise 2.8
colData(sum_exp)$age_category = ifelse((sum_exp)$age_at_diagnosis < 50, "Young", 
                                       "Old")

#exercise 2.9 
head(rowData(sum_exp))
dim(rowData(sum_exp))

#exercise 2.10
"APC" %in% rowData(sum_exp)$external_gene_name
"BRAF" %in% rowData(sum_exp)$external_gene_name
print(rowData(sum_exp)$external_gene_name)

#exercise 2.11
assays(sum_exp)$"HTSeq - Counts"[20:25, 30:35]
rowData(sum_exp)[rowData(sum_exp)$external_gene_name == "APC",]
rowData(sum_exp)[rowData(sum_exp)$external_gene_name == "BRAF",]

#exercise 2.12 
geneA_id_mask = (rowData(sum_exp)$external_gene_name == "APC")
geneB_id_mask = (rowData(sum_exp)$external_gene_name == "BRAF")
sum(geneA_id_mask) #sum is 1 b/c there is only one of that gene
sum(geneB_id_mask) #sum is 1 b/c there is only one of that gene

##Is this right? 
ensembl_geneA = rowData(sum_exp)$ensembl_gene_id[geneA_id_mask]
ensembl_geneB = rowData(sum_exp)$ensembl_gene_id[geneB_id_mask]
print(ensembl_geneA)
print(ensembl_geneB)

#exercise 2.13 
dim(assays(sum_exp)$"HTSeq - Counts")
#Ensemble ID its a row in the assays table. 

#exercise 2.14
min(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA,])
summary(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB,])

#exercise 2.15 
data_for_geneA = assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA,]
data_for_geneB = assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB,]
plot(data_for_geneA,
     data_for_geneB,
     xlab = "APC Gene", # remember to rename axes!
     ylab = "BRAF Gene"
)


#exercise 2.16 
bool_age_na = is.na(colData(sum_exp)$age_category)
num_na = sum(bool_age_na)
num_na
#I do not foresee an issue, b/c there are 0 NAs. 

#exercise 2.17 
age_cat_no_NAs = colData(sum_exp)$age_category[bool_age_na]
print(age_cat_no_NAs)

#exercise 2.18 
length(age_cat_no_NAs)
dim( colData(sum_exp) ) #gives number of rows then number of columns
dim( colData(sum_exp) )[1] #gives number of rows
dim( colData(sum_exp) )[2]

dim( colData(sum_exp) )[1] + length(age_cat_no_NAs) == num_na

#exercise 2.19 
dim(assays(sum_exp)$"HTSeq - Counts")
length(assays(sum_exp)$"HTSeq - Counts"[1,])

#exercise 2.20 
identical(rownames(colData(sum_exp)), colnames(assays(sum_exp)$"HTSeq - Counts")
          )
gene_counts = assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, !bool_age_na]
print(gene_counts)

#exercise 2.21 
length(age_cat_no_NAs) == length(gene_counts)

#exercise 2.22
boxplot(gene_counts ~ age_cat_no_NAs, 
        xlab = "Age Category", 
        ylab = "APC Gene Counts")

#exercise 3.1 
## Access using assays(sum_exp)$"HTSeq - Counts". Rows are gene IDs and columns 
## are patient IDs. 

#exercise 3.2 
## Access using rowData(sum_exp). Rows are ensembl_gene_id. Rows of count df and 
## rows of rowData are the same thing. 

#exercise 3.3 
## Access using colData(sum_exp). Columns are patient IDs. Additional, informati
## on about patients is stored here. The columns should be the same between 
## count df and the colData df. 

