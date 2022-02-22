library(TCGAbiolinks)
library(DESeq2)
library(SummarizedExperiment)


#exercise 1.1 
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
sum_exp <- GDCprepare(query)

#exercise 1.2
assays(sum_exp)$"HTSeq - Counts"
counts = assays(sum_exp)$"HTSeq - Counts"
bool_age_na = is.na(colData(sum_exp)$age_at_index)
patient_data = colData(sum_exp)[!bool_age_na, ]
counts = counts[,!bool_age_na]
dim(bool_age_na)
dim(patient_data)
dim(counts)
length(bool_age_na)

patient_data$age_category = ifelse(patient_data$
                                     age_at_index < 50,
                                   "Young", "Old")
patient_data$age_category = factor(patient_data$age_category, levels = 
                                     c("Young", "Old"))

#exercise 1.3 

##First step does not work. 
if (all(rownames(counts) == names(rowRanges(sum_exp)))){
  print("Changing row names!")
  rownames(counts) = rowRanges(sum_exp)$external_gene_name
}
count_row_sums = rowSums(counts)
counts = counts[count_row_sums >= 10,]
dim(counts)
dim(patient_data)

#exercise 2.1 
dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = patient_data,
                             design = ~age_category)

dds_obj = DESeq(dds)
resultsNames(dds_obj)  # see what comparisons got run

# get the young vs. old comparison
results = results(dds_obj, format = "DataFrame", contrast = c("age_category", "Young", "Old"))

#exercise 2.2
dim(results)
str(results)

#exercise 2.3
row_order = order(results$padj)
results = results [row_order,]
head(results, 20)

## ENSG00000202198 is 7SKRNA. It is more expressed in older populations as seen 
## by the negative log2Foldchange number. The gene is studied for miRNA and 
## inhibitory RNAs. 

#exercise 2.4 
log2FoldChange_threshold = 1
padj_threshold = 0.05
results_log2FCg = results[results$log2FoldChange > log2FoldChange_threshold,]
results_log2FCl = results[results$log2FoldChange < -log2FoldChange_threshold,]
results_padjl = results[results$padj < padj_threshold,]

#exercise 2.5
fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05

# fill in your plot code here!
# be sure to relabel the axes!
# note: you can perform the log transformation directly in the plot function
plot(x = results$log2FoldChange,
     y = -log10(results$padj),
     xlab = "Fold Change (Young over Old)", # be sure the specify that it's young over old!
     ylab = "Negative Log 10 of P-Adjusted Value",
     pch = 20) # smaller solid circles
abline(v=c(-log2(fc_threshold), log2(fc_threshold)), h= c(-log10(p_threshold)), 
       col="green")

# these lines put the lines on the plot
# abline() plots straight lines on an R plot.
# v argument is for a vertical line, h argument is for a horizontal line, col 
# argument is color
fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05


fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05


library(ggplot2)

volcano_plot = ggplot(data = data.frame(results), 
                      aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = ifelse(log2FoldChange < -1 & padj < 0.05, "lower in young",
                                ifelse(log2FoldChange > 1 & padj < 0.05, "higher in young", "NS"))),
             size = 0.5) + 
  theme_minimal() + # make things pretty +
  theme(legend.title = element_blank()) + 
  # next 2 lines draw lines at the thresholds
  geom_vline(xintercept=c(-log2(fc_threshold), log2(fc_threshold)), color="green") + 
  geom_hline(yintercept=-log10(p_threshold), color="green") + 
  scale_color_discrete(type=c("red", "blue", "black")) +
  labs(x = "log2 Fold Change (Young/Old)",
       y = "-log10 Adjusted p-value")


volcano_plot


#exercise 2.6
write.csv(x = results,
          file = '/Users/vbchavali/Documents/qbio_data_analysis_chavaliv/week6_DESeq2/results.csv',row.names = FALSE)


