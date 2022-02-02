if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")
if(!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")


clin_query <- GDCquery(project = "TCGA-COAD", 
                       data.category = "Clinical",
                       file.type = "xml")
##Only need this once. 
GDCdownload(clin_query) 

clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
print (clinic)

#exercise 1.1
str(clinic)
head(clinic)

#exercise 1.2 
colnames(clinic)
clinic$bcr_patient_barcode

#exercise 2.1
plot(clinic$age_at_initial_pathologic_diagnosis, clinic$weight, xlab = "Age",
     ylab = "Weight")

#exercise 2.2
boxplot(clinic$age_at_initial_pathologic_diagnosis ~ clinic$race_list)
unique(clinic$race_list)
par(mar=c(10,1,1,1))
boxplot(clinic$age_at_initial_pathologic_diagnosis ~ clinic$race_list, las = 2, 
        cex.axis = 0.5)

#exercise 2.3
clinic$race_list
clinic$race_list = ifelse(clinic$race_list == "", "No Data", clinic$race_list)
clinic$race_list

#exercise 2.4
min(clinic$age_at_initial_pathologic_diagnosis)
max(clinic$age_at_initial_pathologic_diagnosis)
mean(clinic$age_at_initial_pathologic_diagnosis)
median(clinic$age_at_initial_pathologic_diagnosis)
summary(clinic$age_at_initial_pathologic_diagnosis)

#exercise 2.5
young = clinic[clinic$age_at_initial_pathologic_diagnosis < 50,]
old = clinic[clinic$age_at_initial_pathologic_diagnosis >= 50,]
dim(young)
dim(old)
dim(clinic)
#The rows do total up. 

#exercise 2.6 
young_pt_ids = young$patient_id
length(young_pt_ids)
old_pt_ids = old$patient_id
length(old_pt_ids)

#exercise 2.7
clinic$age_category = ifelse(clinic$age_at_initial_pathologic_diagnosis < 50, 
                             "Young", "Old")
head(clinic$age_category)

#exercise 2.8
clinic[1,1] #"TCGA-3L-AA1B"
clinic[1,] #Accesses first row and all columns
clinic[2:5,] #Accesses rows 2-5
clinic[,3] #Accesses all data in column 3

#exercise 2.9 
young_clinic = clinic[clinic$age_category == 'Young',]
print(young_clinic)
old_clinic = clinic[clinic$age_category == 'Old',]

#exercise 2.10
young_clinic_one_line = clinic[clinic$age_at_initial_pathologic_diagnosis < 50,]
identical(dim(young_clinic), dim(young_clinic_one_line))

#exercise 3
install.packages("survival")
install.packages("survminer")
library("survival")
library("survminer")

#exercise 3.1
clinic$survival_time = ifelse(is.na(clinic$days_to_death),
                               clinic$days_to_last_follow_up, 
                               clinic$days_to_death )
#exercise 3.2
clinic$death_event = ifelse(clinic$vital_status == "Alive" , 0, 1 )

#exercise 3.3
# We initialize a 'survival' object first, which contains the data we need.
surv_object <- Surv(time = clinic$days_to_death, 
                    event = clinic$death_event)

# We then create a fit object
race_fit <- surv_fit( surv_object ~ clinic$race_list, data = clinic )

#the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot = ggsurvplot(race_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")

p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p

# save the plot as a png
# if you want to present the data, change the strata labels!
ggsave("../week4_clinical/kmplot_by_race.png", plot = p, width = 12, height = 9)
##It could be improved by making it more human-readable. The shape of the graph
##definitely makes it difficult.

#exercise 4.1
write.csv(clinic, "/Users/vbchavali/Desktop/coad_clinical_data.csv", 
          row.names = F)
clinic_read_in <- read.csv("/Users/vbchavali/Desktop/coad_clinical_data.csv")
