
#### Set Working directory ####
setwd("/projects/sunlab/Students_Work/Caspar_work/Thesis/UKB")

#### Loading package ####
library(pacman)
# Load pacman to import the function p_load()

p_load(tidyverse, ggrepel, lubridate, scales, matrixStats, reshape2, data.table, gtreg, gt, ggfortify, ggpubr, msigdbr, fgsea, Hmisc, VennDiagram)
# Load essential package for data management and GSEA analysis

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.20")
# Load Bioconductor

#### Enrollment dataset ####
data1 <- data.table::fread("/projects/sunlab/UKB/Data_clean/ukb27656.csv")
data1 <- as.data.frame(data1)
# Read enrollment dataset

data1 <- rename(data1, "ID"           = "eid",
                       "enrol_date"   = "53-0.0",
                       "lost_fu_date" = "191-0.0",
                       "age_enrolled" = "21022-0.0",
                       "birth_year"   = "34-0.0",
                       "birth_month"  = "52-0.0",
                       "smoking"      = "20116-0.0",
		       "sex"          = "31-0.0",
                       "ethnic"	      = "21000-0.0",
		       "BMI"	      = "21001-0.0")
# Rename variable name

data1$birthday_art <- as.Date(paste0(data1$birth_year, "-", data1$birth_month, "-", "15"), origin = "1900-01-01")
data1$age <- data1$age_enrolled
data1$ethnic[data1$ethnic %in% c(1,1001,1002,1003)]<- "White"
data1$ethnic[data1$ethnic %in% c(2,2001,2002,2003,2004)]<- "Mixed"
data1$ethnic[data1$ethnic %in% c(3001,3002,3003)]<- "SouthAsian"
data1$ethnic[data1$ethnic %in% c(3,3004,5)]<- "OtherAsian"
data1$ethnic[data1$ethnic %in% c(4,4001,4002,4003)]<- "Black"
data1$ethnic[data1$ethnic %in% c(6)]<- "Other ethnic group"
# Recategorize ethnicity into subtypes
 
data1$white <- ifelse(data1$ethnic == "White",1,0)
# Dicotomize ethnicity into white and non-white

data1$smokings <- ifelse((data1$smoking == 1|data1$smoking == 2), 1, 0)
# Dichotomize smoking history into ever or never smoke

data1 <- data1[, c("ID", "enrol_date", "lost_fu_date", "age_enrolled",
                   "birthday_art", "smokings","age", "sex", "white", "BMI", "birthday_art")]
# Subset enrollment dataset

print(dim(data1))
# 502536*11

#### Disease dataset ####
disease<-data.table::fread("/projects/sunlab/UKB/Core_Data/ukb675990.tab",stringsAsFactors = F,header=T)
disease<-as.data.frame(disease)
# Load disease dataset

print(dim(disease))
# 502316*2197

##Death_Date ####
field_Death <- ("40000")
for (i in 1:length(field_Death)){print(grep(field_Death[i], names(disease), value = T))}
# Find all death related events(columns) from the column names of the dataset
# "f.40000.0.0" "f.40000.1.0" 
# "f.40000.0.0 codes for death date
# "f.40000.1.0 codes for death cause
 
field_Death <- paste0("f.", field_Death,c(".0.0", ".1.0"))
print(field_Death)
# Create column names for death events
# "f.40000.0.0" "f.40000.1.0"

Death <- disease[, c("f.eid", field_Death)]
# Subset Death dataset

for (i in 2:ncol(Death)){Death[, i] <- as.numeric(as.Date(Death[, i], "%Y-%m-%d", origin = "1900-01-01"))}
# Convert death date into numeric form

Death[(Death$f.40000.0.0 != Death$f.40000.1.0) & (!is.na(Death$f.40000.0.0) & !is.na(Death$f.40000.1.0)), ]
# Subset the participants (rows) with discrepency in death date

Death$Death_date <- matrixStats::rowMins(as.matrix(Death[, -1], na.rm = T))
# Aggregate death dates by creating a single column with the earliest death date from the two field number

Death$Death_date <- as.Date(Death$Death_date, "%Y-%m-%d", origin = "1970-01-01")
# Convert the death date back to date form

Death <- (Death[Death$Death_date!= Inf, c("f.eid", "Death_date")])
# Subset the Death dateset with only participant number and death date

print(dim(Death)) ## 502316*2

#### First Occurence Pneumonia ####
field_pneumonia <- c("131442", #Date J11 first reported (influenza, virus not identified)
                     "131444", #Date J12 first reported (viral pneumonia, not elsewhere classified)
                     "131446", #Date J13 first reported (pneumonia due to streptococcus pneumoniae)
                     "131448", #Date J14 first reported (pneumonia due to haemophilus influenzae)
                     "131450", #Date J15 first reported (bacterial pneumonia, not elsewhere classified)
                     "131452", #Date J16 first reported (pneumonia due to other infectious organisms, not elsewhere classified)
                     "131454", #Date J17 first reported (pneumonia in diseases classified elsewhere)
                     "131456"  #Date J18 first reported (pneumonia, organism unspecified)
                    )
# Create pnuemonia case definition with ICD10 code

omit_code <- as.Date(c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03","1909-09-09", "2037-07-07"))
# Omit the incidence precede the date of birth or at date of birth
# "1901-01-01" codes for event date before participant's date of birth
# "1902-02-02" codes for event date matching participant's date of birth
# "1903-03-03" codes for event date after participant's date of birth and falls in the same calendar year as date of birth
# "1909-09-09" codes for missing or invalid events
# "2037-07-07" codes for event date in the future

for (i in 1:length(field_pneumonia)){print(grep(field_pneumonia[i], names(disease), value = T))}
# Find out all columns matching the case definition of pneumonia

field_pneumonia <- paste0("f.",field_pneumonia,".0.0")
# Create column name list

pneumonia <- disease[, c("f.eid", field_pneumonia)]
# Subset disease dataset with participants ID and columns coding pneumonia only

for (i in 2:ncol(pneumonia)){pneumonia[,i] <- as.numeric(as.Date(pneumonia[,i],"%Y-%m-%d",origin = "1900-01-01"))}
# Convert the date that pneumonia cases were first reported into numeric form

pneumonia$pneumonia_onset <- matrixStats::rowMins(as.matrix(pneumonia[,-1]),na.rm = T)
# Create a single column to store the earliest date of pnuemonia onset
 
pneumonia$pneumonia_onset <- as.Date(pneumonia$pneumonia_onset,"%Y-%m-%d",origin = "1970-01-01")
# Convert it back to date form

pneumonia[sapply(pneumonia, is.infinite)] <- NA
# Assign NA to infinite value

pneumonia <- pneumonia[!pneumonia$pneumonia_onset %in% omit_code, ]
# Remove the participants(rows) with the omit code

pneumonia <-(pneumonia[, c("f.eid","pneumonia_onset")])
# Create a new dateset with only participant ID and pneumonia onset date

print(dim(pneumonia)) 
# 502280*2

#### Discharge dataset ####
exclude <- read.csv("/projects/sunlab/UKB/General/Withdraw_w34031_2023-04-25.csv",stringsAsFactors=F,header=F)
# Load discharge dateset

#### Merge dataset ####
analysis <- merge(Death, pneumonia, by = "f.eid", all.y = T)
analysis <- merge(data1, analysis, by.x = "ID", by.y = "f.eid", all.x = T)
# Merge Death, pneumonia, and enrollment dateset

analysis <- analysis[!analysis$ID%in%exclude[,1],]
# Exclude participants discharged from the UKB study

analysis <- analysis[!is.na(analysis$age), ]
analysis <- analysis[!is.na(analysis$sex), ]
analysis <- analysis[!is.na(analysis$white), ]
analysis <- analysis[!is.na(analysis$smokings), ]
# Exclude the participatns with missing values in the covariates of the model

print(dim(analysis))
# 501468*13

#### Create prevalence_pneumonia ####
date <- c("enrol_date","lost_fu_date","Death_date")
# Define baseline date

for (i in 1:length(date)){analysis[,date[i]] <- as.Date(analysis[,date[i]],"%Y-%m-%d",origin = "1900-01-01")}
# Convert them to date form

table(analysis$pneumonia_onset > analysis$enrol_date)
# Define prevelant pneumonia cases as having pneumonia onset earlier than the enrollment date
# Pneumonia onset after baseline date will be defined as incident cases
# TRUE: 33668 (incident case) 
# FALSE: 20075 (prevalent case) 

length(na.omit(analysis$pneumonia_onset))
# Omit missing value of pneumonia onset
# 53743 total cases

analysis$pneumonia_exposure <- ifelse(is.na(analysis$pneumonia_onset)|analysis$pneumonia_onset > analysis$enrol_date,0,1)
# pneumonia exposure is defined as prevalent pneumonia cases

print(table(analysis$pneumonia_exposure))
# 20075 prevalent cases

print(table(analysis$pneumonia_exposure)/nrow(analysis)) 
# prevalent rate 0.04003296

#### Create pneumonia incident to enrollment (ite) and ite group ####
analysis$pneumonia_ite <- difftime(as.Date(analysis$pneumonia_onset), as.Date(analysis$enrol_date), unit = "days")
# Calculate pneumonia incidence time to enrollment

analysis$ite_group <- ifelse(analysis$pneumonia_exposure == 0, "control",
		      ifelse(analysis$pneumonia_exposure == 1 & analysis$pneumonia_ite < -3650, "long", "short"))
# Define pneumonia ITE group with the threshhold of 10 years
# Pneumonia ITE > 10 yr (< -3650 days) is defined as long, else short pneumonia ITE

analysis$ite_short <- ifelse(analysis$pneumonia_exposure == 0, "control",
                      ifelse(analysis$pneumonia_exposure == 1 & analysis$pneumonia_ite < -1825, "long", "short")) 
# Alternative defination of pneumonia short with ITE < 5 yr(> -1825 days)
# This definition was not applied in the following analysis

analysis$pneumonia_aai <- difftime(as.Date(analysis$pneumonia_onset), as.Date(analysis$birthday_art), unit = "days") / 365
# Create age at pneumonia incidence

#### Load Proteomics Dataset ####
load("/projects/sunlab/UKB/Proteomics/olink_data_675033_WideFormat.rda")
# Load .rda file

print(dim(data_wide))
# Check dataset dimension
# 55313*2924

#### Read dataset containing batch number and plate ID ####
sample_data <- data.table::fread("/projects/sunlab/UKB/Proteomics/ukb675033.tab")
# Read sample_data dataset

batch_number <- read.table("/projects/sunlab/Students_Work/Caspar_work/Thesis/UKB/olink_batch_number.dat", header = T)
# Read batch number dataset

batch_number$PlateID <- as.character(batch_number$PlateID)
# Convert plate ID into character to match the class of PlateID in sample_data dataset

print(head(batch_number$PlateID))
# "890000000001" "890000000002" "890000000003" "890000000004" "890000000005"
 
#### Modify sample dataset #####
sample_data <- sample_data %>% rename("ID" = "f.eid", "PlateID" = "f.30901.0.0")
# Rename the column of participant ID and plate ID

sample_data <- sample_data[, c("ID", "PlateID")]
# Subset sample_data dataset

sample_data$ID <- as.character(sample_data$ID)
# Convert participant ID into character form
 
print(dim(sample_data))
# Check dataset dimension
# 502357*2

#### Proteomic data QC: filter (protein missing rate > 0.2, individual call rate < 0.8) ####
data_wide$index <- gsub("_0", "", data_wide$index)
# Clean participant ID by replacing _0 with space

data_wide <- data_wide[!is.na(data_wide$index), ]
# Exclude rows with missing values

protein <- grep("result", colnames(data_wide), value = T)
# Find all protein columns with the prefix "result"
# Total 2923 proteins

exclude <- rowSums(is.na(data_wide[, protein])) / length(protein)    
# Create a vector to store individual call rate across all proteins
 
data_wide_1 <- data_wide[exclude <= 0.2, ]
# Create a new dataset to store all rows > 80% individual call rate
# Exclude all rows with a missing rate > 20%

protein <- grep("result", colnames(data_wide_1), value = T)
# Create a vector of all proteins' column name in data_wide_1 dataset

ex_pro <- colMeans(is.na(data_wide_1[,protein]))
# Calcualte the missing rate of each protein

data_pro <- data_wide_1[, protein %in% names(ex_pro[ex_pro <= 0.2])]
# Create a new dataset with proteins < 20% missing rate

protein <- grep("result", colnames(data_pro), value = T)
# Replace the vector after excluding proteins with a missing rate

print(length(protein))
# Check the number of protein remain
# Total 2919 protein

print(head(protein)) 
# "result.729" "result.730" "result.732" "result.733" "result.734"

index <- NA
# Create an empty vector to store protein code

for (i in 1 :length(protein)) {list <- strsplit(protein[i], split = "[.]")
			       index[i] <- list[[1]][2]}
# Create a for loop that extract the protein code from each protein name

print(head(index))
# "729" "730" "732" "733" "734" "737"

colnames(data_pro)[1] <- "ID"
# Rename the first column of data_pro as ID

data_pro <- merge(data_pro, analysis, by = "ID")
data_pro <- merge(data_pro, sample_data, by = "ID", all.x = T)
# Merge data_pro with analysis and sample_data dataset by ID

match <- match(data_pro$PlateID, batch_number$PlateID)
# Create a vector containing the row position match with the plate ID of each participant in batch_number dataset 

match_1 <- c()
# Create an empty vector to store the batch number corresponding to the plate ID

for (i in 1:length(match)) {match_1 <- c(match_1, batch_number$Batch[match[i]])}
# Create a for loop to extract the batch number corresponding to the plate ID

data_pro$batch <- match_1
# Create a new vector in data_pro to store the value in match_1 vector

print(dim(data_pro))
# Check the dimension of data_pro dataset
# 44134*2939

#### Check infection time to enrollment ####
ite_sum <- data_pro %>% group_by(pneumonia_exposure) %>%
                        summarise(N_case   = table(pneumonia_exposure),
                                  Mean_ite = as.numeric(mean(pneumonia_ite, na.rm = TRUE)),
                                  Med_ite  = as.numeric(median(pneumonia_ite, na.rm = TRUE), units = "days"),
                                  SD_ite   = as.numeric(sd(pneumonia_ite, na.rm = TRUE), units = "days"),
                                  IQR_ite  = as.numeric(IQR(pneumonia_ite, na.rm = TRUE), units = "days"),
                                  Min_ite  = as.numeric(min(pneumonia_ite), units = "days"),
                                  Max_ite  = as.numeric(max(pneumonia_ite), units = "days"),
				  Qtl1_ite = as.numeric(quantile(pneumonia_ite, probs = 0.25, na.rm = TRUE), units = "days"),
				  Qtl3_ite = as.numeric(quantile(pneumonia_ite, probs = 0.75, na.rm = TRUE), units = "days"))

# Create a summary table stratified by pneumonia_exposure(1|0) with the following parameters:
# Number of cases
# Mean of incidence time to enrollment
# Median of incidence time to enrollment
# Standard deviation of incidence time to enrollment
# Interquartile range of incidence time to enrollment
# Minimum of incidence time to enrollment
# Maximum of incidence time to enrollment
# 1st quartile of incidence time to enrollment
# 3rd quartile of incidence time to enrollment
 
print(ite_sum)
# print the summary table

#### Tabel shell ####
Table_pro <- data_pro %>%
	     select(age, sex, white, smokings, pneumonia_ite, pneumonia_exposure) %>%
	     gtreg::tbl_reg_summary(by = pneumonia_exposure, include = c(age, sex, white, smokings, pneumonia_ite)) %>%
	     as_gt() %>%
             gtsave("Table_pro.docx")
# Create a table with the following variables: age, sex, white, smoking status, pneumonia_ite, pneumonia_exposure
# Stratify these variables with pneumonia_exposure
# Save Table_pro

#### Plot incident time to enrollment ####
data_plot <- subset(data_pro, ite_group == "long" | ite_group == "short")
# Create a new dataset only with the rows with ite_group equals to long or short
# Excldue control (participants without pneumonia_exposure)

pneumonia_Hist_1 <- ggplot(data_plot, aes(x = as.numeric(pneumonia_ite), fill = ite_group)) + 
		    geom_histogram(binwidth = 200, position = "dodge") + scale_fill_manual(values = c("lightsalmon", "orangered")) + 
		    theme_bw() + labs(x = "Incident time to enrollment (Days)", y = "Frequency") + 
		    theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), 
                          panel.background = element_blank())
# Plot a histogram of pneumonia incident time to enrollment stratified by ite_group

ggsave("pneumonia_Hist_1.png") 
# Save the plot
 
## Plot age at pneumonia incidence
pneumonia_Hist_2 <- ggplot(data_plot, aes(x = as.numeric(pneumonia_aai), fill = ite_group)) +
                    geom_histogram(binwidth = 1, position = "dodge") + scale_fill_manual(values = c("lightsalmon", "orangered")) + 
                    theme_bw() + labs(x = "Age at pneumonia incidence (Years)", y = "Frequency") +
                    theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank())

# Plot a hidtogram of age at pneumonia incidence stratified by ite_group

ggsave("pneumonia_Hist_2.png")
# Save the plot

pneumonia_His_3 <- ggarrange(pneumonia_Hist_1, pneumonia_Hist_2, 
			     legend = "right", common.legend = T,
			     labels = c("a", "b"))
# Combine two histograms with labels and common legend

ggsave("pneumonia_Hist_3.png")
# Save the plot

## Modify Annotation dataset
Annotation <- data.table::fread("/projects/sunlab/UKB/Proteomics/UKB_coding143_Olink.tsv")
# Read annotation dataset

Annotation$Feature <- NA
Annotation$Name    <- NA
# Create two empty variables to store protein name and code

for (i in 1:nrow(Annotation)) {list <- strsplit(Annotation$meaning[i], ";")
                               Annotation$Feature[i] <- list[[1]][1]
                               Annotation$Name[i] <- list[[1]][2]}
# Create a for loop to extract the protein name and code from the column of meaning

head(rownames(Annotation))
# Check the row names of Annotation dataset
# "result.729" "result.730" "result.732" "result.733" "result.734"

Annotation$coding <- as.character(Annotation$coding)
# Convert the coding column into characters

Annotation <- Annotation[Annotation$coding %in% index, c("Feature", "Name", "coding")]
# Subset Annotation dataset containing the proteins (rows) match the vector of index (protein code)
 
Annotation <- Annotation[order(match(Annotation$coding, index)), , drop = F]
# Arrange Annotation dataset with the order of protein code in index

for (i in 1:nrow(Annotation)) {rownames(Annotation)[i] <- paste0("result.", Annotation$coding[i])}
# Create a for loop to replace the row name of each protein with the same format in protein vector

#### Inverse Normal Transformation (#Phenotypes with missing values have rank of NA; smallest value rank=1) ####
for (i in 1:length(protein)){
 data_pro[,protein[i]][is.na(data_pro[,protein[i]])] <- min(data_pro[,protein[i]],na.rm = T)
 data_pro$rank<-rank(data_pro[,protein[i]],na.last="keep")
 data_pro$rank_P<-(data_pro$rank-0.5)/length(na.omit(data_pro[,protein[i]]))
 data_pro[,protein[i]]=qnorm(data_pro$rank_P)
 data_pro$rank<-NULL
 data_pro$rank_P<-NULL}
# Create a for loop to transform each protein
# Impute missing value in each protein as the minimum
# Rank the rows within each protein individually and retain all missing value
# Calculate the percentile of the values in each protein (subtract the rank with 0.5 to avoid infinite Z value)
# Apply INT to the values in each protein based on their p-values
# Remove rank and rank_P column 

#### Variable Modification ####
data_pro$sex <- factor(data_pro$sex)
# Convert sex into factor variable 
# All character variables will be converted into factor variable

data_pro$white <- factor(data_pro$white)
# Convert white into factor variable

data_pro$smokings <- factor(data_pro$smokings)
# Convert smoking status into factor variable

data_pro$pneumonia_exposure <- factor(data_pro$pneumonia_exposure)
# Convert pneumonia_expousre into factor variable

data_pro$ite_group <- factor(data_pro$ite_group)
# Convert ite_group into factor variable

data_pro$ite_short <- factor(data_pro$ite_short)
# Convert ite_short into factor variable

data_pro$agesq <- data_pro$age^2
# Create age square variable

save(data_pro, file = "/projects/sunlab/Students_Work/Caspar_work/Thesis/UKB/data_pro.rda")
# Save data_pro dataset

#### Create short and long pneumonia ITE dataset for sensitivity analysis ####
data_pro_short <- subset(data_pro, ite_group == "control" | ite_group == "short")
# Create a subset of data_pro dataset containing control and short pneumonia ITE cases

data_pro_long  <- subset(data_pro, ite_group == "control" | ite_group == "long")
# Create a subset of data_pro dataset containing control and long pneumonia ITE cases

data_pro_short$ite_group <- factor(data_pro_short$ite_group)
# Convert ite_group in data_pro_short into factor variable

data_pro_long$ite_group  <- factor(data_pro_long$ite_group)
# Convert ite_group in data_pro_long into factor variable

for (i in 1:length(protein)){
 data_pro_short[,protein[i]][is.na(data_pro_short[,protein[i]])] <- min(data_pro_short[,protein[i]],na.rm = T)
 data_pro_short$rank<-rank(data_pro_short[,protein[i]],na.last="keep")
 data_pro_short$rank_P<-(data_pro_short$rank-0.5)/length(na.omit(data_pro_short[,protein[i]]))
 data_pro_short[,protein[i]]=qnorm(data_pro_short$rank_P)
 data_pro_short$rank<-NULL
 data_pro_short$rank_P<-NULL}
# INT for the proteins in data_pro_short

for (i in 1:length(protein)){
 data_pro_long[,protein[i]][is.na(data_pro_long[,protein[i]])] <- min(data_pro_long[,protein[i]],na.rm = T)
 data_pro_long$rank<-rank(data_pro_long[,protein[i]],na.last="keep")
 data_pro_long$rank_P<-(data_pro_long$rank-0.5)/length(na.omit(data_pro_long[,protein[i]]))
 data_pro_long[,protein[i]]=qnorm(data_pro_long$rank_P)
 data_pro_long$rank<-NULL
 data_pro_long$rank_P<-NULL}
# INT for the proteins in  data_pro_long

Table_pro_short <- data_pro_short %>%
          select(age, sex, white, smokings, ite_group) %>%
          gtreg::tbl_reg_summary(by = ite_group, include = c(age, sex, white, smokings)) %>%
          as_gt() %>%
          gtsave("Table_pro_short.docx")
# Create table with same variables in Table_pro for data_pro_short and save it  

Table_pro_long <- data_pro_long %>%
          select(age, sex, white, smokings, ite_group) %>%
          gtreg::tbl_reg_summary(by = ite_group, include = c(age, sex, white, smokings)) %>%
          as_gt() %>%
          gtsave("Table_pro_long.docx")
# Create table with same variables in Table_pro for data_pro_short and save it

save(data_pro_short, file = "/projects/sunlab/Students_Work/Caspar_work/Thesis/UKB/data_pro_short.rda")
save(data_pro_long , file = "/projects/sunlab/Students_Work/Caspar_work/Thesis/UKB/data_pro_long.rda")
# Save two datasets

#### Linear Model Construction ####
# lm model
lm_pro <- function(data, protein, mode, file) {
# Create a function to run linear model with different datasets

table <- tibble(protein = character(),
                beta    = numeric(),
                SE      = numeric(),
                pval    = numeric())
# Create an empty datatable to store summary statistics
# protein name
# beta coefficient of pneumonia exposure, ite_group for short, ite_group for long
# standard error of beta coefficient
# p-value of beta coefficient

if (mode == "full") {pneumonia_exposure <-  "pneumonia_exposure"}
# pneumonia_exposure for full model

else if (mode == "short") {pneumonia_exposure <- "ite_group"}
# Change variable to ite_group for short pneumonia ITE model

else if (mode == "long")  {pneumonia_exposure <- "ite_group"}
# Change variable to ite_group for long pneumonia ITE model

else {stop("Invalid")}
# Terminate function

for (i in 1:length(protein)) {
formula <- paste0(protein[i], " ~ ", pneumonia_exposure, " + age + sex + age:sex + agesq + agesq:sex + white + smokings + batch")
model <- lm(data = data, as.formula(formula))
coeff <- coef(summary(model))
table <- rbind(table, tibble(Protein = Annotation$Feature[i],
                             beta    = coeff[2, 1],
                             SE      = coeff[2, 2],
                             pval    = coeff[2, 4]))}
# Create a for loop to execute lm function to analyze each protein in different dataset
# paste formuula with variables correspond to each dataset
# Conduct linear regression with lm() function for each protein
# Extract coefficients from the summary table of linear regression for each protein
# Store all parameters in the empty table created above
 
table <- data.frame(table) %>%
# Convert the table into a data frame

mutate(bonf = p.adjust(pval, method = "bonf"))
# Create a new variable to store Bonferroni-corrected p-value

table <- table[order(table$bonf), ]
# Arrange the rows based on adjusted p-values

print(head(table, 20))
# Print the top associations

write_csv(table, paste0(file,".csv"))}
# Save the results in .csv file

lm_pro(data = data_pro, protein = protein, mode = "full", file = "PWAS_Full")
lm_pro(data = data_pro_short, protein = protein, mode = "short", file = "PWAS_Short")
lm_pro(data = data_pro_long, protein = protein, mode = "long",  file = "PWAS_Long")
# Execute lm_pro function for full, short and long pneumonia ITE model

Full  <- read.csv("PWAS_Full.csv")
Short <- read.csv("PWAS_Short.csv")
Long  <- read.csv("PWAS_Long.csv")
# Read the summary statistics generated from lm_pro()

#### Volcano Plot ####
Volcano <- function(data, filename) {
# Create a function for volcano plot

data <- data %>%
mutate(effect = ifelse(bonf >= 0.05, "NS", 
                ifelse(bonf < 0.05 & beta > 0, "Up", "Down")),
       neg.log.p = -log10(pval), 
       neg.log.bonf = -log10(bonf))
# Create a new column to classify protein association
# Adjusted p-value >= 0.05, non-significant (NS)
# Adjusted p-value >= 0.05 & beta coefficient > 0, Up-regulation (Up)
# Adjusted p-value >= 0.05 & beta coefficient < 0, Down-regulation (Down)
# Calculate -log10 tra# Statistically significantnsformed raw p-value
# Calculate -log10 transformed adjusted p-value

ggplot(data, aes(x = beta, y = neg.log.p, color = factor(effect))) + 
      geom_point() + scale_color_manual(name = "Effect", values = c("Up" = "red", "NS" = "grey", "Down" = "blue")) +
      ggrepel::geom_text_repel(data = data[data$bonf < 0.05, ], 
                               aes(label = Protein), max.overlaps = 20,force = 2) +
      geom_hline(yintercept = -log10(0.05/2919) , linetype = "dashed")  + 
      geom_vline(xintercept = 0, linetype = "dashed") + 
      labs(x = "beta", y = "-log10(P)") + 
      theme_classic() + theme(plot.title = element_text(hjust = 0.5),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), 
                              panel.background = element_blank()) + 
                              scale_x_continuous(limits = c(-0.32, 0.32)) 
ggsave(paste0(filename, ".png"))}
# Generate a scater plot through ggplot
# Set the corresponding colors for effect group
# Add annotations of proteins labels for significant associations (bonf < 0.05)
# Add horizontal line for -log10 transformed significance threshold (-log10(0.05/2919))
# Add vertical line for beta coefficient = 0
# Add labels
# Adjust plot theme and x axis limits
# Save the plot 

Volcano(Full, filename = "Vol_plot_Full")
Volcano(Short, filename = "Vol_plot_Short")
Volcano(Long, filename = "Vol_plot_Long")
# Plot volcano plot for Full, Short, and Long pneumonia ITE model

#### GSEA Analysis ####
GSEA <- function(data, msig_collection, subcollection, database, range) {
# Create a function for GSEA analysis

set.seed(42)
# Set the random number starting point to ensure data reproducibility

data_GSEA <- data %>%
  mutate(rank = sign(beta)*-log10(pval)) %>%
  arrange(desc(rank))
# Create a new variable rank as the product of -log10 transformed p-value and the sign of beta coefficient
# Arrange rank variable in a descending order

rankings <- data_GSEA$rank
# Assign the values in rank to a new vector rankings

names(rankings) <- data_GSEA$Protein
# Assign the protein name in data_GSEA to rankings

print(head(rankings))
# Check the top rankings

list <- msigdbr(species = "Homo sapiens", collection = paste0(msig_collection), subcollection = paste0(subcollection))
# Get msigdbr database of designated species, collection, and subcollection

list_split <- split(x = list$gene_symbol, f = list$gs_name)
# Create a list of the genes splitted by gene set

print(head(list_split))
# Print the first few gene set

if (database == "KEGG") {
gsea <- fgseaMultilevel(pathway = list_split,
                        stats = rankings,
                        minSize = 15,
                        maxSize = 500,
			nPermSimple = 100000) %>%
	arrange(padj) %>%
        filter(!is.na(pval)) %>%
        mutate(pathway = str_remove_all(pathway, "KEGG"),
               pathway = str_replace_all(pathway, "_", " "),
               pathway = str_to_title(pathway),
	       pathway = str_trim(pathway))

# Create an if loop for KEGG database
# Conduct GSEa analysis with fgseaMultilevel()
# Gene set list is the input of pathway
# Rankings from PWAS is the input of stats
# minSize was set as 15 to exclude gene set with limited genes
# maxSize was set as 500 to balance computational efficiency, and biological relevancy
# nPermSimple was set as 100000 to increase randomness
# Arrange the table based on ajusted p-value
# Filter missing values
# Change gene set names

print(dim(gsea))
print(head(gsea))
# Check the output of fgseaMultilevel

plot_1 <- ggplot(gsea[1:range, ], aes(reorder(pathway, -log10(padj)), -log10(padj))) +
	  geom_col(aes(fill = padj < 0.05)) + 
	  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + 
	  coord_flip() + theme_classic() + 
	  theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
		legend.position = "bottom") +
	  scale_fill_manual(values = c(`TRUE` = "orangered", `FALSE` = "lightsalmon")) +
	  labs(x = "Pathway", y = "-log10(padj)",
               title="KEGG pathways from GSEA")
# Create a vertical bar chart through ggplot
# Pathways were reordered based on -log10 transformed adjusted p-value
# Classify the pathway by the threshhold of 0.05 for adjusted p-value
# Add a horizontal line that marks the threshold (-log10(0.05))
# Flip the plot vertically
# Adjust plot theme and legend position
# Assign color for the pathways
# Adjust labels and title

write.csv(gsea[, 1:7], file = paste0("GSEA_", deparse(substitute(data)), "_", database, ".csv"))
ggsave(paste0("GSEA_", deparse(substitute(data)), "_", database, ".png"), plot = plot_1)}
# Save the table and plot individually

else if (database == "GO") {
gsea <- fgseaMultilevel(pathway = list_split,
                        stats = rankings,
                        minSize = 15,
                        maxSize = 500,
	                nPermSimple = 100000) %>%
        arrange(padj) %>%
        filter(!is.na(pval)) %>%
	mutate(pathway = str_remove_all(pathway, "GOBP"),
               pathway = str_replace_all(pathway, "_", " "),
               pathway = str_to_title(pathway),
	       pathway = str_trim(pathway))
# Same parameters as above but different gene set database

print("Adaptive Immune Response Based On Somatic Recombination Of Immune Receptors Built From Immunoglobulin Superfamily Domains" %in% gsea$pathway)
# Check the existence of a specific pathway before modify its name 

gsea[gsea$pathway == "Adaptive Immune Response Based On Somatic Recombination Of Immune Receptors Built From Immunoglobulin Superfamily Domains", "pathway"] <-
                     "Adaptive Immunity by Somatic Recombination of Ig Domains"
# Modify pathway name to make it more concise 

print(dim(gsea))
print(head(gsea))
# Check the output of fgseaMultilevel

plot_2 <- ggplot(gsea[1:range, ], aes(reorder(pathway,-log(padj)), -log(padj))) +
	  geom_col(aes(fill = padj < 0.05)) + 
	  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + 
	  coord_flip() + theme_classic() + 
	  theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
		legend.position = "bottom") +
	  scale_fill_manual(values = c(`TRUE` = "orangered", `FALSE` = "lightsalmon")) + 
          labs(x = "Pathway", y = "-log10(padj)",
               title = "GO pathways from GSEA")
# Same ggplot layout
# Adjust plot title
 
write.csv(gsea[, 1:7], file = paste0("GSEA_", deparse(substitute(data)), "_", database, ".csv"))
ggsave(filename = paste0("GSEA_", deparse(substitute(data)), "_", database, ".png"), plot = plot_2)}}
# Save the table and plot individually

GSEA(Full, "C2", "CP:KEGG_LEGACY", "KEGG", 20)
GSEA(Short, "C2", "CP:KEGG_LEGACY", "KEGG", 20)
GSEA(Full, "C5", "GO:BP", "GO", 20)
GSEA(Short, "C5", "GO:BP", "GO", 20)
## Execute GSEA function to generate GSEA tables and plots for full and short pneumonia ITE model

#### Venn diagram ####
gsea_full <- read.csv("GSEA_Full_GO.csv")
# Read gsea_full table

gsea_short <- read.csv("GSEA_Short_GO.csv")
# Read gsea_short table

venn_full <- unlist(gsea_full[gsea_full$padj < 0.05,"pathway"])
# Creat a list of all significant pathways in gsea_full

venn_short <- unlist(gsea_short[gsea_short$padj < 0.05, "pathway"])
# Creat a list of all significant pathways in gsea_short

diff_short <- list(venn_short[!venn_short %in% venn_full])
# Create a list of unique pathways in gsea_short

diff_full <- list(venn_full[!venn_full %in% venn_short])
# Create a list of unique pathways in gsea_full
 
print(diff_short)
print(diff_full)
# Check unique pathways in both gsea_full and gsea_short

venn.diagram(
  x = list(venn_full, venn_short),
  category.names = c("Full" , "Short"),
  filename = 'venn_diagramm_GSEA_GO.png',
  col = c("navyblue", "yellow"),
  fill = c(alpha("navyblue", 0.3), alpha("yellow", 0.3)),
  cat.col = c("navyblue", "yellow"),
  output= TRUE)
# plot a Venn diagram to show the number of shared and unique pathways in two models

#### Correlation plot ####
merge <- left_join(x = Full, y = Short, by = "Protein") %>%
  mutate(Significance = ifelse(bonf.x < 0.05 & bonf.y < 0.05, "Both", 
                        ifelse(bonf.x < 0.05 & bonf.y > 0.05, "Full",
                        ifelse(bonf.x > 0.05 & bonf.y < 0.05, "Short",
                        "NS")))) %>%
  arrange(Protein)
# Merge PWAS_Full (Full) and PWAS_Short(Short) dataset
# Create a new variable indicating the model associated with proteins' statistical significance
# Statistically significant in both Full and Short, "Both"
# Statistically significant in Full not Short, "Full"
# Statistically significant in Short not Full, "Short"
# Statistically insignificant in all model, "NS"
# Arrange the protein
 
cor <- rcorr(x = merge$beta.x, y = merge$beta.y, type = "pearson")
# Calculate pearson correlation coefficient between Full and Short

coeff <- cor$r[1,2]
# Extract correlation coefficient

coeff <- round(coeff, 4)
# Round to four decimals

cor_plot <- ggplot(data = merge, aes(x = beta.y, y = beta.x)) + 
            geom_point(aes(shape = Significance,
                           colour = Significance)) + 
            stat_smooth(method = lm, color = "red3") +
            labs(y = "Fully adjusted model (beta)",
                 x = "Short-term ITE model (beta)") + theme_classic() +
            scale_shape_manual(name = "Significance", values = c("NS" = 1, "Full" = 19, 
                                                                 "Short" = 19, "Both" = 19)) +
            scale_color_manual(name = "Significance", values = c("NS" = "black", 
					                         "Full" = "#55C667FF", 
                                                                 "Short" = "#FDE725FF", 
                                                                 "Both" = "#FF4500")) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
            annotate(geom = "text", x = 0.05, y = 0.25, label = paste0("R = ", coeff)) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), 
                  panel.background = element_blank())
# Generate a scatter plot through ggplot
# Plot the points' shape and color based on Significance
# Add linear regression of Full and Short to the plot
# Adjust labels and plot theme
# Set color and shape for each Significant group
# Add a reference line with a slope of 1
# Add text annotation for correlation coefficient derived from cor
# Adjust plot theme 

ggsave("cor_plot.png", plot = cor_plot)
# Save the correlation plot

################################################################################################################################################
