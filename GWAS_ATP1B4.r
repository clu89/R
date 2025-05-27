 

## **Load library and dataset** ##
install.packages("pacman")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.20")
library(pacman)
p_load(tidyverse, GENESIS, SeqArray, SeqVarTools, Biobase, SPAtest, qqman)

## Set Working Directory
setwd("C:/Users/caspe/OneDrive/Desktop/Casper/Emory University/Global Health/Fall 2/EPI-556 Applied Genome Epidemiology/Final Project")

## Load dataset
data <- SeqArray::seqOpen("GWAS.gds")
data

## Create annotation dataset
variant.id <- SeqArray::seqGetData(data, "variant.id")
sample.id <- SeqArray::seqGetData(data, "sample.id")
rs <- SeqArray::seqGetData(data, "annotation/id")
chr <- SeqArray::seqGetData(data, "chromosome")
pos <- SeqArray::seqGetData(data, "position")
allele <- SeqArray::seqGetData(data, "allele")
alt <- SeqArray::seqGetData(data, "$alt")
annot <- data.frame(variant.id=variant.id,rs=rs,chromosome=chr,
                  position=pos,allele=allele,alt=alt,stringsAsFactors=F)

## Genotype dataset
geno <- SeqArray::seqGetData(data, "genotype")

## **Quality Control** ##
## Missing Rate < 0.05
## HWE > 1*10^(-4)
## MAF > 0.1
## Missing Rate Filter
Mis_Var <- SeqVarTools::missingGenotypeRate(data, margin = "by.variant")
Mis_Var <- Mis_Var > 0.05
length(Mis_Var) ## Total 434985
Mis_Ex <- variant.id[Mis_Var]
length(Mis_Ex) ## Exclude 108501

## HWE Filter
HWE <- hwe(data)
dim(HWE) ## Total 434985 
HWE_low <- !is.na(HWE$p) & HWE$p<1e-4
HWE_Ex <- variant.id[HWE_low]
length(HWE_Ex) ## Exclude 503

## Minor allele Frequency Filter
MAF <- seqAlleleFreq(data, minor = T)
length(MAF) ## Total 434985
rare <- !is.na(MAF) & MAF < 0.1
MAF_Ex <- variant.id[rare]
length(MAF_Ex) ## Exclude 148044

## Combine Filter and SNP selection
Ex_Gen <- c(Mis_Ex, HWE_Ex, MAF_Ex)
keepSNP <- setdiff(variant.id, Ex_Gen)
length(keepSNP) ## Total exclude variant 213650

## Individual Call rate filter for samples
Mis_Sam <- missingGenotypeRate(data, margin = "by.sample")
Sam_Mis <- Mis_Sam > 0.05
Sam_Ex <- sample.id[Sam_Mis]
length(Sam_Ex) ## Remain 28

## Samples selection
keepSample <- setdiff(sample.id, Sam_Ex)
length(keepSample) ## Remain 150

## **Incorporate phenotype dataset** ##
Phe <- read.csv("GENESIS_final_pheno.csv")
dim(Phe)
names(Phe)
head(Phe)
Phe <- as(Phe, "AnnotatedDataFrame")
Biobase::varMetadata(Phe)
table(is.na(Phe$ATP1B4_1)) ## No missing ATP1B4_1
table(is.na(Phe$inferred_population)) ## No missing inferred race data
Cauc_ID <- Biobase::pData(Phe)[pData(Phe)$inferred_population =="Cauc" & 
                    !is.na(pData(Phe)$inferred_population),]$sample.id
Trig_ID <- Biobase::pData(Phe)[!is.na(pData(Phe)$Trig),]$sample.id
keepSample <- intersect(keepSample, Cauc_ID)
keepSample <- intersect(keepSample, Trig_ID)
table(Phe$inferred_population) ## 178 Cauc
summary(Phe$Trig)

## Single-variant analysis_1: ATP1B4_1 ~ SNP
SeqArray::seqSetFilter(data, variant.id = keepSNP, sample.id = keepSample)
seqData <- SeqVarTools::SeqVarData(data, sampleData = Phe)
iterators <- SeqVarTools::SeqVarBlockIterator(seqData)
nullmod1 <- GENESIS::fitNullModel(iterators,outcome = "ATP1B4_1",covars = NULL ,family = "gaussian")
assoc_1 <- GENESIS::assocTestSingle(iterators, nullmod1, test = "Score")
summary(assoc_1$Score.pval)
assoc_1 <- merge(annot, assoc_1, by = "variant.id", all.y = T)
assoc_1 <- assoc_1[order(assoc_1$Score.pval),]

## Single-variant analysis_2: ATP1B4_1 ~ SNP + Age + sex
SeqArray::seqSetFilter(data, variant.id = keepSNP, sample.id = keepSample)
seqData <- SeqVarTools::SeqVarData(data, sampleData = Phe)
iterators <- SeqVarTools::SeqVarBlockIterator(seqData)
nullmod2 <- GENESIS::fitNullModel(iterators,outcome = "ATP1B4_1",covars = c("Age", "sex") ,family = "gaussian")
assoc_2 <- GENESIS::assocTestSingle(iterators, nullmod2, test = "Score")
summary(assoc_2$Score.pval)
assoc_2 <- merge(annot, assoc_2, by = "variant.id", all.y = T)
assoc_2 <- assoc_2[order(assoc_2$Score.pval),]

## QQ plot
qq(assoc_1$Score.pval)
qq(assoc_2$Score.pval)

## Manhattan plot
assoc_1$chromosome <- as.numeric(assoc_1$chromosome)
manhattan(assoc_1,chr="chromosome",bp="pos.x",p="Score.pval",snp="rs.x",suggestiveline=-log10(1e-2),genomewideline=-log10(5e-08))

assoc_2$chromosome <- as.numeric(assoc_2$chromosome)
manhattan(assoc_2,chr="chromosome",bp="pos.x",p="Score.pval",snp="rs.x",suggestiveline=-log10(1e-2),genomewideline=-log10(5e-08))
