
##set work directory
getwd()
setwd("C:\\Users\\caspe\\OneDrive\\Desktop\\Casper\\Emory University\\Global Health\\APE_REAL\\Singh Lab\\Projects\\Phenylalanine Level Monitor")

##Package installation
install.packages("pacman")
install.packages("readxl")

library(pacman)
p_load(rio, tidyverse, janitor, table1, BlandAltmanLeh, ggplot2, readxl, 
       meantables, summarytools, rempsyc, flextable, viridis, epiR, SimplyAgree, 
       cowplot, mcr, Cairo, flextable, webshot, htmlwidgets, DescTools)

##Import Camp phe monitor spreadsheet and calculate mean and difference
Camp_Phe <- read_excel('C:\\Users\\caspe\\OneDrive\\Desktop\\Casper\\Emory University\\Global Health\\APE_REAL\\Singh Lab\\Projects\\Phenylalanine Level Monitor\\Camp Dataset\\Modified dataset\\Camp dataset_De-identify.xlsx') %>%
            mutate(PAA_Phe = as.numeric(PAA_Phe), 
                   DBS_Phe = as.numeric(DBS_Phe),
                   avg_Phe = (PAA_Phe+DBS_Phe)/2,
                   Phe_dif = (PAA_Phe - DBS_Phe),
                   Phe_dif_prop = (Phe_dif/avg_Phe)*100,
                   Treatment = replace(Treatment, 
                   Treatment == "Palynziq", "Pegvaliase"),
                   Treatment = replace(Treatment, 
                   Treatment == "Kuvan", "Sapropterin Dihydrochloride"),
                   Treatment = replace_na(Treatment, "Diet only"))

table(Camp_Phe$Treatment)
##Examination of data normality with Shapiro-wilk test
s1 <- shapiro.test(Camp_Phe$PAA_Phe)
s2 <- shapiro.test(Camp_Phe$DBS_Phe)
s1
s2

##Summarize percent mean difference, upper and lower limits
sum_Phe <- summarize(Camp_Phe, `mean difference (%)` = mean(Camp_Phe$Phe_dif_prop),
                   `lower limit` = `mean difference (%)`-1.96*sd(Phe_dif_prop),
                   `upper limit` = `mean difference (%)`+1.96*sd(Phe_dif_prop))

##Combining vectors to form new dataframe for annotating the ggplot for Bland-Altman
text_data <- data.frame(label = c("lower", "mean", "upper"),
                        y_position = c(sum_Phe$`lower limit`,
                                       sum_Phe$`mean difference (%)`,
                                       sum_Phe$`upper limit`))

##Create new dataset for slope chart
Camp_Phe_table <- Camp_Phe %>%
  pivot_longer(cols = c(DBS_Phe, PAA_Phe), 
               names_to = "Sample", 
               values_to = "Phe_level_blood")


##Descriptive statistical table with Wilcoxon signed rank test 
Phe_t <- t.test(Camp_Phe$DBS_Phe, Camp_Phe$PAA_Phe, paired = TRUE)
Phe_t

Phe_wil <- wilcox.test(Camp_Phe$DBS_Phe, Camp_Phe$PAA_Phe, paired = TRUE)

Phe_descr <- as.data.frame(t(descr(Camp_Phe[5:6]))) %>%
             rownames_to_column("Sample") %>%
             select('Sample', 'N.Valid', 'Mean', 'Std.Dev', 'Median', 'Min', 'Max') %>%
             mutate(`p` = Phe_wil$p.value) %>%
             rename('N' = N.Valid, 'SD' = Std.Dev)

Phe_descr$Sample[Phe_descr$Sample == c("DBS_Phe","PAA_Phe")] <- 
                                     c("Dried Blood Spots (DBS)", 
                                       "Plasma Amino Acid (PAA)")            

Phe_table <- nice_table(Phe_descr, title =  c("Table 1.", "Baseline characteristics of Phe level in plasma and DBS"),
             note = c("Statistical significance was calculated with Wilcoxon signed-rank test.",
                  "* p < .05, ** p < .01, *** p < .001",
                  "Phe concentration was noted as (uM)"))
print(Phe_table)

save_as_docx(Phe_table, path = "Phe_table.docx")



##Generate bland altman plot
Plot1 <- bland.altman.plot(Camp_Phe$PAA_Phe, 
                           Camp_Phe$DBS_Phe, mode = 1, 
                           graph.sys = "ggplot2") + 
         labs(x= "Phe average of plasma and DBS (uM)", 
              y = "Plasma Phe - DBS Phe (uM)")

Plot1 + aes(ymax = 800, ymin = -800) + 
        geom_text(data = data.frame(y = Plot1$layers[[2]]$data$yintercept) %>%
                  rownames_to_column(), 
        aes(x = 1000, y = y, label = paste0(rowname, " ", round(y,2))),
                  hjust = 0.4, vjust = -0.8) + 
        ggtitle("Bland-Altman plot for Phe detection in PAA and DBS") +
        theme_bw() + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(),
              plot.title = element_text(hjust = 0.5, vjust = 0.8, face = "bold"))

## plot size 700*420
                   
##Bland Altman plot with percent mean difference
Plot2 <- ggplot(data = Camp_Phe, 
         aes(x= avg_Phe, y = Phe_dif_prop, ymax = 200, ymin =-200)) + 
         geom_point(aes(color = Treatment)) + 
         geom_hline(yintercept = text_data$y_position, linetype = c(2,1,2)) + 
         labs(x= "Phe average of plasma amino acid and dry blood spots (μM)", y = "PAA Phe - DBS Phe (%)")

Plot2 + geom_text(data = text_data, aes(x = 1000, y = y_position, 
        label = paste0(label,":", round(y_position,2)), vjust = -0.8, hjust = 0.3)) +
        ggtitle("Bland-Altman plot for Phe detection in PAA and DBS") +
        theme_bw() + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(),
              plot.title = element_text(hjust = 0.5, vjust = 0.8, face = "bold")) +
        theme(legend.position = c(0.5, 0.01),
              legend.justification = c(0.5, 0),
              legend.direction = "horizontal")

## plot size 700*420

##Slope chart
Camp_Phe_table_1 <- Camp_Phe_table

Camp_Phe_table_1$Sample[Camp_Phe_table_1$Sample == c("DBS_Phe","PAA_Phe")] <- 
                                                   c("Dried Blood Spots (DBS)", 
                                                   "Plasma Amino Acid (PAA)")
Phe_slope <- ggplot(data = Camp_Phe_table_1,
             aes(x = Sample,
             y = Phe_level_blood,
             group = ID)) +
             geom_line(linewidth = 0.5, alpha = 1) +
             geom_point(size = 1.5) + 
             labs(y = "Phe level (uM)") + 
             ggtitle("Slope plot for Phe detection in PAA and DBS") + 
             theme(plot.title = element_text(hjust = 0.5, vjust = 0.8, face = "bold"))+
             theme_bw() + 
             theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), 
                  panel.background = element_blank(), 
                  axis.title.y = element_blank())
  
print(Phe_slope)

##Regression and Correlation
##Concordance Correlation Coefficient
Phe_ccc <- epi.ccc(Camp_Phe$DBS_Phe, Camp_Phe$PAA_Phe, ci = "z-transform", 
                   conf.level = 0.95)
Phe_ccc

##Simple Linear Regression (SLR)
Phe_lm = lm(PAA_Phe ~ DBS_Phe, data = Camp_Phe)
Phe_SLR <- ggplot(data = Camp_Phe, aes(x = DBS_Phe, y = PAA_Phe)) + 
           geom_point() + 
           geom_smooth(method = "lm", se = TRUE, lwd = 0.5, colour = "black") + 
           labs(x = "Phe in dry blood spots (uM)", y = "Phe in plasma amino acid (uM)") + 
           annotate("text", label = "y = 1.39x + 37.46", x = 900, y = 300) + 
           ggtitle("Simple linear regression for PAA and DBS") + 
           theme(plot.title = element_text(hjust = 0.5, vjust = 0.8, face = "bold"))

data_lm <- data.frame(Phe_lm$coefficients)

Phe_lm
print(Phe_SLR)

##Deming Regression
Phe_DR <- dem_reg(x = "DBS_Phe", y = "PAA_Phe", data = Camp_Phe, weighted = FALSE)

Phe_DR

check(Phe_DR)

data_DR <- data.frame(Phe_DR$model)


Phe_DR_plot <- ggplot(data = Camp_Phe, aes(x = DBS_Phe, y = PAA_Phe)) + 
               geom_point() + 
               geom_abline(slope = c(data_DR[2,1], 1, data_lm[2,1]), 
               intercept = c(data_DR[1,1], 0, data_lm[1,1]), 
               linetype = c(1, 2, 1), color = c("steelblue", "black", "navy")) + 
               labs(x = "Phe in dry blood spots (uM)", 
               y = "Phe in plasma amino acid (uM)") + 
               ggtitle("Deming regression for PAA and DBS") + 
               theme(plot.title = element_text(hjust = 0.5, vjust = 0.8, face = "bold")) + 
               annotate("text", label = "y = 1.56x - 32.74", x = 900, y = 300)
                
Phe_DR_plot 

##Passing-Bablock Regression
Phe_PBR <- mcreg(Camp_Phe$DBS_Phe, Camp_Phe$PAA_Phe, method.reg = "PaBa")

summary(Phe_PBR)

data_PBR <- data.frame(Phe_PBR@para)

Phe_PBR_plot <- ggplot(data = Camp_Phe, aes(x = DBS_Phe, y = PAA_Phe)) + 
                geom_point() + 
                geom_abline(slope = c(data_PBR[2,1], 1, data_lm[2,1]), 
                intercept = c(data_PBR[1,1], 0, data_lm[1,1]), 
                linetype = c(1, 2, 1), 
                color = c("springgreen", "black", "navy")) + 
                labs(x = "Phe in dry blood spots (uM)", y = "Phe in plasma amino acid (uM)") + 
                ggtitle("Passing-Bablok regression for PAA and DBS") + 
                theme(plot.title = element_text(hjust = 0.5, vjust = 0.8, face = "bold")) + 
                annotate("text", label = "y = 1.48x - 10.98", x = 900, y = 300)

print(Phe_PBR_plot)

##Summrized plot
Model <- c("Reference line", "Simple linear regression", 
           "Deming regression", "Passing-Bablok regression")
Slope <- c(1, data_lm[2,1], data_DR[2,1], data_PBR[2,1])
Intercept <- c(0, data_lm[1,1], data_DR[1,1], data_PBR[1,1])
data_sum <- data.frame(t(rbind(Model, Slope, Intercept)))
line <- as.factor(c(2, 1, 1, 1))
  
Phe_Sum_plot <- ggplot(data = Camp_Phe, aes(x = DBS_Phe, y = PAA_Phe)) + 
  geom_point() + 
  geom_abline(data = data_sum,
              aes(slope = as.numeric(Slope), 
                  intercept = as.numeric(Intercept), 
                  linetype = Model,
                  colour = Model)) + 
  scale_color_manual(name = "Model", 
                     values = c("Reference line" = "black",
                               "Simple linear regression" = "turquoise",
                               "Unweighted Deming regression" = "violet",
                               "Passing-Bablok regression" = "orangered")) +
  scale_linetype_manual(name = "",
                        guide = "none",
                        values = c("Reference line" = "dashed",
                                   "Simple linear regression" = "solid",
                                   "Unweighted Deming regression" = "solid",
                                   "Passing-Bablok regression" = "solid")) +
  labs(x = "Phe in Dry Blood Spots (uM)", 
       y = "Phe in Plasma Amino Acid (uM)") + 
  ggtitle("Linear regression for PAA and DBS") + 
  theme(plot.title = element_text(hjust = 0.5, 
                                  vjust = 0.8,
                                  face = "bold")) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position = c(0.8, 0.2),
        plot.title = element_text(hjust = 0.5, vjust = 0.8, face = "bold"))

print(Phe_Sum_plot)
##Plot size 700*420

##Model Examination with ccc correlation
Camp_Phe_model <- Camp_Phe %>%
  mutate(DBS_Dem = 1.56*DBS_Phe - 32.74,
         DBS_PBR = 1.476*DBS_Phe - 10.983,
         DBS_SLR = 1.388*DBS_Phe + 36.461)

Dem_cal <- epi.ccc(Camp_Phe_model$DBS_Dem, 
           Camp_Phe_model$PAA_Phe, 
           ci = "z-transform", 
           conf.level = 0.95)

PBR_cal <- epi.ccc(Camp_Phe_model$DBS_PBR, 
           Camp_Phe_model$PAA_Phe, 
           ci = "z-transform", 
           conf.level = 0.95)

SLR_cal <- epi.ccc(Camp_Phe_model$DBS_SLR, 
           Camp_Phe_model$PAA_Phe, 
           ci = "z-transform", 
           conf.level = 0.95)

##Conclusion table
Phe_Sum <- data.frame(rbind(Phe_ccc$rho.c,
                            SLR_cal$rho.c,
                            Dem_cal$rho.c,
                            PBR_cal$rho.c))

Phe_Sum$Model <- c("Raw data",
                   "Simple Linear",
                   "Deming Unweighted",
                   "Passing Bablok") 

Phe_Sum <- Phe_Sum[, c(4,1,2,3)]

colnames(Phe_Sum) <- c("Model",
                       "Estimate",
                       "Lower C.I",
                       "Upper C.I")

round_4 <- function(x)
                   {paste(round(x, 4))}

Sum_table <- nice_table(Phe_Sum,
                        title =  "Table 2. Agreement rates of Phe level in DBS and PAA with model adjustment",
                        note = c("Agreement rates was calculated with Concordance Correlation Coefficient.",
                                 "Coefficients were transformed with Fisher's Z transformation.",
                                 "C.I indicates 95% confidence interval for Z distribution."),
                        col.format.custom = 2:4,
                        format.custom = "round_4")
Sum_table

save_as_docx(Sum_table, path = "Phe_table_2.docx")
