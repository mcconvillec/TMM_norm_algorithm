# load libraries
library("tidyverse")
library("psych")
library("gridExtra")
library("scales")
library("limma") 
library("edgeR")
library("readr")
library("dplyr")
library("limma")
library("Glimma")
library("manhattanly")

# read the grouped protein summary file
MQ_raw <- read_tsv("./proteinGroups.txt", col_names = TRUE, guess_max = 5427)
MQ_tmt <- MQ_raw %>% select("Gene names", starts_with("Reporter intensity corrected 1"), starts_with("Reporter intensity corrected 5"), starts_with("Reporter intensity corrected 6"), starts_with("Reporter intensity corrected 7"))
MQ_tmt <- MQ_tmt %>% select(-"Reporter intensity corrected 1",-"Reporter intensity corrected 5", -"Reporter intensity corrected 6", -"Reporter intensity corrected 7", -starts_with("Reporter intensity corrected 10"), -starts_with("Reporter intensity corrected 11"))

#separate gene names from data
accession <- as.data.frame(MQ_tmt$`Gene names`)
MQ_tmt <- MQ_tmt %>% select(-"Gene names")

#rename columns and create groups/conditions
colnames(MQ_tmt) <- c("eHap_1", "eHap_2", "eHap_3",
                      "COMMD6_1", "COMMD6_2", "COMMD6_3",
                      "COMMD9_1", "COMMD9_2", "COMMD9_3",
                      "COMMD1+6_1", "COMMD1+6_2", "COMMD1+6_3")

#Load data into DGEList object
group <- c(rep("wildtype", 3), rep("COMMD6KO", 3), rep("COMMD9KO", 3), rep("COMMD1+6KO", 3))
y <- DGEList(counts=MQ_tmt, group=group, genes=accession)


#RUN TMM NORMALISATION

#---set colors for all future boxplots
color <- c(rep("red", 3), rep("blue", 3))

#---compute normalised intensities 
apply_tmm_factors <- function(y, color = NULL, plot = TRUE) {
  
  #compute grand total (library size scalings)
  #i.e. what is the ratio between the lib.size of each sample and the mean lib.size
  lib_facs <- mean(y$samples$lib.size) / y$samples$lib.size
  
  #TMM factors are library adjustment factors (so divide by them)
  norm_facs <- lib_facs / y$samples$norm.factors
  cat("Overall Factors (lib.size+TMM): \n", sprintf("%-5s -> %f\n",colnames(y$counts), norm_facs))
  
  #compute the normalized data as a new data frame
  tmt_tmm <- as.data.frame(sweep(y$counts, 2, norm_facs, FUN="*"))
  colnames(tmt_tmm) <- str_c(colnames(y$counts), "_tmm")
  
  #
  y.name = deparse(substitute(y))
  #visualise results and return data frame
  if(plot == TRUE){
    boxplot(log2(tmt_tmm), col=color,notch=TRUE, main=sprintf("TMM Normalised data %s", y.name))
  }
  tmt_tmm #returns the dataframe
}

#TMM normalise data
MQ_tmt_tmm <- apply_tmm_factors(y, color)
MQ_tmt_tmm <- bind_cols(MQ_tmt_tmm, accession)
write.table(MQ_tmt_tmm, row.names = FALSE, 
            file = "./11plex_proteinGroups_TMM.txt", sep="\t", eol="\r\n")

