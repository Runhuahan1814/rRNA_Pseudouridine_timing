# libraries ---- 
library(here)
library(vroom)
library(tidyverse)
library(ape)
library(GenomicAlignments)

# set folder as current
setwd("/input_directory")

# calculate read identity
read_in_bam <- function(file){
  init <- readGAlignments(file, use.names = T, param = ScanBamParam(tag=c("NM"), what="mapq"))
  init_t <- GenomicAlignments::as.data.frame(init) %>%
    dplyr::mutate(minion_read_name = names(init),
                  mapped_gene = seqnames) 
  init_t_final <- init_t %>%
    dplyr::mutate(identity = (1 - NM/aligned_reads)*100) %>%
    dplyr::group_by(minion_read_name) %>%
    dplyr::filter(identity == max(identity),
                  aligned_reads == max(aligned_reads)) %>%
    dplyr::distinct(minion_read_name, .keep_all = T) %>%
    dplyr::mutate(gene = str_split_fixed(mapped_gene,"-",2)[,2])
  return(init_t_final)
}

# Read in desired .bam file
pre_rRNA <- read_in_bam("./input_directory/input.sorted.bam")

# Filter bam file for reads with identity of 80% or higher
pre_rRNA_good <- pre_rRNA %>%
  dplyr::filter(identity >= 80) 

# Extract read name
dput(pre_rRNA_good$minion_read_name, "pre_rRNA_good.txt")

# Classify the reads in the bam file with the processing stage, then group into processing stages
# The first several mutate commands add a column for a processing stage (i.e. 20S, 27SA2)
# The if_else command places a 1 (yes) if the read has the desired qualities, and a 0 (no) if it does not
# For example, 20S is a yes if the read starts at position 701 or more and ends between 2501 and 2712
# These qualities can be adjusted or enhanced to fine tune the clustering
# The final mutate command adds a new column called group, and assigns the group based on the values of the processing states
pre_rRNA_good_sort <- pre_rRNA_good  %>%
  dplyr::mutate(twentyS = if_else(start >= 701 && between(end, 2501, 2712), 1, 0)) %>%
  dplyr::mutate(twentysevenSAtwo = if_else(between(start, 2713, 2784) && end <= 6661, 1, 0)) %>%
  dplyr::mutate(twentysevenSAthree = if_else(between(start, 2785, 2856) && end <= 6647, 1, 0)) %>%
  dplyr::mutate(twentysevenSBS = if_else(start >=2862 && end <= 6647, 1, 0)) %>%
  dplyr::mutate(twentysevenSBL = if_else(between(start, 2856, 2861) && end <= 6647, 1, 0)) %>%
  dplyr::mutate(sevenSS = if_else(start >=2862 && end <= 3158, 1, 0)) %>%
  dplyr::mutate(sevenSL = if_else(between(start, 2856, 2861) && end <= 3158, 1, 0)) %>%
  dplyr::mutate(twentyfivepointfiveS = if_else(between(start, 3159, 3252) && end <= 6647, 1, 0)) %>%
  dplyr::mutate(fiftyeightSL = if_else(between(start, 2856, 2861) && end <= 3158, 1, 0)) %>%
  dplyr::mutate(eighteenS = if_else(start >= 701 && end <= 2500, 1, 0)) %>%
  dplyr::mutate(twentyfiveS = if_else(start >= 3252 && end <= 6647, 1, 0)) %>%
  dplyr::mutate(fiftyeightS = if_else(start >= 2862 && end <= 3019, 1, 0)) %>%
  dplyr::mutate(thirtyfiveS = if_else(between(start, 0, 609) | end >= 6662, 1, 0)) %>%
  dplyr::mutate(thirtythreeS = if_else(between(start, 609, 700) | between(start, 609, 700) && between(end, 6647, 6662) | between(start, 609, 700) && end >=2713, 1, 0)) %>%
  dplyr::mutate(thirtytwoS = if_else(between(start, 701, 2712) && end >=2713 | between(start, 701, 2712) && between(end, 6647, 6662), 1, 0)) %>%
  dplyr::mutate(group = case_when((twentyS == 1 | twentysevenSAtwo == 1) ~ "intermediate",
                                  (eighteenS == 1 | twentyfiveS == 1 | fiftyeightS == 1) ~ "mature",
                                  (twentysevenSAthree == 1 |
                                     twentysevenSBS == 1 |
                                     twentysevenSBL == 1 | 
                                     sevenSS == 1 |
                                     sevenSL == 1 |
                                     twentyfivepointfiveS == 1 |
                                     fiftyeightSL == 1) ~ "late",
                                  (thirtyfiveS == 1) ~ "primary",
                                  (thirtytwoS == 1 | thirtythreeS == 1) ~ "early"))

# Break classified reads into individual tables for the different processing stages                                  
pre_rRNA_good_primary <- pre_rRNA_good_sort %>%
  dplyr::filter(group == "primary")
pre_rRNA_good_early <- pre_rRNA_good_sort %>%
  dplyr::filter(group == "early")
pre_rRNA_good_intermediate <- pre_rRNA_good_sort %>%
  dplyr::filter(group == "intermediate")
pre_rRNA_good_late <- mpre_rRNA_good_sort %>%
  dplyr::filter(group == "late")
pre_rRNA_good_mature <- pre_rRNA_good_sort %>%
  dplyr::filter(group == "mature")

# Output the list of minion read names from the different processing stages.
# This file can then be formatted and used with samtools to output a bam file containing only these reads.
dput(pre_rRNA_good_primary$minion_read_name, "pre_rRNA_good_primary.txt")
dput(pre_rRNA_good_early$minion_read_name, "pre_rRNA_good_early.txt")
dput(pre_rRNA_good_intermediate$minion_read_name, "pre_rRNA_good_intermediate.txt")
dput(pre_rRNA_good_late$minion_read_name, "pre_rRNA_good_late.txt")
dput(pre_rRNA_good_mature$minion_read_name, "pre_rRNA_good_mature.txt")

#classification of the rRNA reads into every single intermediate in the rRNA maturation pathway
pre_rRNA_good_sort_each <- pre_rRNA_good %>%
  dplyr::mutate(twentyS = if_else(start >= 701 && between(end, 2501, 2712), 1, 0)) %>%
  dplyr::mutate(twentysevenSAtwo = if_else(between(start, 2713, 2784) && end <= 6661, 1, 0)) %>%
  dplyr::mutate(twentysevenSAthree = if_else(between(start, 2785, 2856) && end <= 6647, 1, 0)) %>%
  dplyr::mutate(twentysevenSBS = if_else(start >=2862 && end <= 6647, 1, 0)) %>%
  dplyr::mutate(twentysevenSBL = if_else(between(start, 2856, 2861) && end <= 6647, 1, 0)) %>%
  dplyr::mutate(sevenSS = if_else(start >=2862 && end <= 3158, 1, 0)) %>%
  dplyr::mutate(sevenSL = if_else(between(start, 2856, 2861) && end <= 3158, 1, 0)) %>%
  dplyr::mutate(twentysixS = if_else(between(start, 3159, 3252) && end <= 6647, 1, 0)) %>%
  dplyr::mutate(fiftyeightSL = if_else(between(start, 2856, 2861) && end <= 3158, 1, 0)) %>%
  dplyr::mutate(eighteenS = if_else(start >= 701 && end <= 2500, 1, 0)) %>%
  dplyr::mutate(twentyfiveS = if_else(start >= 3252 && end <= 6647, 1, 0)) %>%
  dplyr::mutate(fiftyeightS = if_else(start >= 2862 && end <= 3019, 1, 0)) %>%
  dplyr::mutate(thirtyfiveS = if_else(between(start, 0, 609) | end >= 6662, 1, 0)) %>%
  dplyr::mutate(thirtythreeS = if_else(between(start, 609, 700) | between(start, 609, 700) && between(end, 6647, 6662) | between(start, 609, 700) && end >=2713, 1, 0)) %>%
  dplyr::mutate(thirtytwoS = if_else(between(start, 701, 2712) && end >=2713 | between(start, 701, 2712) && between(end, 6647, 6662), 1, 0)) %>%
  dplyr::mutate(group = case_when((twentyS == 1) ~ "20S",
                                  (twentysevenSAtwo == 1) ~ "27SA2",
                                  (eighteenS == 1) ~ "18S",
                                  (twentyfiveS == 1) ~ "25S",
                                  (fiftyeightS == 1) ~ "5.8S",
                                  (twentysevenSAthree == 1) ~ "27SA3",
                                  (twentysevenSBS == 1) ~ "27SBS",
                                  (twentysevenSBL == 1) ~ "27SBL",
                                  (sevenSS == 1) ~ "7SS",
                                  (sevenSL == 1) ~ "7SL",
                                  (twentyfivepointfiveS == 1) ~ "25.5S",
                                  (fiftyeightSL == 1) ~ "5.8SL"
                                  (thirtyfiveS == 1) ~ "35S",
                                  (thirtythreeS == 1) ~ "33S",
                                  (thritytwoS == 1) ~ "32S"))

pre_rRNA_good_20S <- pre_rRNA_good_sort_each %>%
  dplyr::filter(group == "20S")
pre_rRNA_good_27SA2 <- pre_rRNA_good_sort_each %>%
  dplyr::filter(group == "27SA2")
pre_rRNA_good_18S <- pre_rRNA_good_sort_each %>%
  dplyr::filter(group == "18S")
pre_rRNA_good_25S <- pre_rRNA_good_sort_each %>%
  dplyr::filter(group == "25S")
pre_rRNA_good_5.8S <- pre_rRNA_good_sort_each %>%
  dplyr::filter(group == "5.8S")
pre_rRNA_good_27SA3 <- pre_rRNA_good_sort_each %>%
  dplyr::filter(group == "27SA3")
pre_rRNA_good_27SBS <- pre_rRNA_good_each %>%
  dplyr::filter(group == "27SBS")
pre_rRNA_good_27SBL <- pre_rRNA_good_sort_each %>%
  dplyr::filter(group == "27SBL")
pre_rRNA_good_7SS <- pre_rRNA_good_sort_each %>%
  dplyr::filter(group == "7SS")
pre_rRNA_good_7SL <- pre_rRNA_good_sort_each %>%
  dplyr::filter(group == "7SL")
pre_rRNA_good_26S <- pre_rRNA_good_sort_each %>%
  dplyr::filter(group == "25.5S")
pre_rRNA_good_5.8SL <- pre_rRNA_good_sort_each %>%
  dplyr::filter(group == "5.8SL")
pre_rRNA_good_35S <- pre_rRNA_good_sort_each %>%
  dplyr::filter(group == "35S")
pre_rRNA_good_33S <- pre_rRNA_good_sort_each %>%
  dplyr::filter(group == "33S")
pre_rRNA_good_32S <- pre_rRNA_good_sort_each %>%
  dplyr::filter(group == "32S")


