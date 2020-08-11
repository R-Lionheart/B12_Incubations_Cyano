source("B12_Inc_Functions.R")
# source("src/biostats.R")
options(scipen = 999)

library(tidyverse) 
library(vegan)

# User data
pattern = "BMIS"
percentMissing = 0.5

# Functions
makeWide <- function(df) {
  df.wide <- df %>%
    ungroup() %>%
    tidyr::spread(Replicate.Name, Adjusted.Area) %>%
    as.data.frame()
  
  df.rownames <- df.wide[,-1]
  rownames(df.rownames) <- df.wide[,1]
  
  df.rownames[is.na(df.rownames)] <- NA
  
  #df.noNA <- na.omit(df.rownames)
  
  return(df.rownames)
}

RemoveCsv <- function(full.filepaths) {
  # Remove a .csv file extension and obtain basename from a given list of filepaths.
  #
  # Args
  #   Character strings of filepaths in a directory.
  #
  # Returns
  #   Character strings of file basenames, without a csv extension.
  #
  no.path <- substr(full.filepaths, 1, nchar(full.filepaths)-4)
  no.ID <-   gsub("\\_.*","", no.path)
  
  return(no.path)
}


# Data import and first filtering of unnecessary SampIDs --------------------------------------------
filename <- RemoveCsv(list.files(path = 'data_processed/', pattern = pattern))
filepath <- file.path('data_processed', paste(filename, ".csv", sep = ""))

Cyano_all <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE)) %>%
  select(Mass.Feature, runDate:replicate, Adjusted.Area) %>%
  unite(Replicate.Name, runDate:replicate, sep = "_") %>%
  filter(!str_detect(Replicate.Name, "Sept29QC|TruePoo|TruePooWeek2|TruePooWeek3|TruePooWeek4|DSW700m")) %>%
  filter(!Mass.Feature == "Inj_vol") %>%
  filter(!str_detect(Mass.Feature, ","))

Cyano_fixed <- Cyano_all %>%
  mutate(Replicate.Name = recode(Replicate.Name, 
                                 "171013_Smp_IT0_1" ="171013_Smp_IL1IT0_1", 
                                 "171013_Smp_IT0_2" = "171013_Smp_IL1IT0_2",
                                 "171013_Smp_IT0_3" = "171013_Smp_IL1IT0_3",
                                 "171013_Smp_IT05um_1" = "171013_Smp_IL1IT05um_1",
                                 "171013_Smp_IT05um_2" = "171013_Smp_IL1IT05um_2",
                                 "171013_Smp_IT05um_3" = "171013_Smp_IL1IT05um_3",
                                 "171030_Smp_IT0_1" = "171030_Smp_IL2IT0_1",
                                 "171030_Smp_IT0_2" = "171030_Smp_IL2IT0_2",
                                 "171030_Smp_IT0_3" = "171030_Smp_IL2IT0_3",
                                 "171030_Smp_IT05um_1" = "171030_Smp_IL2IT05um_1",
                                 "171030_Smp_IT05um_2" = "171030_Smp_IL2IT05um_2",
                                 "171030_Smp_IT05um_3" = "171030_Smp_IL2IT05um_3"))

#Then turn it back into a factor with the levels in the correct order
Cyano_fixed$Mass.Feature <- factor(Cyano_fixed$Mass.Feature, levels=unique(Cyano_fixed$Mass.Feature))
all.cyano <- ggplot(Cyano_fixed, aes(x = Mass.Feature, y = Adjusted.Area)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("B12 Incubations: Raw Cyano Area")
print(all.cyano)

Cyano_filtered <- Cyano_fixed %>%
  group_by(Mass.Feature) %>%
  mutate(Missing = sum(is.na(Adjusted.Area))) %>%
  mutate(MF.Count = n()) %>%
  filter(!Missing > (percentMissing*MF.Count)) %>%
  select(-c("Missing", "MF.Count"))
