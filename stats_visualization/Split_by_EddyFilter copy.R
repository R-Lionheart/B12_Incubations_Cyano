library(tidyverse)
options(scipen = 999)

source("src/B12_Functions.R")

## Import and filtering for stats analysis 

# User data
BMIS.pattern = "BMIS_Output" # To upload the BMISd data
percentMissing = 0.5 # Use as percent cutoff for dropping compounds missing this percent of peaks.

# BMIS import and filtering of unnecessary SampIDs --------------------------------------------
filename <- RemoveCsv(list.files(path = "data_processed/", pattern = BMIS.pattern))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))

all.HILIC <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE)) %>%
  select(Mass.Feature, runDate:replicate, Adjusted.Area) %>%
  unite(Replicate.Name, runDate:replicate, sep = "_") %>%
  filter(!str_detect(Replicate.Name, "Sept29QC|TruePooWeek1|TruePooWeek2|TruePooWeek3|TruePooWeek4|DSW700m"),
         Mass.Feature != "Inj_vol",
         !str_detect(Mass.Feature, ",")) # This line drops Internal Standards. Temporary fix!

# Adjust for T0 naming issues
all.HILIC <- all.HILIC %>%
  mutate(Replicate.Name = recode(Replicate.Name, 
                                 "171002_Smp_IT0_1" ="171002_Smp_IL1IT0_1", 
                                 "171002_Smp_IT0_2" = "171002_Smp_IL1IT0_2",
                                 "171002_Smp_IT0_3" = "171002_Smp_IL1IT0_3",
                                 "171009_Smp_IT05um_1" = "171009_Smp_IL1IT05um_1",
                                 "171009_Smp_IT05um_2" = "171009_Smp_IL1IT05um_2",
                                 "171009_Smp_IT05um_3" = "171009_Smp_IL1IT05um_3",
                                 "171016_Smp_IT0_1" = "171016_Smp_IL2IT0_1",
                                 "171016_Smp_IT0_2" = "171016_Smp_IL2IT0_2",
                                 "171016_Smp_IT0_3" = "171016_Smp_IL2IT0_3",
                                 "171023_Smp_IT05um_1" = "171023_Smp_IL2IT05um_1",
                                 "171023_Smp_IT05um_2" = "171023_Smp_IL2IT05um_2",
                                 "171023_Smp_IT05um_3" = "171023_Smp_IL2IT05um_3"))
currentDate <- Sys.Date()
csvFileName <- paste("data_processed/BMISd_Time0_Fixed_", currentDate, ".csv", sep = "")
write.csv(all.HILIC, csvFileName, row.names = FALSE)


# Filter out compounds that are missing too many peaks --------------------------------------------
HILIC_filtered <- all.HILIC %>%
  group_by(Mass.Feature) %>%
  mutate(Missing = sum(is.na(Adjusted.Area))) %>%
  mutate(MF.Count = n()) %>%
  filter(!Missing > (percentMissing*MF.Count)) %>%
  select(-c("Missing", "MF.Count"))

# Separate groups into analysis and save to data_processed/ subdirectory.
IsoLagran1 <- HILIC_filtered %>%
  filter(str_detect(Replicate.Name, "IL1"))
IsoLagran2 <- HILIC_filtered %>%
  filter(str_detect(Replicate.Name, "IL2"))

IsoLagran1_0.2 <- IsoLagran1 %>%
  filter(!str_detect(Replicate.Name, "5um"))
IsoLagran1_5 <- IsoLagran1 %>%
  filter(str_detect(Replicate.Name, "5um"))
IsoLagran2_0.2 <- IsoLagran2 %>%
  filter(!str_detect(Replicate.Name, "5um"))
IsoLagran2_5 <- IsoLagran2 %>%
  filter(str_detect(Replicate.Name, "5um"))

write.csv(IsoLagran1_0.2, "data_processed/IsoLagran1_0.2_notstd.csv")
write.csv(IsoLagran1_5, "data_processed/IsoLagran1_5_notstd.csv")
write.csv(IsoLagran2_0.2, "data_processed/IsoLagran2_0.2_notstd.csv")
write.csv(IsoLagran2_5, "data_processed/IsoLagran2_5_notstd.csv")

