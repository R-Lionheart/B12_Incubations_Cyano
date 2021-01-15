library(tidyverse)
options(scipen = 999)

replace_nonvalues <- function(x) (gsub(NaN, NA, x))

## Import and combine relevant Vitamins and Incubation files.

Vitamins_1500QC <- read.csv("data_processed/Skyline_QE_QC_Output_Vitamins_1500QC.csv") %>%
  slice(-1:-9) %>%
  select(-c(Description, Value)) %>%
  filter(!str_detect(Replicate.Name, 
                     "Blk|Std|TruePoo|CultureMED4|Sept29QC|TruePooWeek1|TruePooWeek2|
                     TruePooWeek3|TruePooWeek4|DSW700m|Process")) %>%
  mutate(QC.Level = "1500",
         Dataset = "Vitamins") %>%
  mutate(Area.with.QC = ifelse(str_detect(all.Flags, "Blank.Flag"), 
                               NA, Area.with.QC)) %>%
  select(Replicate.Name, Precursor.Ion.Name, Area, Area.with.QC, QC.Level, Dataset)

# Munge full Incubations dataset
Incubations <- read.csv("data_processed/MSDial_QE_QC_Output_B12-Incubations.csv",
                        stringsAsFactors = FALSE) %>%
  slice(-1:-6) %>%
  select(-c(Description, Value)) %>%
  filter(!str_detect(Replicate.Name, 
                     "Blk|Std|TruePoo|CultureMED4|Sept29QC|TruePooWeek1|TruePooWeek2|
                     TruePooWeek3|TruePooWeek4|DSW700m|Process")) %>%
  rename(Precursor.Ion.Name = Metabolite.Name,
         Area = Area.Value) %>%
  mutate(QC.Level = "10000",
         Dataset = "Incubations") %>%
  mutate(Area.with.QC = ifelse(str_detect(all.Flags, "Blank.Flag"), 
                               NA, Area.with.QC)) %>%
  select(Replicate.Name, Precursor.Ion.Name, Area, Area.with.QC, QC.Level, Dataset)


# Join together
Complete.Dataset <- Vitamins_1500QC %>%
  rbind(Incubations) %>%
  filter(!str_detect(Precursor.Ion.Name, ","),
         !str_detect(Replicate.Name, "Inj")) %>%
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
                                 "171030_Smp_IT05um_3" = "171030_Smp_IL2IT05um_3")) %>%
  mutate(Size.Fraction = ifelse(str_detect(Replicate.Name, "5um"), "Large.Filter", "Small.Filter"),
         Eddy = ifelse(str_detect(Replicate.Name, "IL1"), "Cyclonic", "Anticyclonic")) %>%
  group_by(Precursor.Ion.Name, Eddy, Size.Fraction) %>%
  mutate(Binned.Group = ifelse(str_detect(Replicate.Name, "DSW"), "DeepSeaWater",
                               ifelse(str_detect(Replicate.Name, "T0"), "TimeZero",
                                      ifelse(str_detect(Replicate.Name, "Control"), "Control",
                                             ifelse(str_detect(Replicate.Name, "DMBnoBT|noBT"), "Nutrients", "NoNutrients"))))) %>%
  mutate(Binned.Group = ifelse(str_detect(Replicate.Name, "5um"), paste(Binned.Group, "LargeFilter", sep = "_"),
                               paste(Binned.Group, "SmallFilter", sep = "_"))) %>%
  mutate(Binned.Group = ifelse(Eddy == "Cyclonic", paste(Binned.Group, "Cyclonic", sep = "_"), 
                               paste(Binned.Group, "Anticyclonic", sep = "_")))

write.csv(Complete.Dataset, "data_processed/Vitamins_Incubations_CompleteDataset.csv")

# Filter Vitamins B2-IS and DMB from total dataset averages
Complete.Dataset.Avg <- Complete.Dataset %>%
  filter(!Precursor.Ion.Name == "B2-IS|DMB") %>%
  separate(Replicate.Name, into = c("Date", "runtype", "Supergroup", "replicate"), remove = FALSE) %>%
  group_by(Precursor.Ion.Name, Binned.Group, Size.Fraction, Eddy) %>%
  mutate(Area.mean = mean(Area, na.rm = TRUE),
         Area.with.QC.mean = mean(Area.with.QC, na.rm = TRUE)) %>%
  mutate_at(c("Area.with.QC.mean"), replace_nonvalues) %>%
  mutate(Area.with.QC.mean = as.numeric(Area.with.QC.mean)) %>%
  select(Precursor.Ion.Name, Size.Fraction, Eddy,
         Binned.Group, Area.with.QC.mean, Dataset) %>%
  unique() %>%
  arrange(Binned.Group) 
write.csv(Complete.Dataset.Avg, "data_processed/Vitamins_Incubations_AvgCompleteDataset.csv")
