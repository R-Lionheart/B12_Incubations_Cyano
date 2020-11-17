library(tidyverse)
library(vegan)
options(scipen = 999)

source("src/B12_Functions.R")

# Upload required files: BMISd and split by eddy and fraction.
Isolagran.pattern <- "IsoLagran"
filenames <- RemoveCsv(list.files(path = "data_processed/", pattern = Isolagran.pattern))
for (i in filenames) {
  filepath <- file.path("data_processed", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE))
}

# Upload Chlorophyll data
Chlorophyll.data <- read.csv("data_raw/Dyhrman_MS_Chla.csv", stringsAsFactors = FALSE) %>%
  select(1:6) %>%
  rename(Filter.size = Filter..µm.,
         Chla = Chl.a..µg.L.) %>%
  mutate(Date = as.character(Date),
         Eddy = ifelse(str_detect(Date, "13|9"), "IL2", "IL1")) %>%
  unite("Replicate.Name", SAMPLE, Eddy, sep = "_", remove = TRUE) %>%
  filter(Filter.size == 5.0) %>%
  mutate(Chla = as.numeric(Chla))


# Fix Replicate Name labels 
Chlorophyll_fixed <- Chlorophyll.data %>%
  mutate(Replicate.Name = recode(Replicate.Name,
                                 "C-1_IL1" ="171009_Smp_IL1Control5um_1",
                                 "C-1_IL2" = "171009_Smp_IL2Control5um_1",
                                 "C-2_IL1" = "171009_Smp_IL1Control5um_2",
                                 "C-2_IL2" = "171009_Smp_IL2Control5um_2",
                                 "C-3_IL1" = "171009_Smp_IL1Control5um_3",
                                 "C-3_IL2" = "171009_Smp_IL2Control5um_3",
                                 "ØBT-1_IL1" = "171009_Smp_IL1noBT5um_1",
                                 "ØBT-1_IL2" = "171009_Smp_IL2noBT5um_1",
                                 "ØBT-2_IL1" = "171009_Smp_IL1noBT5um_2",
                                 "ØBT-2_IL2" = "171009_Smp_IL2noBT5um_2",
                                 "ØBT-3_IL1" = "171009_Smp_IL1noBT5um_3",
                                 "ØBT-3_IL2" = "171009_Smp_IL2noBT5um_3",
                                 "DMB-1_IL1" = "171009_Smp_IL1DMB5um_1",
                                 "DMB-1_IL2" = "171009_Smp_IL2DMB5um_1",
                                 "DMB-2_IL1" = "171009_Smp_IL1DMB5um_2",
                                 "DMB-2_IL2" = "171009_Smp_IL2DMB5um_2",
                                 "DMB-3_IL1" = "171009_Smp_IL1DMB5um_3",
                                 "DMB-3_IL2" = "171009_Smp_IL2DMB5um_3",
                                 "DMBØBT-1_IL1" = "171009_Smp_IL1DMBnoBT5um_1",
                                 "DMBØBT-1_IL2" = "171009_Smp_IL2DMBnoBT5um_1",
                                 "DMBØBT-2_IL1" = "171009_Smp_IL1DMBnoBT5um_2",
                                 "DMBØBT-2_IL2" = "171009_Smp_IL2DMBnoBT5um_2",
                                 "DMBØBT-3_IL1" = "171009_Smp_IL1DMBnoBT5um_3",
                                 "DMBØBT-3_IL2" = "171009_Smp_IL2DMBnoBT5um_3",
                                 "DSW-1_IL1" = "171009_Smp_IL1DSW5um_1",
                                 "DSW-1_IL2" = "171009_Smp_IL2DSW5um_1",
                                 "DSW-2_IL1" = "171009_Smp_IL1DSW5um_2",
                                 "DSW-2_IL2" = "171009_Smp_IL2DSW5um_2",
                                 "DSW-3_IL1" = "171009_Smp_IL1DSW5um_3",
                                 "DSW-3_IL2" = "171009_Smp_IL2DSW5um_3",
                                 "IS-1_IL1" = "171009_Smp_IL1IT05um_1",
                                 "IS-1_IL2" = "171009_Smp_IL2IT05um_1",
                                 "IS-2_IL1" = "171009_Smp_IL1IT05um_2",
                                 "IS-2_IL2" = "171009_Smp_IL2IT05um_2",
                                 "IS-3_IL1" = "171009_Smp_IL1IT05um_3",
                                 "IS-3_IL2" = "171009_Smp_IL2IT05um_3",
                                 "WBT-1_IL1" = "171009_Smp_IL1WBT5um_1",
                                 "WBT-1_IL2" = "171009_Smp_IL2WBT5um_1",
                                 "WBT-2_IL1" = "171009_Smp_IL1WBT5um_2",
                                 "WBT-2_IL2" = "171009_Smp_IL2WBT5um_2",
                                 "WBT-3_IL1" = "171009_Smp_IL1WBT5um_3",
                                 "WBT-3_IL2" = "171009_Smp_IL2WBT5um_3")) %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "replicate")) %>%
  select(SampID, replicate, Chla) %>%
  unite(SampID, replicate, col = "Replicate.Name")

write.csv(Chlorophyll_fixed, "data_processed/Renamed_Chlorophyll.csv")

# Fix BMISd names to match new chlorophyll names
BMIS_fixed_IL15 <- FixBMISNames(IsoLagran1_5_notstd)
BMIS_fixed_IL25 <- FixBMISNames(IsoLagran2_5_notstd)

# Split Chlorophyll by Eddy
Chlorophyll.data_IL15 <- Chlorophyll.data %>% filter(str_detect(Replicate.Name, "IL1"))
Chlorophyll.data_IL25 <- Chlorophyll.data %>% filter(str_detect(Replicate.Name, "IL2"))


# Normalize to chlorophyll ------------------------------------------------------------
NormalizeToChla <- function(df) {
  complete_df <- df %>%
    left_join(Chlorophyll_fixed %>% select(Replicate.Name, Chla)) %>%
    mutate(Chla = as.numeric(Chla)) %>%
    mutate(Normalized.by.Chla = Adjusted.Area/Chla) %>%
    filter(!Chla < 0) %>% # Remove negative Chl-A values
    select(Mass.Feature, Replicate.Name, Normalized.by.Chla) %>%
    na.omit()
  return(complete_df)
}

complete.set.IL15 <- NormalizeToChla(BMIS_fixed_IL15)
complete.set.IL25 <- NormalizeToChla(BMIS_fixed_IL25)


## Standardize
complete.wide.IL15 <- MakeWide(complete.set.IL15, "Normalized.by.Chla")
complete.wide.IL25 <- MakeWide(complete.set.IL25, "Normalized.by.Chla")

complete.wide.IL15[is.na(complete.wide.IL15)] <- 1000
complete.wide.IL25[is.na(complete.wide.IL25)] <- 1000

complete.wideT.IL15 <- t(complete.wide.IL15)
complete.wideT.IL25 <- t(complete.wide.IL25)

complete.wide.normalizedT.IL15 <- decostand(complete.wideT.IL15, method = "standardize", na.rm = TRUE)
complete.wide.normalizedT.IL25 <- decostand(complete.wideT.IL25, method = "standardize", na.rm = TRUE) 


## Save Chlorophyll files to data_processed/
write.csv(complete.wide.normalizedT.IL15, paste("data_processed/IL1_5um_ChlAnormd_std.csv", sep = ""))
write.csv(complete.set.IL15, paste("data_processed/IL1_5um_ChlAnormd_notstd.csv", sep = ""))

write.csv(complete.wide.normalizedT.IL25, paste("data_processed/IL2_5um_ChlAnormd_std.csv", sep = ""))
write.csv(complete.set.IL25, paste("data_processed/IL2_5um_ChlAnormd_notstd.csv", sep = ""))

