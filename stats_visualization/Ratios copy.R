library(ggplot2)
library(tidyverse)

source("src/B12_Functions.R")
## Specific compound ratios

Glutamic.Acid.RF <- 1.83
Glutamine.RF <- 0.912

# Import + clean BMIS files
dataset.pattern <- "Time0"

filename <- RemoveCsv(list.files(path = "data_processed/", pattern = dataset.pattern))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))

BMISd <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, 
                                               check.names = FALSE)) %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "four"), fill = "right", remove = FALSE) %>%
  select(Mass.Feature, SampID, Adjusted.Area) %>%
  drop_na() %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(Averages = mean(Adjusted.Area, na.rm = TRUE)) %>%
  select(Mass.Feature, SampID, Averages) %>%
  unique()


# # 1_0.2
# Dataset$SampID <- factor(Dataset$SampID,
#                               levels = c("IL1IT0", "IL1Control", "IL1DMB", "IL1WBT", "IL1DSW",
#                                          "IL1DMBnoBT", "IL1noBT"))
# 
# # 1_5
# Dataset$SampID <- factor(Dataset$SampID,
#                               levels = c("IL1IT05um", "IL1Control5um", "IL1DMB5um", "IL1WBT5um", "IL1DSW5um",
#                                          "IL1DMBnoBT5um", "IL1noBT5um"))
# # 2_0.2
# Dataset$SampID <- factor(Dataset$SampID,
#                               levels = c("IL2IT0", "IL2Control", "IL2DMB", "IL2WBT", "IL2DSW",
#                                          "IL2DMBnoBT", "IL2noBT"))
# # # 2_5
# Dataset$SampID <- factor(Dataset$SampID,
#                               levels = c("IL2IT05um", "IL2Control5um", "IL2DMB5um", "IL2WBT5um",
#                                          "IL2DSW5um","IL2DMBnoBT5um", "IL2noBT5um"))


### SAM SAH
makeRatiodf <- function(df, filtersize, eddy, mycompounds) {
  RatioDataframe <- df %>%
    filter(!str_detect(SampID, filtersize)) %>%
    filter(str_detect(SampID, eddy)) %>%
    filter(Mass.Feature %in% mycompounds) %>%
    group_by(SampID) %>%
    add_tally() %>%
    filter(!n < 2) %>%
    mutate(Ratios =
             Averages[Mass.Feature == mycompounds[1]] /
             Averages[Mass.Feature == mycompounds[2]])
  return(RatioDataframe)
}
makeRatioPlot <- function(df, mytitle) {
  Ratio.plot <- ggplot(df, aes(x = SampID, y = Averages, fill = Mass.Feature)) +
    geom_bar(stat = "identity", position = "fill") +
    geom_hline(yintercept = 0) +
    coord_flip() +
    ggtitle(paste(mytitle))
  print(Ratio.plot)
}

SAM.SAH_IL15um_df <- makeRatiodf(BMISd, filtersize = "5um", eddy = "IL1", 
                                  mycompounds = c("Adenosyl Methionine", "Adenosyl Homocysteine"))
SAM.SAH_IL15um_plot <- makeRatioPlot(SAM.SAH_IL15um_df, mytitle = "Cyclonic, 5um, SAM:SAH")

GG_IL15um_df <- makeRatiodf(BMISd, filtersize = "5um", eddy = "IL1",
                            mycompounds = c("Glutamic acid", "Glutamine"))
GG_IL15um_plot <- makeRatioPlot(GG_IL15um_df, mytitle = "Cyclonic, 5um, Glutamate:Glutamine")


# RF Normalization, other compound comparisons ----------------------------

# Glutamine : Glutamic acid Normed to RF
Dataset <- BMISd %>%
  filter(!str_detect(SampID, "5um")) %>%
  filter(str_detect(SampID, "IL1"))

glutamine.glutamate <- Dataset %>%
  filter(Mass.Feature %in% c("Glutamic acid", "Glutamine")) %>%
  mutate(Normd.to.RF = ifelse(Mass.Feature == "Glutamic acid", 
                              Averages / Glutamic.Acid.RF, Averages / Glutamine.RF)) %>%
  group_by(SampID) %>%
  add_tally() %>%
  filter(!n < 2) %>%
  mutate(Ratios =
           Normd.to.RF[Mass.Feature == "Glutamic acid"] /
           Normd.to.RF[Mass.Feature == "Glutamine"])

Plot2_5 <- ggplot(glutamine.glutamate, aes(x = SampID, y = Normd.to.RF, fill = Mass.Feature)) +
  geom_bar(stat = "identity", position = "fill") +
  #geom_text(aes(label = Ratios)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  ggtitle("Anticyclonic Eddy, 5um. Glutamine:Glutamate, RF-Normalized")
Plot2_5

# ketoglutaric Acid
KGA <- Dataset %>%
  filter(Mass.Feature == "Ketoglutaric Acid")

# Glutathione : GSSG
glutathione.gssg <- Dataset %>%
  filter(Mass.Feature %in% c("Glutathione Disulfide"))


