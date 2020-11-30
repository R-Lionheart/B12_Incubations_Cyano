library(tidyverse)

options(scipen = 999)
options(digits = 5)
source("src/Functions.R")

## Recreating Angie's graphs for the Cyano full dataset w/ Vitamins

# The mass features that are significantly different between treatments in 
# all SF, Eddy variations (four total).

# y axis is peak, x axis is treatment, facet wrap is eddy/sf.

test <- read.csv("data_processed/Vitamins_Incubations_AvgCompleteDataset.csv") %>%
  select(Precursor.Ion.Name, Area.with.QC.mean, Eddy, Size.Fraction, Replicate.Bin.Group) %>%
  filter(Precursor.Ion.Name == "Acetyl-L-carnitine") %>%
  mutate(Replicate.Bin.Group = substr(Replicate.Bin.Group, 1, nchar(Replicate.Bin.Group)-2)) %>%
  unique() %>%
  unite("GroupingID", Eddy:Size.Fraction, sep = "_")

ggplot(test, aes(x = GroupingID, y = Area.with.QC.mean)) +
  geom_boxplot()