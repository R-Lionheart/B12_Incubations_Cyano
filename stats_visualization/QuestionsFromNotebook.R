## Answering questions from notebook

library(forcats)
library(tidyverse)

options(scipen = 999)
options(digits = 5)
source("src/Functions.R")

replace_nonvalues <- function(x) (gsub(NaN, NA, x))

## QC threshold is lowered to 1500 for vitamins, compared with ordinary 5000 level

# NUTRIENTS: [DMBnoB12, noB12], both containing a f/2 spike 

# NONUTRIENTS: [WB12, WDMB], no additional nutrients added.

# DEEPSEAWATER: [DSW], 10% Deep Sea Water (DSW) taken from 700 m, 
# filtered and evenly mixed with 90% in situ water from 25 m. This mixture was then incubated and filtered. 

# TIME 0: [T0] In situ condition at the start of the experiment, aka time 0. 
# Water that was filtered, notincubated, and frozen. 

# CONTROL: [Control], 100% in situ water from 25 m, incubated and then filtered.

# Import files --------------------------------------------------
Standards <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv",
                               stringsAsFactors = FALSE, header = TRUE) %>%
  select(Compound.Name_old, Compound.Type) %>%
  rename(Precursor.Ion.Name = Compound.Name_old)

Complete.Dataset <- read.csv("data_processed/Vitamins_Incubations_CompleteDataset.csv")
Complete.Dataset.Avg <- read.csv("data_processed/Vitamins_Incubations_AvgCompleteDataset.csv")
# Notebook questions --------------------------------------------------------
# Which compounds exist in which size fraction?
Which.Cmpds <- Complete.Dataset %>%
  group_by(Precursor.Ion.Name, Size.Fraction) %>%
  mutate(
    Number.Of.Samples = n(),
    Number.Missing = sum(is.na(Area.with.QC)),
    Percent.Present = 1 - (sum(is.na(Area.with.QC))/Number.Of.Samples),
    Over.Half.Missing = ifelse(Percent.Present >= 0.5, FALSE, TRUE)
  ) %>%
  ungroup()

ggplot(Which.Cmpds, aes(x = reorder(Precursor.Ion.Name, -Percent.Present),
                        y = Percent.Present)) +
  facet_wrap(~Size.Fraction + Eddy) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_discrete(guide = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25)) +
  scale_x_discrete(name ="Precursor.Ion.Name") +
  scale_y_continuous(name = "Percent of Samples Present in Dataset") +
  ggtitle("Which Compounds Exist in Which Size Fractions?")
#ggsave("figures/Which.Compounds.Present.png")

# Which compounds were completely in one or the other size fraction?
WhichCompounds <- function(df, EddyDirection) {
  Table.DF <- df %>%
    filter(Eddy == EddyDirection) %>%
    group_by(Precursor.Ion.Name) %>%
    mutate(Number.Of.Samples = n(),
           Number.Missing = sum(is.na(Area.with.QC.mean))) %>%
    filter(Number.Of.Samples != Number.Missing) %>%
    group_by(Precursor.Ion.Name, Size.Fraction) %>%
    mutate(Is.Present = ifelse(is.na(Area.with.QC.mean), TRUE, FALSE)) %>%
    group_by(Precursor.Ion.Name, Size.Fraction) %>%
    mutate(Group.Present = ifelse(any(Is.Present == FALSE), TRUE, FALSE)) %>%
    ungroup() %>%
    group_by(Precursor.Ion.Name) %>%
    mutate(Filter.Check = ifelse(any(Group.Present == FALSE), TRUE, FALSE)) %>%
    filter(Filter.Check == TRUE) %>%
    filter(Group.Present == TRUE) %>%
    select(Precursor.Ion.Name, Size.Fraction, Eddy) %>%
    unique()
  return(Table.DF)
}

Table.DF.Cyclonic <- WhichCompounds(Complete.Dataset.Avg, "Cyclonic")
Table.DF.Anticyclonic <- WhichCompounds(Complete.Dataset.Avg, "Anticyclonic")
write.csv(Table.DF.Cyclonic, "data_intermediate/WhichCompoundsSizeFraction_Cyclonic.csv")
write.csv(Table.DF.Anticyclonic, "data_intermediate/WhichCompoundsSizeFraction_Anticyclonic.csv")

## Do nutrients change with additions/size fractions?
Nutrient.Change <- Complete.Dataset.Avg %>%
  filter(Eddy == "Anticyclonic") %>%
  select(Precursor.Ion.Name, Area.with.QC.mean, Binned.Group, Size.Fraction, Eddy, Dataset) %>%
  unique() %>%
  mutate(SampID = str_extract(Binned.Group, "[^_]+")) %>%
  filter(!is.na(Area.with.QC.mean)) %>%
  mutate(Area.with.QC.mean = as.numeric(Area.with.QC.mean)) %>%
  mutate(Dataset.Position = ifelse(Area.with.QC.mean > 633000,
                                   "Higher.Average.Area", "Lower.Average.Area")) %>%
  arrange(desc(Area.with.QC.mean)) 

ggplot(Nutrient.Change, aes(x = reorder(Precursor.Ion.Name, -Area.with.QC.mean),
                            y = Area.with.QC.mean, fill = SampID)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_wrap(~Dataset.Position + Size.Fraction, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous("QC'd Average Area") +
  scale_x_discrete("Precursor Ion Name") +
  ggtitle("Anticyclonic: Do nutrients change with additions/size fractions?")
# ggsave("figures/DoNutrientsChange.png")


## Do compound proportions stay the same w/ DSW, based on compound type?
DSWproportions <- Complete.Dataset.Avg %>%
  mutate(Precursor.Ion.Name = as.character(Precursor.Ion.Name)) %>%
  left_join(Standards) %>%
  filter(!is.na(Compound.Type)) %>%
  filter(Eddy == "Anticyclonic",
         Size.Fraction == "Large.Filter")%>%
  group_by(Binned.Group, Compound.Type) %>%
  summarize(Sum.per.CT = sum(Area.with.QC.mean, na.rm = TRUE)) %>%
  filter(Sum.per.CT != 0)

Compound.Type.Levels <- DSWproportions %>%
  group_by(Compound.Type) %>%
  summarize(Total.Sum = sum(Sum.per.CT)) %>%
  arrange(desc(Total.Sum)) %>%
  pull(Compound.Type)

DSWproportions.with.Factor <- DSWproportions %>%
  mutate(Compound.Type = factor(Compound.Type, levels = Compound.Type.Levels)) 
DSWproportions.with.Factor$Binned.Group <- factor(DSWproportions.with.Factor$Binned.Group,
                                                  levels=c("DeepSeaWater_LargeFilter_Anticyclonic",
                                                           "Control_LargeFilter_Anticyclonic",
                                                           "TimeZero_LargeFilter_Anticyclonic",
                                                           "Nutrients_LargeFilter_Anticyclonic",
                                                           "NoNutrients_LargeFilter_Anticyclonic"))

ggplot(DSWproportions.with.Factor, aes(x = Binned.Group, y = Sum.per.CT, fill = Compound.Type)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_brewer(palette = "Dark2") +
  ggtitle("Deep Sea Water Proportions")

#ggsave("figures/DoCompoundProportionsChange_CmpdType.png")


## Ask a compound: how much of all of the compound are in large/small size fraction? 
# Assume n peak area 0.2
# + n peak area 5 = 100%.

plotchart <- function(data) {
  # Quick plotting function for size fraction ratios.
  plot <- ggplot(data) +
    aes(x = reorder(Precursor.Ion.Name, Order), 
        y = Percent,
        fill = Size.Fraction) +
    geom_bar(stat = "identity", position = "fill") +
    coord_flip() +
    ggtitle(paste(unique(data[[1]][1])))
  
  print(plot)
}

test <- Complete.Dataset.Avg %>%
  mutate(Binned.Group = str_extract(Binned.Group, "[^_]+"))
  

Size.Fraction.Ratios <- test %>%
  ungroup() %>%
  filter(Eddy == "Anticyclonic") %>%
  mutate_at(c("Area.with.QC.mean"), replace_nonvalues) %>%
  mutate(Area.with.QC.mean = as.numeric(Area.with.QC.mean)) %>%
  group_by(Precursor.Ion.Name, Binned.Group, Size.Fraction) %>%
  mutate(Sum.per.SF = sum(Area.with.QC.mean, na.rm = TRUE)) %>%
  select(Precursor.Ion.Name, Eddy, Size.Fraction, Sum.per.SF, Binned.Group) %>%
  unique() %>%
  filter(Sum.per.SF != 0) %>%
  group_by(Precursor.Ion.Name) %>%
  add_tally() %>%
  filter(n != 5) %>%
  select(-c(Eddy, n)) %>%
  ungroup()
  
SF.Ratios.Plot <- Size.Fraction.Ratios %>%
  complete(nesting(Binned.Group, Precursor.Ion.Name), Size.Fraction, fill = list(Sum.per.SF = 0)) %>%
  group_by(Binned.Group, Precursor.Ion.Name) %>%
  mutate(Percent = Sum.per.SF / sum(Sum.per.SF)) %>%
  mutate(Large.Filter.Percent = Percent[Size.Fraction == "Large.Filter"]) %>%
  arrange(Binned.Group, desc(Large.Filter.Percent), desc(Precursor.Ion.Name), Size.Fraction) %>%
  ungroup() %>%
  mutate(Order = row_number()) %>%
  group_by(Binned.Group) %>%
  group_split() 

names(SF.Ratios.Plot) <- unique(Size.Fraction.Ratios$Binned.Group)

plotlist <- lapply(SF.Ratios.Plot, plotchart)


## Box and whisker for how these plots are behaving similarly. Which compounds are behaving the same in all treatments?
## Add line to that to show moving directions of compounds.
## Cluster them to see whos moving together. Stuff is always going up in deep sea water? Always going down in Time0?
## Also just do a set of amino acids for all the graphs. 

## Venn diagram for overlapping compounds between experiments
## 3D intensity visualizations

## Good ionizers will be more abundant. 
## Does the pathway change under different scenarios?
## Ask Sonya (for example) who makes ectoine? 