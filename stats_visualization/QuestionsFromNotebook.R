## Answering questions from notebook

library(gt)
library(kableExtra)
library(paletteer)
library(patchwork)
library(superheat)
library(tidyverse)

options(scipen = 999)
options(digits = 5)
source("src/Functions.R")

replace_nonvalues <- function(x) (gsub(NaN, NA, x))

## QC threshold is lowered to 1500 for vitamins, compared with ordinary 5000 level
## S/N flags need to be incorporated
## DMB should be removed from averaged group, otherwise spikes the overall. 


## Combine with bigger project, but keep the QC separate

# NUTRIENTS: [DMBnoB12, noB12], both containing a f/2 spike 

# NONUTRIENTS: [WB12, WDMB], no additional nutrients added.

# DEEPSEAWATER: [DSW], 10% Deep Sea Water (DSW) taken from 700 m, 
# filtered and evenly mixed with 90% in situ water from 25 m. This mixture was then incubated and filtered. 

# TIME 0: [T0] In situ condition at the start of the experiment, aka time 0. 
# Water that was filtered, notincubated, and frozen. 

# CONTROL: [Control], 100% in situ water from 25 m, incubated and then filtered.

replace_nonvalues <- function(x) (gsub(NaN, NA, x))

# Import files --------------------------------------------------
# Munge vitamin files
Vitamins_5000QC <- read.csv("data_processed/Skyline_QE_QC_Output_Vitamins_5000QC.csv", stringsAsFactors = FALSE) %>%
  slice(-1:-9) %>%
  select(-c(Description, Value)) %>%
  filter(!str_detect(Replicate.Name, 
                     "Blk|Std|TruePoo|CultureMED4|Sept29QC|TruePooWeek1|TruePooWeek2|
                     TruePooWeek3|TruePooWeek4|DSW700m|Process")) %>%
  separate(Replicate.Name, into = c("Date", "runtype", "Supergroup", "replicate"), remove = FALSE) %>%
  mutate(QC.Level = "5000") %>%
  select(Supergroup, Precursor.Ion.Name, Area, Area.with.QC, all.Flags, QC.Level) %>%
  mutate(Has.Flag = ifelse(is.na(all.Flags), "No Flag", "Flagged Peak"),
         Dataset = "Vitamins5000")

Vitamins_1500QC <- read.csv("data_processed/Skyline_QE_QC_Output_Vitamins_1500QC.csv") %>%
  slice(-1:-9) %>%
  select(-c(Description, Value)) %>%
  filter(!str_detect(Replicate.Name, 
                     "Blk|Std|TruePoo|CultureMED4|Sept29QC|TruePooWeek1|TruePooWeek2|
                     TruePooWeek3|TruePooWeek4|DSW700m|Process")) %>%
  separate(Replicate.Name, into = c("Date", "runtype", "Supergroup", "replicate"), remove = FALSE) %>%
  mutate(QC.Level = "1500") %>%
  select(Supergroup, Precursor.Ion.Name, Area, Area.with.QC, all.Flags, QC.Level) %>%
  mutate(Has.Flag = ifelse(is.na(all.Flags), "No Flag", "Flagged Peak"),
         Dataset = "Vitamins1500")

# Munge full Incubations dataset
Incubations <- read.csv("data_processed/MSDial_QE_QC_Output_B12-Incubations_2020-10-12.csv",
                        stringsAsFactors = FALSE) %>%
  slice(-1:-6) %>%
  select(-c(Description, Value)) %>%
  filter(!str_detect(Replicate.Name, 
                     "Blk|Std|TruePoo|CultureMED4|Sept29QC|TruePooWeek1|TruePooWeek2|
                     TruePooWeek3|TruePooWeek4|DSW700m|Process")) %>%
  separate(Replicate.Name, into = c("Date", "runtype", "Supergroup", "replicate"), remove = FALSE) %>%
  rename(Precursor.Ion.Name = Metabolite.Name,
         Area = Area.Value) %>%
  mutate(QC.Level = "5000") %>%
  select(Supergroup, Precursor.Ion.Name, Area, Area.with.QC, all.Flags, QC.Level) %>%
  mutate(Dataset = "Incubations")

## QC Comparisons  --------------------------------------------------
All.Vitamins <- Vitamins_1500QC %>%
  rbind(Vitamins_5000QC) %>%
  group_by(QC.Level, Has.Flag) %>%
  add_tally() %>%
  mutate(QC.Level = as.numeric(QC.Level))

## Compare QC filtering
ggplot(All.Vitamins, aes(x = QC.Level, y = n, fill = Has.Flag)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_text(aes(label=n), position=position_dodge(width = 0.9), 
            vjust = -0.25, check_overlap = TRUE) +
  annotate("text", x = 1500, y = 600, label = "510 + 414 = 924") +
  annotate("text", x = 6000, y = 600, label = "651 + 273 = 924") +
  ggtitle("Comparison of Area Min Flags in QC Levels: Vitamins")

ggsave("figures/QC_Comparison.png")


## Combine with full Cyano dataset --------------------------------------------------
Vitamins.for.Analysis <- All.Vitamins %>%
  ungroup() %>%
  filter(QC.Level == 1500) %>%
  select(Supergroup:Dataset, -Has.Flag)

# Join together
Complete.Dataset <- Vitamins.for.Analysis %>%
  rbind(Incubations) %>%
  filter(!str_detect(Precursor.Ion.Name, ",")) %>%
  mutate(Size.Fraction = ifelse(str_detect(Supergroup, "5um"), "Large.Filter", "Small.Filter"),
         Eddy = ifelse(str_detect(Supergroup, "IL1"), "Cyclonic", "Anticyclonic"),
         Dataset = "Vitamins") %>%
  group_by(Precursor.Ion.Name, Eddy, Size.Fraction) %>%
  mutate(Binned.Group = ifelse(str_detect(Supergroup, "DSW"), "DeepSeaWater",
                               ifelse(str_detect(Supergroup, "T0"), "TimeZero",
                                      ifelse(str_detect(Supergroup, "Control"), "Control",
                                             ifelse(str_detect(Supergroup, "DMBnoBT|noBT"), "Nutrients", "NoNutrients")))))
write.csv(Complete.Dataset, "data_intermediate/Vitamins_Incubations_CompleteDataset.csv")

# Filter Vitamins B2-IS and DMB from total dataset averages
Complete.Dataset.Avg <- Complete.Dataset %>%
filter(!Precursor.Ion.Name == "B2-IS" && Dataset == "Vitamins",
       !Precursor.Ion.Name == "DMB" && Dataset == "Vitamins") %>%
  group_by(Precursor.Ion.Name, Binned.Group, Size.Fraction, Eddy) %>%
  mutate(Area.mean = mean(Area, na.rm = TRUE),
         Area.with.QC.mean = mean(Area.with.QC, na.rm = TRUE)) %>%
  mutate_at(c("Area.with.QC.mean"), replace_nonvalues) %>%
  mutate(Area.with.QC.mean = as.numeric(Area.with.QC.mean)) %>%
  select(Precursor.Ion.Name, Size.Fraction, Eddy,
         Binned.Group, Area.mean, Area.with.QC.mean, Dataset) %>%
  unique() %>%
  arrange(Binned.Group)
write.csv(Complete.Dataset.Avg, "data_intermediate/Vitamins_Incubations_CompleteDatasetAvg.csv")


 
# Notebook questions --------------------------------------------------------
# Which compounds exist in which size fraction?
Which.Cmpds <- Complete.Dataset %>%
  group_by(Precursor.Ion.Name, Size.Fraction) %>%
  mutate(
    NumberOfSamples = n(),
    Number.Missing = sum(is.na(Area.with.QC)),
    Percent.Present = 1 - (sum(is.na(Area.with.QC))/NumberOfSamples),
    Over.Half.Missing = ifelse(Percent.Present >= 0.5, FALSE, TRUE)
  ) %>%
  ungroup()

ggplot(Which.Cmpds, aes(x = reorder(Precursor.Ion.Name, -Percent.Present), 
                        y = Percent.Present)) +
  facet_wrap(~Size.Fraction) +
  geom_bar(position = "dodge", stat = "identity", 
           aes(fill = as.numeric(row.names(Which.Cmpds)) %% 2 == 0)) +
  scale_fill_discrete(guide = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25)) +
  scale_x_discrete(name ="Precursor.Ion.Name") +
  scale_y_continuous(name = "Percent of Samples Present in Dataset") +
  ggtitle("Which Compounds Exist in Which Size Fractions?")
ggsave("figures/Which.Compounds.Present.png")

Table.DF <- Which.Cmpds %>%
  select(Precursor.Ion.Name, Percent.Present) %>%
  filter(Percent.Present == 1 | Percent.Present == 0) %>%
  unique() %>%
  arrange(Percent.Present)
write.csv(Table.DF, "data_intermediate/WhichCompoundsTable.csv")

## Do nutrients change with additions/size fractions?
Nutrient.Change <- Complete.Dataset.Avg %>%
  filter(Eddy == "Anticyclonic") %>%
  select(Precursor.Ion.Name, Area.with.QC.mean, Binned.Group, Size.Fraction, Eddy, Dataset) %>%
  unique() %>%
  filter(!is.na(Area.with.QC.mean)) %>%
  mutate(Area.with.QC.mean = as.numeric(Area.with.QC.mean)) %>%
  mutate(Dataset.Position = ifelse(Area.with.QC.mean > 633000, 
                                   "Higher.Average.Area", "Lower.Average.Area"))

ggplot(Nutrient.Change, aes(x = reorder(Precursor.Ion.Name, -Area.with.QC.mean), 
                            y = Area.with.QC.mean, fill = Binned.Group)) +
  geom_bar(position = "dodge", stat = "identity") +
  coord_flip() +
  facet_wrap(~Size.Fraction + Dataset.Position, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_y_continuous("QC'd Average Area") +
  scale_x_discrete("Precursor Ion Name") +
  ggtitle("Anticyclonic: Do nutrients change with additions/size fractions?")
ggsave("figures/DoNutrientsChange.png")


## Do compound proportions stay the same w/ DSW?
DSWproportions <- Complete.Dataset.Avg %>%
  filter(Eddy == "Anticyclonic",
         Size.Fraction == "Large.Filter") %>%
  select(Precursor.Ion.Name, Binned.Group, Area.with.QC.mean, Size.Fraction, Eddy) %>%
  unique() 

ggplot(DSWproportions, aes(x = Binned.Group, y = Area.with.QC.mean, fill = Precursor.Ion.Name)) +
  geom_bar(position = "fill", stat = "identity") +
  coord_flip()
ggsave("figures/DoCompoundProportionsChange_BetaineCholine.png")

# Hopefully will work at some point...
# ggplot(data=DSWproportions, aes(x=Area.with.QC.mean, group=Precursor.Ion.Name, fill=Precursor.Ion.Name)) +
#   facet_wrap(~Binned.Group) +
#   geom_density(adjust=1.5, alpha=.4) 
# ggsave("figures/DensityProportions_BetaineCholine.png")


## Ask a compound: how much of all of the compound are in large/small size fraction? Assume n peak area 0.2
# + n peak area 5 = 100%.

SF.ratios <- Complete.Dataset.Avg %>%
  ungroup() %>%
  filter(Eddy == "Anticyclonic") %>%
  mutate_at(c("Area.with.QC.mean"), replace_nonvalues) %>%
  mutate(Area.with.QC.mean = as.numeric(Area.with.QC.mean)) %>%
  group_by(Precursor.Ion.Name, Size.Fraction) %>%
  mutate(Sum.per.SF = sum(Area.with.QC.mean, na.rm = TRUE)) %>%
  select(Precursor.Ion.Name, Eddy, Size.Fraction, Sum.per.SF, Binned.Group) %>%
  unique() %>%
  filter(!Sum.per.SF == 0) %>%
  group_by(Precursor.Ion.Name) %>%
  add_tally() %>%
  filter(!n == 5) 

SF.ratios$Precursor.Ion.Name <- factor(SF.ratios$Precursor.Ion.Name, levels=rev(unique(SF.ratios$Precursor.Ion.Name)))

ggplot(SF.ratios, aes(x = Sum.per.SF, y = Precursor.Ion.Name, fill = Size.Fraction)) +
  facet_wrap(~Binned.Group) +
  geom_bar(position = "fill", stat = "identity") +
  ggtitle("Anticyclonic: How much of the compound is in the large and small size fraction?")
ggsave("figures/HowMuchinSizeFractions.png")

## Which ones were completely in one or the other
# Check these graphs!! last one is artifice
## Box and whisker for how these plots are behaving similarly. Which compounds are behaving the same in all treatments?
## Add line to that to show moving directions of compounds.
## Cluster them to see whos moving together. Stuff is always going up in deep sea water? Always going down in Time0?
## Also just do a set of amino acids for all the graphs. 

## Good ionizers will be more abundant. 
## Ask Sonya (for example) who makes ectoine? 
