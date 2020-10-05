## Answering questions from notebook

library(gt)
library(kableExtra)
library(paletteer)
library(patchwork)
library(tidyverse)

options(scipen = 999)
options(digits = 1)
source("src/Functions.R")

## QC threshold is lowered to 1500 for vitamins
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


file.pattern <- "Vitamins|MSDial"

replace_nonvalues <- function(x) (gsub(NaN, NA, x))

# Import files --------------------------------------------------
filenames <- RemoveCsv(list.files(path = "data_processed/", pattern = file.pattern))

for (i in filenames) {
  filepath <- file.path("data_processed/", paste(i,".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE, check.names = TRUE))
}

# Munge vitamins file
Vitamins <- Skyline_QE_QC_Output_Vitamins_2020.10.05 %>%
  slice(-1:-9) %>%
  select(-c(Description, Value)) %>%
  filter(!str_detect(Replicate.Name, 
                     "Blk|Std|TruePoo|CultureMED4|Sept29QC|TruePooWeek1|TruePooWeek2|
                     TruePooWeek3|TruePooWeek4|DSW700m|Process")) %>%
  separate(Replicate.Name, into = c("Date", "runtype", "Supergroup", "replicate"), remove = FALSE) %>%
  select(Supergroup, Precursor.Ion.Name, Area, Area.with.QC, all.Flags)

# Munge full Incubations dataset
Incubations <- MSDial_QE_QC_Output_B12.Incubations_2020.10.05 %>%
  slice(-1:-6) %>%
  select(-c(Description, Value)) %>%
  filter(!str_detect(Replicate.Name, 
                     "Blk|Std|TruePoo|CultureMED4|Sept29QC|TruePooWeek1|TruePooWeek2|
                     TruePooWeek3|TruePooWeek4|DSW700m|Process")) %>%
  separate(Replicate.Name, into = c("Date", "runtype", "Supergroup", "replicate"), remove = FALSE) %>%
  rename(Precursor.Ion.Name = Metabolite.Name,
         Area = Area.Value) %>%
  select(Supergroup, Precursor.Ion.Name, Area, Area.with.QC, all.Flags) 


# Join together
Complete.Dataset <- Incubations %>%
  rbind(Vitamins) %>%
  mutate(SizeFraction = ifelse(str_detect(Supergroup, "5um"), "5um", "0.2um"),
         Eddy = ifelse(str_detect(Supergroup, "IL1"), "Cyclonic", "Anticyclonic"),
         Dataset = "Vitamins") %>%
  group_by(Precursor.Ion.Name, Eddy, SizeFraction) %>%
  mutate(Binned.Group = ifelse(str_detect(Supergroup, "DSW"), "DeepSeaWater",
                               ifelse(str_detect(Supergroup, "T0"), "TimeZero",
                                      ifelse(str_detect(Supergroup, "Control"), "Control",
                                             ifelse(str_detect(Supergroup, "DMBnoBT|noBT"), "Nutrients", "NoNutrients")))))

  
  

# Plot Preparation --------------------------------------------------------
Combined <- Groups.Binned %>%
  filter(!Precursor.Ion.Name == "B2-IS",
         !Precursor.Ion.Name == "DMB") %>%
  group_by(Precursor.Ion.Name, Binned.Group, SizeFraction, Eddy) %>%
  mutate(Area.mean = mean(Area, na.rm = TRUE),
         Area.with.QC.mean = mean(Area.with.QC, na.rm = TRUE)) %>%
  select(Supergroup, Precursor.Ion.Name, SizeFraction, Eddy, 
         Binned.Group, Area.with.QC, Area.mean, Area.with.QC.mean)

## Which compounds exist in which size fraction? 
Which.Cmpds <- Combined %>%
  mutate(
    NumberOfSamples = n(),
    Percent.Missing = (sum(is.na(Area.with.QC))/NumberOfSamples)
  )

MakeTable <- function(df, EddyVorticity) {
  myTable <- df %>%
    filter(Eddy == EddyVorticity) %>%
    mutate(Binned.Group = case_when(
      Binned.Group == "DeepSeaWater" ~ "Deep Sea Water. No added nutrients: only deep sea water from 700m",
      Binned.Group == "TimeZero" ~ "Time Zero. In situ condition at the start of the experiment.",
      Binned.Group == "noNutrients" ~ "No Nutrients. Consists of DMBnoB12 and noB12 samples.",
      Binned.Group == "Nutrients" ~ "Nutrients. Spikes added after collection: B12, DMB, and Control",
    )) %>% 
    mutate(Supergroup = factor(Supergroup, 
                               levels = c("IL1T0", "IL1Control", "IL1DMB", "IL1WBT", "IL1DSW", "IL1DMBnoBT", "IL1noBT",
                                          "IL1T05um", "IL1Control5um", "IL1DMB5um", "IL1WBT5um", "IL1DSW5um","IL1DMBnoBT5um", "IL1noBT5um",
                                          "IL2T0", "IL2Control", "IL2DMB", "IL2WBT", "IL2DSW", "IL2DMBnoBT", "IL2noBT",
                                          "IL2T05um", "IL2Control5um", "IL2DMB5um", "IL2WBT5um", "IL2DSW5um","IL2DMBnoBT5um", "IL2noBT5um"))) %>%
    dplyr::group_by(Eddy) %>%
    dplyr::group_by(SizeFraction, Precursor.Ion.Name) %>% 
    dplyr::summarize(
      NumberOfSamples = n(),
      PercentMissing = (sum(is.na(Area.with.QC))/NumberOfSamples)
    ) %>% 
    gt(rowname_col = "Supergroup") %>% 
    tab_header(title = md(paste(EddyVorticity, "Eddy")),
               subtitle = "Which Compounds Exist in Which Size Fractions?") %>%
    cols_align(align = "right", columns = TRUE) %>%
    data_color(
      columns = vars(PercentMissing),
      colors = scales::col_numeric(
        palette = paletteer::paletteer_d(
          palette = "ggsci::red_material"
        ) %>% as.character(),
        domain = NULL
      ),
      alpha = 0.5
    )
  print(myTable)
}

CyclonicTable <- MakeTable(Which.Cmpds, EddyVorticity = "Cyclonic")
AnticyclonicTable <- MakeTable(Which.Cmpds, EddyVorticity = "Anticyclonic")


# Assign factor levels to SampID and generate heatmap plot
Which.Cmpds$Binned.Group <- factor(Which.Cmpds$Binned.Group, 
  levels = c("Control", "DeepSeaWater", "TimeZero", "NoNutrients", "Nutrients"))

heatmap <- ggplot(data = Which.Cmpds, aes(x = Precursor.Ion.Name, y = Binned.Group, 
                                          fill = Percent.Missing)) + 
  facet_wrap(~SizeFraction + Eddy) +
  geom_tile(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 10)) +
  scale_y_discrete(limits = rev(levels(as.factor(Which.Cmpds$Binned.Group)))) +
  ggtitle("Where are Samples Missing?")
print(heatmap)


## Do nutrients change with additions/size fractions?
NutrientChange <- Combined %>%
  select(Precursor.Ion.Name, Area.with.QC.mean, Binned.Group, SizeFraction) %>%
  unique()

ggplot(NutrientChange, aes(x = reorder(Precursor.Ion.Name, -Area.with.QC.mean), y = Area.with.QC.mean,
                     fill = Binned.Group)) +
  geom_bar(position = "dodge", stat = "identity") +
  coord_flip() +
  facet_wrap(~SizeFraction) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_text(aes(label=round(Area.with.QC.mean)),
            position=position_dodge(width = 0.9), vjust = -1, hjust = -2) +
  scale_y_continuous("QC'd Average Area") +
  scale_x_discrete("Precursor Ion Name") +
  ggtitle("Anticyclonic: Do nutrients change with additions/size fractions?")

## Do compound proportions stay the same w/ DSW?
DSWproportions <- Combined %>%
  filter(Eddy == "Anticyclonic") %>%
  select(Precursor.Ion.Name, Binned.Group, Area.with.QC.mean, SizeFraction, Eddy) %>%
  unique()

ggplot(DSWproportions, aes(x = Binned.Group, y = Area.with.QC.mean, fill = Precursor.Ion.Name)) +
  geom_bar(position = "fill", stat = "identity") +
  coord_flip()

## Ask a compound: how much of all of the compound are in large/small size fraction? Assume x peak area 0.2
# + y peak area 5 = 100%.

## remove treatments where DMB and B12 are added
SF.ratios <- Combined %>%
  ungroup() %>%
  filter(Eddy == "Anticyclonic") %>%
  mutate_at(c("Area.with.QC.mean"), replace_nonvalues) %>%
  mutate(Area.with.QC.mean = as.numeric(Area.with.QC.mean)) %>%
  group_by(Precursor.Ion.Name, SizeFraction) %>%
  mutate(Sum.per.SF = sum(Area.with.QC.mean, na.rm = TRUE)) %>%
  select(Precursor.Ion.Name, Eddy, SizeFraction, Sum.per.SF) %>%
  unique() 

ggplot(SF.ratios, aes(x = Sum.per.SF, y = Precursor.Ion.Name, fill = SizeFraction)) +
  geom_bar(position = "fill", stat = "identity") +
  ggtitle("Anticyclonic: How much of the compound is in the large and small size fraction?")

## Find size fraction ratios for each compound. Do they change?
SF.ratios.all <- Combined %>%
  ungroup() %>%
  filter(Eddy == "Anticyclonic") %>%
  mutate_at(c("Area.with.QC.mean"), replace_nonvalues) %>%
  mutate(Area.with.QC.mean = as.numeric(Area.with.QC.mean)) %>%
  group_by(Precursor.Ion.Name, SizeFraction, Binned.Group) %>%
  mutate(Sum.per.SF = sum(Area.with.QC.mean, na.rm = TRUE)) %>%
  select(Precursor.Ion.Name, Eddy, SizeFraction, Binned.Group, Sum.per.SF) %>%
  unique() 

ggplot(SF.ratios.all, aes(x = Sum.per.SF, y = Precursor.Ion.Name, fill = SizeFraction)) +
  facet_wrap(~Binned.Group) +
  geom_bar(position = "fill", stat = "identity") +
  ggtitle("Anticyclonic: Do ratios change with treatments?")
  

## Good ionizers will be more abundant. 
## Ask Sonya (for example) who makes ectoine? 
