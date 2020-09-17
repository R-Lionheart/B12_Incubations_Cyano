## Answering questions from notebook

## Lower QC threshold again and incorporate S/N flags

## Combine with bigger project, but keep the QC separate

# NUTRIENTS: [DMBnoB12, noB12], both containing a f/2 spike 
# NONUTRIENTS: [WB12, WDMB], no additional nutrients added.
# DEEPSEAWATER: [DSW]
# TIME 0: [T0]
# CONTROL: [Control]

## DMB should be removed from averaged group, otherwise spikes the overall. 
# treatments that got b12 and dmb is lumped into control, reconsider. Maybe can bin together in full
## metabolites treatment.

library(gt)
library(paletteer)
library(tidyverse)
options(scipen = 999)
options(digits = 1)
source("src/Functions.R")

file.pattern <- "Vitamins"

replace_nonvalues <- function(x) (gsub(NaN, NA, x))

mytheme <- theme(panel.grid.major = element_line(colour="black", size = (0.03)),
                 panel.grid.minor = element_line(size = (0.2), colour="grey"))

# Import files --------------------------------------------------
filename <- RemoveCsv(list.files(path = "data_processed/", pattern = file.pattern))

filepath <- file.path("data_processed/", paste(filename,".csv", sep = ""))
assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, check.names = TRUE) %>%
  slice(-1:-9) %>%
  select(-c(Description, Value)) %>%
  filter(!str_detect(Replicate.Name, 
                     "Blk|Std|TruePoo|CultureMED4|Sept29QC|TruePooWeek1|TruePooWeek2|
                     TruePooWeek3|TruePooWeek4|DSW700m|Process")) %>%
  separate(Replicate.Name, into = c("Date", "runtype", "Supergroup", "replicate"), remove = FALSE) %>%
  mutate(SizeFraction = ifelse(str_detect(Supergroup, "5um"), "5um", "0.2um")) %>%
  mutate(Eddy = ifelse(str_detect(Supergroup, "IL1"), "Cyclonic", "Anticyclonic")))

Groups.Binned <- Skyline_QE_QC_Output_Vitamins_2020.09.10 %>%
  group_by(Precursor.Ion.Name, Eddy, SizeFraction) %>%
  mutate(Binned.Group = ifelse(str_detect(Supergroup, "DSW"), "DeepSeaWater",
                               ifelse(str_detect(Supergroup, "T0"), "TimeZero",
                                      ifelse(str_detect(Supergroup, "Control"), "Control",
                                        ifelse(str_detect(Supergroup, "DMBnoBT|noBT"), "Nutrients", "NoNutrients"))))) %>%
  select(Supergroup, Precursor.Ion.Name, Area, Area.with.QC, all.Flags, SizeFraction, Eddy, Binned.Group) %>%
  ungroup()


# Plot Preparation --------------------------------------------------------
Combined <- Groups.Binned %>%
  group_by(Precursor.Ion.Name, Binned.Group, SizeFraction, Eddy) %>%
  mutate(Area.mean = mean(Area, na.rm = TRUE),
         Area.with.QC.mean = mean(Area.with.QC, na.rm = TRUE)) %>%
  filter(!Precursor.Ion.Name == "B2-IS")

## Which compounds exist in which size fraction? 
## AT THE MOMENT THIS IS ONLY ANTICYCLONIC ##
Which.Cmpds <- Combined %>%
  select(Supergroup, Precursor.Ion.Name, SizeFraction, Eddy, 
         Binned.Group, Area.with.QC, Area.mean, Area.with.QC.mean) %>%
  filter(Eddy == "Anticyclonic")

Which.Cmpds %>%
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
  dplyr::group_by(SizeFraction, Precursor.Ion.Name) %>% 
  dplyr::summarize(
    NumberOfSamples = n(),
    PercentMissing = (sum(is.na(Area.with.QC))/NumberOfSamples)
  ) %>% 
  gt(rowname_col = "Supergroup") %>% 
  tab_header(title = md("Anticyclonic Eddy"),
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

## Do nutrients change with additions/size fractions?
ggplot(Combined, aes(x = reorder(Precursor.Ion.Name, -Area.with.QC.mean), y = Area.with.QC.mean,
                     fill = Binned.Group)) +
  geom_bar(position = "dodge", stat = "identity") +
  coord_flip() +
  facet_wrap(~SizeFraction) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
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
