library(tidyverse)
options(scipen = 999)

source("src/Functions.R")

dataset.pattern <- "IsoLagran|wide"
BMIS.pattern <- "Time0"

## Import your datasets. This will import a lot of information.
filenames <- RemoveCsv(list.files(path = "data_processed/", pattern = dataset.pattern))
filepath <- file.path("data_processed", paste(filenames, ".csv", sep = ""))
for (i in filenames) {
  filepath <- file.path("data_processed/", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE))
}


BMISd_1_0.2_long <- MakeLong(IsoLagran1_Cyclonic_0.2um_wide_std)
BMISd_1_5_long <- MakeLong(IsoLagran1_Cyclonic_5um_wide_std)
BMISd_2_0.2_long <- MakeLong(IsoLagran2_Anticyclonic_0.2um_wide_std)
BMISd_2_5_long <- MakeLong(IsoLagran2_Anticyclonic_5um_wide_std)

BMISD.normd <- BMISd_1_0.2_long %>%
  rbind(BMISd_1_5_long) %>%
  rbind(BMISd_2_0.2_long) %>%
  rbind(BMISd_2_5_long)

rm(list=ls(pattern="long|wide"))


# Complete BMISd not normd

filenames <- RemoveCsv(list.files(path = "data_processed/", pattern = BMIS.pattern))
filepath <- file.path("data_processed", paste(filenames, ".csv", sep = ""))
for (i in filenames) {
  filepath <- file.path("data_processed/", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE))
}

BMISd <- BMISd_Time0_Fixed_2020.07.15
grouped.BMISd <- BMISd %>%
  left_join(BMISD.normd) %>%
  separate(Replicate.Name, into = c("date", "runtype", "SampID", "replicate")) %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(Avg.Area = mean(Adjusted.Area)) %>%
  mutate(Count = n()) %>%
  filter(!Count < 3) %>%
  ungroup() %>%
  mutate(SampID = factor(SampID, ordered = TRUE)) 

glimpse(grouped.BMISd)
levels(grouped.BMISd$SampID)

ggplot(grouped.BMISd, aes(x = Mass.Feature, y = Avg.Area, fill = SampID)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_blank())

# Apply ANOVA to dataframe, summarize and check significance
AnovaList <- lapply(split(grouped.BMISd, grouped.BMISd$Mass.Feature), function(i) {
  aov(lm(Adjusted.Area ~ SampID, data = i))
})
AnovaListSummary <- lapply(AnovaList, function(i) {
  summary(i)
})
  
# Summarize ANOVA and create dataframe of significance
AnovaDF <- as.data.frame(do.call(rbind, lapply(AnovaListSummary, function(x) {temp <- unlist(x)})))
colnames(AnovaDF)[9] <- "AnovaP"
AnovaDF$AnovaQ <- p.adjust(AnovaDF$AnovaP, method = "fdr")

AnovaDF <- AnovaDF %>%
  rownames_to_column(var = "Mass.Feature") %>%
  mutate(AnovaSig = ifelse(AnovaQ < 0.1, TRUE, FALSE)) %>%
  select(Mass.Feature, AnovaP, AnovaQ, AnovaSig) %>%
  arrange(Mass.Feature)

TukeyList <- lapply(AnovaList, function(x) TukeyHSD(x))
TukeyDF <- as.data.frame(do.call(rbind, lapply(TukeyList, function(x) {temp <- unlist(x)}))) %>%
  #select(SampID10:12) %>%
  rownames_to_column("Mass.Feature") %>%
  arrange(Mass.Feature)

toPlot <- grouped.BMISd %>%
  left_join(AnovaDF) %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(AveSmp = mean(Area.BMISd.Normd, na.rm = TRUE)) %>%
  select(Mass.Feature, SampID, AveSmp, AnovaSig) %>%
  unique() %>%
  arrange(Mass.Feature) %>%
  drop_na() %>%
  group_by(Mass.Feature) %>%
  mutate(TotalAve = mean(AveSmp))


a <- ggplot(toPlot, aes(x = Mass.Feature, y = AveSmp, fill = AnovaSig)) +
  geom_point(size = 2, shape = 21) +  
  scale_fill_manual(values = c("grey", "royalblue4")) +
  ggtitle("oops") +
  theme(plot.title = element_text(size = 15),
        legend.position="none",
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text=element_text(size=9, angle = 90)) +
  theme(legend.position="right") 
a
