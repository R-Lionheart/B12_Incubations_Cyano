library(tidyverse)
options(scipen = 999)

source("src/Functions.R")

dataset.pattern <- "_CompleteDataset"


## Import your datasets.
filenames <- RemoveCsv(list.files(path = "data_processed/", pattern = dataset.pattern))
filepath <- file.path("data_processed", paste(filenames, ".csv", sep = ""))
for (i in filenames) {
  filepath <- file.path("data_processed", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE))
}

AnovaData.Cyano <- Vitamins_Incubations_CompleteDataset %>%
  select(Precursor.Ion.Name, Area.with.QC, Binned.Group) %>%
  mutate(Binned.Group = factor(Binned.Group, ordered = TRUE)) %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(CountVals = sum(!is.na(Area.with.QC))) %>%
  filter(CountVals > 2) %>%
  ungroup()

glimpse(AnovaData.Cyano)
levels(AnovaData.Cyano$Binned.Group)

AnovaList.Cyano <- lapply(split(AnovaData.Cyano, AnovaData.Cyano$Precursor.Ion.Name), function(i) { 
  aov(lm(Area.with.QC ~ Binned.Group, data = i))
})
AnovaListSummary.Cyano <- lapply(AnovaList.Cyano, function(i) {
  summary(i)
})

# Summarize ANOVA and create dataframe of significance
AnovaDF <- as.data.frame(do.call(rbind, lapply(AnovaListSummary.Cyano, function(x) {temp <- unlist(x)})))
colnames(AnovaDF)[9] <- "AnovaP"
AnovaDF$AnovaQ <- p.adjust(AnovaDF$AnovaP, method = "fdr")

AnovaDF <- AnovaDF %>%
  rownames_to_column(var = "Mass.Feature") %>%
  mutate(AnovaSig = ifelse(AnovaQ < 0.1, TRUE, FALSE)) %>%
  select(Mass.Feature, AnovaP, AnovaQ, AnovaSig) %>%
  arrange(Mass.Feature)

TukeyList <- lapply(AnovaList.Cyano, function(x) TukeyHSD(x))
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
