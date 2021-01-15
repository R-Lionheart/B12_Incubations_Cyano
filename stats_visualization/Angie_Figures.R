library(gridExtra)
library(patchwork)
library(tidyverse)

options(scipen = 999)
options(digits = 5)
source("src/Functions.R")

## Recreating Angie's graphs for the Cyano full dataset w/ Vitamins

# The mass features that are significantly different between treatments in 
# all SF, Eddy variations (four total).

# y axis is peak, x axis is treatment, facet wrap is eddy/sf.

Vitamins.Complete.Set <-  read.csv("data_processed/Vitamins_Incubations_CompleteDataset.csv") %>% 
  select(Precursor.Ion.Name, Area.with.QC, Binned.Group) %>%
  mutate(Binned.Group = factor(Binned.Group, ordered = TRUE)) %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(CountVals = sum(!is.na(Area.with.QC))) %>%
  filter(CountVals > 2) %>%
  ungroup() %>%
  separate(Binned.Group, c("SampID", "B", "C")) %>%
  unite("Grouping.ID", B:C)

# glimpse(Vitamins.Complete.Set)
# levels(Vitamins.Complete.Set$Binned.Group)
# 
# VCS.Anova <- lapply(split(Vitamins.Complete.Set, Vitamins.Complete.Set$Precursor.Ion.Name), function(i) {
#   aov(lm(Area.with.QC ~ Binned.Group, data = i))
# })
# VCS.Summary <- lapply(VCS.Anova, function(i) {
#   summary(i)
# })
# VCS.AnovaDF <- as.data.frame(do.call(rbind, lapply(VCS.Summary, function(x) {temp <- unlist(x)})))
# 
# ANOVAByGroup <- function(df, GroupID) {
#   df.split <- df %>%
#     filter(Grouping.ID == GroupID) %>%
#     mutate(SampID = factor(SampID, ordered = TRUE)) %>%
#     group_by(Precursor.Ion.Name) %>%
#     mutate(CountVals = sum(!is.na(Area.with.QC))) %>%
#     filter(CountVals > 6) %>%
#     ungroup()
#   
#   df.anova <- lapply(split(df.split, df.split$Precursor.Ion.Name), function(i) {
#     aov(lm(Area.with.QC ~ SampID, data = i))
#   })
#   df.anova.summary <- lapply(df.anova, function(i) {
#     summary(i)
#   })
# 
#   # Summarize ANOVA and create dataframe of significance
#   AnovaDF <- as.data.frame(do.call(rbind, lapply(df.anova.summary, function(x) {temp <- unlist(x)})))
#   colnames(AnovaDF)[9] <- "AnovaP"
#   AnovaDF$AnovaQ <- p.adjust(AnovaDF$AnovaP, method = "fdr")
# 
#   AnovaDF <- AnovaDF %>%
#     rownames_to_column(var = "Mass.Feature") %>%
#     mutate(AnovaSig = ifelse(AnovaQ < 0.1, TRUE, FALSE)) %>%
#     select(Mass.Feature, AnovaP, AnovaQ, AnovaSig) %>%
#     arrange(Mass.Feature)
#   
#   return(AnovaDF)
# }
# 
# Cyc_Small <- ANOVAByGroup(Vitamins.Complete.Set, "SmallFilter_Cyclonic") %>%
#   filter(AnovaSig == TRUE) 
# Cyc_Large <- ANOVAByGroup(Vitamins.Complete.Set, "LargeFilter_Cyclonic") %>%
#   filter(AnovaSig == TRUE) 
# Anti_Small <- ANOVAByGroup(Vitamins.Complete.Set, "SmallFilter_Anticyclonic") %>%
#   filter(AnovaSig == TRUE) 
# Anti_Large <- ANOVAByGroup(Vitamins.Complete.Set, "LargeFilter_Anticyclonic") %>%
#   filter(AnovaSig == TRUE) 
# 
# common <- Reduce(intersect, list(Cyc_Large$Mass.Feature, Cyc_Small$Mass.Feature,
#                                  Anti_Large$Mass.Feature, Anti_Small$Mass.Feature))
# 
# Tryptophan <- Vitamins.Complete.Set %>%
#   filter(Precursor.Ion.Name %in% common)
# 
# 
# trypto <- ggplot(Tryptophan, aes(x = SampID, y = Area.with.QC, fill = SampID)) +
#   facet_wrap(~ Grouping.ID, scales = "free") +
#   geom_boxplot() +
#   scale_fill_grey() +
#   ggtitle("Tryptophan")
# print(trypto)

PlotAllCompounds <- function (df) {
  all.plots <- ggplot(df, aes(x = SampID, y = Area.with.QC, fill = SampID)) +
    facet_wrap(~ Grouping.ID, scales = "free") +
    geom_boxplot() +
    scale_fill_grey() +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(paste(unique(df[[1]])))
  
  return(all.plots)
}

to.plot <- Vitamins.Complete.Set %>%
  group_by(Precursor.Ion.Name) %>%
  group_split()

plotlist <- lapply(to.plot, PlotAllCompounds)

pdf("~/Downloads/B12Incubations_Cyano_AllCompounds.pdf")
lapply(plotlist, print)
dev.off()