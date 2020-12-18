source("src/Functions.R")
source("src/biostats.R")
options(scipen = 999)

library(tidyverse) 
library(vegan)

file.pattern = "_CompleteDataset"

# BMIS import and filtering of unnecessary SampIDs --------------------------------------------
filename <- RemoveCsv(list.files(path = "data_processed/", pattern = file.pattern))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))

vitamins.all <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE))

all.vitamins.mean <- vitamins.all %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(Total.Average = mean(Area.with.QC))

## Stacked graphs
top.15 <- all.vitamins.mean %>%
  arrange(desc(Total.Average)) %>%
  select(Precursor.Ion.Name) %>%
  unique() %>%
  head(15)
  

stacked.vitamins.data <- vitamins.all %>%
  mutate(Full.Total = sum(Area.with.QC, na.rm = TRUE)) %>%
  filter(Precursor.Ion.Name %in% top.15$Precursor.Ion.Name) %>%
  #filter(str_detect(Replicate.Name, Treatments)) %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "four")) %>%
  group_by(Precursor.Ion.Name, SampID) %>%
  mutate(My.Average = mean(Area.with.QC, na.rm = TRUE)) %>%
  mutate(Percent.Total = (My.Average / Full.Total)) %>%
  select(Precursor.Ion.Name, SampID, My.Average, Percent.Total) %>%
  unique()

# Stacked Cyano
stacked.vitamins.data$SampID <- factor(stacked.vitamins.data$SampID, levels = 
                                       c("IL1IT0", "IL1Control", "IL1DMB", "IL1WBT", "IL1DSW", "IL1DMBnoBT", "IL1noBT", 
                                         "IL2IT0", "IL2Control", "IL2DMB", "IL2WBT", "IL2DSW", "IL2DMBnoBT", "IL2noBT",
                                         "IL1IT05um", "IL1Control5um", "IL1DMB5um", "IL1WBT5um", "IL1DSW5um","IL1DMBnoBT5um", "IL1noBT5um",
                                         "IL2IT05um", "IL2Control5um", "IL2DMB5um", "IL2WBT5um", "IL2DSW5um","IL2DMBnoBT5um", "IL2noBT5um"))

ggplot(stacked.vitamins.data, aes(fill=Precursor.Ion.Name, y=My.Average, x=SampID)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, size = 15, vjust = .55),
        legend.text=element_text(size=15)) +
  ggtitle("Top 15 Most Abundant Compounds")


# All Cyano plotted, no filtering
all.vitamins.mean <- vitamins.all %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(Total.Average = mean(Area.with.QC, na.rm = TRUE)) %>%
  select(Precursor.Ion.Name, Total.Average) %>%
  unique()

all.hilics <- ggplot(all.vitamins.mean, aes(x = reorder(Precursor.Ion.Name, -Total.Average), 
                                          y = Total.Average)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(all.hilics)


