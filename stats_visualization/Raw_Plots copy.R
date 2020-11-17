source("src/B12_Functions.R")
source("src/biostats.R")
options(scipen = 999)

library(tidyverse) 
library(vegan)

BMIS.pattern = "Time0_Fixed"

# BMIS import and filtering of unnecessary SampIDs --------------------------------------------
filename <- RemoveCsv(list.files(path = "data_processed/", pattern = BMIS.pattern))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))

HILIC_all <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE))

all.hilics.data <- HILIC_all %>%
  group_by(Mass.Feature) %>%
  mutate(Total.Average = mean(Adjusted.Area))

## Stacked graphs
top.15 <- all.hilics.data %>%
  arrange(desc(Total.Average)) %>%
  select(Mass.Feature) %>%
  unique() %>%
  head(15)
  

stacked.hilics.data <- HILIC_all %>%
  mutate(Full.Total = sum(Adjusted.Area, na.rm = TRUE)) %>%
  filter(Mass.Feature %in% top.15$Mass.Feature) %>%
  #filter(str_detect(Replicate.Name, Treatments)) %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "four")) %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(My.Average = mean(Adjusted.Area, na.rm = TRUE)) %>%
  mutate(Percent.Total = (My.Average / Full.Total)) %>%
  select(Mass.Feature, SampID, My.Average, Percent.Total) %>%
  unique()

# Stacked HILICS
stacked.hilics.data$SampID <- factor(stacked.hilics.data$SampID, levels = 
                                       c("IL1IT0", "IL1Control", "IL1DMB", "IL1WBT", "IL1DSW", "IL1DMBnoBT", "IL1noBT", 
                                         "IL2IT0", "IL2Control", "IL2DMB", "IL2WBT", "IL2DSW", "IL2DMBnoBT", "IL2noBT",
                                         "IL1IT05um", "IL1Control5um", "IL1DMB5um", "IL1WBT5um", "IL1DSW5um","IL1DMBnoBT5um", "IL1noBT5um",
                                         "IL2IT05um", "IL2Control5um", "IL2DMB5um", "IL2WBT5um", "IL2DSW5um","IL2DMBnoBT5um", "IL2noBT5um"))

ggplot(stacked.hilics.data, aes(fill=Mass.Feature, y=My.Average, x=SampID)) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(angle = 90, size = 15, vjust = .55),
        legend.text=element_text(size=15)) +
  ggtitle("Top 15 Most Abundant Compounds")


# All HILICS plotted, no filtering
all.hilics.data <- HILIC_all %>%
  group_by(Mass.Feature) %>%
  mutate(Total.Average = mean(Adjusted.Area, na.rm = TRUE)) %>%
  select(Mass.Feature, Total.Average) %>%
  unique()

all.hilics <- ggplot(all.hilics.data, aes(x = reorder(Mass.Feature, -Total.Average), 
                                          y = Total.Average)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(all.hilics)


