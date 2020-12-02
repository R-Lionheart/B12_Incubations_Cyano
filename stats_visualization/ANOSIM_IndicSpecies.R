source("src/Functions.R")
source("src/biostats.R")
options(scipen = 999)

library(grid)
library(gridExtra)
library(indicspecies)
library(tidyverse) 
library(vegan)

dataset.pattern <- "IsoLagran"

## Import your datasets. This will import all datasets split by eddy and filter, non-standardized long formats,
## and those that have been pivoted wider and standardized (aka gone through the NMDS distance matrix production).
filenames <- RemoveCsv(list.files(path = "data_processed/", pattern = dataset.pattern))
filepath <- file.path("data_processed", paste(filenames, ".csv", sep = ""))
for (i in filenames) {
  filepath <- file.path("data_processed/", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE))
}

## Set information for ANOSIM test
EddyInformation <- "2_Anticyclonic_5um" # EddyNumber_EddyDirection_FilterSize. e.g. "1_Cyclonic_0.2um"
mydf <- IsoLagran2_5_notstd # Non-standardized, long format.
mydf.wide.std <- IsoLagran2_Anticyclonic_5um_wide_std # Standardized, wide format.
hasChlorophyll <- "no" # yes or no

## Create Treatment dataframe for ANOSIM analysis
Treatment <- mydf %>%
  ungroup() %>%
  select(Replicate.Name) %>%
  unique() 

if (str_detect(hasChlorophyll, regex("no", ignore_case = TRUE))) {
  Treatment <- Treatment %>%
    separate(Replicate.Name, into = c("Date", "runtype", "Supergroup", "replicate"), remove = FALSE) 
  
} else if (str_detect(hasChlorophyll, regex("yes", ignore_case = TRUE))) {
  Treatment <- Treatment %>%
    separate(Replicate.Name, into = c("Supergroup", "replicate"), remove = FALSE)
}

Treatment <- Treatment %>%
  mutate(Control.Status = ifelse(str_detect(Supergroup, "IT0"),
                                 "Incubation", ifelse(str_detect(Supergroup, "DSW"), "DeepSeaWater", 
                                                      ifelse(str_detect(Supergroup, "Control"), "Control", "Treatments")))) %>%
  mutate(Treatment.Status = ifelse(Control.Status == "Control", "Control",
                                   ifelse(Control.Status == "DeepSeaWater", "DeepSeaWater",
                                          ifelse(Control.Status == "Incubation", "TimeZero",
                                                 ifelse(str_detect(Supergroup, "DMBnoBT"), "DMBnoB12",
                                                        ifelse(str_detect(Supergroup, "WBT"), "B12",
                                                               ifelse(str_detect(Supergroup, "DMB"), "DMB", "noB12"))))))) %>%
  select(Replicate.Name, Control.Status, Treatment.Status, Supergroup)

## Transform 
# The ANOSIM statistic “R” compares the mean of ranked dissimilarities between groups 
# to the mean of ranked dissimilarities within groups. An R value close to “1.0” suggests 
# dissimilarity between groups while an R value close to “0” suggests an even distribution
# of high and low ranks within and between groups

# A Significance value less than 0.05 is generally considered to be statistically significant, 
# and means the null hypothesis can be rejected. 
# Euclidean distance vs gower distance produces very similar results, both confirming statistical
# significance.

ano.Treatment.Status = anosim(mydf.wide.std[,-1], Treatment$Treatment.Status, 
                              distance = "euclidean", permutations = 9999)
ano.Treatment.Status

ano.Control.Status = anosim(mydf.wide.std[,-1], Treatment$Control.Status, 
                            distance = "euclidean", permutations = 9999)
ano.Control.Status


# Indicator Species Analysis ----------------------------------------------
# Indicator species are “A species whose status provides information on the 
# overall condition of the ecosystem and of other species in that ecosystem. 
# They reflect the quality and changes in environmental conditions as well 
# as aspects of community composition.”

# You can use the same groups used for the ANOSIM test. 
# Then you can check to see which species are most responsible for the 
# differences in microbial community composition between groups.

which.treatment <- Treatment$Treatment.Status

# The mulitpatt command results in lists of species that are associated
# with a particular group of samples. In our case, which compounds are associated with which treatments.
# If your group has more than 2 categories, multipatt will also identify species that are 
# statistically more abundant in combinations of categories.

species.indication = multipatt(mydf.wide.std[,-1], which.treatment,
                               func = "r.g", control = how(nperm=9999))
summary(species.indication)

# The first list contains the species found significantly more often in the “DSW” grouping. 
# The #sps shows the number of species that were identified as indicators for this group. 
# The first column contains species names, the next column contains the stat value 
# (higher means the OTU is more strongly associated). The p.value column contains 
# the statistical p values for the species association (lower means stronger significance). 
# The final column shows the significance level, which is explained by the Signif. codes
# at the bottom of the output.

individual.species <- as.data.frame(species.indication$sign) %>%
  rownames_to_column() %>%
  mutate(Significant = ifelse(p.value <= 0.05, "Significant", "NotSignificant")) %>%
  rename(Mass.Feature = rowname) %>%
  pivot_longer(cols = s.B12:s.TimeZero, names_to = "SampID") %>%
  mutate(SampID = substr(SampID, 3, nchar(SampID)))

individual.graph <- individual.species %>%
  group_by(Mass.Feature) %>%
  mutate(GroupName = paste(SampID[value != 0], collapse = "_")) %>%
  select(-SampID, -value) %>%
  unique() 


tbd.layout <- as.data.frame(species.indication$comb) %>%
  cbind(mydf.wide.std[1, ])

heatmap.data <- individual.graph %>%
  filter(Significant == "Significant")
heatmap <- ggplot(data = heatmap.data, aes(x = Mass.Feature, y = GroupName, fill = stat)) + 
  geom_tile(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 10)) +
  ggtitle("Species Indicator Analysis")
print(heatmap)


# Table layout
d1 <- individual.graph[1:3, 1, drop=FALSE]
d2 <- individual.graph[1:2,1:2]

g1 <- tableGrob(d1)
g2 <- tableGrob(d2)

haligned <- gtable_combine(g1,g2, along=1)
valigned <- gtable_combine(g1,g2, along=2)
grid.newpage()
grid.arrange(haligned, valigned, ncol=2)

# Heatmap layout
ggplot(heatmap.data, aes(Significant, stat, label = Mass.Feature)) +
  geom_point(mapping = aes(color = Significant)) +
  geom_text() + 
  facet_wrap(~GroupName) +
  theme(axis.text.x = element_blank())


# Boxplot layout (shitty)
gg = ggplot(individual.species, aes(x = Mass.Feature, y = stat, fill = SampID)) + 
  geom_boxplot(colour = "black", position = position_dodge(0.5)) +
  theme(legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), legend.position = "right", 
        axis.text.x = element_text(colour = "black", size = 12,
                                   angle = 90), 
        axis.text.y = element_text(size = 12, colour = "black"), 
        axis.title.y = element_text(size = 14, colour = "black"), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        legend.key=element_blank()) 
gg

# Point layout
ggplot(individual.species, aes(x=stat, y=p.value, size = value, color = SampID)) +
  geom_point(alpha=0.7) +
  ggtitle("Species Indicator Analysis")

# Polar point layout
ggplot(individual.species, aes(x=p.value, y=stat)) +
  geom_point(mapping = aes(color = value), alpha = 1/5) + 
  scale_color_gradient(low="blue", high="orange") +
  coord_polar() +
  facet_wrap(~SampID)