library(tidyverse)
options(scipen = 999)

source("src/B12_Functions.R")

# Setup -------------------------------------------------------------------
# Specific size/eddy

# Set treatment filters
myTreatments <- c("IL1DSW|IL1noBT|IL1Control")
myTitle <- "Anticyclonic Eddy (#1), 5um size fraction"
myTreatmentTitle <- "Deep Sea Water v no B12 v Control"
myTreatmentsSplit <- unlist(strsplit(myTreatments, split = '|', fixed = TRUE))

dataset.pattern <- "notstd|wide"

## Import your datasets. This will import a lot of information.
filenames <- RemoveCsv(list.files(path = "data_processed/", pattern = dataset.pattern))
filepath <- file.path("data_processed", paste(filenames, ".csv", sep = ""))
for (i in filenames) {
  filepath <- file.path("data_processed/", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE))
}

# BMISd_1_0.2_notnormd <- read.csv("data_processed/IsoLagran1_0.2_notnormd.csv", stringsAsFactors = FALSE) %>%
#   select(Mass.Feature:Adjusted.Area)
# BMISd_1_5_notnormd <- read.csv("data_processed/IsoLagran1_5_notnormd.csv", stringsAsFactors = FALSE) %>%
#   select(Mass.Feature:Adjusted.Area)
# BMISd_2_0.2_notnormd <- read.csv("data_processed/IsoLagran2_0.2_notnormd.csv", stringsAsFactors = FALSE) %>%
#   select(Mass.Feature:Adjusted.Area)
# BMISd_2_5_notnormd <- read.csv("data_processed/IsoLagran2_5_notnormd.csv", stringsAsFactors = FALSE) %>%
#   select(Mass.Feature:Adjusted.Area)
# 
# BMISd_1_0.2_wide <- read.csv("data_processed/IsoLagran1_0.2_normd.csv", stringsAsFactors = FALSE, check.names = FALSE)
# BMISd_1_5_wide <- read.csv("data_processed/IsoLagran1_5_normd.csv", stringsAsFactors = FALSE, check.names = FALSE)
# BMISd_2_0.2_wide <- read.csv("data_processed/IsoLagran2_0.2_normd.csv", stringsAsFactors = FALSE, check.names = FALSE)
# BMISd_2_5_wide <- read.csv("data_processed/IsoLagran2_5_normd.csv", stringsAsFactors = FALSE, check.names = FALSE)

# Assign names
BMISd.long.notstd <- IsoLagran1_0.2_notstd
BMISd.wide.std <- IsoLagran1_Cyclonic_0.2um_wide_std

positions <- c(2:4)
colnames(BMISd.wide.std)[1] <- "Replicate.Name"
names(BMISd.long.notstd) <- make.names(names(BMISd.long.notstd))
BMISd.wide.notstd <- BMISd.long.notstd %>%
  select(-X) %>%
  pivot_wider(names_from = Replicate.Name, 
              values_from = Adjusted.Area)

# Change all sets to long format
BMISd.long.std <- BMISd.wide.std %>%
  pivot_longer(-Replicate.Name,
               names_to = "Mass.Feature",
               values_to = "Area.BMISd.Normd") %>%
  select(Mass.Feature, Replicate.Name, Area.BMISd.Normd)

BMISd.wide.std <- BMISd.long.std %>%
  spread(key = "Replicate.Name", value = "Area.BMISd.Normd")


# Combine all data and rearrange
full.BMISd <- BMISd.long.std %>%
  left_join(BMISd.long.notstd) %>%
  filter(str_detect(Replicate.Name, myTreatments)) %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "four"), remove = FALSE) %>% 
  select(Mass.Feature, Replicate.Name, SampID, Area.BMISd.Normd, Adjusted.Area) %>%
  unique() %>%
  arrange(Mass.Feature)
full.BMISd <- full.BMISd[complete.cases(full.BMISd), ]

# Quick fold change test
FC.test <- full.BMISd %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(Average.Area = mean(Adjusted.Area, na.rm = TRUE)) %>%
  select(Mass.Feature, SampID, Average.Area) %>%
  filter(Average.Area > 10000000) %>% # Filter out small peaks for visibility
  unique()

FC.plot <- ggplot(FC.test, aes(x = Mass.Feature, y = Average.Area)) +
  geom_bar(position = "dodge", stat="identity", aes(fill = SampID)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ggtitle(paste("Areas between samples:", myTitle))
FC.plot 

# Start analysis ----------------------------------------------------------

# Calculate fold changes between treatments
Analysis <- BMISd.wide.notstd[complete.cases(BMISd.wide.notstd), ]
mySamps <- colnames(Analysis)

myTreat1 <- mySamps[grepl(myTreatmentsSplit[1], mySamps)] # Treat1 = DSW
myTreat2 <- mySamps[grepl(myTreatmentsSplit[2], mySamps)] # Treat2 = noB12
myTreat3 <- mySamps[grepl(myTreatmentsSplit[3], mySamps)] # Treat3 = Control
myTreatsdf <- Analysis[, c(myTreat1, myTreat2, myTreat3)]

Analysis <- Analysis %>%
  mutate(!!paste(myTreatmentsSplit[1], "v", myTreatmentsSplit[2], "_FC", sep = "") := 
           log2(rowMeans(Analysis[, myTreat1]) / rowMeans(Analysis[, myTreat2]))) %>% 
  mutate(!!paste(myTreatmentsSplit[2], "v", myTreatmentsSplit[3], "_FC", sep = "") := 
           log2(rowMeans(Analysis[, myTreat2]) / rowMeans(Analysis[, myTreat3]))) %>% 
  mutate(!!paste(myTreatmentsSplit[1], "v", myTreatmentsSplit[3], "_FC", sep = "") :=
           log2(rowMeans(Analysis[, myTreat1]) / rowMeans(Analysis[, myTreat3]))) %>%
  select(matches("Mass|FC")) %>%
  mutate(PlotA_Yaxis = .[[2]],
         PlotB_Yaxis = .[[3]],
         PlotC_Yaxis = .[[4]])

# Set up data for ANOVA
AnovaDSW <- full.BMISd %>%
  select(-Replicate.Name) %>%
  mutate(SampID = factor(SampID, ordered = TRUE)) %>%
  group_by(Mass.Feature) %>%
  mutate(Count = n()) %>%
  filter(!Count < 8) %>% # DSW so 8 instead of 9
  arrange(Mass.Feature) %>%
  select(-Count) 


# Graph normalized areas for reference
Normd.Areas <- ggplot(AnovaDSW, aes(x = SampID, y = Area.BMISd.Normd, fill = SampID)) +
  geom_boxplot() +
  facet_wrap(~Mass.Feature) +
  theme(axis.text.x = element_blank()) +
  geom_jitter(shape = 15,
              color = "steelblue",
              position = position_jitter(0.21)) +
  ggtitle(myTitle) 
Normd.Areas

# Apply ANOVA to dataframe, summarize and check significance
AnovaList <- lapply(split(AnovaDSW, AnovaDSW$Mass.Feature), function(i) {
  aov(lm(Adjusted.Area ~ SampID, data = i))
}) 
AnovaListSummary <- lapply(AnovaList, function(i) {
  summary(i)
})


# Apply Tukey to dataframe and check significance
TukeyList <- lapply(AnovaList, function(x) TukeyHSD(x))

# Summarize ANOVA and create dataframe of significance
AnovaDF <- as.data.frame(do.call(rbind, lapply(AnovaListSummary, function(x) {temp <- unlist(x)})))
colnames(AnovaDF)[9] <- "AnovaP"
AnovaDF$AnovaQ <- p.adjust(AnovaDF$AnovaP, method = "fdr")

AnovaDF <- AnovaDF %>%
  rownames_to_column(var = "Mass.Feature") %>%
  mutate(AnovaSig_P = ifelse(AnovaP < 0.1, TRUE, FALSE)) %>%
  mutate(AnovaSig_Q = ifelse(AnovaQ < 0.1, TRUE, FALSE)) %>%
  select(Mass.Feature, AnovaP, AnovaSig_P, AnovaQ, AnovaSig_Q) %>%
  arrange(Mass.Feature)

# Summarize Tukey HSD and create dataframe of significance
TukeyDF <- as.data.frame(do.call(rbind, lapply(TukeyList, function(x) {temp <- unlist(x)}))) %>%
  select(SampID10:12) %>%
  rownames_to_column("Mass.Feature") %>%
  mutate(Sig_1 = ifelse(SampID10 < 0.1, "Significant", "NotSig"),
         Sig_2 = ifelse(SampID11 < 0.1, "Significant", "NotSig"),
         Sig_3 = ifelse(SampID12 < 0.1, "Significant", "NotSig")) %>%
  rename(!!paste(myTreatmentsSplit[1], "v", myTreatmentsSplit[2], sep = "") := SampID10,
         !!paste(myTreatmentsSplit[2], "v", myTreatmentsSplit[3], sep = "") := SampID11,
         !!paste(myTreatmentsSplit[1], "v", myTreatmentsSplit[3], sep = "") := SampID12) %>%
  rename(!!paste(myTreatmentsSplit[1], "v", myTreatmentsSplit[2], "_Sig", sep = "") := Sig_1,
         !!paste(myTreatmentsSplit[2], "v", myTreatmentsSplit[3], "_Sig", sep = "") := Sig_2,
         !!paste(myTreatmentsSplit[1], "v", myTreatmentsSplit[3], "_Sig", sep = "") := Sig_3) %>%
  # mutate(SampID10_Q = p.adjust(SampID10, method = "fdr"),
  #        SampID11_Q = p.adjust(SampID11, method = "fdr"),
  #        SampID12_Q = p.adjust(SampID12, method = "fdr")) %>%
  # mutate(Sig_1Q = ifelse(SampID10_Q < 0.1, TRUE, FALSE),
  #        Sig_2Q = ifelse(SampID11_Q < 0.1, TRUE, FALSE),
  #        Sig3_Q = ifelse(SampID12_Q < 0.1, TRUE, FALSE)) %>%
  arrange(Mass.Feature)

# Join with original data to plot
toPlot <- full.BMISd %>%
  left_join(AnovaDF) %>%
  group_by(Mass.Feature, SampID) %>%
  mutate(AveSmp = mean(Area.BMISd.Normd, na.rm = TRUE)) %>%
  select(Mass.Feature, SampID, AveSmp, AnovaSig_P, AnovaSig_Q) %>%
  unique() %>%
  left_join(TukeyDF) %>% # for tukey-determined significance, p only (no fdr)
  left_join(Analysis) %>%
  arrange(Mass.Feature) %>%
  drop_na() %>%
  group_by(Mass.Feature) %>%
  mutate(TotalAve = mean(AveSmp))

# DSW v noB12
a <- ggplot(toPlot, aes_string(x = names(toPlot)[18], y = (names(toPlot)[12]), fill = names(toPlot)[5],
                        label = names(toPlot)[1])) +
  geom_point(size = 3, shape = 21, stroke=0) +  
  scale_fill_manual(values = c("grey", "royalblue4")) +
  scale_alpha_manual(values = c(1, 0.5)) +
  ggtitle(paste("DSW v noB12", myTitle)) +
  theme(plot.title = element_text(size = 15),
        legend.position="none",
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text=element_text(size=9)) +
  labs(x="Average normalized peak size", y=expression(paste(Log[2], "DSW v Control", sep = ""))) +
  theme(legend.position="right") +
  geom_text(data = subset(toPlot, toPlot[5] == TRUE),
            nudge_y = -0.06, nudge_x = 0.02, check_overlap = FALSE)
a

# noB12 v Control
b <- ggplot(toPlot, aes_string(x = names(toPlot)[18], y = (names(toPlot)[13]), fill = names(toPlot)[5],
                               label = names(toPlot)[1])) +
  geom_point(size = 3, shape = 21, stroke=0) +  
  scale_fill_manual(values = c("grey", "royalblue4")) +
  scale_alpha_manual(values = c(1, 0.5)) +
  ggtitle(paste("noB12 v Control", myTitle)) +
  theme(plot.title = element_text(size = 15),
        legend.position="none",
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text=element_text(size=9)) +
  labs(x="Average normalized peak size", y=expression(paste(Log[2], "DSW v Control", sep = ""))) +
  theme(legend.position="right") +
  geom_text(data = subset(toPlot, toPlot[5] == TRUE), 
            nudge_y = -0.06, nudge_x = 0.02, check_overlap = FALSE)
b

# DSW v Control
c <- ggplot(toPlot, aes_string(x = names(toPlot)[18], y = (names(toPlot)[14]), fill = names(toPlot)[5],
                               label = names(toPlot)[1])) +
  geom_point(size = 3, shape = 21, stroke=0) +  
  scale_fill_manual(values = c("grey", "royalblue4")) +
  scale_alpha_manual(values = c(1, 0.5)) +
  ggtitle(paste("DSW v Control", myTitle)) +
  theme(plot.title = element_text(size = 15),
        legend.position="none",
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text=element_text(size=9)) +
  labs(x="Average normalized peak size", y=expression(paste(Log[2], "DSW v Control", sep = ""))) +
  theme(legend.position="right") +
  geom_text(data = subset(toPlot, toPlot[5] == TRUE), 
            nudge_y = -0.06, nudge_x = 0.02, check_overlap = FALSE)
c
