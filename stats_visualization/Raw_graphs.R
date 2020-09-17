library(ggrepel)
library(patchwork)
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
assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, check.names = TRUE))

QCd <- Skyline_QE_QC_Output_Vitamins_2020.09.10 %>%
  slice(-1:-6) %>%
  select(-c(Description, Value)) %>%
  select(Replicate.Name, Precursor.Ion.Name, Area, Area.with.QC) %>%
  filter(!str_detect(Replicate.Name, "Blk|Std")) %>%
  mutate(Replicate.Name = Replicate.Name %>%
           str_replace("-",".")) %>%
  filter(!str_detect(Replicate.Name, 
                     "TruePoo|CultureMED4|Sept29QC|TruePooWeek1|TruePooWeek2|TruePooWeek3|TruePooWeek4|DSW700m|Process|Std")) %>%
  separate(Replicate.Name, into = c("Date", "runtype", "Supergroup", "replicate"), remove = FALSE) %>%
  group_by(Precursor.Ion.Name, Supergroup) %>%
  mutate(Area.mean = mean(Area, na.rm = TRUE)) %>%
  mutate(Area.with.QC.mean = mean(Area.with.QC, na.rm = TRUE)) %>%
  select(Supergroup, Precursor.Ion.Name, Area.mean, Area.with.QC.mean) %>%
  unique() %>%
  mutate_at(c("Area.with.QC.mean"), replace_nonvalues) %>%
  mutate(Area.with.QC.mean = as.numeric(Area.with.QC.mean)) %>%
  mutate(Control.Status = ifelse(str_detect(Supergroup, "T0"),
                                 "Incubation", ifelse(str_detect(Supergroup, "DSW"), "DeepSeaWater", 
                                                      ifelse(str_detect(Supergroup, "Control"), "Control", "Treatments")))) %>%
  mutate(Treatment.Status = ifelse(Control.Status == "Control", "Control",
                                   ifelse(Control.Status == "DeepSeaWater", "DeepSeaWater",
                                          ifelse(Control.Status == "Incubation", "TimeZero",
                                                 ifelse(str_detect(Supergroup, "DMBnoBT"), "DMBnoB12",
                                                        ifelse(str_detect(Supergroup, "WBT"), "B12",
                                                               ifelse(str_detect(Supergroup, "DMB"), "DMB", "noB12"))))))) %>%
  mutate(SizeFraction = ifelse(str_detect(Supergroup, "5um"), "5um", "0.2um")) %>%
  mutate(Eddy = ifelse(str_detect(Supergroup, "IL1"), "Cyclonic", "Anticyclonic")) %>%
  select(Supergroup, Precursor.Ion.Name, Area.mean, Area.with.QC.mean, Treatment.Status, SizeFraction, Eddy)


QCd$Supergroup <- factor(QCd$Supergroup, 
                         levels = c("IL1T0", "IL1Control", "IL1DMB", "IL1WBT", "IL1DSW", "IL1DMBnoBT", "IL1noBT",
                                    "IL1T05um", "IL1Control5um", "IL1DMB5um", "IL1WBT5um", "IL1DSW5um","IL1DMBnoBT5um", "IL1noBT5um",
                                    "IL2T0", "IL2Control", "IL2DMB", "IL2WBT", "IL2DSW", "IL2DMBnoBT", "IL2noBT",
                                    "IL2T05um", "IL2Control5um", "IL2DMB5um", "IL2WBT5um", "IL2DSW5um","IL2DMBnoBT5um", "IL2noBT5um"))


all.plot <- ggplot(QCd, aes(x=Supergroup, y=Area.mean, fill = Treatment.Status)) +
  facet_wrap(~Precursor.Ion.Name, scales = "free") +
  geom_bar(stat = "identity") +
  ggtitle("Vitamins: B12 Incubation, all peaks, averaged") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
print(all.plot + mytheme)


QCd.plot <- ggplot(QCd, aes(x=Supergroup, y=Area.with.QC.mean, fill = Treatment.Status)) +
  facet_wrap(~Precursor.Ion.Name, scales = "free") +
  geom_bar(stat = "identity") +
  ggtitle("Vitamins: B12 Incubation, QC'd, averaged") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) 
print(QCd.plot + mytheme)

NonQCd <- QCd %>%
  filter(is.na(Area.with.QC.mean)) 

NonQCd.plot <- ggplot(NonQCd, aes(x=Supergroup, y=Area.mean, fill=Treatment.Status)) +
  facet_wrap(~Precursor.Ion.Name) +
  geom_bar(stat = "identity") +
  ggtitle("Vitamins: B12 Incubation, ONLY nonQC PEAKS, averaged") +
  theme(axis.text.x = element_text(angle = 90)) 
print(NonQCd.plot + mytheme)

print(paste("Total number of peaks:", length(QCd$Area)))
print(paste("Number of NAs thrown out due to area minimum (less than 5000):", sum(is.na(QCd$Area.with.QC))))


# DMB section -------------------------------------------------------------

### CYCLONIC 

DMB_B12_0.2 <- QCd %>%
  filter(Precursor.Ion.Name %in% c("Ado-B12", "CN-B12", "CN-pB12", "Me-B12",
                                   "Me-pB12", "OH-B12", "OH-pB12")) %>%
  mutate(HasDMB = ifelse(str_detect(Supergroup, "DMB"), "HasDMB", "noDMB")) %>%
  filter(str_detect(Supergroup, "IL1")) %>%
  filter(!str_detect(Supergroup, "5um"))

dmbplot_0.2 <- ggplot(DMB_B12_0.2, aes(x=Precursor.Ion.Name, 
                                       y=round(Area.with.QC.mean,digits = 0),
                                       fill = HasDMB)) +
  facet_wrap(~Supergroup) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_text(aes(label=Area.with.QC.mean),
            position=position_dodge(width=0.9), vjust=-0.25) +
  ggtitle("B12 Incubations, Cyclonic 0.2um: With and Without Added DMB")
print(dmbplot_0.2)

DMB_B12_5 <- QCd %>%
  filter(Precursor.Ion.Name %in% c("Ado-B12", "CN-B12", "CN-pB12", "Me-B12",
                                   "Me-pB12", "OH-B12", "OH-pB12")) %>%
  mutate(HasDMB = ifelse(str_detect(Supergroup, "DMB"), "HasDMB", "noDMB")) %>%
  filter(str_detect(Supergroup, "IL1")) %>%
  filter(str_detect(Supergroup, "5um"))

dmbplot_5 <- ggplot(DMB_B12_5, aes(x=Precursor.Ion.Name, 
                                   y=round(Area.with.QC.mean,digits = 0),
                                   fill = HasDMB)) +
  facet_wrap(~Supergroup) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_text(aes(label=Area.with.QC.mean), 
            position=position_dodge(width=0.9), vjust=-0.25) +
  ggtitle("B12 Incubations, Cyclonic 5um: With and Without Added DMB")

dmbplot_0.2 | dmbplot_5

### CYCLONIC 

DMB_B12_anti_0.2 <- QCd %>%
  filter(Precursor.Ion.Name %in% c("Ado-B12", "CN-B12", "CN-pB12", "Me-B12",
                                   "Me-pB12", "OH-B12", "OH-pB12")) %>%
  mutate(HasDMB = ifelse(str_detect(Supergroup, "DMB"), "HasDMB", "noDMB")) %>%
  filter(str_detect(Supergroup, "IL2")) %>%
  filter(!str_detect(Supergroup, "5um"))

dmbplot_anti_0.2 <- ggplot(DMB_B12_anti_0.2, aes(x=Precursor.Ion.Name, 
                                                 y=round(Area.with.QC.mean, digits = 0),
                                                 fill = HasDMB)) +
  facet_wrap(~Supergroup) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_text(aes(label=Area.with.QC.mean), 
            position=position_dodge(width=0.9), vjust=-0.25) +
  ggtitle("B12 Incubations, AntiCyclonic 0.2um: With and Without Added DMB")

DMB_B12_anti_5 <- QCd %>%
  filter(Precursor.Ion.Name %in% c("Ado-B12", "CN-B12", "CN-pB12", "Me-B12",
                                   "Me-pB12", "OH-B12", "OH-pB12")) %>%
  mutate(HasDMB = ifelse(str_detect(Supergroup, "DMB"), "HasDMB", "noDMB")) %>%
  filter(str_detect(Supergroup, "IL2")) %>%
  filter(str_detect(Supergroup, "5um"))

dmbplot_anti_5 <- ggplot(DMB_B12_anti_5, aes(x=Precursor.Ion.Name, 
                                             y=round(Area.with.QC.mean, digits = 0),
                                             fill = HasDMB)) +
  facet_wrap(~Supergroup) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_text(aes(label=Area.with.QC.mean), 
            position=position_dodge(width=0.9), vjust=-0.25) +
  ggtitle("B12 Incubations, Antiyclonic 5um: With and Without Added DMB")

dmbplot_anti_0.2 | dmbplot_anti_5

# print(paste("Peak area for", unique(DMB_B12$Supergroup), ":", 
#                                     DMB_B12$Area.with.QC.mean))






