source("src/B12_Functions.R")

########################################################
## TEMP METABOLITE FILTER FOR WEI AND ANITRA 6/20/20
mystandards <- c("L-Glutamic acid", "L-Glutamine", "2-Ketoglutaric acid", "N6-Acetyl-L-lysine")
## I matched Glutamic acid and Glutamine to their standard stereoisomers: is that ok?
## There are inconsistent standard runs for some compounds. Why?
########################################################

# Import standards and filter NAs ---------------------------------------------------------------
Raw.Standards <- getURL("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv")
Ingalls.Standards <- read.csv(text = Raw.Standards, stringsAsFactors = FALSE) %>%
  filter(Column == "HILIC") %>%
  rename(Metabolite.Name = Compound.Name) %>%
  select(Metabolite.Name,Compound.Type, QE.RF.ratio, Conc..uM, HILICMix, Emperical.Formula) %>%
  filter(!is.na(Conc..uM),
        Metabolite.Name %in% mystandards)
Ingalls.Standards$Metabolite.Name <- TrimWhitespace(Ingalls.Standards$Metabolite.Name)

# Import BMIS'd sample file ---------------------------------------------------------------
filename <- RemoveCsv(list.files(path = "data_processed/", pattern = "BMIS_Output"))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))
BMISd.data <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE)) 

BMISd.data$Mass.Feature[BMISd.data$Mass.Feature == "Glutamic acid"] <- "L-Glutamic acid"
BMISd.data$Mass.Feature[BMISd.data$Mass.Feature == "Glutamine"] <- "L-Glutamine"
BMISd.data$Mass.Feature[BMISd.data$Mass.Feature == "Ketoglutaric Acid"] <- "2-Ketoglutaric acid"
BMISd.data$Mass.Feature[BMISd.data$Mass.Feature == "N(e)-Acetyl-Lysine"] <- "N6-Acetyl-L-lysine"

BMISd.data <- BMISd.data %>%
  filter(Mass.Feature %in% mystandards) 

# Import QC'd files and remove parameter data ------------------------------
filename <- RemoveCsv(list.files(path = "data_processed/", pattern = "QC_Output"))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))

QCd.data <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  slice(-1:-6) %>%
  select(-c(Description, Value, Mz.Value, RT.Value, SN.Value, Alignment.ID, contains("Flag"))) %>%
  select(Replicate.Name, Metabolite.name, Area.with.QC, Area.Value, Run.Type, everything()) 

test.glutamine <- QCd.data %>% filter(Metabolite.name == "Glutamine",
                                      Run.Type == "std")
test.keto <- QCd.data %>% filter(Metabolite.name == "Ketoglutaric Acid",
                                 Run.Type == "std")

QCd.data$Metabolite.name[QCd.data$Metabolite.name == "Glutamic acid"] <- "L-Glutamic acid"
QCd.data$Metabolite.name[QCd.data$Metabolite.name == "Glutamine"] <- "L-Glutamine"
QCd.data$Metabolite.name[QCd.data$Metabolite.name == "Ketoglutaric Acid"] <- "2-Ketoglutaric acid"
QCd.data$Metabolite.name[QCd.data$Metabolite.name == "N(e)-Acetyl-Lysine"] <- "N6-Acetyl-L-lysine"

QCd.data <- QCd.data %>%
  filter(Metabolite.name %in% mystandards)

QCd.data <- data.frame(lapply(QCd.data, function(x) {gsub("171002_Std_4uMStdinH20_1a", "171002_Std_4uMStdinH20_1", x)}),
                              stringsAsFactors = FALSE)
# Adjust for ILT0 naming issues --------------------------------------------
QCd.data <- QCd.data %>%
  mutate(Replicate.Name = recode(Replicate.Name, 
                                 "171002_Smp_IT0_1" ="171002_Smp_IL1IT0_1", 
                                 "171002_Smp_IT0_2" = "171002_Smp_IL1IT0_2",
                                 "171002_Smp_IT0_3" = "171002_Smp_IL1IT0_3",
                                 "171009_Smp_IT05um_1" = "171009_Smp_IL1IT05um_1",
                                 "171009_Smp_IT05um_2" = "171009_Smp_IL1IT05um_2",
                                 "171009_Smp_IT05um_3" = "171009_Smp_IL1IT05um_3",
                                 "171016_Smp_IT0_1" = "171016_Smp_IL2IT0_1",
                                 "171016_Smp_IT0_2" = "171016_Smp_IL2IT0_2",
                                 "171016_Smp_IT0_3" = "171016_Smp_IL2IT0_3",
                                 "171023_Smp_IT05um_1" = "171023_Smp_IL2IT05um_1",
                                 "171023_Smp_IT05um_2" = "171023_Smp_IL2IT05um_2",
                                 "171023_Smp_IT05um_3" = "171023_Smp_IL2IT05um_3"))

# Apply appropriate filters and isolate standards ---------------------------------------------------------------
QCd.w.Standards <- QCd.data %>%
  rename(Metabolite.Name = Metabolite.name) %>%
  filter(str_detect(Replicate.Name, "Std")) %>%
  left_join(Ingalls.Standards, by = "Metabolite.Name") %>%
  select(Replicate.Name, Metabolite.Name, Compound.Type, everything()) %>%
  unique()

# Get response factors for compounds ----------------------------------
Full.data.RF <- QCd.w.Standards %>%
  mutate(Area.with.QC = as.numeric(Area.with.QC),
         Area.Value = as.numeric(Area.Value),
         Conc..uM = as.numeric(Conc..uM)) %>%
  mutate(RF = Area.with.QC/Conc..uM) %>%
  mutate(Replicate.Name = substr(Replicate.Name, 1, nchar(Replicate.Name)-2))


# Calculate RF max and min using only standards in water.
Full.data.RF.dimensions <- Full.data.RF %>%
  filter(str_detect(Replicate.Name, "StdinH20")) %>%
  group_by(Metabolite.Name) %>%
  mutate(RF.max = max(RF, na.rm = TRUE),
         RF.min = min(RF, na.rm = TRUE),
         RF.diff = RF.max/RF.min) %>% 
  select(Replicate.Name, Metabolite.Name, Area.with.QC, Conc..uM, QE.RF.ratio, RF:RF.diff) %>%
  unique()

# Quantify samples for the BMIS'd dataset ---------------------------------
BMISd.data.filtered <- BMISd.data %>%
  separate(Run.Cmpd, c("Sample.Name"), extra = "drop", fill = "right") %>%
  mutate(Metabolite.Name = Mass.Feature) %>%
  left_join(Full.data.RF.dimensions %>% select(Metabolite.Name, RF.max, RF.min, QE.RF.ratio) %>% unique(), by = "Metabolite.Name") %>%
  select(Metabolite.Name, FinalBMIS, Sample.Name, Adjusted.Area, everything())

Quantitative.data <- BMISd.data.filtered %>%
  mutate(RF.ave = as.numeric(rowMeans(BMISd.data.filtered[, c("RF.min", "RF.max")]))) %>%
  mutate(QE.RF.ratio = as.numeric(QE.RF.ratio)) %>%
  mutate(umol.in.vial.ave = Adjusted.Area/RF.ave/QE.RF.ratio,
         umol.in.vial.max = Adjusted.Area/RF.min/QE.RF.ratio,
         umol.in.vial.min = Adjusted.Area/RF.max/QE.RF.ratio) %>%
  filter(type == "Smp") %>%
  select(Metabolite.Name, FinalBMIS, Sample.Name, type, SampID, replicate, Adjusted.Area, Area.with.QC, RF.max:umol.in.vial.min) %>%
  unite(Replicate.Name, c(Sample.Name:replicate))

Quantitative.data <- Quantitative.data %>%
  mutate(Replicate.Name = recode(Replicate.Name, 
                                 "171002_Smp_IT0_1" ="171002_Smp_IL1IT0_1", 
                                 "171002_Smp_IT0_2" = "171002_Smp_IL1IT0_2",
                                 "171002_Smp_IT0_3" = "171002_Smp_IL1IT0_3",
                                 "171009_Smp_IT05um_1" = "171009_Smp_IL1IT05um_1",
                                 "171009_Smp_IT05um_2" = "171009_Smp_IL1IT05um_2",
                                 "171009_Smp_IT05um_3" = "171009_Smp_IL1IT05um_3",
                                 "171016_Smp_IT0_1" = "171016_Smp_IL2IT0_1",
                                 "171016_Smp_IT0_2" = "171016_Smp_IL2IT0_2",
                                 "171016_Smp_IT0_3" = "171016_Smp_IL2IT0_3",
                                 "171023_Smp_IT05um_1" = "171023_Smp_IL2IT05um_1",
                                 "171023_Smp_IT05um_2" = "171023_Smp_IL2IT05um_2",
                                 "171023_Smp_IT05um_3" = "171023_Smp_IL2IT05um_3")) %>%
  separate(Replicate.Name, into = c("Sample.Name", "type", "SampID", "replicate"))


# Add in dilution factor and filtered volume --------------------------------------------------
All.Info.Quantitative <- Quantitative.data %>%
  #filter(str_detect(SampID, "IL1IT0|IL1IT05um|IL2IT0|IL2IT05um")) %>%
  #mutate(nmol.in.vial.ave = (umol.in.vial.ave*10^-6*Injection.Volume/Volume.Filtered*1000*Dilution.Factor)) %>%
  left_join(QCd.w.Standards %>% select(Metabolite.Name, Emperical.Formula)) %>%
  unique()

ggplot(All.Info.Quantitative, aes(x = SampID, y = umol.in.vial.ave, fill = Metabolite.Name)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x="Sample ID", 
       y="Approximate umol in vial (ave)",
       title="HILIC B12 Incubations: ESTIMATED QUANTIFICATION",
       subtitle="All sizes and eddies, unsorted") + 
  theme(plot.subtitle=element_text(size=10, color="black"),
        axis.text.x = element_text(angle = 90))


# Isolate Gln/Glu ratios to check Wei's work ------------------------------
Gln_Glu <- All.Info.Quantitative %>%
  select(Metabolite.Name, SampID, umol.in.vial.ave) %>%
  filter(str_detect(Metabolite.Name, "Glutamic|Glutamine")) %>%
  filter(str_detect(SampID, "T0")) %>%
  group_by(Metabolite.Name, SampID) %>%
  mutate(Total.umol.Ave = mean(umol.in.vial.ave, na.rm = TRUE)) %>%
  mutate(Eddy = ifelse(str_detect(SampID, "IL1"), "Cyclonic", "Anticyclonic")) %>%
  mutate(Filter.Size = ifelse(str_detect(SampID, "5um"), "5um_Filter", "0.2um_Filter")) %>%
  select(-umol.in.vial.ave) %>%
  unique() %>%
  ungroup() %>%
  group_by(SampID) %>%
  mutate(Ratios = (Total.umol.Ave[Metabolite.Name == "L-Glutamine"]) 
         / (Total.umol.Ave[Metabolite.Name == "L-Glutamic acid"]))

ggplot(data = Gln_Glu, aes(Eddy, Ratios, fill = Eddy)) +
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("B12 HILIC Incubations Cyclonic + Anticyclonic")
