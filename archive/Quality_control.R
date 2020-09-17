# Quality control script
source("B12_Inc_Functions.R")

area.min   <- 5000 # QE suggestion: HILIC - 1000, Cyano - 5000
RT.flex    <- 0.2 # QE suggestion: +/- 0.4 min for HILIC, +/- 0.2 min for Cyano
blk.thresh <- 0.2 # QE suggestion: +/- 0.2
SN.min     <- 5 # QE suggestion: 5 for Cyano, 4 for HILIC
ppm.flex   <- 7 # QE suggestion: 7

pattern = "combined"


# Import QC'd files and clean parameter data.
filename <- RemoveCsv(list.files(path = 'data_processed/', pattern = pattern))
filepath <- file.path('data_processed', paste(filename, ".csv", sep = ""))

combined <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
  select(Replicate.Name:Alignment.ID, Metabolite.name) %>%
  mutate(Run.Type = (tolower(str_extract(Replicate.Name, "(?<=_)[^_]+(?=_)")))) %>%
  rename(Metabolite.Name = Metabolite.name)

# # Quick KRH compound analysis
KRH_compounds <- read.csv("data_extras/KRH_compounds.csv", stringsAsFactors = FALSE) %>%
  select(X) %>%
  rename(Metabolite.Name = X)

for_anitra <- combined %>% filter(Metabolite.Name == "Coenzyme B12")
ggplot(for_anitra, aes(x = Replicate.Name, y = Area.Value)) +
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Coenzyme B12 in Cyano Incubation")

for_anitra_all <- combined %>% filter(Metabolite.Name %in% KRH_compounds$Metabolite.Name)
ggplot(for_anitra_all, aes(x = Replicate.Name, y = Area.Value)) +
  geom_bar(position = "dodge", stat="identity", aes(color = Metabolite.Name)) +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(Metabolite.Name ~ ., scales = "free_y") + theme(legend.position = "none") +
  ggtitle("KRH Compounds in Cyano Incubation")



## Continue analysis
msdial.runtypes <- IdentifyRunTypes(combined)

RT.table <- combined %>%
  filter(Run.Type == "std") %>%
  mutate(RT.Value = na_if(RT.Value, 0)) %>%
  arrange(Metabolite.Name) %>%
  group_by(Metabolite.Name) %>%
  mutate(RT.min = min(RT.Value, na.rm = TRUE)) %>%
  mutate(RT.max = max(RT.Value, na.rm = TRUE)) %>%
  mutate(RT.diff = abs(RT.max - RT.min)) %>%
  select(Metabolite.Name:RT.diff) %>%
  unique()

blank.table <- combined %>%
  filter(Run.Type == "blk") %>%
  mutate(Blk.Area = Area.Value) %>%
  arrange(Metabolite.Name) %>%
  group_by(Metabolite.Name) %>%
  mutate(Blk.min = min(Area.Value)) %>%
  mutate(Blk.max = max(Area.Value)) %>%
  select(Metabolite.Name:Blk.max) %>%
  select(-Blk.Area) %>%
  unique()

# Create signal to noise (SN) and area minimum flags --------------------------------
SN.Area.Flags <- combined %>%
  arrange(Metabolite.Name) %>%
  mutate(SN.Flag       = ifelse(((SN.Value) < SN.min), "SN.Flag", NA)) %>%
  mutate(Area.Min.Flag = ifelse((Area.Value < area.min), "Area.Min.Flag", NA))

# Create retention time flags ---------------------------------------
add.RT.Flag <- SN.Area.Flags %>%
  left_join(RT.table %>% select(-Run.Type)) %>%
  mutate(RT.Flag = ifelse((RT.Value >= (RT.max + RT.flex) | RT.Value <= (RT.min - RT.flex)), "RT.Flag", NA)) %>%
  select(-c("RT.max", "RT.min", "RT.diff"))

# Create blank flags ---------------------------------------
add.blk.Flag <- add.RT.Flag %>%
  left_join(blank.table %>% select(-Run.Type)) %>%
  mutate(Blank.Flag = ifelse((Area.Value / Blk.max) < blk.thresh, "Blank.Flag", NA)) %>%
  select(-c("Blk.min", "Blk.max"))

# Combine all the flags ---------------------------------------------------
final.table <- add.blk.Flag %>%
  mutate(all.Flags      = paste(SN.Flag, Area.Min.Flag, RT.Flag, Blank.Flag, sep = ", ")) %>%
  mutate(all.Flags      =all.Flags %>% str_remove_all("NA, ") %>% str_remove_all("NA")) %>%
  mutate(all.Flags      = ifelse(all.Flags == "", NA, all.Flags)) %>%
  mutate(Area.with.QC   = ifelse(is.na(Area.Min.Flag), Area.Value, NA)) %>%
  select(Replicate.Name:Area.Value, Area.with.QC, everything())


test <- final.table %>% filter(!Metabolite.Name == "Betaine")
ggplot(test, aes(x = reorder(Metabolite.Name, -Area.Value), 
                        y = Area.Value)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Pre-QC Cyano Compounds")

# Print to file with comments and a new name ------------------------------
Description <- c(as.character(anydate(Sys.Date())),
                 "Hello! Welcome to the world of MSDIAL QE Quality Control! ",
                 "Minimum area for a real peak: ",
                 "RT flexibility: ",
                 "Blank can be this fraction of a sample: ",
                 "S/N ratio: ")
Value <- as.character(c(NA, NA, area.min, RT.flex, blk.thresh, SN.min))

df <- data.frame(Description, Value)
final.table <- bind_rows(df, final.table)

currentDate <- Sys.Date()
csvFileName <- paste("data_processed/QC_Cyano_Output_", currentDate, ".csv", sep = "")

write.csv(final.table, csvFileName, row.names = FALSE)

rm(list = ls()[!ls() %in% c("final.table", lsf.str())])
