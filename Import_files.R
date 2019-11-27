# Actual script goes here

source("B12_Inc_Functions.R")

# Import all MSDial files --------------------------------------------------

filenames <- RemoveCsv(list.files(path = 'data_raw', pattern = '*.csv'))

for (i in filenames) {
  filepath <- file.path('data_raw', paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE))
}

matching.variable <- "cyano"


columns.to.drop <- c('Average.Rt.min.', 'Formula', 'Ontology', 'INCHIKEY', 
                     'SMILES', 'Isotope.tracking.parent.ID', 'Isotope.tracking.weight.number',
                     'MS1.isotopic.spectrum', 'MS.MS.spectrum', 'Average.Mz', 'Post.curation.result', 
                     'Fill..', 'Annotation.tag..VS1.0.', 'RT.matched',
                     'm.z.matched', 'MS.MS.matched', 'Manually.modified', 'Total.score', 
                     'RT.similarity', 'Dot.product', 'Reverse.dot.product', 'Fragment.presence..')

# Set header, filter unknowns ---------------------------------------

runs <- grep(matching.variable, names(.GlobalEnv), value = TRUE, ignore.case = TRUE)
runlist <- do.call("list", mget(runs))


headers.set <- lapply(names(runlist), function(x) SetHeader(runlist[[x]]))
names(headers.set) <- runs

for (df in seq_along(headers.set)) { 
  headers.set[[df]] <- headers.set[[df]] %>% filter(!Metabolite.name == "Unknown")
  headers.set[[df]] <- headers.set[[df]] %>% select(-one_of(columns.to.drop))
}

# Change variable classes -------------------------------------------------

classes.changed <- lapply(names(headers.set), function(x) ChangeClasses(headers.set[[x]]))
names(classes.changed) <- runs


list2env(classes.changed, globalenv())


# Rearrange datasets ------------------------------------------------------

SN_Cyano_B12.Incubations <- SN_Cyano_B12.Incubations %>%
  tidyr::gather(
    key = "Replicate.Name",
    value = "SN.Value",
    starts_with("X")) %>%
  select(Replicate.Name, SN.Value, everything())

RT_Cyano_B12.Incubations <- RT_Cyano_B12.Incubations %>%
  tidyr::gather(
    key = "Replicate.Name",
    value = "RT.Value",
    starts_with("X")) %>%
  select(Replicate.Name, RT.Value, everything())

Area_Cyano_B12.Incubations <- Area_Cyano_B12.Incubations %>%
  tidyr::gather(
    key = "Replicate.Name",
    value = "Area.Value",
    starts_with("X")) %>%
  select(Replicate.Name, Area.Value, everything())

Mz_Cyano_B12.Incubations <- Mz_Cyano_B12.Incubations %>%
  tidyr::gather(
    key = "Replicate.Name",
    value = "MZ.Value",
    starts_with("X")) %>%
  select(Replicate.Name, MZ.Value, everything())


# Combine to one dataset --------------------------------------------------
combined <- Area_Cyano_B12.Incubations %>%
  left_join(Mz_Cyano_B12.Incubations) %>%
  left_join(SN_Cyano_B12.Incubations) %>%
  left_join(RT_Cyano_B12.Incubations) %>%
  select(Replicate.Name, Area.Value, MZ.Value, RT.Value, SN.Value, everything()) %>%
  mutate(Metabolite.name = ifelse(str_detect(Metabolite.name, "Ingalls_"), sapply(strsplit(Metabolite.name, "_"), `[`, 2), Metabolite.name)) 

combined$Replicate.Name <- gsub("^.{0,1}", "", combined$Replicate.Name)

currentDate <- Sys.Date()
csvFileName <- paste("data_processed/MSDial_combined_", currentDate, ".csv", sep = "")

write.csv(combined, csvFileName, row.names = FALSE)

rm(list = ls())