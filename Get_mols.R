# # Information
# # This code retrieves mol/L from peak areas of targeted compounds.
# source("B12_Inc_Functions.R")
# 
# 
# # Get information from standards that can be quantified.
# 
# # Import standards, filter NAs --------------------------------------------
# Ingalls.Standards <- read.csv("data_extras/Ingalls_Lab_Standards.csv", stringsAsFactors = FALSE) %>%
#   filter(Column == "HILIC") %>%
#   select(Compound.Name, QE.RF.ratio, HILICMix, Conc..uM, Emperical.Formula) %>%
#   filter(!is.na(Conc..uM)) 
# 
# 
# # Import QC'd files and clean parameter data.
# filename <- RemoveCsv(list.files(path = 'data_processed/', pattern = 'QC'))
# filepath <- file.path('data_processed', paste(filename, ".csv", sep = ""))
# 
# HILICS <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, header = TRUE)) %>%
#   slice(-1:-6) %>%
#   select(-c(Description, Value)) %>%
#   select(Replicate.Name, Metabolite.name, Column, Area.with.QC, Area.Value, Run.Type) 
# 
# 
# #********************************************
# #  *Manual removal of duplicate compounds*
# #  * THIS NEEDS TO BE FIXED!!!*
# #  Currently, they are just being removed. 
# HILICS.duplicates <- IdentifyDuplicates(HILICS)
# 
# HILICS <- HILICS %>%
#   filter(!Metabolite.name %in% HILICS.duplicates$Metabolite.name)
# #********************************************
# 
# 
# # Filter for compounds detected in Ingalls Standards.
# HILICS.raw <- HILICS %>%
#   filter(Metabolite.name %in% Ingalls.Standards$Compound.Name) %>%
#   rename(Compound.Name = Metabolite.name)
# 
# # Isolate standards.
# HILICS <- HILICS %>%
#   filter(Metabolite.name %in% Ingalls.Standards$Compound.Name) %>%
#   filter(str_detect(Replicate.Name, "Std")) %>%
#   rename(Compound.Name = Metabolite.name) %>%
#   left_join(Ingalls.Standards, by = "Compound.Name") %>%
#   select(Replicate.Name, Compound.Name, everything()) %>%
#   unique()
# 
# # Get response factors for compounds using standards in water.
# response.factors <- HILICS %>%
#   mutate(RF = as.numeric(as.character(Area.with.QC))/as.numeric(Conc..uM)) %>%
#   mutate(Type = ifelse(str_detect(Replicate.Name, "H20"), "Standards_Water", "Standards_Matrix")) %>%
#   filter(!str_detect(Compound.Name, ",")) %>%
#   mutate(Replicate.Name = gsub("(.*)_.*", "\\1", Replicate.Name))
# 
# # Calculate RF max and min using only standards in water.
# RFs.dimensions <- response.factors %>%
#   filter(Type == "Standards_Water") %>%
#   group_by(Compound.Name) %>%
#   mutate(RF.max = max(RF, na.rm = TRUE),
#          RF.min = min(RF, na.rm = TRUE))
# 
# RFs.dimensions$RF.max[is.infinite(RFs.dimensions$RF.max) | is.nan(RFs.dimensions$RF.max) ] <- NA
# RFs.dimensions$RF.min[is.infinite(RFs.dimensions$RF.min) | is.nan(RFs.dimensions$RF.min) ] <- NA
# 
# RFs.dimensions <- RFs.dimensions %>%
#   mutate(RF.diff = RF.max/RF.min) %>%
#   unique()
# 
# # Calculate the response factor ratios using (Standards in Matrix - Water in Matrix) / (Standards in Water) for each replicate.
# RF.ratios <- response.factors %>%
#   group_by(Compound.Name, Type) %>%
#   mutate(RF.mean.per_sampleID = mean(RF, na.rm = TRUE)) %>%
#   select(Replicate.Name, Compound.Name, Type, RF.mean.per_sampleID) %>%
#   unique() %>%
#   group_by(Compound.Name) %>% filter(n() >= 3) %>%
#   mutate(RF.ratio = 
#            ((RF.mean.per_sampleID[Type == "Standards_Matrix"] - RF.mean.per_sampleID[Type == "Water_Matrix"]) / 
#               RF.mean.per_sampleID[Type == "Standards_Water"])) %>%
#   select(Compound.Name, RF.ratio) %>%
#   unique()
# 
# RF.ratios$RF.ratio[RF.ratios$Compound.Name == "Choline"] <- 1
# RF.ratios$RF.ratio[RF.ratios$Compound.Name == "Trimethyl-L-lysine"] <- NA ## DROPPED DUE TO SUPER WEIRD NUMBERS
# RF.ratios$RF.ratio[is.nan(RF.ratios$RF.ratio)] <- NA
# 
# rm(list = c("RFs.transect", "RFs2.transect"))
# 
# 
# # Replace NA RF ratios with QE RF ratios from Ingalls Standards.
# 
# test.RFratios <- transect.RFratios %>%
#   filter(is.na(RF.ratio))
# 
# test.standards <- Ingalls.Standards %>%
#   filter(Compound.Name %in% test.RFratios$Compound.Name) %>%
#   rename(RF.ratio = QE.RF.ratio) %>%
#   select(Compound.Name, RF.ratio) %>%
#   mutate(Compound.Name = as.character(Compound.Name)) %>%
#   mutate(RF.ratio = as.character(RF.ratio))
# 
# transect.RFratios <- transect.RFratios %>%
#   as.data.frame() %>%
#   filter(!is.na(RF.ratio)) %>%
#   rbind(test.standards)
# 
# 
# # Switch to BMIS'd transect dataset.
# 
# Samp.Data.transect <- read.csv("data_processed/BMIS_Output.csv") %>%
#   separate(Run.Cmpd, sep = " ", into = c("Sample.Name")) %>%
#   mutate(Compound.Name = MassFeature) %>%
#   filter(Compound.Name %in% transect.RFratios$Compound.Name) %>%
#   left_join(transect.RFratios) %>%
#   left_join(transect.RFs.dimensions %>% select(Compound.Name, Column, RF.max, RF.min) %>% unique(), by = "Compound.Name") %>%
#   mutate(RF.ratio = as.numeric(RF.ratio)) %>%
#   select(Compound.Name, FinalBMIS, Sample.Name, Adjusted_Area, everything())  
# 
# 
# 
# # *Volume of seawater filtered for all TRANSECT samples was 5 L*
# 
# # Calculate umol in vial for transect compounds without an internal standard.
# 
# Quan.Dat.transect <- Samp.Data.transect %>%
#   mutate(RF.ave = as.numeric(rowMeans(Samp.Data.transect[, c("RF.min", "RF.max")]))) %>%
#   mutate(umol.in.vial.ave = Adjusted_Area/RF.ave/RF.ratio,
#          umol.in.vial.max = Adjusted_Area/RF.min/RF.ratio,
#          umol.in.vial.min = Adjusted_Area/RF.max/RF.ratio) %>%
#   select(Compound.Name:Adjusted_Area, everything())
# 
# 
# # Pull out data for eddy center Internal Standard and matched ones.
# 
# original.IS.key <- read.csv("data_extras/InternalStandardNames.csv") %>%
#   rename(FinalBMIS = Internal_Standards) %>%
#   mutate(FinalBMIS = as.character(FinalBMIS))
# 
# IS.key <- Samp.Data.transect %>%
#   select(FinalBMIS, MassFeature) %>%
#   rename(Compound.Name = MassFeature) %>%
#   unique() %>%
#   left_join(original.IS.key %>% select(FinalBMIS, Concentration_nM)) %>%
#   filter(str_detect(FinalBMIS, as.character(Compound.Name)))
# 
# rm(original.IS.key)
# 
# 
# #*This section is a little hairy and needs to be checked and cleaned of redundancies. Also check the Isethionic acid that has that has 186 umol/vial. Another result of that low standard area.*
# 
# # Calculate umol in vial for eddycenter compounds that have a matched internal standard.
# 
# IS.data.transect <- HILICS.transect %>%
#   filter(as.character(Compound.Name) %in% IS.key$FinalBMIS) %>%
#   mutate(IS_Area = Area.with.QC,
#          FinalBMIS = Compound.Name) %>%
#   select(IS_Area, FinalBMIS, ReplicateName) %>%
#   left_join(IS.key %>% select(FinalBMIS, Concentration_nM))
# 
# IS.names <- data.frame(Compounds = c(IS.key[ ,"FinalBMIS"], as.character(IS.key[ ,"Compound.Name"])))
# 
# IS.smp.data.transect <- HILICS.raw.transect %>%
#   left_join(IS.data.transect %>% select(Compound.Name, Concentration_nM)) %>%
#   unique() %>%
#   filter(Compound.Name %in% IS.names$Compounds) %>%
#   filter(!str_detect(ReplicateName, "Std")) %>%
#   mutate(Std.Type = ifelse(str_detect(Compound.Name, ","), "Internal_std", "Standard")) %>%
#   mutate(testcol1 = ifelse(str_detect(Compound.Name, ","), sapply(strsplit(Compound.Name, ","), `[`, 1), Compound.Name)) %>%
#   mutate(Names = ifelse(str_detect(testcol1, "-"), sapply(strsplit(testcol1, "-"), `[`, 2), testcol1)) %>%
#   mutate(Pairs = ifelse(!str_detect(Compound.Name, ","), Compound.Name, paste(Names, "IS", sep = "_"))) %>%
#   select(-c("Pairs", "testcol1", "Run.Type")) %>%
#   mutate(Area.with.QC = as.numeric(as.character(Area.with.QC))) %>%
#   arrange(ReplicateName) %>%
#   group_by(Names) %>%
#   group_split()
# 
# 
# IS.mid_frame <- lapply(IS.smp.data.transect, function(x) group_by(x, ReplicateName))
# 
# IS.mid_frame2 <- lapply(IS.mid_frame,
#                         function(x)
#                           mutate(x,
#                                  umol.in.vial_IS = (Area.with.QC[Std.Type == "Standard"] / Area.with.QC[Std.Type == "Internal_std"]) * (Concentration_nM[Std.Type == "Internal_std"]/1000)))
# 
# IS.smp.data.transect <- do.call(rbind, IS.mid_frame2) %>%
#   filter(!str_detect(Compound.Name, ",")) %>%
#   rename(Sample.Name = ReplicateName) %>%
#   select(Sample.Name:Area.with.QC, Concentration_nM, umol.in.vial_IS)
# 
# rm(list = c("IS.names", "HILICS.raw.eddycenter", "IS.mid_frame", "IS.mid_frame2"))
# 
# 
# # Add IS_smp info back into main frame.
# # *Add a test to make sure the sample names are what they should be*
# all.info <- Quan.Dat.transect %>%
#   left_join(IS.smp.data.transect %>% select(Sample.Name, Compound.Name, umol.in.vial_IS)) %>%
#   mutate(umol.in.vial.ave = ifelse(is.na(umol.in.vial_IS), umol.in.vial.ave, umol.in.vial_IS),
#          umol.in.vial.max = ifelse(is.na(umol.in.vial_IS), umol.in.vial.max, NA),
#          umol.in.vial.min = ifelse(is.na(umol.in.vial_IS), umol.in.vial.min, NA)) %>%
#   rename(ReplicateName = Sample.Name) %>%
#   filter(!str_detect(ReplicateName, "DDA"))
# 
# 
# # Add in dilution factor from standard.info.
# 
# # (nmol.in.Enviro.ave = umol.in.vial.ave * conversion factor * injection volume / Volume filtered * conversion factor * Dilution.Factor)
# Dilution.Factor = 2
# quanDat2 <- all.info %>%
#   mutate(nmol.in.Enviro.ave = (umol.in.vial.ave*10^-6*400/5*1000 * Dilution.Factor)) %>% 
#   left_join(HILICS.transect %>% select(Compound.Name, Emperical.Formula)) %>%
#   select(Compound.Name, ReplicateName, Adjusted_Area, Orig_RSD:Emperical.Formula) %>%
#   unique()
# 
# 
# # Obtain numbers for Carbon and Nitrogen according to their emperical formulas
# 
# # Okay dokay, go get how many Carbons and Nitrogens there are here.
# 
# quanDat3 <- quanDat2  %>%
#   mutate(C = ifelse(is.na(str_extract(Emperical.Formula, "^C\\d\\d")),
#                     str_extract(Emperical.Formula, "^C\\d"),
#                     str_extract(Emperical.Formula, "^C\\d\\d"))) %>%
#   mutate(C = as.numeric(str_replace_all(C, "C", ""))) %>%
#   mutate(N = ifelse(str_detect(Emperical.Formula, "N\\D"),
#                     1,
#                     str_extract(Emperical.Formula, "N\\d"))) %>%
#   mutate(N = as.numeric(str_replace_all(N, "N", ""))) %>%
#   mutate(nmol.C.ave = nmol.in.Enviro.ave*C,
#          nmol.N.ave = nmol.in.Enviro.ave*N ) %>%
#   select(Compound.Name, SampID, ReplicateName, everything())
# 
# quanDatSum <- quanDat3 %>%
#   group_by(Compound.Name) %>%
#   summarise(nmol.Enviro.med = median(nmol.in.Enviro.ave, na.rm  = T),
#             nmol.Enviro.min = min(nmol.in.Enviro.ave, na.rm  = T),
#             nmol.Enviro.max = max(nmol.in.Enviro.ave, na.rm  = T),
#             nmol.C.med = median(nmol.C.ave, na.rm  = T),
#             nmol.C.min = min(nmol.C.ave, na.rm  = T),
#             nmol.C.max = max(nmol.C.ave, na.rm  = T)) %>%
#   arrange(desc(nmol.Enviro.med))
# 
# 
# # Organize final outputs and download.
# 
# # Cacluate mole fractions of each compound.
# TotalMoles <- quanDat3  %>%
#   select(SampID, nmol.C.ave, nmol.N.ave) %>%
#   group_by(SampID) %>%
#   summarise(totalCmeasured_nM_perID = sum(as.numeric(nmol.C.ave), na.rm = TRUE),
#             totalNmeasured_nM_perID = sum(as.numeric(nmol.N.ave), na.rm = TRUE))
# 
# quanDat4 <- quanDat3 %>%
#   unique() %>%
#   left_join(TotalMoles) %>%
#   mutate(ratioCN = totalCmeasured_nM_perID / totalNmeasured_nM_perID) %>%
#   mutate(molFractionC = nmol.C.ave/totalCmeasured_nM_perID,
#          molFractionN = nmol.N.ave/totalNmeasured_nM_perID) %>%
#   select(Compound.Name, ReplicateName, Adjusted_Area, Area.with.QC, RF.ratio:molFractionN) %>%
#   unique()
# 
# quanDatSum <- quanDat4 %>%
#   group_by(Compound.Name) %>%
#   summarise(nmol.Enviro.med = median(nmol.in.Enviro.ave, na.rm  = T),
#             nmol.Enviro.min = min(nmol.in.Enviro.ave, na.rm  = T),
#             nmol.Enviro.max = max(nmol.in.Enviro.ave, na.rm  = T),
#             nmol.C.med = median(nmol.C.ave, na.rm  = T),
#             nmol.C.min = min(nmol.C.ave, na.rm  = T),
#             nmol.C.max = max(nmol.C.ave, na.rm  = T),
#             mol.C.Fraction.med = median(molFractionC, na.rm = T),
#             mol.C.Fraction.min = min(molFractionC, na.rm = T),
#             mol.C.Fraction.max = max(molFractionC, na.rm = T)) %>%
#   arrange(desc(Compound.Name))
# 
# 
# write.csv(quanDatSum, "data_processed/Quantified_Summary.csv")
# write.csv(quanDat4, "data_processed/Quantified_Full.csv")
# write.csv(TotalMoles, "data_processed/Quantified_per_SampID.csv")
