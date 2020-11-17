library(patchwork)
library(tidyverse)
source("src/B12_Functions.R")

createHeatmap <- function(data, mytitle, factor_levels) {
  
  # Rearrange to required format
  colnames(data)[1] <- "Replicate.Name"
  data2 <- t(data)
  colnames(data2) <- as.character(unlist(data2[1,]))
  data2 <- as.data.frame(data2[-1,])
  
  data3 <- data.frame(row.names(data2), data2, row.names = NULL)
  colnames(data3)[1] <- "Mass.Feature"
  
  for (i in c(2:ncol(data3))) {
    data3[, i] <- as.numeric(as.character(data3[, i]))
  }
  
  data3 <- data3 %>%
    mutate(Mass.Feature = as.character(Mass.Feature))
  
  # HEATMAPS for Chl AND non-Chl standardized data -------------------------------------------------------------------
  if (str_detect(file.pattern, "IsoLagran")) {
    
    data4 <- data3 %>%
      pivot_longer(cols = starts_with("X"),
                   names_to = "Replicate.Name",
                   values_to = "Area.BMISd.Normd")
    
    data4$Replicate.Name <- gsub("^.{0,1}", "", data4$Replicate.Name)
    
    heatmap.data <- data4 %>%
      separate(Replicate.Name, into = c("Date", "runtype", "SampID", "replicate")) %>%
      select(Mass.Feature, SampID, Area.BMISd.Normd) %>%
      group_by(Mass.Feature, SampID) %>%
      mutate(Average.Area.Normd = mean(Area.BMISd.Normd)) %>%
      select(Mass.Feature, SampID, Average.Area.Normd) %>%
      unique()
    
  } else if (str_detect(file.pattern, "IL")) {
    
    data4 <- data3 %>%
      pivot_longer(cols = starts_with("IL"),
                   names_to = "Replicate.Name",
                   values_to = "Area.BMISd.Normd")
    
    heatmap.data <- data4 %>%
      separate(Replicate.Name, into = c("SampID", "replicate")) %>%
      select(Mass.Feature, SampID, Area.BMISd.Normd) %>%
      group_by(Mass.Feature, SampID) %>%
      mutate(Average.Area.Normd = mean(Area.BMISd.Normd)) %>%
      select(Mass.Feature, SampID, Average.Area.Normd) %>%
      unique()
  }
  
  # Assign factor levels to SampID and generate heatmap plot
  heatmap.data$SampID <- factor(heatmap.data$SampID, levels = factor_levels)
  
  heatmap <- ggplot(data = heatmap.data, aes(x = Mass.Feature, y = SampID, fill = Average.Area.Normd)) + 
    geom_tile(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
          #axis.text.x = element_blank(),
          axis.text.y = element_text(size = 10),
          strip.text = element_text(size = 10)) +
    scale_y_discrete(limits = rev(levels(as.factor(heatmap.data$SampID)))) +
    ggtitle(paste(mytitle, Norm.Type))
  print(heatmap)
  
}
filterForOsmolytes <- function(df) {
  osmolyte.locations <- which(colnames(df) %in% osmocolumns)
  df.with.osmolytes <- df[, c(osmolyte.locations)]
  full.osmolyte.df <- cbind(df[1], df.with.osmolytes)
  
  return(full.osmolyte.df)
}

# Upload files and enter user data ----------------------------------------
file.pattern = "IsoLagran" # "IL" or Chl-normalized data. "IsoLagran" for plain BMIS-adjusted data
Norm.Type = "BMISd" # For titles & labeling. Add "with ChlA" when appropriate

filenames <- RemoveCsv(list.files(path = "data_processed/", pattern = file.pattern))
filepath <- file.path("data_processed", paste(filenames, ".csv", sep = ""))

for (i in filenames) {
  filepath <- file.path("data_processed/", paste(i, ".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE))
}

# Create heatmaps for all non-Chl datasets -----------------------
IL1_0.2_noChl <- createHeatmap(IsoLagran1_Cyclonic_0.2um_wide_std, mytitle = "Cyclonic Eddy, 0.2um",
                               c("IL1IT0", "IL1Control", "IL1DMB", "IL1WBT", "IL1DSW", "IL1DMBnoBT", "IL1noBT"))
IL1_5_noChl <- createHeatmap(IsoLagran1_Cyclonic_5um_wide_std, mytitle = "Cyclonic Eddy, 5um",
                             c("IL1IT05um", "IL1Control5um", "IL1DMB5um", "IL1WBT5um", "IL1DSW5um","IL1DMBnoBT5um", "IL1noBT5um"))
IL2_0.2_noChl <- createHeatmap(IsoLagran2_Anticyclonic_0.2um_wide_std, mytitle = "Anticyclonic Eddy, 0.2um", 
                               c("IL2IT0", "IL2Control", "IL2DMB", "IL2WBT", "IL2DSW", "IL2DMBnoBT", "IL2noBT"))
IL2_5_noChl <- createHeatmap(IsoLagran2_Anticyclonic_5um_wide_std, mytitle = "Anticyclonic Eddy, 5um",
                             c("IL2IT05um", "IL2Control5um", "IL2DMB5um", "IL2WBT5um", "IL2DSW5um","IL2DMBnoBT5um", "IL2noBT5um"))

# Create heatmaps for all Chl-normd datasets -----------------------
IL1_5_wChl <- createHeatmap(IL1_5um_ChlAnormd_std, mytitle = "Cyclonic Eddy, 5um",
                             c("IL1IT05um", "IL1Control5um", "IL1DMB5um", "IL1WBT5um", "IL1DSW5um","IL1DMBnoBT5um", "IL1noBT5um"))
IL2_5_wChl <- createHeatmap(IL2_5um_ChlAnormd_std, mytitle = "Anticyclonic Eddy, 5um",
                             c("IL2IT05um", "IL2Control5um", "IL2DMB5um", "IL2WBT5um", "IL2DSW5um","IL2DMBnoBT5um", "IL2noBT5um"))

# Create heatmaps for Osmolyte-filtered datasets -----------------------
osmocolumns <- Osmolytes_edited[[2]]

## Chl-adjusted osmolytes
IL1_5um_ChlAnormd_std_osmolytes<- filterForOsmolytes(IL1_5um_ChlAnormd_std)
IL1_5um_ChlAnormd_std_osmolytes_heatmap <- createHeatmap(IL1_5um_ChlAnormd_std_osmolytes, 
                                                         mytitle = "Cyclonic Eddy, 5um, Chl-Normalized, Osmolytes",
                      c("IL1IT05um", "IL1Control5um", "IL1DMB5um", "IL1WBT5um", "IL1DSW5um","IL1DMBnoBT5um", "IL1noBT5um"))

IL2_5um_ChlAnormd_std_osmolytes <- filterForOsmolytes(IL2_5um_ChlAnormd_std)
IL2_5um_ChlAnormd_std_osmolytes_heatmap <- createHeatmap(IL2_5um_ChlAnormd_std_osmolytes, 
                                                         mytitle = "Anticyclonic Eddy, 5um, Chl-Normalized, Osmolytes",
                                                         c("IL2IT05um", "IL2Control5um", "IL2DMB5um", "IL2WBT5um", "IL2DSW5um","IL2DMBnoBT5um", "IL2noBT5um"))


# Plots
IL1_5_noChl / IL1_5_wChl
IL2_5_noChl / IL2_5_wChl


