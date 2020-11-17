library(tidyverse)
options(scipen = 999)

source("src/Functions.R")
currentDate <- Sys.Date()


# Import non-standardized files ------------------------------------------------------

Complete.Dataset <- read.csv("data_intermediate/Vitamins_Incubations_CompleteDataset.csv", stringsAsFactors = FALSE) 

# Set filtering conditions that correspond to the treatments you are comparing.
Condition1 <- "Control_SmallFilter_Cyclonic" 
Condition2 <- "DeepSeaWater_SmallFilter_Cyclonic"
myConditions <- c(Condition1, Condition2)
SigValue <- "pvalue" # alternative is "qvalue", when using fdr-corrected values.
file.pattern <- "Cyclonic_5um" # will be used as a search ID and title for graphs. 
SigNumber <- 0.1 # Pvalue cutoff
myDataFrame <- Complete.Dataset # Assign correct dataframe for analysis. Should be non-standardized.


# Create table for analysis -----------------------------------------------
Wide.myDataFrame <- myDataFrame %>%
  select(Precursor.Ion.Name, Area.with.QC, Replicate.Bin.Group) %>%
  pivot_wider(names_from = "Replicate.Bin.Group", values_from = "Area.with.QC") %>%
  select(Precursor.Ion.Name, matches(paste(myConditions, collapse="|"))) %>%
  drop_na()

mySamps <- colnames(Wide.myDataFrame)

# Add Condition1 vs Condition2 stats
myTreat1 <- mySamps[grepl(Condition1, mySamps)]
myTreat2 <- mySamps[grepl(Condition2, mySamps)]
myTreatsdf <- Wide.myDataFrame[, c(myTreat1, myTreat2)] %>%
  drop_na()


# Add a Pvalue for between the two treatments for QC
Wide.myDataFrame[, paste(Condition1, "v", Condition2, "_pvalue", sep = "")] <- apply(myTreatsdf, 1, function(x) 
{t.test(x[myTreat1], x[myTreat2])$p.value}) 
# Add a false-discovery-rate-corrected q value
Wide.myDataFrame[, paste(Condition1, "v", Condition2, "_qvalue", sep = "")] <- p.adjust(Wide.myDataFrame[, ncol(Wide.myDataFrame)], method = "fdr") 
# Calculate fold change: Condition 1 / Condition 2
Wide.myDataFrame[, paste(Condition1, "v", Condition2, "_FC", sep = "")] <- log2(rowMeans(Wide.myDataFrame[, myTreat1]) / rowMeans(Wide.myDataFrame[, myTreat2]))
# Calculate Condition 1 Average
Wide.myDataFrame[, paste(Condition1, "_Ave", sep = "")] <- rowMeans(Wide.myDataFrame[, myTreat1])  
# Calculate Condition 2 Average
Wide.myDataFrame[, paste(Condition2, "_Ave", sep = "")] <- rowMeans(Wide.myDataFrame[, myTreat2])
# Calculate complete row means
Wide.myDataFrame$AveSmp <- rowMeans(Wide.myDataFrame[, c(myTreat1, myTreat2)])
# Organize columns and assign significance
Wide.myDataFrame <- Wide.myDataFrame %>%
  select(Precursor.Ion.Name, contains(SigValue), everything()) %>%
  # mutate(Significance = ifelse(.[[2]] < SigNumber, "Significant",
  #                       ifelse(between(.[[2]], SigNumber, 0.5), "CloseSig", "NotSig")))
  mutate(Significance = ifelse(.[[2]] < SigNumber, "Significant", "NotSig"))

# Adjust fold change axis
FC_Yaxis <- Wide.myDataFrame %>%
  select(Precursor.Ion.Name, contains("FC")) %>%
  mutate(FC_Yaxis = .[[2]]) # Multiply by -1 here to reverse y axis

# Combine for plot.
dataToPlot <- Wide.myDataFrame %>%
  left_join(FC_Yaxis) %>%
  select(Precursor.Ion.Name, Significance, AveSmp, contains(Condition1), contains(Condition2), contains("FC"))

## Sanity Check for fold change ratios
sanitycheck <- myDataFrame %>%
  filter(Binned.Group == Condition1 | Binned.Group == Condition2) %>%
  group_by(Precursor.Ion.Name, Binned.Group) %>%
  mutate(myave = mean(Area.with.QC, na.rm = TRUE)) %>%
  drop_na()

# Plotting section --------------------------------------------------------

# Plot non-transformed data
ggplot(sanitycheck, aes(x = reorder(Precursor.Ion.Name, -myave), 
                        y = myave, fill = Binned.Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, size = 11, vjust = 0.5)) +
  ggtitle(paste(file.pattern, Condition1, "v", Condition2))

# Plot square root data
squareRoot <- sanitycheck %>%
  mutate(square.root = sqrt(myave))

ggplot(squareRoot, aes(x = reorder(Precursor.Ion.Name, -square.root), 
                       y = square.root, fill = Binned.Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, size = 11, vjust = 0.5)) +
  ggtitle(paste(file.pattern, Condition1, "v", Condition2, "_SquareRoot"))


# Condition1 v Condition 2 Significance
SignificancePlot <- ggplot(dataToPlot, aes(x = AveSmp, y = FC_Yaxis, fill = Significance, 
                                           label = Precursor.Ion.Name)) +
  geom_point(size = 3, shape = 21, stroke=0) +
  scale_fill_manual(values = c("grey", "royalblue")) +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_x_log10() +
  ggtitle(paste(file.pattern, Condition1, "v", Condition2)) +
  theme(plot.title = element_text(size = 15),
        legend.position="left",
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.text=element_text(size=10)) +
  labs(x="Average peak size", y=paste("Log2", Condition1, "/", Condition2, sep = "")) +
  theme(legend.position="right") +
  scale_y_continuous(limits=c(-3, 3)) +
  geom_text(data = subset(dataToPlot, Significance == "Significant"),
            hjust = "inward", nudge_x = 0.05, check_overlap = TRUE, size = 4) 
# geom_text(data = subset(dataToPlot, Significance == "CloseSig"),
#           hjust = "inward", nudge_x = 0.05, check_overlap = TRUE,
#           color = "lightskyblue3" )
SignificancePlot

figureFileName <- paste("Figure_", file.pattern, "_", 
                        Condition1, "v", Condition2, "_", 
                        currentDate, ".png", sep = "")

ggsave(filename = figureFileName, plot = SignificancePlot, path = "figures/")



