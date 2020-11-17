library(cluster)
library(Rtsne)
library(tidyverse)
options(scipen = 999)

source("src/B12_Functions.R")

# Import BMIS files
dataset.pattern <- "Time0"

filename <- RemoveCsv(list.files(path = "data_processed/", pattern = dataset.pattern))
filepath <- file.path("data_processed", paste(filename, ".csv", sep = ""))

BMISd <- assign(make.names(filename), read.csv(filepath, stringsAsFactors = FALSE, check.names = FALSE))

# Combine BMIS and standards. 
BMISd.tidied <- BMISd %>%
  separate(Replicate.Name, into = c("one", "two", "SampID", "four")) %>%
  select(Mass.Feature, SampID, Adjusted.Area) %>%
  mutate(Mass.Feature = as.factor(Mass.Feature),
         SampID = as.factor(SampID),
         Adjusted.Area = as.numeric(Adjusted.Area)) 

# Calculating Distance
# In order for a yet-to-be-chosen algorithm to group observations together, 
# we first need to define some notion of (dis)similarity between observations. 
# A popular choice for clustering is Euclidean distance. 
# However, Euclidean distance is only valid for continuous variables, 
# and thus is not applicable here. In order for a clustering algorithm to yield sensible results, 
# we have to use a distance metric that can handle mixed data types. 
# In this case, we will use something called Gower distance.

# My version, clustering mixed data types --------------------------------------------------------
full.data <- BMISd.tidied %>%
  mutate(isControl = as.factor(ifelse(str_detect(SampID, "Control"), "Control", "NonControl")),
         eddy = ifelse(str_detect(SampID, "IL1"), "Cyclonic", "Anticyclonic"),
         size = ifelse(str_detect(SampID, "5um"), "5um", "0.2um")) %>%
  unite(eddy_size, c("eddy", "size"), sep = "_") %>%
  mutate(eddy_size = as.factor(eddy_size))
glimpse(full.data)

# Gower distance
# The concept of Gower distance is actually quite simple. For each variable type, 
# a particular distance metric that works well for that type is used and scaled 
# to fall between 0 and 1. Then, a linear combination using user-specified weights 
# (most simply an average) is calculated to create the final distance matrix. 
# The metrics used for each data type are described below:
  
# quantitative (interval): range-normalized Manhattan distance
# ordinal: variable is first ranked, then Manhattan distance is used with a special adjustment
# for ties
# nominal: variables of k categories are first converted into k binary columns 
# and then the Dice coefficient is used

# pros: Intuitive to understand and straightforward to calculate
# cons: Sensitive to non-normality and outliers present in continuous variables, 
# so transformations as a pre-processing step might be necessary. 
# Also requires an NxN distance matrix to be calculated, which is computationally 
# intensive to keep in-memory for large samples

full.data <- full.data %>%
  mutate(Adjusted.Area = sqrt(Adjusted.Area))

gower.distance <- daisy(full.data, # Remove mass.feature
                      metric = "gower",
                      type = list(logratio = 3))

summary(gower.distance) # check data types: mix of N and I
gower.matrix <- as.matrix(gower.distance)


# Sanity Check
# similar (can this be expanded to include more grouping?)
full.data[
  which(gower.matrix == min(gower.matrix[gower.matrix != min(gower.matrix)]),
        arr.ind = TRUE)[1, ], ]
# dissimilar
full.data[
  which(gower.matrix == max(gower.matrix[gower.matrix != max(gower.matrix)]),
        arr.ind = TRUE)[1, ], ]

# Calculate silhouette width for many k using PAM, and plot sihouette width (higher is better)
# If you know the k-means algorithm, partioning around medoids (PAM) might look very familiar. 
# In fact, both approaches are identical, except k-means has cluster centers defined by 
# Euclidean distance (i.e., centroids), while cluster centers for PAM are restricted 
# to be the observations themselves (i.e., medoids).

# pros: Easy to understand, more robust to noise and outliers when compared to k-means, 
# and has the added benefit of having an observation serve as the exemplar for each cluster
# cons: Both run time and memory are quadratic (i.e., $O(n^2)$)

# mine
# We will use silhouette width, an internal validation metric which is an aggregated 
# measure of how similar an observation is to its own cluster compared its closest 
# neighboring cluster. The metric can range from -1 to 1, where higher values are better. 

silhouette.width <- c(NA)
for(i in 2:10){
  PAM.fit <- pam(gower.distance,
                 diss = TRUE,
                 k = i)
  silhouette.width[i] <- PAM.fit$silinfo$avg.width
}
plot(1:10, silhouette.width,
     xlab = "My Number of clusters",
     ylab = "My Silhouette Width")
lines(1:10, silhouette.width)


# mine
PAM.fit <- pam(gower.distance, diss = TRUE, k = 8)
PAM.results <- full.data %>%
  dplyr::select(-Mass.Feature) %>%
  mutate(cluster = PAM.fit$clustering) %>%
  group_by(cluster) %>%
  do(the_summary = summary(.))
PAM.results$the_summary
full.data[PAM.fit$medoids, ]

# One way to visualize many variables in a lower dimensional space is with 
# t-distributed stochastic neighborhood embedding, or t-SNE. 
# This method is a dimension reduction technique that tries to preserve local structure 
# so as to make clusters visible in a 2D or 3D visualization. 
# While it typically utilizes Euclidean distance, it has the ability to handle 
# a custom distance metric like the one we created above. In this case, 
# the plot shows the three well-separated clusters that PAM was able to detect. 

TSNE.obj <- Rtsne(gower.distance, is_distance = TRUE)
TSNE.data <- TSNE.obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(PAM.fit$clustering),
         name = full.data$Mass.Feature,
         SampID = full.data$SampID) %>%
  unique()

ggplot(aes(x = X, y = Y), data = TSNE.data) +
  geom_point(aes(color = cluster)) +
  ggtitle("PAM K-means")

TSNE.data.filter <- TSNE.data %>%
  filter(X > -30 & X < 0,
         Y > 30)
