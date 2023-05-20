###########################################################################################
#Task 1
## Load Libraries
library(tidyverse)
library(readxl)
library(factoextra)
library(tidyverse)
# Importing data set

Cancer_data <- read_xlsx("D:/breast_cancer_data.xlsx")
str(Cancer_data)

# Examine the distribution of the target variable
ggplot(Cancer_data, aes(x = diagnosis, fill = diagnosis)) +
  geom_bar() +
  ggtitle("Distribution of Diagnosis")

# Select only the columns for mean values
Cancer_data_means <- Cancer_data %>% 
  select(diagnosis, radius_mean:fractal_dimension_mean)

# Melt the dataset for plotting
Cancer_data_means_melted <- Cancer_data_means %>% 
  pivot_longer(cols = -diagnosis, names_to = "Measurement", values_to = "Value")

# Plot the distribution of different measurements
ggplot(Cancer_data_means_melted, aes(x = Value, fill = diagnosis)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Measurement, scales = "free") +
  ggtitle("Distribution of Mean Measurements")
##1
# Extract the first 5 and first 10 columns for analysis
data_pca5 <- Cancer_data[, 1:5]
data_pca10 <- Cancer_data[, 1:10]
# Convert diagnosis column to binary 0's and 1's
data_pca5$diagnosis <- ifelse(data_pca5$diagnosis == "M", 1, 0)
data_pca10$diagnosis <- ifelse(data_pca10$diagnosis == "M", 1, 0)
# Perform PCA on the first 5 columns
pca5 <- prcomp(data_pca5, center = TRUE, scale. = TRUE)
# Scree plot
plot(pca5, type = "l")
abline(h = 1, col = "red")
# Loadings plot for PC1 and PC2
biplot(pca5, scale = 0)
# Biplot using PC1 and PC2
biplot(pca5, scale = 0, choices = c(1, 2))
# Biplot using PC2 and PC3
biplot(pca5, scale = 0, choices = c(2, 3))
# Perform PCA on the first 10 columns
pca10 <- prcomp(data_pca10, center = TRUE, scale. = TRUE)
# Scree plot
plot(pca10, type = "l")
abline(h = 1, col = "red")
# Loadings plot for PC1 and PC2
biplot(pca10, scale = 0)

# Biplot using PC1 and PC2
biplot(pca10, scale = 0, choices = c(1, 2))

# Biplot using PC2 and PC3
biplot(pca10, scale = 0, choices = c(2, 3))

##2
# Load Libraries
library(tidyverse)
library(readxl)
library(cluster)
library(factoextra)
library(dendextend)

# Convert diagnosis column to binary 0's and 1's
Cancer_data$diagnosis <- ifelse(Cancer_data$diagnosis == "M", 1, 0)

# Extract the first 5 and first 10 columns for analysis
data_cluster_all <- Cancer_data[, 3:32]
data_cluster_10 <- Cancer_data[, 3:12]

# Cluster Analysis on all cell nucleus variables
# Using different distance metrics and hierarchical clustering methods
set.seed(123)
dist_methods <- c("euclidean", "manhattan", "maximum", "canberra")
hclust_methods <- c("single", "complete", "ward")
results_all <- list()

for (i in dist_methods) {
  for (j in hclust_methods) {
    hc <- hclust(dist(data_cluster_all, method = i), method = j)
    results_all <- c(results_all, list(list(dist_method = i, hclust_method = j, hc = hc)))
  }
}

# Dendrogram 1: Euclidean distance and Ward's method
dend1 <- results_all[[1]]$hc %>%
  as.dendrogram() %>%
  color_branches(k = 4)
plot(dend1)
# Dendrogram 2: Manhattan distance and Complete linkage
dend2 <- results_all[[7]]$hc %>%
  as.dendrogram() %>%
  set("branches_k_color", k = 4)
plot(dend2)

# Cluster Analysis on first 10 cell nucleus measurements
# Using different distance metrics and hierarchical clustering methods
results_10 <- list()
for (i in dist_methods) {
  for (j in hclust_methods) {
    hc <- hclust(dist(data_cluster_10, method = i), method = j)
    results_10 <- c(results_10, list(list(dist_method = i, hclust_method = j, hc = hc)))
  }
}

# Dendrogram 3: Cosine distance and Complete linkage
dend3 <- results_10[[12]]$hc %>%
  as.dendrogram() %>%
  set("branches_k_color", k = 3)
plot(dend3)

# Cluster Analysis on cell nucleus measurements (columns)
# Using the best combination of distance metric and hierarchical clustering method
set.seed(123)
best_hc <- hclust(dist(t(data_cluster_all), method = "euclidean"), method = "ward.D2")

# Dendrogram 4: Clustering of columns using Euclidean distance and Ward's method
fviz_dend(as.dendrogram(best_hc), rect = TRUE, rect_border = "red")
##################################################################################
#Task 2
# Load the required libraries
library(ggplot2)
library(GGally)
library(MASS)

# Load the dataset
data <- read_xlsx("D:/breast_cancer_data.xlsx")
# Create a new variable with diagnosis values (either M=malignant or B=benign)
data$diagnosis_value <- ifelse(data$diagnosis == "M", "Malignant", "Benign")
# Build a scatter matrix using ggpairs() and include the categorical variable
ggpairs(data[, c("radius_mean", "texture_mean", "perimeter_mean", "area_mean", "diagnosis_value")], 
        mapping = aes(color = diagnosis_value), 
        lower = list(continuous = "smooth"),
        diag = list(continuous = "density"))
# Identify strongly correlated variables
correlations <- cor(data[, c("radius_mean", "texture_mean", "perimeter_mean", "area_mean")])
print(correlations)
# The variables radius_mean, perimeter_mean and area_mean are strongly correlated with each other.
# Fit a linear model to predict the texture variable
lm_model <- lm(texture_mean ~ radius_mean + perimeter_mean + area_mean, data = data)
summary(lm_model)
# Based on the summary, the coefficients of the linear model indicate that the variable area_mean has the highest impact on texture_mean. Therefore, a single-predictor linear model with only the variable area_mean would "best" predict texture_mean.
# Identify high leverage points
leverage <- hatvalues(lm_model)
high_leverage_points <- which(leverage > (2 * ncol(data) / nrow(data)))
print(high_leverage_points)
# Remove high leverage points
data_cleaned <- data[-high_leverage_points,]
###B
# Fit Model #1 using radius and perimeter as predictors
model1 <- lm(texture_mean ~ radius_mean + perimeter_mean, data=data)
summary(model1)
cat("Model #1 AIC:", AIC(model1), "\n")
# Fit Model #2 using the best two-predictor linear model
model2 <- lm(texture_mean ~ concavity_mean + compactness_se, data = data)
summary(model2)
AIC(model2)
# Fit Model #3 using the best four-predictor linear model
model3 <- lm(texture_mean ~ concavity_mean + compactness_se + symmetry_se + smoothness_se, data = data)
summary(model3)
AIC(model3)
# Fit Model #4 using any subset of measurements
# Drop the ID column
data <- data[, -1]
# Split the data into training and test sets
set.seed(123)
train_index <- sample(nrow(data), nrow(data) * 0.7)
train_data <- data[train_index, ]
test_data <- data[-train_index, ]
# Fit the initial linear model with all predictors
model <- lm(texture_mean ~., data = train_data)
# Use stepwise selection to find the best subset of predictors
best_model <- stepAIC(model, direction = "both", trace = FALSE)
# Print the best model
summary(best_model)
# Calculate AIC for the best model
AIC(best_model)
# Create a vector of AIC scores
aic_scores <- c(AIC(model1), AIC(model2), AIC(model3), AIC(best_model))

# Create a bar plot of AIC scores
barplot(aic_scores, names.arg = c("Model #1", "Model #2", "Model #3", "Model #4"), 
        ylab = "AIC Score", main = "Comparison of AIC Scores")

##2
model1 <- lm(texture_mean ~ radius_mean + perimeter_mean, data=data)
# Residuals vs Fitted Values
plot(model1$fitted.values, model1$residuals, main="Residuals vs Fitted Values", xlab="Fitted Values", ylab="Residuals")
# Normal Q-Q Plot
qqnorm(model1$residuals, main="Normal Q-Q Plot")
qqline(model1$residuals)
# Scale-Location Plot
plot(model1$fitted.values, sqrt(abs(model1$residuals)), main="Scale-Location Plot", xlab="Fitted Values", ylab="Sqrt(|Residuals|)")
# Residuals vs Leverage
plot(model1, which=5, main="Residuals vs Leverage")
#3 
# Subset the data into two separate datasets for malignant and benign cases
malignant_data <- subset(data, diagnosis == "M")
benign_data <- subset(data, diagnosis == "B")

# Fit linear model for malignant cases
model_malignant <- lm(texture_mean ~ radius_mean + perimeter_mean + concavity_mean, data = malignant_data)
summary(model_malignant)

# Fit linear model for benign cases
model_benign <- lm(texture_mean ~ radius_mean + perimeter_mean + concavity_mean, data = benign_data)
summary(model_benign)

# Create residuals vs fitted values plot for both models
par(mfrow = c(1,2))
plot(model_malignant, which = c(1,3), main = "Malignant Cases")
plot(model_benign, which = c(1,3), main = "Benign Cases")


