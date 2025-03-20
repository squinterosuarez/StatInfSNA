######## LISS2108 Statistical Inference for Social Networks Analysis ##########
# Author: Santiago Quintero (KCL) - 2025

# In this workshop, we'll use the class's running examples and will run 
# permutation tests and autoregressive models and will compare them to
# running standard regressions. We'll create our own data to keep track of how 
# the models capture the "real" data generating process.


# Clear workspace and set seed for reproducibility
rm(list = ls())
set.seed(123)
options(scipen = 100) # turn off scientific notation

## Libraries we'll work with
packages <- c("tidyr", "lmPerm", "igraph", "spdep", "Matrix", 
              "ggplot2", "dplyr", "spatialreg")

# lmPerm  - Permutation tests in linear models 
# igraph  - Network creation 
# spdep   - Weight matrix for the autoregression
# spatialreg - Spatial Autoregressive models
# Matrix  - Matrix operations
# ggplot2 - Plotting 
# dplyr   - Data manipulation 
# tidyr   - Data manipulation

# If the packages are not installed, run:
# install.packages(packages)


# Call the libraries
lapply(packages, require, character.only = TRUE)


######## Permutation tests for group level analysis ######## 
# Running example: Does a more centralized friendship network improve the effects 
# of health interventions in a classroom? 
#
# In this case, we want to analyse if the structure of the friendship network 
# in a classroom is a moderating factor of a smoking-prevention campaing 
# (e.g., if the more “popular” central nodes changed their behaviour, it is likely 
# that other classmates did too). 
#
# For this analysis, we have a dataset of 20 classrooms, and we have data on the 
# intervention´s outcome (# of students that quit smoking). We also have data
# on the network centralisation (highest degree in the friendship network).
# We also collect data on the size of the classroom (total # of students).

# Let's generate our synthetic data
n_classrooms <- 20 # Number of classrooms
centralisation <- sample(0:15, n_classrooms, replace = TRUE) # Centralisation measure
error <- rt(n_classrooms, df = 20) * 50 #  effect with error term following a tailed distribution
smoking_reduction <- 5 + 4 * centralisation + error 

plot(error)

classroom_data <- data.frame(
  Classroom = 1:n_classrooms,
  Centralisation = centralisation,
  SmokingReduction = smoking_reduction)

### Let's first fit a standard Linear Regression (OLS) and run the usual checks to
### see the behaviour of the residuals.
lm_model <- lm(SmokingReduction ~ Centralisation, 
                data = classroom_data)
summary_lm_model <- summary(lm_model)
print(summary_lm_model)

beta_lm <- coef(summary_lm_model)["Centralisation", "Estimate"]
p_lm <- coef(summary_lm_model)["Centralisation", "Pr(>|t|)"] 
se_lm <- coef(summary_lm_model)["Centralisation", "Std. Error"]


## Usual visual inspections
# Plot residuals to check normality
ggplot(data.frame(Residuals = residuals(lm_model)), aes(x = Residuals)) +
  geom_density() +
  labs(title = "Residual Distribution of Standard Linear Regression",
       subtitle = "Heavy-tailed residuals violate normality assumption")

# And fitted vs residuals
ggplot(data.frame(Fitted = fitted(lm_model), Residuals = residuals(lm_model)), 
       aes(x = Fitted, y = Residuals)) +
  geom_point(color = "blue", size = 3, alpha = 0.7) +  # Residuals as points
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Zero reference line
  labs(title = "Residuals vs. Fitted Values",
       subtitle = "Points should be randomly scattered for a well-fitting model",
       x = "Fitted Values",
       y = "Residuals")

## Can we assume homoscedasticity? 


## Now let's run the permutation-based test manually
nperm <- 5000   # Number of permutations
perm_betas <- numeric(nperm)

for (i in 1:nperm) {
  permuted_outcome <- sample(classroom_data$SmokingReduction)  # Shuffle outcome
  perm_model <- lm(permuted_outcome ~ classroom_data$Centralisation)  # Fit permuted model
  perm_betas[i] <- coef(perm_model)[2]  # Extract beta for centralisation
}

# Compute permutation-based standard error
se_perm <- sd(perm_betas)  # Standard deviation of permuted coefficients

# Compute Permutation-based p-value
p_perm <- sum(abs(perm_betas) >= abs(beta_lm))/nperm # Two-tailed p-value

# --- Compare Results ---
comparison <- data.frame(
  Method = c("Standard Regression", "Permutation Test"),
  StdError = c(se_lm, se_perm),
  Pvalue = c(p_lm, p_perm))
print(comparison)


# Let's visualization the distribution of permuted coefficients
ggplot(data.frame(PermutedBeta = perm_betas), aes(x = PermutedBeta)) +
  geom_histogram(bins = 50, fill = "gray20", alpha = 0.6) +
  geom_vline(xintercept = beta_lm, color = "red", linetype = "dashed") +
  labs(title = "Permutation Test: Null Distribution of Coefficient",
       subtitle = "Red dashed line = observed coefficient",
       x = "Permuted Centralisation Effect",
       y = "Frequency")


## Very conveniently, there are several packages that conduct the permutation
## and calculate the significance test more efficiently
perm_lm <- lmp(SmokingReduction ~ Centralisation,
               data = classroom_data,
               perm = "Prob") # the "Prob" parameter lets the model decides the number of permutations automatically
summary_perm_lm <- summary(perm_lm)
print(summary_perm_lm) ## The specific value of the statistic test might vary because of the number of permutations



#### Permutation-test task #######
## First, run the "manual" permutation procedure varying the number of iterations
## Is there any difference?

## Now let's assume we have a much smaller sample, say 10 classrooms. Run the 
## whole analysis again using these synthetic data:
n_classrooms2 <- 10 # New number of classrooms
centralisation2 <- sample(0:15, n_classrooms2, replace = TRUE) # Centralisation measure
error2 <- rt(n_classrooms2, df = 20) * 50 #  effect with error term following a heavy tailed distribution
smoking_reduction2 <- 5 + 4 * centralisation2 + error2

# New dataset - Use this for the analysis
classroom_data2 <- data.frame(
  Classroom = 1:n_classrooms2,
  Centralisation = centralisation2,
  SmokingReduction = smoking_reduction2)

## How are the results different? Why?


## Now, we'll assume that the error is actually normally distributed.
## Run the analysis with this data
n_classrooms3 <- 10 # New number of classrooms
centralisation3 <- sample(0:15, n_classrooms3, replace = TRUE) # Centralisation measure
smoking_reduction3 <- 5 + 4 * centralisation3 + rnorm(n_classrooms, mean = 0, sd = 2)

# New dataset - Use this for the analysis
classroom_data3 <- data.frame(
  Classroom = 1:n_classrooms3,
  Centralisation = centralisation3,
  SmokingReduction = smoking_reduction3)

## What are the differences and why? Is it easier to assume homoscedasticity
## in this case?




######## Network autoregressive model for node-level analysis ######## 
# Running example: Do countries with more trading partners have higher GDP?
#
# In this case, we want to understand the impact of trade on GDP, particularly 
# the number of trading partners a country has. But we suspect that there might
# be important dependencies to account for. For instance, it might be the case
# that the GDP of the trading partners affect the impact of trading with them
# on our own GDP (countries with a higher GDP might produce more added-value 
# products).


### Simulate a trade network
n_countries <- 50  # Number of countries

# Create a random trade network - this is an Igraph network object
trade <- sample_gnp(n_countries, p = 0.1, directed = FALSE)  

## Let's plot the network for the sake of fun
plot(trade,
     vertex.label = NA,
     vertex.size = 7)

dev.off()

# Convert to adjacency matrix - We'll use it to calculate our weighting matrix
W <- as.matrix(as_adjacency_matrix(trade))
diag(W) <- 0  # Remove self-loops; In this case, we only care about *international* trade


## Let's take a look at the adjacency matrix
View(W)


# To calculate the degree of each node (the number of trading partners),
# we simply take the sum of that node's row in the adjacency matrix
trade_degree <- rowSums(W)

# We can also calculate the degree of a node apling the Igraph function degree()
# to the network object
degree(trade)


# Let's simulate an additional exogenous variable (eg, investment)
investment <- rnorm(n_countries, mean = 100000, sd = 3000)

# Base GDP
base_GDP <- 500 * investment + 5000 * trade_degree + rnorm(n_countries, mean = 10000, sd = 1000)

# Generate GDP with spatial dependence (true SAR process)
# Spatial dependence strength - this is the dependence we want to model
rho <- 0.5
w_scale <- W / rowSums(W) # scale influence matrix
GDP <- solve(diag(n_countries) - rho * w_scale) %*% base_GDP


# Create DataFrame
inter_trade <- data.frame(Country = 1:n_countries, 
                          GDP = as.vector(GDP), 
                          trade_degree, 
                          investment)

# Let's take a look at the trading partners - GDP relationship
ggplot(inter_trade, aes(x = trade_degree, y = GDP)) +
  geom_point(color = "blue", size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  labs(title = "GDP vs. Number of Trade Partners",
       x = "Number of Trade Partners",
       y = "GDP (in $)")


## Let's first run a standard OLS regression
lm_trade <- lm(GDP ~ trade_degree + investment, 
               data = inter_trade)
summary(lm_trade)

## Usual visual inspections
# Plot residuals to check normality
ggplot(data.frame(Residuals = residuals(lm_trade)), aes(x = Residuals)) +
  geom_density() +
  labs(title = "Residual Distribution of Standard Linear Regression",
       subtitle = "Heavy-tailed residuals violate normality assumption")

# And fitted vs residuals
ggplot(data.frame(Fitted = fitted(lm_trade), Residuals = residuals(lm_trade)), 
       aes(x = Fitted, y = Residuals)) +
  geom_point(color = "blue", size = 3, alpha = 0.7) +  # Residuals as points
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Zero reference line
  labs(title = "Residuals vs. Fitted Values",
       subtitle = "Points should be randomly scattered for a well-fitting model",
       x = "Fitted Values",
       y = "Residuals")


# Plot Actual GDP vs Fitted GDP for OLS model
plot(inter_trade$GDP, fitted(lm_trade), 
     main = "OLS: Actual vs Fitted GDP", 
     xlab = "Actual GDP", 
     ylab = "Fitted GDP", 
     col = "darkblue", pch = 16)
abline(0, 1, col = "darkred")  # 45-degree line (perfect fit)


## Now we run the spatial Autoregressive Model (SAR)

# First, we create spatial weights matrix (Row-standardized). We 
# need the parameter "zero.policy = TRUE" because we might have 
# "isolates" in the network, ie, countries without partners.
W_list <- mat2listw(W, style = "W", zero.policy = TRUE)


# Now, let's run the SAR Model (Using Maximum Likelihood Estimation)
sar_trade <- lagsarlm(GDP ~ trade_degree + investment, 
                      data = inter_trade, 
                      listw = W_list,
                      zero.policy = TRUE, 
                      tol.solve = 1e-60)

summary(sar_trade)


## Usual visual inspections
# Plot residuals to check normality
ggplot(data.frame(Residuals = residuals(sar_trade)), aes(x = Residuals)) +
  geom_density() +
  labs(title = "Residual Distribution of SAR model")

# And fitted vs residuals
ggplot(data.frame(Fitted = fitted(sar_trade), Residuals = residuals(sar_trade)), 
       aes(x = Fitted, y = Residuals)) +
  geom_point(color = "blue", size = 3, alpha = 0.7) +  # Residuals as points
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Zero reference line
  labs(title = "Residuals vs. Fitted Values",
       subtitle = "Points should be randomly scattered for a well-fitting model",
       x = "Fitted Values",
       y = "Residuals")

# Plot Actual GDP vs Fitted GDP for SAR model
plot(inter_trade$GDP, fitted(sar_trade), 
     main = "SAR: Actual vs Fitted GDP", 
     xlab = "Actual GDP", 
     ylab = "Fitted GDP", 
     col = "darkgreen", pch = 16)
abline(0, 1, col = "darkred")

### What happened to the sign of the number of trading partners term?
### Why? 
### Could we identify the spatial dependence? Value for Rho
### Does this model offer a better fit? Compare the AIC values (lower is better)





