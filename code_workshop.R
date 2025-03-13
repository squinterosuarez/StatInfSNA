######## LISS2108 Statistical Inference for Social Networks Analysis ##########
# Workshop: Inference in Social Network Analysis
# Level 1: Group-Level Analysis using Permutation Tests
# Level 2: Node-Level Analysis using Autoregressive (Network Lag) Models
#
# In the group-level example, we ask:
#   "Does a more centralized friendship network improve the effects of health interventions in a classroom?"
#
# In the node-level example, we ask:
#   "Do countries with more trading partners have higher GDP?"
#
# In this workshop we:
#   1. Create simulated datasets for each level.
#   2. Run a standard linear regression and then a specialized method:
#        - Permutation test for the group-level (small, nonrandom sample)
#        - A network autoregressive model for the node-level, which uses an autoregressive term
#          defined as the weighted sum of the GDPs of the *other* countries.
#   3. Explicitly set the number of permutations in the permutation test.
#   4. For the node-level, include the degree (i.e. the number of trading partners) as a regressor
#      in addition to the network lag variable.
#   5. Print a side-by-side comparison of p-values from the standard regression vs. the network AR model.
#   6. Provide two exercises with incomplete code for further practice.
#
# ---------------------------
# Preliminary Setup

# Clear workspace and set seed for reproducibility
rm(list = ls())
set.seed(123)

# Load required libraries
library(ggplot2)      # For plotting
library(dplyr)        # For data manipulation
library(lmPerm)       # For permutation tests in linear models

# ---------------------------
# Group-Level Analysis: Permutation Test
# Example: Does a more centralized friendship network improve the effects of health interventions?
#
# Simulate data for 20 classrooms:
#  - 'Centralization': a network centralization score between 0.2 and 0.8.
#  - 'SmokingReduction': outcome, modeled as:
#       SmokingReduction = 5 + 10 * Centralization + error

n_classrooms <- 20
centralization <- runif(n_classrooms, min = 0.2, max = 0.8)
smoking_reduction <- 5 + 10 * centralization + rnorm(n_classrooms, mean = 0, sd = 2)
classroom_data <- data.frame(
  Classroom = 1:n_classrooms,
  Centralization = centralization,
  SmokingReduction = smoking_reduction
)

# Standard linear regression (assuming independent classrooms):
lm_group <- lm(SmokingReduction ~ Centralization, data = classroom_data)
cat("Standard Linear Regression (Group-Level) Summary:\n")
print(summary(lm_group))
se_lm <- coef(summary(lm_group))[, "Std. Error"]

# Permutation test: for small samples with nonrandom sampling,
# we use lmPerm's lmp() function.
# Explicitly set the number of permutations:
nperm <- 10000  # (Adjust this value as needed)
perm_model <- lmp(SmokingReduction ~ Centralization, data = classroom_data, nperm = nperm)
cat("\nPermutation Test (lmPerm) Summary (with", nperm, "permutations):\n")
print(summary(perm_model))
se_perm <- coef(summary(perm_model))[, "Std. Error"]

cat("\nComparison of Standard Errors (lm() minus lmp()):\n")
print(se_lm - se_perm)
# (Differences may occur because the permutation approach uses the randomization distribution.)

# Plot the relationship:
p1 <- ggplot(classroom_data, aes(x = Centralization, y = SmokingReduction)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Group-Level Analysis: Centralization vs Smoking Reduction",
       subtitle = paste("n =", n_classrooms, "classrooms"),
       x = "Network Centralization",
       y = "Reduction in Smoking")
print(p1)

# ---------------------------
# Node-Level Analysis: Network Autoregressive Model
# Example: Do countries with more trading partners have higher GDP?
#
# We simulate data for 50 countries.
# We also simulate a random binary trading network (adjacency matrix) among countries.
#
# Here we define:
#   - TradingPartners: the degree (number of trading partners) of a country.
#   - NetworkGDP: for each country, the weighted sum of the GDPs of all other countries,
#       computed via a row-normalized trading matrix (with zeros on the diagonal).
#
# The standard regression will be: GDP ~ TradingPartners
# The network autoregressive model will be: GDP ~ NetworkGDP + TradingPartners
# (Thus, the degree is explicitly included as a regressor.)

n_countries <- 50
# Simulate a random binary trading network (sparse)
trade_matrix <- matrix(rbinom(n_countries * n_countries, 1, prob = 0.1), nrow = n_countries)
diag(trade_matrix) <- 0  # Remove self-ties

# Simulate baseline GDP:
# Let the number of trading partners positively influence GDP.
trading_partners <- rowSums(trade_matrix)
gdp <- 1000 + 500 * trading_partners + rnorm(n_countries, mean = 0, sd = 1000)
country_data <- data.frame(
  Country = 1:n_countries,
  GDP = gdp,
  TradingPartners = trading_partners
)

# Standard regression: GDP ~ TradingPartners
lm_country <- lm(GDP ~ TradingPartners, data = country_data)
cat("\nStandard Regression (Node-Level) Summary:\n")
print(summary(lm_country))
se_lm_country <- coef(summary(lm_country))[, "Std. Error"]

# Construct network weight matrix P by row-normalizing the trade matrix.
# (Rows with zero sum become all 0.)
W <- trade_matrix / rowSums(trade_matrix)
W[is.na(W)] <- 0

# Compute the network autoregressive term:
# For each country, NetworkGDP = weighted sum of the GDPs of the other countries.
NetworkGDP <- as.numeric(W %*% country_data$GDP)
country_data$NetworkGDP <- NetworkGDP

# The "degree" variable (number of trading partners) is already in TradingPartners.
# Now fit the network autoregressive model:
# GDP is modeled as a function of:
#    (i) NetworkGDP: the network lag term,
#   (ii) TradingPartners: the degree.
nar_model <- lm(GDP ~ NetworkGDP + TradingPartners, data = country_data)
cat("\nNetwork Autoregressive Model Summary:\n")
print(summary(nar_model))
se_nar <- coef(summary(nar_model))[, "Std. Error"]

# Compare p-values between the two models.
# Standard regression (only TradingPartners) vs. network AR model (NetworkGDP + TradingPartners):
p_lm_country <- coef(summary(lm_country))[, "Pr(>|t|)"]
p_nar <- coef(summary(nar_model))[, "Pr(>|t|)"]

# For clarity, create a comparison table.
# Note: lm_country has only two rows (Intercept and TradingPartners)
# whereas nar_model has three rows (Intercept, NetworkGDP, TradingPartners).
# We can compare the p-value for TradingPartners across models.
pvals_comparison <- data.frame(
  Predictor = c("Intercept", "TradingPartners"),
  StandardRegression = c(coef(summary(lm_country))["(Intercept)", "Pr(>|t|)"],
                         coef(summary(lm_country))["TradingPartners", "Pr(>|t|)"]),
  NetworkAR = c(coef(summary(nar_model))["(Intercept)", "Pr(>|t|)"],
                coef(summary(nar_model))["TradingPartners", "Pr(>|t|)"])
)
cat("\nComparison of p-values for TradingPartners and the Intercept:\n")
print(pvals_comparison)

# We may also report the p-value for the autoregressive term from the network AR model:
cat("\nP-value for the NetworkGDP coefficient (autoregressive term):\n")
print(coef(summary(nar_model))["NetworkGDP", "Pr(>|t|)"])

# Plot the standard regression relationship:
p2 <- ggplot(country_data, aes(x = TradingPartners, y = GDP)) +
  geom_point(size = 3, color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(title = "Standard Regression: Trading Partners vs GDP",
       x = "Number of Trading Partners",
       y = "GDP")
print(p2)

# ---------------------------
# Student Exercises: Incomplete Code for Further Practice
#
# Exercise 1: Group-Level Analysis
# Research Question: Does the density (rather than centralization) of a classroom's friendship network
# influence the effectiveness of a health intervention?
#
# Simulate data for 20 classrooms:
set.seed(456)
density <- runif(n_classrooms, min = 0.1, max = 0.9)  # Network density
intervention_effect <- 3 + 7 * density + rnorm(n_classrooms, mean = 0, sd = 1.5)
exercise1_data <- data.frame(
  Classroom = 1:n_classrooms,
  Density = density,
  InterventionEffect = intervention_effect
)

# Incomplete code for students to complete:
# Use lmPerm::lmp to perform a permutation test:
# (Fill in the number of permutations, e.g., nperm = 5000)
# perm_ex1 <- lmp(InterventionEffect ~ Density, data = exercise1_data, nperm = ________)
# summary(perm_ex1)

#
# Exercise 2: Node-Level Analysis
# Research Question: Do countries with higher centrality in the trade network have higher GDP growth rates?
#
# Simulate data for 50 countries:
set.seed(789)
gdp_growth <- runif(n_countries, min = -2, max = 10)  # GDP growth rates (in percent)
# For network centrality, we use the same trade_matrix as above.
degree_centrality <- rowSums(trade_matrix)
exercise2_data <- data.frame(
  Country = 1:n_countries,
  GDPGrowth = gdp_growth,
  DegreeCentrality = degree_centrality
)

# Incomplete code for students to complete:
# (1) Compute the network lag of GDPGrowth using the weight matrix W.
# network_lag_growth <- as.numeric(W %*% exercise2_data$GDPGrowth)
#
# (2) Fit a network autoregressive model:
# model_ex2 <- lm(GDPGrowth ~ network_lag_growth + DegreeCentrality, data = exercise2_data)
# summary(model_ex2)
# (Hint: Compare the results with a simple regression of GDPGrowth on DegreeCentrality.)

# ---------------------------
# End of Workshop Script