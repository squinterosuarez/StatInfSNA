######## LISS2108 Statistical Inference for Social Networks Analysis ##########
# Author: Santiago Quintero (KCL) - 2025

# In this workshop, we'll use a canonical network dataset, the Florentine marriages 
# from Padgett and Ansell (1993), to practice the application of Quadratic Assignment
# Procedure to predict network ties (dyad-level analysis)

# Clear workspace and set seed for reproducibility
rm(list = ls())
set.seed(123)
options(scipen = 100) # turn off scientific notation

## Libraries we'll work with
packages <- c("sna", "network", "ergm")

# sna      - For the QAP tests
# network  - For network manipulations
# ergm     - For the Florentine network

# If the packages are not installed, run:
# install.packages(packages)

# Call the libraries
lapply(packages, require, character.only = TRUE)



######## Florentine families networks ######## 
## First, we call the main datasets we'll work with, the marriage and business ties
## network from the most important families in Florence in the 15th century
# ?florentine  - Run this to learn more about the datase
data(florentine) 


## These data contain two network objects we'll use, the marriages network (flomarriage),
## and the business network (flobusiness)

## Let's first inspect both networks

### Marriages network ###
print(flomarriage)

plot(flomarriage,
     displaylabels = TRUE,    # Show node labels
     label.pos = 5,           # Center labels on nodes
     label.cex = 0.8,         # Label size
     vertex.cex = 3,          # Node size
     vertex.col = "gold",     # Node color
     edge.col = "gray20",     # Edge color
     main = "Florentine Marriages Network")

## Resize the nodes by wealth
plot(flomarriage,
     displaylabels = TRUE,    
     label.pos = 5,           
     label.cex = 0.8,         
     vertex.cex = flomarriage %v% "wealth"/15,  # Node size by wealth
     vertex.col = "gold",     
     edge.col = "gray20",     
     main = "Florentine Marriages Network - Node size by wealth")


### Business network ###
print(flobusiness)

# Plot the business network
plot(flobusiness,
     displaylabels = TRUE,    
     label.pos = 5,           
     label.cex = 0.8,         
     vertex.cex = 3,          
     vertex.col = "pink",     
     edge.col = "gray20",     
     main = "Florentine Business Network")


## Resize this one by wealth too
plot(flobusiness,
     displaylabels = TRUE,    
     label.pos = 5,           
     label.cex = 0.8,         
     vertex.cex = flomarriage %v% "wealth"/15,  
     vertex.col = "pink",     
     edge.col = "gray20",     
     main = "Florentine Business Network - Node size by wealth")


## How are the two networks different?


########  Bivariate correlation using QAP ######## 

## First, we'll check whether there is a correlation between marriage and business
## connections: Are families with marriage ties more likely to do business?

# Let's run it as a simple correlation
marriage_mat <- as.matrix(flomarriage)
business_mat <- as.matrix(flobusiness)

# We need to vectorise the matrices to be able to compute a simple correlation
cor_flo <- cor.test(as.vector(marriage_mat), as.vector(business_mat))

print(cor_flo)     # How correlated are the two networks? 
                   # Is this correlation statistically significant?


### Now let's compute the correlation using QAP. We use the 'qaptest' function
### from the sna package

qap_flo <- qaptest(marriage_mat, gcor,  business_mat, reps = 1000)  # Make sure to use the matrices and not the network objects
summary(qap_flo)  # Is the correlation different? Why? 
plot(qap_flo)     # Is it statistically significant?


########  Multivariate Regression QAP (MRQAP) ######## 

## First, let's say we want to predict "strength" of the relationship between
## the families. For that we'll create a new network that adds the business
## and marriage ties.

# Create a combined network by summing the two networks
strength_network <- marriage_mat + business_mat

### The dataset contains some node attributes that we can use to further
### control for  other factors that can be affecting the strength of the relationship
### between families. We have access to data about families' wealth and number of 
### priorates (seats on the civic council - a proxy for political power). 
### We'll calculate a dyadic variable that reflect the difference in these
### attributes between families.

# Extract the 'wealth' and priorates attributes from one of the networks
wealth <- flomarriage %v% "wealth"
priorates <- flomarriage %v% "priorates"

# Create a dyadic matrix of absolute differences in wealth and priorates
wealth_diff <- abs(outer(wealth, wealth, "-"))
priorates_diff <- abs(outer(priorates, priorates, "-"))


## Now we have a valued network, with values of 2 for families that have both
## marriage and business connections. We'll use the wealth and prioirates 
## differentials networks we created  to predict the strength of the
## family ties. As the ties are not binary, we use an MRQAP.

mrqap_result <- netlm(strength_network, # matrix to be predicted - make sure to use the matrices and not the network objects
                      list(wealth_diff, priorates_diff),  ## predictors
                      nullhyp = "qap",                                      
                      mode = "graph", # specify that it is an undirected network
                      reps = 1000)                        

summary(mrqap_result)
## What is the effect of wealth and political power differentials?
## Are these effects significant?



########  Logistic Regression QAP (LRQAP) ######## 

### Now, we want to further understand what factors might lead to families
### engaging in business ties. For this, we'll predict the business network
### using the marriage network and the attributes differentials. Because 
### the business network is binary, we are better off by directly estimating
### a probability model, so we use a LRQAP.


## Run the logit QAP 
lrqap_flo <- netlogit(business_mat, 
                      list(marriage_mat, wealth_diff, priorates_diff),
                      nullhyp = "qap",
                      mode = "graph",
                      reps = 1000)

summary(lrqap_flo) 

## What do the results tell us about the factors driving business relationships?


######## TASK ########  

### Now, let's say our measure of political power relationship is not very interesting.
### What if we want to understand the political power from a homophily perspective,
### i.e., if two families are very powerful, they might be more likely to form ties.
### Usually, for this approach, we use the *product* of the attributes we want to analyse.
### Calculate a dyadic matrix multiplying the number of priorates (instead of the difference).
### Additionally, there is another nodal attribute in the datasets, 'totalties ', which
### refers to the total number of marriages and business connections each family 
### has with a wider number of families in Florence (including many more than the
### ones in the main dataset).With these new variables --the product of priorates and
### the dyadic transformation of your choice for the total ties variable 
### (justify why using the absolute difference or the product)-- run the 
### adequate analysis to predict the emergence of *marriage ties*. 
### Also, use the business network as a predictor.
















