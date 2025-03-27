######## LISS2108 Statistical Inference for Social Networks Analysis ##########
# Author: Santiago Quintero (KCL) - 2025

# In this workshop, we'll work with the "General Relativity and Quantum Cosmology 
# collaboration network" collected by Leskovec et al (2007) (http://doi.acm.org/10.1145/1217299.1217301).
# This network contains all academic collaborations in the domain of GR&QC between
# 1993-2003. It is an undirected network that connects nodes (authors) if they 
# wrote a paper together. We'll practice fitting ERGMs with this network.

# Clear workspace and set seed for reproducibility
rm(list = ls())
set.seed(123)
options(scipen = 100) # turn off scientific notation


## Libraries we'll work with
packages <- c("igraph", "network", "dplyr", "ergm", "intergraph")

# If the packages are not installed, run:
# install.packages(packages)

# Call the libraries
lapply(packages, require, character.only = TRUE)


# Download the file and save it
download.file(url = "https://snap.stanford.edu/data/ca-GrQc.txt.gz",
              destfile = paste0(getwd(), "arxiv-net.txt.gz"))

# Read dataset
collab <- read.table(gzfile(paste0(getwd(), "arxiv-net.txt.gz")))
names(collab) <- c("FromNodeId", "ToNodeId")


# Delete self-loops
collab <- collab[collab$FromNodeId != collab$ToNodeId, ]

## This is a pretty large network, so we'll work with a subset 
## for the sake of simplicity and to facilitate ERGM running times
collab_small <- sample_n(collab, 1000)




# Create igraph object
net.igraph <- graph_from_edgelist(as.matrix(collab_small),
                           directed = FALSE)

## We'll only consider the nodes that have at least 3 collaborations (to simplify and reduce running times)
net.igraph <- delete_vertices(net.igraph, V(net.igraph)[degree(net.igraph) < 3]) 

## Delete remaining isolates
net.igraph <- delete_vertices(net.igraph, V(net.igraph)[degree(net.igraph) == 0])

## The dataset does not have node-level attributes, so let's create one
## randomly to use as a variable later. Let's assume that we add values for
## the prestige of the university each author is in. And let's assume that that
## value is correlated with the degree centrality of the author

degree_centrality <- degree(net.igraph) # centrality
random_variable <- rnorm(vcount(net.igraph), mean = 5, sd = 2) # add some noice

prestige <- degree_centrality * random_variable # create the prestige variable

V(net.igraph)$prestige <- prestige ## Add the variable to the network


## Let's first have a look at the network
summary(net.igraph)

## And graphically
plot(net.igraph,
     edge.width = 0.5,
     edge.col = "gray10",
     vertex.label = NA,
     vertex.size = V(net.igraph)$prestige/5,
     vertex.col = "gold",
     main = "GR & QC Co-Authorship Network")


## Transform the Igraph object into a Network object to run the ERGMs
net <- asNetwork(net.igraph)

detach(package:igraph) # sometimes the Network and Igraph packages clash, so better to keep them separate



##### Fitiing ERGMs ######

## First, we fit the simplest model possible: only considering the number of ties
ergm.base <- ergm(net ~ edges)  # Number of edges in the network (equal to
                                # the intercept term in a regression)
summary(ergm.base)  


## Now let's try to model the structural properties of the network directly.
## From the literature, we suspect that some important social forces
## behind academic collaboration have to do with popularity (more prolific)
## academics might collaborate more and with triadic closure (it is more
## likely to collaborate with collaborators of my collaborators). We first
## model these network dependencies using an "old" approach.

ergm.simple <- ergm(net ~ 
                      edges +
                      triangle +       # Number of triangles
                      kstar(2))        # Number of 2-stars 

summary(ergm.simple)

## How do we interpret these parameters? What is the unit of the coefficients


## Now, let's model these dependencies with a more "contemporary" approach
## to calculate these dependencies.
ergm.structure <- ergm(net ~ 
                    edges +                         
                    gwesp(0.5, fixed = TRUE) +      # Geometrically weighted edgewise shared partner distribution
                    gwdegree(0.5, fixed = TRUE))     # Geometrically weighted degree distribution

summary(ergm.structure)

## Are these parameters different? How do we interpret these?


## Now, let's add the node-level attribute we created before to the model
ergm.node.att <- ergm(net ~ 
                         edges +                         
                         gwesp(0.5, fixed = TRUE) +     
                         gwdegree(0.5, fixed = TRUE)+
                         nodecov("prestige"))    ## Add node attribute

summary(ergm.node.att)

# What if we exponentiate the coefficients
OR.ergm.node.att <- exp(coef(ergm.node.att))
print(OR.ergm.node.att)


## We can go even more sophisticated. What if authors prefer to collaborate with
## people in equally good institutions (i.e., prestige is similar)
ergm.node.att.diff <- ergm(net ~ 
                        edges +                         
                        gwesp(0.5, fixed = TRUE) +     
                        gwdegree(0.5, fixed = TRUE)+
                        nodecov("prestige")+
                        absdiff("prestige")) 

summary(ergm.node.att.diff)

# Exponentiate the coefficients
OR.ergm.node.att2 <- exp(coef(ergm.node.att.diff))
print(OR.ergm.node.att2)





####### Goodness of Fit checks #########

## Now, how do we know if our model is a good representation of the 
## observed network?
par(mfrow = c(2, 2))     # To setup the GOF display 

## GOF for the base model

gof_base <- gof(ergm.base)
plot(gof_base)


## MCMC diagnostics and GOF for the simple model
mcmc.diagnostics(ergm.simple)
gof_simple <- gof(ergm.simple)
plot(gof_simple)


## MCMC diagnostics and GOF for the more sophisticated model
mcmc.diagnostics(ergm.structure)
gof_structure <- gof(ergm.structure)
plot(gof_structure)


## MCMC diagnostics and GOF for the model with attributes
mcmc.diagnostics(ergm.node.att)
gof_node.att <- gof(ergm.node.att)
plot(gof_node.att)


## MCMC diagnostics and GOF for the model with attributes
mcmc.diagnostics(ergm.node.att.diff)
gof_node.att.diff <- gof(ergm.node.att.diff)
plot(gof_node.att.diff)


## It seems we are not capturing degree=2 very well. What
## if we simply add a term doing that?
ergm.node.att.d2 <- ergm(net ~ 
                         edges +                         
                         gwesp(0.5, fixed = TRUE) +     
                         gwdegree(0.5, fixed = TRUE)+
                         nodecov("prestige")+
                         absdiff("prestige")+
                         degree(2)) 

summary(ergm.node.att.d2)
## MCMC diagnostics and GOF for the model with attributes
mcmc.diagnostics(ergm.node.att.d2)
gof_node.att.d2 <- gof(ergm.node.att.d2)
plot(gof_node.att.d2)


## How do we select the best model? AIC, BIC?

###### TASK ######
## Now, you will run your own ERGMs using the Florentine network we worked
## with last class.
data(florentine) 

## Remeber there are two networks in this dataset 
print(flomarriage) ## The marriage network and 
plot(flomarriage,
     displaylabels = TRUE,    
     label.pos = 5,           
     label.cex = 0.8,         
     vertex.cex = flomarriage %v% "wealth"/15,  # Node size by wealth
     vertex.col = "gold",     
     edge.col = "gray20",     
     main = "Florentine Marriages Network - Node size by wealth")


print(flobusiness) ## The business network
plot(flobusiness,
     displaylabels = TRUE,    
     label.pos = 5,           
     label.cex = 0.8,         
     vertex.cex = flomarriage %v% "wealth"/15,  
     vertex.col = "pink",     
     edge.col = "gray20",     
     main = "Florentine Business Network - Node size by wealth")

## The networks are already set to run ERGMs. How would you model these networks?
## What parameters would you include? 
## Spend some time looking at the possibilities here: https://zalmquist.github.io/ERGM_Lab/ergm-terms.html



