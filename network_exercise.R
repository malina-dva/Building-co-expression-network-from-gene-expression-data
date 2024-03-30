# Chapter 1: Install and load the required packages 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install('mixOmics')
# BiocManager::install('igraph')
# BiocManager::install('RColorBrewer')

library(mixOmics)
library(igraph)
library(RColorBrewer)

# Chapter 2: Set working directory  
# Directory for all output files
setwd("C://Users//Malina//Desktop//igraphs//")
""

# Chapter 3: Load the data 
data(breast.TCGA)
mRNA <- breast.TCGA$data.train$mrna
dim(mRNA)  # Display dimensions of the data
head(mRNA[, 1:20])  # Display first 20 rows of mRNA data

# Chapter 4: Calculate pairwise correlation between genes and create an adjacency matrix
# Setting the right correlation threshold is very important
correlation_matrix <- cor(mRNA, method = "spearman")
head(correlation_matrix[,1:30])
threshold <- 0.5

# Create an adjacency matrix based on the correlation threshold
adjacency_matrix <- ifelse(abs(correlation_matrix) > threshold, 1, 0)
head(adjacency_matrix[,1:20])

# Remove the diagonal values 
diag(adjacency_matrix) <- 0
head(adjacency_matrix[,1:20])

# Chapter 5: Build undirected, unweighted graph (network)
graph <- graph_from_adjacency_matrix(adjacency_matrix , 
                                     mode = "undirected")
graph

# Chapter 6: Get the network topology 
# Cluster based on edge betweenness
ceb <- cluster_edge_betweenness(graph)

# Extract community memberships from the clustering
community_membership <- membership(ceb)
head(community_membership)
length(ceb)

# Count the number of nodes in each cluster
cluster_sizes <- table(community_membership)
cluster_sizes

# Chapter 7: Add attributes to the vertices
# Add cluster memberships as a vertex attribute
igraph::V(graph)$cluster_membership <- community_membership

# Map cluster memberships to colors using a color palette function
cluster_colors <- rainbow(max(community_membership))
head(cluster_colors)

# Create a color variable based on cluster memberships
V(graph)$color_variable <- cluster_colors[community_membership]

# vetrex degree
vertex_degrees <- igraph::degree(graph)
head(vertex_degrees)

# Add vertex degree as a vertex attribute
igraph::V(graph)$degree <- vertex_degrees

# Chapter 8: Combine the attributes of the vertices in a list
vertex.attr <- list(
  cluster_membership = V(graph)$cluster_membership,
  name = V(graph)$name,
  degrees = V(graph)$degree,
  color = V(graph)$color_variable)

# Chapter 9: Save the GraphML file
# Specify the file name for saving the GraphML file
graphml_file <- "graph_mRNA_with_metadata_unweighted_undirected.graphml"

# Write the igraph object to GraphML format with metadata
write_graph(
  graph,
  graphml_file,
  format  = "graphml")

# Chapter 10: Save the vertices attributes as a dataframe
# Combine vertex attributes into a data frame
vertex_data <- data.frame(vertex.attr)

head(vertex_data)
write.csv(vertex_data, "vertex_data_mRNA_with_metadata_unweighted_undirected.csv", row.names = FALSE)

# Chapter 11: Bonus Analyses 
# Section 1: Additional measures of network topology 

# Calculate degree centrality
degree_centrality <- degree(graph)
mean(degree_centrality)

# Average path length
average.path.length(graph)

# Calculate betweenness centrality
betweenness_centrality <- betweenness(graph)
mean(betweenness_centrality)

# Add the betweenness centrality as an attribute of the vertices
V(graph)$betweenness = betweenness_centrality

# Calculate closeness centrality and add as an attribute of the vertices
closeness_centrality <- igraph::closeness(graph)
V(graph)$closeness = closeness_centrality

# Calculate eigenvector and add as an attribute of the vertices
Eigenvectors  = eigen_centrality(graph)$vector
V(graph)$eigen = Eigenvectors

# Add the additional node properties to an extended dataframe
vertex.attr.ext <- list(
  cluster_membership = V(graph)$cluster_membership,
  name = V(graph)$name,
  degrees = V(graph)$degree,
  color = V(graph)$color_variable, 
  betweeness = V(graph)$betweenness,
  closeness = V(graph)$closeness,
  eigenvectors = V(graph)$eigen)

vertex_data_ext <- data.frame(vertex.attr.ext)
head(vertex_data_ext)

# In case you'd like to save the extended vertex attr list you can do 
write.csv(vertex_data_ext, "vertex_data_mRNA_with_metadata_EXT_unweighted_undirected.csv", row.names = FALSE)

# Calculate clustering coefficient
clustering_coefficient <- transitivity(graph, type = "global")
clustering_coefficient

# Compute and plot degree distribution
ddist <- igraph::degree.distribution(graph)

# Data frame needed for ggplot2
df <- data.frame(Degree = as.factor((seq_along(ddist)) - 1),
                 Fraction = ddist)

ggplot(data = df, aes(x = Degree, y = Fraction, group = 1)) +
  geom_line() +
  geom_point() +
  theme_bw()

# Section 2: Additional method to obtain community membership 

# Alternative 2: Using cluster_walktrap
cwt <- cluster_walktrap(graph)
membership2 <- membership(cwt)
table(membership2)

# Section 3: Creating network with edge weights

data(breast.TCGA)
mRNA <- breast.TCGA$data.train$mrna

# Calculate pairwise correlation between genes
correlation_matrix <- cor(mRNA,method = "spearman")

# Setting the right correlation threshold is very important
threshold <- 0.5

# Create an adjacency matrix based on the correlation threshold
# This matrix and graph will only serve to get directional correlation scores
# Weights have to be positive values for the clustering to work 
adjacency_matrix <- ifelse(abs(correlation_matrix) > threshold, as.numeric(correlation_matrix), 0)
head(adjacency_matrix[,1:40])
diag(adjacency_matrix) <- 0

graph1 <- graph_from_adjacency_matrix(adjacency_matrix , 
                                      mode = "undirected", weighted = TRUE)
graph1

# Create the graph with weights (absolute (positive) correlation values)

adjacency_matrix <- ifelse(abs(correlation_matrix) > threshold, as.numeric(abs(correlation_matrix)), 0)
head(adjacency_matrix[,1:40])
diag(adjacency_matrix) <- 0

# Now the weighted option is going to be set to TRUE
graph2 <- graph_from_adjacency_matrix(adjacency_matrix , 
                                      mode = "undirected", weighted = TRUE)
graph2

# We can have a quick sneak peak at the values of the weights
# They are assigned automatically to the edge graph property (E) as weight
# e.g E(graph2)$weight
head(E(graph2)$weight)

# Now we will proceed to get the vertices attributes (short version)

# Calculate the degree of each vertex
vertex_degrees <- igraph::degree(graph2)

# Cluster based on cluster_fast_greedy, which is a function that considers the effect of the weights when defining clusters

cfg <- cluster_fast_greedy(graph2)
membership2 <- membership(cfg)
table(membership2)

# Cluster 1,2,3 have the highest number of members

# Extract community memberships from the clustering
community_membership <- membership(cfg)

# Add cluster memberships as a vertex attribute
igraph::V(graph2)$cluster_membership <- community_membership

# Calculate vetrex degree
vertex_degrees <- igraph::degree(graph2)

# Add degree as a vertex attribute
igraph::V(graph2)$degree <- vertex_degrees


# Map cluster memberships to colors using a color palette function
# Using ColorBrewer palettes
# The number of colors in Set1 is limited, so not all clusters will have color, but in reality we only need colors for the first several clusters that had most of members

cluster_colors <- brewer.pal(max(community_membership), "Set1")

# Create a color variable based on cluster memberships

V(graph2)$color_variable <- cluster_colors[community_membership]

# Weight is already added as an attr of the edges
# Now we are going to add the corr scores from graph1 as corr_scores attr to graph2 
# Corr scores from graph1 are actually exactly the same weights as the weights of graph2 but with sign (pos or neg) respectively

E(graph2)$corr_scores = E(graph1)$weight

# Add the node properties to the vertex attr list

vertex.attr = list(
  cluster_membership = V(graph2)$cluster_membership,
  name = V(graph2)$name,
  degrees = V(graph2)$degree,
  color = V(graph2)$color_variable,
  size  = V(graph2)$degree)

# Add the edge property to the edge attr list
edge.attr = list(
  corr_scores = round(E(graph2)$corr_scores, 3 ))

# Specify the file path for saving the GraphML file
graphml_file <- "graph_mRNA_with_metadata_undirected_WEIGHTED.graphml"

# Write the igraph object to GraphML format with metadata
write_graph(
  graph2,
  graphml_file,
  format  = "graphml")
