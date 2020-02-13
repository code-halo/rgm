library(data.table)
library(igraph)
library(ggplot2)
library(Matrix)

setwd('.')

colnames <- read.table('../Datasets/Real/TCGA/ColNames.csv',header=FALSE);
colnames_tf <- read.table('../Datasets/Real/TCGA/ColNames_TF.csv',header=FALSE)

mutant_graph <- read.table('../Datasets/Real/TCGA/Gene_Mutation_A.csv',header=FALSE)
mutant_matrix <- matrix(0,nrow=length(colnames_tf$V1),ncol=length(colnames$V1))
for (i in 1:nrow(mutant_graph))
{
  src <- mutant_graph[i,1]
  tgt <- mutant_graph[i,2]
  mutant_matrix[src,tgt] <- mutant_graph[i,3]
}

wildtype_graph <- read.table('../Datasets/Real/TCGA/Gene_NonMutation_B.csv',header=FALSE)
wildtype_matrix <- matrix(0,nrow=length(colnames_tf$V1),ncol=length(colnames$V1))
for (i in 1:nrow(wildtype_graph))
{
  src <- wildtype_graph[i,1]
  tgt <- wildtype_graph[i,2]
  wildtype_matrix[src,tgt] <- wildtype_graph[i,3]
}

rownames(mutant_matrix) <- as.character(colnames_tf$V1)
colnames(mutant_matrix) <- as.character(colnames$x);
rownames(wildtype_matrix) <- as.character(colnames_tf$V1)
colnames(wildtype_matrix) <- as.character(colnames$x);

nodes_connected_components <- read.table('results/TCGA/output_node_cluster.csv',sep=",",header=FALSE);
nodes_connected_components$labels <- as.character(colnames_tf$V1[nodes_connected_components$V1]);

list_of_clusters <- unique(nodes_connected_components$V2);
for (i in 1:length(list_of_clusters))
{
  nodes_in_cluster_i <- nodes_connected_components[nodes_connected_components$V2==list_of_clusters[i],]$labels;
  list_mutant <- NULL;
  list_wildtype <- NULL;
  edgelist_mutant <- NULL;
  edgelist_wildtype <- NULL;
  mutant_subgraph_matrix <- mutant_matrix[nodes_in_cluster_i,];
  for(j in 1:nrow(mutant_subgraph_matrix))
  {
    connections <- which(mutant_subgraph_matrix[j,]>0);
    names(connections) <- rownames(mutant_subgraph_matrix)[which(mutant_subgraph_matrix[j,]>0)]
    if (length(connections)>0)
    {
      source_node <- rownames(mutant_subgraph_matrix)[j];
      list_mutant <- c(list_mutant,c(source_node,names(connections)));
      for(k in 1:length(connections))
      {
        edgelist_mutant <- rbind(edgelist_mutant,cbind(source_node,names(connections)[k],mutant_subgraph_matrix[source_node,connections[k]]));
      }
    }
  }
  edgelist_mutant <- as.data.frame(edgelist_mutant)
  colnames(edgelist_mutant) <- c("Source","Target","Weight");
  edgelist_mutant$Source <- as.character(as.vector(edgelist_mutant$Source))
  edgelist_mutant$Target <- as.character(as.vector(edgelist_mutant$Target))
  edgelist_mutant$Weight <- as.numeric(as.vector(edgelist_mutant$Weight))
  list_mutant <- unique(list_mutant);
  
  
  wildtype_subgraph_matrix <- wildtype_matrix[nodes_in_cluster_i,];
  for(j in 1:nrow(wildtype_subgraph_matrix))
  {
    connections <- which(wildtype_subgraph_matrix[j,]>0);
    names(connections) <- rownames(wildtype_subgraph_matrix)[which(wildtype_subgraph_matrix[j,]>0)];
    if (length(connections)>0)
    {
      source_node <- rownames(wildtype_subgraph_matrix)[j];
      list_wildtype <- c(list_wildtype,c(source_node,names(connections)));
      for(k in 1:length(connections))
      {
        edgelist_wildtype <- rbind(edgelist_wildtype,cbind(source_node,names(connections)[k],wildtype_subgraph_matrix[source_node,connections[k]]));
      }
    }
  }
  edgelist_wildtype <- as.data.frame(edgelist_wildtype)
  colnames(edgelist_wildtype) <- c("Source","Target","Weight");
  edgelist_wildtype$Source <- as.character(as.vector(edgelist_wildtype$Source))
  edgelist_wildtype$Target <- as.character(as.vector(edgelist_wildtype$Target))
  edgelist_wildtype$Weight <- as.numeric(as.vector(edgelist_wildtype$Weight))
  list_wildtype <- unique(list_wildtype);
  
  write.csv(list_mutant,file=paste0("results/TCGA/nodelist_mutant_cluster_",list_of_clusters[i],".csv"),row.names=F);
  write.csv(list_wildtype,file=paste0("results/TCGA/nodelist_wildtype_cluster_",list_of_clusters[i],".csv"),row.names=F);
  write.table(edgelist_mutant,file=paste0("results/TCGA/edgelist_mutant_cluster_",list_of_clusters[i],".csv"),row.names=F,col.names=T);
  write.table(edgelist_wildtype,file=paste0("results/TCGA/edgelist_wildtype_cluster_",list_of_clusters[i],".csv"),row.names=F,col.names=T);
}



