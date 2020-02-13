library(data.table)
library(igraph)
library(dplyr)

setwd('.')

mutant_cluster_1 <- read.table('edgelist_mutant_cluster_1.csv',header=TRUE);
wildtype_cluster_1 <- read.table('edgelist_wildtype_cluster_1.csv',header=TRUE);
mutant_cluster_1$Source <- as.character(as.vector(mutant_cluster_1$Source))
mutant_cluster_1$Target <- as.character(as.vector(mutant_cluster_1$Target))
mutant_cluster_1$Weight <- as.numeric(as.vector(mutant_cluster_1$Weight))
wildtype_cluster_1$Source <- as.character(as.vector(wildtype_cluster_1$Source))
wildtype_cluster_1$Target <- as.character(as.vector(wildtype_cluster_1$Target))
wildtype_cluster_1$Weight <- as.numeric(as.vector(wildtype_cluster_1$Weight))



new_edgelist <- NULL;
for (i in 1:nrow(wildtype_cluster_1))
{
  source <- as.character(wildtype_cluster_1$Source[i]);
  target <- as.character(wildtype_cluster_1$Target[i]);
  if (source!=target)
  {
    if (length(which(mutant_cluster_1$Source==source & mutant_cluster_1$Target==target))>0)
    {
      weight <- 0
      
      new_edgelist <- rbind(new_edgelist,cbind(source,target,weight));
    }
    else{
      weight <- 1;
      new_edgelist <- rbind(new_edgelist,cbind(source,target,weight));
    }
  }
}
for (i in 1:nrow(mutant_cluster_1))
{
  source <- mutant_cluster_1$Source[i]
  target <- mutant_cluster_1$Target[i]
  if (source!=target)
  {
    if (length(which(wildtype_cluster_1$Source==source & wildtype_cluster_1$Target==target))>0)
    {}
    else{
      weight <- -1
      new_edgelist <- rbind(new_edgelist,cbind(source,target,weight))
    }
  }
}
new_edgelist <- as.data.frame(new_edgelist);
colnames(new_edgelist) <- c("Source","Target","Type")
new_edgelist$Source <- as.character(as.vector(new_edgelist$Source));
new_edgelist$Target <- as.character(as.vector(new_edgelist$Target));
new_edgelist$Type <- as.numeric(as.vector(new_edgelist$Type));

write.table(new_edgelist,"edgelist_info_1.csv",col.names=T,row.names=F)

#Get node list
node_list_mutant_cluster_2 <- read.table("nodelist_mutant_cluster_2.csv",header=TRUE);
node_list_wildtype_cluster_2 <- read.table("nodelist_wildtype_cluster_2.csv",header=TRUE);


#GO TERM Analysis
go_terms_mutant_cluster_2 <- read.table("enriched_goterms_mutant_cluster_2.csv",header=TRUE,sep="\t");
go_terms_wildtype_cluster_2 <- read.table("enriched_goterms_wildtype_cluster_2.csv",header=TRUE,sep="\t")

common_goterms <- intersect(go_terms_mutant_cluster_2$term_goid,go_terms_wildtype_cluster_2$term_goid);
go_terms_common_cluster_2 <- go_terms_mutant_cluster_2[go_terms_mutant_cluster_2$term_goid %in% common_goterms,]

go_terms_mutant_cluster_2 <- go_terms_mutant_cluster_2[!go_terms_mutant_cluster_2$term_goid %in% common_goterms,]
go_terms_wildtype_cluster_2 <- go_terms_wildtype_cluster_2[!go_terms_wildtype_cluster_2$term_goid %in% common_goterms,]


table(go_terms_common_cluster_2$term_category);
table(go_terms_mutant_cluster_2$term_category)
table(go_terms_wildtype_cluster_2$term_category)

#Pathway Analysis
go_pathways_mutant_cluster_2 <- read.table("enriched_pathways_mutant_cluster_2.csv",header=TRUE,sep="\t");
go_pathways_wildtype_cluster_2 <- read.table("enriched_pathways_wildtype_cluster_2.csv",header=TRUE,sep="\t");

common_pathways <- intersect(go_pathways_mutant_cluster_2$external_id,go_pathways_wildtype_cluster_2$external_id)
go_pathways_common_cluster_2 <- go_pathways_mutant_cluster_2[go_pathways_mutant_cluster_2$external_id%in% common_pathways,]

go_pathways_mutant_cluster_2 <- go_pathways_mutant_cluster_2[!go_pathways_mutant_cluster_2$external_id%in%common_pathways,]
go_pathways_wildtype_cluster_2 <- go_pathways_wildtype_cluster_2[!go_pathways_wildtype_cluster_2$external_id%in%common_pathways,]