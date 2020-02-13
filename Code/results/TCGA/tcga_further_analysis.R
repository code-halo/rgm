library(data.table)
library(lattice)
library(Matrix)
library(ggplot2)
library(ggpubr)

setwd(".")

file1 <- read.table(file = "../../../Datasets/Real/TCGA/X.txt",header=FALSE)
file2 <- read.table(file = "../../../Datasets/Real/TCGA/Y.txt",header=FALSE)

tf_names <- read.csv("../../../Datasets/Real/TCGA/ColNames_TF.csv",header=FALSE)
tf_names <- as.character(as.vector(tf_names$V1))

pspi_result_P <- read.table(file="Pinf/p_nu_1.000000_mu_1.000000.txt",sep=",")
nips_result_P <- read.table(file="Pinf/P_NIPS.txt",sep=",")
pspi_bijection <- NULL
nips_bijection <- NULL
for (i in 1:nrow(pspi_result_P))
{
  pspi_bijection <- c(pspi_bijection,which(pspi_result_P[i,]>0))
  nips_bijection <- c(nips_bijection,which(nips_result_P[i,]>0))
}
df_nodes <- cbind(c(1:length(pspi_bijection)),pspi_bijection,nips_bijection)
df_nodes <- as.data.frame(df_nodes)
colnames(df_nodes) <- c("Original","New_PSPI_Ones","New_NIPS_Ones")
df_nodes$Mapping_PSPI <- df_nodes$Original-df_nodes$New_PSPI_Ones
df_nodes$Mapping_PSPI <- as.numeric(df_nodes$Mapping_PSPI==0)
df_nodes$Mapping_NIPS <- df_nodes$Original-df_nodes$New_NIPS_Ones
df_nodes$Mapping_NIPS <- as.numeric(df_nodes$Mapping_NIPS==0)


df_nodes$X_PSPI_nodes <- tf_names[df_nodes$New_PSPI_Ones]
df_nodes$X_NIPS_nodes <- tf_names[df_nodes$New_NIPS_Ones]
df_nodes$Y_nodes <- tf_names[df_nodes$Original]

X <- file1
Y <- file2
rownames(X) <- tf_names
colnames(X) <- tf_names
rownames(Y) <- tf_names
colnames(Y) <- tf_names

new_bijection_match <- NULL
for (i in 1:nrow(df_nodes))
{
  node_pspi_x <- df_nodes$X_PSPI_nodes[i]
  node_nips_x <- df_nodes$X_NIPS_nodes[i]
  node_name_y <- df_nodes$Y_nodes[i]
  x_pspi_genes <- names(X[node_pspi_x,which(X[node_pspi_x,]>0)])
  x_nips_genes <- names(X[node_nips_x,which(X[node_nips_x,]>0)])
  y_genes <- names(Y[node_name_y,which(Y[node_name_y,]>0)])
  pspi_y_genes <- names(Y[node_pspi_x,which(Y[node_pspi_x,]>0)])
  nips_y_genes <- names(Y[node_nips_x,which(Y[node_nips_x,]>0)])
  intersect_pspi_genes <- intersect(x_pspi_genes,y_genes)
  intersect_nips_genes <- intersect(x_nips_genes,y_genes)
  original_pspi_intersect <- intersect(x_pspi_genes,pspi_y_genes)
  original_nips_intersect <- intersect(x_nips_genes,nips_y_genes)
  union_pspi <- union(x_pspi_genes,y_genes)
  union_nips <- union(x_nips_genes,y_genes)
  temp <- cbind(node_pspi_x,node_nips_x,node_name_y,
                length(x_pspi_genes),length(x_nips_genes),
                length(y_genes),length(intersect_pspi_genes),
                length(intersect_nips_genes),length(original_pspi_intersect),
                length(original_nips_intersect),length(union_pspi),length(union_nips))
  new_bijection_match <- rbind(new_bijection_match,temp)
}
new_bijection_match <- as.data.frame(new_bijection_match)
colnames(new_bijection_match) <- c("Gene_PSPI_X","Gene_NIPS_X","Gene_Y",
                                   "No_Edges_PSPI_X","No_Edges_NIPS_X","No_Edges_Y",
                                   "No_Matched_Edges_PSPI","No_Matched_Edges_NIPS",
                                   "No_Orig_Matched_Edges_PSPI","No_Orig_Matched_Edges_NIPS",
                                   "Union_PSPI_Y","Union_NIPS_Y")
new_bijection_match$Gene_PSPI_X <- as.character(as.vector(new_bijection_match$Gene_PSPI_X))
new_bijection_match$Gene_NIPS_X <- as.character(as.vector(new_bijection_match$Gene_NIPS_X))
new_bijection_match$Gene_Y <- as.character(as.vector(new_bijection_match$Gene_Y))
new_bijection_match$No_Edges_PSPI_X <- as.numeric(as.vector(new_bijection_match$No_Edges_PSPI_X))
new_bijection_match$No_Edges_NIPS_X <- as.numeric(as.vector(new_bijection_match$No_Edges_NIPS_X))
new_bijection_match$No_Edges_Y <- as.numeric(as.vector(new_bijection_match$No_Edges_Y))
new_bijection_match$No_Matched_Edges_PSPI <- as.numeric(as.vector(new_bijection_match$No_Matched_Edges_PSPI))
new_bijection_match$No_Matched_Edges_NIPS <- as.numeric(as.vector(new_bijection_match$No_Matched_Edges_NIPS))
new_bijection_match$No_Orig_Matched_Edges_PSPI <- as.numeric(as.vector(new_bijection_match$No_Orig_Matched_Edges_PSPI))
new_bijection_match$No_Orig_Matched_Edges_NIPS <- as.numeric(as.vector(new_bijection_match$No_Orig_Matched_Edges_NIPS))
new_bijection_match$Union_PSPI_Y <- as.numeric(as.vector(new_bijection_match$Union_PSPI_Y))
new_bijection_match$Union_NIPS_Y <- as.numeric(as.vector(new_bijection_match$Union_NIPS_Y))
new_bijection_match <- new_bijection_match[order(new_bijection_match$No_Matched_Edges_PSPI,decreasing = T),]
write.table(new_bijection_match,"Matching_of_P_matrices.csv",row.names=F,col.names = T,quote=F)

indices_PSPI <- which(new_bijection_match$No_Edges_PSPI_X>0)
indices_NIPS <- which(new_bijection_match$No_Edges_NIPS_X>0)
#important_PSPI <- new_bijection_match[new_bijection_match$No_Edges_PSPI_X>0,]$No_Matched_Edges_PSPI
#important_NIPS <- new_bijection_match[new_bijection_match$No_Edges_NIPS_X>0,]$No_Matched_Edges_NIPS
#important_NIPS <- sort(important_NIPS,decreasing = T)
bijection_pspi <- new_bijection_match[indices_PSPI,]
bijection_nips <- new_bijection_match[indices_NIPS,]
jacc_PSPI <- bijection_pspi$No_Matched_Edges_PSPI/bijection_pspi$Union_PSPI_Y
jacc_NIPS <- bijection_nips$No_Matched_Edges_NIPS/bijection_nips$Union_NIPS_Y
N <- length(jacc_PSPI)
#important_df <- cbind(c(important_PSPI,important_NIPS),c(rep("PSPI",N),rep("MGM",N)))
important_df <- cbind(c(jacc_PSPI,jacc_NIPS),c(rep("STEPD",N),rep("MGM",N)))
important_df <- as.data.frame(important_df)
colnames(important_df) <- c("Density_Info","Method")
important_df$Density_Info <- as.numeric(as.vector(important_df$Density_Info))
important_df$Method <- as.character(as.vector(important_df$Method))

#Sample data
g <- ggplot(important_df, aes(x = Density_Info, fill = Method)) + 
     geom_density(aes(fill=Method), alpha = 0.35) + 
     geom_vline(data=important_df[important_df$Method=="STEPD",], aes(xintercept = mean(Density_Info)), color = "blue", 
                linetype = "dashed") + theme_bw() +
     geom_vline(data=important_df[important_df$Method=="MGM",], aes(xintercept = mean(Density_Info)), color = "red", 
             linetype = "dashed") + theme_bw() + ylab("Percentage of total TFs") +
     xlab("Distribution of Jaccard coefficient") + ggtitle("On Avg Differnce of 0.03 in Jaccard Coefficient by PSPI vs MGM") +
     theme(plot.title = element_text(size=20,hjust = 0.5)) + theme(plot.title = element_text(face = "bold")) +
     theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))+
  theme(legend.text=element_text(size=18))

important_df_PSPI <- important_df[important_df$Method=="STEPD",]
g1 <- ggdensity(data=important_df_PSPI,x="Density_Info",color="blue",ylab="Percentage of total TFS",
                xlab="Distribution of Jaccard Coefficient",ylim=c(0,28), title = "Jaccard Coefficient Distribution for STEPD Method", 
                font.label = list(size = 16, type="bold", color = "black"),fill="cyan") +
  theme(plot.title = element_text(size=20,hjust = 0.5)) + theme(plot.title = element_text(face = "bold")) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))+
  theme(legend.text=element_text(size=18))

important_df_MGM <- important_df[important_df$Method=="MGM",]
g2 <- ggdensity(data=important_df_MGM,x="Density_Info",color="purple",ylab="Percentage of total TFS",
                xlab="Distribution of Jaccard Coefficient",ylim=c(0,28), title = "Jaccard Coefficient Distribution for MGM Method", 
                font.label = list(size = 16, type="bold", color = "black"),fill="pink") +
  theme(plot.title = element_text(size=20,hjust = 0.5)) + theme(plot.title = element_text(face = "bold")) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))+
  theme(legend.text=element_text(size=18))


g <- ggarrange(g1,g2,
          ncol = 2, nrow = 1)
ggsave(filename = "Topology_Diff_PSPI_vs_MGM.pdf", plot = g, device = pdf(), width = 15, height = 6, units = "in", dpi = 300)     
dev.off()

#Load the graphs as edgelist
idh_mutant_X <- read.table("../../../Datasets/Real/TCGA/Gene_Mutation_A.csv",header=F)
idh_wildtype_Y <- read.table("../../../Datasets/Real/TCGA/Gene_NonMutation_B.csv",header=F)
gene_names_all <- read.csv("../../../Datasets/Real/TCGA/ColNames.csv",header=F)
gene_names_all <- as.character(as.vector(gene_names_all$V1))
idh_mutant_X$V1 <- gene_names_all[idh_mutant_X$V1]
idh_mutant_X$V2 <- gene_names_all[idh_mutant_X$V2]
rev_idh_mutant_X <- NULL
for (i in 1:nrow(idh_mutant_X))
{
  if (idh_mutant_X$V1[i]!=idh_mutant_X$V2[i])
  {
    rev_idh_mutant_X <- rbind(rev_idh_mutant_X,idh_mutant_X[i,])
  }
}
rev_idh_mutant_X <- as.data.frame(rev_idh_mutant_X)
colnames(rev_idh_mutant_X) <- c("TF","Target","Weight")
rev_idh_mutant_X$TF <- as.character(as.vector(rev_idh_mutant_X$TF))
rev_idh_mutant_X$Target <- as.character(as.vector(rev_idh_mutant_X$Target))
rev_idh_mutant_X$Weight <- as.numeric(as.vector(rev_idh_mutant_X$Weight))

idh_wildtype_Y$V1 <- gene_names_all[idh_wildtype_Y$V1]
idh_wildtype_Y$V2 <- gene_names_all[idh_wildtype_Y$V2]
rev_idh_wildtype_Y <- NULL
for (i in 1:nrow(idh_wildtype_Y))
{
  if (idh_wildtype_Y$V1[i]!=idh_wildtype_Y$V2[i])
  {
    rev_idh_wildtype_Y <- rbind(rev_idh_wildtype_Y,idh_wildtype_Y[i,])
  }
}
rev_idh_wildtype_Y <- as.data.frame(rev_idh_wildtype_Y)
colnames(rev_idh_wildtype_Y) <- c("TF","Target","Weight")
rev_idh_wildtype_Y$TF <- as.character(as.vector(rev_idh_wildtype_Y$TF))
rev_idh_wildtype_Y$Target <- as.character(as.vector(rev_idh_wildtype_Y$Target))
rev_idh_wildtype_Y$Weight <- as.numeric(as.vector(rev_idh_wildtype_Y$Weight))

shared_targets_pspi <- NULL
for (i in 1:nrow(bijection_pspi))
{
  if ((bijection_pspi$Gene_PSPI_X[i] %in% unique(rev_idh_mutant_X$TF)) &&
      (bijection_pspi$Gene_Y[i] %in% unique(rev_idh_wildtype_Y$TF)))
  {
    v1 <- rev_idh_mutant_X[rev_idh_mutant_X$TF==bijection_pspi$Gene_PSPI_X[i],]$Target
    v2 <- rev_idh_wildtype_Y[rev_idh_wildtype_Y$TF==bijection_pspi$Gene_Y[i],]$Target
    v <- intersect(v1,v2)
    shared_targets_pspi <- rbind(shared_targets_pspi,cbind(bijection_pspi$Gene_PSPI_X[i],
                                                           bijection_pspi$Gene_Y[i],
                                                           length(v1),length(v2),
                                                           length(v)))
  }
}
shared_targets_pspi <- as.data.frame(shared_targets_pspi)
colnames(shared_targets_pspi) <- c("TF_PSPI_X","TF_Y","Targets_PSPI_X","Targets_Y","Matched_Targets")
shared_targets_pspi$TF_PSPI_X <- as.character(as.vector(shared_targets_pspi$TF_PSPI_X))
shared_targets_pspi$TF_Y <- as.character(as.vector(shared_targets_pspi$TF_Y))
shared_targets_pspi$Targets_PSPI_X <- as.numeric(as.vector(shared_targets_pspi$Targets_PSPI_X))
shared_targets_pspi$Targets_Y <- as.numeric(as.vector(shared_targets_pspi$Targets_Y))
shared_targets_pspi$Matched_Targets <- as.numeric(as.vector(shared_targets_pspi$Matched_Targets))


edgelist_gsx1_x_and_e2f7_y <- NULL
x_names <- names(X["GSX1",which(X["GSX1",]>0)])
y_names <- names(Y["E2F7",which(Y["E2F7",]>0)])
edgelist_gsx1_x_and_e2f7_y <- cbind(rep("GSX1",length(x_names)),x_names)
edgelist_gsx1_x_and_e2f7_y <- rbind(edgelist_gsx1_x_and_e2f7_y,
                                    cbind(rep("E2F7",length(y_names)),y_names))
intersect_genes <- intersect(x_names,y_names)
label_vector <- c(rep(1,length(x_names)),rep(3,length(y_names)))
edgelist_gsx1_x_and_e2f7_y <- cbind(edgelist_gsx1_x_and_e2f7_y,label_vector)
edgelist_gsx1_x_and_e2f7_y <- as.data.frame(edgelist_gsx1_x_and_e2f7_y)
colnames(edgelist_gsx1_x_and_e2f7_y) <- c("TF","Target","Weight")
edgelist_gsx1_x_and_e2f7_y$TF <- as.character(as.vector(edgelist_gsx1_x_and_e2f7_y$TF))
edgelist_gsx1_x_and_e2f7_y$Target <- as.character(as.vector(edgelist_gsx1_x_and_e2f7_y$Target))
edgelist_gsx1_x_and_e2f7_y$Weight <- as.numeric(as.vector(edgelist_gsx1_x_and_e2f7_y$Weight))
edgelist_gsx1_x_and_e2f7_y[edgelist_gsx1_x_and_e2f7_y$Target %in% intersect_genes,]$Weight <- 2
write.table(edgelist_gsx1_x_and_e2f7_y,"GSX1_X_vs_E2F7_Y.csv",row.names = F,col.names = T,quote=F)

all_nodes_gsx1_e2f7 <- union(c(c("GSX1","E2F7"),x_names),y_names)
all_nodes_table <- cbind(all_nodes_gsx1_e2f7,rep(1,length(all_nodes_gsx1_e2f7)))
all_nodes_table <- as.data.frame(all_nodes_table)
colnames(all_nodes_table) <- c("Nodes","Color")
all_nodes_table$Nodes <- as.character(as.vector(all_nodes_table$Nodes))
all_nodes_table$Color <- as.numeric(as.vector(all_nodes_table$Color))
all_nodes_table[all_nodes_table$Nodes %in% x_names, ]$Color <- 2
all_nodes_table[all_nodes_table$Nodes %in% y_names, ]$Color <- 3
all_nodes_table[all_nodes_table$Nodes %in% intersect_genes, ]$Color <- 4
all_nodes_table[c(1:2),]$Color <- 1
write.table(all_nodes_table,"GSX1_E2F7_Node_Info.csv",row.names=F,col.names = T,quote=F)



edgelist_esr1_x_and_ventx_y <- NULL
x_names <- names(X["ESR1",which(X["ESR1",]>0)])
y_names <- names(Y["VENTX",which(Y["VENTX",]>0)])
edgelist_x_and_y <- cbind(rep("ESR1",length(x_names)),x_names)
edgelist_x_and_y <- rbind(edgelist_x_and_y,
                                    cbind(rep("VENTX",length(y_names)),y_names))
intersect_genes <- intersect(x_names,y_names)
label_vector <- c(rep(1,length(x_names)),rep(3,length(y_names)))
edgelist_x_and_y <- cbind(edgelist_x_and_y,label_vector)
edgelist_x_and_y <- as.data.frame(edgelist_x_and_y)
colnames(edgelist_x_and_y) <- c("TF","Target","Weight")
edgelist_x_and_y$TF <- as.character(as.vector(edgelist_x_and_y$TF))
edgelist_x_and_y$Target <- as.character(as.vector(edgelist_x_and_y$Target))
edgelist_x_and_y$Weight <- as.numeric(as.vector(edgelist_x_and_y$Weight))
edgelist_x_and_y[edgelist_x_and_y$Target %in% intersect_genes,]$Weight <- 2
write.table(edgelist_x_and_y,"ESR1_X_vs_VENTX_Y.csv",row.names = F,col.names = T,quote=F)

all_nodes_x_y <- union(c(c("ESR1","VENTX"),x_names),y_names)
all_nodes_table <- cbind(all_nodes_x_y,rep(1,length(all_nodes_x_y)))
all_nodes_table <- as.data.frame(all_nodes_table)
colnames(all_nodes_table) <- c("Nodes","Color")
all_nodes_table$Nodes <- as.character(as.vector(all_nodes_table$Nodes))
all_nodes_table$Color <- as.numeric(as.vector(all_nodes_table$Color))
all_nodes_table[all_nodes_table$Nodes %in% x_names, ]$Color <- 2
all_nodes_table[all_nodes_table$Nodes %in% y_names, ]$Color <- 3
all_nodes_table[all_nodes_table$Nodes %in% intersect_genes, ]$Color <- 4
all_nodes_table[c(1:2),]$Color <- 1
write.table(all_nodes_table,"ESR1_VENTX_Node_Info.csv",row.names=F,col.names = T,quote=F)


edgelist_x_and_y <- NULL
x_names <- names(X["IRF7",which(X["IRF7",]>0)])
y_names <- names(Y["PLAG1",which(Y["PLAG1",]>0)])
edgelist_x_and_y <- cbind(rep("IRF7",length(x_names)),x_names)
edgelist_x_and_y <- rbind(edgelist_x_and_y,
                          cbind(rep("PLAG1",length(y_names)),y_names))
intersect_genes <- intersect(x_names,y_names)
label_vector <- c(rep(1,length(x_names)),rep(3,length(y_names)))
edgelist_x_and_y <- cbind(edgelist_x_and_y,label_vector)
edgelist_x_and_y <- as.data.frame(edgelist_x_and_y)
colnames(edgelist_x_and_y) <- c("TF","Target","Weight")
edgelist_x_and_y$TF <- as.character(as.vector(edgelist_x_and_y$TF))
edgelist_x_and_y$Target <- as.character(as.vector(edgelist_x_and_y$Target))
edgelist_x_and_y$Weight <- as.numeric(as.vector(edgelist_x_and_y$Weight))
edgelist_x_and_y[edgelist_x_and_y$Target %in% intersect_genes,]$Weight <- 2
write.table(edgelist_x_and_y,"IRF7_X_vs_PLAG1_Y.csv",row.names = F,col.names = T,quote=F)

all_nodes_x_y <- union(c(c("IRF7","PLAG1"),x_names),y_names)
all_nodes_table <- cbind(all_nodes_x_y,rep(1,length(all_nodes_x_y)))
all_nodes_table <- as.data.frame(all_nodes_table)
colnames(all_nodes_table) <- c("Nodes","Color")
all_nodes_table$Nodes <- as.character(as.vector(all_nodes_table$Nodes))
all_nodes_table$Color <- as.numeric(as.vector(all_nodes_table$Color))
all_nodes_table[all_nodes_table$Nodes %in% x_names, ]$Color <- 2
all_nodes_table[all_nodes_table$Nodes %in% y_names, ]$Color <- 3
all_nodes_table[all_nodes_table$Nodes %in% intersect_genes, ]$Color <- 4
all_nodes_table[c(1:2),]$Color <- 1
write.table(all_nodes_table,"IRF7_PLAG1_Node_Info.csv",row.names=F,col.names = T,quote=F)
