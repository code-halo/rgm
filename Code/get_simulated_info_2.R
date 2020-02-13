library(data.table)
library(stringr)
library(plotly)

setwd('.')

get_pr_metrics <- function(folder_name,experimentid)
{
  #Show the precision recall metric
  pr_info <- read.table(paste0('results/',folder_name,'/',experimentid,'_PR.info'),header=FALSE,sep=" ");
  pspi_info <- data.frame(str_split_fixed(pr_info$V1,fixed("_"),9));
  pspi_info <- data.frame(pspi_info[,c(1:7)],str_split_fixed(pspi_info$X8,fixed(".txt"),2))
  colnames(pspi_info) <- c(paste0("X",c(1:ncol(pspi_info))))
  pspi_info$X4 <- as.numeric(as.vector(pspi_info$X6));
  pspi_info$X6 <- as.numeric(as.vector(pspi_info$X8))
  
  nu_range <- c(0.03125,0.0625,0.125,0.25,0.5,1);
  mu_range <- c(0.03125,0.0625,0.125,0.25,0.5,1)
  
  avg_precision <- NULL;
  avg_recall <- NULL;
  avg_fscore <- NULL;
  metric_df <- NULL;
  for (i in 1:length(nu_range))
  {
    for(j in 1:length(mu_range))
    {
      nu = nu_range[i];
      mu = mu_range[j];
      id = which(pspi_info$X4==nu & pspi_info$X6==mu);
      if (length(id)>0)
      {
        val_prec <- pr_info[id,]$V9;
        val_rec <- pr_info[id,]$V11
        val_fmeasure <- (2*val_prec*val_rec)/(val_rec+val_prec+.Machine$double.eps);
      }
      else{
        val_prec <- c(0,0,0)
        val_rec <- c(0,0,0)
        val_fmeasure <- c(0,0,0)
      }
      metric_df <- rbind(metric_df,cbind(rep(nu,3),rep(mu,3),val_prec,val_rec,val_fmeasure,c(1,2,3)));
    }
  }
  
  metric_df <- as.data.frame(metric_df);
  colnames(metric_df) <- c("Nu","Mu","Precision","Recall","Fscore","Exp_Id")
  metric_df$Nu <- as.numeric(as.vector(metric_df$Nu));
  metric_df$Mu <- as.numeric(as.vector(metric_df$Mu));
  metric_df[,3] <- as.numeric(as.vector(metric_df[,3]));
  metric_df[,4] <- as.numeric(as.vector(metric_df[,4]));
  metric_df[,5] <- as.factor(as.vector(metric_df[,5]));
  nu_mu_tuple <- paste0(as.character(as.vector(metric_df$Nu)),",",as.character(as.vector(metric_df$Mu)));
  metric_df$nu_mu_tuple <- nu_mu_tuple;
  metric_df2 <- melt(metric_df[,3:ncol(metric_df)],id.vars = c("nu_mu_tuple","Exp_Id"))
  return(metric_df2)
}

get_qualitative_metrics <- function(folder_name,experimentid)
{
  fmax_info <- read.table(paste0('results/',folder_name,'/',experimentid,'_Fmax_inferred.info'),header=FALSE,sep=" ");
  fmax_true_info <- read.table(paste0('results/',folder_name,'/',experimentid,'_Fmax_gd.info'),header=FALSE,sep=" ");
  me_info <- read.table(paste0('results/',folder_name,'/',experimentid,'_ME.info'),header=FALSE,sep=" ")
  
  fmax_info[,1] <- as.character(as.vector(fmax_info[,1]))
  fmax_true_info[,1] <- as.character(as.vector(fmax_true_info[,1]))
  me_info[,1] <- as.character(as.vector(me_info[,1]))
  
  newdf <- data.frame(V1=fmax_info$V1,V2=fmax_info$V4,V3=fmax_true_info$V4,V4=fmax_info$V8,V5=fmax_true_info$V8)
  newdf[,1] <- as.character(as.vector(newdf[,1]))
  
  nu_range <- c("0.031250","0.062500","0.125000","0.250000","0.500000","1.000000");
  mu_range <- c("0.031250","0.062500","0.125000","0.250000","0.500000","1.000000");
  
  fmax_metric_df <- NULL;
  me_metric_df <- NULL
  for (i in 1:length(nu_range))
  {
    for(j in 1:length(mu_range))
    {
      nu = nu_range[i];
      mu = mu_range[j];
      fmax_id = which(grepl(pattern=paste0("nu_",nu,"_mu_",mu),x=newdf$V1));
      me_id <- which(grepl(pattern=paste0("nu_",nu,"_mu_",mu),x=me_info$V1))
      if (length(fmax_id)==3)
      {
        val_infer_fmax <- newdf[fmax_id,]$V3
        val_diff_fmax <- abs(newdf[fmax_id,]$V3-newdf[fmax_id,]$V2)
        val_diff_k <- abs(newdf[fmax_id,]$V5-newdf[fmax_id,]$V4)
        val_opt_k <- newdf[fmax_id,]$V5
        val_infer_k <- newdf[fmax_id,]$V4
      }
      else{
        val_infer_fmax <- c(0,0,0)
        val_diff_fmax <- c(1,1,1)
        val_diff_k <- unique(newdf$V5)
        val_opt_k <- unique(newdf$V5)
        val_infer_k <- c(0,0,0)
      }
      val_me <- me_info[me_id,]$V4
      fmax_metric_df <- rbind(fmax_metric_df,cbind(rep(nu,3),rep(mu,3),val_diff_fmax,val_diff_k,val_infer_fmax,val_opt_k,
                                                   val_infer_k,c("Exp_1","Exp_2","Exp_3")));
      me_metric_df <- rbind(me_metric_df,cbind(rep(nu,3),rep(mu,3),val_me,c("Exp_1","Exp_2","Exp_3")))
    }
  }
  fmax_metric_df <- as.data.frame(fmax_metric_df);
  colnames(fmax_metric_df) <- c("Nu","Mu","Diff_Fmeasure","Diff_Opt_k","Fmeasure_Infer","Opt k","Infer k",experimentid)
  fmax_metric_df$Nu <- as.numeric(as.vector(fmax_metric_df$Nu));
  fmax_metric_df$Mu <- as.numeric(as.vector(fmax_metric_df$Mu));
  fmax_metric_df[,3] <- as.numeric(as.vector(fmax_metric_df[,3]));
  fmax_metric_df[,4] <- as.numeric(as.vector(fmax_metric_df[,4]));
  fmax_metric_df[,5] <- as.numeric(as.vector(fmax_metric_df[,5]))
  fmax_metric_df[,6] <- as.numeric(as.vector(fmax_metric_df[,6]))
  fmax_metric_df[,7] <- as.numeric(as.vector(fmax_metric_df[,7]))
  fmax_metric_df[,8] <- as.factor(as.vector(fmax_metric_df[,8]))
  
  nu_mu_tuple <- paste0(as.character(as.vector(fmax_metric_df$Nu)),",",as.character(as.vector(fmax_metric_df$Mu)));
  fmax_metric_df$nu_mu_tuple <- nu_mu_tuple;
  fmax_metric_df2 <- melt(fmax_metric_df[,3:ncol(fmax_metric_df)],id.vars = c("nu_mu_tuple",experimentid))
  
  me_metric_df <- as.data.frame(me_metric_df);
  colnames(me_metric_df) <- c("Nu","Mu","Matching_Error",experimentid);
  me_metric_df$Nu <- as.numeric(as.vector(me_metric_df$Nu));
  me_metric_df$Mu <- as.numeric(as.vector(me_metric_df$Mu));
  me_metric_df[,3] <- as.numeric(as.vector(me_metric_df[,3]));
  me_metric_df[,4] <- as.factor(as.vector(me_metric_df[,4]));
  nu_mu_tuple <- paste0(as.character(as.vector(me_metric_df$Nu)),",",as.character(as.vector(me_metric_df$Mu)));
  me_metric_df$nu_mu_tuple <- nu_mu_tuple;
  me_metric_df2 <- melt(me_metric_df[,3:ncol(me_metric_df)],id.vars = c("nu_mu_tuple",experimentid))
  
  return(list(fmax_metric_df2,me_metric_df2))
}

#EXPERIMENT FOR Precision-Recall
df_exp1 <- get_pr_metrics("Exp_1_Exhaustive","Exp_1");
df_exp2 <- get_pr_metrics("Exp_2_Exhaustive","Exp_2");
df_exp3 <- get_pr_metrics("Exp_3_Exhaustive","Exp_3");
df_exp4 <- get_pr_metrics("Exp_4_Exhaustive","Exp_4");
df_exp5 <- get_pr_metrics("Exp_5_Exhaustive","Exp_5");
 
out_value_df2 <- as.data.frame(cbind(df_exp1$value,df_exp2$value,df_exp3$value,df_exp4$value,df_exp5$value))
for (i in 1:ncol(out_value_df2))
{
   out_value_df2[,i] <- as.numeric(as.vector(out_value_df2[,i]))
}
value <- apply(out_value_df2,1,function(x) {mean(x,na.rm = TRUE)})
final_metric_df2 <- as.data.frame(cbind(df_exp1[,c(1:(ncol(df_exp1)-1))],value))
colnames(final_metric_df2) <- c("nu_mu_tuple","Exp_Id","Variable","Value")
final_metric_df2$nu_mu_tuple <- as.character(as.vector(final_metric_df2$nu_mu_tuple))
final_metric_df2$Exp_Id <- as.factor(as.vector(final_metric_df2$Exp_Id))
final_metric_df2$Variable <- as.factor(as.vector(final_metric_df2$Variable))
final_metric_df2$Value <- as.numeric(as.vector(final_metric_df2$Value))
 
g <- ggplot(data=final_metric_df2, aes(x=nu_mu_tuple, y=Value, shape=Variable, color=Exp_Id)) +
      geom_point(aes(shape=Variable),size=3) + 
   xlab("Parameter Tuples (Nu,Mu)") + ylab("Metric Values") + theme_bw() + ylim(0,1) +
   theme(text = element_text(size=18), axis.text.x = element_text(angle = -90, hjust = 1)) + 
   ggtitle("Graph Discrimination experiments") + 
   theme(plot.title = element_text(hjust = 0.5))

#Experiment for Unsupervised metric
out_exp1 <- get_qualitative_metrics("Exp_1_Exhaustive","Exp_1");
fmax_df_exp1 <- out_exp1[[1]]
me_df_exp1 <- out_exp1[[2]]

out_exp2 <- get_qualitative_metrics("Exp_2_Exhaustive","Exp_2");
fmax_df_exp2 <- out_exp2[[1]]
me_df_exp2 <- out_exp2[[2]]

out_exp3 <- get_qualitative_metrics("Exp_3_Exhaustive","Exp_3");
fmax_df_exp3 <- out_exp3[[1]]
me_df_exp3 <- out_exp3[[2]]

out_exp4 <- get_qualitative_metrics("Exp_4_Exhaustive","Exp_4");
fmax_df_exp4 <- out_exp4[[1]]
me_df_exp4 <- out_exp4[[2]]

out_exp5 <- get_qualitative_metrics("Exp_5_Exhaustive","Exp_5");
fmax_df_exp5 <- out_exp5[[1]]
me_df_exp5 <- out_exp5[[2]]

out_value_df <- as.data.frame(cbind(fmax_df_exp1$value,fmax_df_exp2$value,fmax_df_exp3$value,fmax_df_exp4$value,fmax_df_exp5$value))
for (i in 1:ncol(out_value_df))
{
   out_value_df[,i] <- as.numeric(as.vector(out_value_df[,i]))
}
value <- apply(out_value_df,1,function(x) {mean(x,na.rm = TRUE)})
final_metric_df <- as.data.frame(cbind(fmax_df_exp1[,c(1:(ncol(fmax_df_exp1)-1))],value))
colnames(final_metric_df) <- c("nu_mu_tuple","Exp_Id","variable","value")
final_metric_df$nu_mu_tuple <- as.character(as.vector(final_metric_df$nu_mu_tuple))
final_metric_df$Exp_Id <- as.factor(as.vector(final_metric_df$Exp_Id))
final_metric_df$variable <- as.factor(as.vector(final_metric_df$variable))
final_metric_df$value <- as.numeric(as.vector(final_metric_df$value))

p <- ggplot(data=final_metric_df, aes(x=nu_mu_tuple, y=value, group=variable, color=Exp_Id)) +
  geom_point(aes(shape=variable),size=3) +   scale_shape_manual(values=c(2, 3, 4, 5, 6))+ scale_color_manual(values=c("red","blue","green"))+
  xlab("Parameter Tuples (Nu,Mu)") + ylab("Metric Values") + theme_bw() + ylim(0,7) +
  theme(text = element_text(size=15),axis.text.x = element_text(angle = -90, hjust = 1)) + 
  ggtitle("Graph Discrimination Experiments") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0('results/Plots/Qualitative_Comparison_Graph_Discrimination_Task.pdf'), plot = p, height=10, width=15)

me_value_df <- as.data.frame(cbind(me_df_exp1$value,me_df_exp2$value,me_df_exp3$value,me_df_exp4$value,me_df_exp5$value))
for (i in 1:ncol(out_value_df))
{
  me_value_df[,i] <- as.numeric(as.vector(me_value_df[,i]))
}
value <- apply(me_value_df,1,function(x) {mean(x,na.rm = TRUE)})
final_metric_df3 <- as.data.frame(cbind(me_df_exp1[,c(1:(ncol(me_df_exp1)-1))],value))
colnames(final_metric_df3) <- c("nu_mu_tuple","Exp_Id","variable","value")
final_metric_df3$nu_mu_tuple <- as.character(as.vector(final_metric_df3$nu_mu_tuple))
final_metric_df3$Exp_Id <- as.factor(as.vector(final_metric_df3$Exp_Id))
final_metric_df3$variable <- as.factor(as.vector(final_metric_df3$variable))
final_metric_df3$value <- as.numeric(as.vector(final_metric_df3$value))


p1 <- ggplot(data=final_metric_df3, aes(x=nu_mu_tuple, y=value, group=Exp_Id, color=variable)) +
  geom_point(aes(shape=Exp_Id),size=3) +   scale_shape_manual(values=c(2, 3, 4))+ scale_color_manual(values=c("red","blue","green"))+
  xlab("Parameter Tuples (Nu,Mu)") + ylab("Metric Values") + theme_bw() + ylim(17,21) +
  theme(text = element_text(size=15),axis.text.x = element_text(angle = -90, hjust = 1)) + 
  ggtitle("Graph Discrimination Experiments") + 
  theme(plot.title = element_text(hjust = 0.5))
