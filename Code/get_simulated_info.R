library(data.table)
library(stringr)
library(plotly)

setwd('.')

get_pr_metrics <- function(folder_name,experimentid)
{
  #Show the precision recall metric
  pr_info <- read.table(paste0('results/',folder_name,'/PR.info'),header=FALSE,sep=" ");
  pspi_info <- data.frame(str_split_fixed(pr_info$V1,fixed("_"),7));
  pspi_info$X4 <- as.numeric(as.vector(pspi_info$X4));
  pspi_info$X6 <- as.numeric(as.vector(pspi_info$X6))
  
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
        val_prec <- mean(pr_info[id,]$V8);
        avg_precision <- c(avg_precision,val_prec);
        val_rec <- mean(pr_info[id,]$V10);
        avg_recall <- c(avg_recall,val_rec);
      }
      else{
        val_prec <- 0;
        val_rec <- 0;
        avg_precision <- c(avg_precision,0);
        avg_recall <- c(avg_recall,0);
      }
      val_fmeasure <- (2*val_prec*val_rec)/(val_rec+val_prec+.Machine$double.eps);
      avg_fscore <- c(avg_fscore,val_fmeasure);
      metric_df <- rbind(metric_df,cbind(nu,mu,val_prec,val_rec,val_fmeasure));
    }
  }
  
  metric_df <- as.data.frame(metric_df);
  colnames(metric_df) <- c("Nu","Mu",paste0("Precision_",experimentid),paste0("Recall_",experimentid),
                           paste0("Fscore_",experimentid));
  metric_df$Nu <- as.numeric(as.vector(metric_df$Nu));
  metric_df$Mu <- as.numeric(as.vector(metric_df$Mu));
  metric_df[,3] <- as.numeric(as.vector(metric_df[,3]));
  metric_df[,4] <- as.numeric(as.vector(metric_df[,4]));
  metric_df[,5] <- as.numeric(as.vector(metric_df[,5]));
  nu_mu_tuple <- paste0(as.character(as.vector(metric_df$Nu)),",",as.character(as.vector(metric_df$Mu)));
  metric_df$nu_mu_tuple <- nu_mu_tuple;
  metric_df2 <- melt(metric_df[,3:ncol(metric_df)],id.vars = "nu_mu_tuple")
  return(metric_df2)
}

get_fmax_metrics <- function(folder_name,experimentid)
{
  fmax_info <- read.table(paste0('results/',folder_name,'/Fmax_inferred.info'),header=FALSE,sep=" ");
  fmax_true_info <- read.table(paste0('results/',folder_name,'/Fmax_gd.info'),header=FALSE,sep=" ");
  pspi_info <- data.frame(str_split_fixed(fmax_info$V1,fixed("_"),7));
  pspi_info$X4 <- as.numeric(as.vector(pspi_info$X4));
  pspi_info$X6 <- as.numeric(as.vector(pspi_info$X6))
  
  nu_range <- c(0.03125,0.0625,0.125,0.25,0.5,1);
  mu_range <- c(0.03125,0.0625,0.125,0.25,0.5,1)
  
  metric_df <- NULL;
  for (i in 1:length(nu_range))
  {
    for(j in 1:length(mu_range))
    {
      nu = nu_range[i];
      mu = mu_range[j];
      id = which(pspi_info$X4==nu & pspi_info$X6==mu);
      val_fmax_true <- mean(fmax_true_info[id,]$V3);
      if (length(id)>0)
      {
        val_fmax <- mean(fmax_info[id,]$V3);
      }
      else{
        val_fmax <- 0;
      }
      metric_df <- rbind(metric_df,cbind(nu,mu,val_fmax,val_fmax_true));
    }
  }
  metric_df <- as.data.frame(metric_df);
  colnames(metric_df) <- c("Nu","Mu",paste0("Fmeasure(I)_",experimentid),paste0("Fmeasure(O)_",experimentid));
  metric_df$Nu <- as.numeric(as.vector(metric_df$Nu));
  metric_df$Mu <- as.numeric(as.vector(metric_df$Mu));
  metric_df[,3] <- as.numeric(as.vector(metric_df[,3]));
  metric_df[,4] <- as.numeric(as.vector(metric_df[,4]));
  nu_mu_tuple <- paste0(as.character(as.vector(metric_df$Nu)),",",as.character(as.vector(metric_df$Mu)));
  metric_df$nu_mu_tuple <- nu_mu_tuple;
  metric_df2 <- melt(metric_df[,3:ncol(metric_df)],id.vars = "nu_mu_tuple")
  return(metric_df2)
}

#EXPERIMENT FOR Precision-Recall
df_exp1 <- get_pr_metrics("r015_N_500_m1_50_m2_100_noise_30","ER_NP_30");
df_exp2 <- get_pr_metrics("r015_N_500_m1_50_m2_100_noise_30_both","ER_P_30");
df_exp3 <- get_pr_metrics("r015_N_500_m1_50_m2_100_noise_50","ER_NP_50");
df_exp4 <- get_pr_metrics("r015_N_500_m1_50_m2_100_noise_50_both","ER_P_50");

metric_df_exp1 <- as.data.frame(rbind(df_exp1,df_exp2,df_exp3,df_exp4));
metric_df_exp1[,1] <- as.character(as.vector(metric_df_exp1[,1]));
metric_df_exp1[,2] <- as.character(as.vector(metric_df_exp1[,2]));
metric_df_exp1[,3] <- as.numeric(as.vector(metric_df_exp1[,3]));

val <- c('Fscore','Fscore','Fscore','Fscore','Precision','Precision','Precision','Precision',
         'Recall','Recall','Recall','Recall')
val <- val[as.factor(metric_df_exp1$variable)];
metric_df_exp1$shape <- as.factor(val);
metric_df_exp1 <- metric_df_exp1[order(metric_df_exp1$nu_mu_tuple),];

g <- ggplot(data=metric_df_exp1, aes(x=nu_mu_tuple, y=value, group=variable, color=variable)) +
     geom_line(aes(linetype=variable),lwd=0.75) + geom_point(aes(shape=shape),size=4) + 
     xlab("Parameter Tuples (Nu,Mu)") + ylab("Metric Values") + theme_bw() + ylim(0,1) +
     scale_color_manual(values=c("red", "pink", "red", "pink", "blue","cyan", "blue", "cyan", "green","darkgreen",
                                 "green","darkgreen"))+
     theme(text = element_text(size=18), axis.text.x = element_text(angle = -90, hjust = 1)) + 
     ggtitle("Effect of Model Parameters on ER graph experiments") + 
     theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0('results/Plots/PRF_r015_Score.pdf'), plot = g, height=6, width=12)

#Experiment for Unsupervised metric
df2_exp1 <- get_fmax_metrics("a015_N_500_m1_50_m2_100_noise_30","Exp1");
df2_exp2 <- get_fmax_metrics("a015_N_500_m1_50_m2_100_noise_30_both","Exp2");
df2_exp3 <- get_fmax_metrics("a015_N_500_m1_50_m2_100_noise_50","Exp3");
df2_exp4 <- get_fmax_metrics("a015_N_500_m1_50_m2_100_noise_50_both","Exp4");

metric_df2_exp1 <- as.data.frame(rbind(df2_exp1,df2_exp2,df2_exp3,df2_exp4));
metric_df2_exp1[,1] <- as.character(as.vector(metric_df2_exp1[,1]));
metric_df2_exp1[,2] <- as.character(as.vector(metric_df2_exp1[,2]));
metric_df2_exp1[,3] <- as.numeric(as.vector(metric_df2_exp1[,3]));

val <- c('Fmeasure(I)','Fmeasure(I)','Fmeasure(I)','Fmeasure(I)',
         'Fmeasure(O)','Fmeasure(O)','Fmeasure(O)','Fmeasure(O)');
val <- val[as.factor(metric_df2_exp1$variable)];
metric_df2_exp1$shape <- as.factor(val);
metric_df2_exp1 <- metric_df2_exp1[order(metric_df2_exp1$nu_mu_tuple),];

p <- ggplot(data=metric_df2_exp1, aes(x=nu_mu_tuple, y=value, group=variable, color=variable)) +
  geom_line(aes(linetype=variable),size=1.5) + geom_point(aes(shape=shape),size=3) + 
  xlab("Parameter Tuples (Nu,Mu)") + ylab("Metric Values") + theme_bw() + ylim(0,0.5) +
  scale_color_manual(values=c("red", "green", "blue", "orange", "pink","darkgreen","cyan", "yellow")) + 
  theme(text = element_text(size=15),axis.text.x = element_text(angle = -90, hjust = 1)) + 
  ggtitle("Effect of Model Parameters on ER graph experiments") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0('results/Plots/Fmeasure_a015_Score.pdf'), plot = p, height=10, width=15)