#load packages
library(ggplot2, lib.loc = "/scratch/user/rbovio/R_LIBS_USER/3.6.2-intel-2019b-recommended-mt")
library(reshape2, lib.loc = "/scratch/user/rbovio/R_LIBS_USER/3.6.2-intel-2019b-recommended-mt")
library(DescTools, lib.loc = "/scratch/user/rbovio/R_LIBS_USER/3.6.2-intel-2019b-recommended-mt")
library(ggpubr, lib.loc = "/scratch/user/rbovio/R_LIBS_USER/3.6.2-intel-2019b-recommended-mt")
library(cowplot, lib.loc = "/scratch/user/rbovio/R_LIBS_USER/3.6.2-intel-2019b-recommended-mt")

#store command line arguements as variables
arrArgs <- commandArgs(trailingOnly = TRUE)
parameters<-as.character(arrArgs[1])
p1_af<-as.character(arrArgs[2])
p2_af<-as.character(arrArgs[3])

parameters_list = strsplit(parameters,"_")
linkage<- parameters_list[[1]][1]
gene_action = parameters_list[[1]][2]
esize<- parameters_list[[1]][3]
ss<- parameters_list[[1]][4]
ns<- parameters_list[[1]][5]


#create a dataframe to store the mean/sd/cv phenotype for each iteration/replicate
df = data.frame(NULL)

#create parameter counters
proportion_overlap_mean_p1 = 0
proportion_overlap_mean_p2 = 0
novelty = 0 #95% CI doesn't overlap with p1 or p2
num_failed_sims = 0

#loop through each directory (ie iteration)
for(i in 1:100){
  setwd(paste0("/scratch/user/rbovio/admixem3.0/admixem-master3.0/bin/F1_",parameters,"_",i))
  f = 'hybrids00'
  if(file.exists(f)==FALSE){
    num_failed_sims = num_failed_sims + 1
    next
  }
  #get parental and hybrid phenotype values
  p1m_data = read.delim(file = "parents00", header = TRUE)
  p2f_data = read.delim(file = "parents01", header = TRUE)
  hm_data = read.delim(file = "hybrids00", header = TRUE)
  hf_data = read.delim(file = "hybrids01", header = TRUE)
  
  
  #append to the dataframe
  df = rbind(df, c(i,
                   mean(p1m_data$Trait1),
                   mean(p2f_data$Trait1),
                   mean(c(hm_data$Trait1,hf_data$Trait1)), # EDIT HERE mean(ns$x)
                   sd(p1m_data$Trait1),
                   sd(p2f_data$Trait1),
                   sd(c(hm_data$Trait1,hf_data$Trait1)),
                   sd(p1m_data$Trait1)/mean(p1m_data$Trait1),
                   sd(p2f_data$Trait1)/mean(p2f_data$Trait1),
                   sd(c(hm_data$Trait1,hf_data$Trait1))/mean(c(hm_data$Trait1,hf_data$Trait1))))
  
  # Summary statistics for each replicate in order to determine proportion of replicates that show novel phenotypes
  p1_mean = mean(p1m_data$Trait1)
  p1_sd = sd(p1m_data$Trait1)
  p1_cv = p1_sd/p1_mean
  p1_n = nrow(p1m_data)
  p1_error <- qnorm(0.975)*p1_sd/sqrt(p1_n)
  p1_left  <- p1_mean - p1_error
  p1_right <- p1_mean + p1_error
  p1_sum_stat <- cbind(p1_left, p1_mean, p1_right)
  
  p2_mean = mean(p2f_data$Trait1)
  p2_sd = sd(p2f_data$Trait1)
  p2_cv = p2_sd/p2_mean
  p2_n = nrow(p2f_data)
  p2_error <- qnorm(0.975)*p2_sd/sqrt(p2_n)
  p2_left  <- p2_mean - p2_error
  p2_right <- p2_mean + p2_error
  p2_sum_stat <- cbind(p2_left, p2_mean, p2_right)
  
  h_mean = mean(c(hm_data$Trait1,hf_data$Trait1))
  h_sd = sd(c(hm_data$Trait1,hf_data$Trait1))
  h_cv = h_sd/h_mean
  h_n = nrow(hm_data)+nrow(hf_data)
  h_error <- qnorm(0.975)*h_sd/sqrt(h_n)
  h_left  <- h_mean - h_error
  h_right <- h_mean + h_error
  h_sum_stat <- cbind(h_left, h_mean, h_right)
  
  #check range overlap
  p1_range = c(p1_sum_stat[1,1],p1_sum_stat[1,3])
  p2_range = c(p2_sum_stat[1,1],p2_sum_stat[1,3])
  h_range = c(h_sum_stat[1,1],h_sum_stat[1,3])
  if(p1_range %overlaps% h_range == TRUE){
    proportion_overlap_mean_p1 = proportion_overlap_mean_p1 + 1
    # print(p1_range)
    # print(h_range)
  }
  if(p2_range %overlaps% h_range == TRUE){
    proportion_overlap_mean_p2 = proportion_overlap_mean_p2 + 1
    # print(p2_range)
    # print(h_range)
  }
  if(p1_range %overlaps% h_range == FALSE & p2_range %overlaps% h_range == FALSE){
    novelty = novelty + 1
  }
  colnames(df) = c("replicate", 
                   "p1_mean","p2_mean","h_mean",
                   "p1_sd","p2_sd","h_sd",
                   "p1_cv","p2_cv","h_cv")  
     
  setwd("/scratch/user/rbovio/admixem3.0/admixem-master3.0/bin")
}
#add all parameters to df dataframe for writing
df$linkage = linkage
df$gene_action =  gene_action
df$p1_af = p1_af
df$p2_af = p2_af
df$esize = esize
df$ss = ss
df$ns = ns
write.csv(df, paste0("df_",parameters,".csv"), quote = FALSE, row.names = FALSE)

#remove unneccesary columns for remaining analysis
df = df[,2:10]

p_o_mean = matrix(ncol = 3, nrow = 1)
p_o_mean[1,1] = proportion_overlap_mean_p1
p_o_mean[1,2] = proportion_overlap_mean_p2
p_o_mean[1,3] = novelty
write.table(p_o_mean, paste0("prop_overlap_mean_",parameters,".txt"), quote = FALSE, row.names = FALSE)
write.table(num_failed_sims, paste0("num_failed_sims_",parameters,".txt"), quote = FALSE, row.names = FALSE)

df_mean = df[,1:3]
colnames(df_mean) = c("p1","p2","h")
df_mean = melt(df_mean, measure.vars = c(1,2,3))
df_mean$variable <- factor(df_mean$variable, levels = c("p1", "h", "p2"))

df_sd = df[,4:6]
colnames(df_sd) = c("p1","p2","h")
df_sd = melt(df_sd, measure.vars = c(1,2,3))
df_sd$variable <- factor(df_sd$variable, levels = c("p1", "h", "p2"))

df_cv = df[,7:9]
colnames(df_cv) = c("p1","p2","h")
df_cv = melt(df_cv, measure.vars = c(1,2,3))
df_cv$variable <- factor(df_cv$variable, levels = c("p1", "h", "p2"))
df_cv[is.na(df_cv)] <- 0

#summary statistics over all replicates to get 95% CI
mean_vec <- colMeans(df)
sd_vec   <- apply(df, 2, sd)
n        <- nrow(df) #number of replicates

error <- qnorm(0.975)*sd_vec/sqrt(n)
left  <- mean_vec - error
right <- mean_vec + error

sum_stat <- cbind(left, mean_vec, right)

#mean
sum_stat_mean = sum_stat[1:3,]

sum_stat_mean = as.data.frame(sum_stat_mean)
sum_stat_mean$y = 0
sum_stat_mean$variable = factor(c("p1","p2","h"), levels = c("p1","p2","h"), ordered = TRUE)
sum_stat_mean$variable = factor(sum_stat_mean$variable, levels = c("p1","h","p2"))

l <- density(df_mean$value)
p <- ggplot(data =df_mean, aes(x=value, group = variable, fill = variable)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  scale_fill_manual(values=c("steelblue1", "springgreen1","tomato1")) +
  xlim(c(0,1)) +
  xlab("mean") +
  labs(fill = "population") +
  geom_point(inherit.aes = FALSE, data = sum_stat_mean, 
             aes(x= mean_vec, y = y, color = variable, size = 1.5),show.legend = FALSE)+
  geom_errorbarh(inherit.aes = FALSE, data = sum_stat_mean, 
                 aes(xmin = left, xmax = right, y = y, color = variable), 
                 height = 0.1, show.legend = FALSE, size = 1.5) +
  scale_color_manual(values=c( "steelblue4","springgreen4", "tomato4")) +
  theme_bw()

# #sd
sum_stat_sd = sum_stat[4:6,]

sum_stat_sd = as.data.frame(sum_stat_sd)
sum_stat_sd$y = 0
sum_stat_sd$variable = factor(c("p1","p2","h"), levels = c("p1","p2","h"), ordered = TRUE)
sum_stat_sd$variable = factor(sum_stat_sd$variable, levels = c("p1","h","p2"))

l <- density(df_sd$value)
pp <- ggplot(data =df_sd, aes(x=value, group = variable, fill = variable)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values=c("steelblue1", "springgreen1","tomato1")) +
  xlim(range(l$x)) +
  xlab("standard deviation") +
  labs(fill = "population") +
  geom_point(inherit.aes = FALSE, data = sum_stat_sd,
             aes(x= mean_vec, y = y, color = variable, size = 1.5),show.legend = FALSE)+
  geom_errorbarh(inherit.aes = FALSE, data = sum_stat_sd,
                 aes(xmin = left, xmax = right, y = y, color = variable),
                 height = 0.1, show.legend = FALSE, size = 1.5) +
  scale_color_manual(values=c( "steelblue4", "springgreen4","tomato4")) +
  theme_bw()

#cv
sum_stat_cv = sum_stat[7:9,]

sum_stat_cv = as.data.frame(sum_stat_cv)
sum_stat_cv$y = 0
sum_stat_cv$variable = factor(c("p1","p2","h"), levels = c("p1","p2","h"), ordered = TRUE)
sum_stat_cv$variable = factor(sum_stat_cv$variable, levels = c("p1","h","p2"))

l <- density(df_cv$value)
ppp <- ggplot(data =df_cv, aes(x=value, group = variable, fill = variable)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  scale_fill_manual(values=c("steelblue1", "springgreen1","tomato1")) +
  xlim(range(l$x)) +
  xlab("coefficient of variation") +
  labs(fill = "population") +
  geom_point(inherit.aes = FALSE, data = sum_stat_cv, 
             aes(x= mean_vec, y = y, color = variable, size = 1.5),show.legend = FALSE)+
  geom_errorbarh(inherit.aes = FALSE, data = sum_stat_cv, 
                 aes(xmin = left, xmax = right, y = y, color = variable), 
                 height = 0.1, show.legend = FALSE, size = 1.5) +
  scale_color_manual(values=c( "steelblue4", "springgreen4","tomato4")) +
  theme_bw()

#PLOT
legend <- cowplot::get_legend(pp)
# text <- paste0("p1 = ",p1_af,", p2 =",p2_af)
# # Create a text grob
# tgrob <- text_grob(text,size = 20)
# # Draw the text
# plot_0 <- as_ggplot(tgrob) + theme(plot.margin = margin(0,3,0,0, "cm"))
# save plot as jpeg
options(bitmapType='cairo')
jpeg(file=paste0("variation_mean_and_cv_density_",parameters,"_p1-",p1_af,"_p2-",p2_af,".jpeg"), width = 1000, height = 600)
# ggarrange(plot_0, NULL, NULL,p,ppp,legend,nrow=2,ncol=3, widths = c(3,3,1), heights = c(1,10))
ggarrange(p,ppp,legend,nrow=1,ncol=3, widths = c(3,3,1))
dev.off()

