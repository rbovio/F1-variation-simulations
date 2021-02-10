library(ggplot2)
library(reshape2)

arrArgs <- commandArgs(trailingOnly = TRUE)
parameters<-as.character(arrArgs[1])

path = getwd()
data<-data.frame(NULL)
file.names <- dir(path, pattern =paste0("prop_overlap_mean_",parameters))
for(i in 1:length(file.names)){
  file <- read.table(file.names[i],header=TRUE, sep=" ")
  data <- rbind(data, file)
}

data$initial_af_p1 = seq(1,0.5,-0.1)
data$initial_af_p2 = seq(0,0.5,0.1)
colnames(data) = c('p1-hybrid overlap', "p2-hybrid overlap",'novelty', 'initial_af_p1', "initial_af_p2")
data = melt(data, id.vars = c('initial_af_p1', "initial_af_p2"))
data$value = data$value #if you do 1000 replicates, divide by 10 to make it a proportion between 0-100.

p = ggplot(data, aes(x=interaction(initial_af_p1, initial_af_p2), y = value, group = variable, color = variable)) +
    theme_bw() +
    geom_point() +
    geom_line() + 
    ggtitle("proportion of replicates where the mean phenotype between\nhybrid and parental populations overlap") +
    annotate(geom = "text", x = seq(1,6), y = -10, label = unique(data$initial_af_p1), size = 3) +
    annotate(geom = "text", x = seq(1,6), y = -15, label = unique(data$initial_af_p2), size = 3) +
    annotate(geom = "text", x = 7, y = -10, label = "p1 initial allele frequency", size = 3, hjust = 0) +
    annotate(geom = "text", x = 7, y = -15, label = "p2 initial allele frequency", size = 3, hjust = 0) +
    coord_cartesian(ylim = c(-5, 105), expand = FALSE, clip = "off") +

    theme(plot.margin = unit(c(1, 1, 4, 1), "lines"),
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())



options(bitmapType='cairo')
jpeg(file=paste0("proportion_overlap_mean_",parameters,".jpeg"))
p
dev.off()



