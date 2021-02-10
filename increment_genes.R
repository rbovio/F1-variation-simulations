arrArgs <- commandArgs(trailingOnly = TRUE)
genes_file<-as.character(arrArgs[1])
increment_value<-as.numeric(arrArgs[2])

my_func = function(increment_value){
genes = read.delim(paste(genes_file), sep = "\t")
genes$Dominant.Freq.Pop1[2:11] = genes$Dominant.Freq.Pop1[2:11] - increment_value
genes$Dominant.Freq.Pop2[2:11] = genes$Dominant.Freq.Pop2[2:11] + increment_value
write.table(genes, paste(genes_file), sep = '\t', quote = FALSE, row.names = FALSE)
}

my_func(increment_value)
