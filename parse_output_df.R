dir = c("output_unlinked_adddom_esize0.1_noss_nons",
        "output_unlinked_adddom_esize0.1_noss_ns",
        "output_unlinked_adddom_esize0.1_ssds_nons",
        "output_unlinked_adddom_esize0.1_ssfd_nons",
        "output_unlinked_adddom_esize0.9_noss_nons",
        "output_unlinked_adddom_esize0.9_noss_ns",
        "output_unlinked_adddom_esize0.9_ssds_nons",
        "output_unlinked_adddom_esize0.9_ssfd_nons",
        "output_unlinked_additive_esize0.1_noss_nons",
        "output_unlinked_additive_esize0.1_noss_ns",
        "output_unlinked_additive_esize0.1_ssds_nons",
        "output_unlinked_additive_esize0.1_ssfd_nons",
        "output_unlinked_additive_esize0.9_noss_nons",
        "output_unlinked_additive_esize0.9_noss_ns",
        "output_unlinked_additive_esize0.9_ssds_nons",
        "output_unlinked_additive_esize0.9_ssfd_nons",
        "output_unlinked_dominant_esize0.1_noss_nons",
        "output_unlinked_dominant_esize0.1_noss_ns",
        "output_unlinked_dominant_esize0.1_ssds_nons",
        "output_unlinked_dominant_esize0.1_ssfd_nons",
        "output_unlinked_dominant_esize0.9_noss_nons",
        "output_unlinked_dominant_esize0.9_noss_ns",
        "output_unlinked_dominant_esize0.9_ssds_nons",
        "output_unlinked_dominant_esize0.9_ssfd_nons",
        "output_linked_adddom_esize0.1_noss_nons",
        "output_linked_adddom_esize0.1_noss_ns",
        "output_linked_adddom_esize0.1_ssds_nons",
        "output_linked_adddom_esize0.1_ssfd_nons",
        "output_linked_adddom_esize0.9_noss_nons",
        "output_linked_adddom_esize0.9_noss_ns",
        "output_linked_adddom_esize0.9_ssds_nons",
        "output_linked_adddom_esize0.9_ssfd_nons",
        "output_linked_additive_esize0.1_noss_nons",
        "output_linked_additive_esize0.1_noss_ns",
        "output_linked_additive_esize0.1_ssds_nons",
        "output_linked_additive_esize0.1_ssfd_nons",
        "output_linked_additive_esize0.9_noss_nons",
        "output_linked_additive_esize0.9_noss_ns",
        "output_linked_additive_esize0.9_ssds_nons",
        "output_linked_additive_esize0.9_ssfd_nons",
        "output_linked_dominant_esize0.1_noss_nons",
        "output_linked_dominant_esize0.1_noss_ns",
        "output_linked_dominant_esize0.1_ssds_nons",
        "output_linked_dominant_esize0.1_ssfd_nons",
        "output_linked_dominant_esize0.9_noss_nons",
        "output_linked_dominant_esize0.9_noss_ns",
        "output_linked_dominant_esize0.9_ssds_nons",
        "output_linked_dominant_esize0.9_ssfd_nons")

#create blank dataframe what we'll combine all the other dataframes to. 
df = data.frame(NULL)

for(i in dir){
  #loop into each output directory
  setwd(paste0("/scratch/user/rbovio/admixem3.0/admixem-master3.0/bin/",i))
  
  parameters_list = strsplit(i,"_")
  linkage<- parameters_list[[1]][2]
  gene_action = parameters_list[[1]][3]
  esize<- parameters_list[[1]][4]
  ss<- parameters_list[[1]][5]
  ns<- parameters_list[[1]][6]
  
  #check to see if you're in a natural selection folder
  if(ns == "ns"){
    #load dataframe
    tmp1 = read.csv(file = paste0("df_",linkage,"_",gene_action,"_",esize,"_",ss,"_nsfd_all.csv"))
    tmp2 = read.csv(file = paste0("df_",linkage,"_",gene_action,"_",esize,"_",ss,"_nsnfd_all.csv"))
    tmp3 = read.csv(file = paste0("df_",linkage,"_",gene_action,"_",esize,"_",ss,"_nsds_all.csv"))
    #append to combined dataframe
    df = rbind(df, tmp1, tmp2, tmp3)
    setwd("/scratch/user/rbovio/admixem3.0/admixem-master3.0/bin/")
    next
  }
  else{
    #load dataframe
    tmp = read.csv(file = paste0("df_",linkage,"_",gene_action,"_",esize,"_",ss,"_",ns,"_all.csv"))
    #append to combined dataframe
    df = rbind(df, tmp)
    setwd("/scratch/user/rbovio/admixem3.0/admixem-master3.0/bin/")
  }
}

write.csv(df, "all_parameters.csv", quote = FALSE, row.names = FALSE)
