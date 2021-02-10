# F1-variation-simulations

This project uses Admixem to simulate F1 hybrids under various genetic architectures and selection parameters and compares the coefficient of variation to that of the parental species. See: https://github.com/melop/admixem for instructions on installing and running Admixem. 
Next, a decision tree (https://www.guru99.com/r-decision-trees.html) is used to predict which genetic architectures and selection parameters are more likely to cause F1 hybrids to have a smaller coefficient of variaton than the parental species.

The files in this repository only simulate a single set of genetic architectures. Specifically, one in which the phenotype of interest is determined by 10 additive unlinked loci, distributed among 10 indepedent chromosomes, each with an effect size of 0.1 and no natural or sexual selection.  
See parameters_simulated.xlsx for the various parameter combinations that were simulated and F1_variation_markdown.Rmd for results. 

#Begin by installing Admixem

Follow the instuctions on https://github.com/melop/admixem for installing Admixem

#Configure input files

Follow the instuctions on https://github.com/melop/admixem for configuring the input files

Admixem requires a number of input files. These include:

sexualsel_unlinked_additive_esize0.1_noss_nons.txt
naturalsel_unlinked_additive_esize0.1_noss_nons.txt
admixsimul_unlinked_additive_esize0.1_noss_nons.txt
phenotypes_unlinked_additive_esize0.1_noss_nons.txt
genes_unlinked_additive_esize0.1_noss_nons.txt
samplegens.txt
markers.txt
recombination_hotspots.txt

The samplesgens.txt, markers.txt, and recombination_hotspots.txt files do not need to be editted between each set of simulation parameters. The remaining files with the genetic architecture and selection parameter suffixes do need to editted for each set of simulations. 

#Run loop_admixem_unlinked_additive_esize0.1_noss_nons.bash

This job file will simulate 100 replicates of a given genetic architecture (e.g. unlinked_additive_esize0.1_noss_nons) for 6 different inital allele frequencies and output several files.

This job file requires the following Rscripts:

increment_genes.R
parse_phenotypes3.R
parse_phenotypes4.R
prop_overlap.R

The output is a directory named 'output_unlinked_additive_esize0.1_noss_nons' with the following files:

#summary statistics for each set of initial allele frequencies

df_unlinked_additive_esize0.1_noss_nons.1.csv
df_unlinked_additive_esize0.1_noss_nons.2.csv
df_unlinked_additive_esize0.1_noss_nons.3.csv
df_unlinked_additive_esize0.1_noss_nons.4.csv
df_unlinked_additive_esize0.1_noss_nons.5.csv
df_unlinked_additive_esize0.1_noss_nons.6.csv

#summary statistics for each set of initial allele frequencies combined

df_unlinked_additive_esize0.1_noss_nons_all.csv

#number of replicates in which the cross between parent species 1 and parent species 2 did not produce any offspring 

num_failed_sims_unlinked_additive_esize0.1_noss_nons.1.txt
num_failed_sims_unlinked_additive_esize0.1_noss_nons.2.txt
num_failed_sims_unlinked_additive_esize0.1_noss_nons.3.txt
num_failed_sims_unlinked_additive_esize0.1_noss_nons.4.txt
num_failed_sims_unlinked_additive_esize0.1_noss_nons.5.txt
num_failed_sims_unlinked_additive_esize0.1_noss_nons.6.txt

#proportion of replicates in which the mean hybrid phenotype overlaps with either of the parental species' mean phenotype

proportion_overlap_mean_unlinked_additive_esize0.1_noss_nons.jpeg
prop_overlap_mean_unlinked_additive_esize0.1_noss_nons.1.txt
prop_overlap_mean_unlinked_additive_esize0.1_noss_nons.2.txt
prop_overlap_mean_unlinked_additive_esize0.1_noss_nons.3.txt
prop_overlap_mean_unlinked_additive_esize0.1_noss_nons.4.txt
prop_overlap_mean_unlinked_additive_esize0.1_noss_nons.5.txt
prop_overlap_mean_unlinked_additive_esize0.1_noss_nons.6.txt

#standard output

unlinked_additive_esize0.1_noss_nons-6753949.err
unlinked_additive_esize0.1_noss_nons-6753949.out

#density plots of mean phenotype and coefficient of variation for each set of initial allele frequencies 

variation_mean_and_cv_density_unlinked_additive_esize0.1_noss_nons_p1-0.5_p2-0.5.jpeg
variation_mean_and_cv_density_unlinked_additive_esize0.1_noss_nons_p1-0.6_p2-0.4.jpeg
variation_mean_and_cv_density_unlinked_additive_esize0.1_noss_nons_p1-0.7_p2-0.3.jpeg
variation_mean_and_cv_density_unlinked_additive_esize0.1_noss_nons_p1-0.8_p2-0.2.jpeg
variation_mean_and_cv_density_unlinked_additive_esize0.1_noss_nons_p1-0.9_p2-0.1.jpeg
variation_mean_and_cv_density_unlinked_additive_esize0.1_noss_nons_p1-1_p2-0.jpeg

#Run parse_output_df.R

After all parameter combinations have been simulated run: Rscript parse_output_df.R
This script will loop through each output directory and combine all the summary statistics files (e.g. df_unlinked_additive_esize0.1_noss_nons_all.csv) into a sinlge .csv.

#Run decision_tree.R
