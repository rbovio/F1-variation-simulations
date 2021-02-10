#!/bin/sh
#SBATCH -J unlinked_additive_esize0.1_noss_nons
#SBATCH --time=06:00:00
#SBATCH -n1
#SBATCH --mem=32000
#SBATCH --mail-user=rbovio@tamu.edu
#SBATCH -o unlinked_additive_esize0.1_noss_nons-%j.out
#SBATCH -e unlinked_additive_esize0.1_noss_nons-%j.err



#######################
#EDIT THESE PARAMETERS#
#######################
parameters=unlinked_additive_esize0.1_noss_nons;
genes_file=genes_$parameters.txt
admixem_file=admixsimul_$parameters.txt

#make and empty directory so admixem starts making directories for each replicate that have the prefix number
mkdir F1_$parameters

export LANG=en_US.UTF-8  
#load R module
module load R_tamu/3.6.2-intel-2019b-recommended-mt

#iterate through initial allele frequencies
#change k for the number of initial allele frequencies you want to simulate
for((k=1; k<=6; k++)){
	#value to increment intital allele frequencies by
	increment_value=0.10
	Rscript increment_genes.R $genes_file $increment_value

	#Run job_file.slrm 1000 times
	echo "starting at date on hostname"
	#change i for number of replicates you want to simulate
	for ((i=1;i<=100;i++))
	do
	echo "replicate $i"

	#perl -i.bak -lpe 'BEGIN { sub inc { my ($num) = @_; ++$num } } s/(RandomSeed=)(\d+)/$1 . (inc($2))/eg' $admixem_file
	#rm $admixem_file.bak

	./admixemp $admixem_file

	rm /scratch/user/rbovio/admixem3.0/admixem-master3.0/bin/F1_$parameters\_$i/Gen0_genes.txt
	rm /scratch/user/rbovio/admixem3.0/admixem-master3.0/bin/F1_$parameters\_$i/Gen0_markers.txt
	rm /scratch/user/rbovio/admixem3.0/admixem-master3.0/bin/F1_$parameters\_$i/Gen0_natselprobdump.txt
	rm /scratch/user/rbovio/admixem3.0/admixem-master3.0/bin/F1_$parameters\_$i/Gen0_phenostats.txt
	rm /scratch/user/rbovio/admixem3.0/admixem-master3.0/bin/F1_$parameters\_$i/Gen1_genes.txt
	rm /scratch/user/rbovio/admixem3.0/admixem-master3.0/bin/F1_$parameters\_$i/Gen1_markers.txt
	rm /scratch/user/rbovio/admixem3.0/admixem-master3.0/bin/F1_$parameters\_$i/Gen1_natselprobdump.txt
	rm /scratch/user/rbovio/admixem3.0/admixem-master3.0/bin/F1_$parameters\_$i/Gen1_phenostats.txt
	
	sleep 1s
	#rm slurm*
	done
	echo "ended at date on hostname"
	#exit 0

	#iterate through each output directory
	for dir in F1_$parameters\_*/; do
		(cd $dir
		#split phenotypes file
		csplit -z -f parents Gen0_phenotypes.txt '/id/' '{*}'
		csplit -z -f hybrids Gen1_phenotypes.txt '/id/' '{*}'
		rm Gen0_phenotypes.txt
		rm Gen1_phenotypes.txt
		cd ..)
	done
	
	#set command line variables
	af1=$(awk 'NR == 4 {print $9}' $genes_file); echo $af1
	af2=$(awk 'NR == 4 {print $10}' $genes_file); echo $af2
	
	#Run phenotype_parser - plots mean/sd/cv and calculates proportion of overlap between hybrid and parent populations among all replicates.
	Rscript phenotype_parser3.R $parameters $af1 $af2
	mv prop_overlap_mean_$parameters.txt prop_overlap_mean_$parameters.$k.txt
	mv num_failed_sims_$parameters.txt num_failed_sims_$parameters.$k.txt 
	mv df_$parameters.csv df_$parameters.$k.csv
	rm -r F1_$parameters\_*

}

#Rscript to plot prop_overlap x initial allele frequencies
Rscript prop_overlap.R $parameters

#combine dataframe output for each initial allele freq into a single spreadsheet
head -1 df_$parameters.1.csv > df_$parameters\_all.csv; tail -n +2 -q df_$parameters.* >> df_$parameters\_all.csv

mkdir output_$parameters
mv variation_mean_and_cv_density_$parameters* output_$parameters
mv prop_overlap_mean_$parameters* output_$parameters
mv proportion_overlap_mean_$parameters* output_$parameters
mv num_failed_sims_$parameters* output_$parameters
mv df_$parameters* output_$parameters
mv unlinked_additive_esize0.1_noss_nons-* output_$parameters
