#### Change path to directory

cd /path/to/Genomic_Data/


## Select SNP Subset for regenie Step 1

for Anc in Chinese Indian Malay;
	do
		/bin/plink2 --bfile HELIOS_${Anc}_QCrem_miss02_maf005_hwe1e06 \
				--maf 0.05 \
				--exclude range rangefile \
				--indep-pairwise 200 100 0.5 \
				--out HELIOS_${Anc}
	done


### Regenie Step 1

for ${Anc} in Chinese Indian Malay
	do
		/bin/regenie --step 1 \
			--bed HELIOS_${Anc}_QCrem_miss02_maf005_hwe1e06 \
			--extract HELIOS_${Anc}.prune.in \
			--covarFile ${Anc}_metab_covar.txt \
			--catCovarList Gender \
			--phenoFile Phenotypes/${Anc}_Phenotype_File_N.txt \ ## replace N with number of the Phenotype file
			--loocv \
			--bsize 1000 \
			--force-step1 \
			--lowmem \
			--out Results/HELIOS_${Anc}
	done



### Regenie Step 2

for ${Anc} in Chinese Indian Malay
	do
		/bin/regenie --step 2 \
			--bed HELIOS_${Anc}_QCrem_miss02_maf005_hwe1e06 \
			--covarFile ${Anc}_metab_covar.txt \
			--catCovarList Gender \
			--phenoFile Phenotypes/${Anc}_Phenotype_File_N.txt \ ## replace N with number of the Phenotype file
			--bsize 200 \
			--minMAC 0.5 \
			--pred Results/HELIOS_${Anc}_pred.list \
			--out Results/HELIOS_${Anc}_QCrem_miss02_maf005_hwe1e06
	done


### Output files : *regenie summary statistics for each ancestry and each phenotype




######################################################################################

## Run Meta-analysis

while read F
do
	/bin/metal metal_inputfile_${F}.txt
done <Pheno_list.txt
## Clean metal results to select final list of variants (filtering for present in 2 studies and hetereogeneity)


while read F
	do
		awk '{if($14>=1) print $0}' Results/HELIOS_All_QCrem_miss02_maf005_hwe1e06_${F}_1.sumstats.txt | \ ## Present in 2 studies (Het_df>1)
		awk '{if($15>0.05) print $0}' | \ ## Het_p > 0.05
		awk '{if($7-$6 < 0.5) print $0}' \ ## AF difference between any 2 populations  < 0.5
		>Results/Final_Sumstats_filtered/HELIOS_All_${F}_filt2studies_hetp05_AFdiff05.summary.txt 
	done <Pheno_list.txt


###################################


## Extract Significant Associations for all Metabolites


while read F
	do
		awk '{if(12 <=5e-08) print $0}' Results/Final_Sumstats_filtered/HELIOS_All_${F}_filt2studies_hetp05_AFdiff05.summary.txt \
		>Results/Final_Sumstats_filtered/${F}.significant.associations ## Used for SMR and Plotting

	done <Pheno_list.txt


### Clumping to identify top independent variants


while read F
	do
		/bin/plink 
				--bfile HELIOS_Chinese_QCrem_miss02_maf005_hwe1e06 \ ## Using largest population as reference
				--clump Results/Final_Sumstats_filtered/HELIOS_All_${F}_filt2studies_hetp05_AFdiff05.summary.txt \
        		--clump-p1 5e-08 \
        		--clump-r2 0.01 \
        		--clump-kb 500 \
        		--clump-range ~/hg38_genelist \ ##To find genes close top SNP
        		--clump-field P-value \
        		--clump-snp-field MarkerName \
        		--out Methods_Paper_Analysis/Signif_Assoc/${F}
    done <Pheno_list.txt




