#### Change path to directory

cd /path/to/Genomic_Data/

#### Seperate files based on Ancestry and perform GWAS pre-processing QC steps


for i in {1..22};
do
	for Anc in Chinese Indian Malay
	do
		/bin/plink2 --bfile HELIOS_All_Chr${i} \
			--keep ${Anc}.ids.txt \
			--make-bed \
			--out HELIOS_${Anc}_Chr{i}
	done
done


#### Merge all chromosomes into 1 file

for Anc in Chinese Indian Malay
do
	/bin/plink2 --bfile HELIOS_${Anc}_Chr1 \
		--pmerge-list ${Anc}_mergelist.txt \
		--make-bed \
		--out HELIOS_${Anc}
done



#### Quality Control

for Anc in Chinese Indian Malay
do
# Variant QC
        /bin/plink --bfile HELIOS_${Anc} \
                --geno 0.02 \
                --mind 0.02 \
                --maf 0.005 \
		--hwe 1e-06 \
		--make-bed \
		--out HELIOS_${Anc}_QC1		

# Individual QC
	/bin/plink --bfile HELIOS_${Anc}_QC1 \
		--maf 0.05 \
		--exclude range rangefile.txt #MHC and INV \
		--indep-pairwise 200 100 0.1 \
		--out Indep_SNP

	/bin/plink --bfile HELIOS_${Anc}_QC1 \
		--extract Indep_SNP.prune.in \
		--het \
		--genome \
		--out HELIOS_${Anc}_reports

### extract list of samples with IBC >|0.2| and pi_hat >0.75 [${Anc}.exclude.list.txt]

	/bin/plink --bfile HELIOS_4{Anc}_QC1 \
		--remove ${Anc}.exclude.list.txt \
		--maf 0.005 \
		--hwe 1e-06 \
		--geno 0.02 \
		--make-bed \
		--out HELIOS_${Anc}_QCrem_miss02_maf005_hwe1e06


		rm HELIOS_${Anc}_QC1*
done



####### Generate top 20 PCs for each clean dataset ( to be used as Covariates in GWAS)

for Anc in Chinese Indian Malay
do
	/bin/plink --bfile HELIOS_${Anc}_QCrem_miss02_maf005_hwe1e06 \
                --maf 0.05 \
                --exclude range rangefile.txt #remove MHC and INV \
                --indep-pairwise 200 100 0.1 \ 
                --out Indep_SNP
                        
        /bin/plink --bfile HELIOS_${Anc}_QCrem_miss02_maf005_hwe1e06 \
                --extract Indep_SNP.prune.in \
                --pca 20 header
                --out HELIOS_${Anc}_QCrem_miss02_maf005_hwe1e06    #output_files: HELIOS_${Anc}_QCrem_miss02_maf005_hwe1e06.eigenvec
done











