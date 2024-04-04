#### Change path to directory

cd /path/to/Genomic_Data/

#### QC of Raw PLINK Binary files

for i in {1..22};
do 
	/bin/plink2 --bfile HELIOS_All_Chr${i} \
		--remove WGS_QC_failed_samples_others_and_dups.txt \ #list of samples failing WGS QC, Samples with Non-Southeast Asian reported Ancestry and 1 sample out of each duplicate pair [Columns: FID, IID]
		--mind 0.02 \
		--make-bed \
		--out Tmp_file_Chr${i}
	
	/bin/plink2 --bfile Tmp_file_Chr${i} \
		--geno 0.02 \
		--maf 0.05 \
		--hwe 1e-06 \
		--exclude range rangefile.txt \ [Remove MHC and INV Regions]
		--make-bed \
		--out HELIOS_All_QCed_Chr${i}
	
	rm Tmp_file_Chr${i}*

#### Merge QCed files

for i in {2..22}; echo HELIOS_All_QCed_Chr${i} >>mergelist.txt; done

/bin/plink2 --bfile HELIOS_All_QCed_Chr1 \
	--pmerge-list mergelist.txt \
	--make-bed \
	--out HELIOS_All_QCed

#### Pruning to select Independent Variants

/bin/plink2--bfile HELIOS_All_QCed \
	--indep-pairwise 200 100 0.1 \
	--out Prunelist
	
/bin/plink2 --bfile HELIOS_All_Qced \
	--extract Prunelist.prune.in \
	--make-bed \
	--out HELIOS_All_QCed_pruned


	
