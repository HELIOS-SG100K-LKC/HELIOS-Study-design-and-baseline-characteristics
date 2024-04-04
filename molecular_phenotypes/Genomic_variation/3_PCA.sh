#### Change path to directory

cd /path/to/Genomic_Data/

#### PCA of QCed and Pruned File

/bin/plink --bfile HELIOS_All_QCed_pruned \
	--pca 20 header \
	--out HELIOS_All

## Output files : HELIOS_All.eigenvec and HELIOS_All.eigenval

#### PCA with 1000 Genomes Phase 3 data

awk 'NR==FNR{a[$2];next}$2 in a' HELIOS_All_QCed_pruned.bim 1kg_GRCh38.bim >common_variant_list.txt

/bin/plink --bfile HELIOS_All_QCed_pruned \
        --extract common_variant_list.txt \
        --make-bed \
        --out HELIOS_filtered 

/bin/plink --bfile 1kg_GRCh38_All \
	--keep Sample_list.txt \     #Only EUR, AFR, SAS and EAS individuals
	--extract common_variant_list.txt \
	--make-bed \
	--out 1KG_filtered

/bin/plink --bfile HELIOS_filtered \
	--bmerge 1KG_filtered \
	--make-bed \
	--out HELIOS_1KG

/bin/plink --bfile HELIOS_1KG \
	--pca 20 header \
	--out HELIOS_1KG

## Output files : HELIOS_1KG.eigenvec and HELIOS_1KG.eigenval
