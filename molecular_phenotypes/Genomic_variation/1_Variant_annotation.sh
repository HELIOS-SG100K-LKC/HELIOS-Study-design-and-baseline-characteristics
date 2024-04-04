#### Change path to directory with plink binaries

cd /path/to/Genomic_data/


#### Prepare combined file of all autosomal variants [239.6 million]

for i in {1..22}
do
	awk '{print $1,$4,$4,$6,$5}' HELIOS_All_Chr${i}.bim >>HELIOS_all_variants.bed  #Column IDs : Chr, BP_start, BP_end, Ref_Allele, Alt_Allele [BP start and end were same values]
done

## Make sure that for INDEL the reference allele is always the single Allele and the Alt is the INDEL


#### Annotate variants using Annovar

/bin/annovar/annotate_variation.pl -geneanno \
				-dbtype refGene \
				-buildver hg38 \
				HELIOS_all_variants.bed \
				/bin/annovar/humandb/ \
				-out HELIOS_all_variants


## Output files :  HELIOS_all_variants.variant_function & HELIOS_all_variants.exonic_variant_function


#### Get count of different variant_function

awk '{print $1}' HELIOS_all_variants.variant_function | sort | uniq -c >HELIOS_variant_function_counts.txt
awk '{print $1}' HELIOS_all_variants.exonic_variant_function | sort | uniq -c >HELIOS_exonic_variant_function_counts.txt



#### Determine number of unreported variants in our dataset

wget https://ftp.ncbi.nlm.nih.gov/snp/archive/b156/VCF/GCF_000001405.40.gz   #download dbSNP data [requires post-processing to update Chromosome IDs present in first 37 lines]
gunzip GCF_000001405.40.gz
awk 'NR>37' GCF_000001405.40 | awk '{print $1,$2,$3,$4,$5}' >dbsnp_b156_snps.txt  #Columns IDs: Chr, Pos, rsID, Ref_Allele, Alt_Allele


## prep files for matching
for i in {1..22}
	do
		awk '{if($1=="'${i}'") print $3,"chr"$1":"$2":"$4":"$5}' >dbSNP_Chr_${i}.snps.txt #Column IDs: rsID, Matching_ID
	done

## prep HELIOS_variant_files

for i in {1..22};
	do
		awk '{print $2,"chr"$1":"$4":"$5":"$6}' HELIOS_All_Chr${i}.bim >HELIOS_chr${i}_all_variants.txt ## Column IDs: SNP_ID, Matching_ID
		awk '{print $2,"chr"$1":"$4":"$6":"$5}' HELIOS_All_Chr${i}.bim >>HELIOS_chr${i}_all_variants.txt ## both orientations in case the alleles are switched in dbSNP database

	done

## compare HELIOS variants to dbSNP

for i in {1..22};
	do
		awk 'NR==FNR{a[$2];next} !($2 in a)' dbSNP_Chr_${i}.snps.txt HELIOS_chr${i}_all_variants.txt | \
		awk '{print $1}' | \
		sort | \
		uniq -c | \
		if($1==2) print $2}' >HELIOS_chr${i}_novel_variants.txt
 









				
