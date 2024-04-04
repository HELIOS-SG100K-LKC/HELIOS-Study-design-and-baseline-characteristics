cd /Path/to/Genomic_data/Results

mkdir smr_analysis_sumstats
mkdir SMR_reults


### Replace PHENO with metabolite of interest
### Downloaded eqtlgen cis-eQTL data


echo "SNP A1 A2 freq b se p N" >smr_analysis_sumstats/PHENO.smr.sumstats
awk -F "\t" 'NR>1{if($9<=1e-05) print "chr"$1":"$2,$4,$5,$6,$7,$8,$9,$10}' Final_Sumstats_filtered/HELIOS_All_PHENO_filt2studies_hetp05_AFdiff05.summary.txt >>smr_analysis_sumstats/PHENO.smr.sumstats
grep -w PHENO ../Metab_all_signif_5e08.txt | awk '{print "chr"$2":"$3}' > PHENO.target.snps ## SNPs with p <5e-08
~/bin/smr-1.3.1 \
	--bfile HELIOS_Chinese_QCrem_miss02_maf005_hwe1e06 \
	--beqtl-summary ~/eQTLgen_Data/grch38_cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212_besd-dense \
	--gwas-summary smr_analysis_sumstats/PHENO.smr.sumstats \
	--peqtl-cis 5e-08 \
	--extract-snp PHENO.target.snps  \
	--heidi-off \
	--diff-freq 1 \
	--out SMR_Results/PHENO