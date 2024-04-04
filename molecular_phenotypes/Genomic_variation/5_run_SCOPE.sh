#### Change path to directory

cd /path/to/Genomic_Data/

#### Run SCOPE using pruned data with supervised information about allele frequency of variants within each of the 3 homogenous clusters

/bin/plink --bfile HELIOS_All_QCed_pruned \
	--freq \
	--within Homogenous_ancestry_info.txt \
	--out HELIOS_All

## Output file : HELIOS_All.frq.strat

awk 'NR>1{if($3!=4) print $0}' HELIOS_All_frq.strat >HELIOS_All_CLST3.frq.strat 


## run Supervised SCOPE

/bin/scope -g HELIOS_All_QCed_pruned \
	-freq HELIOS_All_CLST3.frq.strat \
	-output HELIOS_All_Supervised_PC1_5_

## Output files : HELIOS_All_Supervised_PC1_5_V.txt; HELIOS_All_Supervised_PC1_5_Qhat.txt (Ancestry proportions); HELIOS_All_Supervised_PC1_5_Phat.txt (Allele frequencies per population)

  
