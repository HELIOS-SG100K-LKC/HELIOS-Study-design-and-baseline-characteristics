#### Change path to directory

cd /path/to/Genomic_Data/


#### Sumstats downloaded from PGS Catalog [harmonised] 
#### Code shown for 1 disorder [same followed for all traits]


awk 'NR>20{print $9":"$10,$4,$5,$6}' PGS_T2D.txt >T2D.sumstats.txt #Column IDs: Chr:Pos(ID),A1,A2,Effect_size

#### Seperate HELIOS data into Ancestry groups based on predicted ancestry

awk '{if($4=="Chinese") print $1,$2}' Final_Clustering_Information.txt >Chinese.ids.txt
awk '{if($4=="Indian") print $1,$2}' Final_Clustering_Information.txt >Indian.ids.txt
awk '{if($4=="Malay") print $1,$2}' Final_Clustering_Information.txt >Malay.ids.txt

###### Split files

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

####### Create merge lists

for i in {2..22}; do
	for Anc in Chinese Indian Malay; do
	echo HELIOS_${j}_Chr{i} >>${Anc}_mergelist.txt
done
done

####### Merge all chromosomes into 1 file

for Anc in Chinese Indian Malay
do
	/bin/plink2 --bfile HELIOS_${Anc}_Chr1 \
		--pmerge-list ${Anc}_mergelist.txt \
		--make-bed \
		--out HELIOS_${Anc}
done



##### PRS esimtation [Make sure the SNP IDs are in the same format in the Sumstats and the PLINK bim files]

for Anc in Chinese Indian Malay; 
do
	/bin/plink --bfile HELIOS_${Anc} \
		--score T2D.sumstats.txt 1 2 4 header \
		--out ${Anc}_T2D
done



## output : HELIOS_"Ancestry"_T2D.profile


