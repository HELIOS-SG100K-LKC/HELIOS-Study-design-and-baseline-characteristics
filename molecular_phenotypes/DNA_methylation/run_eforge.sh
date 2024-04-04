#### Change path to directory

cd /path/to/Epigenomic_Data/


#### Split 16,444 Sentinel CpGs into 22 files of 750 CpGs

split Sentinel_CpGs.txt CpGs_List_ -d -a2 -l 750

#### Run eFORGE_V2 for each CpG List

mkdir Results

for F in {00..21}
do
	for i in H3K4me1 H3K4me3 H3K9me3 H3K27me3 H3K36me3 DHS chromatin15state-all
	do
		perl /bin/eforge.pl --label CPG_list_${F} --array 850k --f CPGs_List_${F} --dataset erc2-${i} --format probeid --noproxy --reps 1000  --noplot
		mv 0* Results/CPG_list_${F}_${i}
	done
done

#### Post-processing in R



