#### Change path to directory

cd /path/to/Epigenomic_Data/


#### Download Information on TF Binding Sites : ReMAP 2022 Database - CRM used

filename : remap2022_crm_macs2_hg38_v1_0.bed


#### Prepare Seperate bed files for each TF

## prepare list of all TF - TF_all.txt

while read F
do
        awk -F',' '{ for (i = 1; i <= NF; i++) if ($i == "'${F}'") { print; break; } }' remap2022_crm_macs2_hg38_v1_0.bed | awk '{print $1,$2,$3}' >TF_bedfiles/${F}.crm.bed #Column IDs: Chr,Start_Pos, End_Pos
done <TF_all.txt


######## Map all 850K CpGs to each TFBS to get background count of all CpGs binding to the specific TFBS

mkdir TFBS_CpG_Overlap

## Divided into 84 sets of 10,000 CpGs to parallelize

split CpG_list_850K.txt CpG_Sets/Probe_set_ -d -a2 -l 10000

## Determine CpGs falling in different TFBS

while read F
do
                                echo "Starting" ${F} 
                                awk '{print $1,$2,$3}' TF_bedfiles/${F}.crm.bed >Match_file
                                for r in {00..84}
                                        do
                                                cat CpG_Sets/Probe_set_${r} >probe_set
                                                file1="probe_set"
                                                file2="Match_file"

# Loop over each line of file1
                                                while read -r col1 col2 col3; do
# Match col2 to file2 and select the corresponding rows using col3
                                                        selected_rows=$(awk -v col2="$col2" -v col3="$col3" '$1 == col2 && $2 <= col3 && $3 >= col3 {print $1, $2, $3}' "$file2")

# Output the selected rows with all three columns from file1
                                                        if [[ -n "$selected_rows" ]]; then
                                                                echo "$selected_rows" | while read -r col1_file2 col2_file2 col3_file2; do
                                                                echo "$col1 $col2 $col3 $col1_file2 $col2_file2 $col3_file2"
                                                        done
                                                        fi
                                                done < "$file1" >TF_CpG_Overlap/${F}.all.overlap
                                        done
                                echo "done" ${F}


# Get list of Sentinel CpGs falling in different TFBS

				awk 'NR==FNR{a[$1];next}$1 in a' Sentinel_CpGs.txt TF_CpG_Overlap/${F}.all.overlap >TF_CpG_Overlap/${F}.sentinel.overlap
				rm Match_file
				rm probe_set
done <TF_all.txt



####### get the background and sentinel count for each TF

while read F; 
do
	background_count=$(wc -l TF_CpG_Overlap/${F}.all.overlap)
	sentinel_count=$(wc -l TF_CpG_Overlap/${F}.sentinel.overlap)
	echo ${F}" "${background_count}" "${sentinel_count} >>TF_Count_file.txt
done

########## Post_processing and enrichment analysis in R







