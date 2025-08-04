#!/usr/bin/env nextflow

// I/O Directories
prjDir="nf_mRNA_seq"
rawFASTQ="RawFASTQ"
resultDir="mRNA_seq_results"
prjPath="rna_seq_pipeline/"

userDir="rna_seq/"
mRnaSeqDir="nf_rna_seq_pipeline"
toolsDir="rna_seq/mRNA_seq_tools/"

params.outdir="${prjPath}/${userDir}/${mRnaSeqDir}/${prjDir}"

FASTQC_HOME="${prjPath}/${toolsDir}/opt/FastQC"
TRIMGALORE_PATH="${prjPath}/${toolsDir}/opt/TrimGalore-0.6.6"
STAR_PATH="${prjPath}/${toolsDir}/opt/STAR"
SORTMERNA_PATH="${prjPath}/${toolsDir}/opt/sortMeRNA/sortmerna-2.1b"
RSEM_PATH="${prjPath}/${toolsDir}/opt/rsem/bin"
PRESEQ_PATH="${prjPath}/${toolsDir}/opt/preseq/preseq_v2.0"
BEDTOOLS_PATH="${prjPath}/${toolsDir}/opt/bedtools/usr/local/bin/"

CUTADAPT_PATH="~/.local/bin/cutadapt"
PICARD="${prjPath}/${toolsDir}/opt/picard/picard.jar"
GENOME_INDEX="${prjPath}/${toolsDir}/ref_genome/star_index"
GTF_FILE="${prjPath}/${toolsDir}/ref_genome/gtf/gencode.v38.primary_assembly.annotation.gtf"
RSEM_REF="${prjPath}/${toolsDir}/ref_genome/rsem_ref/"
PICARD_RNAREF="${prjPath}/${toolsDir}/ref_genome/picard_rnaREF/refFlat.txt"
GENOME_BED="${prjPath}/${toolsDir}/ref_genome/bedfile/genome_GRCh38_bed_file_final.bed"

// Samples to be analyzed
read_pairs_ch = Channel.fromFilePairs("${baseDir}/$rawFASTQ/*_{1,2}.fq.gz")
read_pairs_ch2 = Channel.fromFilePairs("${baseDir}/$rawFASTQ/*_{1,2}.fq.gz")  


/******* STEP1 - FastQC_raw *********
***************************
*/
process raw_fastqc {

    tag "$sample_id"
    echo true
    
    publishDir "${params.outdir}/${resultDir}/${sample_id}/fastqc_raw", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) from read_pairs_ch

    output:
    tuple val(sample_id),file('*_fastqc.{zip,html}') into ch_fastqc_results_for_multiqc1

    
    script:
    """
    module load python/2.7.12

    ${FASTQC_HOME}/fastqc -t 2 ${reads[0]} ${reads[1]}

    """
}


/******* STEP2 - Trim Galore! *********
***************************
*/
process trim_galore {

    tag "$sample_id"
    echo true
    
    publishDir "${params.outdir}/${resultDir}/${sample_id}/trimmed_reads", pattern: "*trimming_report.txt", mode:'copy'
   
    input:
    tuple val(sample_id), path(reads) from read_pairs_ch2
   
    output:
    tuple val(sample_id), file('*_val_1.fq.gz'), file('*_val_2.fq.gz') into ch_trimmed_reads_for_sortMeRNA, ch_trimmed_reads_fastQC
    tuple val(sample_id), file("*trimming_report.txt") into ch_trim_galore_results_for_multiqc2

    script:

    """             
    module rm python
    module load python/3.5.1

    $TRIMGALORE_PATH/trim_galore --cores 4 --path_to_cutadapt ${CUTADAPT_PATH} --phred33 --paired ${reads[0]} ${reads[1]}
	    
    """
}

/******* STEP3 - FastQC for trimmed reads *********
***************************
*/
process trim_fastqc {

    tag "$sample_id"
    echo true
    
    publishDir "${params.outdir}/${resultDir}/${sample_id}/fastqc_trimmed", mode: 'copy'

    input:
    tuple val(sample_id), file(trimmed_reads_r1), file(trimmed_reads_r2) from ch_trimmed_reads_fastQC

    output:
    tuple val(sample_id), file('*_fastqc.{zip,html}') into ch_trimmed_fastqc_results_for_multiqc3
    
    script:
    """
    module load python/2.7.12

    ${FASTQC_HOME}/fastqc -t 2 ${trimmed_reads_r1} ${trimmed_reads_r2}

    """
}

/******* STEP4 - Detect rRNA *********
***************************
*/
process sortmeRNA {

    tag "$sample_id"
    echo true
    
    publishDir "${params.outdir}/${resultDir}/${sample_id}/sortMeRNA", pattern: "*.log", mode: 'copy'

    input:
    tuple val(sample_id), file(trimmed_reads_r1), file(trimmed_reads_r2) from ch_trimmed_reads_for_sortMeRNA

    output:
    tuple val(sample_id), file('*_rm_rRNA_R1.fastq'), file('*_rm_rRNA_R2.fastq') into ch_non_rRNA_reads_for_alignment, ch_non_rRNA_reads_for_fastqc
    tuple val(sample_id), file('*.log') into ch_star_results_for_multiqc4
    
    

    script:
    """
    
    gunzip -f < ${trimmed_reads_r1} > ${sample_id}_trimmed_R1.fastq
    gunzip -f < ${trimmed_reads_r2} > ${sample_id}_trimmed_R2.fastq

    rm -rf ${params.outdir}/sortMeRNA_tmp/${sample_id}/

    ${SORTMERNA_PATH}/sortmerna --ref /home/projects/12002085/darwin/mRNA_seq_tools/opt/sortMeRNA/sortmerna-2.1b/rRNA_databases/rfam-5.8s-database-id98.fasta \
    --ref /home/projects/12002085/darwin/mRNA_seq_tools/opt/sortMeRNA/sortmerna-2.1b/rRNA_databases/rfam-5s-database-id98.fasta \
    --ref /home/projects/12002085/darwin/mRNA_seq_tools/opt/sortMeRNA/sortmerna-2.1b/rRNA_databases/silva-arc-16s-id95.fasta \
    --ref /home/projects/12002085/darwin/mRNA_seq_tools/opt/sortMeRNA/sortmerna-2.1b/rRNA_databases/silva-arc-23s-id98.fasta \
    --ref /home/projects/12002085/darwin/mRNA_seq_tools/opt/sortMeRNA/sortmerna-2.1b/rRNA_databases/silva-bac-16s-id90.fasta \
    --ref /home/projects/12002085/darwin/mRNA_seq_tools/opt/sortMeRNA/sortmerna-2.1b/rRNA_databases/silva-bac-23s-id98.fasta \
    --ref /home/projects/12002085/darwin/mRNA_seq_tools/opt/sortMeRNA/sortmerna-2.1b/rRNA_databases/silva-euk-18s-id95.fasta \
    --ref /home/projects/12002085/darwin/mRNA_seq_tools/opt/sortMeRNA/sortmerna-2.1b/rRNA_databases/silva-euk-28s-id98.fasta \
    --reads ${sample_id}_trimmed_R1.fastq --reads ${sample_id}_trimmed_R2.fastq \
    --workdir ${params.outdir}/sortMeRNA_tmp/${sample_id}/ \
    --num_alignments 1 \
    --aligned rRNA_reads --other non_rRNA_reads --fastx --threads 10 --paired_in --out2


    mv non_rRNA_reads_fwd.fq ${sample_id}_rm_rRNA_R1.fastq
    mv non_rRNA_reads_rev.fq ${sample_id}_rm_rRNA_R2.fastq

    rm -rf ${params.outdir}/sortMeRNA_tmp/${sample_id}/
    rm ${sample_id}_trimmed_R1.fastq ${sample_id}_trimmed_R2.fastq

    """
}

/******* STEP5 - FastQC after rRNA removal *********
***************************
*/
process rRNA_fastqc {

    tag "$sample_id"
    echo true
    
    publishDir "${params.outdir}/${resultDir}/${sample_id}/fastqc_sortMeRNA", mode: 'copy'

    input:
    tuple val(sample_id), file(non_rRNA_reads_r1), file(non_rRNA_reads_r2) from ch_non_rRNA_reads_for_fastqc

    output:
    tuple val(sample_id), file('*_fastqc.{zip,html}') into ch_trimmed_fastqc_results_for_multiqc5
    
    script:
    """
    module load python/2.7.12

    ${FASTQC_HOME}/fastqc -t 2 ${non_rRNA_reads_r1} ${non_rRNA_reads_r2}

    """
}

/******* STEP6 - mapping of reads to GRCh38 *********
***************************
*/
process starAlign {

    tag "$sample_id"
    echo true
    
    publishDir "${params.outdir}/${resultDir}/${sample_id}/star_align", mode: 'copy'

    input:
    tuple val(sample_id), file(non_rRNA_R1), file(non_rRNA_R2) from ch_non_rRNA_reads_for_alignment

    output:
    tuple val(sample_id), file('*toTranscriptome.out.bam') into ch_star_results_for_rsem
    tuple val(sample_id), file('*Aligned.sortedByCoord.out.bam') into ch_sorted_results_for_picard_rnaMetrics, ch_sorted_results_for_picard_insert_size, ch_sorted_results_for_preseq, ch_sorted_results_for_bedtools
    tuple val(sample_id), file('*') into ch_star_results_for_multiqc6
    

    script:
    """
  
    ${STAR_PATH}/STAR-2.7.9a/source/STAR  \
    --genomeDir ${GENOME_INDEX} \
    --readFilesIn ${non_rRNA_R1} ${non_rRNA_R2} \
    --runThreadN 10 --outFileNamePrefix ${sample_id}_ \
    --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
    --twopassMode Basic --sjdbGTFfile ${GTF_FILE}

    rm ${non_rRNA_R1} ${non_rRNA_R2}

#    --runRNGseed 0 --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 \
#    --outSAMattributes NH HI AS NM MD --quantTranscriptomeBan Singleend \

    """
}

/******* STEP7 - Quantify transcript abundances *********
***************************
*/
process rsem {

    tag "$sample_id"
    echo true
    
    publishDir "${params.outdir}/${resultDir}/${sample_id}/rsem", mode: 'copy'

    input:
    tuple val(sample_id), file(bamTranscriptome) from ch_star_results_for_rsem

    output:
    tuple val(sample_id), file('*') into ch_trimmed_fastqc_results_for_multiqc6

    script:
    """
    
    rm -rf ${params.outdir}/rsem_tmp/${sample_id}/
    mkdir -p ${params.outdir}/rsem_tmp/${sample_id}/
 
    ${RSEM_PATH}/rsem-calculate-expression -p 4 \
    --temporary-folder ${params.outdir}/rsem_tmp/${sample_id}/ \
    --alignments --paired-end ${bamTranscriptome} \
    ${RSEM_REF}/GRCh38 \
    ${sample_id}


    
    """
}

/******* STEP8 - Compute RNA-seq metrics *********
***************************
*/
process rnaMetrics {
    tag "$sample_id"
    echo true
   
    publishDir "${params.outdir}/${resultDir}/${sample_id}/picard_rnaMetrics",  mode:"copy"

    input:
    tuple val(sample_id), path(sorted_bam) from ch_sorted_results_for_picard_rnaMetrics

    output:
    
    tuple val(sample_id),file("*") into ch_trimmed_fastqc_results_for_multiqc7

    script:
    """   

    java -jar $PICARD CollectRnaSeqMetrics \
      I=${sorted_bam} \
      O=${sample_id}_rna_metrics.txt \
      REF_FLAT=${PICARD_RNAREF} \
      STRAND=NONE \
      CHART_OUTPUT=${sample_id}_chart.pdf


    """

}


/******* STEP9 - Determine insert size *********
***************************
*/
process insert_size {
    tag "$sample_id"
    echo true
    
    publishDir "${params.outdir}/${resultDir}/${sample_id}/insert_size_files", mode:"copy"

    input:
    tuple val(sample_id), path(sorted_bam) from ch_sorted_results_for_picard_insert_size

    output:
    tuple val(sample_id),file("*") into ch_trimmed_fastqc_results_for_multiqc8

    script:
    """    

    java -jar $PICARD CollectInsertSizeMetrics \
        I=${sorted_bam} \
        O=${sample_id}_insert_size_metrics.txt \
        H=${sample_id}_insert_size_histogram.pdf \
	
    """

}

/******* STEP10 - Compute coverage *********
***************************
*/
process bedTools {
    tag "$sample_id"
    echo true
   
    publishDir "${params.outdir}/${resultDir}/${sample_id}/bedTools",  mode:"copy"

    input:
    tuple val(sample_id), path(sorted_bam) from ch_sorted_results_for_bedtools

    output:
    
    tuple val(sample_id),file("*") into ch_trimmed_fastqc_results_for_multiqc9

    script:
    """   

    ${BEDTOOLS_PATH}/bedtools genomecov -ibam ${sorted_bam} \
    -bga -split | ${BEDTOOLS_PATH}/bedtools sort > ${sample_id}_locus_coverage.bedGraph
    
    ${BEDTOOLS_PATH}/bedtools coverage -hist -a ${GENOME_BED} \
    -b ${sorted_bam} > ${sample_id}_gene_coverage.txt


    """

}





