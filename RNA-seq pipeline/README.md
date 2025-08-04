# mRNA-Seq Pipeline
 This nextflow pipeline is used for processing mRNA sequencing fastq data and has been used in the TOAST manuscript [ref]. 

# Processing steps
 This nextflow pipeline contains the following key processing steps: 
 1) Quality check using [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
 2) Trim adaptors and low quality bases using [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/). 
 3) Quality check (using [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) of trimmed reads from step 2.
 4) Removal of ribosomal RNA (rRNA) using [sortMeRNA](https://github.com/sortmerna/sortmerna).
 5) Quality check (using [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) of RNA-seq reads (with rRNA removed) from step 4. 
 6) Map QC-ed reads to GRCh38 reference genome using [STAR algorithm](https://github.com/alexdobin/STAR/releases/tag/2.7.9a).
 7) Quantification of transcript abundances using [RSEM](https://github.com/deweylab/RSEM).
 8) Collect RNA-seq metrics using [Picard algorithm](https://broadinstitute.github.io/picard/).
 9) Determine the library insert size using [Picard algorithm](https://broadinstitute.github.io/picard/).
 10) Determine coverage using [bedtools](https://bedtools.readthedocs.io/en/latest/).

<br/>
