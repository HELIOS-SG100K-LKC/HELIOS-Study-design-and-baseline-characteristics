### This script contains the code used to create the circosplot ###

# install.packages("circlize")
library(circlize)

circos.initializeWithIdeogram(chromosome.index = paste0("chr", c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)))


### INNERMOST RING - T2D ###
# read file containing the 314 T2D sentinels - file consists of 5 columns (cpg, chr, start, end, minuslog10P)
# minuslog10P contains the -log10(P) for the p-values generated from association testing with T2D as outcome in TOAST study
# we set the 'start' and 'end' position for each CpG as a +/-2Mb region around the CpG to ensure bars displayig minuslog10P will be thick enough to be visible
# in order to not compress the vertical axes too much, minuslog10P >22 are all set to 22 
T2D_sentinel <- read.delim("<insert filename>", header = TRUE, stringsAsFactors = FALSE)

#set track height for T2D ring
circos.par("track.height" = 0.1)

#plot innermost ring - light yellow for T2D association results
circos.genomicTrack(T2D_sentinel[,2:5], panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, ytop.column = 1, col="cyan", ...)
},bg.col = "lightyellow")


### 2ND INNERMOST RING - BMI ###
# read file containing BMI association results for the 314 T2D sentinels - file consists of 5 columns (cpg, chr, start, end, minuslog10P)
# minuslog10P contains the -log10(P) for the p-values generated from association testing with BMI as outcome in HELIOS study
# as above, we set the 'start' and 'end' position for each CpG as a +/-2Mb region around the CpG to ensure bars displayig minuslog10P will be thick enough to be visible
T2D_sentinel_BMI <- read.delim("<insert filename>", header=T, stringsAsFactors = FALSE)

#set track height for BMI ring
circos.par("track.height" = 0.1)

#plot 2nd innermost ring - light blue for BMI association results 
circos.genomicTrack(T2D_sentinel_BMI[,2:5], panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, ytop.column = 1, col="pink", ...)
},bg.col = "lightblue")


### 3RD INNERMOST RING - METHYLATION FOR T2D LOCI BY ETHNICITY ###
# read file stating which ethnicity has the highest risk for T2D according to average methylation values for the ethnic group - file consists of 5 columns (cpg, chr, start, end, ethnicity)
# ethnicity column indicates 1, 2 or 3 (1: Chinese, 2: Malay, 3: Indian)
ethnic_worst_T2D <- read.delim("<insert filename>", header=T, stringsAsFactors = FALSE)

#set track height for T2D risk by ethnicity ring
circos.par("track.height" = 0.1)

#plot 3rd innermost ring - colored by which ethnicity has the highest risk for T2D 
#heatmap color - Chinese, Malay, Indian
col_fun = colorRamp2(breaks = c(1, 2, 3), colors = c("green", "yellow", "red"))
circos.genomicTrack(ethnic_worst_T2D[,2:5],  panel.fun = function(region, value, ...) {
        circos.genomicRect(region, value, col = col_fun(value[1]), border = NA, ...)
})


### CORRELATION BETWEEN PAIRS OF CPG ###
# read first file 
# file consists of pairs of CpGs with absolute (Pearson) correlation between CpG sites that lies between 0.5 and 0.6 - file consists of 4 columns (chr of CpG1, bp of CpG1, chr of CpG2, bp of CpG2)
cor1 <- read.delim("<insert filename>", header=T, stringsAsFactors = FALSE)

# Each loop plots 1 GRAY link between a pair of CpG sites with absolute correlation between 0.5 and 0.6
for (i in 1:nrow(cor1)){
	circos.link(cor1[i,1],cor1[i,2],cor1[i,3],cor1[i,4],col="gray")
	}

# read second file
# file consists of pairs of CpGs with absolute (Pearson) correlation between CpG sites greater than 0.6 - file consists of 4 columns (chr of CpG1, bp of CpG1, chr of CpG2, bp of CpG2)
cor2 <- read.delim("<insert filename>", header=T, stringsAsFactors = FALSE)
# blue links indicate absolute (Pearson) correlation between CpG sites greater than 0.6

# Each loop plots 1 BLUE link between a pair of CpG sites with absolute correlation greater than 0.6
for (i in 1:nrow(cor2)){
	circos.link(cor2[i,1],cor2[i,2],cor2[i,3],cor2[i,4],col="blue")
	}
    
circos.clear()