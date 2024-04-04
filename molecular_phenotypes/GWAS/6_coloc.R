setwd("Path/to/Genomic_Data/Results")

##### Prepare regional summary statistics for each metabolite Gene pair based on SMR results (1MB on each side of top snp)
## Columns : SNP, A1, A2, Freq, Beta, SE, P, N


library(coloc)
## read regional sumstats
dat1 <- read.table("Metab_sumstats/METAB.regional.sumstats",header=T)
dat2 <- read.table("Gene_sumstats/GENE_ID.regional.sumstats",header=T)
## harmonize sumstats
Common_table <- merge.data.frame(dat1,dat2,by="SNP")
Common_table$b_gene <- ifelse(Common_table$A1==Common_table$A1_gene & Common_table$A2==Common_table$A2_gene,Common_table$b_gene,
				ifelse(Common_table$A1==Common_table$A2_gene & Common_table$A2==Common_table$A1_gene,-1*Common_table$b_gene,NA))
Common_table <- na.omit(Common_table)
Common_table$varbeta <- (Common_table$se)^2
Common_table$varbeta_gene <- (Common_table$se_gene)^2

## Prepare coloc files
D1 <- as.list(Common_table[,c(1,2,5,6,17,9)])
D1$type <- "quant"
names(D1) <- c("snp","position","MAF","beta","varbeta","N","type")
D2 <- as.list(Common_table[,c(1,2,12,13,18,16)])
D2$type <- "quant"
names(D2) <- c("snp","position","MAF","beta","varbeta","N","type")

## Perform coloc
my_res <- coloc.abf(dataset1=D1,dataset2=D2)

## save P-value results
output <- cbind(my_res$summary[[1]],my_res$summary[[2]],my_res$summary[[3]],my_res$summary[[4]],my_res$summary[[5]],my_res$summary[[6]])

## save SNP_Count,H0,H1,H2,H3,H4 posterior p-values

write.table(output,file="Coloc_results/coloc_METAB_GENE_ID.results.txt",row.names=F,col.names=F,quote=F,sep=" ")