## Study: Rare Variant Analyses in Ancestrally Diverse Cohorts Reveal Novel ADHD Risk Genes
## Analysis: heterogeneity analyses using Fisher’s exact test between ADHD and ASD
## Author: Seulgi Jung

rm(list=ls())
library(dplyr)

## set work directory
setwd("C:\\Users\\jsg79\\Dropbox\\MSSM\\Seulgi\\ADHD\\Het_test") 


adhd <- read.table("SPARK_Olfson_de_novo_variant_merge.txt", header=T)
asd <- read.table("Fu_de_novo_probands.txt", header=T)
genes <- read.table("15_all_genes_with_q_0.1.txt", header=T)
total_genes <- read.table("18128genes_in_Fu.txt", header=T)

head(adhd)
head(asd)

adhd <- adhd[which(adhd$Gene %in% total_genes$gene),]
asd <- asd[which(asd$gene %in% total_genes$gene),]

## Keep protein-truncating variants (PTVs) and missense variants with MPC score >= 1 within 15 ADHD risk genes
adhd_remain1 <- adhd[which(adhd$isPTV == "TRUE"),]
adhd_remain2 <- adhd[which((adhd$isMis == "TRUE") & (adhd$MPC >= 1)),]
adhd_remain_pre <- rbind(adhd_remain1, adhd_remain2)
adhd_remain <- adhd_remain_pre[-which(adhd_remain_pre$Gene %in% genes$Gene),]
nrow(adhd_remain) ## 240


asd_remain1 <- asd[which(asd$isPTV == "TRUE"),]
asd_remain2 <- asd[which((asd$isMis == "TRUE") & (asd$MPC >= 1)),]
asd_remain_pre <- rbind(asd_remain1, asd_remain2)
asd_remain <- asd_remain_pre[-which(asd_remain_pre$gene %in% genes$Gene),]
nrow(asd_remain) ## 5050


## Make input file
adhd_ptv <- adhd[which(adhd$isPTV == "TRUE"),]
adhd_mis <- adhd[which((adhd$isMis == "TRUE") & (adhd$MPC >= 1)),]
adhd_input_pre <- rbind(adhd_ptv, adhd_mis)
adhd_input <- adhd_input_pre[which(adhd_input_pre$Gene %in% genes$Gene),]

# Count the number of occurrences of each gene_id and make a table
adhd_count <- adhd_input %>%
  group_by(Gene) %>%
  summarise(adhd_count = n()) %>%
  as.data.frame()


# View the table
adhd_count


asd_ptv <- asd[which(asd$isPTV == "TRUE"),]
asd_mis <- asd[which((asd$isMis == "TRUE") & (asd$MPC >= 1)),]
asd_input_pre <- rbind(asd_ptv, asd_mis)
asd_input <- asd_input_pre[which(asd_input_pre$gene %in% genes$Gene),]

# Count the number of occurrences of each gene_id and make a table
asd_count <- asd_input %>%
  group_by(gene) %>%
  summarise(asd_count = n()) %>%
  as.data.frame()

names(asd_count)[1] <- c("Gene")
# View the table
asd_count

#########################

## merge

m <- merge(adhd_count, asd_count, by="Gene", all=TRUE)
head(m)

m[is.na(m)] <- 0

m$adhd_remain <- nrow(adhd_remain)
m$asd_remain <- nrow(asd_remain)

m

#############################################

## Fisher

# Create a new column for Fisher's exact p-value
m$fisher_p_value <- apply(m, 1, function(row) {
  # Create a 2x2 contingency table
  contingency_table <- matrix(c(as.numeric(row["adhd_count"]), 
                                as.numeric(row["adhd_remain"]),
                                as.numeric(row["asd_count"]),
                                as.numeric(row["asd_remain"])),
                              nrow = 2)
  # Calculate Fisher's exact test
  fisher_test <- fisher.test(contingency_table)
  # Return the p-value
  return(fisher_test$p.value)
})


m$OR <- apply(m, 1, function(row) {
  # Create a 2x2 contingency table
  contingency_table <- matrix(c(as.numeric(row["adhd_count"]), 
                                as.numeric(row["adhd_remain"]),
                                as.numeric(row["asd_count"]),
                                as.numeric(row["asd_remain"])),
                              nrow = 2)
  # Calculate Fisher's exact test
  fisher_test <- fisher.test(contingency_table)
  # Return the p-value
  return(fisher_test$estimate)
})



# View the updated table with Fisher's p-values
m
m[which(m$fisher_p_value < (0.05/15)),]

m$fdr_BH <- p.adjust(m$fisher_p_value, method = "fdr")

write.table(m, "Run_Fisher_test_ADHD_ASD_total_15_denovo_genes_with_q_0.1_result.txt", col.names=T, row.names=F, quote=F, sep="\t")


################################################################################
################################################################################
## Case-control

adhd <- read.table("C:\\Users\\jsg79\\Dropbox\\MSSM\\Seulgi\\ADHD\\TADA\\Fu\\tada_template_spark_902cases_3508controls_olfson_allofus_kyle_ID_alphamis_input_current_final.txt", header=T)
asd <- read.table("Fu_case-control.txt", header=T)
genes <- read.table("15_all_genes_with_q_0.1.txt", header=T)
total_genes <- read.table("18128genes_in_Fu.txt", header=T)

head(adhd)
head(asd)

adhd <- adhd[which(adhd$gene %in% total_genes$gene),]
asd <- asd[which(asd$gene %in% total_genes$gene),]


## remaining variants
adhd_remain <- adhd[-which(adhd$gene %in% genes$Gene),]
asd_remain <- asd[-which(asd$gene %in% genes$Gene),]

adhd_remain_count <- sum(adhd_remain$cc.ptv.case, adhd_remain$cc.misb.case)
asd_remain_count <- sum(asd_remain$case.ptv, asd_remain$case.misb)


## deleterious variants
adhd_risk <- adhd[which(adhd$gene %in% genes$Gene),]
asd_risk <- asd[which(asd$gene %in% genes$Gene),]

adhd_risk$adhd_count <- adhd_risk$cc.ptv.case + adhd_risk$cc.misb.case
asd_risk$asd_count <- asd_risk$case.ptv + asd_risk$case.misb


## merge

m <- merge(adhd_risk, asd_risk, by="gene", all=TRUE)
head(m)

m$adhd_remain <- adhd_remain_count
m$asd_remain <- asd_remain_count

m[is.na(m)] <- 0
m
#############################################

## Fisher

# Create a new column for Fisher's exact p-value
m$fisher_p_value <- apply(m, 1, function(row) {
  # Create a 2x2 contingency table
  contingency_table <- matrix(c(as.numeric(row["adhd_count"]), 
                                as.numeric(row["adhd_remain"]),
                                as.numeric(row["asd_count"]),
                                as.numeric(row["asd_remain"])),
                              nrow = 2)
  # Calculate Fisher's exact test
  fisher_test <- fisher.test(contingency_table)
  # Return the p-value
  return(fisher_test$p.value)
})




# View the updated table with Fisher's p-values
m
m[which(m$fisher_p_value < (0.05/15)),]

cc_result <- select(m, c(gene, adhd_count,	asd_count,	adhd_remain,	asd_remain,	fisher_p_value))

cc_result$fdr_BH <- p.adjust(cc_result$fisher_p_value, method = "fdr")

write.table(cc_result, "Run_Fisher_test_ADHD_ASD_total_result_case-control.txt", col.names=T, row.names=F, quote=F, sep="\t")


