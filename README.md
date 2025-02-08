# DR_Analysis
##Differential Expression Analysis of DR-Related Genes in Breast Cancer
counts01A <- read.table("counts01A.txt", header = T, sep = "\t", row.names = 1)
counts01A <- t(counts01A)
counts01A <- as.data.frame(counts01A)
counts01A$Group <- "Disease"
clinical <- read.csv("Clinical.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
clinical$A0_Samples <- paste0(clinical$A0_Samples, ".01A")
clinical$A0_Samples <- gsub("-", ".", clinical$A0_Samples)

er_negative_samples <- clinical[clinical$breast_carcinoma_estrogen_receptor_status == "Negative", "A0_Samples"]
counts01A_er_negative <- counts01A[rownames(counts01A) %in% er_negative_samples, ]

counts11A <- read.table("counts11A.txt", header = T, sep = "\t", row.names = 1)
counts11A <- t(counts11A)
counts11A <- as.data.frame(counts11A)
counts11A$Group <- "Control"
countsA <- rbind(counts01A_er_negative, counts11A)
print(head(countsA))
countsA <- rbind(counts01A, counts11A)
colnames(countsA)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
countsA <- read.table("countsA.txt", header = T, sep = "\t", row.names = 1)
gene_names <- colnames(countsA)[-which(colnames(countsA) == "Group")]
count_data <- as.matrix(countsA[, gene_names])
rownames(count_data) <- countsA$Group
group <- factor(countsA$Group)
y <- DGEList(counts = t(count_data), genes = gene_names, group = group)
y <- calcNormFactors(y)
y <- estimateDisp(y)
fit <- glmQLFit(y)
qlf <- glmQLFTest(fit, coef=2)
res <- topTags(qlf, n = nrow(y$counts))
res_df <- as.data.frame(res) %>%
  rownames_to_column("SYMBOL") %>%
  mutate(logFC = logFC,
         padj = PValue) 
write.csv(res_df, "res_df_er.csv", quote = F, row.names = F)
df <- res_df %>%
  mutate(significant = case_when(logFC > 1 & padj < 0.05 ~ "Up",
                                 abs(logFC) < 1 | padj > 0.05 ~ "None",
                                 logFC < -1 & padj < 0.05 ~ "Down"))

df$significant <- as.factor(df$significant)
write.csv(df, "df_significant_counts_er-.csv", quote = F, row.names = F)
library(dplyr)
rm(countsA)
df_significant <- read.csv("df_significant_counts_er-.csv", header = TRUE)
significant_data <- df_significant %>% filter(significant != "None")
DR_genes <- read.csv("DR-基因-火山图.csv", header = TRUE)
colnames(DR_genes)
DR_gene_symbols <- DR_genes$gene
significant_genes <- significant_data$SYMBOL
common_genes <- intersect(DR_gene_symbols, significant_genes)
df_common_genes <- df_significant %>% filter(SYMBOL %in% common_genes)
write.csv(df_common_genes, "common_genes_counts_er-.csv", row.names = FALSE)
df_common_genes<-read.csv("common_genes_counts_er-.csv")
down_count <- df_common_genes %>% filter(significant == "Down") %>% nrow()
down_count
up_count <- df_common_genes %>% filter(significant == "Up") %>% nrow()
up_count

##Mendelian Randomization Analysis of DR-Related Genes and Breast Cancer Risk
library(TwoSampleMR)
library(ggplot2)
#set up your work directory
setwd("/Users/caroline/Desktop/孟德尔/DR")
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
expo_rt<-extract_instruments(outcome="ieu-a-300",clump = FALSE)
expo_rt<-clump_data(expo_rt,clump_kb = 100,
                    clump_r2 = 0.3,
                    clump_p1 = 5e-08,
                    clump_p2 = 5e-08,
)
PSMB8_snp<-subset(expo_rt,chr.exposure==6 & pos.exposure>32808494-100000 & pos.exposure<32812456+100000)
#load outcome data 
outc_rt <- extract_outcome_data(
  snps = PSMB8_snp$SNP,
  outcomes = 'ieu-a-1126'
)
#harmonise and merge data
harm_rt <- harmonise_data(
  exposure_dat = PSMB8_snp, 
  outcome_dat = outc_rt,action=1)
get_gene_snp <- function(chr, start_pos, end_pos, buffer = 100000) {
  subset(exposure_clumped_data, chr.exposure == chr & pos.exposure > start_pos - buffer & pos.exposure < end_pos + buffer)
}

library(readxl)
library(openxlsx)
file_path <- "gene.xlsx"
gene_data <- read_excel(file_path)
gene_data$SNP_results <- NA
for (i in 1:nrow(gene_data)) {
  chr <- gene_data$CHr[i]
  start_pos <- gene_data$`Location L`[i]
  end_pos <- gene_data$`Location R`[i]
  
  cat("Processing gene:", gene_data$`Positive-related`[i], "\n")
  cat("Chromosome:", chr, "Start Position:", start_pos, "End Position:", end_pos, "\n")
  
  gene_snp <- get_gene_snp(chr, start_pos, end_pos)
  
  cat("Number of SNPs found:", nrow(gene_snp), "\n\n")
 
  gene_data$SNP_results[i] <- nrow(gene_snp)
}


write.xlsx(gene_data, "gene_results.xlsx", rowNames = FALSE)
filtered_gene_data <- gene_data[gene_data$SNP_results != 0, ]

write.xlsx(filtered_gene_data, "gene_results_filtered2.xlsx", rowNames = FALSE)

library(TwoSampleMR)
library(ggplot2)
rm(BMI_data)
rm(clump_dat)
rm(exposure_data)
rm(state)
library(vroom)
breast_gwas_data<-VariantAnnotation::readVcf("ieu-a-1126.vcf.gz")
outcome_data<-gwasvcf_to_TwoSampleMR(vcf = breast_gwas_data)
rm(breast_gwas_data)
colnames(outcome_data)
outcome1_data<-outcome_data[!is.na(outcome_data$SNP),]
rm(outcome1_data)
outcome_data<-as.data.frame(outcome_data)
outcome_breast<-format_data(dat = outcome_data,
                            type = "outcome",
                            snp_col = "SNP",
                            beta_col = "beta.exposure",
                            se_col = "se.exposure",
                            eaf_col = "eaf.exposure",
                            effect_allele_col = "effect_allele.exposure",
                            other_allele_col = "other_allele.exposure",
                            pval_col = "pval.exposure",
                            chr_col = "chr.exposure",
                            pos_col = "pos.exposure")
write.csv(outcome_breast,file = "outcome_breast.csv")
rm(outcome_data)
outcome_breast<-vroom("outcome_breast.csv")
colnames(outcome_breast)
###re
RAB5B_snp<-subset(exposure_clumped_data,chr.exposure==12 & pos.exposure>56367733-100000 & pos.exposure<56390467+100000)
outc_rt <- read_outcome_data(
  snps =  RAB5B_snp$SNP,
  filename = "outcome_breast.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  eaf_col = "eaf.outcome",
  pval_col = "pval_origin.outcome")
harm_rt <- harmonise_data(
  exposure_dat = RAB5B_snp, 
  outcome_dat = outc_rt,action=1)
mr_result<- mr(harm_rt)
View(mr_result)
library(readxl)
genes_data <- read_excel("gene_results_filtered.xlsx", sheet = "Sheet 1")
results_list <- list()
for (i in 1:nrow(genes_data)) {
  gene_name <- genes_data$`Positive-related`[i]
  chr <- genes_data$CHr[i]
  start_pos <- genes_data$`Location L`[i]
  end_pos <- genes_data$`Location R`[i]
 
  gene_snp <- subset(exposure_clumped_data, chr.exposure == chr & pos.exposure > start_pos - 100000 & pos.exposure < end_pos + 100000)

  if (nrow(gene_snp) == 0) {
    next
  }
  outc_rt <- read_outcome_data(
    snps = gene_snp$SNP,
    filename = "outcome_breast.csv",
    sep = ",",
    snp_col = "SNP",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    effect_allele_col = "effect_allele.outcome",
    other_allele_col = "other_allele.outcome",
    eaf_col = "eaf.outcome",
    pval_col = "pval_origin.outcome"
  )

  if (nrow(outc_rt) == 0) {
    next
  }
  harm_rt <- harmonise_data(
    exposure_dat = gene_snp, 
    outcome_dat = outc_rt, action = 1
  )

  mr_result <- mr(harm_rt)
  mr_result$Gene <- gene_name
  results_list[[gene_name]] <- mr_result
}
final_results <- do.call(rbind, results_list)
library(openxlsx)
write.xlsx(final_results, "mr_results1.xlsx", rowNames = FALSE)
View(all_results)
write.csv(all_results, "mr_results.csv", row.names = FALSE)


