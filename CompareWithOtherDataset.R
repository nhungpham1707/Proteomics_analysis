## THe same pathway activities was carried out to compare with the tissue specific proteomic
# data was collected from blood in COVID-19 patients with different severity. control group is healthy people 
# https://www.nature.com/articles/s41598-021-85877-0#data-availability 

library(rstatix)
library(org.Hs.eg.db)
library(rWikiPathways)
library(RCy3)
library(readr)

source("MultiOrganProteomicFunction.R")
# read the edit data

blood_data <- data.frame(read_csv("Data/S-BSST416(1)/blood_proteomics_edit.csv"))
blood_data
# gene_list <- colnames(blood_data)[2:ncol(blood_data)] #368 genes

blood_data2 <- t(blood_data)
blood_data3 <- blood_data2[2:nrow(blood_data2),] # remove the group row
blood_data4 <- data.frame(rownames(blood_data3), blood_data3)
colnames(blood_data4)[1] <- "uniprot"
blood_data4

mild_index <- which(blood_data2[1,] == 1)
moderate_index <- which(blood_data2[1,] == 2)
severe_index <- which(blood_data2[1,]==3)
control_index <- which(blood_data2[1,] == 0)

# map to pathway 
wp_file <- "Data/wikipathways-20210910-gmt-Homo_sapiens.gmt"
dm_file <- "Data/COVID19_DiseaseMap_June2021.gmt"
pwys <- c("WP4936","WP5020","WP5021","WP5035") # extra pathways in wikipathways covid portal but not in the gmt file
unip2entrez <- clusterProfiler::bitr(blood_data4$uniprot, fromType = "UNIPROT",toType = c("ENTREZID","SYMBOL","ENSEMBL"), OrgDb = org.Hs.eg.db)
data <- merge(blood_data4, unip2entrez, by.x="uniprot", by.y="UNIPROT", all.x = TRUE)

# Generate pathway gene list from wikipathway and diseasemap 
pwy2gene <- combineWP_DM(wp_file, dm_file,pwys)
# merge gene in dataset with pathway-gene list in wikipathways and disease map through entrezid
pathway <- merge(data,pwy2gene, by.x = "ENTREZID", by.y = 'gene', all.y = TRUE)

unique_pathway <- unique(pathway$wpid) #726 pathways

control_data <- data.frame(pathway[,1:2],pathway[,control_index+2], pathway$SYMBOL, pathway$ENSEMBL, pathway$name, pathway$wpid)
control_data

mild_data <- data.frame(pathway[,1:2],pathway[,mild_index+2], pathway$SYMBOL, pathway$ENSEMBL, pathway$name, pathway$wpid)
mild_data

moderate_data <- data.frame(pathway[,1:2],pathway[,moderate_index+2], pathway$SYMBOL, pathway$ENSEMBL, pathway$name, pathway$wpid)
moderate_data

severe_data <- data.frame(pathway[,1:2],pathway[,severe_index+2], pathway$SYMBOL, pathway$ENSEMBL, pathway$name, pathway$wpid)
severe_data

# calculate wilcoxon 
pathway_to_compare <- c("WP4868", "WP366", "WP2332", "WP364")
  wilcoxon_mild_res <- data.frame((matrix(NA, ncol = 5, nrow = length(pathway_to_compare))))
  colnames(wilcoxon_mild_res) <- c('pathway', 'p_value', 'effec_size', 'magnitude', 'difference')
  wilcoxon_mild_res$pathway <- pathway_to_compare
  
  wilcoxon_moderate_res <- data.frame((matrix(NA, ncol = 5, nrow = length(pathway_to_compare))))
  colnames(wilcoxon_moderate_res) <- c('pathway', 'p_value', 'effec_size', 'magnitude', 'difference')
  wilcoxon_moderate_res$pathway <- pathway_to_compare
  
  wilcoxon_severe_res <- data.frame((matrix(NA, ncol = 5, nrow = length(pathway_to_compare))))
  colnames(wilcoxon_severe_res) <- c('pathway', 'p_value', 'effec_size', 'magnitude', 'difference')
  wilcoxon_severe_res$pathway <- pathway_to_compare
  
  
i = 1

    print(paste(pathway_to_compare[i])
          control_pathway1 <- control_data[which(control_data$pathway.wpid == pathway_to_compare[i]),]
          control_pathway1 
          control_col <- ncol(control_pathway1)-4
          control_stack <- stack(control_pathway1[,3:control_col])
          control_stack
          
          mild_pathway1 <- mild_data[which(mild_data$pathway.wpid == pathway_to_compare[i]),]
          mild_pathway1
          index_col <- ncol(mild_pathway1)-4
          mild_p1_stack <- stack(mild_pathway1[,3:index_col])
          mild_p1_stack
          
          moderate_p1 <- moderate_data[which(moderate_data$pathway.wpid == pathway_to_compare[i]),]
          moderate_p1
          index_col <- ncol(moderate_p1)-4
          moderate_p1_stack <- stack(moderate_p1[,3:index_col])
          moderate_p1_stack
          
          severe_p1 <- severe_data[which(severe_data$pathway.wpid == pathway_to_compare[i]),]
          severe_p1
          index_col <- ncol(severe_p1)-4
          severe_p1_stack <- stack(severe_p1[,3:index_col])
          severe_p1_stack
          
          
          dat_all_gene <- data.frame(
            Group = as.factor(c(rep("mild", length(mild_p1_stack$values)), rep("noncovid", length(control_stack$values)))),
            protein_value = c(as.numeric(mild_p1_stack$values), as.numeric(control_stack$values)
            )
          )
          test <- wilcox.test(protein_value ~ Group, data = dat_all_gene, conf.int = TRUE, paired = FALSE)
          wilcoxon_mild_res$p_value [i] <- test$p.value
          wilcoxon_mild_res$difference [i] <- test$estimate #the median of the difference between a sample from x and a sample from y.
          effec_res <- dat_all_gene %>% wilcox_effsize(protein_value ~ Group, paired = FALSE )
          wilcoxon_mild_res$effec_size [i] <- effec_res$effsize # The effect size r is calculated as Z statistic divided by square root of the sample size (N) (\(Z/\sqrt{N}\))
          wilcoxon_mild_res$magnitude [i] <- effec_res$magnitude
          
          dat_moderate_gene <- data.frame(
            Group = as.factor(c(rep("moderate", length(moderate_p1_stack$values)), rep("noncovid", length(control_stack$values)))),
            protein_value = c(as.numeric(moderate_p1_stack$values), as.numeric(control_stack$values)
            )
          )
          test <- wilcox.test(protein_value ~ Group, data = dat_moderate_gene, conf.int = TRUE, paired = FALSE)
          wilcoxon_moderate_res$p_value [i] <- test$p.value
          wilcoxon_moderate_res$difference [i] <- test$estimate #the median of the difference between a sample from x and a sample from y.
          effec_res <- dat_moderate_gene %>% wilcox_effsize(protein_value ~ Group, paired = FALSE )
          wilcoxon_moderate_res$effec_size [i] <- effec_res$effsize # The effect size r is calculated as Z statistic divided by square root of the sample size (N) (\(Z/\sqrt{N}\))
          wilcoxon_moderate_res$magnitude [i] <- effec_res$magnitude
          
          
          dat_severe_gene <- data.frame(
            Group = as.factor(c(rep("severe", length(severe_p1_stack$values)), rep("noncovid", length(control_stack$values)))),
            protein_value = c(as.numeric(severe_p1_stack$values), as.numeric(control_stack$values)
            )
          )
          test <- wilcox.test(protein_value ~ Group, data = dat_severe_gene, conf.int = TRUE, paired = FALSE)
          wilcoxon_severe_res$p_value [i] <- test$p.value
          wilcoxon_severe_res$difference [i] <- test$estimate #the median of the difference between a sample from x and a sample from y.
          effec_res <- dat_all_gene %>% wilcox_effsize(protein_value ~ Group, paired = FALSE )
          wilcoxon_severe_res$effec_size [i] <- effec_res$effsize # The effect size r is calculated as Z statistic divided by square root of the sample size (N) (\(Z/\sqrt{N}\))
          wilcoxon_severe_res$magnitude [i] <- effec_res$magnitude
    


# result: i= 2,3,4 do not have gene measured in control, mild. only some in modrate and severe
          # i = 1 has 5 genes (15%) 
          # WP4868 reduce in mild and moderate, increase in severe

## Cytoscape
pathway_to_check <- "WP4868"
          
severe_p <- filter(severe_data, pathway.wpid %in% pathway_to_check)
pathway_data <- data.frame(matrix(NA,ncol = 6 , nrow = length(severe_p$ENTREZID)))
colnames(pathway_data) <- c("gene", "Ensemble", "control", "mild", "moderate", "severe")
pathway_data$gene <- severe_p$uniprot
ncol <- ncol(severe_data) - 4
severe_p_data <- severe_p[,3:ncol]
pathway_data$severe <- apply(severe_p_data,1, function(x){ ifelse(sum(!is.na(x))/ncol(severe_p_data) > 0.1, log10(mean(x, na.rm = TRUE)+1), NA)})

control_p <- filter(control_data, pathway.wpid %in% pathway_to_check)
ncontrol <- ncol(control_p) - 4
control_p_data <- control_p[, 3:ncontrol]
 pathway_data$control <- apply(control_p_data,1, function(x){ ifelse(sum(!is.na(x))/ncol(control_p_data) > 0.1, log10(mean(x, na.rm = TRUE)+1), NA)})

mild_p <- filter(mild_data, pathway.wpid %in% pathway_to_check)
nmild <- ncol(mild_p) - 4
mild_p_data <- mild_p[, 3:nmild]
pathway_data$mild <- apply(mild_p_data,1, function(x){ ifelse(sum(!is.na(x))/ncol(mild_p_data) > 0.1, log10(mean(x, na.rm = TRUE)+1), NA)})

moderate_p <- filter(moderate_data, pathway.wpid %in% pathway_to_check)
nmoderate <- ncol(moderate_p) - 4
moderate_p_data <- moderate_p[, 3:nmoderate]
pathway_data$moderate <- apply(moderate_p_data,1, function(x){ ifelse(sum(!is.na(x))/ncol(moderate_p_data) > 0.1, log10(mean(x, na.rm = TRUE)+1), NA)})

pathway_data$Ensemble[2] <- 'ENSG00000213928'
pathway_data$Ensemble[3] <- 'ENSG00000213928'
pathway_data$Ensemble[4] <- 'ENSG00000107201'
pathway_data$Ensemble[16] <- 'ENSG00000198001'
pathway_data$Ensemble[21] <- 'ENSG00000130234'

           
RCy3::commandsRun(paste('wikipathways import-as-pathway id=', pathway_to_check))
toggleGraphicsDetails()
loadTableData(pathway_data, data.key.column = "Ensemble", table.key.column = "Ensembl")
# setNodeCustomBarChart(check_columns, type = "GROUPED", colors = c("red","blue"), orientation = "HORIZONTAL", style.name = "WikiPathways")
# Saving output