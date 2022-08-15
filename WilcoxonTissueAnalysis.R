# to conduct wilcoxon rank test to compare pathway activities between COVID-19 and non-COVID-19 patients in different tissues
# proteomics data was downloaded from the supplementary files from Nie. X et al (2021)
# Nhung, 12 April 2022

library(clusterProfiler)
library(readr)
library(data.table)
library(ggplot2)
library(dplyr)
library(tibble)
library(data.table)
library(ggridges)
library(org.Hs.eg.db)
library(rWikiPathways)
library(RCy3)
library(RColorBrewer)
library(ggpubr)
library(msigdbr)
library(purrr)
library(gplots)
library(ggstatsplot)
library(coin)
library(rstatix)
library(ggplot2)

source("MultiOrganProteomicFunction.R")

# create tissue data frame to be used in the whole analysis. All results will be in this order 
tissue_df <- data.frame(tissue = c("lung","spleen","liver","heart","kidney","testi","thyroid")) #renal medulla and cortex are grouped to kidney

## All COVID patients: load data, map to WikiPathways
dataname <- "Data/mmc2_covid.csv" 
df_label <- read.csv(dataname) # to extract patient labels
covid_label <- df_label[,2]
df_data <- read.csv('Data/mmc3.csv') # protein matrix of all covid and noncovid patients
covid_data <- subset(df_data[,c('Uniprot.ID', "Gene.name", covid_label)])

# Map to pathway databases
# Load pathway database: gmt file from wikipathways and the disease map
wp_file <- "Data/wikipathways-20210910-gmt-Homo_sapiens.gmt"
dm_file <- "Data/COVID19_DiseaseMap_June2021.gmt"
pwys <- c("WP4936","WP5020","WP5021","WP5035") # extra pathways in the wikipathways covid portal but not in the gmt file
all_pathway <- pathway_ID(covid_data, wp_file, dm_file, pwys)

## nonCOVID patients: load data, map to wikipathways

dataname <- "Data/mmc2_noncovid.csv"
noncovid_df_label <- read.csv(dataname) # to extract patient labels
non_covid_label <- noncovid_df_label[,2]
non_covid_data <- subset(df_data[,c('Uniprot.ID', "Gene.name", non_covid_label)])
nc_all_pathway <- pathway_ID(non_covid_data, wp_file, dm_file, pwys)

# extract pathways with more than 3 genes and 30% gene measured 

ProteinMean <- protein_mean(covid_data,tissue_df)
all_pathway_mean <- pathway_ID(ProteinMean[1], wp_file, dm_file, pwys)

NC_ProteinMean <- protein_mean(non_covid_data,tissue_df)
NC_all_pathway_mean <- pathway_ID(NC_ProteinMean[1], wp_file, dm_file, pwys)

n <-  length(tissue_df$tissue)
m1 <- 4 # (the first tissue column in all_pathway_mean)
m2 <- 10 # (the last tissue column in all_pathway_mean)
tol1 <- 30 # at least 30% gene measured per pathway
tol2 <- 3 # at least 3 genes measured per pathway
measured_pathway_mean <- tol_measured_pathway(all_pathway_mean,n,m1,m2,tol1,tol2)
NC_measured_pathways_mean <- tol_measured_pathway(NC_all_pathway_mean,n,m1,m2,tol1,tol2)
length(unique(measured_pathway_mean$wpid))
length(unique(NC_measured_pathways_mean$wpid))
overlap_pathway <- intersect(measured_pathway_mean$wpid, NC_measured_pathways_mean$wpid)

# Calculate Wilcoxon rank test for each pathway in COVID-19 and non-COVID-19
patients <- data.frame(matrix(NA, nrow = length(tissue_df$tissue), ncol = 3))
colnames(patients) <- c('tissue', 'covid', 'noncovid')

for (t in 1:length(tissue_df$tissue)){
  organ <- tissue_df$tissue[t]
  # organ_data <- all_pathway %>% dplyr::select(contains(organ))
  # dim(organ_data) # 30 patients in lung 
  # patients$covid[t] <- length(unique(colnames(organ_data)))
  # compare gene expression distribution for noncovid and covid for all the same pathway
  wilcoxon_all_gene_res <- data.frame((matrix(NA, ncol = 5, nrow = length(overlap_pathway))))
  colnames(wilcoxon_all_gene_res) <- c('pathway', 'p_value', 'effec_size', 'magnitude', 'difference')
  wilcoxon_all_gene_res$pathway <- overlap_pathway
  for (i in 1:length(overlap_pathway)){
    tryCatch({
      print(paste(overlap_pathway[i], 'in', tissue_df$tissue[t] ))
      nc_all_gene_pathway1 <- nc_all_pathway[which(nc_all_pathway$wpid == overlap_pathway[i]),]
      nc_tissue_data <- nc_all_gene_pathway1 %>% dplyr::select(contains(organ))
      nc_tissue_data_stack <- stack(nc_tissue_data)
      
      all_gene_pathway1 <- all_pathway[which(all_pathway$wpid == overlap_pathway[i]),]
      tissue_data <- all_gene_pathway1 %>% dplyr :: select(contains(organ))
      tissue_data_stack <- stack(tissue_data)
      
      patients$tissue [t] <- tissue_df$tissue[t]
      patients$covid [t] <- dim(tissue_data)[2]
      patients$noncovid [t] <- dim(nc_tissue_data)[2]
      
      dat_all_gene <- data.frame(
        Group = as.factor(c(rep("covid", length(tissue_data_stack$values)), rep("noncovid", length(nc_tissue_data_stack$values)))),
        protein_value = c(as.numeric(tissue_data_stack$values), as.numeric(nc_tissue_data_stack$values)
        )
      )
      test <- wilcox.test(protein_value ~ Group, data = dat_all_gene, conf.int = TRUE, paired = FALSE)
      wilcoxon_all_gene_res$p_value [i] <- test$p.value
      wilcoxon_all_gene_res$difference [i] <- test$estimate #the median of the difference between a sample from x and a sample from y.
      effec_res <- dat_all_gene %>% wilcox_effsize(protein_value ~ Group, paired = FALSE )
      wilcoxon_all_gene_res$effec_size [i] <- effec_res$effsize # The effect size r is calculated as Z statistic divided by square root of the sample size (N) (\(Z/\sqrt{N}\))
      wilcoxon_all_gene_res$magnitude [i] <- effec_res$magnitude
      
      ggbetweenstats( # independent samples
        data = dat_all_gene,
        x = Group,
        y = protein_value,
        plot.type = "box", # for boxplot
        type = "nonparametric", # for wilcoxon , func used wilcox.test(x, ...)
        centrality.plotting = FALSE # remove median
      )
      
      path <- gsub(" ","",paste("Wilcoxon_Fig/",tissue_df$tissue[t]))
      dir.create(path)
      filename <- gsub(" ","",paste(path, "/", overlap_pathway[i], ".png"))
      # ggsave(paste("Wilcoxon_Fig/lung/all_gene/", overlap_pathway[i], ".png"))
      ggsave(filename)
    }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
  }
  
  write_csv(wilcoxon_all_gene_res, paste(gsub(" ","", paste(path, "/",tissue_df$tissue[t], 'all_gene_wilcoxon.csv'))))
  
  sig_p <- wilcoxon_all_gene_res[which(wilcoxon_all_gene_res$p_value < 0.05),]
  length(sig_p$pathway)
  sig_big_ef <- sig_p[which(sig_p$magnitude > 2),]
  sig_big_ef 
  sig_rank <- sig_p[order(-sig_p$effec_size),]
  sig_rank
  write_csv(sig_rank, paste(gsub(" ","", paste(path, "/",tissue_df$tissue[t], 'significant_all_gene.csv'))))
} 

## Analyze result
heart <- read.csv('Wilcoxon_Fig/heart/heartsignificant_all_gene.csv')
lung <- read.csv("Wilcoxon_Fig/lung/lungsignificant_all_gene.csv")
spleen <- read.csv("Wilcoxon_Fig/spleen/spleensignificant_all_gene.csv")
thyroid <- read.csv("Wilcoxon_Fig/thyroid/thyroidsignificant_all_gene.csv")
testi <- read.csv("Wilcoxon_Fig/testi/testisignificant_all_gene.csv")
liver <- read.csv("Wilcoxon_Fig/liver/liversignificant_all_gene.csv")
kidney <- read.csv("Wilcoxon_Fig/kidney/kidneysignificant_all_gene.csv")

# only keep pathways with more than 30% genes and 3 genes measured per tissue
n <-  length(tissue_df$tissue)
m1 <- 4 # (the first tissue column in all_pathway_mean)
m2 <- 10 # (the last tissue column in all_pathway_mean)
tol1 <- 30 # at least 30% gene measured per pathway
tol2 <- 3 # at least 3 genes measured per pathway
table <- prot_measured_per_pathway(all_pathway_mean, n,m1,m2)

number_prot <- data.frame(table[[1]])
percent_prot <- data.frame(table[[2]])

heart_tol <- percent_prot$wpid[which(as.numeric(number_prot$heart) >3 & as.numeric(percent_prot$heart) >0.3)]
heart_limit <- filter(heart, pathway %in% heart_tol)

lung_tol <- percent_prot$wpid[which(as.numeric(number_prot$lung) >3 & as.numeric(percent_prot$lung) >0.3)]
lung_limit <- filter(lung, pathway %in% lung_tol)

kidney_tol <- percent_prot$wpid[which(as.numeric(number_prot$kidney) >3 & as.numeric(percent_prot$kidney) >0.3)]
kidney_limit <- filter(kidney, pathway %in% kidney_tol)

thyroid_tol <- percent_prot$wpid[which(as.numeric(number_prot$thyroid) >3 & as.numeric(percent_prot$thyroid) >0.3)]
thyroid_limit <- filter(thyroid, pathway %in% thyroid_tol)

spleen_tol <- percent_prot$wpid[which(as.numeric(number_prot$spleen) >3 & as.numeric(percent_prot$spleen) >0.3)]
spleen_limit <- filter(spleen, pathway %in% spleen_tol)

liver_tol <- percent_prot$wpid[which(as.numeric(number_prot$liver) >3 & as.numeric(percent_prot$liver) >0.3)]
liver_limit <- filter(liver, pathway %in% liver_tol)

testi_tol <- percent_prot$wpid[which(as.numeric(number_prot$testi) >3 & as.numeric(percent_prot$testi) >0.3)]
testi_limit <- filter(testi, pathway %in% testi_tol)

# extract the most significant changed pathways per tissue. pval < 0.05 and effect size > 0.3 or magnitude > 1
top_kidney <- kidney_limit[which(kidney_limit$magnitude >1),]
top_lung <- lung_limit[which(lung_limit$magnitude >1),]
top_liver <- liver_limit[which(liver_limit$magnitude >1),]
top_spleen <- spleen_limit[which(spleen_limit$magnitude >1),]
top_heart <- heart_limit[which(heart_limit$magnitude >1),]
top_testi <- testi_limit[which(testi_limit$magnitude >1),]
top_thyroid <- thyroid_limit[which(thyroid_limit$magnitude >1),]

# create a matrix with all significant changed pathways with their values (p_value, effect size, maginitude, difference) in each tissue
all_sig_p <- cbind(top_kidney$pathway, top_liver$pathway, top_lung$pathway, top_thyroid$pathway,
                   top_testi$pathway, top_spleen$pathway, top_heart$pathway)
all_sig_p <- stack(all_sig_p)
all_sig_p <- unique(all_sig_p$value) # 343 significant pathways, 55 pathways with size effect > 0.3
matrix <- data.frame(matrix(NA, ncol = 7, nrow = length(all_sig_p)))
colnames(matrix) <- tissue_df$tissue

for (i in 1:length(all_sig_p)){
  # lung
  lung_index <- which(lung_limit$pathway == all_sig_p[i])
  if (is_empty(lung_index)) { matrix$lung [i] <- NA} else {
    if (lung_limit$difference[lung_index] > 0) { 
      matrix$lung [i] <- lung_limit$effec_size[lung_index]
    } else { matrix$lung [i] <- lung_limit$effec_size[lung_index] * -1}
  }
  
  # liver
  liver_index <- which(liver_limit$pathway == all_sig_p[i])
  if (is_empty(liver_index)) { matrix$liver [i] <- NA} else {
    if (liver_limit$difference[liver_index] > 0) { 
      matrix$liver [i] <- liver_limit$effec_size[liver_index]
    } else { matrix$liver [i] <- liver_limit$effec_size[liver_index] * -1}
  }
  
  # spleen
  spleen_index <- which(spleen_limit$pathway == all_sig_p[i])
  if (is_empty(spleen_index)) { matrix$spleen [i] <- NA} else {
    if (spleen_limit$difference[spleen_index] > 0) { 
      matrix$spleen [i] <- spleen_limit$effec_size[spleen_index]
    } else { matrix$spleen [i] <- spleen_limit$effec_size[spleen_index] * -1}
  }
  
  # heart
  heart_index <- which(heart_limit$pathway == all_sig_p[i])
  if (is_empty(heart_index)) { matrix$heart [i] <- NA} else {
    if (heart_limit$difference[heart_index] > 0) { 
      matrix$heart [i] <- heart_limit$effec_size[heart_index]
    } else { matrix$heart [i] <- heart_limit$effec_size[heart_index] * -1}
  }
  
  # testi
  testi_index <- which(testi_limit$pathway == all_sig_p[i])
  if (is_empty(testi_index)) { matrix$testi [i] <- NA} else {
    if (testi_limit$difference[testi_index] > 0) { 
      matrix$testi [i] <- testi_limit$effec_size[testi_index]
    } else { matrix$testi [i] <- testi_limit$effec_size[testi_index] * -1}
  }
  
  #thyroid
  thyroid_index <- which(thyroid_limit$pathway == all_sig_p[i])
  if (is_empty(thyroid_index)) { matrix$thyroid [i] <- NA} else {
    if (thyroid_limit$difference[thyroid_index] > 0) { 
      matrix$thyroid [i] <- thyroid_limit$effec_size[thyroid_index]
    } else { matrix$thyroid [i] <- thyroid_limit$effec_size[thyroid_index] * -1}
  }
  
  # kidney
  kidney_index <- which(kidney_limit$pathway == all_sig_p[i])
  if (is_empty(kidney_index)) { matrix$kidney [i] <- NA} else {
    if (kidney_limit$difference[kidney_index] > 0) { 
      matrix$kidney [i] <- kidney_limit$effec_size[kidney_index]
    } else { matrix$kidney [i] <- kidney_limit$effec_size[kidney_index] * -1}
  }
}
rownames(matrix) <- all_sig_p

# heatmap
# gradient colors
stack_data3 <- stack(matrix)
rank_data3 <- stack_data3[order(-stack_data3$values),]
rank_data3 <- na.omit(rank_data3)
breaks = seq(min(rank_data3$values), max(rank_data3$values), length.out= length(rank_data3$values))

gradient1 <- colorpanel(sum(breaks[-1] < -0.3), "#003993", "#51a7db" ) # down regulated pathways
gradient2 <- colorpanel(sum(breaks[-1] >= -0.3 & breaks[-1] <= 0), "#51a7db", "#ffffff") # no change pathways
gradient3 <- colorpanel(sum(breaks[-1] > 0 & breaks[-1] <= 0.3), "#ffffff" , "#db8151") # no change pathways
gradient4 <- colorpanel(sum(breaks[-1] > 0.3), "#db8151", '#930000') # upregulated pathways

hm.colors = c(gradient1, gradient2, gradient3, gradient4)

matrix2 <- matrix(unlist(matrix),nrow= nrow(matrix), ncol = 7)
colnames(matrix2) <- tissue_df$tissue

heatmap.2(matrix2, na.color = 'grey', distfun = dist_no_na, col= hm.colors, breaks = breaks, trace = "none")

# Analyze heat map result
group_up <- c()
group_down <- c()
top_red <- all_sig_p[group_up]
up_pathways <- filter(all_pathway,wpid %in% top_red)
unique(up_pathways$name)

bottom_blue <- all_sig_p[group_down]
down_blue <- filter(all_pathway,wpid %in%bottom_blue)
unique(down_blue$name)

## Analyze overlap among tissues 

all_rank <- data.frame(matrix(NA, ncol = 5, nrow = length(all_sig_p)))
colnames(all_rank) <- c('sig_pathway', 'all_tissues','up', 'down', 'abs_effect_size')
all_rank[,1] <- all_sig_p

for (i in 1:length(all_sig_p)){
  all_rank [i,2] <- rowSums(!is.na(matrix[i,])) # number of tissue this pathway change
  all_rank [i,3] <- rowSums(matrix[i,] >0, na.rm = TRUE) # number of tissue this pathway is up
  all_rank [i,4] <- rowSums(matrix[i,] <0, na.rm = TRUE) # number of tissue this pathway is down
  all_rank [i,5] <- rowSums(abs(matrix[i,]) > 0.3, na.rm = TRUE) # number of tissue this pathway has size effect > 0.3
}

all_rank_order <- all_rank[order(-all_rank$abs_effect_size),]
all_rank_order
length(all_sig_p) # 73 out of 736 pathways that change significantly in at least 1 tissue # 69?

# Pathways that change in more than 2 tissues with effect size > 0.3
top_p <- all_rank$sig_pathway[which(all_rank$abs_effect_size > 1)]
top_change_p <- filter(all_pathway, wpid%in% top_p)
unique(top_change_p$name)

# write to table
multi_tissue_p <- data.frame(matrix(NA, ncol = 9, nrow = length(top_p)))
colnames(multi_tissue_p) <- c('wpid', 'name', tissue_df$tissue)
for (i in 1:length(top_p)){
  index <- which(all_sig_p == top_p[i])
  multi_tissue_p[i, 3:9] <- matrix[index,]
  multi_tissue_p[i,1] <- top_p[i]
  multi_tissue_p[i,2] <- unique(all_pathway$name[which(all_pathway$wpid == top_p[i])])
}

write_csv(multi_tissue_p, 'Results/15_pathways_altered_in_more_than_1_tissue.csv')
# bar plot for number of changed pathways per tissue

result <- matrix(NA, ncol = 7, nrow = 2)
result [1,1] <- length(which(top_lung$difference > 0))
result [2,1] <- length(which(top_lung$difference < 0))
result [1,2] <- length(which(top_spleen$difference > 0))
result [2,2] <- length(which(top_spleen$difference < 0))
result [1,3] <- length(which(top_liver$difference > 0))
result [2,3] <- length(which(top_liver$difference < 0))
result [1,4] <- length(which(top_heart$difference > 0))
result[2,4] <- length(which(top_heart$difference < 0))
result [1,5] <- length(which(top_kidney$difference> 0))
result [2,5] <- length(which(top_kidney$difference < 0))
result [1,6] <- length(which(top_testi$difference > 0))
result [2,6] <- length(which(top_testi$difference < 0))
result [1,7] <- length(which(top_thyroid$difference > 0))
result [2,7] <- length(which(top_thyroid$difference < 0))

rownames(result) <- c("up", "down")
colnames(result) <- tissue_df$tissue
result

cut_off <- 0.3
layout(matrix(c(1,1,2,3,4,5,6,7,8), nrow = 3, ncol = 3, byrow = TRUE))
barplot(result, main = "The number of altered pathways in COVID-19 patients", names.arg =  colnames(result),
        col= c('#920000', '#2588b3'),beside = TRUE, border= NA, cex.lab = 1.5,cex.axis = 1.5, cex.names = 1.5, cex.main = 1.5,
        ylab = " Number of pathways")

legend("topleft",
       c("More active in COVID-19","Less active in COVID-19"), cex = 1.5,bty = 'n',
       fill = c('#920000','#2588b3'))
for (i in 1: length(tissue_df$tissue)){
  tissue <- tissue_df$tissue[i]
  plot(matrix[,i], col = ifelse(matrix[,i] >0, '#920000', '#2588b3'),ylim =c(-0.7, 0.7) , pch =20, cex = 2.5,
       cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5,
       main = tissue, xlab ='Pathway index', ylab='Effect size')
  abline(h= cut_off, col="black", lty = 2, lwd = 1.5)
  abline(h= -cut_off, col='black',lty = 2, lwd = 1.5)
}            # save image as sgv and add legend using inkscape

## ANalyze result for each tissue

# testi
for (i in 1:length(testi_limit$pathway)){
  testi_limit$name [i] <- unique(all_pathway$name[which(all_pathway$wpid ==testi_limit$pathway[i])])
}

testi_limit


# liver
for (i in 1:length(liver_limit$pathway)){
  liver_limit$name [i] <- unique(all_pathway$name[which(all_pathway$wpid ==liver_limit$pathway[i])])
}

liver_limit

# thyroid
for (i in 1:length(thyroid_limit$pathway)){
  thyroid_limit$name [i] <- unique(all_pathway$name[which(all_pathway$wpid ==thyroid_limit$pathway[i])])
}

thyroid_limit

# lung
for (i in 1:length(lung_limit$pathway)){
  lung_limit$name [i] <- unique(all_pathway$name[which(all_pathway$wpid ==lung_limit$pathway[i])])
}

lung_limit

# heart
for (i in 1:length(heart_limit$pathway)){
  heart_limit$name [i] <- unique(all_pathway$name[which(all_pathway$wpid ==heart_limit$pathway[i])])
}

heart_limit

# kidney
for (i in 1:length(kidney_limit$pathway)){
  kidney_limit$name [i] <- unique(all_pathway$name[which(all_pathway$wpid ==kidney_limit$pathway[i])])
}

kidney_limit

# spleen 
for (i in 1:length(spleen_limit$pathway)){
  spleen_limit$name [i] <- unique(all_pathway$name[which(all_pathway$wpid ==spleen_limit$pathway[i])])
}

spleen_limit
top_check <- spleen_limit[1:17,]
top_check[which(top_check$difference <0),]

# write to file

write_csv(testi_limit,"Results/testi_Changed.csv")
write_csv(liver_limit,"Results/liver_Changed.csv")
write_csv(lung_limit,"Results/lung_Changed.csv")
write_csv(heart_limit,"Results/heart_Changed.csv")
write_csv(thyroid_limit,"Results/thyroid_Changed.csv")
write_csv(spleen_limit,"Results/spleen_Changed.csv")
write_csv(kidney_limit,"Results/kidney_Changed.csv")

## Analyze thyroid

thyroid_name <- filter(all_pathway, wpid %in% thyroid$pathway[1:73])
unique(thyroid_name$name)

# create visualization for all sig  pathways in thyroid
# i = 7 # thyroid location in the tissue_df 
# for (j in 6: length(top_thyroid$pathway)) {
#   pathway_to_check <- top_thyroid$pathway[j]
#   if (startsWith(pathway_to_check, "WP")){
#     print (pathway_to_check)
#     unique(all_pathway$name[which(all_pathway_mean$wpid == pathway_to_check)])
#     pathway_data <- data.frame(all_pathway_mean[which(all_pathway_mean$wpid == pathway_to_check),],
#                                NC_all_pathway_mean[which(NC_all_pathway_mean$wpid == pathway_to_check),])
#     check_columns <- c(tissue_df$tissue[i],gsub(" ","",paste(tissue_df$tissue[i], ".1")))
#     RCy3::commandsRun(paste('wikipathways import-as-pathway id=', pathway_to_check))
#     toggleGraphicsDetails()
#     loadTableData(pathway_data, data.key.column = "ENSEMBL", table.key.column = "Ensembl")
#     
#     setNodeCustomBarChart(check_columns, type = "GROUPED", colors = c("red","blue"), orientation = "HORIZONTAL", style.name = "WikiPathways")
#     # Saving output
#     path <- gsub(" ","",paste("Wilcoxon_res_fig/",tissue_df$tissue[i]))
#     dir.create(path)
#     filename <- gsub(" ","",paste(path, "/",pathway_to_check))
#     exportImage(filename,'SVG')
#     exportImage(filename,'PNG', zoom = 500)
#     saveSession(filename)
#     RCy3::closeSession(save.before.closing = F)
#   }
# }

## Visualize most altered pathways per tissue. Make sure to have cytoscape open in advance
dir.create("Wilcoxon_res_fig")
for (i in 1: length(tissue_df$tissue)){
  print (tissue_df$tissue[i])
  index <- which(abs(matrix[,i]) > 0.3)
  for (j in 1: length(index)) {
    pathway_to_check <- all_sig_p[index[j]]
    if (startsWith(pathway_to_check, "WP")){
      print (pathway_to_check)
      unique(all_pathway$name[which(all_pathway_mean$wpid == pathway_to_check)])
      pathway_data <- data.frame(all_pathway_mean[which(all_pathway_mean$wpid == pathway_to_check),],
                                 NC_all_pathway_mean[which(NC_all_pathway_mean$wpid == pathway_to_check),])
      check_columns <- c(tissue_df$tissue[i],gsub(" ","",paste(tissue_df$tissue[i], ".1")))
      RCy3::commandsRun(paste('wikipathways import-as-pathway id=', pathway_to_check))
      toggleGraphicsDetails()
      loadTableData(pathway_data, data.key.column = "ENSEMBL", table.key.column = "Ensembl")
      setNodeCustomBarChart(check_columns, type = "GROUPED", colors = c("red","blue"), orientation = "HORIZONTAL", style.name = "WikiPathways")
      # Saving output
      path <- gsub(" ","",paste("Wilcoxon_res_fig/",tissue_df$tissue[i]))
      dir.create(path)
      filename <- gsub(" ","",paste(path, "/",pathway_to_check))
      exportImage(filename,'SVG')
      exportImage(filename,'PNG', zoom = 500)
      saveSession(filename)
      RCy3::closeSession(save.before.closing = F)
    }
  }
}

# write change pathways into file for each tissue
dir.create("Wilcoxon_Fig")

for (i in 1:length(tissue_df$tissue)){
  print(tissue_df$tissue)
  index <- which(abs(matrix[,i]) > 0.3)
  table <- data.frame()
  for (j in 1:length(index)){
    name <- all_pathway$name[which(all_pathway$wpid == all_sig_p[index[j]])]
    table [j,1] <- all_sig_p[index[j]]
    table [j,2] <- name[1]
    table [j,3] <- matrix[index[j],i]
  }
  colnames(table) <- c("id", "name", "difference")
  write_csv(table,paste("Wilcoxon_Fig/",Sys.Date(),"changed pathways in", tissue_df$tissue[i], ".csv"))
  
}

## Create network of all changed pathways 

changed_pathways <- filter(all_pathway,wpid %in% all_sig_p)
dir.create("Results")
index <- !is.na(changed_pathways$Uniprot.ID)
node_list <- data.frame(changed_pathways$Uniprot.ID[index], changed_pathways$wpid[index])
write_csv(node_list, "Results/node_list2") # create network in cytoscape from this file

# prepare table data to load to the network in cytoscape
pathway_list <- unique(changed_pathways$wpid) 
database <- data.frame(matrix(,ncol = 1, nrow = length(pathway_list)))
change_dif <- data.frame(matrix(,ncol = 7, nrow = length(pathway_list)))
colnames(change_dif) <- c("plung","pspleen","pliver","pheart","pkidney","ptesti","pthyroid")

for (i in 1:length(pathway_list)){
  index <- which(all_sig_p == pathway_list[i])
  change_dif$plung[i] <- matrix[index,1]
  change_dif$pspleen[i] <- matrix[index,2]
  change_dif$pliver[i] <- matrix[index,3]
  change_dif$pheart[i] <- matrix[index,4]
  change_dif$pkidney[i] <- matrix[index,5]
  change_dif$ptesti[i] <- matrix[index,6]
  change_dif$pthyroid[i] <- matrix[index,7]
  if (startsWith(as.character(pathway_list[i]), "WP")){
    database$matrix...ncol...1..nrow...length.pathway_list..[i] <- "wikipathways"
  } else if (startsWith(as.character(pathway_list[i]),"MINERVA")){
    database$matrix...ncol...1..nrow...length.pathway_list..[i] <- "DiseaseMap"}
}

all_table <- data.frame(pathway_list, database,change_dif)
colnames(all_table)[1:2] <- c("wpid","database")
write_csv(all_table,"Results/all_pathways_with_dif2.csv") # load this table data to the network and adjust visualization

## pathways that chnaged in literature

lit <- c('WP4868', 'WP366', 'WP3333', 'WP3858', 'WP30360', 'WP400', 'WP364')
for (i in 1:length(lit)){
print(which(all_sig_p == lit[i])) # there is none
}
# inspect pathway of interest with differential expression data
# load differential expression data
dif_data <- read.csv('Data/mmc4_covidvsnoncovid.csv')
dif_data <-  dif_data[2:nrow(dif_data),] # remove header line
# data order: log2FC, pvalue, adjusted pvalue  (significant adjusted p_value < 0.05, |log2FC| > log2[1.2])
lung <- c(3,4,5) # the tissue columns in the differential expression data
spleen <- c(6,7,8)
liver <- c(9,10,11)
heart <- c(12,13,14)
renal_cortex <- c(15,16,17)
renal_medulla <- c(18,19,20)
kidney <- c(15:20)
testi <- c(21,22,23)
thyroid <- c(24,25,26) 

pathway_to_check <- "WP4586"
all_pathway$name[which(all_pathway$wpid == pathway_to_check)][1]
matrix[which(all_sig_p == pathway_to_check),]
i <- 7 # tissue position in tissue_df
tissue_col <- thyroid
inspect_pathway(number_prot, pathway_to_check, i, tissue_col, all_pathway, nc_all_pathway)





