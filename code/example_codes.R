
library(dplyr)

### Data

all.drugs <- c("Asparaginase", "Bortezomib", "CHZ868", "Cytarabine", "Dasatinib", "Daunorubicin", "Dexamethasone", "Ibrutinib", "Mercaptopurine", "Nelarabine", "Panobinostat", "Prednisolone", "Ruxolitinib", "Thioguanine", "Trametinib", "Venetoclax", "Vincristine", "Vorinostat")

lc50.data <- read.csv('phenotypes.csv', check.names=FALSE)

### clean up drug names, removing units
names(lc50.data) <- gsub(' +\\(.*\\)$', '', names(lc50.data))
table(all.drugs %in% names(lc50.data))

drug.concentrations <- read.csv('data/drug_concentration.csv')

table(drug.concentrations$drug %in% names(lc50.data))

## Normalization of LC50

# lc50 = (log(lc50) - log(min_lc50)) / (log(max_lc50) - log(min_lc50))

### remove the existing columns with normalized data
lc50.data <- lc50.data[, !grepl('_normalized', names(lc50.data))]

for (curr.drug in all.drugs) {
    print(curr.drug)
    curr.lc50 <- lc50.data[[curr.drug]]
    curr.conc.range <- subset(drug.concentrations, drug==curr.drug)
    max_conc <- curr.conc.range$max_conc[1]
    min_conc <- curr.conc.range$min_conc[1]
    
    curr.lc50.normalized <- (log(curr.lc50) - log(curr.conc.range$min_conc)) / (log(max_conc) - log(min_conc))
    curr.lc50.normalized[curr.lc50.normalized <= 0]  <- 0
    curr.lc50.normalized[curr.lc50.normalized >= 1]  <- 1
    lc50.data[[paste0(curr.drug, '_normalized')]] <- curr.lc50.normalized
}

## Multiple imputation

### prepare data for imputation, use protocol to adjust for potential batch effect

lc50.data$subjectid <- 1:nrow(lc50.data)

lc50.data.for.imputation <- lc50.data[, c('subjectid', paste0(all.drugs, '_normalized'))]

lc50.data.for.imputation <- cbind(lc50.data.for.imputation, model.matrix(~Protocol, data=lc50.data)[,-1])


library(mice)

n.impute <- 10

### first round of imputation using MICE

## linear model
imputed.lc50.data.orig <- mice(lc50.data.for.imputation[, -1], m=n.impute, maxit = 5, printFlag=FALSE, ntree=10)

imputed.lc50.data <- complete(imputed.lc50.data.orig, "long")

## ----mice_r2, include=FALSE, cache=TRUE---------------------------------------

### continue to impute the variables excluded in r1 due to collinearity

failed.drugs.r1 <- names(imputed.lc50.data)[colSums(is.na(imputed.lc50.data)) > 100]

failed.drugs.r1

for (i in 1:n.impute) {
  temp.data <- imputed.lc50.data %>% filter(.imp==i)
  
  ### linear model
  
  imputed.lc50.data.orig <- mice(temp.data[, -c(1,2)], m=5, maxit = 5, printFlag=FALSE, ntree=10)
  
  temp.long <- complete(imputed.lc50.data.orig)
  
  temp.long$.imp <- i
  temp.long$.id <- 1:nrow(temp.long)
  
  imputed.lc50.data <- imputed.lc50.data %>% filter(.imp!=i) %>% bind_rows(temp.long)
}

failed.drugs.r2 <- names(imputed.lc50.data)[colSums(is.na(imputed.lc50.data)) > 100]

print(failed.drugs.r2)

### Clustering of patients based on the imputed lc50s

library(ComplexHeatmap)
library(circlize)
library(dendsort)

names(imputed.lc50.data) <- gsub('_normalized', '', names(imputed.lc50.data))

mat2 <- as.matrix(imputed.lc50.data[, all.drugs])
dim(mat2)

tmp.clust <- hclust(dist(mat2, method='manhattan'), method='ward.D2')

tmp.clust.dendro <- tmp.clust %>% as.dendrogram %>% dendsort

n.stacking.clust <- 6

ht_list = Heatmap(mat2, col = colorRamp2(c(-0.2, 0.5, 1.2), c("blue", "white", "red")), 
    name = "lc50", column_title = paste0("LC50 for ", nrow(mat2)/n.impute, " Patients (stacking)"),
    show_column_names = TRUE, show_row_names=FALSE, width = unit(9, "cm"),
    clustering_distance_columns = "pearson",  clustering_method_columns='ward.D',
    cluster_rows=tmp.clust.dendro,
    heatmap_legend_param = list(title = "LC50"), show_heatmap_legend=FALSE,
    left_annotation = rowAnnotation(block = anno_block(gp = gpar(fill = 2:6, col = NA)),width=unit(3,'mm')),
    row_split = n.stacking.clust)

draw(ht_list)

###

sessionInfo()




