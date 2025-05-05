#######
# EQTL MAPPING
#######

# Generate datasets for eqtl mapping, including subsets for saturation curve
# makes covariate files for each subset
# prints cov-free matrices to be used in modularity analysis and to plot
# diff. expression results


#### BODY
{   
  #load raw counts of expressed genes
  body.info <- read.table("AllDatasets/Info_Samples_eqtlCatalog_Aug3021.txt",h=T)
  body.info <- subset(body.info, body.info$part=="body")
  body.raw <- read.table("AllDatasets/Rawcounts_CPM1_body_eqtCatalog_Aug3021.txt",h=T, check.names = F)
  
  #subsample 200, 500, 700, 939(all) individuals randomly
  bodyall <- body.raw
  set.seed(2345)
  body200 <- body.raw[, colnames(body.raw) %in% sample(body.info$id, 200, replace=F)]
  set.seed(2345)
  body400 <- body.raw[, colnames(body.raw) %in% 
                        union(sample(setdiff(body.info$id, colnames(body200)), 200, replace=F), colnames(body200))]
  set.seed(2345)
  body600 <- body.raw[, colnames(body.raw) %in% 
                        union(sample(setdiff(body.info$id, colnames(body400)), 200, replace=F), colnames(body400))]
  set.seed(2345)
  body800 <- body.raw[, colnames(body.raw) %in% 
                        union(sample(setdiff(body.info$id, colnames(body600)), 200, replace=F), colnames(body600))]
  
  #filter by average CPM>1 and >CPM1 in at least 20% flies, and just keep genes that overlap all datasets
  
  tmp = list(body200, body400, body600, body800, bodyall)
  tmp2 <- lapply(tmp, function(x) DGEList(counts=x))
  tmp2 <- lapply(tmp2, function(x) calcNormFactors(object = x))
  tmp2 <- lapply(tmp2, function(x) cpm(x, log = F) )
  
  tokeep <- lapply(tmp2, function(x) apply(x,1, function(y) (sum(y)/length(y))>1 & (sum(y>1)/length(y))>0.2))
  tokeep <- lapply(tokeep, function(x) subset(x, x==TRUE))
  
  #subset raw counts based on cpm filter
  test <- list()
  for (i in 1:5){
    test[[i]] <- subset(tmp[[i]], rownames(tmp[[i]]) %in% names(tokeep[[i]]))
  }
  
  #transform voom
  hvoom <- lapply(test, function(x) DGEList(counts=x))
  hvoom <- lapply(hvoom, function(x) calcNormFactors(object = x))
  hvoom <- lapply(hvoom, function(x) voom(x, design = NULL) )
  
  #pca to check for batches
  pcavoom <- lapply(hvoom, function(x) prcomp(t(x$E)))
  par(mfrow=c(2,3))   
  for (i in 1:5){
    plot(pcavoom[[i]]$x[,2] ~ pcavoom[[i]]$x[,1])
  }
  
  #remove plate effect
  voommatrix <- lapply(hvoom,function(x) as.matrix(x$E))
  #have to create a body.info file for each matrix
  bodyinfolist <- list()
  for (i in 1:5){
    bodyinfolist[[i]] <- subset(body.info, body.info$id %in% colnames(voommatrix[[i]]))
  }
  
  voomnoplate <- list()
  for (i in 1:5){
    voomnoplate[[i]] <- removeBatchEffect(voommatrix[[i]], batch=bodyinfolist[[i]]$plate)
  }
  
  voomnoplate.pca <- lapply(voomnoplate, function(x) prcomp(t(x)))
  
  par(mfrow=c(2,3))   
  for (i in 1:5){
    plot(voomnoplate.pca[[i]]$x[,2] ~ voomnoplate.pca[[i]]$x[,1])
  }
  
  pcaweight <- list()
  for (i in 1:5){
    for (j in 1:length(voomnoplate.pca[[i]]$sdev)) { 
      pcaweight[[i]][j]<-(voomnoplate.pca[[i]]$sdev[j])^2/sum(voomnoplate.pca[[i]]$sdev^2) } 
    plot(pcaweight[[i]][1:10])
  }
  
  
  #write voom tables and cov tables
  write.table(hvoom[[1]]$E, "AllDatasets/VOOMCounts_CPM1_body_ctrl_200ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F )
  write.table(hvoom[[2]]$E, "AllDatasets/VOOMCounts_CPM1_body_ctrl_400ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F )
  write.table(hvoom[[3]]$E, "AllDatasets/VOOMCounts_CPM1_body_ctrl_600ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F )
  write.table(hvoom[[4]]$E, "AllDatasets/VOOMCounts_CPM1_body_ctrl_800ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F )
  write.table(hvoom[[5]]$E, "AllDatasets/VOOMCounts_CPM1_body_ctrl_939ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F )
  
  #Print matrix of covariates for matrixeqtl 
  cov <- list()
  for (i in 1:5){
    cov[[i]] <- as.data.frame(t(bodyinfolist[[i]][,3]))
    colnames(cov[[i]]) <- bodyinfolist[[i]]$id
    rownames(cov[[i]]) <- "plate"
  }   
  
  write.table(cov[[1]], "AllDatasets/Cov_forMapping_body_ctrl_200ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F)
  write.table(cov[[2]], "AllDatasets/Cov_forMapping_body_ctrl_400ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F)
  write.table(cov[[3]], "AllDatasets/Cov_forMapping_body_ctrl_600ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F)
  write.table(cov[[4]], "AllDatasets/Cov_forMapping_body_ctrl_800ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F)
  write.table(cov[[5]], "AllDatasets/Cov_forMapping_body_ctrl_939ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F)
  
  #write covfree table that includes only samples used for mapping, 
  #939indv, 8651 genes
  write.table(voomnoplate[[5]], "AllDatasets/VOOMCounts_CPM1_body_ctrl_939ind_covfree_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F )
  
}

####  HEAD
{   
  #load raw counts of expressed genes
  head.info <- read.table("AllDatasets/Info_Samples_eqtlCatalog_Aug3021.txt",h=T)
  head.info <- subset(head.info, head.info$part=="head")
  head.raw <- read.table("AllDatasets/Rawcounts_CPM1_head_eqtCatalog_Aug3021.txt",h=T, check.names = F)
  
  #subsample 200, 500, 700, 940(all) individuals randomly
  headall <- head.raw
  set.seed(2345)
  head200 <- head.raw[, colnames(head.raw) %in% sample(head.info$id, 200, replace=F)]
  set.seed(2345)
  head400 <- head.raw[, colnames(head.raw) %in% 
                        union(sample(setdiff(head.info$id, colnames(head200)), 200, replace=F), colnames(head200))]
  set.seed(2345)
  head600 <- head.raw[, colnames(head.raw) %in% 
                        union(sample(setdiff(head.info$id, colnames(head400)), 200, replace=F), colnames(head400))]
  set.seed(2345)
  head800 <- head.raw[, colnames(head.raw) %in% 
                        union(sample(setdiff(head.info$id, colnames(head600)), 200, replace=F), colnames(head600))]
  
  #filter by CPM>1 again, and just keep genes that overlap all datasets
  
  tmp = list(head200, head400, head600, head800, headall)
  tmp2 <- lapply(tmp, function(x) DGEList(counts=x))
  tmp2 <- lapply(tmp2, function(x) calcNormFactors(object = x))
  tmp2 <- lapply(tmp2, function(x) cpm(x, log = F) )
  
  tokeep <- lapply(tmp2, function(x) apply(x,1, function(y) (sum(y)/length(y))>1 & (sum(y>1)/length(y))>0.2))
  tokeep <- lapply(tokeep, function(x) subset(x, x==TRUE))
  
  #subset raw counts based on cpm filter
  test <- list()
  for (i in 1:5){
    test[[i]] <- subset(tmp[[i]], rownames(tmp[[i]]) %in% names(tokeep[[i]]))
  }
  
  #transform voom
  hvoom <- lapply(test, function(x) DGEList(counts=x))
  hvoom <- lapply(hvoom, function(x) calcNormFactors(object = x))
  hvoom <- lapply(hvoom, function(x) voom(x, design = NULL) )
  
  #pca to check for batches
  pcavoom <- lapply(hvoom, function(x) prcomp(t(x$E)))
  par(mfrow=c(2,3))   
  for (i in 1:5){
    plot(pcavoom[[i]]$x[,2] ~ pcavoom[[i]]$x[,1])
  }
  
  #remove plate effect
  voommatrix <- lapply(hvoom,function(x) as.matrix(x$E))
  #have to create a head.info file for each matrix
  headinfolist <- list()
  for (i in 1:5){
    headinfolist[[i]] <- subset(head.info, head.info$id %in% colnames(voommatrix[[i]]))
  }
  
  voomnoplate <- list()
  for (i in 1:5){
    voomnoplate[[i]] <- removeBatchEffect(voommatrix[[i]], batch=headinfolist[[i]]$plate)
  }
  
  voomnoplate.pca <- lapply(voomnoplate, function(x) prcomp(t(x)))
  
  par(mfrow=c(2,3))   
  for (i in 1:5){
    plot(voomnoplate.pca[[i]]$x[,2] ~ voomnoplate.pca[[i]]$x[,1])
  }
  
  #estimate SVs
  #SVA
  mod <-list() ; mod0<-list() ; counts<-list(); svavoom <-list()
  for (i in 1:5){
    mod[[i]] <- model.matrix(~ plate, data = headinfolist[[i]]) # Plate info is known, i just want to identify other hidden batches
    mod0[[i]] <- model.matrix(~ 1,  data = headinfolist[[i]])
    counts[[i]] <- voommatrix[[i]]
    svavoom[[i]] <- sva(counts[[i]], mod[[i]], mod0[[i]], n.sv=10)
  }
  
  #remove svas + plate
  voomnoplatesva <- list();voomnoplate.pca<-list()
  
  for (i in 1:5){
    voomnoplatesva[[i]] <- removeBatchEffect(voommatrix[[i]], batch=headinfolist[[i]]$plate,
                                             covariates=svavoom[[i]]$sv[,c(1:3)])
    voomnoplate.pca[[i]] <- prcomp(t(voomnoplatesva[[i]]))
  }
  
  par(mfrow=c(2,3))   
  for (i in 1:5){
    plot(voomnoplate.pca[[i]]$x[,2] ~ voomnoplate.pca[[i]]$x[,1])
  }
  
  pcaweight <- vector("list", 5)
  for (i in 1:5){
    for (j in 1:length(voomnoplate.pca[[i]]$sdev)) { 
      pcaweight[[i]][j]<-(voomnoplate.pca[[i]]$sdev[j])^2/sum(voomnoplate.pca[[i]]$sdev^2)
    } 
    plot(pcaweight[[i]][1:10])
  }
  
  
  #write voom tables and cov tables
  write.table(hvoom[[1]]$E, "AllDatasets/VOOMCounts_CPM1_head_ctrl_200ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F )
  write.table(hvoom[[2]]$E, "AllDatasets/VOOMCounts_CPM1_head_ctrl_400ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F )
  write.table(hvoom[[3]]$E, "AllDatasets/VOOMCounts_CPM1_head_ctrl_600ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F )
  write.table(hvoom[[4]]$E, "AllDatasets/VOOMCounts_CPM1_head_ctrl_800ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F )
  write.table(hvoom[[5]]$E, "AllDatasets/VOOMCounts_CPM1_head_ctrl_940ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F )
  
  #Print matrix of covariates for matrixeqtl , only 3 svas are necesary
  cov <- list()
  for (i in 1:5){
    cov[[i]] <- as.data.frame(t(cbind.data.frame(as.factor(headinfolist[[i]]$plate), svavoom[[i]]$sv[,1], svavoom[[i]]$sv[,2], svavoom[[i]]$sv[,3])))
    colnames(cov[[i]]) <- headinfolist[[i]]$id
    rownames(cov[[i]]) <- c("plate", "sv1","sv2", "sv3")
  }   
  
  
  write.table(cov[[1]], "AllDatasets/Cov_forMapping_head_ctrl_200ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F)
  write.table(cov[[2]], "AllDatasets/Cov_forMapping_head_ctrl_400ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F)
  write.table(cov[[3]], "AllDatasets/Cov_forMapping_head_ctrl_600ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F)
  write.table(cov[[4]], "AllDatasets/Cov_forMapping_head_ctrl_800ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F)
  write.table(cov[[5]], "AllDatasets/Cov_forMapping_head_ctrl_940ind_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F)
  
  #write covfree table 
  #940indv, 8877 genes
  write.table(voomnoplatesva[[5]], "AllDatasets/VOOMCounts_CPM1_head_ctrl_940ind_covfree_Aug3121.txt", row.names = T, col.names = T, sep="\t", quote=F )
  
  
}

#### GEMMA 
{
  #to be run in the server for each subset
  #example of code is available in this folder to:
  #estimate GRM and make fam files for the subsets generated in this R script
  #run GEMMA and extract subset of pvalues for FDR correction in R
  
  #all files needed to run the full data set are in Edmond https://doi.org/10.17617/3.8JLETW
 
  
  
}


