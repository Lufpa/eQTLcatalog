#########
# DIFF EXPRESSION
#########
# generate joint matrix of body/head shared genes and covs
# 

#get genes overlapping between both tissues
#using raw data because the transformations are goona be done directly in the
#diff expression analysis modelling
{     
  bodyinfo <- read.table("AllDatasets/Info_Samples_eqtlCatalog_Aug3021.txt",h=T)
  bodyinfo <- subset(body.info, part=="body")
  body.raw <- read.table("AllDatasets/Rawcounts_CPM1_body_eqtCatalog_Aug3021.txt",h=T, check.names = F)
  
  head.raw <- read.table("AllDatasets/Rawcounts_CPM1_head_eqtCatalog_Aug3021.txt",h=T, check.names = F)
  headinfo <- read.table("AllDatasets/Info_Samples_eqtlCatalog_Aug3021.txt",h=T)
  headinfo <- subset(head.info, head.info$part=="head")
  
  #Using all the QCed RNAseq data, there are 593 paired samples
  #347 unique head samples
  #346 unique body samples
  commonind <- intersect(bodyinfo$id, headinfo$id) ; length(commonind)
  uniqind.h <- setdiff(headinfo$id, bodyinfo$id) ; length(uniqind.h)
  uniqind.b <- setdiff(bodyinfo$id, headinfo$id) ; length(uniqind.b)
  
  #Recode the batches to differentiate body from head 
  bodyinfo$RNAseqBatch[bodyinfo$RNAseqBatch==1] <- 6 #to differentiate from head batches
  bodyinfo$RNAseqBatch[bodyinfo$RNAseqBatch==2] <- 3 #to match head batch3 that this body2 batch is part of
  
  bodyinfo$RNAlibBatch[bodyinfo$RNAlibBatch==1] <- 4 #to differentiate from head batches
  bodyinfo$RNAlibBatch[bodyinfo$RNAlibBatch==2] <- 5 #to differentiate from head batches
  
  #identify head and body samples in the counts matrix
  colnames(body.raw) <- paste(colnames(body.raw),"_b",sep="")
  colnames(head.raw) <- paste(colnames(head.raw),"_h", sep="")
  
  #7800 overlapping genes, 1079 head unique,579
  commongenes <- intersect(rownames(head.raw), rownames(body.raw)) ; length(commongenes)
  length(setdiff(rownames(head.raw), rownames(body.raw)))
  length(setdiff(rownames(body.raw), rownames(head.raw)))
  
  bothcounts <- cbind(body.raw[rownames(body.raw) %in% commongenes,], head.raw[rownames(head.raw) %in% commongenes,])
  write.table(bothcounts, "AllDatasets/RawCounts_CPM1_headbody_Sep1621.txt",
              col.names = T, row.names = T, quote=F, sep="\t")
  
}  

#Get covariates for diff. expression analysis
#this will be estimated in the joint table

{
  #voom transform the raw counts since that matrix is what is gonna be used in DE analysis
  countdata <- read.table("AllDatasets/RawCounts_CPM1_heabody_Sep1621.txt",h=T,check.names = F)
  normcounts <- DGEList(counts=countdata)
  normcounts <- calcNormFactors(normcounts) #gets norm. factors based on TMM (controlling not only for lib size, but also composition)
  normcounts <- voom(normcounts, design = NULL, plot=F) #the E matrix of normalized counts is the same as in boomWithQualityWeights
  # log2(counts+0.5/normalized lib.size in Millions), norm library size = totalreads*norm.factors
  normcounts <- normcounts$E
  
  coldata <- read.table("AllDatasets/Info_Samples_eqtlCatalog_Aug3021.txt",h=T)
  coldata$well <- as.factor(coldata$well)
  coldata$plate <- as.factor(coldata$plate)
  coldata$treatment <- as.factor(coldata$treatment)
  coldata$RNAlibBatch <- as.factor(coldata$RNAlibBatch)
  coldata$RNAseqBatch <- as.factor(coldata$RNAseqBatch)
  coldata$egglayBatch <- as.factor(coldata$egglayBatch)
  coldata$eclosionBatch <- as.factor(coldata$eclosionBatch)
  coldata$platingBatch <- as.factor(coldata$platingBatch)
  coldata$part <- as.factor(coldata$part)
  
  #sort coldata to match countmatrix
  coldata <- coldata[order(match(rownames(coldata), colnames(normcounts))),]
  all.equal(rownames(coldata), colnames(normcounts))
  
  #to easily plot pcas w/o Deseq
  #install.packages("ggfortify")
  library("ggfortify")
  #install.packages("vctrs")
  library(vctrs)
  pca<-prcomp(t(normcounts))
  autoplot(pca, data=coldata, colour="egglayBatch",x=1, y=2)
  pcaw <- NULL ; for (i in 1:length(pca$sdev)) { pcaw[i]= (pca$sdev[i])^2/sum(pca$sdev^2) } 
  plot(pcaw); plot(pcaw,xlim=c(0,10)); sum(pcaw[1:100]) #first 100 pcs account for 66.6% of var
  
  #exploring effect of covs
  table(coldata$part, coldata$plate) #body and head have same plate id, so 
  #i could again just use palte to accoutn for egglay, eclossion, plating date
  #even for processing
  
  table(coldata$part, coldata$RNAlibBatch) #totally confunded by bodypart - can't use
  table(coldata$part, coldata$RNAseqBatch) #6batches, only one has both bodyparts
  table(coldata$part, coldata$egglayBatch) #plate = egglay = eclosion = plating
  table(coldata$part, coldata$eclosionBatch) #ojo plate 108 has two eclosion batches 2 and 5, this is a mistake
  table(coldata$plate, coldata$platingBatch) #plate 108 eclosion 22 and 23 days. mistake. eclosion ~= plating batch use only one
  
  
  car::Anova(lm(pca$x[,1:1000]~ coldata$part + coldata$RNAseqBatch + coldata$egglayBatch + coldata$platingBatch))
  car::Anova(lm(pca$x[,1:1000]~ coldata$part + coldata$egglayBatch + coldata$platingBatch))
  car::Anova(lm(pca$x[,1:100]~ coldata$part + coldata$plate))
  
  #how does the pca look like after removing all known covs
  nullmodel <- model.matrix(~coldata$part)
  counts.new <- limma::removeBatchEffect(normcounts, batch = coldata$plate,  design = nullmodel,) 
  
  pca1<-prcomp(t(counts.new))
  autoplot(pca1, data=coldata, colour="plate",x=1, y=2)
  autoplot(pca1, data=coldata, colour="plate",x=3, y=4)
  autoplot(pca1, data=coldata, colour="plate",x=5, y=6)
  pca1w <- NULL ; for (i in 1:length(pca1$sdev)) { pca1w[i]= (pca1$sdev[i])^2/sum(pca1$sdev^2) } 
  plot(pca1w); plot(pca1w[1:10]); sum(pca1w[1:100]) #first 100 pcs account for 65.4% of var
  
  car::Anova(lm(pca1$x[,1:1000]~ coldata$part + coldata$plate))
  
  #do I need to add SVAs? check if they identify clear batches, or make disapear the head pattern.
  library(sva)
  
  mod  <- model.matrix(~ plate+part, data = coldata) #full model, known batches + variable of interest
  mod0 <- model.matrix(~ plate, data = coldata) #null model, only known batches
  
  # svaseq is recommended for RNAseq data, while sva for microarray-like data (symmetrically distributed)
  # but, voom-transformed data makes RNAseq count data behave normally distributed, so I'll use sva https://support.bioconductor.org/p/54723/
  
  normcountsm <- as.matrix(normcounts) #needs counts to be in a matrix
  sva.sva = sva(normcountsm, mod, mod0, n.sv = 10)
  par(mfrow=c(2,2))
  plot(sva.sva$sv[,1], sva.sva$sv[,2], col=coldata$part, pch=16) 
  plot(sva.sva$sv[,3], sva.sva$sv[,4], col=coldata$part, pch=16)
  plot(sva.sva$sv[,5], sva.sva$sv[,6], col=coldata$part, pch=16); 
  plot(sva.sva$sv[,1], pca1$x[,1]) #it seems that sv1 and sv2 might be capturing the weird head pattern
  plot(sva.sva$sv[,2], pca1$x[,2])
  plot(sva.sva$sv[,3], pca1$x[,1])
  
  #how important are those sVs in explaining gene expression variation
  
  car::Anova(lm(pca$x[,1:1000]~ coldata$part +coldata$plate +
                  sva.sva$sv[,1] + sva.sva$sv[,2] + sva.sva$sv[,3] + sva.sva$sv[,4]+ sva.sva$sv[,5]+ sva.sva$sv[,6] ))
  
  sv1<- sva.sva$sv[,1]; sv2<- sva.sva$sv[,2]; sv3<- sva.sva$sv[,3]
  sv4 <- sva.sva$sv[,4] ; sv5 <- sva.sva$sv[,5] ; sv6 <- sva.sva$sv[,6]
  
  #how does the PCA look after adding SVAs to the known batches
  plate<-coldata$plate
  contrasts(plate) <- contr.sum(levels(plate))
  
  covariates <- model.matrix(~plate+sv1)
  counts.new2 <- limma::removeBatchEffect(normcounts, covariates=covariates[,-1], design = nullmodel) 
  covariates <- model.matrix(~plate+sv1+sv2)
  counts.new3 <- limma::removeBatchEffect(normcounts, covariates=covariates[,-1], design = nullmodel) 
  covariates <- model.matrix(~plate+sv1+sv2+sv3)
  counts.new4 <- limma::removeBatchEffect(normcounts, covariates=covariates[,-1], design = nullmodel) 
  covariates <- model.matrix(~plate+sv1+sv2+sv3+sv4)
  counts.new5 <- limma::removeBatchEffect(normcounts, covariates=covariates[,-1], design = nullmodel) 
  
  pca2 <- prcomp(t(counts.new2))
  pca3 <- prcomp(t(counts.new3))
  pca4 <- prcomp(t(counts.new4))
  pca5 <- prcomp(t(counts.new5))
  
  
  par(mfrow=c(3,2))
  autoplot(pca2, data=coldata, colour="part",x=1, y=2);autoplot(pca2, data=coldata, colour="part",x=3, y=4)
  autoplot(pca3, data=coldata, colour="part",x=1, y=2);autoplot(pca3, data=coldata, colour="part",x=3, y=4)
  autoplot(pca4, data=coldata, colour="part",x=1, y=2);autoplot(pca4, data=coldata, colour="part",x=3, y=4)
  autoplot(pca5, data=coldata, colour="part",x=1, y=2);autoplot(pca5, data=coldata, colour="part",x=3, y=4)
  
  #write svs into the covariate file
  coldata2 <- cbind(coldata, sv1, sv2, sv3, sv4)
  write.table(coldata2, "AllDatasets/Covariates_forDEA_bodyhead_Sep1621.txt",col.names = T, row.names = T, sep="\t", quote=F)
  
  #write covfree voom table
  write.table(counts.new5, "AllDatasets/VOOMCounts_CPM1_bodyhead_covfree_Sep1621.txt",
              col.names = T, row.names = T, quote=F, sep="\t")
  
  #write the voom table raw
  write.table(normcounts, "AllDatasets/VOOMCounts_CPM1_bodyhead_Sep1621.txt",
              col.names = T, row.names = T, quote=F, sep="\t")
  
  #make sure all printed tables are in the right order
  all.equal(colnames(counts.new5), rownames(coldata2))
  all.equal(colnames(normcounts), rownames(coldata2))
  
} 


##run differential expression analyses in the server
{
  #code available in this folder: DEA_limma_weights_onlymainchr.R  

  
  
}


