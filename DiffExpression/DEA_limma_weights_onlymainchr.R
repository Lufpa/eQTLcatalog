#activate conda env that has limma installed
#conda activate limma
#check that R version is 4.

#only samples used in mapping will be used here (head = 940, body=939)
#only genes in main chromosomes will be used here (n= 7800 shared between head and body)

library(limma)
library(edgeR) 

start=Sys.time()
Sys.time()

cat("loading data \n")
countdata <- read.table("RawCounts_CPM1_headbody_Sep1621.txt",h=T, check.names = F)
coldata <- read.table("Covariates_forDEA_bodyhead_Sep1621.txt",h=T)
coldata$well <- as.factor(coldata$well)
coldata$plate <- as.factor(coldata$plate)
coldata$treatment <- as.factor(coldata$treatment)
coldata$RNAlibBatch <- as.factor(coldata$RNAlibBatch)
coldata$RNAseqBatch <- as.factor(coldata$RNAseqBatch)
coldata$egglayBatch <- as.factor(coldata$egglayBatch)
coldata$eclosionBatch <- as.factor(coldata$eclosionBatch)
coldata$platingBatch <- as.factor(coldata$platingBatch)
coldata$part <- as.factor(coldata$part)
cat ("data loaded \n")

design <- model.matrix(~0+coldata$part+
                         coldata$sv1+coldata$sv2+coldata$sv3+coldata$sv4+
                         coldata$plate) #body=0, head=1 (resutls are head-body)
colnames(design) <- c("body","head", "sv1", "sv2", "sv3", "sv4", "plate25", "plate31", "plate82", "plate86", "plate88","plate91","plate106","plate108","plate109","plate111","plate112","plate113","plate115","plate137", "plate161","plate206")
contrast <- makeContrasts(HvsB=head-body, levels = design)

tmp1<-Sys.time()
cat("normalizing and estimating mean-variance weights \n")
countdata.list <- DGEList(counts=countdata)
countdata.norm <- calcNormFactors(countdata.list) #gets norm. factors based on TMM (controlling not only for lib size, but also composition)
countdata.voom <- voom(countdata.norm, design = design, plot=F) #the E matrix of normalized counts is the same as in boomWithQualityWeights

cat("estimating correlation between same individuals \n")
corfit <- duplicateCorrelation(countdata.voom,design,block=coldata$id)
tmp2<-Sys.time()
cat("average correlation = ", corfit$consensus, "\n")
cat("est. corr = ", tmp2-tmp1, "\n")

cat("running analysis \n")
fit <- lmFit(countdata.voom,design,block=coldata$id,correlation=corfit$consensus)
fitc <- contrasts.fit(fit, contrasts = contrast)
fit2 <- eBayes(fitc)
summary(decideTests(fit2))
res<- topTable(fit2, number = Inf)
tmp3<-Sys.time()
cat("analysis done", tmp3-tmp2, "\n")

cat("writing results table and RDS object with all info from fit2 object \n")
write.table(res, "DiffExpression_results_headbody_Sep1621.txt", col.names = T, row.names = T, quote = F, sep="\t")
saveRDS(fit2, "DiffExpression_results_headbody_Sep1621.rds")
saveRDS(fit, "DiffExpression_lmFit_headbody_Sep1621.rds") #from here I get the mean expression logcpm per group fit$coefficients
#and the stderror fit$stdev.unscaled * fit$sigma

Sys.time()
end=Sys.time()
cat("total time = ", end-start)


