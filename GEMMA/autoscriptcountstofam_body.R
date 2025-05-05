#Script to make 2) counts into fam and 2) cov into cov file for GEMMA
### BODY (path to files is to body folders)


###
# Counts matrix into fam file
###

workingdir <- "/Genomics/ayroleslab2/lamaya/bigProject/Genotypes_feb2020/headbody_notsameindv/ForGemma"

files <- c("200", "400", "600", "800", "939")

	for (j in 1:5){

	i <- files[j]
	cat("subset being processed: ",i)

	#subset binary files
	#get samples to keep
	cov <- read.table(paste("/Genomics/ayroleslab2/lamaya/bigProject/eQTLcatalog/body/GEMMA/subsets/Cov_forMapping_body_ctrl_", i, "ind_Aug3121.txt", sep=""), check.names=F)
	cov.t <- t(cov)
	write.table(cbind(cov.t[,1], rownames(cov.t)), "tokeep", col.names=F, row.names=F, quote=F)
	system(paste("plink --bfile ", workingdir, "/ctrl_headbody.final.bodyonlynorelatedness1 --keep ", "tokeep --make-bed --out ", workingdir, "/ctrl_headbody.final.bodyonlynorelatedness1.", i, "ind.Sep221", sep="")) 

	#load gene counts matrix
	#each matrix is unique, not a subsample of the entire ~950indv matrix because the VOOM
	#normalization changes depending on the individuals included, so the normalization
	#was runned for each subset separately

	counts <- read.table(paste("/Genomics/ayroleslab2/lamaya/bigProject/eQTLcatalog/body/GEMMA/subsets/VOOMCounts_CPM1_body_ctrl_",i, "ind_Aug3121.txt", sep=""),h=T, row.names=1, check.names=F)

	#counts to fam

	fam <- read.table(paste(workingdir, "/ctrl_headbody.final.bodyonlynorelatedness1.", i, "ind.Sep221.fam", sep=""))
	counts.t <- t(counts)
	test <- as.data.frame(cbind("V2"= rownames(counts.t), counts.t))
	newfam <- merge(fam[,-6], test, by="V2") #add gene expression info

	 if(isTRUE(all.equal(newfam$V2, fam$V2))){
      		cat("all samples are in the right order\n")
	    } else {
		cat("ERROR: something went wrong with the sample order in fam and new fam\n")
	    }

	newfam2 <- newfam[ , c(2,1,3:ncol(newfam))]
	geneorder <- colnames(newfam2)
	geneorder <- geneorder[-1:-5]

	write.table(geneorder, paste(workingdir,"/geneorderFAMfile_body_",i,"ind.Sep221.txt", sep=""), quote=F, col.names=F, row.names=F)
	#remove fam subseted above in plink, it contained gene counts from the original fam file.
	system(paste("rm ", workingdir, "/ctrl_headbody.final.bodyonlynorelatedness1.",i,"ind.Sep221.fam",sep=""))

	#write new fam with same name as removed before
	write.table(newfam2, paste(workingdir,"/ctrl_headbody.final.bodyonlynorelatedness1.",i,"ind.Sep221.fam",sep=""), quote=F, col.names=F, row.names=F)


	#cov into cov GEMMA file
	#load covariate file #done above

	cov.t <- cbind(1, cov.t)
	cov2 <- cov.t[order(match(rownames(cov.t), newfam2$V2)),]

	 if(isTRUE(all.equal(as.character(newfam$V2), rownames(cov2)))){
	      cat("all samples are in the right order\n")
	    } else {
	      cat("ERROR: something went wrong with the sample order in cov file\n")
	    }

	write.table(cov2, paste(workingdir,"/ctrl_headbody.final.body.covs.forGemma.",i,"ind.Sep221", sep=""), col.names=F, row.names=F, sep="\t", quote=F)

cat("done, files created for subset", i, "path: /Genomics/ayroleslab2/lamaya/bigProject/Genotypes_feb2020/headbody_notsameindv/ForGemma/ subset.fam, subset.bed, subset.bim, subset.cov, geneorder")

}

