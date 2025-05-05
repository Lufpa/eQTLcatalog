Example of code to run differential expression analyses with partially overlapping samples (e.g. head and body RNAseq data comes from the same individual). 

- DiffExpression.R Describes how the matrices for this analyses were processed and how covariates including Surrogate Variables were estimated

- DEA_limma_weights.R Describes how to run the DEA in limma

Data needed to use this code:
- Rawcounts_CPM1_body_eqtCatalog_Aug3021.txt (The raw gene count tables are available in Edmond: https://edmond.mpg.de/dataset.xhtml?persistentId=doi:10.17617/3.8JLETW)
- Rawcounts_CPM1_head_eqtCatalog_Aug3021.txt (The raw gene count tables are available in Edmond: https://edmond.mpg.de/dataset.xhtml?persistentId=doi:10.17617/3.8JLETW)
- Info_Samples_eqtlCatalog_Aug3021.txt (This table is deposited here)