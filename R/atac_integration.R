library(ExperimentHub)
library(bsseq)
library(data.table)
library(ggplot2)
library(ComplexHeatmap)

eh <- ExperimentHub()
# methylation
BS.obj.ex.fit.combined <- eh[["EH1073"]]
# rna-seq
se_rpkms <- eh[["EH1074"]]

regions <- GRanges(
  seqnames = c("X"),
  ranges = IRanges(start = c(7579676),
                   end = c(7595243)
  )
)

# helper function to convert ids
resolveGenesymbolToEntrez <- function(genesysmbol = NA) {
  # browser()
  library(biomaRt)
  mart <- useMart("ensembl")
  mart <- useDataset("mmusculus_gene_ensembl", mart = mart)
  genes.table <- getBM(filters= "ensembl_gene_id",
                       attributes= c("ensembl_gene_id", "entrezgene",
                                     "external_gene_name"), values=genesysmbol, mart=mart)
  genesysmbol <- data.frame(ensembl_gene_id = genesysmbol)
  results <- merge(x=as.data.frame(genesysmbol), y=genes.table, by.x="ensembl_gene_id", by.y="ensembl_gene_id", all.x = T)
  # cut results further, only return unique results
  results <- results[isUnique(results$external_gene_name),]
  return(results)
}

pData <- pData(BS.obj.ex.fit.combined)
pData$col <- c("red", "blue", "yellow", "orange", "green")
pData(BS.obj.ex.fit.combined) <- pData
plotRegion(BS.obj.ex.fit.combined, region=atac_regions[3,], extend = 2000)
### ATAC seq data: periphery tregs
tissetregs_periphery_peaks.coverage_qnorm <- read.table(file = "/Users/charlie/Dropbox/tissuetregs_Charles/tissetregs_periphery_peaks.coverage_qnorm.annotated.csv", sep = ",", header = TRUE)
periphery_differential.deseq_result <- read.table(file = "/Users/charlie/Dropbox/tissuetregs_Charles/periphery_differential.deseq_result.all_comparisons.csv", sep = ",", header = TRUE)
#
samples_of_interest <- c("Fat_Treg_1", "Fat_Treg_2", "Fat_Treg_3")
# TODO: select Peaks differently of course
atac_seq <- subset(tissetregs_periphery_peaks.coverage_qnorm, distance <= 5000)
# construct the actual genomic regions
atac_regions <- GRanges(
  seqnames = gsub("chr", "", atac_seq$chrom),
  ranges = IRanges(start = atac_seq$start,
                   end = atac_seq$end
  )
)
# methylation heatmap
meth_subset <- as.matrix(getMeth(BS.obj.ex.fit.combined, regions = atac_regions, type = "smooth", what = "perRegion"))
meth_heatmap <- Heatmap(meth_subset[complete.cases(meth_subset),],
                        show_column_names = T,
                        show_row_names = F,
                        cluster_columns = F
)
# atac
# remove regions where there is no data from methylation assay
atac_heatmap <- Heatmap(atac_seq[complete.cases(meth_subset),samples_of_interest],
                        show_column_names = T,
                        show_row_names = F,
                        cluster_columns = F
)
# convert ids from RNA-seq to make compatible
rna_seq_mapping <- resolveGenesymbolToEntrez(rownames(assay(se_rpkms)))
# get expression of genes
genes_external_gene_name <- as.character(atac_seq[complete.cases(meth_subset),]$gene_name)
genes_ensembl_gene_id <- rna_seq_mapping[genes_external_gene_name %in% rna_seq_mapping$external_gene_name,]$ensembl_gene_id
# TODO: finish, integrate somehow

# integration ATAC + Meth

atac_meth_integrated_fat <- data.frame(atac = rowMeans(atac_seq[complete.cases(meth_subset),samples_of_interest]),
                                   meth = meth_subset[complete.cases(meth_subset),c("Fat-Treg")])
plot(atac_meth_integrated_fat, xlim=c(-6,+10))

atac_meth_integrated_skin <- data.frame(atac = rowMeans(atac_seq[complete.cases(meth_subset), grepl("Skin", colnames(atac_seq))]),
                                       meth = meth_subset[complete.cases(meth_subset),c("Skin-Treg")])
plot(atac_meth_integrated_skin, xlim=c(-6,+10))


# plot heatmap
# atac_heatmap + meth_heatmap

# test for atac seq
mean_tissue <- data.frame(
  Fat = rowMeans(tissetregs_periphery_peaks.coverage_qnorm[grepl("Fat",colnames(tissetregs_periphery_peaks.coverage_qnorm))]),
  Lung = rowMeans(tissetregs_periphery_peaks.coverage_qnorm[grepl("Lung",colnames(tissetregs_periphery_peaks.coverage_qnorm))]),
  Colon = rowMeans(tissetregs_periphery_peaks.coverage_qnorm[grepl("Colon",colnames(tissetregs_periphery_peaks.coverage_qnorm))]),
  Skin = rowMeans(tissetregs_periphery_peaks.coverage_qnorm[grepl("Skin",colnames(tissetregs_periphery_peaks.coverage_qnorm))])
)
mean_tissue_melted <- melt(mean_tissue)
ggplot(mean_tissue_melted, aes(x=value, group=variable, col=variable)) +
  geom_freqpoly(binwidth = 0.05)
# get methylation values just to play with first
atac_seq <- subset(tissetregs_periphery_peaks.coverage_qnorm)
# construct the actual genomic regions
atac_regions <- GRanges(
  seqnames = gsub("chr", "", atac_seq$chrom),
  ranges = IRanges(start = atac_seq$start,
                   end = atac_seq$end
  )
)
atac_regions_meth <- getMeth(BS.obj.ex.fit.combined, regions = atac_regions, type = "smooth", what = "perRegion")
atac_regions_meth_melted <- melt(as.matrix(atac_regions_meth))
# plot "histogram" (frequency polygons) of beta values for each sample 
ggplot(atac_regions_meth_melted, aes(x=value, group=Var2, col=Var2)) +
  geom_freqpoly(binwidth = 0.01)
