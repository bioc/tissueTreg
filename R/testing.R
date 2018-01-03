# all testing

# loading
readMetadataFromCsv("/home/imbusch/git/tissueTreg/")
# local testing

# using SummarizedExperiment but doesn't keep the smoothed values
# library(SummarizedExperiment)
# mytest <- loadHDF5SummarizedExperiment(dir = "/icgc/dkfzlsdf/analysis/G200/imbusch/treg/treg_experimentHub/SummarizedExperiment/")
# pheatmap::pheatmap(getCoverage(mytest[1000000:1000100,]))

# directly using bsseq values
mytest <- readRDS(file = "/icgc/dkfzlsdf/analysis/G200/imbusch/treg/treg_experimentHub/treg_twgbs.rds")
pheatmap::pheatmap(getMeth(mytest)[11000:11100,], cluster_rows = F, cluster_cols = F)
