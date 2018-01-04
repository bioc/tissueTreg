# all testing

# loading
library(ExperimentHub)
readMetadataFromCsv("/home/imbusch/git/tissueTreg/")
# local testing

# using SummarizedExperiment but doesn't keep the smoothed values
# library(SummarizedExperiment)
# mytest <- loadHDF5SummarizedExperiment(dir = "/icgc/dkfzlsdf/analysis/G200/imbusch/treg/treg_experimentHub/SummarizedExperiment/")
# pheatmap::pheatmap(getCoverage(mytest[1000000:1000100,]))

# directly using bsseq values
BS.obj.ex.fit <- readRDS(file = "/icgc/dkfzlsdf/analysis/G200/imbusch/treg/treg_experimentHub/treg_twgbs.rds")
# example plot with FOXP3 using the bsseq object and smoothed values
# add annotation and colors ...
regions <- GRanges(
  seqnames = c("X"),
  ranges = IRanges(start = c(7579676),
                   end = c(7595243)
  )
)
plotRegion(BS.obj.ex.fit, region=regions[1,], extend = 2000)
