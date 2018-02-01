# loading
library(ExperimentHub)
readMetadataFromCsv("/home/imbusch/git/tissueTreg/")

### local testing of methylation data
# directly using bsseq values
BS.obj.ex.fit <- readRDS(file = "/icgc/dkfzlsdf/analysis/G200/imbusch/treg/treg_experimentHub/treg_twgbs.rds")
BS.obj.ex.fit <- readRDS(file = "/icgc/dkfzlsdf/analysis/G200/imbusch/treg/treg_experimentHub/treg_twgbs_per_group.rds")
# example plot with FOXP3 using the bsseq object and smoothed values
# add annotation and colors ...
regions <- GRanges(
  seqnames = c("X"),
  ranges = IRanges(start = c(7579676),
                   end = c(7595243)
  )
)
plotRegion(BS.obj.ex.fit, region=regions[1,], extend = 2000)

### local testing of RNA-seq data
se_rpkms <- readRDS(file = "/icgc/dkfzlsdf/analysis/G200/imbusch/treg/treg_experimentHub/RPKM_table.rds")
