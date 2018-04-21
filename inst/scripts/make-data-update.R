# convert data to new bsseq object which has a newer DelayedArray version

library(ExperimentHub)
library(bsseq)

eh <- ExperimentHub()
BS.obj.ex.fit <- eh[["EH1072"]]
BS.obj.ex.fit2 <- BSseq(M = assays(BS.obj.ex.fit)[["M"]], Cov = assays(BS.obj.ex.fit)[["Cov"]], gr = rowRanges(BS.obj.ex.fit))
BS.obj.ex.fit2_smoothed <- BSmooth(BS.obj.ex.fit2, mc.cores=8, verbose=TRUE, parallelBy = "chr")
saveRDS(object = BS.obj.ex.fit2_smoothed, file = "treg_twgbs_per_sample.rds")

BS.obj.ex.fit.combined <- eh[["EH1073"]]
BS.obj.ex.fit.combined2 <- BSseq(M = assays(BS.obj.ex.fit.combined)[["M"]], Cov = assays(BS.obj.ex.fit.combined)[["Cov"]], gr = rowRanges(BS.obj.ex.fit.combined))
BS.obj.ex.fit.combined_smoothed <- BSmooth(BS.obj.ex.fit.combined2, mc.cores=8, verbose=TRUE, parallelBy = "chr")
saveRDS(object = BS.obj.ex.fit.combined_smoothed, file = "treg_twgbs_per_group.rds")
