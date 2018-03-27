library(testthat)
library(tissueTreg)

# init
library(ExperimentHub)
library(bsseq)
eh <- ExperimentHub()
### check on object classes in the package
# dna methylation
expect_is(eh[["EH1072"]], "BSseq")
expect_is(eh[["EH1073"]], "BSseq")
# # rna-seq
expect_is(eh[["EH1074"]], "SummarizedExperiment")
expect_is(eh[["EH1075"]], "SummarizedExperiment")
# check on number of samples in each object
expect_equal(dim(eh[["EH1072"]])[2], 15)
expect_equal(dim(eh[["EH1073"]])[2], 5)
expect_equal(dim(eh[["EH1074"]])[2], 15)
expect_equal(dim(eh[["EH1075"]])[2], 15)
