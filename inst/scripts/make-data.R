library(bsseq)
library(data.table)

# creates a bsseq object and stores it as a rds file
# files are read from a custom in-house bed file with CpG context only and read counts from
# + and - strand have prior to this point already been merged

chromosomes <- c(seq(1,19), "X", "Y")
replicates <- c("R1", "R2", "R3")
tissue_cells <- c("Fat-Treg", "Liver-Treg", "Skin-Treg", "Lymph-N-Tcon", "Lymph-N-Treg")
# two lists for coverage and methylation info
list_cov <- list()
list_ureads <- list()
# keep track of chromosomes and positions seperately
chrs <- NULL
positions <- NULL
created_flag <- FALSE

# load data
for (tissue_cell in tissue_cells) {
  for (replicate in replicates) {
    dm_per_replicate_cov <- NULL
    dm_per_replicate_ureads <- NULL
    for (chromosome in chromosomes) {
      filename <- paste0("/icgc/dkfzlsdf/analysis/G200/imbusch/treg/treg_experimentHub/bed_files/", tissue_cell, "-", replicate, "-chr", chromosome, ".bed")
      message(paste("Loading", filename))
      dm <- fread(input = filename, data.table=F, col.names = c("chr", "start", "stop", "strand", "context", "ureads", "creads"))
      # coverage
      dm_cov <- as.numeric(dm$ureads) + as.numeric(dm$creads)
      dm_per_replicate_cov <- c(dm_per_replicate_cov, dm_cov)
      # ureads
      dm_ureds <- as.numeric(dm$ureads)
      dm_per_replicate_ureads <- c(dm_per_replicate_ureads, dm_ureds)
      # save chromosome and position as well
      if( !created_flag ) {
        chrs <- c(chrs, dm$chr)
        positions <- c(positions, dm$start)
      }
    }
    # set flag
    created_flag <- TRUE
    # sample name
    sample_name <- paste(tissue_cell, replicate, sep="-")
    # save
    list_cov[[sample_name]] <- dm_per_replicate_cov
    list_ureads[[sample_name]] <- dm_per_replicate_ureads
    # size check
    sl <- object.size(list_cov)
    print(sl, units = "auto", standard = "SI")
  }
}
# combine them
dm_big_cov <- do.call(cbind, list_cov)
dm_big_ureads <- do.call(cbind, list_ureads)
# sanity check
stopifnot(all(colnames(dm_big_cov) == colnames(dm_big_ureads)))
### create bsseq object and smoothing
BS.obj.ex <- BSseq(M=dm_big_ureads, Cov=dm_big_cov, chr=chrs, pos=positions)
BS.obj.ex.fit <- BSmooth(BS.obj.ex, mc.cores=24, verbose=TRUE, parallelBy = "chr")
# TODO: remove this step? saveHDF5SummarizedExperiment does not seem to preserve the smoothing part from bsseq
BS.obj.ex_hdf5 <- saveHDF5SummarizedExperiment(x = BS.obj.ex, dir = "/icgc/dkfzlsdf/analysis/G200/imbusch/treg/treg_experimentHub/SummarizedExperiment/", replace=TRUE)
# save as RDS file
saveRDS(object = BS.obj.ex.fit, file = "/icgc/dkfzlsdf/analysis/G200/imbusch/treg/treg_experimentHub/treg_twgbs.rds")
