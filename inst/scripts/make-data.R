library(bsseq)
library(data.table)

# creates a bsseq object and stores it as a rds file
# files are read from a custom in-house bed file with CpG context only and read counts from
# + and - strand have prior to this point already been merged

maindir <- "/icgc/dkfzlsdf/analysis/G200/imbusch/treg/treg_experimentHub/"
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
# save as RDS file
saveRDS(object = BS.obj.ex.fit, file = paste0(maindir, "treg_twgbs_per_sample.rds"))
# collapse the samples into tissue groups, perform smoothing and save
group_sample_labels <- rep(tissue_cells, rep(3, 5))
names(group_sample_labels) <- paste(rep(tissue_cells, rep(3, 5)) , replicates, sep="-")
BS.obj.ex.collapsed <- collapseBSseq(BS.obj.ex, columns = group_sample_labels)
BS.obj.ex.collapsed.fit <- BSmooth(BS.obj.ex, mc.cores=24, verbose=TRUE, parallelBy = "chr")
saveRDS(object = BS.obj.ex.collapsed.fit, file = paste0(maindir, "treg_twgbs_per_group.rds"))

### compiles RNA-seq data

### process rpkms
# load and prepare the data
filePath <- "/home/imbusch/git/treg_experimentHub/RNA-seq/D100-batch1_rpkm_table.tsv"
rna_seq <- as.data.frame(read.csv(filePath, header = TRUE, sep = "\t"))
colnames(rna_seq) <- c("symbol",
                       "Fat-Treg-R1", "Fat-Treg-R2", "Fat-Treg-R3",
                       "Liver-Treg-R1", "Liver-Treg-R2", "Liver-Treg-R3",
                       "Lymph-N-Tcon-R1", "Lymph-N-Tcon-R2", "Lymph-N-Tcon-R3",
                       "Lymph-N-Treg-R1", "Lymph-N-Treg-R2", "Lymph-N-Treg-R3",
                       "Skin-Treg-R1", "Skin-Treg-R2", "Skin-Treg-R3"
)
# remove last column of NAs
rna_seq <- rna_seq[,1:16]
rna_seq$symbol <- as.character(rna_seq$symbol)
# clean data
exprData <- as.matrix(rna_seq[, -1])
# add annotations
colData <- DataFrame(tissue_cell=rep(c("Fat-Treg", "Liver-Treg", "Lymph-N-Tcon", "Lymph-N-Treg", "Skin-Treg"), 3),
                     row.names = colnames(exprData)
)
rowData <- DataFrame(id_symbol=rna_seq$symbol)
#
se <- SummarizedExperiment(assays=list(counts=exprData), colData=colData, rowData=rowData)
saveRDS(object = se, file = paste0(maindir, "RPKM_table.rds"))
rm(rna_seq)
### process htseq counts
filePath <- "/home/imbusch/git/treg_experimentHub/RNA-seq/D100_htseq_counts_header.txt"
rna_seq <- as.data.frame(read.csv(filePath, header = TRUE, sep = "\t"))
colnames(rna_seq) <- c("EnsemblID",
                       "Fat-Treg-R1", "Fat-Treg-R2", "Fat-Treg-R3",
                       "Liver-Treg-R1", "Liver-Treg-R2", "Liver-Treg-R3",
                       "Lymph-N-Tcon-R1", "Lymph-N-Tcon-R2", "Lymph-N-Tcon-R3",
                       "Lymph-N-Treg-R1", "Lymph-N-Treg-R2", "Lymph-N-Treg-R3",
                       "Skin-Treg-R1", "Skin-Treg-R2", "Skin-Treg-R3"
)
exprData <- as.matrix(rna_seq[, -1])
# add annotations
colData <- DataFrame(tissue_cell=rep(c("Fat-Treg", "Liver-Treg", "Lymph-N-Tcon", "Lymph-N-Treg", "Skin-Treg"), 3),
                     row.names = colnames(exprData)
)
se <- SummarizedExperiment(assays=list(counts=exprData), colData=colData)
saveRDS(object = se, file = paste0(maindir, "htseq_table.rds"))
