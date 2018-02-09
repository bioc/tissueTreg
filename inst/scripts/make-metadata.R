
## Each row of the metadata file represents a resource added to one of
## the 'Hubs'. This example creates a metadata.csv file for a single resource.
## In the case of multiple resources, the arguments below would be character
## vectors that produced multiple rows in the data.frame.

meta_bsseq_per_sample <- data.frame(
  Title = "Bisulfite sequencing data from tissue Tregs (per sample)",
  Description = paste0("Bisulfite sequencing data (PMID: 28783152) containing regulatory T cells from four different tissues ",
                       "and T conventional cells from a single tissue represented as a ",
                       "bsseq object"),
  BiocVersion = "3.7",
  Genome = "mm10",
  SourceType = "BAM",
  SourceUrl = "NA",
  SourceVersion = "Feb 08 2018",
  Species = "Mus musculus",
  TaxonomyId = 10090,
  Coordinate_1_based = TRUE,
  DataProvider = "DKFZ",
  Maintainer = "c.imbusch@dkfz.de",
  RDataClass = "bsseq",
  DispatchClass = "RDS",
  RDataPath = "tissueTreg/treg_twgbs_per_sample.rds"
)

meta_bsseq_per_group <- data.frame(
  Title = "Bisulfite sequencing data from tissue Tregs (per tissue/cell type)",
  Description = paste0("Bisulfite sequencing data (PMID: 28783152) containing regulatory T cells from four different tissues ",
                       "and T conventional cells from a single tissue represented as a ",
                       "bsseq object"),
  BiocVersion = "3.7",
  Genome = "mm10",
  SourceType = "BAM",
  SourceUrl = "NA",
  SourceVersion = "Feb 08 2018",
  Species = "Mus musculus",
  TaxonomyId = 10090,
  Coordinate_1_based = TRUE,
  DataProvider = "DKFZ",
  Maintainer = "c.imbusch@dkfz.de",
  RDataClass = "bsseq",
  DispatchClass = "RDS",
  RDataPath = "tissueTreg/treg_twgbs_per_group.rds"
)

meta_rna_seq_per_per_sample_rpkm <- data.frame(
  Title = "RNA-seq data from tissue Tregs (RPKM values)",
  Description = paste0("RNA-seq data (PMID: 28783152) containing regulatory T cells from four different tissues ",
                       "and T conventional cells from a single tissue represented as a ",
                       "SummarizedExperiment"),
  BiocVersion = "3.7",
  Genome = "mm10",
  SourceType = "BAM",
  SourceUrl = "NA",
  SourceVersion = "Feb 08 2018",
  Species = "Mus musculus",
  TaxonomyId = 10090,
  Coordinate_1_based = TRUE,
  DataProvider = "DKFZ",
  Maintainer = "c.imbusch@dkfz.de",
  RDataClass = "SummarizedExperiment",
  DispatchClass = "RDS",
  RDataPath = "tissueTreg/RPKM_table.rds"
)

meta_rna_seq_per_per_sample_htseq <- data.frame(
  Title = "RNA-seq data from tissue Tregs (htseq values)",
  Description = paste0("RNA-seq data (PMID: 28783152) containing regulatory T cells from four different tissues ",
                       "and T conventional cells from a single tissue represented as a ",
                       "SummarizedExperiment"),
  BiocVersion = "3.7",
  Genome = "mm10",
  SourceType = "BAM",
  SourceUrl = "NA",
  SourceVersion = "Feb 08 2018",
  Species = "Mus musculus",
  TaxonomyId = 10090,
  Coordinate_1_based = TRUE,
  DataProvider = "DKFZ",
  Maintainer = "c.imbusch@dkfz.de",
  RDataClass = "SummarizedExperiment",
  DispatchClass = "RDS",
  RDataPath = "tissueTreg/htseq_table.rds"
)

## Write the data out and put in the inst/extdata directory.
write.csv(rbind(meta_bsseq_per_sample,
                meta_bsseq_per_group,
                meta_rna_seq_per_per_sample_rpkm,
                meta_rna_seq_per_per_sample_htseq),
          file="inst/extdata/metadata.csv",
          row.names=FALSE)
