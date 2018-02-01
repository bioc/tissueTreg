
## Each row of the metadata file represents a resource added to one of
## the 'Hubs'. This example creates a metadata.csv file for a single resource.
## In the case of multiple resources, the arguments below would be character
## vectors that produced multiple rows in the data.frame.

meta_bsseq_per_sample <- data.frame(
  Title = "TWGBS dataset from study PMID: 28783152",
  Description = paste0("TWGBS data from study PMID: 28783152 containing regulatory T cells from four different tissues ",
                       "and T conventional cells from a single tissue represented as a ",
                       "SummarizedExperiment"),
  BiocVersion = "3.4",
  Genome = "mm10",
  SourceType = "BAM",
  SourceUrl = "",
  SourceVersion = "Jan 03 2018",
  Species = "Mus musculus",
  TaxonomyId = 10090,
  Coordinate_1_based = TRUE,
  DataProvider = "",
  Maintainer = "c.imbusch@dkfz.de",
  RDataClass = "bsseq",
  DispatchClass = "rds",
  ResourceName = "treg_twgbs_per_sample.rds",
  RDataPath = "treg/SummarizedExperiment/",
  Location_Prefix = "http://computational-biology-public.imbusch.net/"
)

meta_bsseq_per_group <- data.frame(
  Title = "TWGBS dataset from study PMID: 28783152",
  Description = paste0("TWGBS data from study PMID: 28783152 containing regulatory T cells from four different tissues ",
                       "and T conventional cells from a single tissue represented as a ",
                       "SummarizedExperiment"),
  BiocVersion = "3.4",
  Genome = "mm10",
  SourceType = "BAM",
  SourceUrl = "",
  SourceVersion = "Jan 03 2018",
  Species = "Mus musculus",
  TaxonomyId = 10090,
  Coordinate_1_based = TRUE,
  DataProvider = "",
  Maintainer = "c.imbusch@dkfz.de",
  RDataClass = "bsseq",
  DispatchClass = "rds",
  ResourceName = "treg_twgbs_per_group.rds",
  RDataPath = "treg/SummarizedExperiment/",
  Location_Prefix = "http://computational-biology-public.imbusch.net/"
)

meta_rna_seq_per_per_sample_rpkm <- data.frame(
  Title = "RNA-seq dataset from study PMID: 28783152",
  Description = paste0("RNA-seq data from study PMID: 28783152 containing regulatory T cells from four different tissues ",
                       "and T conventional cells from a single tissue represented as a ",
                       "SummarizedExperiment"),
  BiocVersion = "3.4",
  Genome = "mm10",
  SourceType = "BAM",
  SourceUrl = "",
  SourceVersion = "Jan 03 2018",
  Species = "Mus musculus",
  TaxonomyId = 10090,
  Coordinate_1_based = TRUE,
  DataProvider = "",
  Maintainer = "c.imbusch@dkfz.de",
  RDataClass = "SummarizedExperiment",
  DispatchClass = "rds",
  ResourceName = "RPKM_table.rds",
  RDataPath = "treg/SummarizedExperiment/",
  Location_Prefix = "http://computational-biology-public.imbusch.net/"
)

meta_rna_seq_per_per_sample_htseq <- data.frame(
  Title = "RNA-seq dataset from study PMID: 28783152",
  Description = paste0("RNA-seq data from study PMID: 28783152 containing regulatory T cells from four different tissues ",
                       "and T conventional cells from a single tissue represented as a ",
                       "SummarizedExperiment"),
  BiocVersion = "3.4",
  Genome = "mm10",
  SourceType = "BAM",
  SourceUrl = "",
  SourceVersion = "Jan 03 2018",
  Species = "Mus musculus",
  TaxonomyId = 10090,
  Coordinate_1_based = TRUE,
  DataProvider = "",
  Maintainer = "c.imbusch@dkfz.de",
  RDataClass = "SummarizedExperiment",
  DispatchClass = "rds",
  ResourceName = "htseq_table.rds",
  RDataPath = "treg/SummarizedExperiment/",
  Location_Prefix = "http://computational-biology-public.imbusch.net/"
)

## Write the data out and put in the inst/extdata directory.
write.csv(rbind(meta_bsseq_per_sample,
                meta_bsseq_per_group,
                meta_rna_seq_per_per_sample_rpkm),
          file="inst/extdata/metadata.csv",
          row.names=FALSE)
