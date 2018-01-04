
## Each row of the metadata file represents a resource added to one of
## the 'Hubs'. This example creates a metadata.csv file for a single resource.
## In the case of multiple resources, the arguments below would be character
## vectors that produced multiple rows in the data.frame.

meta <- data.frame(
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
  ResourceName = "treg_twgbs.rds",
  RDataPath = "treg/SummarizedExperiment/",
  Location_Prefix = "http://computational-biology-public.imbusch.net/"
)

# meta <- data.frame(
#   Title = "TWGBS dataset from study PMID: 28783152",
#   Description = paste0("TWGBS data from study PMID: 28783152 containing regulatory T cells from four different tissues ",
#                        "and T conventional cells from a single tissue represented as a ",
#                        "SummarizedExperiment"),
#   BiocVersion = "3.4",
#   Genome = "mm10",
#   SourceType = "BAM",
#   SourceUrl = "",
#   SourceVersion = "Jan 03 2018",
#   Species = "Mus musculus",
#   TaxonomyId = 10090,
#   Coordinate_1_based = TRUE,
#   DataProvider = "",
#   Maintainer = "c.imbusch@dkfz.de",
#   RDataClass = "SummarizedExperiment",
#   DispatchClass = "HDF5Array",
#   ResourceName = "assays.h5",
#   RDataPath = "SummarizedExperiment_test/",
#   Location_Prefix = "http://computational-biology-public.imbusch.net/"
# )
#
# meta_append <- data.frame(
#   Title = "TWGBS dataset from study PMID: 28783152",
#   Description = paste0("TWGBS data from study PMID: 28783152 containing regulatory T cells from four different tissues ",
#                        "and T conventional cells from a single tissue represented as a ",
#                        "SummarizedExperiment"),
#   BiocVersion = "3.4",
#   Genome = "mm10",
#   SourceType = "BAM",
#   SourceUrl = "",
#   SourceVersion = "Jan 03 2018",
#   Species = "Mus musculus",
#   TaxonomyId = 10090,
#   Coordinate_1_based = TRUE,
#   DataProvider = "",
#   Maintainer = "c.imbusch@dkfz.de",
#   RDataClass = "SummarizedExperiment",
#   DispatchClass = "HDF5Array",
#   ResourceName = "se.rds",
#   RDataPath = "SummarizedExperiment_test/",
#   Location_Prefix = "http://computational-biology-public.imbusch.net/"
# )

## Write the data out and put in the inst/extdata directory.
write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
