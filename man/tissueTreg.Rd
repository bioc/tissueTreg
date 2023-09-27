\name{tissueTreg}
\alias{tissueTreg}
\docType{data}
\title{
Epigenomes and transcriptomes of tissue resident regulatory T cells
}
\description{
The package provides ready to use epigenomes (obtained from TWGBS)
and transcriptomes (RNA-seq) from various tissues.
Regulatory T cells (Treg cells) perform two distinct functions: they maintain
self-tolerance, and they support organ homeostasis by differentiating into
specialized tissue Treg cells. The underlying dataset characterises the epigenetic
and transcriptomic modifications for specialized tissue Treg cells.
}
\source{
Delacher, M, Imbusch, CD, Weichenhan, D, Breiling, A, Hotz-Wagenblatt, A, Träger, U, Hofer, AC, Kägebein, D, Wang, Q, Frauhammer, F, Mallm, JP, Bauer, K, Herrmann, C, Lang, PA, Brors, B, Plass, C, Feuerer, M (2017). Genome-wide DNA-methylation landscape defines specialization of regulatory T cells in tissues. Nat. Immunol., 18, 10:1160-1172.
}
\examples{
eh <- ExperimentHub::ExperimentHub()

# RNA-seq RPKM data:
se_rpkms <- eh[["EH1074"]]

# Whole genome bisulfite sequencing data as bsseq objects:
tregs_per_sample <- eh[["EH1072"]]
tregs_per_tissue <- eh[["EH1073"]]
}
\keyword{datasets}
