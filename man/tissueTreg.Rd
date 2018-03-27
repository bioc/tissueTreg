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
\usage{data("tissueTreg")}
\format{
  A data frame with 0 observations on the following 2 variables.
  \describe{
    \item{\code{x}}{a numeric vector}
    \item{\code{y}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Delacher and Imbusch 2017, PMID: 28783152
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
library(ExperimentHub)
library(bsseq)
eh <- ExperimentHub()
BS.obj.ex.fit <- eh[["EH1072"]]
}
\keyword{datasets}
