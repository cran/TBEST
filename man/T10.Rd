\name{T10}
\alias{T10}
\docType{data}
\title{
Breast tumor single cells data
}
\description{
This data set summarizes DNA copy number variation in 100 individual cancer cells harvested from a breast tumor. The cells belong to four subtypes, differing by ploidy.There are 47 Diploid+Pseudo-diploid, 24 Hypo-diploid, 4 Aneuploid B and 25 Aneuploid A cells. Their copy number profiles are summarized in terms of 354 amplification and deletion "cores", are computed by the CORE package.
}
\usage{data(T10)}
\format{
   A data frame with 101 rows and 354 columns. Each column corresponds to a core. The first row contains the chromosome number for each core. Then from second row, each row represent a single cancer cell.
}
\details{
Please remove the first row before computing the distance matrix.
}
\source{
Alexander Krasnitz, Guoli Sun, Peter Andrews, and Michael Wigler(2013)
Target inference from collections of genomic intervals
}
\references{
Alexander Krasnitz, Guoli Sun, Peter Andrews, and Michael Wigler(2013)
Target inference from collections of genomic intervals		
}
\examples{
data(T10)
dim(T10)
}
