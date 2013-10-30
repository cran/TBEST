\name{leukemia}
\alias{leukemia}
\docType{data}
\title{
Leukemia data
}
\description{
This data set represents mRNA expression of 500 genes in 38 patient cases of leukemia. These 38 cases fall into 3 subtypes: AML (11), T-lineage ALL (8) and B-lineage ALL (19). The set was obtained by removing 499 genes from Golub's leukemia data, to facilitate the execution of examples for this package.
}
\usage{data(leukemia)}
\format{
A data frame with 38 observations (rows) of 500 variables (columns).
}
\details{
Bone marrow samples obtained from acute leukemia patients at the time of diagnosis.
}
\source{
http://www.broadinstitute.org/cgi-bin/cancer/datasets.cgi
}
\references{
T.R. Golub, D.K. Slonim et al(1999)
Molecular Classification of Cancer: Class Discovery and Class Prediction by Gene Expression;				

Stefano Monti, Pablo Tamayo, Jill Mesirov, and Todd Golub(2003)
Consensus Clustering: A resampling-based method for class discovery and visualization of gene expression microarray data
}
\examples{
data(leukemia)
dim(leukemia)
}
