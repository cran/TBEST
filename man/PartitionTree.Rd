\name{PartitionTree}
\alias{PartitionTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the most detailed partition of a tree into tight branches.
}
\description{
Description: The function finds the most detailed partition of a hierartchical tree into tight branches, given a level of significance for tightness.
}
\usage{
PartitionTree(x,siglevel=0.05,statname="fldc",
        sigtype=c("raw","corrected","fdr"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{\link{best}}, such as computed by function \code{\link{SigTree}}.
}
  \item{siglevel}{
Threshold of significance for tightness of branches. Default is 0.05.
}
  \item{statname}{
A character string specifying the name of measure of tighness whose is significance is to be used for partition. The choices are \code{"fldc"}(default),\code{"bldc"},\code{"fldcc"}.
}
  \item{sigtype}{
A character string specifying how the significance threshold \code{siglevel} should be interpreted. If \code{"raw"}, the threshold will be applied directly to the p-values tabulated for each tree node in \code{x}. With \code{"corrected"} chosen, the threshold will be applied to the p-values corrected for multiplicity: \cr \code{p_cor = 1 - (1 - p)^(N - 2)}, where \code{N} is the number of leaves in the tree. of significance. If \code{"fdr"}, \code{siglevel} is interpreted as a threshold on false discovery rate.
}
}
\value{
An object of class \code{\link{partition}}. See \code{?partition} for details.
}
\author{
Guoli Sun, Alex Krasnitz
}
%\note{
%When using false discovery detection, package \code{qvalue} should be installed in advance.
%}
\seealso{
\code{\link{SigTree}},\code{\link{partition}},\code{\link{best}}
}
\examples{
\dontrun{
data(leukemia)
mytable<-SigTree(data.matrix(leukemia),mystat="all",
        mymethod="ward",mymetric="euclidean",rand.fun="shuffle.column",
        distrib="Rparallel",njobs=2,Ptail=TRUE,tailmethod="ML")
mypartition<-PartitionTree(x=mytable,siglevel=0.001,statname="fldc",
        sigtype="raw")
partition1<-mypartition$partition
sigmatrix1<-mypartition$sigvalue
fix(partition1)
fix(sigmatrix1)
}
}
