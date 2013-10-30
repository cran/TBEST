\name{PartitionTree}
\alias{PartitionTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluate branches with given significance level
}
\description{
Description: The evaluation process finds the most detailed partition of the tree into significantly tight branches, given the significance level, and a choice of a statistic for tightness. This function returns an object of class "partition".
}
\usage{
PartitionTree(x,sigp=0.05,statname="fldc",fdr=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{\link{best}}, such as computed by function \code{\link{SigTree}}.
}
  \item{sigp}{
Significance level for tightness of branches, default is 0.05.
}
  \item{statname}{
A character string specifying the name of statistic for which the significance of tightness is computed. The choices are \code{"fldc"}(default),\code{"bldc"},\code{"fldcc"}.
}
  \item{fdr}{
A logical value. If \code{fdr=TRUE}, \code{sigp} is interpreted as false discovery rate.
}
}
\value{
Output returns an object of class \code{\link{partition}}. See \code{?partition} for details.
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
        mymethod="ward",mymetric="euclidean",rand.fun="shuffle.factor",
        distrib="Rparallel",njobs=2,Ptail=TRUE)
mypartition<-PartitionTree(x=mytable,sigp=0.0001,statname="fldc",fdr=FALSE)
partition1<-mypartition$partition
fix(partition1)
mypartition<-PartitionTree(x=mytable,sigp=0.05,statname="fldc",fdr=FALSE)
}
}
