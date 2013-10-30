\name{partition}
\alias{partition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An object of class "partition"
}
\description{
Description: This object contains four variables, including info on detailed partition of tree.
}
\value{
\item{best}{
An object of class "best", see \code{\link{best}} for more info.
}
\item{statistic}{
A character string suggesting which statistic is used to evaluate the branches.
}
\item{pvalue}{
A numeric value specifying the significance level.
}
\item{partition}{
If a most detailed partition of the tree exists, this is a two-column index table. First column is the sample label, and the second column is the corresponding subtype index, in terms of branch number. If there's no such partition, this is NA.
}
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Guoli Sun
}
\seealso{
\code{\link{PartitionTree}},\code{\link{best}},\code{\link{SigTree}}
}
\examples{
\dontrun{
data(leukemia)
mytable<-SigTree(data.matrix(leukemia),mystat="all",
       mymethod="ward",mymetric="euclidean",rand.fun="shuffle.factor",
       distrib="Rparallel",njobs=2,Ptail=TRUE)
class(mytable)
mypartition<-PartitionTree(x=mytable,sigp=0.0001,statname="fldc",fdr=FALSE)
class(mypartition)
names(mypartition)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
