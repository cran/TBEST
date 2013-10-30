\name{partition}
\alias{partition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An object of class "partition"
}
\description{
Description: This object is a list of four items, which jointly specify a detailed partition of a hierarchical tree into tight branches.
}
\value{
\item{Call}{
An object of class Call, specifying the function call which generated the list.
}
\item{best}{
An object of class "best", see \code{\link{best}} for more info.
}
\item{sigvalue}{
A two-column matrix, with one row per each internal node of the tree. The first column enumerates the nodes. The second column profides the significance estimate for the tightness of the node.
}
\item{partition}{
A two-column data frame specifying the partitition. The first column is a character vector with the names of the leaves. The second column provides the number of the part to which the leaf belongs.
}
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Guoli Sun, Alex Krasnitz
}
\seealso{
\code{\link{PartitionTree}},\code{\link{best}},\code{\link{SigTree}}
}
\examples{
\dontrun{
data(leukemia)
mytable<-SigTree(data.matrix(leukemia),mystat="all",
       mymethod="ward",mymetric="euclidean",rand.fun="shuffle.column",
       distrib="Rparallel",njobs=2,Ptail=TRUE,tailmethod="ML")
class(mytable)
mypartition<-PartitionTree(x=mytable,siglevel=0.001,statname="fldc",
       sigtype="raw")
class(mypartition)
names(mypartition)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
