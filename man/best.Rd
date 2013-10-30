\name{best}
\alias{best}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An object of class "best"
}
\description{
Description: This object is a list of four items. It contains a statistical assessment of the tightness of branches in a hierarchical tree.
}
\value{
\item{data}{
A matrix from which the distance matrix for growing the tree is computed, with the rows corresponding to the items being clustered.
}
\item{method}{
A character string specifying the linkage used to grow the tree.
}
\item{metric}{
A character string specifying the method used to compute distances among the items.
}
\item{indextable}{
A matrix with X columns and the number of rows one less than the number of items being clustered. The columns are as follows: First four columns specify internal tree structure, \code{"index1"},\code{"index2"},\code{"height"},\code{"clustersize"}; The following columns are chosen statistics and corresponding p-values. For example, if we choose statistic \code{"fldc"}, the following columns are \code{"fldc"}, \code{"pfldc"}.
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
\code{\link{SigTree}},\code{\link{plot.best}}
}
\examples{
\dontrun{
data(leukemia)
mytable<-SigTree(data.matrix(leukemia),mystat="all",
       mymethod="ward",mymetric="euclidean",rand.fun="shuffle.column",
       distrib="Rparallel",njobs=2,Ptail=TRUE)
class(mytable)
names(mytable)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
