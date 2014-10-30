\name{best}
\alias{best}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An object of class "best"
}
\description{
Description: This object is a list of three items. It contains a statistical assessment of the tightness of branches in a hierarchical tree.
}
\value{
\item{Call}{
An object of class Call, specifying the parameters used.
}
\item{data}{
A matrix from which the distance matrix used for growing the tree is computed, with the rows corresponding to the items being clustered.
}
\item{indextable}{
If measure of tightness is not \code{"slb"}, this is a matrix with the number of rows one less than the number of items being clustered.  Each row corresponds to an internal node in the tree. The columns are as follows. First two columns specify the merging order of the tree, as in the \code{merge} component of the class \code{hclust}. The third column contains the node heights, as in the \code{height} component of \code{hclust}. The fourth column provides the number of leaves for each node. The corresponding column names are \code{"index1"}, \code{"index2"}, \code{"height"}, \code{"clustersize"}. The remaining columns come in pairs. If the name of the first column in a pair is \code{"x"}, the name of the second one is \code{"px"}. The first column in each pair tabulates a mesure of tightness; the second column provides the corresponding p-value. If measure of tightness is \code{"slb"}, this is a list with two variable, a matrix like above except without column of p-value and a p-value suggesting the significance of two-way split of input data.
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
\code{\link{SigTree}},\code{\link{plot.best}}
}
\examples{
\dontrun{
data(leukemia)
mytable<-SigTree(data.matrix(leukemia),mystat="all",
       mymethod="ward",mymetric="euclidean",rand.fun="shuffle.column",
       distrib="Rparallel",njobs=2,Ptail=TRUE,tailmethod="ML")
class(mytable)
names(mytable)
mytable<-SigTree(data.matrix(leukemia),mystat="slb",
       mymethod="ward",mymetric="euclidean",rand.fun="shuffle.column",
       distrib="Rparallel",njobs=2,Ptail=FALSE)
class(mytable)
names(mytable)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
