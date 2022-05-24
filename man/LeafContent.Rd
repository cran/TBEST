\name{LeafContent}
\alias{LeafContent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find names of leaves belongigng to given branches of a hierarchical tree
}
\description{
Description: find the names of all items comprising one or more branches of a hierarchical tree.
}
\usage{
LeafContent(myinput, mynode=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{myinput}{
An object of class \code{hclust}, \code{\link{best}} or \code{\link{partition}}.
}
  \item{mynode}{
An integer vector of the numbers of branches whose leaf content is desired. The \code{hclust} convention is used for numbering branches and leaves, i.e., the branch numbers can take any value between (-N) and (N-1) excluding 0, where N is the number of leaves in the tree. A negative value refers to an individual leaf whose number is minus that value. If \code{myinput} is of class \code{partition}, this argument is ignored. The function lists the leaf content for each of the branches that form the partition.
}
}
%\details{
%This function can be used to find sample labels in significant sub-clusters detected using function \code{\link{SigTree}} and \code{\link{PartitionTree}}.
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list of items, of the same length as \code{mynode}. Each item corresponds to a branch listed in \code{myneode} and is a character vector containing the names of the leaves in the branch.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Guoli Sun, Alex Krasnitz
}
\examples{
data(leukemia)
hc<-hclust(dist(data.matrix(leukemia)),"ward")
#find the name of leaf 29
LeafContent(hc,mynode=c(-29))
#find the name of leaf 29 and leaves belonging to node 29
LeafContent(hc,mynode=c(-29,29))
\dontrun{
mytable<-SigTree(data.matrix(leukemia),mystat="fldc",
        mymethod="ward",mymetric="euclidean",rand.fun="shuffle.column",
        distrib="Rparallel",njobs=2,Ptail=TRUE,tailmethod="ML")
LeafContent(mytable,mynode=c(-29,29))
mypartition<-PartitionTree(x=mytable,siglevel=0.001,statname="fldc",sigtype="raw")
LeafContent(mypartition)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
%R documentation directory.
