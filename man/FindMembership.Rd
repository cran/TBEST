\name{FindMembership}
\alias{FindMembership}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
List items inside branch/branches in a hierarchical tree
}
\description{
Description: Use this function to find the labels of all items comprising one or more branches of a hierarchical tree.
}
\usage{
FindMembership(myinput, mynode, method, metric)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{myinput}{
Input data table should be numerical and in data class "matrix", with columns of factors and rows of samples.
}
  \item{mynode}{
A numeric vector of the branch numbers you would like to view in the hierarchical tree. The branch number should be less than or equal to Sample Size-1. Negative value means it is an individual object.
}
  \item{method}{
A character string specifying the linkage method for hierarchical clustering, to be used by the \code{hclust} function. See \code{?hclust} argument \code{method} for method options.
}
  \item{metric}{
A character string specifying hierarchical clustering dissimilarity metric.Choices can be correlation based, \code{"pearson"},\code{"kendall"}, or \code{"spearman"}. Choices can also be distance measures to be used by \code{dist} function. See \code{?dist} argument \code{method} for method options.
}
}
\details{
This function can be used to find sample labels in significant sub-clusters detected using function \code{\link{SigTree}}.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of items. Each item is a vector of all members inside the branch, in terms of sample labels.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Guoli Sun
}
\examples{
data(leukemia)
#find an individual item name
FindMembership(data.matrix(leukemia),mynode=c(-29),method="ward",
       metric="euclidean")
#find all members inside a sub-cluster
FindMembership(data.matrix(leukemia),mynode=c(-29,29),
       method="ward",metric="euclidean")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
%R documentation directory.
