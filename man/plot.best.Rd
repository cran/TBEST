\name{plot.best}
\alias{plot.best}
\alias{plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a dendrogam of a hierarchical cluster with branches labeled by their numbers and significance estimates of tightness.
}
\description{
Description: A plot method for the class \code{\link{best}}.
}
\usage{
\method{plot}{best}(x,mystat="fldc",siglevel=0.05,sigtype=c("raw","corrected","fdr"),
        partition=NA,print.num=TRUE,print.lab=TRUE,float=0.01,col.best=c(2,3),
        cex.best=0.8,cex.leaf=0.8,font.best=NULL,main=NULL,sub=NULL,xlab=NULL,
        metric.args=list(),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{\link{best}}, such as computed by the \code{\link{SigTree}} function.
}
\item{mystat}{
A measure of tightness for which p-values are to be shown in the plot. Default is \code{"fldc"}. Other options are \code{"fldcc"} and \code{"bldc"}.
}
\item{siglevel}{
A threshold level of significance for tightness of branches used when \code{partition=NA}. Default is 0.05. If the estimate of significance for a node is below threshold, it will be shown on the plot next to the node.
}
\item{sigtype}{
A character string specifying how the significance threshold \code{siglevel} should be interpreted. If \code{"raw"}, the threshold will be applied directly to the p-values tabulated for each tree node in \code{x}. With \code{"corrected"} chosen, the threshold will be applied to the p-values corrected for multiplicity:\cr \code{p_cor = 1 - (1 - p)^(N - 2)}, where \code{N} is the number of leaves in the tree. of significance. If \code{"fdr"}, \code{siglevel} is interpreted as a threshold on false discovery rate.
}
\item{partition}{
An object of class \code{\link{partition}}, such as computed by the \code{\link{PartitionTree}} function.
} 
\item{print.num}{
Logical. If true, the branch numbers will be indicated.
}
\item{print.lab}{
Logical. If true, the labels will be displayed at the bottom of dendrogram.
}
\item{float}{
A numeric value that can change the vertical location of pvalues.
}
\item{col.best}{
A character vector of length 2, indicating the colors to be used for the p-values and for the numbers of the nodes.
}
\item{cex.best}{
A numeric value for the text size of the branch labels.
}
\item{cex.leaf}{
A numeric value for the text size of the leaf labels.
}
\item{font.best}{
An integer which specifies font choice of text on the plot. See \code{?par} function parameter \code{font} for details.
}
\item{main}{
A character string specifying the title of the plot.
}
\item{sub}{
A character string specifying a subtitle of the plot.
}
\item{xlab}{
A caracter string specifying the label of horizontal axis.
}
\item{metric.args}{
Additional argument from user supplied dissimilariity(distance) function. See \code{details} and \code{examples} below for further explanation.
}
\item{...}{
Further arguments to be passed on to the \code{plot} function.
}
}
\details{
The function plots a dendrogram of the hierarchical tree as specified by the \code{x} argument, an object of class \code{"best"}. When argument \code{partition} is set to an object of class \code{"partition"}, and a partition does exist (see \code{\link{partition}} for description), this plot provides the significance estimates for the nodes that form the partition. Otherwise, this function puts legends on all tight nodes with significance estimates no more than \code{siglevel}. To obtain the leaves descending from a given node, refer to function \code{\link{LeafContent}}.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A plot with all branch numbers and significant pvalues in the hierarchical tree.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Guoli Sun, Alex Krasnitz
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SigTree}}, \code{\link{PartitionTree}},\code{\link{best}},\code{\link{partition}}
}
\examples{
\dontrun{
data(leukemia)
mytable<-SigTree(data.matrix(leukemia),mystat="all",
        mymethod="ward",mymetric="euclidean",rand.fun="shuffle.column",
        distrib="Rparallel",njobs=2,Ptail=TRUE,tailmethod="ML")
plot(x=mytable,mystat="fldc",siglevel=0.001,sigtype="raw")
mypartition<-PartitionTree(x=mytable,siglevel=0.001,statname="fldc",
        sigtype="raw")
plot(x=mytable,mystat="fldc",partition=mypartition)
plot(x=mytable,mystat="fldc",partition=mypartition,print.num=F)
#with user-defined functions
mydist<-function(x,y){return(dist(x)/y)}
myrand<-function(x,z){return(apply(x+z,2,sample))}
mytable<-SigTree(data.matrix(leukemia),mystat="fldc",
	mymethod="ward",mymetric="mydist",rand.fun="myrand",
	distrib="Rparallel",njobs=2,Ptail=TRUE,tailmethod="MOM",metric.args=list(3),
	rand.args=list(2))
plot(mytable,metric.args=list(3))
plot(mytable,metric.args=list(3),cex.leaf=1.5)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
