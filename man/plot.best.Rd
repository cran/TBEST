\name{plot.best}
\alias{plot.best}
\alias{plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a dendrogam of a hierarchical cluster with branches labeled by their numbers and significance estimates of tightness.
}
\description{
Description: This function provides a visual view of the hierarchical structure, together with the computed pvalue for tight branches.
}
\usage{
\method{plot}{best}(x,mystat="fldc",sigp=0.05,partition=NA,print.num=TRUE,
	print.lab=TRUE,float=0.01,col.best=c(2,3),cex.best=0.8,font.best=NULL,
	main=NULL,sub=NULL,xlab=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{\link{best}}, such as computed by the \code{\link{SigTree}} function.
}
\item{mystat}{
A statistic for which the p-values are to be indicated. Default is \code{"fldc"}. Other options are \code{"fldcc"} and \code{"bldc"}.
}
\item{sigp}{
Significance level for tightness of branches when \code{partition=NA}, default is 0.05.
}
\item{partition}{
An object of class \code{partition}, such as computed by the \code{\link{PartitionTree}} function. If this is NA(default), all significant branches are plotted instead of a partition. See below details. 
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
A character vector of length 2, indicating the colors to be used for the p-values and for the branch numbers.
}
\item{cex.best}{
A numeric value for the text size of the branch labels.
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
\item{...}{
Further arguments to be passed on to the \code{plot} function.
}
}
\details{
The plot gives the branch numbers and pvalues in the hierarchical tree. When variable \code{partition} in the input object \code{x} is a two-column index table, this plot gives the pvalues of branches that form the partition. Otherwise, this function puts legends on all tight branches with pvalue \code{<} \code{sigp}. To obtain the members within certain branch, refer to function \code{\link{FindMembership}}.
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
Guoli Sun
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
        distrib="Rparallel",njobs=2,Ptail=TRUE)
plot(x=mytable,mystat="fldc")
mypartition<-PartitionTree(x=mytable,sigp=0.0001,statname="fldc",fdr=FALSE)
plot(x=mytable,mystat="fldc",partition=mypartition)
plot(x=mytable,mystat="fldc",partition=mypartition,print.num=F)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
