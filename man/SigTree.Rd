\name{SigTree}
\alias{SigTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Branches evaluated statistically for tightness with hierarchical clustering.
}
\description{
Description: This function computes tightness statistics, and their corresponding permutation p-values, for branches of a hierarchical tree.
}
\usage{
SigTree(myinput,mystat=c("all","fldc","bldc","fldcc"),
	mymethod="complete",mymetric="euclidean",
	rand.fun=NA,by.block=NA,
	distrib=c("vanilla","Rparallel"),Ptail=TRUE,
	tailmethod="ML",njobs=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{myinput}{
A matrix with rows corresponding to items to be clustered.
}
  \item{mystat}{
A character string specifying the names of statistics to be computed and evaluated for significance. If \code{"all"} is chosen, all statistics and p-values are computed. Otherwise, only the specified statistic and its p-value are computed.
%\code{"fldc"},\code{"bldc"},\code{"fldcc"} are computed based on the parent and children branches' height.
}
  \item{mymethod}{
A character string specifying the linkage method for hierarchical clustering, to be used by the \code{hclust} function. See \code{?hclust} argument \code{method} for method options.
}
  \item{mymetric}{
A character string specifying hierarchical clustering dissimilarity metric.Choices can be correlation based, \code{"pearson"},\code{"kendall"}, or \code{"spearman"}. Choices can also be distance measures to be used by \code{dist} function. See \code{?dist} argument \code{method} for method options. 
}
  \item{rand.fun}{
A character string specifying the permutation method to be applied to \code{myinput}. If NA(default), no permutation is performed. \code{"shuffle.column"} performs a random permutation independently within each column. With \code{"shuffle.block"}, a random permutation is performed independently within each block of columns, as specified by the \code{by.block} argument, and independently from the other blocks. It can also be a string contains any user supplied randomization function for \code{myinput}. For more details, see details and examples below.
}
  \item{by.block}{
A numeric vector of the same length as the column dimension of \code{myinput}, to specify the blocking of columns of \code{myinput}. It is used in conjunction with \code{rand.fun = "shuffle.block"} and optionally with a user-defined permutation method, and is ignored otherwise.  
}
  \item{distrib}{
One of \code{"vanilla", "Rparallel"} to specify the distributed computing option for the cluster assignment step. For \code{"vanilla"} (default)
no distributed computing is performed. For \code{"Rparallel"} the \code{parallel} package of \code{R} core is used for multi-core processing.
}
  \item{Ptail}{
Logical. If \code{Ptail} is TRUE(default), the Generalized Pareto Distribution is used to approximate the null distribution of the chosen statistics. Otherwise, empirical p-values are computed directly from the permutation test.
}
  \item{tailmethod}{
A character string only needed to be specified if the \code{Ptail} is set to TRUE. \code{"ML"} uses maximum likelihood to estimate the parameters of Generalized Pareto Distribution. \code{"MOM"} uses method of moments to estimate those parameters. 
}
  \item{njobs}{
A single integer specifying the number of worker jobs to create in case of distributed computation.
}
}
\details{
When \code{rand.fun} is set to a user supplied randomization function, the first parameter should be used as the input \code{myinput}, and other parameter values should be given in advance. See examples below.
 
The function performs statistical evaluation on a number of different statistics. Given the hierarchical structure, a branch \code{a} has a sibling branch \code{b}, with which share the same parent branch \code{p}. We compute the statistics on measures of branches' height. Assume branch \code{a}, \code{b} and \code{p} have height \code{ha},\code{hb},\code{hp}, respectively. Annotation for statistic of the child branch is \code{Sa}, and the parent branch is \code{Sp}.  

fldc:

\code{Sa = (hp-ha)/hp}

fldcc:

\code{Sa = (hp-(ha-hb)/2)/ha}

bldc:

\code{Sp = (2*hp-ha-hb)/(2*hp)}

The tightness of each branch is computed by chosen statistics, and further evaluated by permutation p-values.
With \code{Ptail} equals TRUE(default), a total number of 1000 permutation tests are run. Otherwise, a number of \code{N}*1000 permutation tests are run (\code{N} is the sample size).
}
\value{
If \code{rand.fun} is set to NA, then this function returns the internal tree structure numeric table with extra columns of specified statistics.
If \code{rand.fun} is set to a specific permutation method, output returns an object of class \code{\link{best}}. See \code{?best} for details.
}
\references{
Theo A. Knijnenburg, Lodewyk F. A. Wessels et al (2009)
Fewer permutations, more accurate P-values
}
\author{
Guoli Sun, Alex Krasnitz
}
\note{
%When using generalized pareto p-value estimation, package \code{signal} should be installed in advance.
If \code{rand.fun} is a customized randomization function, make sure you have read and write permission at your local working environment.
}
\seealso{
\code{\link{best}},\code{\link{plot.best}}
}
\examples{
####leukemia data, with ground truth of three subtypes
data(leukemia)
#output only statistic table
mytable<-SigTree(data.matrix(leukemia),mystat="all",
        mymethod="ward",mymetric="euclidean")
class(mytable)
\dontrun{
#use multicore processing to detect significant sub-clusters
mytable<-SigTree(data.matrix(leukemia),mystat="all",
	mymethod="ward",mymetric="euclidean",rand.fun="shuffle.column",
	distrib="Rparallel",njobs=2,Ptail=TRUE)
class(mytable)
####Breast tumor single cells data, ground truth of four subtypes
data(T10)
#This data set contains chromosome information, 
#then we can perform randomization within chromosomes
chrom<-as.numeric(T10[1,])
mydata<-T10[-1,] 
mytable<-SigTree(data.matrix(mydata),mystat="fldc",        
	mymethod="ward",mymetric="euclidean",rand.fun="shuffle.block",
	by.block=chrom,distrib="Rparallel",njobs=2,Ptail=TRUE,tailmethod="ML")
#using user supplied randomization function
myrand<-function(x,y=2){
	return(apply(x+y,2,sample))
}
mytable<-SigTree(data.matrix(leukemia),mystat="fldc",
        mymethod="ward",mymetric="euclidean",rand.fun="myrand",
	distrib="Rparallel",njobs=2,Ptail=TRUE,tailmethod="MOM")
}
}
