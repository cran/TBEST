\name{SigTree}
\alias{SigTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform statistical analysis of tightness for branches of a hierarchical cluster.
}
\description{
Description: Given data from which a hierarchical tree is grown, compute measures of tightness for each branch, sample from the null distribution of these measures in the randomized data and compute the corresponding p-values.
}
\usage{
SigTree(myinput,mystat=c("all","fldc","bldc","fldcc","slb"),
        mymethod="complete",mymetric="euclidean",rand.fun=NA,
        by.block=NA,distrib=c("vanilla","Rparallel"),Ptail=TRUE,
        tailmethod=c("ML","MOM"),njobs=1,seed=NA,
        Nperm=ifelse(Ptail,1000,1000*nrow(myinput)),
        metric.args=list(),rand.args=list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{myinput}{
A matrix with rows corresponding to items to be clustered.
}
  \item{mystat}{
A character string specifying the measures of tightness to be computed and evaluated for significance of finding. See Details for the definitions of these measures. If \code{"all"} is chosen, all the first three measures, \code{"fldc"}, \code{"bldc"} and \code{"fldcc"}, and the corresponding p-values are computed. Otherwise, only the specified measure and its p-value are computed.
}
  \item{mymethod}{
A character string specifying the linkage method for hierarchical clustering, to be used by the \code{hclust} function. See \code{hclust} argument \code{method} for method options.
}
  \item{mymetric}{
A character string specifying the definition of dissimilarity (distance) among the data items. The options, in addition to those for the argument \code{method} of the \code{dist} functiton, are \code{"pearson"},\code{"kendall"}, and \code{"spearman"}. If one of the latter three is chosen, the distances are computed as \code{as.dist(1 - } 
\cr \code{cor(myinput))}, with the corresponding option for the \code{method} argument of the \code{cor} function.It can also be a character string specifying a user supplied dissimilarity (distance) function for \code{myinput}. See \code{details} and \code{examples} below for further explanation.
}
  \item{rand.fun}{
A character string specifying the permutation method to be applied to \code{myinput}. If NA(default), no permutation is performed. \code{"shuffle.column"} performs a random permutation independently within each column. With \code{"shuffle.block"}, a random permutation is performed independently within each block of columns, as specified by the \code{by.block} argument, and independently from the other blocks. It can also be a character string specifying a user supplied randomization function for \code{myinput}. See \code{details} and \code{examples} below for further explanation.
}
  \item{by.block}{
A vector of the same length as the column dimension of \code{myinput}, to specify the blocking of columns of \code{myinput}. It is used in conjunction with \code{rand.fun}
\cr \code{= "shuffle.block"}, and is ignored otherwise.  
}
  \item{distrib}{
One of \code{"vanilla", "Rparallel"} to specify the distributed computing option for the cluster assignment step. For \code{"vanilla"} (default)
no distributed computing is performed. For \code{"Rparallel"} the \code{parallel} package of \code{R} core is used for multi-core processing.
}
  \item{Ptail}{
Logical. If \code{Ptail} is TRUE(default), the Generalized Pareto Distribution is used to approximate the tail of the null distribution for each of the chosen measures. Otherwise, empirical p-values are computed directly from the corresponding samples.
}
  \item{tailmethod}{
A character string only needed to be specified if the \code{Ptail} is set to TRUE. For \code{"ML"} the parameters of the Generalized Pareto Distribution are estimated by likelihood maximization; for \code{"MOM"} they are estimated by the method of moments.
}
  \item{njobs}{
A single integer specifying the number of worker jobs to create in case of distributed computation if \code{distrib = "Rparallel"}; ignored otherwise.
}
  \item{seed}{
An optional single integer value, to be used to set the random number generator seed (see \code{details}).
}
  \item{Nperm}{
A single integer specifying the size of a sample from the null distribution. See \code{details} for the default sample size.
}
  \item{metric.args}{
Additional arguments for user-supplied dissimilarity (distance) function. See \code{details} and \code{examples} below for further explanation.
}
  \item{rand.args}{
Additional arguments for user-supplied randomization function. See \code{details} and \code{examples} below for further explanation.
}
}
  \details{
When \code{rand.fun} is set to the name of a user supplied randomization function, the first argument of that function should be set to \code{myinput}. See examples below.

The measures of tightness are defined as follows. Denote a node in the tree by \code{a}, its sibling node by \code{b}, and their parent node by \code{p}. Let their respective geights be \code{ha},\code{hb},\code{hp}. Finally, let \code{Sx} mean that the measure \code{S} is computed for the node \code{x}. Then the definitions are

fldc:

\code{Sa = (hp-ha)/hp}

fldcc:

\code{Sa = (hp-(ha-hb)/2)/ha}

bldc:

\code{Sp = (2*hp-ha-hb)/(2*hp)}

slb:

\code{Sp = 2*hp-ha-hb}

The first three measures test tightnss of all internal nodes at the same time, while \code{slb} only tests two-way split of input data. 
The \code{seed} argument is optional. Setting the seed ensures reproducibility of sampling from the null distribution.
}
\value{
If \code{rand.fun} is set to NA, the function returns a matrix whose rows correspond to the internal nodes of the tree and whose columns contain the tree structure as in the \code{merge} component of the class \code{hclust}; the \code{height} component of \code{hclust};and columns tabulating the values of the measures of tightness specified by the \code{mystat} argument.
If \code{rand.fun} is set to a specific randomization method, an object of class \code{\link{best}} is returned. See \code{?best} for details.
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
If \code{mymetric} or \code{rand.fun} is a customized function, make sure you have read and write permission for your working directory.
}
\seealso{
\code{\link{best}},\code{\link{plot.best}}
}
\examples{
####Each column is a gene expression profile for a case of leukemia. 
####Each case belongs to one of three subtypes.
data(leukemia)
#output only statistic table
mytable<-SigTree(data.matrix(leukemia),mystat="all",
        mymethod="ward",mymetric="euclidean")
class(mytable)
\dontrun{
#use multicore processing to detect significant sub-clusters
mytable<-SigTree(data.matrix(leukemia),mystat="all",
	mymethod="ward",mymetric="euclidean",rand.fun="shuffle.column",
	distrib="Rparallel",njobs=2,Ptail=TRUE,tailmethod="ML")
class(mytable)
####Each row after the 1st describes an item belonging to one of four subtypes. 
####Each column corresponds to a genomic location in one of 22 human chromosomes. 
####The 1st row contains the chromosome numbers.
data(T10)
#Perform randomization within each chromosome
chrom<-as.numeric(T10[1,])
mydata<-T10[-1,] 
mytable<-SigTree(data.matrix(mydata),mystat="fldc",        
	mymethod="ward",mymetric="euclidean",rand.fun="shuffle.block",
	by.block=chrom,distrib="Rparallel",njobs=2,Ptail=TRUE,tailmethod="ML")
#Compute dissimilarity using a user-supplied distance function,
#and perform randomization using a user-supplied randomization function, 
#with additional arguments. 
#Both user-supplied functions are only useful as illustration.
mydist<-function(x,y){return(dist(x)/y)}
myrand<-function(x,z){return(apply(x+z,2,sample))}
mytable<-SigTree(data.matrix(leukemia),mystat="fldc",
mymethod="ward",mymetric="mydist",rand.fun="myrand",
distrib="Rparallel",njobs=2,Ptail=TRUE,tailmethod="MOM",metric.args=list(3),
rand.args=list(2))
}
}
