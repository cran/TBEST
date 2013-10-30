plot.best<-function(x,mystat="fldc",sigp=0.05,partition=NA,print.num=TRUE,print.lab=TRUE,float=0.01,
	col.best=c(2,3),cex.best=0.8,font.best=NULL,main=NULL,sub=NULL,xlab=NULL,...){
    	clusterobj<-x
	if (is.null(main)) 
        	main = "Dendrogram with P-values"
	method<-clusterobj$method
	metric<-clusterobj$metric
    	if (is.null(sub)) 
        	sub = paste("Cluster method: ",clusterobj$method)
    	if (is.null(xlab)) 
        	xlab = paste("Distance: ", clusterobj$metric)
	myinput<-clusterobj$data
	sigp<-1-(1-sigp)^(1/nrow(myinput))
	if(data.class(myinput)=="dist")hc<-hclust(myinput,method=method)
        if(data.class(myinput)=="matrix"){
                if(metric!="pearson"&metric!="kendall"&metric!="spearman"){
                        hc<-hclust(dist(myinput,method=clusterobj$metric),method=clusterobj$method)
                }
                if(metric=="pearson"|metric=="kendall"|metric=="spearman"){
                        hc<-hclust(as.dist(1-cor(t(myinput),method=clusterobj$metric,
                                use="pairwise.complete.obs")),method=clusterobj$method)
                }
        }
        if(data.class(myinput)=="hclust")hc<-myinput
	if(print.lab){plot(hc, labels=hc$labels, main = main, sub = sub, xlab = xlab)}
	if(!print.lab){plot(hc, labels=print.lab, main = main, sub = sub, xlab = xlab)}
	pval<-clusterobj$indextable[,paste("p",mystat,sep="")] 
	if(length(partition)==1){
		pval[pval>sigp]<-NA
        	treetext(hc,pval, col = col.best, cex = cex.best, float = float, font=font.best,print.num = print.num)
	}
	if(length(partition)>1){
		if(length(partition$partition)==1){
		pval[pval>sigp]<-NA
                treetext(hc,pval, col = col.best, cex = cex.best, float = float, font=font.best,print.num = print.num)
		}else{
		partitionp<-rep(NA,length(pval))
		partitionp[unique(partition$partition[,2])]<-pval[unique(partition$partition[,2])]
		treetext(hc,partitionp, col = col.best, cex = cex.best, float = float, font=font.best,print.num = print.num)
		}
        }	
}
