FindMembership<-function(myinput,mynode,method,metric){
	indextable<-TreeStat(myinput,method=method,metric=metric)
	if(class(mynode)!="numeric"&class(mynode)!="integer")stop("Inappropriate node number")
        if(max(mynode)>nrow(indextable))stop("Node number should be <= Sample Size - 1")
        if(min(mynode)< -(nrow(indextable)+1))stop("Node number should not be < -Sample Size")
	if(length(which(mynode==0))!=0)stop("Node number can not be 0")
	sigpartition<-mynode
        singleton<--sigpartition[sigpartition<0]
        nodes<-sigpartition[sigpartition>0]
        membership<-vector("list",length(nodes))
	clusters<-vector("list",length(mynode))
	names(clusters)<-paste("branch",mynode)
	pos<-which(sigpartition>0)
	npos<-which(sigpartition<0)
	if(length(nodes)>=1){
        for(i in 1:length(nodes)){
                                myfamily<-nodes[i]
                                nmem<-0
                                while(nmem<indextable[nodes[i],"clustersize"]){
                                        membership[[i]]<-unique(c(myfamily,indextable[myfamily,"index1"],
                                                indextable[myfamily,"index2"]))
                                        myfamily<-membership[[i]][membership[[i]]>0]
                                        nmem<-sum(membership[[i]]<0)
                                }
                        membership[[i]]<-(-membership[[i]][membership[[i]]<0])
                        }
	for(i in 1:length(membership)){
        	clusters[[pos[i]]]<-dimnames(myinput)[[1]][membership[[i]]]
        }
	}
	if(length(npos)!=0){
	for(i in 1:length(npos)){
		clusters[[npos[[i]]]]<-dimnames(myinput)[[1]][singleton[i]]
	}
	}
	return(clusters)
}
