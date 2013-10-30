SigTree <-
function(myinput,mystat=c("all","fldc","bldc","fldcc"),mymethod="complete",mymetric="euclidean",
	rand.fun=NA,by.block=NA,distrib=c("vanilla","Rparallel"),Ptail=TRUE,tailmethod="ML",njobs=1){
        indextable<-TreeStat(myinput,mystat=mystat,method=mymethod,metric=mymetric)
	if(!is.na(rand.fun)){
        if(class(myinput)!="matrix")stop("Inappropriate input data")
	statnames<-c("fldc","bldc","fldcc")
        if(any(mystat=="all"))mystatname<-statnames
	if(!any(mystat=="all")){
		m<-match(mystat,statnames)
		if(length(which(is.na(m)))!=0)stop("Inappropriate statistic name")
		mystatname<-statnames[m]
	}
	size<-ceiling((nrow(indextable)+1)/10)*10
	nperm<-size*1000
	sigp<-1/nperm
        if(nperm%%200!=0) batches<-c(rep(200,floor(nperm/200)),nperm%%200)
        if(nperm%%200==0) batches<-c(rep(200,floor(nperm/200)))
        distrib<-match.arg(distrib)
	if(distrib=="Rparallel"){
                ncores<-min(njobs,length(batches),detectCores())
                cl<-parallel::makeCluster(getOption("cl.cores",ncores))
		if(rand.fun!="shuffle.column"&rand.fun!="shuffle.block"){
			define.function<-get(rand.fun)
			save(define.function,file=paste(getwd(),"/define.function",sep=""))
                	parallel::clusterEvalQ(cl=cl,expr={
			load(paste(getwd(),"/define.function",sep=""))})
			rand.fun<-"define.function"
		}
		}
	if(Ptail){
                nperm<-1000
		batches<-10
                profpack<-vector(mode="list",length=batches)
                for(pn in 1:batches){
                        profpack[[pn]]<-vector(mode="list",length=2)
                        names(profpack[[pn]])<-c("myinput","nperm")
                        profpack[[pn]]$myinput<-myinput
                        profpack[[pn]]$nperm<-100
                }
		processed<-switch(distrib,
                        vanilla=lapply(X=profpack,FUN=RandTail,mystat=mystatname,mymethod=mymethod,
                                mymetric=mymetric,rand.fun=rand.fun,by.block=by.block),
			Rparallel=parLapply(cl,X=profpack,fun=RandTail,mystat=mystatname,mymethod=mymethod,
                                mymetric=mymetric,rand.fun=rand.fun,by.block=by.block))
		#if(distrib=="Rparallel")stopCluster(cl)
                for(pn in 1:length(processed)){
                        if(pn==1)mytailcounts<-processed[[1]]
                        if(pn>1){
			for(i in 1:length(mystatname)){
			mytailcounts[[i]]<-cbind(mytailcounts[[i]],processed[[pn]][[i]])
			}
			}
                }
		names(mytailcounts)<-mystatname
		jointcounts<-matrix(nrow=nrow(indextable),ncol=length(mystatname),
			dimnames=list(c(1:nrow(indextable)),mystatname))
		for(statname in mystatname){
			mystat<-vector(mode="list",length=nrow(indextable))
			for(pn in 1:nrow(indextable)){
                        	mystat[[pn]]<-vector(mode="list",length=2)
                        	names(mystat[[pn]])<-c("x","y")
                        	mystat[[pn]]$x<-indextable[pn,statname]
                        	mystat[[pn]]$y<-mytailcounts[[statname]][pn,]
                	}
			mypval<-switch(distrib,
                        vanilla=lapply(X=mystat,FUN=Pvalue,method=tailmethod,Nexcmax=125),
                        Rparallel=parLapply(cl,X=mystat,fun=Pvalue,method=tailmethod,Nexcmax=125))
			#for(i in 1:nrow(indextable)){
			#	myp<-Pvalue(x=indextable[i,statname],y=mytailcounts[[statname]][i,],method=tailmethod,Nexcmax=125)
			#	jointcounts[i,statname]<-myp[[1]]
			#}
			jointcounts[,statname]<-unlist(mypval)[((1:nrow(indextable))-1)*4+1]
			jointcounts[jointcounts[,statname]==0,statname]<-.Machine$double.eps
		}
		if(distrib=="Rparallel")stopCluster(cl)
        }
	if(!Ptail){
        	profpack<-vector(mode="list",length=length(batches))
        	for(pn in 1:length(batches)){
                	profpack[[pn]]<-vector(mode="list",length=2)
                	names(profpack[[pn]])<-c("myinput","nperm")
                	profpack[[pn]]$myinput<-myinput
                	profpack[[pn]]$nperm<-batches[pn]
        	}
        	processed<-switch(distrib,
                        vanilla=lapply(X=profpack,FUN=RandTree,mystat=mystatname,mymethod=mymethod,
                                mymetric=mymetric,rand.fun=rand.fun,by.block=by.block),
                        Rparallel=parLapply(cl,X=profpack,fun=RandTree,mystat=mystatname,mymethod=mymethod,
                                mymetric=mymetric,rand.fun=rand.fun,by.block=by.block))
        	if(distrib=="Rparallel")stopCluster(cl)
        	for(pn in 1:length(processed)){
                	if(pn==1)jointcounts<-processed[[1]][[1]]
               		if(pn>1)jointcounts<-jointcounts+processed[[pn]][[1]]
        	}
        	jointcounts<-(jointcounts+1)/(nperm+2)
	}
        dimnames(jointcounts)[[2]]<-paste("p",dimnames(jointcounts)[[2]],sep="")
	indextable<-cbind(indextable,jointcounts)
	if(rand.fun=="define.function"){system(paste("rm ",getwd(),"/define.function",sep=""))}
	clusterobj<-vector(mode="list",length=4)
	names(clusterobj)<-c("data","method","metric","indextable")
	clusterobj$data<-myinput
	clusterobj$method<-mymethod
	clusterobj$metric<-mymetric
        clusterobj$indextable<-indextable
	class(clusterobj)<-"best"
        return(clusterobj)
        }
        if(is.na(rand.fun))return(indextable)
}
