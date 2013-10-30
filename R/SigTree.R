SigTree <-
function(myinput,mystat=c("all","fldc","bldc","fldcc"),mymethod="complete",mymetric="euclidean",
	rand.fun=NA,by.block=NA,distrib=c("vanilla","Rparallel"),Ptail=TRUE,tailmethod=c("ML","MOM"),
	njobs=1,seed=NA,Nperm=ifelse(Ptail,1000,1000*nrow(myinput)),...){
        indextable<-TreeStat(myinput,mystat,method=mymethod,metric=mymetric)
	if(!is.na(rand.fun)){
        if(class(myinput)!="matrix")stop("Inappropriate input data")
	mystatname<-match.arg(mystat, several.ok = TRUE)
        if(any(mystatname=="all"))mystatname<-c("fldc","bldc","fldcc")
	#size<-ceiling((nrow(indextable)+1)/10)*10
	size<-nrow(myinput)
	if(Nperm==ifelse(Ptail,1000,1000*nrow(myinput)))nperm<-size*1000
	if(Nperm!=ifelse(Ptail,1000,1000*nrow(myinput)))nperm<-size*Nperm
        if(nperm%%200!=0) batches<-c(rep(200,floor(nperm/200)),nperm%%200)
        if(nperm%%200==0) batches<-c(rep(200,floor(nperm/200)))
        distrib<-match.arg(distrib)
	if(distrib=="Rparallel"){
                ncores<-min(njobs,length(batches),detectCores())
                cl<-parallel::makeCluster(getOption("cl.cores",ncores))
		if(!is.na(seed)) clusterSetRNGStream(cl,iseed=seed)
		if(rand.fun!="shuffle.column"&rand.fun!="shuffle.block"){
			define.function<-get(rand.fun)
			save(define.function,file=paste(getwd(),"/define.function",sep=""))
                	parallel::clusterEvalQ(cl=cl,expr={
			load(paste(getwd(),"/define.function",sep=""))})
			#parallel::clusterExport(cl,"define.function")
			rand.fun<-"define.function"
		}
		}
	if(Ptail){
                if(Nperm==1000|(Nperm*size)>=10000){
			nperm<-1000
			batches<-10
               		profpack<-vector(mode="list",length=batches)
                	for(pn in 1:batches){
                        	profpack[[pn]]<-vector(mode="list",length=2)
                        	names(profpack[[pn]])<-c("myinput","nperm")
                        	profpack[[pn]]$myinput<-myinput
                        	profpack[[pn]]$nperm<-100
                	}
		}else if((Nperm*size)<10000){
			profpack<-vector(mode="list",length=length(batches))
                	for(pn in 1:length(batches)){
                        	profpack[[pn]]<-vector(mode="list",length=2)
                        	names(profpack[[pn]])<-c("myinput","nperm")
                        	profpack[[pn]]$myinput<-myinput
                        	profpack[[pn]]$nperm<-batches[pn]
                	}
		}
		processed<-switch(distrib,
                        vanilla=lapply(X=profpack,FUN=RandTail,mystat=mystatname,mymethod=mymethod,
                                mymetric=mymetric,rand.fun=rand.fun,by.block=by.block,...),
			Rparallel=parLapply(cl,X=profpack,fun=RandTail,mystat=mystatname,mymethod=mymethod,
                                mymetric=mymetric,rand.fun=rand.fun,by.block=by.block,...))
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
		tailmethod<-match.arg(tailmethod)
		for(statname in mystatname){
			mystat<-vector(mode="list",length=nrow(indextable))
			for(pn in 1:nrow(indextable)){
                        	mystat[[pn]]<-vector(mode="list",length=2)
                        	names(mystat[[pn]])<-c("x","y")
                        	mystat[[pn]]$x<-indextable[pn,statname]
                        	mystat[[pn]]$y<-mytailcounts[[statname]][pn,]
                	}
			jointcounts[,statname]<-unlist(switch(distrib,
				vanilla=lapply(X=mystat,FUN=Pvalue,method=tailmethod,Nexcmax=min(150,nperm/4)),
				Rparallel=parLapply(cl,X=mystat,fun=Pvalue,method=tailmethod,Nexcmax=min(150,nperm/4))))
			jointcounts[jointcounts[,statname]==0,statname]<- .Machine$double.eps
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
                                mymetric=mymetric,rand.fun=rand.fun,by.block=by.block,...),
                        Rparallel=parLapply(cl,X=profpack,fun=RandTree,mystat=mystatname,mymethod=mymethod,
                                mymetric=mymetric,rand.fun=rand.fun,by.block=by.block,...))
        	for(pn in 1:length(processed)){
                	if(pn==1)jointcounts<-processed[[1]][[1]]
               		if(pn>1)jointcounts<-jointcounts+processed[[pn]][[1]]
        	}
        	jointcounts<-(jointcounts+1)/(nperm+2)
		if(distrib=="Rparallel")stopCluster(cl)
	}
        dimnames(jointcounts)[[2]]<-paste("p",dimnames(jointcounts)[[2]],sep="")
	indextable<-cbind(indextable,jointcounts)
	if(rand.fun=="define.function"){system(paste("rm ",getwd(),"/define.function",sep=""))}
	clusterobj<-vector(mode="list",length=3)
	names(clusterobj)<-c("Call","data","indextable")
	clusterobj$Call<-match.call()
	clusterobj$data<-myinput
        clusterobj$indextable<-indextable
	class(clusterobj)<-"best"
        return(clusterobj)
        }
        if(is.na(rand.fun))return(indextable)
}
