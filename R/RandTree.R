RandTree <-
function(mydata,mystat,mymethod,mymetric,rand.fun=c("shuffle.column","shuffle.block",
"define.function"),by.block=NA,...){
        myinput<-mydata$myinput
        ntest<-mydata$nperm
        indextable<-TreeStat(myinput,mystat=mystat,method=mymethod,metric=mymetric)
        sizetable<-matrix(ncol=ntest,nrow=nrow(myinput)-1)
        fldctable<-sizetable
        bldctable<-sizetable
        fldcctable<-sizetable
        #randomization
        for(i in 1:ntest){
                if(rand.fun=="shuffle.column")myrdata<-apply(myinput,2,sample)
		#if(rand.fun=="shuffle.column")myrdata<-apply(myinput,2,function(x) runif(nrow(myinput))*(max(x,na.rm=T)-min(x,na.rm=T))+min(x,na.rm=T))
		#if(rand.fun=="shuffle.column")myrdata<-apply(myinput,2,function(x) runif(nrow(myinput),min(x,na.rm=T),max(x,na.rm=T)))
		if(is.na(by.block[1])&rand.fun=="shuffle.block"){
			stop("by.block needs to be specified")
		}		
                if(rand.fun=="shuffle.block"){
                        myrdata<-t(myinput)
			myrlist<-by(myrdata,by.block,FUN=byfactor)
                        for(j in 1:length(myrlist)){
                             	if(j==1){	myrdata<-myrlist[[j]]
				}else{	myrdata<-rbind(myrdata,myrlist[[j]]) 
				}
                        }
                        myrdata<-t(myrdata)
                }
                if(rand.fun=="define.function"){
                        #myrdata<-t(myinput)
			define.function<-match.fun(define.function)
                        #myrlist<-by(myrdata,by.block,FUN=define.function)
			#for(j in 1:length(myrlist)){
                        #        if(j==1){       myrdata<-myrlist[[j]]
                        #        }else{  myrdata<-rbind(myrdata,myrlist[[j]])
                        #        }
                        #}
                        #myrdata<-t(myrdata)
			myrdata<-define.function(myinput,...)
			if(nrow(myrdata)!=nrow(myinput)|ncol(myrdata)!=ncol(myinput)){
				stop("define.function returns wrong dimension")
			}
                }
                rindextable<-TreeStat(myrdata,mystat=mystat,method=mymethod,metric=mymetric)
                sizetable[,i]<-rindextable[,"clustersize"]
		if(!is.na(match("fldc",mystat)))
                fldctable[,i]<-rindextable[,"fldc"]
                if(!is.na(match("bldc",mystat)))
                bldctable[,i]<-rindextable[,"bldc"]
                if(!is.na(match("fldcc",mystat)))
                fldcctable[,i]<-rindextable[,"fldcc"]
        }
        statnames<-mystat
        sortedsize<-apply(sizetable,2,sort)
        allcounts<-matrix(ncol=length(statnames),nrow=nrow(myinput)-1,data=0)
        colnames(allcounts)<-statnames
        for(statname in statnames){
                mystat<-get(paste(statname,"table",sep=""))
                statmax<-max(mystat)
                randomX<-sizetable+0.5*mystat/statmax
                randomX<-apply(randomX,2,sort)
                #compute counts for p value
		rmatch<-apply(sortedsize,2,bestmatch,size=indextable[,"clustersize"])
        	rmatchl<-apply(sortedsize,2,bestmatchl,size=indextable[,"clustersize"])
        	data<-2*statmax*(randomX-sortedsize)[nrow(randomX)*(col(rmatch)-1)+rmatch]
        	datal<-2*statmax*(randomX-sortedsize)[nrow(randomX)*(col(rmatchl)-1)+rmatchl]
        	mydata<-pmax(data,datal)
		nullstat<-matrix(ncol=ncol(rmatch),data=mydata)	
        	allcounts[,statname]<-allcounts[,statname]+
                	rowSums(nullstat>=indextable[,statname])
        }
        return(list(allcounts,ntest))
}
