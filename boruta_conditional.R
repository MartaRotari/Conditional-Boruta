
Condition_Boruta<-function(x,...)
  UseMethod("Boruta")


Boruta.default<-function(x,y,pValue=0.01,mcAdj=TRUE,maxRuns=100,doTrace=0,holdHistory=TRUE,getImp=getImpRfCond,...){
  #Timer starts... now!
  timeStart<-Sys.time()
  
  #Extract the call to store in output
  cl<-match.call()
  cl[[1]]<-as.name('Boruta')
  
  #Convert x into a data.frame
  if(!is.data.frame(x))
    x<-data.frame(x)
  
  ##Some checks on x & y
  if(length(grep('^shadow',names(x)))>0)
    stop('Attributes with names starting from "shadow" are reserved for internal use. Please rename them.')
  if(any(c(is.na(x),is.na(y))))
    stop('Cannot process NAs in input. Please remove them.')
  if(maxRuns<11)
    stop('maxRuns must be greater than 10.')
  
  ##Expands the information system with newly built random attributes and calculates importance
  addShadowsAndGetImp<-function(decReg,runs){
    #xSha is going to be a data frame with shadow attributes; time to init it.
    xSha<-x[,decReg!="Rejected",drop=F]
    while(dim(xSha)[2]<5) xSha<-cbind(xSha,xSha); #There must be at least 5 random attributes.
    
    #Now, we permute values in each attribute
    nSha<-ncol(xSha)
    data.frame(lapply(xSha,sample))->xSha
    names(xSha)<-paste('shadow',1:nSha,sep="")
    
    #Notifying user of our progress
    if(doTrace>1)
      message(sprintf(' %s. run of importance source...',runs))
    
    #Calling importance source; "..." can be used by the user to pass rf attributes (for instance ntree)
    impRaw<-getImp(cbind(x[,decReg!="Rejected"],xSha),y,...)
    if(!is.numeric(impRaw))
      stop("getImp result is not a numeric vector. Please check the given getImp function.")
    if(length(impRaw)!=sum(decReg!="Rejected")+ncol(xSha))
      stop("getImp result has a wrong length. Please check the given getImp function.")
    if(any(is.na(impRaw)|is.nan(impRaw))){
      impRaw[is.na(impRaw)|is.nan(impRaw)]<-0
      warning("getImp result contains NA(s) or NaN(s); replacing with 0(s), yet this is suspicious.")
    }
    
    #Importance must have Rejected attributes put on place and filled with -Infs
    imp<-rep(-Inf,nAtt+nSha);names(imp)<-c(attNames,names(xSha))
    impRaw->imp[c(decReg!="Rejected",rep(TRUE,nSha))]
    shaImp<-imp[(nAtt+1):length(imp)];imp[1:nAtt]->imp
    
    return(list(imp=imp,shaImp=shaImp))
  }
  
  ##Assigns hits
  assignHits<-function(hitReg,curImp){
    curImp$imp>max(curImp$shaImp)->hits
    if(doTrace>2){
      uncMask<-decReg=="Tentative"
      intHits<-sum(hits[uncMask])
      if(intHits>0)
        message(sprintf("Assigned hit to %s attribute%s out of %s undecided.",sum(hits[uncMask]),if(intHits==1) "" else "s",sum(uncMask)))
      else
        message("None of undecided attributes scored a hit.")
    }
    hitReg[hits]<-hitReg[hits]+1
    return(hitReg)
  }
  
  ##Checks whether number of hits is significant
  doTests<-function(decReg,hitReg,runs){
    pAdjMethod<-ifelse(mcAdj[1],'bonferroni','none')
    #If attribute is significantly more frequent better than shadowMax, its claimed Confirmed
    toAccept<-stats::p.adjust(stats::pbinom(hitReg-1,runs,0.5,lower.tail=FALSE),method=pAdjMethod)<pValue
    (decReg=="Tentative" & toAccept)->toAccept
    
    #If attribute is significantly more frequent worse than shadowMax, its claimed Rejected (=irrelevant)
    toReject<-stats::p.adjust(stats::pbinom(hitReg,runs,0.5,lower.tail=TRUE),method=pAdjMethod)<pValue
    (decReg=="Tentative" & toReject)->toReject
    
    #Update decReg
    decReg[toAccept]<-"Confirmed";"Rejected"->decReg[toReject]
    
    #Report progress
    if(doTrace>0){
      nAcc<-sum(toAccept)
      nRej<-sum(toReject)
      nLeft<-sum(decReg=="Tentative")
      if(nAcc+nRej>0)
        message(sprintf("After %s iterations, +%s: ",runs,format(difftime(Sys.time(),timeStart),digits=2)))
      if(nAcc>0)
        message(sprintf(" confirmed %s attribute%s: %s",
                        nAcc,ifelse(nAcc==1,'','s'),.attListPrettyPrint(attNames[toAccept])))
      if(nRej>0)
        message(sprintf(" rejected %s attribute%s: %s",
                        nRej,ifelse(nRej==1,'','s'),.attListPrettyPrint(attNames[toReject])))
      if(nAcc+nRej>0)
        if(nLeft>0){
          message(sprintf(" still have %s attribute%s left.\n",
                          nLeft,ifelse(nLeft==1,'','s')))
        }else{
          if(nAcc+nRej>0) message(" no more attributes left.\n")
        }
    }
    return(decReg)
  }
  
  ##Creating some useful constants
  nAtt<-ncol(x); nrow(x)->nObjects
  attNames<-names(x); c("Tentative","Confirmed","Rejected")->confLevels
  
  ##Initiate state
  decReg<-factor(rep("Tentative",nAtt),levels=confLevels)
  hitReg<-rep(0,nAtt);names(hitReg)<-attNames
  impHistory<-list()
  runs<-0
  
  ##Main loop
  
  while(any(decReg=="Tentative") && (runs+1->runs)<maxRuns){
    curImp<-addShadowsAndGetImp(decReg,runs)
    hitReg<-assignHits(hitReg,curImp)
    decReg<-doTests(decReg,hitReg,runs)
    
    #If needed, update impHistory with scores obtained in this iteration
    if(holdHistory){
      imp<-c(curImp$imp,
             shadowMax=max(curImp$shaImp),
             shadowMean=mean(curImp$shaImp),
             shadowMin=min(curImp$shaImp))
      impHistory<-c(impHistory,list(imp))
    }
  }
  
  ##Building result
  impHistory<-do.call(rbind,impHistory)
  names(decReg)<-attNames
  ans<-list(finalDecision=decReg,ImpHistory=impHistory,
            pValue=pValue,maxRuns=maxRuns,light=TRUE,mcAdj=mcAdj,
            timeTaken=Sys.time()-timeStart,roughfixed=FALSE,call=cl,
            impSource=comment(getImp))
  
  "Boruta"->class(ans)
  return(ans)
}

.attListPrettyPrint<-function(x,limit=5){
  x<-sort(x)
  if(length(x)<limit+1)
    return(sprintf("%s;",paste(x,collapse=", ")))
  sprintf("%s and %s more;",paste(utils::head(x,limit),collapse=", "),length(x)-limit)
}


Boruta.formula<-function(formula,data=.GlobalEnv,...){
  ##Grab and interpret the formula
  stats::terms.formula(formula,data=data)->t
  x<-eval(attr(t,"variables"),data)
  apply(attr(t,"factors"),1,sum)>0->sel
  nam<-rownames(attr(t,"factors"))[sel]
  data.frame(x[sel])->df;names(df)<-nam
  x[[attr(t,"response")]]->dec
  
  ##Run Boruta
  ans<-Boruta.default(df,dec,...)
  ans$call<-match.call()
  ans$call[[1]]<-as.name('Boruta')
  formula->ans$call[["formula"]]
  return(ans)
}


print.Boruta<-function(x,...){
  if(class(x)!='Boruta') stop("This is NOT a Boruta object!")
  cat(paste('Boruta performed ',dim(x$ImpHistory)[1],' iterations in ',format(x$timeTaken),'.\n',sep=''))
  if(x$roughfixed) cat(paste('Tentatives roughfixed over the last ',x$averageOver,' iterations.\n',sep=''))
  if(sum(x$finalDecision=='Confirmed')==0){
    cat(' No attributes deemed important.\n')} else {
      writeLines(strwrap(paste(sum(x$finalDecision=='Confirmed'),' attributes confirmed important: ',
                               .attListPrettyPrint(names(x$finalDecision[x$finalDecision=='Confirmed']))),indent=1))
    }
  if(sum(x$finalDecision=='Rejected')==0){
    cat(' No attributes deemed unimportant.\n')} else {
      writeLines(strwrap(paste(sum(x$finalDecision=='Rejected'),' attributes confirmed unimportant: ',
                               .attListPrettyPrint(names(x$finalDecision[x$finalDecision=='Rejected']))),indent=1))
    }
  if(sum(x$finalDecision=='Tentative')!=0){
    writeLines(strwrap(paste(sum(x$finalDecision=='Tentative'),' tentative attributes left: ',
                             .attListPrettyPrint(names(x$finalDecision[x$finalDecision=='Tentative']))),indent=1))
  }
  invisible(x)
}



getImpLegacyRfcond<-function(x,y,...){
  rff <- party::cforest(y~., data = x, control = cforest_unbiased(mtry = 20, ntree = 500))
  permimp::permimp(rff, conditional = TRUE, progressBar = FALSE)$values
}
comment(getImpLegacyRfcond)<-'randomForest conditional permutation importance'



getImpRfCond<-function(x,y,ntree=500,num.trees=ntree,...){
  x$shadow.Boruta.decision<-y
  rff <- party::cforest(shadow.Boruta.decision~., data = x, control = party::cforest_unbiased(mtry = 20, ntree = num.trees))
  permimp::permimp(rff, conditional = TRUE, progressBar = FALSE)$values
}
comment(getImpRfCond)<-'ranger conditional permutation importance'


getImpLegacyRfZ<-function(x,y,...){
  randomForest::randomForest(x,y,
                             importance=TRUE,keep.forest=FALSE,...)->rf
  randomForest::importance(rf,1,scale=TRUE)[,1]
}
comment(getImpLegacyRfZ)<-'randomForest normalized permutation importance'

#' @rdname getImpLegacyRf
#' @export
getImpLegacyRfRaw<-function(x,y,...){
  randomForest::randomForest(x,y,
                             importance=TRUE,keep.forest=FALSE,...)->rf
  randomForest::importance(rf,1,scale=FALSE)[,1]
}
comment(getImpLegacyRfRaw)<-'randomForest raw permutation importance'


getImpLegacyRfGini<-function(x,y,...){
  randomForest::randomForest(x,y,
                             keep.forest=FALSE,...)->rf
  randomForest::importance(rf,2,scale=FALSE)[,1]
}
comment(getImpLegacyRfGini)<-'randomForest Gini index importance'


getImpRfZ<-function(x,y,ntree=500,num.trees=ntree,...){
  if(inherits(y,"Surv")){
    x$shadow.Boruta.time<-y[,"time"]
    x$shadow.Boruta.status<-y[,"status"]
    return(ranger::ranger(data=x,
                          dependent.variable.name="shadow.Boruta.time",
                          status.variable.name="shadow.Boruta.status",
                          num.trees=num.trees,importance="permutation",
                          scale.permutation.importance=TRUE,
                          write.forest=FALSE,...)$variable.importance)
  }
  #Abusing the fact that Boruta disallows attributes with names
  # starting from "shadow"
  x$shadow.Boruta.decision<-y
  ranger::ranger(data=x,dependent.variable.name="shadow.Boruta.decision",
                 num.trees=num.trees,importance="permutation",
                 scale.permutation.importance=TRUE,
                 write.forest=FALSE,...)$variable.importance
}
comment(getImpRfZ)<-'ranger normalized permutation importance'

getImpRfGini<-function(x,y,ntree=500,num.trees=ntree,...){
  if(inherits(y,"Surv"))
    stop("Ranger cannot produce Gini importance for survival problems.")
  x$shadow.Boruta.decision<-y
  ranger::ranger(data=x,dependent.variable.name="shadow.Boruta.decision",
                 num.trees=num.trees,importance="impurity",
                 scale.permutation.importance=FALSE,
                 write.forest=FALSE,...)$variable.importance
}
comment(getImpRfGini)<-'ranger Gini index importance'

getImpRfRaw<-function(x,y,ntree=500,num.trees=ntree,...){
  if(inherits(y,"Surv")){
    x$shadow.Boruta.time<-y[,"time"]
    x$shadow.Boruta.status<-y[,"status"]
    return(ranger::ranger(data=x,
                          dependent.variable.name="shadow.Boruta.time",
                          status.variable.name="shadow.Boruta.status",
                          num.trees=num.trees,importance="permutation",
                          write.forest=FALSE,...)$variable.importance)
  }
  x$shadow.Boruta.decision<-y
  ranger::ranger(data=x,dependent.variable.name="shadow.Boruta.decision",
                 num.trees=num.trees,importance="permutation",
                 scale.permutation.importance=FALSE,
                 write.forest=FALSE,...)$variable.importance
}
comment(getImpRfRaw)<-'ranger raw permutation importance'


getImpExtraZ<-function(x,y,ntree=500,num.trees=ntree,...)
  getImpRfZ(x,y,ntree=ntree,splitrule="extratrees",...)
comment(getImpExtraZ)<-'ranger normalized permutation importance'


getImpExtraGini<-function(x,y,ntree=500,num.trees=ntree,...)
  getImpRfGini(x,y,ntree=ntree,splitrule="extratrees",...)
comment(getImpExtraGini)<-'ranger extra-trees Gini index importance'

getImpFerns<-function(x,y,...){
  f<-rFerns::rFerns(x,y,
                    saveForest=FALSE,importance=TRUE,...)
  f$importance[,1]
}
comment(getImpFerns)<-'rFerns importance'



attStats<-function(x){
  if(class(x)!='Boruta')
    stop('This function needs Boruta object as an argument.')
  if(is.null(x$ImpHistory))
    stop('Importance history was not stored during the Boruta run.')
  lz<-lapply(1:ncol(x$ImpHistory),function(i) x$ImpHistory[is.finite(x$ImpHistory[,i]),i])
  colnames(x$ImpHistory)->names(lz)
  mr<-lz$shadowMax; lz[1:(length(lz)-3)]->lz
  t(sapply(lz,function(x) c(mean(x),stats::median(x),min(x),max(x),sum(mr[1:length(x)]<x)/length(mr))))->st
  st<-data.frame(st,x$finalDecision)
  names(st)<-c("meanImp","medianImp","minImp","maxImp","normHits","decision")
  return(st)
}


getSelectedAttributes<-function(x,withTentative=FALSE){
  if(class(x)!='Boruta') stop('This function needs Boruta object as an argument.')
  names(x$finalDecision)[
    x$finalDecision%in%(if(!withTentative) "Confirmed" else c("Confirmed","Tentative"))
  ]
}


TentativeRoughFix<-function(x,averageOver=Inf){
  if(!inherits(x,'Boruta'))
    stop('This function needs Boruta object as an argument.')
  if(is.null(x$ImpHistory))
    stop('Importance history was not stored during the Boruta run.')
  if(!is.numeric(averageOver))
    stop('averageOver should be a numeric vector.')
  if(length(averageOver)!=1)
    stop('averageOver should be a one-element vector.')
  if(averageOver<1)
    stop('averageOver should be positive.')
  
  tentIdx<-which(x$finalDecision=='Tentative')
  if(length(tentIdx)==0){
    warning('There are no Tentative attributes! Returning original object.')
    return(x)
  }
  
  nRuns<-dim(x$ImpHistory)[1]
  
  if(averageOver>nRuns)
    averageOver<-nRuns
  
  impHistorySubset<-x$ImpHistory[(nRuns-averageOver+1):nRuns,]
  medianTentImp<-sapply(impHistorySubset[,tentIdx],stats::median)
  medianShaMaxImp<-stats::median(impHistorySubset[,'shadowMax'])
  medianTentImp>medianShaMaxImp->toOrdain
  
  ans<-x
  ans$roughfixed<-TRUE
  ans$averageOver<-averageOver
  ans$originalDecision<-x$finalDecision
  ans$finalDecision[tentIdx[toOrdain]]<-'Confirmed'
  ans$finalDecision[tentIdx[!toOrdain]]<-'Rejected'
  
  return(ans)
}

##generateCol is internally used by plot.Boruta and plotImpHistory
generateCol<-function(x,colCode,col,numShadow){
  #Checking arguments
  if(is.null(col) & length(colCode)!=4)
    stop('colCode should have 4 elements.')
  #Generating col
  if(is.null(col)){
    rep(colCode[4],length(x$finalDecision)+numShadow)->cc
    cc[c(x$finalDecision=='Confirmed',rep(FALSE,numShadow))]<-colCode[1]
    cc[c(x$finalDecision=='Tentative',rep(FALSE,numShadow))]<-colCode[2]
    cc[c(x$finalDecision=='Rejected',rep(FALSE,numShadow))]<-colCode[3]
    col=cc
  }
  return(col)
}


plot.Boruta<-function(x,colCode=c('green','yellow','red','blue'),sort=TRUE,whichShadow=c(TRUE,TRUE,TRUE),
                      col=NULL,xlab='Attributes',ylab='Importance',...){
  #Checking arguments
  if(class(x)!='Boruta')
    stop('This function needs Boruta object as an argument.')
  if(is.null(x$ImpHistory))
    stop('Importance history was not stored during the Boruta run.')
  
  #Removal of -Infs and conversion to a list
  lz<-lapply(1:ncol(x$ImpHistory),function(i) x$ImpHistory[is.finite(x$ImpHistory[,i]),i])
  colnames(x$ImpHistory)->names(lz)
  
  #Selection of shadow meta-attributes
  numShadow<-sum(whichShadow)
  lz[c(rep(TRUE,length(x$finalDecision)),whichShadow)]->lz
  
  #Generating color vector
  col<-generateCol(x,colCode,col,numShadow)
  
  #Ordering boxes due to attribute median importance
  if(sort){
    ii<-order(sapply(lz,stats::median))
    lz[ii]->lz; col<-col[ii]
  }
  
  Labels <- sort(sapply(lz,median))
  
  #Final plotting
  graphics::boxplot(lz,xlab = "", xaxt = "n",ylab=ylab,col=col,...)
  graphics::axis(side = 1,las=2,labels = names(Labels),at = 1:ncol(cb$ImpHistory), cex.axis = 0.7)
  invisible(x)
}


plotImpHistory<-function(x,colCode=c('green','yellow','red','blue'),col=NULL,type="l",lty=1,pch=0,
                         xlab='Classifier run',ylab='Importance',...){
  #Checking arguments
  if(class(x)!='Boruta')
    stop('This function needs Boruta object as an argument.')
  if(is.null(x$ImpHistory))
    stop('Importance history was not stored during the Boruta run.')
  col<-generateCol(x,colCode,col,3)
  
  #Final plotting
  graphics::matplot(0:(nrow(x$ImpHistory)-1),x$ImpHistory,xlab=xlab,ylab=ylab,col=col,type=type,lty=lty,pch=pch,...)
  invisible(x)
}


getConfirmedFormula<-function(x){
  if(!inherits(x,'Boruta'))
    stop('This function needs Boruta object as an argument.')
  if(is.null(x$call[["formula"]]))
    stop('The model for this Boruta run was not a formula.')
  deparse(x$call[["formula"]][[2]])->dec
  preds<-paste(names(x$finalDecision)[x$finalDecision=='Confirmed'],collapse="+")
  return(stats::as.formula(sprintf('%s~%s',dec,preds)))
}

getNonRejectedFormula<-function(x){
  if(!inherits(x,'Boruta'))
    stop('This function needs Boruta object as an argument.')
  if(is.null(x$call[["formula"]]))
    stop('The model for this Boruta run was not a formula.')
  deparse(x$call[["formula"]][[2]])->dec
  preds<-paste(names(x$finalDecision)[x$finalDecision!='Rejected'],collapse="+")
  return(stats::as.formula(sprintf('%s~%s',dec,preds)))
}


