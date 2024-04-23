library(caret)
library(ggplot2)  
library(optmatch)

# This R code is a slight modification of code in Prof. Dylan Small's lecture
# It is used to construct rank based Mahalanobis distance with propensity score caliper


optmatch_caliper<-function(datatemp,nocontrols.per.match,ps.formula,mahal.formula,calipersd=.5){
  # Comment about this code and subsequent matching code
  # There is assumed to be no missing data in the variables that go into the 
  # propensity score model so that there is a propensity score for every variable in 
  # the data frame. 
  # Fit a propensity score using logistic regression with each covariate entering 
  # linearly into the logistic link function
  # Put x=TRUE in order to have model object include design matrix
  propscore.model=glm(ps.formula,family=binomial,data=datatemp)
  
  # This model is to obtain model.matrix for mahalanobis distance.
  mahal.model=glm(mahal.formula,family=binomial,x=TRUE,y=TRUE,data=datatemp) 
  
  datatemp$treated=mahal.model$y
  datatemp$treatment=datatemp$treated
  # Use the caret package to include all categories of categorical variables (i.e.,
  # do not leave out one category) in X matrix
  dmy=dummyVars(mahal.model$formula,data=datatemp)
  Xmat=data.frame(predict(dmy,newdata=datatemp))
  # Matrix of covariates to include in the Mahalanobis distance, for now include all 
  # covariates
  Xmatmahal=Xmat
  
  treated=datatemp$treated
  datatemp$logit.ps=predict(propscore.model) 

  
  # Use Hansen (2009)’s rule for removing subjects who lack overlap 
  logit.propscore=datatemp$logit.ps
  pooled.sd.logit.propscore=sqrt(var(logit.propscore[datatemp$treatment==1])/2+var(logit.propscore[datatemp$treatment==0])/2)
  min.treated.logit.propscore=min(logit.propscore[datatemp$treatment==1])
  max.control.logit.propscore=max(logit.propscore[datatemp$treatment==0])
  # How many treated and control subjects lack overlap by Hansen's criterion
  no.treated.lack.overlap=sum(logit.propscore[datatemp$treatment==1]>(max.control.logit.propscore+.5*pooled.sd.logit.propscore))
  no.control.lack.overlap=sum(logit.propscore[datatemp$treatment==0]<(min.treated.logit.propscore-.5*pooled.sd.logit.propscore))
  # If there are subjects who lack overlap, remove them from the datatemp dataset
  datatemp.original=datatemp
  datatemp.full=datatemp
  Xmat.original=Xmat
  Xmat.full=Xmat
  if(no.treated.lack.overlap+no.control.lack.overlap>0){
    which.remove=which((logit.propscore>(max.control.logit.propscore+.5*pooled.sd.logit.propscore))|(logit.propscore<(min.treated.logit.propscore-.5*pooled.sd.logit.propscore)))
    datatemp=datatemp[-which.remove,]
    datatemp.full=rbind(datatemp,datatemp.original[which.remove,])
    Xmat=Xmat[-which.remove,]
    Xmat.full=rbind(Xmat,Xmat.original[which.remove,])
    Xmatmahal=Xmatmahal[-which.remove,]
  }
  # For the purposes of balance checking later, in datatemp.full, append 
  # the removed rows of datatemp to the end of datatemp
  
  # Make the rownames in datatemp be 1:number of rows
  rownames(datatemp)=seq(1,nrow(datatemp),1) 
  
  # Function for computing 
  # rank based Mahalanobis distance.  Prevents an outlier from
  # inflating the variance for a variable, thereby decreasing its importance.
  # Also, the variances are not permitted to decrease as ties 
  # become more common, so that, for example, it is not more important
  # to match on a rare binary variable than on a common binary variable
  # z is a vector, length(z)=n, with z=1 for treated, z=0 for control
  # X is a matrix with n rows containing variables in the distance
  
  smahal=
    function(z,X){
      X<-as.matrix(X)
      n<-dim(X)[1]
      rownames(X)<-1:n
      k<-dim(X)[2]
      m<-sum(z)
      for (j in 1:k) X[,j]<-rank(X[,j])
      cv<-cov(X)
      vuntied<-var(1:n)
      rat<-sqrt(vuntied/diag(cv))
      cv<-diag(rat)%*%cv%*%diag(rat)
      out<-matrix(NA,m,n-m)
      Xc<-X[z==0,]
      Xt<-X[z==1,]
      rownames(out)<-rownames(X)[z==1]
      colnames(out)<-rownames(X)[z==0]
      library(MASS)
      icov<-ginv(cv)
      for (i in 1:m) out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)
      out
    }
  
  # Function for adding a propensity score caliper to a distance matrix dmat
  # calipersd is the caliper in terms of standard deviation of the logit propensity scoe
  addcaliper=function(dmat,z,logitp,calipersd=.5,penalty=1000){
    # Pooled within group standard devation
    sd.logitp=sqrt((sd(logitp[z==1])^2+sd(logitp[z==0])^2)/2)
    adif=abs(outer(logitp[z==1],logitp[z==0],"-"))
    adif=(adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
    dmat=dmat+adif*penalty
    dmat
  }
  
  
  # Rank based Mahalanobis distance
  distmat=smahal(datatemp$treated,Xmatmahal)
  # Add caliper
  distmat=addcaliper(distmat,datatemp$treated,datatemp$logit.ps,calipersd=.5)
  
  
  # Label the rows and columns of the distance matrix by the rownames in datatemp
  rownames(distmat)=rownames(datatemp)[datatemp$treated==1]
  colnames(distmat)=rownames(datatemp)[datatemp$treated==0]
  
  
  matchvec=pairmatch(distmat,controls=nocontrols.per.match,data=datatemp)
  datatemp$matchvec=matchvec
  ## Create a matrix saying which control units each treated unit is matched to
  ## Create vectors of the subject indices of the treatment units ordered by
  ## their matched set and corresponding control unit
  treated.subject.index=rep(0,sum(treated==1))
  matched.control.subject.index.mat=matrix(rep(0,nocontrols.per.match*length(treated.subject.index)),ncol=nocontrols.per.match)
  matchedset.index=substr(matchvec,start=3,stop=10)
  matchedset.index.numeric=as.numeric(matchedset.index)
  
  for(i in 1:length(treated.subject.index)){
    matched.set.temp=which(matchedset.index.numeric==i)
    treated.temp.index=which(datatemp$treated[matched.set.temp]==1)
    treated.subject.index[i]=matched.set.temp[treated.temp.index]
    matched.control.subject.index.mat[i,]=matched.set.temp[-treated.temp.index]
  }
  matched.control.subject.index=matched.control.subject.index.mat
  
  Xmat.without.missing<-Xmat.full
  treatedmat=Xmat.without.missing[datatemp.full$treated==1,];
  # Standardized differences before matching
  controlmat.before=Xmat.without.missing[datatemp.full$treated==0,];
  controlmean.before=apply(controlmat.before,2,mean,na.rm=TRUE);
  
  treatmean=apply(treatedmat,2,mean,na.rm=TRUE);
  treatvar=apply(treatedmat,2,var,na.rm=TRUE);
  controlvar=apply(controlmat.before,2,var,na.rm=TRUE);
  stand.diff.before=(treatmean-controlmean.before)/sqrt((treatvar+controlvar)/2);
  
  treatmat.after=Xmat.without.missing[treated.subject.index,]
  controlmat.after=Xmat.without.missing[matched.control.subject.index,];
  controlmean.after=apply(controlmat.after,2,mean,na.rm=TRUE);
  treatmean.after=apply(treatmat.after,2,mean,na.rm=TRUE)
  stand.diff.after=(treatmean-controlmean.after)/sqrt((treatvar+controlvar)/2)
  
  res.stand.diff<-cbind(stand.diff.before,stand.diff.after)
  res.mean<-cbind(treatmean.after,controlmean.before,controlmean.after)
  print(round(res.stand.diff,2))
  print(round(res.mean,2))
  
  abs.stand.diff.before=stand.diff.before[-1]
  abs.stand.diff.after=stand.diff.after[-1]
  covariates=names(stand.diff.before[-1])
  plot.dataframe=data.frame(abs.stand.diff=c(abs.stand.diff.before,abs.stand.diff.after),covariates=rep(covariates,2),type=c(rep("Before",length(covariates)),rep("After",length(covariates))))
  p<-ggplot(plot.dataframe,aes(x=abs.stand.diff,y=covariates))+geom_point(size=2,aes(shape=type))+scale_shape_manual(values=c(4,1))+geom_vline(xintercept=c(-.1,.1),lty=2)+xlab("standardized differences in means")+ ylab("")
  return(list(p=p,datatemp=datatemp,treated.subject.index=treated.subject.index,
              matched.control.subject.index=matched.control.subject.index,
              res.stand.diff=res.stand.diff,res.mean=res.mean))
}
