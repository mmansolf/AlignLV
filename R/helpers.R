#' Prepare \code{mirt} estimates for alignment
#'
#' Not generally intended to be used on its own, but exported anyway for didactic purposes.
#'
#' See example for \code{\link{Alignment}} for examples
#'
#' This may differ from what is used in Mplus.
#'
#' @export
getEstimates.mirt=function(fit,SE=T,bifactor.marginal=F){
  #coefficients
  coef.raw=coef(fit,printSE=SE)
  if('lr.betas'%in%names(coef.raw))
    coef.raw=coef.raw[which(names(coef.raw)!='lr.betas')]
  #convert to vector if appropriate; just get item parameters
  if(SE) se.raw=coef.raw%>%map(function(x)if(is.matrix(x))x['SE',] else NULL)
  coef.raw=coef.raw%>%map(function(x)if(is.matrix(x))x['par',] else x)
  #item parameters
  coef.raw=coef.raw[fit@Data$tabdata%>%colnames]
  coef.raw=coef.raw[!is.na(names(coef.raw))]
  if(SE) se.raw=se.raw[fit@Data$tabdata%>%colnames]
  if(SE) se.raw=se.raw[!is.na(names(se.raw))]

  # coef.raw=map(coef.raw,~.[,which(.!=0)])

  #loadings
  a.raw=NULL
  if(SE) se.a.raw=NULL
  for(j in 1:(length(coef.raw))){
    a.raw.j=coef.raw[[j]][which(substr(names(coef.raw[[j]]),1,1)=='a')]
    if(SE) se.a.raw.j=se.raw[[j]][which(substr(names(se.raw[[j]]),1,1)=='a')]
    if(length(a.raw.j)>1 & bifactor.marginal){
      c1=(1/1.7)/sqrt(1+sum((a.raw.j/1.7)^2))
      a.raw.j=a.raw.j*c1 #compare to summary(fit)
      a.sigmasq.i=1-a.raw.j^2
      c2=(1/sqrt(a.sigmasq.i)*1.7)
      a.raw.j=(a.raw.j*c2)[1]
      if(SE)se.a.raw.j=se.a.raw.j[1]*(c1*c2)[1]
    }
    if(j==1){
      a.raw=matrix(a.raw.j)
      if(SE) se.a.raw=matrix(se.a.raw.j)
    } else {
      a.raw=rbind(a.raw,a.raw.j)
      if(SE) se.a.raw=rbind(se.a.raw,se.a.raw.j)
    }
    rownames(a.raw)[j]=names(coef.raw)[j]
    if(SE) rownames(se.a.raw)[j]=names(se.raw)[j]
  }

  #maximum number of categories maxcats
  maxcats=10
  d.raw=matrix(0,1,maxcats)
  if(SE) se.d.raw=matrix(0,1,maxcats)
  for(j in 1:(length(coef.raw))){
    dindices=which(substr(names(coef.raw[[j]]),1,1)=="d")
    #print(length(dindices))
    d.raw=rbind(d.raw,c(coef.raw[[j]][dindices],rep(0,maxcats-length(dindices))))
    if(SE) se.d.raw=rbind(se.d.raw,c(se.raw[[j]][dindices],rep(0,maxcats-length(dindices))))
    rownames(d.raw)[j+1]=names(d.raw)[j]
    if(SE) rownames(se.d.raw)[j+1]=names(se.raw)[j]
  }
  d.raw=d.raw[-1,]
  if(SE) se.d.raw=se.d.raw[-1,]
  d.raw=as.matrix(d.raw[,-which(apply(abs(d.raw),2,sum)==0)])
  if(SE) se.d.raw=as.matrix(se.d.raw[,-which(apply(abs(se.d.raw),2,sum)==0)])
  rownames(d.raw)=names(coef.raw)
  if(SE) rownames(se.d.raw)=names(se.raw)
  if(is.null(colnames(d.raw))) colnames(d.raw)='d.1'
  if(SE) if(is.null(colnames(se.d.raw))) colnames(se.d.raw)='d.1'

  #parameter matrix for alignment input
  pars.raw=data.frame(itemname=rownames(a.raw),a.raw,d.raw)
  names(pars.raw)=c('itemname',paste0("a.",1:ncol(a.raw)),
                    paste0("d.",1:ncol(d.raw)))
  if(SE) ses.raw=data.frame(itemname=rownames(se.a.raw),se.a.raw,se.d.raw)
  if(SE) names(ses.raw)=c('itemname',paste0("a.",1:ncol(se.a.raw)),
                          paste0("d.",1:ncol(se.d.raw)))

  #get model parameters, ready to align
  a.toalign=as.matrix(pars.raw[,c(paste0("a.1"))])
  row.names(a.toalign)=pars.raw[,1]
  d.toalign=as.matrix(pars.raw[,which(substr(names(pars.raw),1,2)=="d.")])
  if(SE) se.a.toalign=as.matrix(ses.raw[,c(paste0("a.1"))])
  if(SE) row.names(se.a.toalign)=ses.raw[,1]
  if(SE) se.d.toalign=as.matrix(ses.raw[,which(substr(names(ses.raw),1,2)=="d.")])
  rownames(d.toalign)=names(coef.raw)
  if(SE) rownames(se.d.toalign)=names(se.raw)
  if(is.null(colnames(d.toalign))) colnames(d.toalign)='d.1'
  if(SE) if(is.null(colnames(se.d.toalign))) colnames(se.d.toalign)='d.1'

  #intercepts
  for(j in 1:ncol(d.toalign)){
    which.na=which(d.toalign[,j]==0)
    if(length(which.na)>0)d.toalign[which.na,j]=NA
    if(length(which.na)>0 & SE)se.d.toalign[which.na,j]=NA
  }

  #return
  if(!SE){
    se.a.toalign=NULL
    se.d.toalign=NULL
  }

  #set column names
  colnames(a.toalign)='G'
  colnames(se.a.toalign)='G'

  return(list(a=a.toalign,d=d.toalign,
              se.a=se.a.toalign,se.d=se.d.toalign))
}

#' Prepare \code{lavaan} estimates for alignment
#'
#' Not generally intended to be used on its own, but exported anyway for didactic purposes.
#'
#' See example for \code{\link{Alignment}} for examples
#'
#' This may differ from what is used in Mplus.
#'
#' @export
getEstimates.lavaan=function(fit,SE=T){
  #coefficients
  coef.raw=lavInspect(fit,'est')
  if(SE)se.raw=lavInspect(fit,'se')
  # if('lr.betas'%in%names(coef.raw))
  #   coef.raw=coef.raw[which(names(coef.raw)!='lr.betas')]
  # #convert to vector if appropriate; just get item parameters
  # if(SE) se.raw=coef.raw%>%map(function(x)if(is.matrix(x))x['SE',] else NULL)
  # coef.raw=coef.raw%>%map(function(x)if(is.matrix(x))x['par',] else x)
  # #item parameters
  # coef.raw=coef.raw[fit@Data$tabdata%>%colnames]
  # coef.raw=coef.raw[!is.na(names(coef.raw))]
  # if(SE) se.raw=se.raw[fit@Data$tabdata%>%colnames]
  # if(SE) se.raw=se.raw[!is.na(names(se.raw))]

  # coef.raw=map(coef.raw,~.[,which(.!=0)])

  #loadings
  lambda.raw=coef.raw$lambda
  #thresholds
  tau.raw=coef.raw$tau%>%as.data.frame%>%rownames_to_column('rowname')%>%
    mutate(Item=strsplit(rowname,'|',fixed=T)%>%map_chr(~.[1]),
           Threshold=strsplit(rowname,'|',fixed=T)%>%map_chr(~.[2]))%>%
    select(-rowname)%>%
    pivot_wider(id_cols=Item,names_from=Threshold,values_from=threshold)%>%
    as.data.frame
  rownames(tau.raw)=tau.raw$Item
  tau.raw=tau.raw%>%select(-Item)%>%as.matrix
  if(SE){
    #loadings
    se.lambda.raw=se.raw$lambda
    #thresholds
    se.tau.raw=se.raw$tau%>%as.data.frame%>%rownames_to_column('rowname')%>%
      mutate(Item=strsplit(rowname,'|',fixed=T)%>%map_chr(~.[1]),
             Threshold=strsplit(rowname,'|',fixed=T)%>%map_chr(~.[2]))%>%
      select(-rowname)%>%
      pivot_wider(id_cols=Item,names_from=Threshold,values_from=threshold)%>%
      as.data.frame
    rownames(se.tau.raw)=se.tau.raw$Item
    se.tau.raw=se.tau.raw%>%select(-Item)%>%as.matrix
  } else {
    se.lambda.raw=NULL
    se.tau.raw=NULL
  }

  #return
  return(list(lambda=lambda.raw,tau=tau.raw,
              se.lambda=se.lambda.raw,se.tau=se.tau.raw,
              parameterization=fit@Options$parameterization))
}

#' Transform \code{mirt} estimates using aligned estimates of latent mean and variance
#'
#' Not generally intended to be used on its own, but exported anyway for didactic purposes.
#'
#' See example for \code{\link{Alignment}} for examples
#'
#' This may differ from what is used in Mplus.
#'
#' @export
transformEstimates.mirt.grm=function(align.mean,align.variance,est){
  # est=out
  # align.mean=means
  # align.variance=variances
  # est=out
  #unpack est
  a=est$a
  d=est$d
  se.a=est$se.a
  se.d=est$se.d
  #transform a and SE's
  if(length(dim(est[[1]]))>2 & length(align.mean)==length(align.variance) & length(align.mean)==dim(est[[1]])[3]){
    #assume 3d arrays
    a=a*array(1/sqrt(align.variance),
              dim(a)[c(3,1,2)])%>%
      aperm(c(2,3,1))
    se.a=se.a*array(1/sqrt(align.variance),
                    dim(a)[c(3,1,2)])%>%
      aperm(c(2,3,1))
    d=d-(array(align.mean,
               dim(d)[c(3,1,2)])%>%
           aperm(c(2,3,1)))*
      (array(a,dim(d)[c(1,3,2)])%>%
         aperm(c(1,3,2)))
  } else if(length(dim(est[[1]]))==2 & length(align.mean)==1 & length(align.variance)==1){
    #assume one group
    a=a*1/sqrt(align.variance)
    se.a=se.a*1/sqrt(align.variance)
    d=d-align.mean*matrix(a,
                          nrow(d),
                          ncol(d),byrow=F)
  } else stop('Either transform one model (scalar mean & variance, 2D arrays in est)
or a set (vector mean & variance with equal lengths, 3D arrays in est, and
length(mean)==dim(est)[3]')
  # a=a*sqrt(1/align.variance)
  # if(!is.null(se.a))se.a=se.a*sqrt(1/align.variance)
  # d=d-align.mean*matrix(a,
  #                       nrow(d),
  #                       ncol(d),byrow=F)
  return(list(a=a,d=d,
              se.a=se.a,se.d=se.d))
}

#' Transform \code{lavaan} estimates using aligned estimates of latent mean and variance
#'
#' Not generally intended to be used on its own, but exported anyway for didactic purposes.
#'
#' See example for \code{\link{Alignment}} for examples
#'
#' This may differ from what is used in Mplus.
#'
#' @export
transformEstimates.lavaan.ordered=function(align.mean,align.variance,est,toCompare=NULL){
  #My current thinking: under the delta parameterization, the transformed
  #estimates (calculate delta, incorporate it into parameters, then
  #transform parameters, BUT don't reverse the delta transformation) do NOT
  #yield an equivalent model, but DO yield a model that can be compared
  #across groups. In order to get an equivalent model, you also need to
  #reverse the delta transformation at the end.
  #To account for this, the extra argument toCompare should be turned on if
  #transformed parameters are to be compared for equivalence across groups.
  #Turning it on results in NOT applying the reverse of the delta transformation
  #at the end.
  #As a result, I should really just use the theta parameterization throughout.

  # est=est.base[[2]]
  # align.mean=-1
  # align.variance=2

  # est=stackEstimates(est.base)
  # align.mean=c(0,1)
  # align.variance=c(1,2)

  #unpack est
  lambda=est$lambda
  tau=est$tau
  se.lambda=est$se.lambda
  se.tau=est$se.tau
  if(length(dim(est[[1]]))>2 & length(align.mean)==length(align.variance) & length(align.mean)==dim(est[[1]])[3]){
    if(all(est$parameterization=='theta')){
      #assume 3d arrays
      lambda=lambda*array(1/sqrt(align.variance),
                          dim(lambda)[c(3,1,2)])%>%
        aperm(c(2,3,1))
      se.lambda=se.lambda*array(1/sqrt(align.variance),
                                dim(lambda)[c(3,1,2)])%>%
        aperm(c(2,3,1))
      tau=tau+(array(align.mean,
                     dim(tau)[c(3,1,2)])%>%
                 aperm(c(2,3,1)))*
        (array(lambda,dim(tau)[c(1,3,2)])%>%
           aperm(c(1,3,2)))
    } else if(all(est$parameterization=='delta') & !is.null(toCompare)){
      #get deltas
      delta=sqrt(1-lambda^2)
      #convert to theta parameterization
      lambda=lambda/delta
      tau=tau/(array(delta,dim(tau)[c(1,3,2)])%>%
                 aperm(c(1,3,2)))
      #now, transform
      lambda=lambda*array(1/sqrt(align.variance),
                          dim(lambda)[c(3,1,2)])%>%
        aperm(c(2,3,1))
      se.lambda=se.lambda*array(1/sqrt(align.variance),
                                dim(lambda)[c(3,1,2)])%>%
        aperm(c(2,3,1))
      tau=tau+(array(align.mean,
                     dim(tau)[c(3,1,2)])%>%
                 aperm(c(2,3,1)))*
        (array(lambda,dim(tau)[c(1,3,2)])%>%
           aperm(c(1,3,2)))
      #convert back if !toCompare
      if(!toCompare){
        lambda=lambda*delta
        tau=tau*(array(delta,dim(tau)[c(1,3,2)])%>%
                   aperm(c(1,3,2)))
      }
    } else stop('Parameterization not found! Must be "delta" or "theta" and must be the same for all models')
  } else if(length(dim(est[[1]]))==2 & length(align.mean)==1 & length(align.variance)==1){
    if(est$parameterization=='theta'){
      #assume one group
      lambda=lambda*1/sqrt(align.variance)
      se.lambda=se.lambda*1/sqrt(align.variance)
      tau=tau+align.mean*matrix(lambda,
                                nrow(tau),
                                ncol(tau),byrow=F)
    } else if(est$parameterization=='delta' & !is.null(toCompare)){
      #get deltas
      delta=sqrt(1-lambda^2)
      #convert to theta parameterization
      lambda=lambda/delta
      tau=tau/matrix(delta,nrow(tau),ncol(tau),byrow=F)
      #now, transform
      lambda=lambda*1/sqrt(align.variance)
      se.lambda=se.lambda*1/sqrt(align.variance)
      tau=tau+align.mean*matrix(lambda,
                                nrow(tau),
                                ncol(tau),byrow=F)
      #convert back?
      if(!toCompare){
        lambda=lambda*delta
        tau=tau*matrix(delta,nrow(tau),ncol(tau),byrow=F)
      }
    } else stop('Parameterization not found! Must be "delta" or "theta"')
  } else stop('Either transform one model (scalar mean & variance, 2D arrays in est)
or a set (vector mean & variance with equal lengths, 3D arrays in est, and
length(mean)==dim(est)[3]. If delta pararameterization is used (the default in lavaan),
toCompare must be set to determine whether you want the parameters of an equivalent
model with incomparable parameters (toCompare=F) or comparable parameters (toCompare=T).')
  return(list(lambda=lambda,tau=tau,
              se.lambda=se.lambda,se.tau=se.tau))
}

#' Estimate \code{mirt} models using aligned parameter estimates
#'
#' Not generally intended to be used on its own, but exported anyway for didactic purposes.
#'
#' See example for \code{\link{Alignment}} for examples
#'
#' This may differ from what is used in Mplus.
#'
#' @export
loadEstimates.mirt.grm=function(fit,align.mean,align.variance,newpars,do.fit=T){
  # fit=fitList[[1]]
  # align.mean=means.vars[[1]][1]
  # align.variance=means.vars[[1]][2]
  # align.variance=means.vars[[1]][2]
  # newpars=test[[1]]
  #get call
  call=fit@Call%>%deparse
  if(length(call)>1)call=paste(call,collapse='')
  #if "fun" made its way in there, replace it with mirt
  call=gsub('fun(','mirt(',call,fixed=T)
  #replace data
  call.split=strsplit(call,'=',fixed=T)
  call.split=map(call.split[[1]],strsplit,',',fixed=T)%>%map(~.[[1]])
  call.split[[2]]=c('fit@Data$data',call.split[[2]][length(call.split[[2]])])
  call=call.split%>%map_chr(paste,collapse=',')%>%paste(collapse='=')
  #run it again with pars='values'
  call.pars=paste0(substr(call,1,nchar(call)-1),
                   ", pars = 'values')")
  pars=eval(parse(text=call.pars))
  #load 'em up: slopes
  for(i in 1:nrow(newpars$a)){
    for(j in 1:ncol(newpars$a)){
      if(!is.na(newpars$a[i]))
        pars$value[pars$item==rownames(newpars$a)[i] &
                     pars$name=='a1']=newpars$a[i]
    }
  }
  #load 'em up: slopes
  for(i in 1:nrow(newpars$d)){
    for(j in 1:ncol(newpars$d)){
      if(!is.na(newpars$d[i,j]))
        pars$value[pars$item==rownames(newpars$d)[i] &
                     pars$name==gsub(".",'',colnames(newpars$d)[j],fixed=T)]=newpars$d[i,j]
    }
  }
  #load 'em up: means and variances
  pars$value[pars$name=='MEAN_1']=align.mean
  pars$value[pars$name=='COV_11']=align.variance
  if(do.fit){
    #add pars
    newcall=paste0(substr(call,1,nchar(call)-1),
                   ", pars = pars)")
    out=eval(parse(text=newcall))
  } else out=pars
  out
}

#' Estimate \code{lavaan} models using aligned parameter estimates
#'
#' Not generally intended to be used on its own, but exported anyway for didactic purposes.
#'
#' See example for \code{\link{Alignment}} for examples
#'
#' This may differ from what is used in Mplus.
#'
#' @export
loadEstimates.lavaan.ordered=function(fit,align.mean,align.variance,newpars,do.fit=T){
  # fit=fit.lavaan
  # align.mean=1
  # align.variance=2
  # newpars=test.lavaan
  # do.fit=T
  #get data and partable
  dat=fit@Data
  pt=fit@ParTable%>%as_tibble
  #lv name
  fname=pt$lhs[pt$op=='=~']%>%unique

  #load 'em up: slopes
  for(i in 1:nrow(newpars$lambda)){
    for(j in 1:ncol(newpars$lambda)){
      pt$start[pt$op=='=~' &
                 pt$rhs==rownames(newpars$lambda)[i]]=newpars$lambda[i]
    }
  }
  #load 'em up: intercepts
  for(i in 1:nrow(newpars$tau)){
    for(j in 1:ncol(newpars$tau)){
      pt$start[pt$lhs==rownames(newpars$tau)[i] &
                 pt$rhs==gsub(".",'',colnames(newpars$tau)[j],fixed=T)]=newpars$tau[i,j]
    }
  }
  #load 'em up: means and variances
  pt$start[pt$lhs==fname & pt$op=='~1']=align.mean
  pt$start[pt$lhs==fname & pt$op=='~~' & pt$rhs==fname]=align.variance
  #just map everything from start to est
  pt$est=pt$start
  if(do.fit){
    out=lavaan(pt%>%as.list,data=dat,parameterization=fit@Options$parameterization)
  } else out=pt
  out
}

#' Stack estimates for optimization
#'
#' Not generally intended to be used on its own, but exported anyway for didactic purposes.
#'
#' See example for \code{\link{Alignment}} for examples
#'
#' This may differ from what is used in Mplus.
#'
#' @export
stackEstimates=function(estList){
  # estList=outByCohort%>%map(~.$fit.mirt.pss)%>%map(getEstimates.mirt,SE=T)
  # estList=est.base
  # estList=list(test.mirt%>%imap(function(x,n){
  #   rownames(x)[2]=paste0(rownames(x)[2],'_ho')
  #   if(!n%in%c('a','se.a'))colnames(x)[2]=paste0(colnames(x)[2],'_ho')
  #   x
  # }),test.mirt%>%imap(function(x,n){
  #   rownames(x)[1]=paste0(rownames(x)[1],'_hi')
  #   if(!n%in%c('a','se.a'))colnames(x)[1]=paste0(colnames(x)[1],'_hi')
  #   x
  # }))
  if(length(estList)==1){
    estList
  } else {
    #transpose
    # estList=sestList%>%transpose
    #output; build iteratively
    out=estList[[1]]
    for(i in 2:length(estList)){
      #get row and column names from each element
      curRows=out%>%map(rownames)
      curCols=out%>%map(colnames)
      newRows=estList[[i]]%>%map(rownames)
      newCols=estList[[i]]%>%map(colnames)
      for(j in 1:length(out)){
        if(names(estList[[i]])[j]=='parameterization'){
          out[[j]]=c(out[[j]],estList[[i]][[j]])
        } else {
          #add empty new rows to out
          if(any(!newRows[[j]]%in%curRows[[j]])){
            for(n in newRows[[j]][!newRows[[j]]%in%curRows[[j]]]){
              dim.out=dim(out[[j]])
              dim.out[1]=1
              out[[j]]=abind(out[[j]],array(NA,dim.out),along=1)
              rownames(out[[j]])[nrow(out[[j]])]=n
            }
          }
          if(any(!newCols[[j]]%in%curCols[[j]])){
            for(n in newCols[[j]][!newCols[[j]]%in%curCols[[j]]]){
              dim.out=dim(out[[j]])
              dim.out[2]=1
              out[[j]]=abind(out[[j]],array(NA,dim.out),along=2)
              colnames(out[[j]])[ncol(out[[j]])]=n
            }
          }
          if(any(!curRows[[j]]%in%newRows[[j]])){
            for(n in curRows[[j]][!curRows[[j]]%in%newRows[[j]]]){
              dim.out=dim(estList[[i]][[j]])
              dim.out[1]=1
              estList[[i]][[j]]=abind(estList[[i]][[j]],array(NA,dim.out),along=1)
              rownames(estList[[i]][[j]])[nrow(estList[[i]][[j]])]=n
            }
          }
          if(any(!curCols[[j]]%in%newCols[[j]])){
            for(n in curCols[[j]][!curCols[[j]]%in%newCols[[j]]]){
              dim.out=dim(estList[[i]][[j]])
              dim.out[2]=1
              estList[[i]][[j]]=abind(estList[[i]][[j]],array(NA,dim.out),along=2)
              colnames(estList[[i]][[j]])[ncol(estList[[i]][[j]])]=n
            }
          }
          #ugh
          what2add=estList[[i]][[j]][rownames(out[[j]]),colnames(out[[j]])]
          if(!is.matrix(what2add)){
            what2add=matrix(what2add,nrow(estList[[i]][[j]]),ncol(estList[[i]][[j]]))
            rownames(what2add)=rownames(estList[[i]][[j]])
            colnames(what2add)=colnames(estList[[i]][[j]])
          }
          out[[j]]=abind(out[[j]],what2add,along=3)
        }
      }
    }
    out
  }
}

#simplicity function
CLF=function(x,e=.01){
  return(sqrt(sqrt(x^2+e)))
}

#weight function
W=function(x,y) return(sqrt(x*y))

#' Simplicity function for alignment
#'
#' Not generally intended to be used on its own, but exported anyway for didactic purposes.
#'
#' See example for \code{\link{Alignment}} for examples
#'
#' This may differ from what is used in Mplus.
#'
#' @export
SF.mplus3D=function(pars,est,comb,nobs,estimator){
  # pars=c(rnorm(ngroups-1,0,5),abs(rnorm(ngroups-1,0,3)))
  # est=stacked
  # comb=comb
  # nobs=n
  # estimator=estimator

  #constants
  # nitems=nrow(mA)
  #unlist sample sizes
  if(is.list(nobs))nobs=unlist(nobs)
  ngroups=length(nobs)

  #extract hyperparameters
  means=pars[1:(ngroups-1)]
  variances=pars[((ngroups-1)+1):(2*(ngroups-1))]
  means=c(0,means)
  variances=c(1/prod(variances),variances)

  #matrices for alignment transformation
  # means1=matrix(means[comb[1,]],nitems,dim(mA)[3],byrow=T)
  # means2=matrix(means[comb[2,]],nitems,dim(mB)[3],byrow=T)
  # variances1=matrix(variances[comb[1,]],nitems,dim(mA)[3],byrow=T)
  # variances2=matrix(variances[comb[2,]],nitems,dim(mB)[3],byrow=T)

  #transform parameters
  if(estimator=='mirt.grm'){
    t.est=transformEstimates.mirt.grm(means,variances,est)
    t.est=t.est[c('a','d')]
  } else if(estimator=='lavaan.ordered'){
    t.est=transformEstimates.lavaan.ordered(means,variances,est,toCompare=T)
    t.est=t.est[c('lambda','tau')]
  }
  # mA[,1,]=mA[,1,]/sqrt(variances1)
  # mB[,1,]=mB[,1,]/sqrt(variances2)
  # dmA2=aperm(array(means1*mA[,1,],dim=c(dim(means1),1)),c(1,3,2))
  # dmA2=do.call(abind,list(rep(list(dmA2),dim(mA)[2]-1),along=2))
  # if(dim(dmA2)[3]==1)dmA2=dmA2[,,1]
  # dmB2=aperm(array(means2*mB[,1,],dim=c(dim(means2),1)),c(1,3,2))
  # dmB2=do.call(abind,list(rep(list(dmB2),dim(mB)[2]-1),along=2))
  # if(dim(dmB2)[3]==1)dmB2=dmB2[,,1]
  # mA[,2:dim(mA)[2],]=mA[,2:dim(mA)[2],]-dmA2
  # mB[,2:dim(mB)[2],]=mB[,2:dim(mB)[2],]-dmB2

  # m=sum(CLF(mA-mB),na.rm=T)
  t.est%>%map(function(x)CLF(x[,,comb[1,]]-x[,,comb[2,]])*
                W(nobs[comb[1,]],nobs[comb[2,]])/
                sum(!is.na(x[,,comb[1,]]-x[,,comb[2,]])))%>%Reduce(c,.)%>%sum(na.rm=T)

  #shrink effect of sample size - log transform
  # nobs=log(nobs)
  # w=W(nobs[comb[1,]],nobs[comb[2,]])
  # # w=rep(1,ncol(comb))
  # x=sum(m*w)
  # x/sum(!is.na(mA-mB))
}

#factory function
factory <- function(fun){
  function(...) {
    warn <- err <- NULL
    res <- withCallingHandlers(
      tryCatch(fun(...), error=function(e) {
        err <<- conditionMessage(e)
        NULL
      }), warning=function(w) {
        warn <<- append(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    list(res, warn=warn, err=err)
  }}
fact.mirt=factory(mirt)

#' Runs alignment optimizer
#'
#' Not generally intended to be used on its own, but exported anyway for didactic purposes.
#'
#' See example for \code{\link{Alignment}} for examples
#'
#' This may differ from what is used in Mplus.
#'
#' @export
align.optim=function(stacked,n,estimator,nstarts=50,ncores=3,parallel=F,center.means){

  # stacked=stackEstimates(est.base)
  # n=c(100,200)
  # nstarts=3
  # ncores=3
  # estimator='lavaan.ordered'
  # center.means=F

  # stacked=test.stack2
  # n=c(100,200)
  # nstarts=3
  # ncores=3
  # estimator='mirt.grm'

  #get sample sizes from mirt objects
  ngroups=length(n)

  if(ngroups>1){
    #matrices to subtract
    comb=combn(1:ngroups,2)

    #############################################################################
    #############################################################################
    #############################################################################
    #############################################################################

    cat("Parameters ready to align, beginning alignment... ")

    #alignment optimization
    parmx=NULL
    #align
    pct=proc.time()
    if(parallel){
      library(doParallel)
      cl=makeCluster(ncores)
      registerDoParallel(cl)
      parmx=foreach(k = 1:nstarts,.export=c("W","SF.mplus3D","CLF",
                                            paste0('transformEstimates.',estimator)),
                    .packages=c("abind",'tidyverse')) %dopar% {
                      set.seed(k)
                      out=optim(c(rnorm(ngroups-1,0,5),abs(rnorm(ngroups-1,0,3))),
                                SF.mplus3D,method="L-BFGS-B",
                                est=stacked,comb=comb,nobs=n,estimator=estimator,
                                control=list(maxit=10000,trace=0),
                                lower=c(rep(-Inf,ngroups-1),rep(1e-6,ngroups-1)))
                      c(value=out$value,convergence=out$convergence,par=out$par)
                    }
      stopCluster(cl)
    } else {
      parmx=list()
      for(k in 1:nstarts) {
        # set.seed(k)
        out=optim(c(rnorm(ngroups-1,0,5),abs(rnorm(ngroups-1,0,3))),
                  SF.mplus3D,method="L-BFGS-B",
                  est=stacked,comb=comb,nobs=n,estimator=estimator,
                  control=list(maxit=10000,trace=0),
                  lower=c(rep(-Inf,ngroups-1),rep(1e-6,ngroups-1)))
        parmx[[k]]=c(value=out$value,convergence=out$convergence,par=out$par)
      }
    }
    proc.time()-pct
    parout=parmx%>%setNames(1:nstarts)%>%bind_cols%>%t

    #reject runs that failed
    failedruns=which(parout[,2]!=0)
    if(length(failedruns)>0)parout=parout[-failedruns,]

    #name parameter list
    print("Printing parout...")
    print(parout)
    print(str(parout))
    print(ngroups)
    colnames(parout)=c("f","convergence",
                       paste0("M.",2:ngroups),
                       paste0("V.",2:ngroups))

    #populate
    align.hyperpars=parout[order(parout[,1],decreasing=F)[1],-c(1,2)]
    align.means=align.hyperpars[1:(ngroups-1)]
    align.variances=align.hyperpars[((ngroups-1)+1):((ngroups-1)*2)]
    #populate means, centering at zero to start
    align.means=c(0,align.means)
    align.variances=c(1/prod(align.variances),align.variances)
    #re-center at weighted averages for grand mean of zero
    if(center.means)align.means=align.means-weighted.mean(align.means,unlist(n))

    #re-center such that weighted product is 1
    weighted.prod=exp(weighted.mean(log(align.variances),unlist(n)))
    align.variances=align.variances/weighted.prod
    weighted.prod=exp(weighted.mean(log(align.variances),unlist(n)))
    weighted.prod #nice

    #means are in standard deviation units for first group; divide by SD of first group
    align.means=align.means/sqrt(align.variances[1])
    #variances just need to be divided by the first one
    align.variances=align.variances/align.variances[1]
  } else {
    align.means=0
    align.variances=1
  }

  #return just means and variances
  return(mapply(function(x,y)c(mean=x,var=y),as.list(align.means),as.list(align.variances),SIMPLIFY=F))
}
