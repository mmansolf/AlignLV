#' Prepare \code{mirt} estimates for alignment
#'
#' Not generally intended to be used on its own, but exported anyway for
#' didactic purposes.
#'
#' See example for \code{\link{Alignment}} for examples
#'
#' This may differ from what is used in Mplus.
#'
#' @param fit A \code{mirt} object compatible with \code{\link{Alignment}}
#' @param SE logical; whether to also obtain standard errors.
#' @param bifactor.marginal See \code{\link{Alignment}} documentation.
#'
#' @returns A `list` of estimates in a format amenable to subsequent alignment
#'
#' @export
getEstimates.mirt=function(fit,SE=FALSE,bifactor.marginal=FALSE){

  #coefficients
  coef.raw=mirt::coef(fit,printSE=SE)
  if('lr.betas'%in%names(coef.raw))
    coef.raw=coef.raw[which(names(coef.raw)!='lr.betas')]
  #convert to vector if appropriate; just get item parameters
  if(SE) se.raw=coef.raw%>%purrr::map(function(x)if(is.matrix(x))x['SE',] else NULL)
  coef.raw=coef.raw%>%purrr::map(function(x)if(is.matrix(x))x['par',] else x)
  #item parameters
  coef.raw=coef.raw[fit@Data$tabdata%>%colnames]
  coef.raw=coef.raw[!is.na(names(coef.raw))]
  if(SE) se.raw=se.raw[fit@Data$tabdata%>%colnames]
  if(SE) se.raw=se.raw[!is.na(names(se.raw))]


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
    if(SE) rownames(se.a.raw)[j]=names(coef.raw)[j]
  }

  #maximum number of categories maxcats
  maxcats=10
  d.raw=matrix(0,1,maxcats)
  if(SE) se.d.raw=matrix(0,1,maxcats)
  for(j in 1:(length(coef.raw))){
    dindices=which(substr(names(coef.raw[[j]]),1,1)=="d")
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
  #shift intercepts according to data mins
  mins=(fit@Data$mins-min(fit@Data$mins))%>%
    stats::setNames(NULL) #fixed min to 0
  d.toalign2=d.toalign
  if(SE) se.d.toalign2=se.d.toalign
  for(r in 1:nrow(d.toalign2)){
    if(mins[r]>0){
      #calculate # of columns to add
      whichExist=which(!is.na(d.toalign2[r,]))
      nCols2add=whichExist[length(whichExist)]-whichExist[1]+mins[r]-ncol(d.toalign2)+1
      # if(!is.na(d.toalign2[r,ncol(d.toalign2)])){
      if(nCols2add>0){
        for(i in 1:nCols2add){
          d.toalign2=cbind(d.toalign2,NA)
          colnames(d.toalign2)[ncol(d.toalign2)]=
            paste0('d.',
                   as.numeric(
                     strsplit(colnames(d.toalign2)[ncol(d.toalign2)-1],'.',fixed=TRUE)%>%purrr::map_chr(~.[2])
                   )+1)
          if(SE){
            se.d.toalign2=cbind(se.d.toalign2,NA)
            colnames(se.d.toalign2)[ncol(se.d.toalign2)]=
              paste0('d.',
                     as.numeric(
                       strsplit(colnames(se.d.toalign2)[ncol(se.d.toalign2)-1],'.',fixed=TRUE)%>%purrr::map_chr(~.[2])
                     )+1)
          }
        }
      }
      d.toalign2[r,(mins[r]+1):(mins[r]+sum(!is.na(d.toalign2[r,])))]=
        d.toalign2[r,1:(sum(!is.na(d.toalign2[r,])))]
      d.toalign2[r,1:mins[r]]=NA
      if(SE){
        se.d.toalign2[r,(mins[r]+1):(mins[r]+sum(!is.na(se.d.toalign2[r,])))]=
          se.d.toalign2[r,1:(sum(!is.na(se.d.toalign2[r,])))]
        se.d.toalign2[r,1:mins[r]]=NA
      }
    }
  }
  d.toalign=d.toalign2
  if(SE) se.d.toalign=se.d.toalign2
  #set column names
  colnames(d.toalign)=paste0('d.',
                             min(mins):(ncol(d.toalign)-min(mins)-1),
                             '.',
                             (min(mins)+1):(ncol(d.toalign)-min(mins)))
  if(SE) colnames(se.d.toalign)=paste0('d.',
                                       min(mins):(ncol(se.d.toalign)-min(mins)-1),
                                       '.',
                                       (min(mins)+1):(ncol(se.d.toalign)-min(mins)))

  #return
  if(!SE){
    se.a.toalign=NULL
    se.d.toalign=NULL
  } else {
    colnames(se.a.toalign)='G'
  }

  #set column names
  colnames(a.toalign)='G'

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
#' @param fit A \code{lavaan} object compatible with \code{\link{Alignment}}
#' @param SE logical; whether to also obtain standard errors.
#'
#' @export
getEstimates.lavaan=function(fit,SE=TRUE){

  #coefficients
  coef.raw=lavaan::lavInspect(fit,'est')
  if(SE)se.raw=lavaan::lavInspect(fit,'se')

  #loadings
  lambda.raw=coef.raw$lambda
  #thresholds
  tau.raw=coef.raw$tau%>%as.data.frame%>%tibble::rownames_to_column('rowname')
  tau.raw=tau.raw%>%
    dplyr::mutate(Item=strsplit(.data$rowname,'|',fixed=T)%>%purrr::map_chr(~.[1]),
           Threshold=strsplit(.data$rowname,'|',fixed=T)%>%purrr::map_chr(~.[2]))%>%
    dplyr::select(-.data$rowname)
  #apply observed levels to thresholds
  obs.thresh=fit@Data@ov$lnam%>%
    strsplit('|',fixed=T)%>%
    purrr::map(~paste0(.[1:(length(.)-1)],'_',.[2:length(.)]))%>%
    purrr::map(~tibble::tibble(Threshold=as.character(1:length(.)),
                boundary=.))%>%
    stats::setNames(fit@Data@ov$name)%>%
    purrr::imap(~dplyr::mutate(.x,Item=.y))%>%dplyr::bind_rows()
  tau.raw=tau.raw%>%
    dplyr::mutate(Threshold=gsub('t','',.data$Threshold,fixed=T))%>%
    dplyr::left_join(obs.thresh,by=c('Item','Threshold'))
  tau.raw=tau.raw%>%
    dplyr::mutate(boundary=paste0('t',.data$boundary))%>%
    dplyr::select(-.data$Threshold)%>%dplyr::rename(Threshold=.data$boundary)
  #finish preparing
  tau.raw=tau.raw%>%
    tidyr::pivot_wider(id_cols=.data$Item,
                       names_from=.data$Threshold,
                       values_from=.data$threshold)%>%
    as.data.frame
  rownames(tau.raw)=tau.raw$Item
  tau.raw=tau.raw%>%dplyr::select(-.data$Item)%>%as.matrix
  if(SE){
    #loadings
    se.lambda.raw=se.raw$lambda
    #thresholds
    se.tau.raw=se.raw$tau%>%as.data.frame%>%
      tibble::rownames_to_column('rowname')%>%
      dplyr::mutate(Item=strsplit(.data$rowname,'|',fixed=T)%>%
                      purrr::map_chr(~.[1]),
             Threshold=strsplit(.data$rowname,'|',fixed=T)%>%
               purrr::map_chr(~.[2]))%>%
      dplyr::select(-.data$rowname)%>%
      dplyr::mutate(Threshold=gsub('t','',.data$Threshold,fixed=T))%>%
      dplyr::left_join(obs.thresh,by=c('Item','Threshold'))%>%
      dplyr::mutate(boundary=paste0('t',.data$boundary))%>%
      dplyr::select(-.data$Threshold)%>%dplyr::rename(Threshold=.data$boundary)%>%
      tidyr::pivot_wider(id_cols=.data$Item,
                         names_from=.data$Threshold,
                         values_from=.data$threshold)%>%
      as.data.frame
    rownames(se.tau.raw)=se.tau.raw$Item
    se.tau.raw=se.tau.raw%>%dplyr::select(-.data$Item)%>%as.matrix
  } else {
    se.lambda.raw=NULL
    se.tau.raw=NULL
  }

  #return
  return(list(lambda=lambda.raw,tau=tau.raw,
              se.lambda=se.lambda.raw,se.tau=se.tau.raw,
              parameterization=fit@Options$parameterization))
}

#' Transform \code{mirt} estimates using aligned estimates of latent mean
#' and variance
#'
#' Not generally intended to be used on its own, but exported anyway for
#' didactic purposes.
#'
#' See example for \code{\link{Alignment}} for examples
#'
#' This may differ from what is used in Mplus.
#'
#' @param align.mean Mean to transform model to.
#' @param align.variance Variance to transform model to.
#' @param est Estimates to transform, from \code{\link{getEstimates.mirt}}
#'
#' @export
transformEstimates.mirt.grm=function(align.mean,align.variance,est){

  #unpack est
  a=est$a
  d=est$d
  se.a=est$se.a
  se.d=est$se.d
  if(is.null(se.a)) SE=FALSE else SE=TRUE
  #transform a and SE's
  if(length(dim(est[[1]]))>2 & length(align.mean)==length(align.variance) &
     length(align.mean)==dim(est[[1]])[3]){
    #assume 3d arrays
    a=a*array(1/sqrt(align.variance),
              dim(a)[c(3,1,2)])%>%
      aperm(c(2,3,1))
    if(SE) se.a=se.a*array(1/sqrt(align.variance),
                           dim(a)[c(3,1,2)])%>%
        aperm(c(2,3,1))
    d=d-(array(align.mean,
               dim(d)[c(3,1,2)])%>%
           aperm(c(2,3,1)))*
      (array(a,dim(d)[c(1,3,2)])%>%
         aperm(c(1,3,2)))
  } else if(length(dim(est[[1]]))==2 &
            length(align.mean)==1 &
            length(align.variance)==1){
    #assume one group
    a=a*1/sqrt(align.variance)
    if(SE) se.a=se.a*1/sqrt(align.variance)
    d=d-align.mean*matrix(a,
                          nrow(d),
                          ncol(d),byrow=FALSE)
  } else stop(paste0('Either transform one model (scalar mean & variance, ',
  '2D arrays in est) or a set (vector mean & variance with equal lengths, ',
  '3D arrays in est, and length(mean)==dim(est)[3]'))
  out=list(a=a,d=d)
  if(SE) out=c(out,list(se.a=se.a,se.d=se.d))
  return(out)
}

#' Transform \code{lavaan} estimates using aligned estimates of latent
#' mean and variance
#'
#' Not generally intended to be used on its own, but exported anyway for
#' didactic purposes.
#'
#' See example for \code{\link{Alignment}} for examples
#'
#' This may differ from what is used in Mplus.
#'
#' @param align.mean Mean to transform model to.
#' @param align.variance Variance to transform model to.
#' @param est Estimates to transform, from \code{\link{getEstimates.lavaan}}
#' @param toCompare Accounts for discrepancies between delta and theta
#' parameterizations; see \code{\link{Alignment}} documentation.
#'
#' @export
transformEstimates.lavaan.ordered=function(align.mean,
                                           align.variance,
                                           est,toCompare=FALSE){

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

  #unpack est
  lambda=est$lambda
  tau=est$tau
  se.lambda=est$se.lambda
  se.tau=est$se.tau
  if(length(dim(est[[1]]))>2 & length(align.mean)==length(align.variance) &
     length(align.mean)==dim(est[[1]])[3]){
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
    } else if(all(est$parameterization=='delta') & toCompare){
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
    } else stop(paste0('Parameterization not found! Must be "delta" or "theta"',
                       'and must be the same for all models'))
  } else if(length(dim(est[[1]]))==2 &
            length(align.mean)==1 &
            length(align.variance)==1){
    if(est$parameterization=='theta'){
      #assume one group
      lambda=lambda*1/sqrt(align.variance)
      se.lambda=se.lambda*1/sqrt(align.variance)
      tau=tau+align.mean*matrix(lambda,
                                nrow(tau),
                                ncol(tau),byrow=FALSE)
    } else if(est$parameterization=='delta' & toCompare){
      #get deltas
      delta=sqrt(1-lambda^2)
      #convert to theta parameterization
      lambda=lambda/delta
      tau=tau/matrix(delta,nrow(tau),ncol(tau),byrow=FALSE)
      #now, transform
      lambda=lambda*1/sqrt(align.variance)
      se.lambda=se.lambda*1/sqrt(align.variance)
      tau=tau+align.mean*matrix(lambda,
                                nrow(tau),
                                ncol(tau),byrow=FALSE)
      #convert back?
      if(!toCompare){
        lambda=lambda*delta
        tau=tau*matrix(delta,nrow(tau),ncol(tau),byrow=FALSE)
      }
    } else stop('Parameterization not found! Must be "delta" or "theta"')
  } else stop(paste0('Either transform one model (scalar mean & variance, ',
                     '2D arrays in est) or a set (vector mean & variance ',
                     'with equal lengths, 3D arrays in est, and ',
                     'length(mean)==dim(est)[3]. If delta pararameterization ',
                     'is used (the default in lavaan), toCompare must be set ',
                     'to determine whether you want the parameters of an ',
                     'equivalent model with incomparable parameters ',
                     '(toCompare=FALSE) or comparable parameters (toCompare=T).'))
  return(list(lambda=lambda,tau=tau,
              se.lambda=se.lambda,se.tau=se.tau))
}

#' Estimate \code{mirt} models using aligned parameter estimates
#'
#' Not generally intended to be used on its own, but exported anyway for
#' didactic purposes.
#'
#' See example for \code{\link{Alignment}} for examples
#'
#' This may differ from what is used in Mplus.
#'
#' @param fit A \code{mirt} object compatible with \code{\link{Alignment}}
#' @param align.mean Mean to transform model to.
#' @param align.variance Variance to transform model to.
#' @param newpars New (transformed) estimates to load into model object.
#' @param do.fit Whether to re-fit the model after loading and fixing estimates.
#'
#' @export
loadEstimates.mirt.grm=function(fit,align.mean,
                                align.variance,newpars,do.fit=TRUE){

  #get call
  call=fit@Call%>%deparse
  if(length(call)>1)call=paste(call,collapse='')
  #if "fun" made its way in there, replace it with mirt
  call=gsub('fun(','mirt(',call,fixed=TRUE)
  call=gsub('..4','fit@Options$technical',call,fixed=TRUE)
  #replace data
  call.split=strsplit(call,'=',fixed=TRUE)
  call.split=purrr::map(call.split[[1]],strsplit,',',fixed=TRUE)%>%purrr::map(~.[[1]])
  call.split[[2]]=c('fit@Data$data',call.split[[2]][length(call.split[[2]])])
  call=call.split%>%purrr::map_chr(paste,collapse=',')%>%paste(collapse='=')
  #run it again with pars='values'
  call.pars=paste0(substr(call,1,nchar(call)-1),
                   ", pars = 'values')")
  # call.pars=gsub('pars = pars,','',call.pars,fixed=TRUE)
  #chatgpt gave me this
  call.pars=gsub("pars = ([^,\\)]+)([,\\)])", "pars = \\1\\2", call.pars, perl = TRUE)

  #bring in mod from @Model$model if model = mod is contained in call.pars
  if(!grepl('model = 1',call.pars,fixed=TRUE)){
    modName=strsplit(call.pars,'model = ',fixed=TRUE)[[1]][2]%>%
      strsplit(',',fixed=TRUE)%>%purrr::map_chr(~.[1])%>%trimws
    call.pars=gsub(modName,'myMod',call.pars,fixed=TRUE)
    call=gsub(modName,'myMod',call,fixed=TRUE)
    myMod=fit@Model$model
  }
  pars=eval(parse(text=call.pars))
  #load 'em up: slopes
  for(i in 1:nrow(newpars$a)){
    for(j in 1:ncol(newpars$a)){
      if(!is.na(newpars$a[i]))
        pars$value[pars$item==rownames(newpars$a)[i] &
                     pars$name=='a1']=newpars$a[i]
    }
  }
  #load 'em up: locations
  #mins
  mins=(fit@Data$mins-min(fit@Data$mins))[
    rownames(newpars$d)
  ]
  for (i in 1:nrow(newpars$d)) {
    for (j in 1:ncol(newpars$d)) {
      if (!is.na(newpars$d[i,j]))
        pars$value[pars$item == rownames(newpars$d)[i] &
                     pars$name == gsub(".","",colnames(newpars$d)[j],
                                       fixed=T)%>%purrr::map_chr(function(x)paste0(
                                         substr(x,1,1),
                                         as.numeric(substr(x,3,3))-mins[i]
                                       ))]=newpars$d[i,j]
    }
  }
  #load 'em up: means and variances
  pars$value[pars$name=='MEAN_1']=align.mean
  pars$value[pars$name=='COV_11']=align.variance
  if(do.fit){
    #add pars, if not there already (?)
    if(!grepl('pars = pars',call,fixed=TRUE)){
      newcall=paste0(substr(call,1,nchar(call)-1),
                     ", pars = pars)")
      out=eval(parse(text=newcall))
    } else out=eval(parse(text=call))
  } else out=pars
  out
}

#' Estimate \code{lavaan} models using aligned parameter estimates
#'
#' Not generally intended to be used on its own, but exported anyway for
#' didactic purposes.
#'
#' See example for \code{\link{Alignment}} for examples
#'
#' This may differ from what is used in Mplus.
#'
#' @param fit A \code{mirt} object compatible with \code{\link{Alignment}}
#' @param align.mean Mean to transform model to.
#' @param align.variance Variance to transform model to.
#' @param newpars New (transformed) estimates to load into model object.
#' @param do.fit Whether to re-fit the model after loading and fixing estimates.
#'
#' @export
loadEstimates.lavaan.ordered=function(fit,align.mean,align.variance,
                                      newpars,do.fit=TRUE){

  #get data and partable
  dat=fit@Data
  pt=fit@ParTable%>%tibble::as_tibble()
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
                 pt$rhs==gsub(".",'',colnames(newpars$tau)[j],fixed=TRUE)]=
        newpars$tau[i,j]
    }
  }
  #load 'em up: means and variances
  pt$start[pt$lhs==fname & pt$op=='~1']=align.mean
  pt$start[pt$lhs==fname & pt$op=='~~' & pt$rhs==fname]=align.variance
  #just map everything from start to est
  pt$est=pt$start
  if(do.fit){
    out=lavaan::lavaan(pt%>%as.list,data=dat,parameterization=
                 fit@Options$parameterization)
  } else out=pt
  out
}

#' Stack estimates for optimization
#'
#' Not generally intended to be used on its own, but exported anyway for
#' didactic purposes.
#'
#' See example for \code{\link{Alignment}} for examples
#'
#' This may differ from what is used in Mplus.
#'
#' @param estList List of estimates from \code{\link{getEstimates.lavaan}} or
#' \code{\link{getEstimates.mirt}} to stack to feed into \code{\link{SF.mplus3D}}
#'
#' @export
stackEstimates=function(estList){

  if(length(estList)==1){
    estList
  } else {
    #output; build iteratively
    out=estList[[1]]
    for(i in 2:length(estList)){
      #get row and column names from each element
      curRows=out%>%purrr::map(rownames)
      curCols=out%>%purrr::map(colnames)
      newRows=estList[[i]]%>%purrr::map(rownames)
      newCols=estList[[i]]%>%purrr::map(colnames)
      for(j in 1:length(out)){
        if(is.null(estList[[i]][[j]])) next
        if(names(estList[[i]])[j]=='parameterization'){
          out[[j]]=c(out[[j]],estList[[i]][[j]])
        } else {
          #add empty new rows to out
          if(any(!newRows[[j]]%in%curRows[[j]])){
            for(n in newRows[[j]][!newRows[[j]]%in%curRows[[j]]]){
              dim.out=dim(out[[j]])
              dim.out[1]=1
              out[[j]]=abind::abind(out[[j]],array(NA,dim.out),along=1)
              rownames(out[[j]])[nrow(out[[j]])]=n
            }
          }
          if(any(!newCols[[j]]%in%curCols[[j]])){
            for(n in newCols[[j]][!newCols[[j]]%in%curCols[[j]]]){
              dim.out=dim(out[[j]])
              dim.out[2]=1
              out[[j]]=abind::abind(out[[j]],array(NA,dim.out),along=2)
              colnames(out[[j]])[ncol(out[[j]])]=n
            }
          }
          if(any(!curRows[[j]]%in%newRows[[j]])){
            for(n in curRows[[j]][!curRows[[j]]%in%newRows[[j]]]){
              dim.out=dim(estList[[i]][[j]])
              dim.out[1]=1
              estList[[i]][[j]]=abind::abind(estList[[i]][[j]],
                                             array(NA,dim.out),along=1)
              rownames(estList[[i]][[j]])[nrow(estList[[i]][[j]])]=n
            }
          }
          if(any(!curCols[[j]]%in%newCols[[j]])){
            for(n in curCols[[j]][!curCols[[j]]%in%newCols[[j]]]){
              dim.out=dim(estList[[i]][[j]])
              dim.out[2]=1
              estList[[i]][[j]]=abind::abind(estList[[i]][[j]],
                                             array(NA,dim.out),along=2)
              colnames(estList[[i]][[j]])[ncol(estList[[i]][[j]])]=n
            }
          }
          #ugh
          what2add=estList[[i]][[j]][rownames(out[[j]]),colnames(out[[j]])]
          if(!is.matrix(what2add)){
            what2add=matrix(what2add,nrow(estList[[i]][[j]]),
                            ncol(estList[[i]][[j]]))
            rownames(what2add)=rownames(estList[[i]][[j]])
            colnames(what2add)=colnames(estList[[i]][[j]])
          }
          out[[j]]=abind::abind(out[[j]],what2add,along=3)
        }
      }
    }
    out
  }
}

#simplicity function
CLF=function(x,e=0.01){
  return(sqrt(sqrt(x^2+e)))
}

#weight function
W=function(x,y) return(sqrt(x*y))

#' Simplicity function for alignment
#'
#' Not generally intended to be used on its own, but exported anyway for
#' didactic purposes.
#'
#' See example for \code{\link{Alignment}} for examples
#'
#' This may differ from what is used in Mplus.
#'
#' @param pars Hyperparameters to feed into optimizer
#' @param est Estimates to transform, from \code{\link{getEstimates.mirt}} or
#' \code{\link{getEstimates.lavaan}}
#' @param comb All combinations of groups from \code{\link[utils]{combn}}
#' @param nobs Sample size in each group
#' @param estimator See \code{\link{Alignment}} documentation.
#' @param eps.alignment See \code{\link{Alignment}} documentation.
#' @param clf.ignore.quantile See \code{\link{Alignment}} documentation.
#' @param hyper Hyperparameter to calculate simplicity function for; see
#' \code{\link{Alignment}} documentation.
#' @param otherHyper Non-included hyperparameter
#'
#' @export
SF.mplus3D=function(pars,est,comb,nobs,estimator,
                    eps.alignment,clf.ignore.quantile,
                    hyper='all',otherHyper=NULL){

  #unlist sample sizes
  if(is.list(nobs))nobs=unlist(nobs)
  ngroups=length(nobs)

  #extract hyperparameters
  if(hyper=='all'){
    means=pars[1:(ngroups-1)]
    variances=pars[((ngroups-1)+1):(2*(ngroups-1))]
  } else if(hyper=='means'){
    means=pars
    variances=otherHyper
  } else if(hyper=='variances'){
    means=otherHyper
    variances=pars
  }
  means=c(0,means)
  variances=c(1/prod(variances),variances)

  #transform parameters
  if(estimator=='mirt.grm'){
    t.est=transformEstimates.mirt.grm(means,variances,est)
    t.est=t.est[c('a','d')]
  } else if(estimator=='lavaan.ordered'){
    t.est=transformEstimates.lavaan.ordered(means,variances,est,toCompare=TRUE)
    t.est=t.est[c('lambda','tau')]
  }
  out=t.est%>%purrr::map(function(x)CLF(x[,,comb[1,]]-x[,,comb[2,]],e=eps.alignment)*
                W(nobs[comb[1,]],nobs[comb[2,]])/
                sum(!is.na(x[,,comb[1,]]-x[,,comb[2,]])))
  out=Reduce(c,out)
  out=ifelse(out<stats::quantile(out,clf.ignore.quantile,na.rm=TRUE),
           0,out)%>%sum(na.rm=TRUE)
}

#' Runs alignment optimizer
#'
#' Not generally intended to be used on its own, but exported anyway
#' for didactic purposes.
#'
#' See example for \code{\link{Alignment}} for examples
#'
#' @param stacked Stacked parameter estimates from \code{\link{stackEstimates}}
#' @param n Sample size in each group
#' @param estimator See \code{\link{Alignment}} documentation.
#' @param nstarts Number of starting values for alignment; default is 10
#' @param ncores See \code{\link{Alignment}} documentation.
#' @param hyper.first See \code{\link{Alignment}} documentation.
#' @param center.means See \code{\link{Alignment}} documentation.
#' @param eps.alignment See \code{\link{Alignment}} documentation.
#' @param clf.ignore.quantile See \code{\link{Alignment}} documentation.
#'
#' @export
align.optim=function(stacked,n,estimator,nstarts=50,ncores=3,
                     hyper.first,center.means,eps.alignment,
                     clf.ignore.quantile){

  #get sample sizes from mirt objects
  ngroups=length(n)

  if(ngroups>1){
    #matrices to subtract
    comb=utils::combn(1:ngroups,2)

    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################

    cat("Parameters ready to align, beginning alignment... ")

    #alignment optimization
    parmx=NULL
    #align
    pct=proc.time()
    if(ncores>1 & requireNamespace("doParallel", quietly = TRUE)){
      cl=parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
      `%dopar%` <- foreach::`%dopar%`
      parmx=foreach::foreach(k = 1:nstarts,.export=c("W","SF.mplus3D",
                                            "CLF",'eps.alignment',
                                            paste0('transformEstimates.',
                                                   estimator),
                                            'clf.ignore.quantile'),
                    .packages=c("abind",'purrr','dplyr','tibble')) %dopar% {
                      set.seed(k)
                      if(hyper.first=='means'){
                        outM=stats::optim(par=stats::rnorm(ngroups-1,0,5),
                                   fn=SF.mplus3D,method="L-BFGS-B",
                                   est=stacked,comb=comb,nobs=n,
                                   estimator=estimator,
                                   eps.alignment=eps.alignment,
                                   clf.ignore.quantile=clf.ignore.quantile,
                                   hyper='means',otherHyper=rep(1,ngroups-1),
                                   control=list(maxit=10000,trace=0),
                                   lower=rep(-Inf,ngroups-1))
                        outV=stats::optim(par=abs(stats::rnorm(ngroups-1,0,3)),
                                   fn=SF.mplus3D,method="L-BFGS-B",
                                   est=stacked,comb=comb,nobs=n,
                                   estimator=estimator,
                                   eps.alignment=eps.alignment,
                                   clf.ignore.quantile=clf.ignore.quantile,
                                   hyper='variances',otherHyper=outM$par,
                                   control=list(maxit=10000,trace=0),
                                   lower=rep(1e-6,ngroups-1))
                        c(value=outM$value+outV$value,
                          convergence=outM$convergence & outV$convergence,
                          par=c(outM$par,outV$par))
                      } else if(hyper.first=='variances') {
                        outV=stats::optim(par=abs(stats::rnorm(ngroups-1,0,3)),
                                   fn=SF.mplus3D,method="L-BFGS-B",
                                   est=stacked,comb=comb,nobs=n,
                                   estimator=estimator,
                                   eps.alignment=eps.alignment,
                                   clf.ignore.quantile=clf.ignore.quantile,
                                   hyper='variances',
                                   otherHyper=rep(0,ngroups-1),
                                   control=list(maxit=10000,trace=0),
                                   lower=rep(1e-6,ngroups-1))
                        outM=stats::optim(par=stats::rnorm(ngroups-1,0,5),
                                   fn=SF.mplus3D,method="L-BFGS-B",
                                   est=stacked,comb=comb,nobs=n,
                                   estimator=estimator,
                                   eps.alignment=eps.alignment,
                                   clf.ignore.quantile=clf.ignore.quantile,
                                   hyper='means',otherHyper=outV$par,
                                   control=list(maxit=10000,trace=0),
                                   lower=rep(-Inf,ngroups-1))

                        c(value=outM$value+outV$value,
                          convergence=outM$convergence & outV$convergence,
                          par=c(outM$par,outV$par))
                      } else if(hyper.first=='no') {
                        out=stats::optim(par=c(stats::rnorm(ngroups-1,0,5),
                                        abs(stats::rnorm(ngroups-1,0,3))),
                                  fn=SF.mplus3D,method="L-BFGS-B",
                                  est=stacked,comb=comb,nobs=n,
                                  estimator=estimator,
                                  eps.alignment=eps.alignment,hyper='both',
                                  control=list(maxit=10000,trace=0),
                                  lower=c(rep(-Inf,ngroups-1),
                                          rep(1e-6,ngroups-1)))
                        c(value=out$value,
                          convergence=out$convergence,
                          par=out$par)
                      }
                    }
      parallel::stopCluster(cl)
    } else {
      if(ncores>1)
        warning(
          paste0('ncores>1 but doParallel is not installed.',
                 'Parallel processing not used'))
      parmx=list()
      for(k in 1:nstarts) {
        if(hyper.first=='means'){
          outM=stats::optim(par=stats::rnorm(ngroups-1,0,5),
                     fn=SF.mplus3D,method="L-BFGS-B",
                     est=stacked,comb=comb,nobs=n,estimator=estimator,
                     eps.alignment=eps.alignment,
                     clf.ignore.quantile=clf.ignore.quantile,
                     hyper='means',otherHyper=rep(1,ngroups-1),
                     control=list(maxit=10000,trace=0),
                     lower=rep(-Inf,ngroups-1))
          outV=stats::optim(par=abs(stats::rnorm(ngroups-1,0,3)),
                     fn=SF.mplus3D,method="L-BFGS-B",
                     est=stacked,comb=comb,nobs=n,estimator=estimator,
                     eps.alignment=eps.alignment,
                     clf.ignore.quantile=clf.ignore.quantile,
                     hyper='variances',otherHyper=outM$par,
                     control=list(maxit=10000,trace=0),
                     lower=rep(1e-6,ngroups-1))
          parmx[[k]]=c(value=outM$value+outV$value,
                       convergence=outM$convergence & outV$convergence,
                       par=c(outM$par,outV$par))
        } else if(hyper.first=='variances') {
          outV=stats::optim(par=abs(stats::rnorm(ngroups-1,0,3)),
                     fn=SF.mplus3D,method="L-BFGS-B",
                     est=stacked,comb=comb,nobs=n,estimator=estimator,
                     eps.alignment=eps.alignment,
                     clf.ignore.quantile=clf.ignore.quantile,
                     hyper='variances',otherHyper=rep(0,ngroups-1),
                     control=list(maxit=10000,trace=0),
                     lower=rep(1e-6,ngroups-1))
          outM=stats::optim(par=stats::rnorm(ngroups-1,0,5),
                     fn=SF.mplus3D,method="L-BFGS-B",
                     est=stacked,comb=comb,nobs=n,estimator=estimator,
                     eps.alignment=eps.alignment,
                     clf.ignore.quantile=clf.ignore.quantile,
                     hyper='means',otherHyper=outV$par,
                     control=list(maxit=10000,trace=0),
                     lower=rep(-Inf,ngroups-1))
          parmx[[k]]=c(value=outM$value+outV$value,
                       convergence=outM$convergence & outV$convergence,
                       par=c(outM$par,outV$par))
        } else if(hyper.first=='no') {
          out=stats::optim(par=c(stats::rnorm(ngroups-1,0,5),abs(stats::rnorm(ngroups-1,0,3))),
                    fn=SF.mplus3D,method="L-BFGS-B",
                    est=stacked,comb=comb,nobs=n,estimator=estimator,
                    eps.alignment=eps.alignment,hyper='all',
                    clf.ignore.quantile=clf.ignore.quantile,
                    control=list(maxit=10000,trace=0),
                    lower=c(rep(-Inf,ngroups-1),rep(1e-6,ngroups-1)))
          parmx[[k]]=c(value=out$value,convergence=out$convergence,par=out$par)
        }
      }
    }
    parout=parmx%>%stats::setNames(1:nstarts)%>%dplyr::bind_cols()%>%t()

    #reject runs that failed
    failedruns=which(parout[,2]!=0)
    nFailedRuns=length(failedruns)
    if(length(failedruns)==nrow(parout))
      stop('No starting values converged. No aligned solution has been found.
Increase starting values or add more common items across groups.')
    if(length(failedruns)>0)parout=parout[-failedruns,]
    if(is.null(nrow(parout)))parout=matrix(parout,nrow=1)

    #name parameter list
    colnames(parout)=c("f","convergence",
                       paste0("M.",2:ngroups),
                       paste0("V.",2:ngroups))

    #populate
    #this line takes the lowest objective function value
    align.hyperpars=parout[order(parout[,1],decreasing=FALSE)[1],-c(1,2)]
    align.means=align.hyperpars[1:(ngroups-1)]
    align.variances=align.hyperpars[((ngroups-1)+1):((ngroups-1)*2)]
    #populate means, centering at zero to start
    align.means=c(0,align.means)
    align.variances=c(1/prod(align.variances),align.variances)
    #re-center at weighted averages for grand mean of zero
    if(center.means)align.means=align.means-stats::weighted.mean(align.means,unlist(n))

    #re-center such that weighted product is 1
    weighted.prod=exp(stats::weighted.mean(log(align.variances),unlist(n)))
    align.variances=align.variances/weighted.prod
    weighted.prod=exp(stats::weighted.mean(log(align.variances),unlist(n)))
    weighted.prod #nice

    #means are in SD units for first group; divide by SD of first group
    align.means=align.means/sqrt(align.variances[1])
    #variances just need to be divided by the first one
    align.variances=align.variances/align.variances[1]
  } else {
    align.means=0
    align.variances=1
    parout=NULL
    nFailedRuns=NA
  }

  #return means and variances, plus parout
  out=list(mv=mapply(function(x,y)c(mean=x,var=y),
                     as.list(align.means),
                     as.list(align.variances),SIMPLIFY=FALSE),
           parout=parout,nFailedRuns=nFailedRuns)
  return(out)
}
