#' Multiple-Group Factor Analysis Alignment from \code{mirt} or \code{lavaan}
#'
#' Performs alignment (\url{https://www.statmodel.com/Alignment.shtml}) using
#' single-group models estimated in mirt or lavaan.
#'
#' @param fitList A \code{list} of fitted model objects. Currently only works
#' for single-group, unidimensional or bifactor models with no covariates
#' estimated in \code{mirt} or \code{lavaan}.
#' @param estimator The model type used, either \code{mirt.grm} for the graded
#' response model estimated in \code{mirt} or \code{lavaan.ordered} for the
#' categorical factor analysis model applied by \code{lavaan} when the
#' \code{ordered} input includes all variables in the model.
#' @param eps.alignment A numeric scalar for the alignment simplicity function,
#' given by (Asparouhov & Muthén, 2014, \emph{Structural Equation Modeling}):
#'
#' \deqn{\sqrt{\sqrt{x^2+\epsilon}}}
#'
#' where $x$ is the difference between corresponding estimates in each pair of
#' aligned models. Lower values may cause numerical instability; default 0.01
#' @param clf.ignore.quantile Another protection from numerical instability;
#' CLF values less than the \code{clf.ignore.quantile} of the full set of CLF
#' values are ignored when calculating the complexity function at each step.
#' Default 0.1 for removing the lowest 10% of CLF values.
#' @param bifactor.marginal A logical scalar indicating whether, for bifactor
#' models, alignment should take place on the marginal, rather than conditional,
#' metric for slopes (Ip, 2010, \emph{Applied Psychological Measurement}).
#' @param hyper.first A string scalar denoting which hyperparameter to align
#' first. Asparouhov & Muthén (2014) align all parameters simultaneously
#' (\code{'no'}); \code{'variances'} (default) performs a two-step process,
#' first aligning variances, then aligning means conditional on variance
#' estimates from the first step. \code{'means'} does the reverse.
#' @param center.means A logical scalar. Alignment fixes the first group's mean
#' to zero to estimate the others. If \code{center.means} is \code{TRUE}
#' (default), aligned means and models are returned after subtracting the
#' weighted mean \code{stats::weighted.mean()} from all mean estimates, yielding
#' a (weighted) grand mean of zero. Variances are automatically rescaled such
#' that their weighted product (i.e., log of weighted mean of
#' \code{e^(variance)}) is 1.
#' @param ncores Number of processor cores to distribute alignment starts
#' across; on systems that support multicore processing, using additional cores
#' can speed up the alignment step by roughly a factor of the number of cores.
#' Defaults to 1 for no parallel processing. Requires the \code{doParallel}
#' package and defaults to parallel processing if not installed.
#'
#' @export
#'
#' @examples
#' #load data
#' library(mirt)
#' library(lavaan)
#' library(purrr)
#' library(tibble)
#' library(magrittr)
#' dat=expand.table(Bock1997)
#' #fit configural models
#' fit.mirt=mirt(dat,1,SE=TRUE)
#' fit.lavaan=cfa(model='G =~ Item.1+Item.2+Item.3',data=dat,
#'                ordered=c('Item.1','Item.2','Item.3'),
#'                std.lv=TRUE,parameterization='delta')
#' (fit.lavaan@ParTable)%>%tibble::as_tibble()%>%print(n=Inf)
#' #test stuff
#' tab=fit.lavaan@ParTable
#' tab$start[23]=3
#' tab$est[23]=3
#' fit.lavaan2=lavaan(tab,data=fit.lavaan@Data)
#'
#' #get estimates
#' est.mirt=getEstimates.mirt(fit.mirt,SE=TRUE,bifactor.marginal=FALSE)
#' est.lavaan=getEstimates.lavaan(fit.lavaan,SE=TRUE)
#'
#' #test transformations
#' newMean=10
#' newVar=2
#' test.mirt=transformEstimates.mirt.grm(newMean,newVar,est.mirt)
#' test.lavaan=transformEstimates.lavaan.ordered(
#'               newMean,newVar,est.lavaan,toCompare=FALSE)
#' #load and test equivalence
#' tfit.mirt=loadEstimates.mirt.grm(fit.mirt,newMean,newVar,newpars=test.mirt)
#' test.mirt=mirt::coef(fit.mirt)
#' test.mirt
#' tfit.lavaan=loadEstimates.lavaan.ordered(
#'               fit.lavaan,newMean,newVar,newpars=test.lavaan)
#' tfit.lavaan@ParTable%>%tibble::as_tibble()%>%print(n=Inf)
#' test.lavaan
#'
#' #now on stacked estimates
#' estList=list(est.mirt%>%purrr::imap(function(x,n){
#'   rownames(x)[2]=paste0(rownames(x)[2],'_ho')
#'   if(!n%in%c('a','se.a'))colnames(x)[2]=paste0(colnames(x)[2],'_ho')
#'   x
#' }),est.mirt%>%purrr::imap(function(x,n){
#'   rownames(x)[1]=paste0(rownames(x)[1],'_hi')
#'   if(!n%in%c('a','se.a'))colnames(x)[1]=paste0(colnames(x)[1],'_hi')
#'   x
#' }))
#' stack=stackEstimates(estList)
#' test.stack=transformEstimates.mirt.grm(c(0,0),c(1,1),stack)
#' sf.stack=SF.mplus3D(c(0,1),stack,combn(1:2,2),c(100,200),'mirt.grm',
#'                                        eps.alignment=0.01,
#'                                        clf.ignore.quantile=0.1)
#' test.stack2=transformEstimates.mirt.grm(c(0,1),c(1,1/2),stack)
#'
#' #try align?
#' #lavaan
#' set.seed(0)
#' sim.base=list(simdata(a=as.numeric(est.mirt$a),d=est.mirt$d,N=5000,
#'                       itemtype='graded',sigma=matrix(1),mu=0),
#'               simdata(a=as.numeric(est.mirt$a),d=est.mirt$d,N=5000,
#'                       itemtype='graded',sigma=matrix(2),mu=1))
#' fit.base=sim.base%>%purrr::map(~cfa(model="G =~ Item_1 + Item_2 + Item_3",
#'                              data=as.data.frame(.),
#'                              ordered=paste0('Item_',1:3),std.lv=TRUE,
#'                              parameterization='delta'))
#' fit.base%>%purrr::map(lavInspect,'est')%>%purrr::transpose()
#' est.base=purrr::map(fit.base,getEstimates.lavaan,SE=TRUE)
#' #not run: using parallel processes with ncores=3
#' # align.stack=align.optim(stackEstimates(est.base),c(100,200),nstarts=3,
#' #                         hyper.first='variances',ncores=3,
#' #                         eps.alignment=0.01,clf.ignore.quantile=0.1,
#' #                         estimator='lavaan.ordered',center.means=FALSE)
#' align.stack=align.optim(stackEstimates(est.base),c(100,200),nstarts=3,
#'                         hyper.first='variances',ncores=1,
#'                         eps.alignment=0.01,clf.ignore.quantile=0.1,
#'                         estimator='lavaan.ordered',center.means=FALSE)
#' align.stack
#' fit.align=Alignment(fit.base,'lavaan.ordered',center.means=FALSE)
#'
#' #mirt
#' fit.base2=list()
#' for(i in 1:length(sim.base)){
#'   fit.base2[[i]]=mirt(sim.base[[i]],1,'graded',SE=TRUE)
#' }
#' est.base2=purrr::map(fit.base2,getEstimates.mirt,SE=TRUE,bifactor.marginal=FALSE)
#' #not run: using parallel processes with ncores=3
#' # align.stack2=align.optim(stackEstimates(est.base2),c(100,200),nstarts=3,
#' #                          hyper.first='variances',ncores=3,
#' #                          eps.alignment=0.01,clf.ignore.quantile=0.1,
#' #                          estimator='mirt.grm',center.means=FALSE)
#' align.stack2=align.optim(stackEstimates(est.base2),c(100,200),nstarts=3,
#'                          hyper.first='variances',ncores=1,
#'                          eps.alignment=0.01,clf.ignore.quantile=0.1,
#'                          estimator='mirt.grm',center.means=FALSE)
#' align.stack2
#' fit.align2=Alignment(fit.base2,'mirt.grm',center.means=FALSE)
#'
#' #did it work?
#' fit.align$hypers
#' fit.align2$hypers
#' fit.align$est%>%purrr::transpose()%>%purrr::map(~mean(.[[1]]-.[[2]]))
#' fit.align2$est%>%purrr::transpose()%>%purrr::map(~mean(.[[1]]-.[[2]]))
#' fit.align$fit
#' fit.align2$fit
#' (fit.align$fit%>%purrr::map(~.@ParTable%>%tibble::as_tibble()%>%dplyr::filter(free!=0))%>%
#'   purrr::transpose())[c('start','est')]%>%purrr::map(~mean(.[[1]]-.[[2]]))
#' (fit.align2$fit%>%purrr::map(coef)%>%
#'     purrr::transpose())[paste0('Item_',1:3)]%>%purrr::map(~mean(.[[1]]-.[[2]]))
#' #appears so!
#'
Alignment=function(fitList,estimator,
                   eps.alignment=0.01,clf.ignore.quantile=0.1,
                   bifactor.marginal=FALSE,
                   hyper.first='variances',center.means=TRUE,
                   ncores=1){
  if(estimator=='mirt.grm'){
    #get all estimates
    est=fitList%>%purrr::map(getEstimates.mirt,SE=FALSE,
                      bifactor.marginal=bifactor.marginal)
    #align
    means.vars.parout=align.optim(est%>%stackEstimates,
                           n=fitList%>%purrr::map_dbl(~.@Data$N),
                           estimator=estimator,
                           eps.alignment=eps.alignment,
                           clf.ignore.quantile=clf.ignore.quantile,
                           hyper.first=hyper.first,
                           center.means=center.means,
                           ncores=ncores)
    means.vars=means.vars.parout$mv
    parout=means.vars.parout$parout
    #get aligned estimates
    test=purrr::map2(est,means.vars,
              ~transformEstimates.mirt.grm(.y[1],.y[2],.x))
    #fitted, aligned models
    tfit=list(fitList,test,means.vars)%>%purrr::pmap(
      function(x,y,z)loadEstimates.mirt.grm(x,z[1],z[2],y,do.fit=TRUE))
  } else if(estimator=='lavaan.ordered'){
    #get all estimates
    est=fitList%>%purrr::map(getEstimates.lavaan,SE=TRUE)
    #align
    means.vars.parout=align.optim(est%>%stackEstimates,
                           n=fitList%>%purrr::map_dbl(~.@Data@nobs[[1]]),
                           estimator=estimator,
                           eps.alignment=eps.alignment,
                           clf.ignore.quantile=clf.ignore.quantile,
                           hyper.first=hyper.first,
                           center.means=center.means,
                           ncores=ncores)
    means.vars=means.vars.parout$mv
    parout=means.vars.parout$parout
    #get aligned estimates
    test=purrr::map2(est,means.vars
              ,~transformEstimates.lavaan.ordered(.y[1],.y[2],.x,
                                                  toCompare=FALSE))
    #fitted, aligned models
    tfit=list(fitList,test,means.vars)%>%purrr::pmap(
      function(x,y,z)loadEstimates.lavaan.ordered(x,z[1],z[2],y,do.fit=TRUE))
  }
  names(means.vars)=names(test)
  #return stuff
  return(list(fit=tfit,est.og=est,est=test,hypers=means.vars,parout=parout))
}
