# source('C:/Users/Max/Dropbox/RPackages/AlignLV/equivalence_functions.R')

#load data
dat=expand.table(Bock1997)
#fit configural models
fit.mirt=mirt(dat,1,SE=T)
fit.lavaan=cfa(model='G =~ Item.1+Item.2+Item.3',data=dat,
               ordered=c('Item.1','Item.2','Item.3'),
               std.lv=T,parameterization='delta')
fit.lavaan@ParTable%>%as_tibble%>%print(n=Inf)
#test stuff
tab=fit.lavaan@ParTable
tab$start[23]=3
tab$est[23]=3
fit.lavaan2=lavaan(tab,data=fit.lavaan@Data)

#get estimates
est.mirt=getEstimates.mirt(fit.mirt,SE=T,bifactor.marginal=F)
est.lavaan=getEstimates.lavaan(fit.lavaan,SE=T)

#test transformations
newMean=10
newVar=2
test.mirt=transformEstimates.mirt.grm(newMean,newVar,est.mirt)
test.lavaan=transformEstimates.lavaan.ordered(newMean,newVar,est.lavaan,toCompare=F)
#load and test equivalence
tfit.mirt=loadEstimates.mirt.grm(fit.mirt,newMean,newVar,newpars=test.mirt)
coef(tfit.mirt)
test.mirt
tfit.lavaan=loadEstimates.lavaan.ordered(fit.lavaan,newMean,newVar,newpars=test.lavaan)
tfit.lavaan@ParTable%>%as_tibble%>%print(n=Inf)
test.lavaan

#now on stacked estimates
estList=list(est.mirt%>%imap(function(x,n){
  rownames(x)[2]=paste0(rownames(x)[2],'_ho')
  if(!n%in%c('a','se.a'))colnames(x)[2]=paste0(colnames(x)[2],'_ho')
  x
}),est.mirt%>%imap(function(x,n){
  rownames(x)[1]=paste0(rownames(x)[1],'_hi')
  if(!n%in%c('a','se.a'))colnames(x)[1]=paste0(colnames(x)[1],'_hi')
  x
}))
stack=stackEstimates(estList)
test.stack=transformEstimates.mirt.grm(c(0,0),c(1,1),stack)
sf.stack=SF.mplus3D(c(0,1),stack,combn(1:2,2),c(100,200),'mirt.grm')
test.stack2=transformEstimates.mirt.grm(c(0,1),c(1,1/2),stack)

#try align?
#lavaan
set.seed(0)
sim.base=list(simdata(a=as.numeric(est.mirt$a),d=est.mirt$d,N=5000,itemtype='graded',sigma=matrix(1),mu=0),
              simdata(a=as.numeric(est.mirt$a),d=est.mirt$d,N=5000,itemtype='graded',sigma=matrix(2),mu=1))
fit.base=sim.base%>%map(~cfa(model="G =~ Item_1 + Item_2 + Item_3",data=as.data.frame(.),
                             ordered=paste0('Item_',1:3),std.lv=T,parameterization='delta'))
fit.base%>%map(lavInspect,'est')%>%transpose
est.base=map(fit.base,getEstimates.lavaan,SE=T)
align.stack=align.optim(stackEstimates(est.base),c(100,200),nstarts=3,ncores=3,estimator='lavaan.ordered',center.means=F)
align.stack
fit.align=Alignment(fit.base,'lavaan.ordered',center.means=F)

#mirt
fit.base2=list()
for(i in 1:length(sim.base)){
  fit.base2[[i]]=mirt(sim.base[[i]],1,'graded',SE=T)
}
est.base2=map(fit.base2,getEstimates.mirt,SE=T,bifactor.marginal=F)
align.stack2=align.optim(stackEstimates(est.base2),c(100,200),nstarts=3,ncores=3,estimator='mirt.grm',center.means=F)
align.stack2
fit.align2=Alignment(fit.base2,'mirt.grm',center.means=F,parallel=F)

#did it work?
fit.align$hypers
fit.align2$hypers
fit.align$est%>%transpose%>%map(~mean(.[[1]]-.[[2]]))
fit.align2$est%>%transpose%>%map(~mean(.[[1]]-.[[2]]))
fit.align$fit
fit.align2$fit
(fit.align$fit%>%map(~.@ParTable%>%as_tibble%>%filter(free!=0))%>%
  transpose)[c('start','est')]%>%map(~mean(.[[1]]-.[[2]]))
(fit.align2$fit%>%map(coef)%>%
    transpose)[paste0('Item_',1:3)]%>%map(~mean(.[[1]]-.[[2]]))
#appears so!
