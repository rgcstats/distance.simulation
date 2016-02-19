##############################################################################################################
##############################################################################################################
# Distance sampling simulations Jan 2016
#   N is defined throughout as population size divided by 1e3
#   This code is owned by Robert Clark (email robertgrahamclark@gmail.com).
#   Any re-use must be non-commercial and fully acknowledged, including a citation to
#   my PLOS-One paper "Statistical Efficiency in Distance Sampling", doi: 10.1371/journal.pone.0149298
##############################################################################################################
##############################################################################################################

#############################################################################################
# Functions
#############################################################################################

library(gtools)
library(numDeriv)
library(MASS)
library(stats)

rdetect <- function(d,g,pcond,num.obs, theta , beta1 , beta2 ,stop=F){
  if(stop) browser()
  probs <- detect.probs(d,g=g,pcond=pcond,num.obs=num.obs, theta=theta , beta1=beta1 , beta2=beta2 )
  randnums <- runif(length(d))
  colnames(probs)[1 + (randnums>probs[,1]) + (randnums>(probs[,1]+probs[,2])) + (randnums>(probs[,1]+probs[,2]+probs[,3])) ]
}


g.hr <- function(d,theta,w){ 1 - exp(-(d/theta[1])^(-theta[2]))}
# function to be used later with true value of theta2
g.hr1 <- function(d,theta,w,theta2){ 1 - exp(-(d/theta)^(-theta2))}
h1.hr <- function(d,theta,w){
  (g.hr(d,theta,w)-1) * (-1) * d / theta[1]^2 * theta[2] * (d/theta[1])^(-theta[2]-1)
}
h2.hr <- function(d,theta,w){
  (g.hr(d,theta,w)-1) * (d/theta[1])^(-theta[2]) * log(d/theta[1])
}


g.unif.cos1 <- function(d,theta,w){
  ( w^(-1) + theta*cos(pi*d/w) ) / ( w^(-1) + theta )
}

g.unif.cos2 <- function(d,theta,w){
  ( w^(-1) + theta[1]*cos(pi*d/w) + theta[2]*cos(2*pi*d/w) ) / ( w^(-1) + theta[1] + theta[2] )
}

g.halfnormal <- function(d,theta,w){exp(-d^2/(2*theta^2))}
g.halfnormal.m2 <- function(d,theta,w){
  y <- d / theta[1]
  exp(-y^2/2) * ( 1 + theta[2]*(y^4-6*y^2+3) ) / ( 1 + theta[2]*3 )
}
h.halfnormal <- function(d,theta,w){
  d^2/theta[1]^3*g.halfnormal(d,theta[1],w)
}

g.exp <- function(d,theta,w){exp(-d/theta)}
g.exp.m2 <- function(d,theta,w){
  exp(-d/theta[1]) * ( 1 + theta[2]*(d/theta[1])^2 )
}

gbar <- function(theta,g,w){
  integrate(f=g,lower=0,upper=w,theta=theta,rel.tol=1e-10,w=w)$value / w
}

l.dN <- function(par,d.s,ff,w,gg,stirling=T,stop=F,reparam=F){
  # if reparam=F, then par=(N,theta1,theta2)
  # if reparam=T, then par=(N*g(w/2),theta1,theta2)
  if(stop) browser()
  if(reparam) N <- par[1] / gg(d=w/2,theta=par[-1],w=w)
  if(!reparam) N <- par[1]
  theta <- par[-1]
  n <- length(d.s)
  gbar <- gbar(theta=theta,g=gg,w=w)
  if(!stirling) out <- lchoose(N*1e3,n) + sum(log(gg(d=d.s,theta=theta,w=w))) + (N*1e3-n)*log(1-ff*gbar)
  if(stirling) out <- N*1e3*log(N*1e3) - (N*1e3-n)*log(N*1e3-n) - n - lfactorial(n) + sum(log(gg(d=d.s,theta=theta,w=w))) + (N*1e3-n)*log(1-ff*gbar)
  out
}

l.d <- function(theta,d.s,ff,w,gg,stirling=T,stop=F){
  if(stop) browser()
  n <- length(d.s)
  gbar.value <- gbar(theta=theta,g=gg,w=w)
  N <- n / gbar.value / ff / 1e3
  if(!stirling) out <- lchoose(N*1e3,n) + sum(log(gg(d=d.s,theta=theta,w=w))) + (N*1e3-n)*log(1-ff*gbar.value)
  if(stirling) out <- N*1e3*log(N*1e3) - (N*1e3-n)*log(N*1e3-n) - n - lfactorial(n) + sum(log(gg(d=d.s,theta=theta,w=w))) + (N*1e3-n)*log(1-ff*gbar.value)
  out
}

l.d.firstpar <- function(theta1,theta2,d,ff,w,gg,stirling=T,stop=F){
  if(stop) browser()
  l.d(theta=c(theta1,theta2),d.s=d,ff=ff,w=w,gg=gg,stirling=stirling,stop=F)
}

mle.dN <- function(theta1min,theta1max,theta2grid,d.s,f,w,gg,stirling=T,stop=F){
  if(stop) browser()
  n <- length(d.s)
  if(missing(theta2grid)) theta <- optimize( f=l.d , interval=c(theta1min,theta1max) , d=d.s ,
                                             ff=f , w=w ,gg=gg , tol=1e-12 , maximum=T )$maximum
  if(!missing(theta2grid)){
    profile.theta2 <- rep(NA,length(theta2grid))
    theta1vals <- rep(NA,length(theta2grid))
    Nvals <- rep(NA,length(theta2grid))
    for(k in c(1:length(theta2grid))){
      optres <- optimize(f=l.d.firstpar,interval=c(theta1min,theta1max),theta2=theta2grid[k] ,
                         d=d.s , ff=f , w=w ,gg=gg , tol=1e-12 , maximum=T )
      profile.theta2[k] <- optres$objective
      theta1vals[k] <- optres$maximum
      Nvals[k] <- n / f / gbar(theta=c(theta1vals[k],theta2grid[k]),g=gg,w=w) / 1e3
    }
    theta2 <- theta2grid[which.max(profile.theta2)]
    theta1 <- theta1vals[which.max(profile.theta2)]
    theta <- c(theta1,theta2)
  }
  profile <- NULL
  if(!missing(theta2grid)) profile <- cbind( N=Nvals , theta1=theta1vals , theta2=theta2grid , loglik=profile.theta2  )
  gbar.value <- gbar(theta=theta,g=gg,w=w)
  N <- n / f / gbar.value / 1e3
  p <- 1 + length(theta)
  if(!missing(theta2grid)) p <- p - (length(theta2grid)<=1)
  AIC <- 2*(1+length(theta)) - 2*l.dN(par=c(N,theta),ff=f,w=w,gg=gg,d.s=d.s)
  varest.par <- matrix(data=NA,nrow=1+length(theta),ncol=1+length(theta))
  inv.varest.par <- -hessian(func=l.dN,x=c(N,theta),ff=f,w=w,gg=gg,d.s=d.s)
  try( varest.par <- solve( inv.varest.par ) , silent=T )
  list(par=c(N,theta),inv.varest.par=inv.varest.par,varest.par=varest.par,gbar=gbar.value,
       sd=sqrt(pmax(0,diag(varest.par))),
       cvpct=sqrt(pmax(0,diag(varest.par)))/c(N,theta)*100,AIC=AIC,profile=profile)
}


gbar.diffsq <- function(theta,gbar.value,g,w){
  ( gbar(theta=theta,g=g,w=w) - gbar.value )^2
}
gbar.diffsq.theta1 <- function(theta1,theta2,gbar.value,g,w){
  ( gbar(theta=c(theta1,theta2),g=g,w=w) - gbar.value )^2
}
gbar.diff <- function(theta1,theta2,gbar.value,g,w){
  ( gbar(theta=c(theta1,theta2),g=g,w=w) - gbar.value )
}

penalty <- function(g,h1,h2,w,f,theta,stop=F){
  if(stop) browser()
  gbar.value <- gbar(theta=theta,g=g,w=w)
  hbar1 <- integrate(f=h1,lower=0,upper=w,theta=theta,w=w)$value/w
  hbar2 <- 0
  if(length(theta)==2) hbar2 <- integrate(f=h2,lower=0,upper=w,theta=theta,w=w)$value/w
  prod <- function(d,g,h.i,h.j,theta,w){
    h.i(d=d,theta=theta,w=w) * h.j(d=d,theta=theta,w=w) /  g(d=d,theta=theta,w=w)
  }
  Delta <- matrix(data=0,nrow=2,ncol=2)
  Delta[1,1] <- integrate( f=prod , lower=0 , upper=w , g=g , h.i=h1 , h.j=h1 , theta=theta , w=w)$value / (w*gbar.value) - (hbar1/gbar.value)^2
  if(length(theta)==2){
    Delta[2,2] <- integrate( f=prod , lower=0 , upper=w , g=g , h.i=h2 , h.j=h2 , theta=theta , w=w)$value / (w*gbar.value) - (hbar2/gbar.value)^2
    Delta[1,2] <- Delta[2,1] <- integrate( f=prod , lower=0 , upper=w , g=g , h.i=h1 , h.j=h2 , theta=theta , w=w)$value / (w*gbar.value) - hbar1*hbar2/gbar.value^2
  }
  hbar <- matrix( c( hbar1 , hbar2 ) , nrow=2,ncol=1)
  penalty <- 1 + t(hbar) %*% ginv(Delta) %*% hbar / gbar.value^2 / (1-f*gbar.value)
  as.numeric(penalty)
}



P.vnarrow <- penalty(g=g.hr,h1=h1.hr,h2=h2.hr,theta=c(0.405,1.25),w=1,f=0.5,stop=F) # HR very narrow shoulder: penalty is 5.4 (very narrow shoulder)
P.narrow <- penalty(g=g.hr,h1=h1.hr,h2=h2.hr,theta=c(0.448,2),w=1,f=0.5,stop=F) # HR narrow shoulder: penalty is 3.1 (narrow shoulder)
P.wide <- penalty(g=g.hr,h1=h1.hr,h2=h2.hr,theta=c(0.484,3),w=1,f=0.5,stop=F) # HR wide shoulder: penalty is 2.2 (wide shoulder)
P.hn <- penalty(g=g.halfnormal,h1=h.halfnormal,h2=NA,theta=0.502,w=1,f=0.5,stop=F) # HN: penalty is 2.0
c(P.vnarrow,P.narrow,P.wide,P.hn)

P.vnarrow.f0.1 <- penalty(g=g.hr,h1=h1.hr,h2=h2.hr,theta=c(0.405,1.25),w=1,f=0.1,stop=F) # HR very narrow shoulder: penalty is 5.4 (very narrow shoulder)
P.narrow.f0.1 <- penalty(g=g.hr,h1=h1.hr,h2=h2.hr,theta=c(0.448,2),w=1,f=0.1,stop=F) # HR narrow shoulder: penalty is 3.1 (narrow shoulder)
P.wide.f0.1 <- penalty(g=g.hr,h1=h1.hr,h2=h2.hr,theta=c(0.484,3),w=1,f=0.1,stop=F) # HR wide shoulder: penalty is 2.2 (wide shoulder)
P.hn.f0.1 <- penalty(g=g.halfnormal,h1=h.halfnormal,h2=NA,theta=0.502,w=1,f=0.1,stop=F) # HN: penalty is 2.0
c(P.vnarrow.f0.1,P.narrow.f0.1,P.wide.f0.1,P.hn.f0.1)


#############################################################################################
# Plot some HR models
#############################################################################################

# calculate parameter values such that gbar is 0.6

optimize(f=gbar.diffsq.theta1,interval=c(0.05,4),gbar.value=0.6,w=1,theta2=2,g=g.hr) # theta1=0.448
optimize(f=gbar.diffsq.theta1,interval=c(0.05,4),gbar.value=0.6,w=1,theta2=3,g=g.hr) # theta1=0.484
optimize(f=gbar.diffsq.theta1,interval=c(0.05,4),gbar.value=0.6,w=1,theta2=0,g=g.halfnormal.m2) # theta1=0.502

dvals <- seq(from=0,to=1,length.out=501)
postscript(file="figure1.eps",width=6,height=6,onefile=FALSE,paper="a4",horizontal=F)
plotcols <- c("black","darkgray","black","black")
plotlty <- c(1,1,5,3)
par(oma=c(6,0,0,0))
plot(dvals,g.hr(dvals,theta=c(0.2,1)),pch="",xlim=c(0,1.06),ylim=c(0,1),xlab="Distance",ylab="Detection Probability",main="")
lines(dvals,g.hr(d=dvals,theta=c(0.405,1.25)),pch="",lty=plotlty[1],col=plotcols[1])
lines(dvals,g.hr(d=dvals,theta=c(0.448,2)),pch="",lty=plotlty[2],col=plotcols[2])
lines(dvals,g.hr(d=dvals,theta=c(0.48404,3)),pch="",lty=plotlty[3],col=plotcols[3])
lines(dvals,g.halfnormal(d=dvals,theta=0.502),pch="",lty=plotlty[4],col=plotcols[4])
legend("bottomleft",lty=plotlty,col=plotcols,
       legend=c(expression(paste("hazard rate very narrow shoulder: ",theta,"=(0.405,1.25)")),
                expression(paste("hazard rate narrow shoulder: ",theta,"=(0.448,2)")),
                expression(paste("hazard rate wide shoulder: ",theta,"=(0.484,3)")),
                expression(paste("half-normal: ",theta,"=0.502"))),cex=0.75)
text(cex=0.75,x=rep(1.05,4),y=c(0.275,0.19,0.14,0.095),
     labels=paste("F=",round(c(P.vnarrow,P.narrow,P.wide,P.hn),digits=1),sep=""))
par(oma=c(0,0,0,0))
dev.off()

# Create a table of asymptotic penalty P for hazard rate

shape.values <- c(1.1,1.25,1.5,2,2.5,3)
fvalues <-c(0.1,0.3,0.6,0.9)
gbar.values <- c(0.3,0.6,0.9)
fvalues.txt <- c("0.1","0.3","0.6","0.9")
gbar.values.txt <- c("0.3","0.6","0.9")
table.F <- matrix(data=NA,nrow=length(fvalues)*length(gbar.values),
                  ncol=length(shape.values)+2)
colnames(table.F) <- c("P","$\bar{g}$",paste("theta2=",shape.values,sep=""))
k <- 1
for(i1 in c(1:length(gbar.values))){
  for(i2 in c(1:length(fvalues))){
    fvalue <- fvalues[i2]
    gbar.value <- gbar.values[i1]
    table.F[k,1:2] <- c(fvalues.txt[i2],gbar.values.txt[i1])
    for(i3 in c(1:length(shape.values))){
      theta2 <- shape.values[i3]
      theta1 <- optimize(f=gbar.diffsq.theta1,interval=c(0.05,4),gbar.value=gbar.value,w=1,theta2=theta2,g=g.hr)$minimum
      Fval <- penalty(g=g.hr,h1=h1.hr,h2=h2.hr,theta=c(theta1,theta2),w=1,f=fvalue,stop=F)
      table.F[k,2+i3] <- format(round(Fval,digits=2),nsmall=2)
    }
    k <- k+1
  }
}
table.F
write.table(table.F,file="tableF.txt",sep=" & " , eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)


# Create a table of asymptotic penalty P for halfnormal

fvalues <- c(0.1,0.3,0.6,0.9)
gbar.values <- c(0.3,0.6,0.9)
fvalues.txt <- c("0.1","0.3","0.6","0.9")
gbar.values.txt <- c("0.3","0.6","0.9")
table.F.hn <- matrix(data=NA,nrow=length(fvalues)*length(gbar.values),
                     ncol=3)
colnames(table.F.hn) <- c("P","$\bar{g}$","Penalty $(F)$")
k <- 1
for(i1 in c(1:length(gbar.values))){
  for(i2 in c(1:length(fvalues))){
    fvalue <- fvalues[i2]
    gbar.value <- gbar.values[i1]
    table.F.hn[k,1:2] <- c(fvalues.txt[i2],gbar.values.txt[i1])
    theta1 <- optimize(f=gbar.diffsq,interval=c(0.025,5),gbar.value=gbar.value,w=1,g=g.halfnormal)$minimum
    Fval <- penalty(g=g.halfnormal,h1=h.halfnormal,theta=theta1,w=1,f=fvalue,stop=F)
    table.F.hn[k,3] <- format(round(Fval,digits=2),nsmall=2)
    k <- k+1
  }
}
table.F.hn
write.table(table.F.hn,file="tableFhn.txt",sep=" & " , eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)

############################################################################
# Simulate Various Options under DS
############################################################################

simulate.ds <- function(N,theta1,theta2,f=0.5,w=1,strip.wvals=c(1:100)/100,R=10,overdispersion=1,g=g.hr,stop=F){
  
  if(stop) browser()
  
  # Generate results object. For each method, record estimate of N and variance of Nhat for each simulation r=1,...,R
  method.names <- c("knowng","hr2","hn1","hn2","exp1","exp2",
                    paste("strip",strip.wvals,sep=""),
                    "model.average")
  quantity.names <- c("Nhat","Nc","Nc.hat","Vhat","AIC","n","gbar","theta1","theta2","posdef","condnum")
  sim.results <- array(data=NA,dim=c(R,length(method.names),length(quantity.names)),
                       dimnames=list(r=paste("rep",1:R,sep=""),
                                     method=method.names,
                                     quantity=quantity.names))
  
  # Now simulate r=1,...,R
  for(r in c(1:R)){
    
    # generate observed distances
    if(overdispersion==1){
      d.U <- runif(n=N,min=0,max=w/f)
    } else{
      shapepar <- (N-1)/(overdispersion-1)
      dvals <- (0.5+c(0:999))/1000 * w/f
      p.dvals <- as.vector( rdirichlet(n=1,alpha=rep(1/1000,1000)*shapepar))
      d.U <- sample(dvals,size=N,prob=p.dvals,replace=TRUE)
    }
    detect.prob.U <- g(d=d.U,theta=c(theta1,theta2)) * (d.U<=w)
    nums <- runif(N)
    detect <- 1*(nums<=detect.prob.U)
    d.s <- d.U[detect==1]
    n.s <- length(d.s)
    sim.results[r,,"n"] <- n.s
    sim.results[r,,"Nc"] <- sum(d.U<=w)
    
    # calculate each estimator of N and estimated variance
    
    for(method in c("hr2","hn1")){
      # nb second parameter of HR model is 1 or higher, following Buckland(1985)
      optim.results <- switch( method , 
                               hr2=optim( par=c(1,2),fn=l.d,d.s=d.s,ff=f,w=w,gg=g.hr,method="L-BFGS-B",
                                          lower=c(0.01,1),upper=c(1.5,6),control=list(fnscale=-1) ),
                               hn1=optim( par=1,fn=l.d,d.s=d.s,ff=f,w=w,gg=g.halfnormal,method="Brent",
                                          lower=0.01,upper=4,control=list(fnscale=-1) ) )
      gfit <- switch( method , hr2=g.hr , hn1=g.halfnormal )
      gbar.est <- gbar(theta=optim.results$par,g=gfit,w=w)
      Nhat <- n.s/f/gbar.est/1e3
      AIC <- 2*(1+length(optim.results$par)) - 2*l.dN(par=c(Nhat,optim.results$par),ff=f,w=w,gg=gfit,d.s=d.s)
      varest.par <- matrix(data=NA,nrow=1+length(optim.results$par),ncol=1+length(optim.results$par))
      inv.varest.par <- -hessian(func=l.dN,x=c(Nhat,optim.results$par),ff=f,w=w,gg=gfit,d.s=d.s)
      try( varest.par <- solve( inv.varest.par ) , silent=T )
      sim.results[r,method,"Nhat"] <- Nhat
      sim.results[r,method,"Vhat"] <- varest.par[1,1]
      sim.results[r,method,"AIC"] <- AIC
      sim.results[r,method,"theta1"] <- optim.results$par[1]
      if(length(optim.results$par)==2) sim.results[r,method,"theta2"] <- optim.results$par[2]
      if(!any(is.na(inv.varest.par))) sim.results[r,method,"posdef"] <- 1*all(eigen(inv.varest.par)$values>0)
      if(!any(is.na(inv.varest.par))) sim.results[r,method,"condnum"] <- kappa(inv.varest.par)
      sim.results[r,method,"gbar"] <- gbar.est
    }
    
    sim.results[r,"knowng","Nhat"] <- n.s / f / gbar(theta=c(theta1,theta2),g=g,w=w) / 1e3
    for(j in c(1:length(strip.wvals))){
      sim.results[r,paste("strip",strip.wvals[j],sep=""),"Nhat"] <- sum(d.s<=strip.wvals[j]) / (f*strip.wvals[j]/w) / 1e3
      sim.results[r,paste("strip",strip.wvals[j],sep=""),"Vhat"] <-
        sum(d.s<=strip.wvals[j]) / (f*strip.wvals[j]/w)^2 * (1-f*strip.wvals[j]/w) / 1e6
    }
    aic.scaled <- sim.results[r,,"AIC"] - min(sim.results[r,,"AIC"],na.rm=T)
    aic.wt <- exp(-0.5*aic.scaled[c("hr2","hn1")]) / sum(exp(-0.5*aic.scaled[c("hr2","hn1")]))
    sim.results[r,"model.average","Nhat"] <- sum( sim.results[r,c("hr2","hn1"),"Nhat"] * aic.wt )
    sim.results[r,"model.average","Vhat"] <- ( sum( aic.wt * sqrt( pmax(0,sim.results[r,c("hr2","hn1"),"Vhat"]) + 
                                                                     (sim.results[r,c("hr2","hn1"),"Nhat"]-sim.results[r,"model.average","Nhat"])^2 ) ) )^2
  }
  sim.results
}

t1 <- proc.time()[3]
R <- 10000
# all.simresults should be a data frame,
#  columns for N, f, sampsize, detfn, overdispersion, method, MSE, bias
all.simresults <- data.frame(N=integer(length=1),f=NA,sampsize=integer(length=1),
                             detfn="",overdispersion=NA,method="",
                             MSE=NA,VAR=NA,bias=NA,stringsAsFactors=F)
thetavals <- rbind( "vnarrow"=c(0.405,1.25) , "narrow"=c(0.448,2) , "wide"=c(0.484,3) , "hn"=c(0.502,0) )
set.seed(72449)
for(sampsize in c(50,100,200,300,400,500,600,700,800,900,1000)){
  if(sampsize %in% c(100,400)){
    fvals <- c(0.1,0.25,0.5) } else{
      fvals <- 0.1
    }
  for(f in fvals){
    for(detfn in c("vnarrow","narrow","wide","hn")){
      for(overdispersion in c(1,2,3)){
        set.seed(72449)
        N <- ceiling(sampsize/0.6/f)
        full.simresults.k <- simulate.ds(N=N,theta1=thetavals[detfn,1],
                                         theta2=thetavals[detfn,2],f=f,w=1,R=R,
                                         overdispersion=overdispersion,
                                         g=switch(detfn,"vnarrow"=g.hr,"narrow"=g.hr,
                                                  "wide"=g.hr,"hn"=g.halfnormal.m2) )
        biases <- apply(full.simresults.k[,,"Nhat"],2,mean) -N/1e3
        MSEs <- apply((full.simresults.k[,,"Nhat"]-N/1e3)^2,2,mean) 
        VARs <- apply(full.simresults.k[,,"Nhat"],2,var) 
        all.simresults <- rbind( all.simresults , 
                                 data.frame(N=N,f=f,sampsize=sampsize,detfn=detfn,
                                            overdispersion=overdispersion,
                                            method=names(biases),MSE=MSEs,VAR=VARs,bias=biases,stringsAsFactors=F) )
      }
    }
  }
}
t2 <- proc.time()[3]
(t2-t1)
(t2-t1)*10000/R/3600

all.simresults <- all.simresults[-1,]
all.simresults$bias.pct <- all.simresults$bias / (all.simresults$N/1e3) * 100
all.simresults$rrmse.pct <- sqrt(all.simresults$MSE) / (all.simresults$N/1e3) * 100
all.simresults$cv.pct <- sqrt(all.simresults$VAR) / (all.simresults$N/1e3) * 100
all.simresults$stripwidth <- as.numeric(substr(all.simresults$method,6,9))
all.simresults$stripwidth[substr(all.simresults$method,1,5)!="strip"] <- NA
all.simresults.knowng <- all.simresults[all.simresults$method=="knowng",]
all.simresults.knowng$MSE.knowng <- all.simresults.knowng$MSE
all.simresults <- merge( all.simresults ,
                         all.simresults.knowng[,c("MSE.knowng","sampsize","f","detfn","overdispersion")] ,
                         by=c("sampsize","f","detfn","overdispersion") )
all.simresults$F <- all.simresults$MSE / all.simresults$MSE.knowng

save.image(file="Post Simulation Image.rdata")
#load(file="Post Simulation Image.rdata")
all.simresults2 <- all.simresults

# Figure 2 (small sample and asymptotic F)

all.simresults3 <- all.simresults2[(all.simresults2$f==0.1)&(all.simresults2$overdispersion==1),]

postscript(file="figure2.eps",width=6,height=6,onefile=FALSE,paper="a4",horizontal=F)

all.simresults4 <- all.simresults3[(all.simresults3$detfn=="vnarrow")&(all.simresults3$method=="hr2"),]
all.simresults4 <- all.simresults4[order(all.simresults4$sampsize),]
plot( c(all.simresults4$sampsize,1250) , c(all.simresults4$F,P.vnarrow.f0.1) ,
      xlab="expected number of detections E(n)" ,
      ylab="Penalty (F)" , type="b" , lty=1 , pch=1 , ylim=c(0,5) ,
      yaxt="n", xaxt="n" )

all.simresults4 <- all.simresults3[(all.simresults3$detfn=="narrow")&(all.simresults3$method=="hr2"),]
all.simresults4 <- all.simresults4[order(all.simresults4$sampsize),]
points( c(all.simresults4$sampsize,1250) , c(all.simresults4$F,P.narrow.f0.1) , type="b" , lty=2 , pch=2 )

all.simresults4 <- all.simresults3[(all.simresults3$detfn=="wide")&(all.simresults3$method=="hr2"),]
all.simresults4 <- all.simresults4[order(all.simresults4$sampsize),]
points( c(all.simresults4$sampsize,1250) , c(all.simresults4$F,P.wide.f0.1) , type="b" , lty=3 , pch=3 )

all.simresults4 <- all.simresults3[(all.simresults3$detfn=="hn")&(all.simresults3$method=="hn1"),]
all.simresults4 <- all.simresults4[order(all.simresults4$sampsize),]
points( c(all.simresults4$sampsize,1250) , c(all.simresults4$F,P.hn.f0.1) , type="b" , lty=4 , pch=4 )

axis(side=1,at=c(c(1:10)*100,1250) , labels=c(c(1:10)*100,"asymptotic") , las=2 )
axis(side=2,at=seq(from=0,to=5,by=0.25),labels=FALSE )
axis(side=2,at=c(0:5),labels=c(0:5) , lwd.ticks=3)
axis(side=4,at=seq(from=0.5,to=6,by=0.25),labels=FALSE )
axis(side=4,at=c(0:5),labels=c(0:5) , lwd.ticks=3)

legend("bottomright",legend=c(expression(paste("hazard rate with very narrow shoulder ",theta,"=(0.405,1.25)")),
                              expression(paste("hazard rate with narrow shoulder ",theta,"=(0.448,2.00)")),
                              expression(paste("hazard rate with wide narrow shoulder ",theta,"=(0.484,3.00)")),
                              expression(paste("half-normal ",theta,"=0.502")) ) ,
       lty=c(1:4),pch=c(1:4) , title="true detection function" ,cex=0.6)

dev.off()

# Figure 3 (small sample F with c=2, no asymptotic available)

all.simresults3 <- all.simresults2[(all.simresults2$f==0.1)&(all.simresults2$overdispersion==2),]

postscript(file="figure3.eps",width=6,height=6,onefile=FALSE,paper="a4",horizontal=F)

all.simresults4 <- all.simresults3[(all.simresults3$detfn=="vnarrow")&(all.simresults3$method=="hr2"),]
all.simresults4 <- all.simresults4[order(all.simresults4$sampsize),]
plot( all.simresults4$sampsize , all.simresults4$F ,
      xlab="expected number of detections E(n)" ,
      ylab="Penalty (F)" , type="b" , lty=1 , pch=1 , ylim=c(0,6) ,
      yaxt="n", xaxt="n" )

all.simresults4 <- all.simresults3[(all.simresults3$detfn=="narrow")&(all.simresults3$method=="hr2"),]
all.simresults4 <- all.simresults4[order(all.simresults4$sampsize),]
points( all.simresults4$sampsize , all.simresults4$F , type="b" , lty=2 , pch=2 )

all.simresults4 <- all.simresults3[(all.simresults3$detfn=="wide")&(all.simresults3$method=="hr2"),]
all.simresults4 <- all.simresults4[order(all.simresults4$sampsize),]
points( all.simresults4$sampsize , all.simresults4$F , type="b" , lty=3 , pch=3 )

all.simresults4 <- all.simresults3[(all.simresults3$detfn=="hn")&(all.simresults3$method=="hn1"),]
all.simresults4 <- all.simresults4[order(all.simresults4$sampsize),]
points( all.simresults4$sampsize , all.simresults4$F , type="b" , lty=4 , pch=4 )

axis(side=1,at=c(1:10)*100 , labels=c(1:10)*100 , las=2 )
axis(side=2,at=seq(from=0,to=6,by=0.25),labels=FALSE )
axis(side=2,at=c(0:6),labels=c(0:6) , lwd.ticks=3)
axis(side=4,at=seq(from=0.5,to=6,by=0.25),labels=FALSE )
axis(side=4,at=c(0:6),labels=c(0:6) , lwd.ticks=3)

legend("bottomright",legend=c(expression(paste("hazard rate with very narrow shoulder ",theta,"=(0.405,1.25)")),
                              expression(paste("hazard rate with narrow shoulder ",theta,"=(0.448,2.00)")),
                              expression(paste("hazard rate with wide narrow shoulder ",theta,"=(0.484,3.00)")),
                              expression(paste("half-normal ",theta,"=0.502")) ) ,
       lty=c(1:4),pch=c(1:4) , title="true detection function" ,cex=0.6)

dev.off()




###############################
# Figures 4 - 7
###############################

# plot MSEs of strip and DS by strip width, for f=0.1, n=100,
#   with (detfn="hn",method="hn1" or "strip")
#        (detfn="vnarrow",method="hr2" or "strip")
#        (detfn="narrow",method="hr2" or "strip")
#        (detfn="wide",method="hr2" or "strip")

range(0.01*c(1:100)[
  all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==50)&(all.simresults$detfn=="narrow")&(all.simresults$overdispersion==2)&!is.na(all.simresults$stripwidth),"MSE"] <
    all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==50)&(all.simresults$detfn=="narrow")&(all.simresults$overdispersion==2)&(all.simresults$method=="hr2"),"MSE"] ] )

range(0.01*c(1:100)[
  all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==100)&(all.simresults$detfn=="narrow")&(all.simresults$overdispersion==2)&!is.na(all.simresults$stripwidth),"MSE"] <
    all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==100)&(all.simresults$detfn=="narrow")&(all.simresults$overdispersion==2)&(all.simresults$method=="hr2"),"MSE"] ] )

range(0.01*c(1:100)[
  all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==200)&(all.simresults$detfn=="narrow")&(all.simresults$overdispersion==2)&!is.na(all.simresults$stripwidth),"MSE"] <
    all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==200)&(all.simresults$detfn=="narrow")&(all.simresults$overdispersion==2)&(all.simresults$method=="hr2"),"MSE"] ] )

range(0.01*c(1:100)[
  all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==400)&(all.simresults$detfn=="narrow")&(all.simresults$overdispersion==2)&!is.na(all.simresults$stripwidth),"MSE"] <
    all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==400)&(all.simresults$detfn=="narrow")&(all.simresults$overdispersion==2)&(all.simresults$method=="hr2"),"MSE"] ] )

which.min(all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==50)&(all.simresults$detfn=="narrow")&(all.simresults$overdispersion==2)&!is.na(all.simresults$stripwidth),"MSE"])/100
which.min(all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==100)&(all.simresults$detfn=="narrow")&(all.simresults$overdispersion==2)&!is.na(all.simresults$stripwidth),"MSE"])/100
which.min(all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==200)&(all.simresults$detfn=="narrow")&(all.simresults$overdispersion==2)&!is.na(all.simresults$stripwidth),"MSE"])/100
which.min(all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==400)&(all.simresults$detfn=="narrow")&(all.simresults$overdispersion==2)&!is.na(all.simresults$stripwidth),"MSE"])/100


postscript(file="strip_comparison_n100_f0.1.eps",width=8,height=8,onefile=FALSE,paper="a4",horizontal=F)

all.simresults.subset <- all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==100),]
par(mfrow=c(2,2))
for(detfn in c("hn","vnarrow","narrow","wide")){
  if(detfn=="hn"){method <- "hn1" } else{method <- "hr2" }
  all.simresults.subset1 <- all.simresults.subset[(all.simresults.subset$detfn==detfn),]
  plot.header <- ""
  if(detfn=="hn") plot.header <- expression(paste("(a) half-normal ",theta,"=0.502"))
  if(detfn=="vnarrow") plot.header <- expression(paste("(b) hazard rate with very narrow shoulder ",theta,"=(0.405,1.25)"))
  if(detfn=="narrow") plot.header <- expression(paste("(c) hazard rate with narrow shoulder ",theta,"=(0.448,2.00)"))
  if(detfn=="wide") plot.header <- expression(paste("(d) hazard rate with wide shoulder ",theta,"=(0.484,3.00)"))
  plot( all.simresults.subset1$stripwidth , all.simresults.subset1$rrmse.pct ,pch="" , xlab="strip width for strip transect estimator" ,
        ylab="RRMSE(%)",ylim=c(0,50),main=plot.header,cex.main=1,cex=1,cex.lab=1,cex.axis=1)
  axis(side=2,at=seq(from=0,to=50,by=5),labels=FALSE,tick=T,tcl=-0.25)
  legend("bottomright",legend=c(expression(paste(hat(N)[ST]," with c=1",sep="")),expression(paste(hat(N)[ST]," with c=2",sep="")),expression(paste(hat(N)[ST]," with c=3",sep="")),
                                expression(paste(hat(N)[CDS]," with c=1",sep="")),expression(paste(hat(N)[CDS]," with c=2",sep="")),expression(paste(hat(N)[CDS]," with c=3",sep=""))),
         lty=c(1,2,3,1,2,3) , col=rep(c("blue","red"),each=3) , cex=0.6 , ncol=2 )
  for(overd in c(1,2,3)){
    all.simresults.subset.overd <- all.simresults.subset1[(all.simresults.subset1$overdispersion==overd),]
    lines(all.simresults.subset.overd$stripwidth,all.simresults.subset.overd$rrmse.pct,lty=overd,col="blue")
    abline(a=all.simresults.subset.overd$rrmse.pct[all.simresults.subset.overd$method==method],b=0,lty=overd,col="red")
  }
}

dev.off()

# now with n=50, f=0.1

all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==50)&(all.simresults$overdispersion==1)&
                 (all.simresults$detfn=="vnarrow")&(is.na(all.simresults$stripwidth)),]
all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==50)&(all.simresults$overdispersion==1)&
                 (all.simresults$method=="hr2"),]

all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==100)&(all.simresults$overdispersion==1)&
                 (all.simresults$method=="hr2"),]


all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==50)&(all.simresults$overdispersion==1)&
                 (all.simresults$method=="hr2"),]




postscript(file="strip_comparison_n50_f0.1.eps",width=8,height=8,onefile=FALSE,paper="a4",horizontal=F)

all.simresults.subset <- all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==50),]
par(mfrow=c(2,2))
for(detfn in c("hn","vnarrow","narrow","wide")){
  if(detfn=="hn"){method <- "hn1" } else{method <- "hr2" }
  all.simresults.subset1 <- all.simresults.subset[(all.simresults.subset$detfn==detfn),]
  plot.header <- ""
  if(detfn=="hn") plot.header <- expression(paste("(a) half-normal ",theta,"=0.502"))
  if(detfn=="vnarrow") plot.header <- expression(paste("(b) hazard rate with very narrow shoulder ",theta,"=(0.405,1.25)"))
  if(detfn=="narrow") plot.header <- expression(paste("(c) hazard rate with narrow shoulder ",theta,"=(0.448,2.00)"))
  if(detfn=="wide") plot.header <- expression(paste("(d) hazard rate with wide shoulder ",theta,"=(0.484,3.00)"))
  plot( all.simresults.subset1$stripwidth , all.simresults.subset1$rrmse.pct ,pch="" , xlab="strip width for strip transect estimator" ,
        ylab="RRMSE(%)",ylim=c(0,90),main=plot.header,cex.main=1,cex=1,cex.lab=1,cex.axis=1)
  axis(side=2,at=seq(from=0,to=90,by=5),labels=FALSE,tick=T,tcl=-0.25)
  legend("bottomright",legend=c(expression(paste(hat(N)[ST]," with c=1",sep="")),expression(paste(hat(N)[ST]," with c=2",sep="")),expression(paste(hat(N)[ST]," with c=3",sep="")),
                                expression(paste(hat(N)[CDS]," with c=1",sep="")),expression(paste(hat(N)[CDS]," with c=2",sep="")),expression(paste(hat(N)[CDS]," with c=3",sep=""))),
         lty=c(1,2,3,1,2,3) , col=rep(c("blue","red"),each=3) , cex=0.6 , ncol=2 )
  for(overd in c(1,2,3)){
    all.simresults.subset.overd <- all.simresults.subset1[(all.simresults.subset1$overdispersion==overd),]
    lines(all.simresults.subset.overd$stripwidth,all.simresults.subset.overd$rrmse.pct,lty=overd,col="blue")
    abline(a=all.simresults.subset.overd$rrmse.pct[all.simresults.subset.overd$method==method],b=0,lty=overd,col="red")
  }
}

dev.off()

# same again, but now with n=200 f=0.1

postscript(file="strip_comparison_n200_f0.1.eps",width=8,height=8,onefile=FALSE,paper="a4",horizontal=F)

all.simresults.subset <- all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==200),]
par(mfrow=c(2,2))
for(detfn in c("hn","vnarrow","narrow","wide")){
  if(detfn=="hn"){method <- "hn1" } else{method <- "hr2" }
  all.simresults.subset1 <- all.simresults.subset[(all.simresults.subset$detfn==detfn),]
  plot.header <- ""
  if(detfn=="hn") plot.header <- expression(paste("(a) half-normal ",theta,"=0.502"))
  if(detfn=="vnarrow") plot.header <- expression(paste("(b) hazard rate with very narrow shoulder ",theta,"=(0.405,1.25)"))
  if(detfn=="narrow") plot.header <- expression(paste("(c) hazard rate with narrow shoulder ",theta,"=(0.448,2.00)"))
  if(detfn=="wide") plot.header <- expression(paste("(d) hazard rate with wide shoulder ",theta,"=(0.484,3.00)"))
  plot( all.simresults.subset1$stripwidth , all.simresults.subset1$rrmse.pct ,pch="" , xlab="strip width for strip transect estimator" ,
        ylab="RRMSE(%)",ylim=c(0,35),main=plot.header,cex.main=1,cex=1,cex.lab=1,cex.axis=1)
  legend("bottomright",legend=c(expression(paste(hat(N)[ST]," with c=1",sep="")),expression(paste(hat(N)[ST]," with c=2",sep="")),expression(paste(hat(N)[ST]," with c=3",sep="")),
                                expression(paste(hat(N)[CDS]," with c=1",sep="")),expression(paste(hat(N)[CDS]," with c=2",sep="")),expression(paste(hat(N)[CDS]," with c=3",sep=""))),
         lty=c(1,2,3,1,2,3) , col=rep(c("blue","red"),each=3) , cex=0.6 , ncol=2 )
  for(overd in c(1,2,3)){
    all.simresults.subset.overd <- all.simresults.subset1[(all.simresults.subset1$overdispersion==overd),]
    lines(all.simresults.subset.overd$stripwidth,all.simresults.subset.overd$rrmse.pct,lty=overd,col="blue")
    abline(a=all.simresults.subset.overd$rrmse.pct[all.simresults.subset.overd$method==method],b=0,lty=overd,col="red")
  }
}

dev.off()

# same again, but now with n=400 f=0.1

postscript(file="strip_comparison_n400_f0.1.eps",width=8,height=8,onefile=FALSE,paper="a4",horizontal=F)

all.simresults.subset <- all.simresults[(all.simresults$f==0.1)&(all.simresults$sampsize==400),]
par(mfrow=c(2,2))
for(detfn in c("hn","vnarrow","narrow","wide")){
  if(detfn=="hn"){method <- "hn1" } else{method <- "hr2" }
  all.simresults.subset1 <- all.simresults.subset[(all.simresults.subset$detfn==detfn),]
  plot.header <- ""
  if(detfn=="hn") plot.header <- expression(paste("(a) half-normal ",theta,"=0.502"))
  if(detfn=="vnarrow") plot.header <- expression(paste("(b) hazard rate with very narrow shoulder ",theta,"=(0.405,1.25)"))
  if(detfn=="narrow") plot.header <- expression(paste("(c) hazard rate with narrow shoulder ",theta,"=(0.448,2.00)"))
  if(detfn=="wide") plot.header <- expression(paste("(d) hazard rate with wide shoulder ",theta,"=(0.484,3.00)"))
  plot( all.simresults.subset1$stripwidth , all.simresults.subset1$rrmse.pct ,pch="" , xlab="strip width for strip transect estimator" ,
        ylab="RRMSE(%)",ylim=c(0,25),main=plot.header,cex.main=1,cex=1,cex.lab=1,cex.axis=1)
  legend("bottomright",legend=c(expression(paste(hat(N)[ST]," with c=1",sep="")),expression(paste(hat(N)[ST]," with c=2",sep="")),expression(paste(hat(N)[ST]," with c=3",sep="")),
                                expression(paste(hat(N)[CDS]," with c=1",sep="")),expression(paste(hat(N)[CDS]," with c=2",sep="")),expression(paste(hat(N)[CDS]," with c=3",sep=""))),
         lty=c(1,2,3,1,2,3) , col=rep(c("blue","red"),each=3) , cex=0.6 , ncol=2 )
  for(overd in c(1,2,3)){
    all.simresults.subset.overd <- all.simresults.subset1[(all.simresults.subset1$overdispersion==overd),]
    lines(all.simresults.subset.overd$stripwidth,all.simresults.subset.overd$rrmse.pct,lty=overd,col="blue")
    abline(a=all.simresults.subset.overd$rrmse.pct[all.simresults.subset.overd$method==method],b=0,lty=overd,col="red")
  }
}

dev.off()

