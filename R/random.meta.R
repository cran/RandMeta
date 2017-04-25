
exact.meta=function(y, v, type="DL", B=500, theta.grd=NULL)
{K=length(y)
xsigntot=matrix(1, 2^K, K)
for(k in 1:K)
{xsigntot[,k]=rep(rep(c(-1,1), c(2^(K-k),2^(K-k))), 2^(k-1))}

if(length(theta.grd)>0)
{B=length(theta.grd)
pvalue=rep(0,B)}

if(length(theta.grd)==0)
{pvalue=rep(0, B)
theta.grd=seq(min(y)-(max(y)-min(y))/2/B, max(y)+(max(y)-min(y))/2/B, length=B)}

theta.grd=unique(sort(c(theta.grd, 0)))
B=length(theta.grd)
index0=(1:B)[theta.grd==0]

for(b in 1:B)
{tau2hat=max(0,   mean((y-theta.grd[b])^2-v))
if(type=="median") weight=1/sqrt(v+tau2hat)
if(type=="DL")     weight=abs(y-theta.grd[b])/(v+tau2hat)
if(type=="wang")   weight=pnorm(abs(y-theta.grd[b])/sqrt(v+tau2hat))-0.5

test=as.vector(rbind(sign(y-theta.grd[b]), xsigntot)%*%weight)
testp.obs=test[1]
testp.null=test[-1]
pvalue[b]=mean(abs(testp.null)>=abs(testp.obs))
}
ci.exact=range(theta.grd[pvalue>=0.05])
theta.exact=mean(theta.grd[which.max(pvalue)])
pvalue.exact=pvalue[index0]

result=list(theta.exact=theta.exact, pvalue.exact=pvalue.exact, ci.exact=ci.exact, pvtot.exact=pvalue, theta.grd=theta.grd)
return(result)
}

##################################################################################################################################
exactwilcox.meta=function(y, v, B=500)
{K=length(y)
xsigntot=matrix(1, 2^K, K)
for(k in 1:K)
{xsigntot[,k]=rep(rep(c(-1,1), c(2^(K-k),2^(K-k))), 2^(k-1))}
weight=1:K
testp.null=as.vector(xsigntot%*%weight)

theta.grd=seq(min(y)-(max(y)-min(y))/2/B, max(y)+(max(y)-min(y))/2/B, length=B)
theta.grd=unique(sort(c(theta.grd, 0)))
B=length(theta.grd)
index0=(1:B)[theta.grd==0]

pvalue=rep(0, B)
for(b in 1:B)
{tau2hat=max(0,   mean((y-theta.grd[b])^2-v))
weight=rank(abs(y-theta.grd[b])/sqrt(v+tau2hat))
testp.obs=sum(sign(y-theta.grd[b])*weight)
pvalue[b]=mean(abs(testp.null)>=abs(testp.obs))
}
ci.exact=range(theta.grd[pvalue>=0.05])
theta.exact=mean(theta.grd[which.max(pvalue)])
pvalue.exact=pvalue[index0]

result=list(theta.exact=theta.exact, pvalue.exact=pvalue.exact, ci.exact=ci.exact)
return(result)
}
###################################################################################################################################

fast.meta=function(y, v, type="median", B=500, theta.grd=NULL)
{K=length(y)
numtot=c(0,0,0,0,0,0,2,4,6,12,18,31,54,87,157,273,453,839,1472,2619)
num=numtot[K]

if(K==7)
  xsigntot=resmat7new
if(K==8)
  xsigntot=resmat8new
if(K==9)
  xsigntot=resmat9new
if(K==10)
  xsigntot=resmat10new
if(K==11)
  xsigntot=resmat11new
if(K==12)
  xsigntot=resmat12new
if(K==13)
  xsigntot=resmat13new
if(K==14)
  xsigntot=resmat14new
if(K==15)
  xsigntot=resmat15new
if(K==16)
  xsigntot=resmat16new
if(K==17)
  xsigntot=resmat17new
if(K==18)
  xsigntot=resmat18new
if(K==19)
  xsigntot=resmat19new
if(K==20)
  xsigntot=resmat20new

if(length(theta.grd)>0)
{B=length(theta.grd)
cov=rep(0,B)}

if(length(theta.grd)==0)
{cov=rep(0, B)
theta.grd=seq(min(y)-(max(y)-min(y))/2/B, max(y)+(max(y)-min(y))/2/B, length=B)}

for(b in 1:B)
{
  tau2hat=max(0, mean((y-theta.grd[b])^2-v))
  if(type=="median") weight=1/sqrt(v+tau2hat)
  if(type=="DL")     weight=abs(y-theta.grd[b])/(v+tau2hat)
  if(type=="wang")   weight=pnorm(abs(y-theta.grd[b])/sqrt(v+tau2hat))-0.5

  id=order(weight)
  weight.sort=weight[id]

  x=sign(y[id]-theta.grd[b])
  test=rbind(x, xsigntot)%*%(weight.sort)

  testp.obs=-abs(test[1])
  testp.null=test[-1]
  testp.null=unique(c(testp.null, testp.obs))
  cov[b]=(sum(testp.null<=testp.obs)+num)/2^K>0.025
}

if(max(cov)==1) ci.exact=range(theta.grd[cov==1])
if(max(cov)==0) ci.exact=NULL
return(list(ci.exact=ci.exact, cov=mean(cov)))
}

###########################################################################################################################

monte.meta=function(y, v, type="median", B=500, N=10000, theta.grd=NULL)
{K=length(y)
xsigntot=matrix(1, N, K)
for(k in 1:K)
{xsigntot[,k]=sample(c(-1,1), N, replace=T)}


if(length(theta.grd)>0)
{B=length(theta.grd)
pvalue=rep(0,B)}

if(length(theta.grd)==0)
{pvalue=rep(0, B)
theta.grd=seq(min(y)-(max(y)-min(y))/2/B, max(y)+(max(y)-min(y))/2/B, length=B)}

theta.grd=unique(sort(c(theta.grd, 0)))
B=length(theta.grd)
index0=(1:B)[theta.grd==0]

for(b in 1:B)
{
  tau2hat=max(0,  mean((y-theta.grd[b])^2-v))
  if(type=="median") weight=1/sqrt(v+tau2hat)
  if(type=="DL")     weight=abs(y-theta.grd[b])/(v+tau2hat)
  if(type=="wang")   weight=pnorm(abs(y-theta.grd[b])/sqrt(v+tau2hat))-0.5

  test=as.vector(rbind(sign(y-theta.grd[b]), xsigntot)%*%weight)
  testp.obs=test[1]
  testp.null=test[-1]
  pvalue[b]=mean(abs(testp.null)>=abs(testp.obs))
}

ci.appx=range(theta.grd[pvalue>=0.05])
theta.appx=mean(theta.grd[which.max(pvalue)])
pvalue.appx=pvalue[index0]

result=list(theta.appx=theta.appx, pvalue.appx=pvalue.appx, ci.appx=ci.appx, pvtot.appx=pvalue, theta.grd=theta.grd)
return(result)
}

##########################################################################################
montewilcox.meta=function(y, v, B=500, N=10000)
{K=length(y)
xsigntot=matrix(1, N, K)
for(k in 1:K)
{xsigntot[,k]=sample(c(-1,1), N, replace=T)}
weight=1:K
testp.null=as.vector(xsigntot%*%weight)

theta.grd=seq(min(y)-(max(y)-min(y))/2/B, max(y)+(max(y)-min(y))/2/B, length=B)
theta.grd=unique(sort(c(theta.grd, 0)))
B=length(theta.grd)
index0=(1:B)[theta.grd==0]

pvalue=rep(0, B)
for(b in 1:B)
{tau2hat=max(0,   mean((y-theta.grd[b])^2-v))
weight=rank(abs(y-theta.grd[b])/sqrt(v+tau2hat))
testp.obs=sum(sign(y-theta.grd[b])*weight)
pvalue[b]=mean(abs(testp.null)>=abs(testp.obs))
}
ci.appx=range(theta.grd[pvalue>=0.05])
theta.appx=mean(theta.grd[which.max(pvalue)])
pvalue.appx=pvalue[index0]

result=list(theta.appx=theta.appx, pvalue.appx=pvalue.appx, ci.appx=ci.appx)
return(result)
}
###########################################################################################################################

random.meta=function(y, v, type="DL", B=500, N=10000, Bstep=5, plot.meta=T)
{K=length(y)

if(type!="median" & type!="DL" & type!="wang" & type!="wilcox")
  print("type must be 'median, DL, wang' or 'wilcox' ")

if(type=="wilcox")
{if(K<6)
{print("This function is used only for meta analysis with more than 5 studies")}
  if(K>=6 && K<=20)
  {result.exact=exactwilcox.meta(y, v, B=B)}
  if(K>=21)
  {result.appx=montewilcox.meta(y, v, B=B, N=N)}
  if(K>=6 && K<=20)
  {print("Exact inference results")
    if(plot.meta==T)
    {plot(c(min(y-1.96*sqrt(v)), max(y+1.96*sqrt(v))), c(-1,K+1), type="n", xlab="theta",  ylab="study")
      for(k in 1:K)
      {lines( c(y[k]-1.96*sqrt(v[k]), y[k]+1.96*sqrt(v[k])), c(k, k))
        points( y[k], k)
      }
      points(result.exact$theta.exact, 0, pch=18, cex=2)
      lines(result.exact$ci.exact, c(0, 0), lwd=3)
      lines(rep(0, 2), c(-1, K+1), lty=2)
      lines(rep(result.exact$theta.exact,2), c(-1, K+1), lwd=2)
    }
    return(list(theta=result.exact$theta.exact, pvalue=result.exact$pvalue.exact, ci95=result.exact$ci.exact))}

  if(K>20)
  {print("Inference results from the Monte Carlo simulation")
    if(plot.meta==T)
    {plot(c(min(y-1.96*sqrt(v)), max(y+1.96*sqrt(v))), c(-1,K),  type="n", xlab="theta",  ylab="study")
      for(k in 1:K)
      {lines( c(y[k]-1.96*sqrt(v[k]), y[k]+1.96*sqrt(v[k])), c(k, k))
        points( y[k], k)
      }
      points(result.appx$theta.appx, 0, pch=18, cex=2)
      lines(result.appx$ci.appx, c(0, 0), lwd=3)
      lines(rep(0, 2), c(-1, K+1), lty=2)
      lines(rep(result.appx$theta.appx,2), c(-1, K+1), lwd=2)
    }
    return(list(theta=result.appx$theta.appx, pvalue=result.appx$pvalue.appx, ci95=result.appx$ci.appx))}
}

if(type=="median" | type=="DL" | type=="wang")
{if(K>=6 && K<=13)
{result.exact=exact.meta(y, v, type=type, B=B)}

  if(K>=14)
  {result.appx=monte.meta(y, v, type=type, B=B, N=N)
  grand.grd=result.appx$theta.grd}

  if(K>=14 &&  K<=20)
  {id.max=(1:B)[which.max(result.appx$pvtot.appx)]
  index.search=max(1, id.max-Bstep):min(B, id.max+Bstep)
  result.exact=exact.meta(y, v, type=type, theta.grd=grand.grd[index.search])
  while(result.exact$theta.exact==min(result.exact$theta.grd))
  {id.max=max(id.max-Bstep, 1)
  index.search=max(1, id.max-Bstep):min(B, id.max+Bstep)
  result.exact=exact.meta(y, v, type=type, theta.grd=grand.grd[index.search])}
  while(result.exact$theta.exact==max(result.exact$theta.grd))
  {id.max=min(id.max+Bstep, B)
  index.search=max(1, id.max-Bstep):min(B, id.max+Bstep)
  result.exact=exact.meta(y, v, type=type, theta.grd=grand.grd[index.search])}

  id.lower=min((1:B)[result.appx$pvtot.appx>=0.05])
  stop=0
  while(stop==0)
  {index.search=max(1, id.lower-Bstep):min(B, id.lower+Bstep)
  fit=fast.meta(y, v, type=type, theta.grd=grand.grd[index.search])
  if(fit$cov==1)
  {id.lower=max(id.lower-Bstep, 1)}
  if(fit$cov==0)
  {id.lower=min(B, id.lower+Bstep)}
  stop=1*(fit$cov>0 & fit$cov<1)
  }
  ci.lower.exact=fit$ci.exact[1]

  id.upper=max((1:B)[result.appx$pvtot.appx>=0.05])
  stop=0
  while(stop==0)
  {index.search=max(1, id.upper-Bstep):min(B, id.upper+Bstep)
  fit=fast.meta(y, v, type=type, B=B, theta.grd=grand.grd[index.search])
  if(fit$cov==1)
  {id.upper=min(id.upper+Bstep, B)}
  if(fit$cov==0)
  {id.upper=max(id.upper-Bstep, 1)}
  stop=1*(fit$cov>0 & fit$cov<1)
  }
  ci.upper.exact=fit$ci.exact[2]

  ci.exact=c(ci.lower.exact, ci.upper.exact)
  }

  if(K<6)
  {print("This function is used only for meta analysis with more than 5 studies")}

  if(K>=6  && K<=13)
  {print("Exact inference results")
    if(plot.meta==T)
    {plot(c(min(y-1.96*sqrt(v)), max(y+1.96*sqrt(v))), c(-1,K+1), type="n", xlab="theta",  ylab="study")
      for(k in 1:K)
      {lines( c(y[k]-1.96*sqrt(v[k]), y[k]+1.96*sqrt(v[k])), c(k, k))
        points( y[k], k)
      }
      points(result.exact$theta.exact, 0, pch=18, cex=2)
      lines(result.exact$ci.exact, c(0, 0), lwd=3)
      lines(rep(0, 2), c(-1, K+1), lty=2)
      lines(rep(result.exact$theta.exact,2), c(-1,K+1), lwd=2)
    }
    return(list(theta=result.exact$theta.exact, pvalue=result.exact$pvalue.exact, ci95=result.exact$ci.exact))}

  if(K>=14 && K<=20)
  {print("Exact inference results")
    if(plot.meta==T)
    {plot(c(min(y-1.96*sqrt(v)), max(y+1.96*sqrt(v))), c(-1,K+1), type="n", xlab="theta",  ylab="study")
      for(k in 1:K)
      {lines( c(y[k]-1.96*sqrt(v[k]), y[k]+1.96*sqrt(v[k])), c(k, k))
        points( y[k], k)
      }
      points(result.exact$theta.exact, 0, pch=18, cex=2)
      lines(ci.exact, c(0, 0), lwd=3)
      lines(rep(0, 2), c(-1, K+1), lty=2)
      lines(rep(result.exact$theta.exact,2), c(-1, K+1), lwd=2)
    }
    return(list(theta=result.exact$theta.exact, pvalue=result.exact$pvalue.exact, ci95=ci.exact))}

  if(K>20)
  {print("Inference results from the Monte Carlo simulation")
    if(plot.meta==T)
    {plot(c(min(y-1.96*sqrt(v)), max(y+1.96*sqrt(v))), c(-1,K),  type="n", xlab="theta",  ylab="study")
      for(k in 1:K)
      {lines( c(y[k]-1.96*sqrt(v[k]), y[k]+1.96*sqrt(v[k])), c(k, k))
        points( y[k], k)
      }
      points(result.appx$theta.appx, 0, pch=18, cex=2)
      lines(result.appx$ci.appx, c(0, 0), lwd=3)
      lines(rep(0, 2), c(-1, K+1), lty=2)
      lines(rep(result.appx$theta.appx,2), c(-1, K+1), lwd=2)
    }
    return(list(theta=result.appx$theta.appx, pvalue=result.appx$pvalue.appx, ci95=result.appx$ci.appx))}
}
}
