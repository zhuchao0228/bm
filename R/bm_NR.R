bm.NR <- function(x,y, vce="oim",mu=NULL,
            start=NULL, noconstant=TRUE,
            control,family,offset=rep(0,nrow(x)),...)
{
  # cat("\nNR",start,"\n")
  # if(noconstant==F) {
  #   sta<<-start
  #   x.t<<-x
  #   y.t<<-y
  # }

  if(!is.null(mu)) {
    res.glm<-suppressWarnings(bm.glm.fit(x=x,y=y, family = family,
                    mustart=mu,control=control,offset=offset,...))
  }
  if(!is.null(start)) {
    res.glm<-suppressWarnings(bm.glm.fit(x=x,y=y, family = family,
                    start= start,control=control,offset=offset,...))
    # sat<<-start
    # if(is.list(res.glm)) print("&&&&&&&")
  }
  if(family$link=="log" & noconstant==F) {
    res.temp<-suppressWarnings(bm.glm.fit(x=x,y=y, family = family,offset=offset,
                    control=control,start= c(-1,rep(0,length(start)-1)),...))
    # print(res.temp$coefficients)
    if(res.temp$deviance<res.glm$deviance) res.glm<-res.temp
    # if(res.temp$deviance<res.glm$deviance) print(1) else print(2)
  }
  beta_old<-res.glm$coefficients
  # if(noconstant==F) cat("\n",res.glm$coefficients,"1\n")
  # if(noconstant==F)res1<<-res.glm
  dev_old<-res.glm$deviance
  # print(11111)
  # print(format(dev_old,nsmall=15))
  # print(beta_old)
  u <- res.glm$fitted.values
  for(i in 1:25)
  {
    res_NR<-tryCatch({
      bm.NR.fit(y=y,x=x,u=u,beta=beta_old,
                link=family$link,offset=offset)
    },
      error=function(err) {
        return(0)
    })
    # print(res_NR)
    if(is.list(res_NR)) {
      beta_new<-res_NR$beta
      u <- res_NR$u
      if(family$validmu(u)) dev_new<-bm.NR.dev(y=y,u=u)
    } else {
      beta_new<-beta_old
      dev_new<-dev_old
    }
    if(!family$validmu(u)) {
      res_retr<-tryCatch({
        bm.retrieve(y=y,x=x,beta_new = beta_new,
              beta_old = beta_old,offset=offset,link=family$link)
      },error=function(err) {
          return(err<-0)
      })
      if(is.list(res_retr)) {
        beta_new<-res_retr$beta
        u<-res_retr$u
        dev_new<-bm.NR.dev(y=y,u=u)
      } else {
        beta_new<-beta_old
        u <- res.glm$fitted.values
        dev_new<-dev_old
      }
    }
    if(abs(dev_new - dev_old)/(0.1 + dev_new)
      <= control$epsilon) {
      coef<-as.vector(beta_new)
      dev<-dev_new
      fitted.values<-u
      names(fitted.values)<-row.names(x)
      break
    } else {
      beta_old<-beta_new
      if(family$link=="log") {
        u <- as.vector(exp(as.matrix(x) %*% beta_new))+offset
      } else {
        u <- as.vector(as.matrix(x) %*% beta_new)+offset
      }
      dev_old<-dev_new
    }
  }
  names(coef)<-colnames(x)
  # print(222222)
  # cat("\nDev_new=",format(dev_new,nsmall=15),"\n")
  vcov<-tryCatch({
    bm.NR.vcov(x=x,y=y,u = u,vce=vce,
                link = family$link)
    }, error=function(err) {
      p <- res.glm$rank
      p1 <- 1L:p
      Qr <- res.glm$qr
      chol2inv(Qr$qr[p1,p1,drop=FALSE])
  })
  dimnames(vcov)<-list(colnames(x),colnames(x))
  # print(fitted.values)
  # print(coef)
  result <- list(coefficients=coef,deviance=dev,
                vcov=vcov, fitted.values=fitted.values,
                linear.predictors=ifelse(family$link=="identity",
                fitted.values,log(fitted.values)),
                null.deviance=res.glm$null.deviance)
  return(result)
}

bm.NR.fit<-function(y=NULL, x=NULL, u=NULL,beta=NULL,
                    link="identity",offset,vce="oim")
{
  if(link=="log") {
    # Establish gradian vactor.
    g_temp<-ifelse((y-u)==0, 0, as.vector((y - u) / (1 - u)))
  } else {
    # Establish gradian vactor.
    g_temp<-ifelse((y-u)==0, 0, as.vector((y - u) / (u*(1 - u))))
  }
  g <-t(x) %*% g_temp
  # Variance-covariance matrix.
  vcov<-bm.NR.vcov(x=x,y=y,u=u,vce=vce,link=link)
  # Update betas by Newton-Raphson.
  beta_new <- beta + as.vector(vcov %*% g)
  if(link=="log") {
    u <- as.vector(exp(as.matrix(x) %*% beta_new))+offset
  } else {
    u <- as.vector(as.matrix(x) %*% beta_new)+offset
  }
  res<-list(beta=beta_new, u=u)
  return(res)
}

# Deviance calculation
bm.NR.dev<-function(y=NULL,u=NULL)
{
  likelihood<-ifelse((1-y)==0, log(u), log(1-u))
  dev<--2*sum(likelihood)
  return(dev)
}

bm.retrieve<-function(y=NULL,x=NULL,link,
  beta_new=NULL,beta_old=NULL,offset)
{
  beta_new_init<-beta_new
  for(m in 1:21)
  {
    for(j in 1:100)
    {
      # print(beta_new)
      beta_new<-step_halving(beta_old = beta_old,beta_new=beta_new)
      if(link=="log") {
        u <- as.vector(exp(as.matrix(x) %*% beta_new))+offset
      } else {
        u <- as.vector(as.matrix(x) %*% beta_new)+offset
      }

      if(family$validmu(u)) break
    }

    if(!(family$validmu(u))) {
      u <- as.vector(as.matrix(x) %*% beta_new_init)+offset
      m1<-ifelse(floor(m/2)==0,1,floor(m/2))
      u[which(u>=1)]<-1-(10^-m1)
      # cat("\n",1-(10^-m1),10^-m1,"\n")
      u[which(u<=0)]<-10^-m1
      res_NR<-bm.NR.fit(y=y,x=x,u=u,beta=beta_new_init,
        offset=offset,link=link)
      beta_new<-res_NR$beta
      # print(beta_new)
      u <- res_NR$u
    }
    if(family$validmu(u)) break
    if(m==21 & !(family$validmu(u))) {
      stop("Fail to converge.")
    }
  }
  res<-list(u=u,beta=beta_new)
  return(res)
}

# Step halving method of R : take average between
# current estimator values and previous admissible
# estimator values.
step_halving <- function(beta_old = NULL, beta_new = NULL)
{
  # Check for fitted values outside domain.
  beta_new <- (beta_old+beta_new )/ 2
  return(beta_new)
}

bm.NR.vcov<-function(x=x,y=y, u=NULL,
                      vce="oim", link)
{
  if(tolower(vce)=="eim") {
    # Expected.
    if(link=="identity") {
      w <- 1/(u*(1 - u))
      # w <- ifelse(u*(1 - u)==0, 0, 1/(u*(1 - u)))
    } else {
      w <- ifelse(1-u==0, 0, u / (1 - u))
    }
  } else {
    # Observed Hessian matrix.
    if(link=="identity") {
      w<-ifelse((1-y)==0, 1/u^2, 1/(1-u)^2)
    } else {
      w<-ifelse(1-y==0, 0, u/(1 - u)^2)
    }
  }
  # Information matrix.
  I <- t(x) %*% diag(as.vector(w)) %*% as.matrix(x)
  # Variance-covariance matrix.
  vcov <- solve(I,tol = 1e-20)
  # vcov<-chol2inv(I)
  # vcov <- solve(I)
  return(vcov)
}
