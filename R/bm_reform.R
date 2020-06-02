bm.reform<-function(x,y, t = NULL, num_bvectors = NULL,family,
                    indicator = NULL,start,vce,control,bound,...)
{
  sub_x<-subset(x,select=-`(Intercept)`)
  rownames.sub_x<-row.names(sub_x)
  # t <- as.data.frame(bvectors)
  rownames.bv<-row.names(t)
  colnames.bv<-colnames(t)
  t.repa<-as.data.frame(matrix(NA,ncol=ncol(t),
                        nrow=nrow(t),byrow=T,
                        dimnames=dimnames(t)))
  for(r in 1:num_bvectors)
  {
    if(r == 1) {
      sub_x<-matrix(unlist(apply(sub_x,1,function(p) p-t[r,])),ncol=ncol(sub_x),byrow=T)
      dimnames(sub_x)<-list(rownames.sub_x,colnames.bv)
      t.repa[r,]<-temp<-t[r, ]
      t<-matrix(unlist(apply(t,1,function(p) p-temp)),ncol=ncol(t),byrow=T)
      dimnames(t)<-list(rownames.bv,colnames.bv)
    } else {
      q <- r-1
      if(t[r,q]!=0) {
        temp <- t[r, ]/t[r, q]
      } else {
        # Find nonzero member in current vector
        nonzero<-which(!t[r,1:ncol(t)]%in%0)
        # print(nonzero)
        nonzero.col<-nonzero[1]
        # Switch the position of current column with the nonzero column
        # print(colnames.bv)
        name.temp<-colnames.bv[q]
        colnames.bv[q]<-colnames.bv[nonzero.col]
        colnames.bv[nonzero.col]<-name.temp
        # print(colnames.bv)
        sub_x<-sub_x[,colnames.bv]
        t<-t[,colnames.bv]

        t.repa<-t.repa[,colnames.bv]
        temp <- t[r, ]/t[r, q]
      }
      t.repa[r,]<-temp
      sub_x<-matrix(unlist(apply(sub_x,1,function(p) p-temp*p[q])),ncol=ncol(sub_x),byrow=T)
      dimnames(sub_x)<-list(rownames.sub_x,colnames.bv)
      t<-matrix(unlist(apply(t,1,function(p) p-temp*p[q])),ncol=ncol(t),byrow=T)
      dimnames(t)<-list(rownames.bv,colnames.bv)
    }
  }
  temp_x <- cbind(indicator=indicator, sub_x)
  temp_x <- temp_x[-which(temp_x[,"indicator"]==1), ]
  y<-y[-which(indicator==1)]
  sub_x<-subset(temp_x,select=-c(1:num_bvectors))
  # cat("\nreform",start,"\n")
  results_temp<-tryCatch({
      res<-bm.NR(x=sub_x,y=y, vce=vce,control=control,
                  family=family,start=start[colnames(sub_x)],
                  offset=rep(ifelse(family$link=="identity" &&
                  bound==1,1,0),nrow(sub_x)),noconstant=T,...)
    },
      error=function(err){
        return(ind=0)
  })
  # print(results_temp$coefficients)
  if(is.list(results_temp)) {
    vcov<-results_temp$vcov
    beta <- results_temp$coefficients
    res <- list(vcov = vcov, beta = beta, t.repa=t.repa,
                deviance=results_temp$deviance,fitted.values=results_temp$fitted.values)
    return(res)
  } else {
    temp_init<-bm.InitValue(x=sub_x,y=y,family=family,
                            offset=rep(ifelse(family$link=="identity" && bound==1,1,0),
                            nrow(sub_x)),noconstant=T,control=control)
    # print(sub_x%*%temp_init+1)
    if(length(temp_init)==nrow(sub_x)) {
      results_temp<-bm.NR(x=sub_x,y=y,mu=temp_init,
                          family=family,vce=vce,control=control,
                          offset=rep(ifelse(family$link=="identity" && bound==1,1,0),
                          nrow(sub_x)),noconstant=T,...)
    } else {
      results_temp<-bm.NR(x=sub_x,y=y,start=temp_init,
                          family=family,vce=vce,control=control,
                          offset=rep(ifelse(family$link=="identity" && bound==1,1,0),
                          nrow(sub_x)),noconstant=T,...)
    }
    # if(is.list(results_temp)) print("$$$$$$$$$")
    vcov <- results_temp$vcov
    beta <- results_temp$coefficients
    # print(results_temp$fitted.values)
    res <- list(vcov = vcov, beta = beta,t.repa=t.repa,
                deviance=results_temp$deviance,fitted.values=results_temp$fitted.values)
    return(res)
  }
}
