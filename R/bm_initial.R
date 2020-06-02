bm.InitValue <-function (x, y, weights = rep(1, nobs), start = NULL,
                          etastart = NULL, mustart = NULL, offset = rep(0, nobs),
                          family, control=bm.control(),noconstant=TRUE)
{
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  # ynames <- if(is.matrix(y)) rownames(y) else names(y)
  # conv <- FALSE
  nobs <- NROW(y)
  nvars <- ncol(x)
  ## define weights and offset if needed
  if (is.null(weights))
    weights <- rep.int(1, nobs)
  if (is.null(offset))
    offset <- rep.int(0, nobs)
  ## Initialize binomial functions:
  variance <- family$variance
  linkinv  <- family$linkinv
  mu.eta <- family$mu.eta
  validmu  <- function(mu) all(is.finite(mu)) && all(mu>0 & mu<1)
  if(is.null(mustart)) {
    ## calculates mustart and may change y and weights and set n (!)
    eval(family$initialize)
  } else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  coefold <- NULL
  eta <-
    if(!is.null(etastart)) etastart
    else if(!is.null(start))
      if (length(start) != nvars)
        stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", nvars, paste(deparse(xnames), collapse=", ")),
            domain = NA)
      else {
        coefold <- start ##Where the starting values were applied.
        offset + as.vector(if (NCOL(x) == 1L) x * start else x %*% start)
      }
    else family$linkfun(mustart)
    mu <- linkinv(eta)
    if (!(validmu(mu)))
      stop("cannot find valid starting values: please specify some", call. = FALSE)
    good <- weights > 0
    varmu <- variance(mu)[good]
    if (anyNA(varmu))
      stop("NAs in V(mu)")
    if (any(varmu == 0))
      stop("0s in V(mu)")
    mu.eta.val <- mu.eta(eta)
    if (any(is.na(mu.eta.val[good])))
      stop("NAs in d(mu)/d(eta)")
    ## drop observations for which w will be zero
    good <- (weights > 0) & (mu.eta.val != 0)
    z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
    # print(z)
    #
    w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
    # print(w)
    fit <- lm.fit(x=x[good, , drop = FALSE]*w, y=z*w, singular.ok=FALSE, tol=min(1e-07, control$epsilon/1000))
    coef<-fit$coefficients
    if(family$link=="identity") mu<-(x %*% fit$coefficients)
      else mu<-exp(x %*% fit$coefficients)
    # print(mu)
    res<-NULL
    if(any(mu<0 | mu>1)) {
      if(family$link=="identity") {
        if(noconstant==F) {
          # print(mu)
          # print(coef)
          coef[1]<-coef[1]-min(mu)
          coef<-0.9*(coef/(max(mu)-min(mu)))
          coef[1]<-coef[1]+0.01
        } else {
          coef<-coef*0.9/max(abs(mu))
        }
      } else {
        if(noconstant==F) {
          coef[1]<-coef[1]-log(max(mu))+log(0.9)
          # print(exp(x %*% coef))
        } else {
          conv<-0
          mu[which(mu>=1)]<-0.9999
          for (p in 1:25)
          {
            if(p==1) {
              ug<-mu
            } else {
              ug<-mu_new
            }
            z<-log(ug)+(y-ug)/ug
            w=ug/(1-ug)
            init_mu_new<-lm.wfit(x=x,y=z,w=w,singular.ok=FALSE, tol=min(1e-07, control$epsilon/1000))
            mu_new<-exp(init_mu_new$fitted.values)
            if(any(mu_new>1 | mu_new<0)) {
              mu_new[which(mu_new>=1)]<-0.9999
            } else {
              mustart<-ug
              conv<-1
              break
            }
            if(p==25 & conv!=1) conv<-0
          }
          res<-mustart
        }
      }
      # print(coef)
    }
    if(is.null(res)) res<-coef
    return(res)
}
