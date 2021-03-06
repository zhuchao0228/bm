
bm.glm.fit <-
  function (x, y, weights = rep(1, nobs), start = NULL,
            etastart = NULL, mustart = NULL, offset = rep(0, nobs),
            family = gaussian(), control = list(), intercept = TRUE)
  {
    control <- do.call("glm.control", control)
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if(is.matrix(y)) rownames(y) else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    ## define weights and offset if needed
    if (is.null(weights))
      weights <- rep.int(1, nobs)
    if (is.null(offset))
      offset <- rep.int(0, nobs)

    ## get family functions:
    variance <- family$variance
    linkinv  <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv) )
      stop("'family' argument seems not to be a valid family object", call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    unless.null <- function(x, if.null) if(is.null(x)) if.null else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu  <- unless.null(family$validmu,  function(mu) TRUE)
    if(is.null(mustart)) {
      ## calculates mustart and may change y and weights and set n (!)
      eval(family$initialize)
    } else {
      mukeep <- mustart
      eval(family$initialize)
      mustart <- mukeep
    }
    if(EMPTY) {
      eta <- rep.int(0, nobs) + offset
      if (!valideta(eta))
        stop("invalid linear predictor values in empty model", call. = FALSE)
      mu <- linkinv(eta)
      ## calculate initial deviance and coefficient
      if (!validmu(mu))
        stop("invalid fitted means in empty model", call. = FALSE)
      dev <- sum(dev.resids(y, mu, weights))
      w <- sqrt((weights * mu.eta(eta)^2)/variance(mu))
      residuals <- (y - mu)/mu.eta(eta)
      good <- rep_len(TRUE, length(residuals))
      boundary <- conv <- TRUE
      coef <- numeric()
      iter <- 0L
    } else {
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
      # cat(eta,"\n")
      if (!(validmu(mu) && valideta(eta)))
        stop("cannot find valid starting values: please specify some", call. = FALSE)
      ## calculate initial deviance and coefficient
      devold <- sum(dev.resids(y, mu, weights))
      boundary <- conv <- FALSE

      ##------------- THE Iteratively Reweighting L.S. iteration -----------
      for (iter in 1L:control$maxit) {
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

        if (all(!good)) {
          conv <- FALSE
          warning(gettextf("no observations informative at iteration %d", iter), domain = NA)
          break
        }
        z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
        w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
        fit <- lm.fit(x=x[good, , drop = FALSE]*w, y=z*w, singular.ok=FALSE, tol=min(1e-07, control$epsilon/1000))
        if (any(!is.finite(fit$coefficients))) {
          conv <- FALSE
          warning(gettextf("non-finite coefficients at iteration %d", iter), domain = NA)
          break
        }
        ## stop if not enough parameters
        if (nobs < fit$rank)
          stop(sprintf(ngettext(nobs,
                                "X matrix has rank %d, but only %d observation",
                                "X matrix has rank %d, but only %d observations"),
                       fit$rank, nobs), domain = NA)
        ## calculate updated values of eta and mu with the new coef:
        start[fit$qr$pivot] <- fit$coefficients
        eta <- drop(x %*% start)
        mu <- linkinv(eta <- eta + offset)
        dev <- sum(dev.resids(y, mu, weights))

        if (control$trace) cat("Deviance =", dev, "Iterations -", iter, "\n")

        ## check for divergence
        boundary <- FALSE
        if (!is.finite(dev)) {
          if(is.null(coefold))
            stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
          warning("step size truncated due to divergence", call. = FALSE)
          ii <- 1
          while (!is.finite(dev)) {
            if (ii > control$maxit)
              stop("inner loop 1; cannot correct step size", call. = FALSE)
            ii <- ii + 1
            start <- round((start + coefold)/2,digits = 10)
            eta <- drop(x %*% start)
            mu <- linkinv(eta <- eta + offset)
            dev <- sum(dev.resids(y, mu, weights))
          }
          boundary <- TRUE

          if (control$trace)
            cat("new deviance =", dev, "\n")
        }
        ## check for fitted values outside domain.
        if (!(valideta(eta) && validmu(mu))) {
          if(is.null(coefold))
            stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
          warning("step size truncated: out of bounds", call. = FALSE)
          ii <- 1
          while (!(valideta(eta) && validmu(mu))) {
            if (ii > control$maxit) {
              etaold <- drop(x %*% coefold)
              muold <- linkinv(etaold <- etaold + offset)
              if(validmu(muold)) start<-coefold
                else stop("inner loop 2; cannot correct step size", call. = FALSE)
            }
            ii <- ii + 1
            start <- round((start + coefold)/2,digits = 10)
            # print(format(start,nsmall=15))
            eta <- drop(x %*% start)
            mu <- linkinv(eta <- eta + offset)
          }
          boundary <- TRUE
          dev <- sum(dev.resids(y, mu, weights))
          if (control$trace)
            cat("new deviance =", dev, "\n")
        }

        ## check for convergence
        if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
          conv <- TRUE
          coef <- start
          break
        }
        else {
          devold <- dev
          coef <- coefold <- start
        }
      } ##-------------- end IRLS iteration -------------------------------

      if (!conv)
        warning("glm.fit: algorithm did not converge", call. = FALSE)
      if (boundary)
        warning("glm.fit: algorithm stopped at boundary value", call. = FALSE)
      eps <- 10*.Machine$double.eps
      if (family$family == "binomial") {
        if (any(mu > 1 - eps) || any(mu < eps))
          warning("glm.fit: fitted probabilities numerically 0 or 1 occurred", call. = FALSE)
      }
      if (family$family == "poisson") {
        if (any(mu < eps))
          warning("glm.fit: fitted rates numerically 0 occurred", call. = FALSE)
      }
      ## If X matrix was not full rank then columns were pivoted,
      ## hence we need to re-label the names ...
      ## Original code changed as suggested by BDR---give NA rather
      ## than 0 for non-estimable parameters
      if (fit$rank < nvars) coef[fit$qr$pivot][seq.int(fit$rank+1, nvars)] <- NA
      xxnames <- xnames[fit$qr$pivot]
      ## update by accurate calculation, including 0-weight cases.
      residuals <-  (y - mu)/mu.eta(eta)
      ##        residuals <- rep.int(NA, nobs)
      ##        residuals[good] <- z - (eta - offset)[good] # z does not have offset in.
      fit$qr$qr <- as.matrix(fit$qr$qr)
      nr <- min(sum(good), nvars)
      if (nr < nvars) {
        Rmat <- diag(nvars)
        Rmat[1L:nr, 1L:nvars] <- fit$qr$qr[1L:nr, 1L:nvars]
      }
      else Rmat <- fit$qr$qr[1L:nvars, 1L:nvars]
      Rmat <- as.matrix(Rmat)
      Rmat[row(Rmat) > col(Rmat)] <- 0
      names(coef) <- xnames
      colnames(fit$qr$qr) <- xxnames
      dimnames(Rmat) <- list(xxnames, xxnames)
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    # for compatibility with lm, which has a full-length weights vector
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    if(!EMPTY)
      names(fit$effects) <-
      c(xxnames[seq_len(fit$rank)], rep.int("", sum(good) - fit$rank))
    ## calculate null deviance -- corrected in glm() if offset and intercept
    wtdmu <-
      if (intercept) sum(weights * y)/sum(weights) else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    ## calculate df
    n.ok <- nobs - sum(weights==0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if(EMPTY) 0 else fit$rank
    resdf  <- n.ok - rank
    ## calculate AIC
    aic.model <-
      aic(y, n, mu, weights, dev) + 2*rank
    ##     ^^ is only initialize()d for "binomial" [yuck!]
    list(coefficients = coef, residuals = residuals, fitted.values = mu,
         effects = if(!EMPTY) fit$effects, R = if(!EMPTY) Rmat, rank = rank,
         qr = if(!EMPTY) structure(fit$qr[c("qr", "rank", "qraux", "pivot", "tol")], class = "qr"),
         family = family,
         linear.predictors = eta, deviance = dev, aic = aic.model,
         null.deviance = nulldev, iter = iter, weights = wt,
         prior.weights = weights, df.residual = resdf, df.null = nulldf,
         y = y, converged = conv, boundary = boundary)
  }
