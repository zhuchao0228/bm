#' @name bm
#' @aliases bm
#' @title bm
#' @import stats utils
#' @keywords log binomial model, boundary, risk ratio/relative risk, convergence
#'
#' @description Fit a binomial regression model (bm) with a link function
#'     by Exact method.The boundary of parameter space in log and identity link
#'     binomial model causes trouble to locate a maximum likelihood solution
#'     and fails to converge when it is close to or on the boundary of space.
#'     Exact method utilizes the property of boundary vectors to re-parametrize
#'     the model without losing any information, and fits the model with no
#'     convergence issues.
#'
#' @usage bm(formula, data, contrasts = NULL, subset,na.action, link="identity",
#'     lfv.u=0.95, lfv.l=0.01, vce = "oim", rescode=NULL, control=bm.control(),\dots)
#'
#' @param formula an object of class "\link[stats]{formula}" (or one that
#'     can be coerced to that class): a symbolic description of the model
#'     to be fitted. The details of model specification are given under \sQuote{Details}.
#' @param data an optional data frame, list or environment (or object
#'     coercible by \code{as.data.frame} to a data frame) containing the variables
#'     in the model. If not found in \code{data}, the variables are taken from
#'     \code{environment(formula)}, typically the environment from which \code{bm}
#'     is called.
#' @param contrasts an optional list. See the \code{contrasts.arg} of
#'     \code{model.matrix.default}.
#' @param subset a specification of the rows to be used: defaults to all rows.
#'     This can be any valid indexing vector (see [\code{.data.frame}] for the
#'     rows of data or if that is not supplied, a data frame made up of the
#'     variables used in formula.
#' @param na.action a function which indicates what should happen when the data
#'     contain \code{NA}s. The default is set by the \code{na.action} setting of
#'     \code{\link[stats]{na.action}}, and is \code{\link[stats]{na.fail}} if that
#'     is unset. The 'factory-fresh' default is \code{\link[stats]{na.omit}}.
#'     Another possible value is \code{NULL}, no action. Value
#'     \code{\link[stats]{na.exclude}} can be useful.
#' @param link a specification for the binomial regression model link function. The
#'     \code{identity}, \code{log} and \code{logit} are supportted as aviliable link.
#'     The default link is \code{identity}.
#' @param lfv.u a probability criterion which decides the range of candidate vectors.
#'     The default value is 0.95, which means a covariate vectors with fitted probability
#'     from 0.95 to 1 included into boundary vector searching system as a candidate vector.
#' @param lfv.l a probability criterion which decides the range of candidate vectors.
#'     The default value is 0.01, which means a covariate vectors with fitted probability
#'     from 0.01 to 0 included into boundary vector searching system as a candidate vector.
#'     (This argument only works in the identity-link binomial regression model because the
#'     log binomial regression model does not have a lower bound of parameter space.)
#' @param vce the variance-covariance matrix to be used in fitting the model.
#'     Two options observed information matrix (OIM) and expected information
#'     matrix (EIM) could be choosen to calculate variance-covariance matrix.
#'     The default vce is \code{"OIM"}. This argument only works in the data with
#'     boundary vector. If there is no boundary vector included in the data, the results
#'     are directly from \code{glm}. In the \code{glm}, the standard error is calculated
#'     by expected information matrix.
#' @param rescode is a option to code the response variable if it is a factor.
#' @param control a list of parameters for controlling the fitting process.
#'     For bm.fit this is passed to bm.control.
#' @param \dots For \code{bm}: arguments to be used to form the default
#'     \code{control} argument if it is not supplied directly.
#'
#' @details A typical predictor has the form \code{response ~ terms} where
#'     response is the (numeric) response vector and terms is a series of
#'     terms which specifies a linear predictor for response. A terms
#'     specification of the form \code{first + second} indicates all the
#'     terms in first together with all the terms in second with any
#'     duplicates removed.

#'     A specification of the form \code{first:second} indicates the set of
#'     terms obtained by taking the interactions of all terms in first with
#'     all terms in second. The specification \code{first*second} indicates
#'     the cross of first and second. This is the same as
#'     \code{first + second + first:second}.

#'     The terms in the formula will be re-ordered so that main effects come
#'     first, followed by the interactions, all second-order, all third-order
#'     and so on: to avoid this pass a terms object as the formula.
#'
#' @return
#' \code{bm} returns an object of class inheriting from \code{"bm"} which
#' inherits from the class \code{"bm"}. The function \code{summary} (i.e.,
#' \code{summary.bm}) can be used to obtain or print a summary of the results
#' and the confidence interval of estimates. The argument \code{CF.lvl} in
#' \code{summary} represents the level of confidence interval claimed in the
#' model. The default value is \code{CF.lvl=0.95}. Optionally, Risk ratio estimates
#' and their related confidence interval are offered as an argument \code{RR} in
#' the \code{summary} for log-link function only. The default \code{RR=FALSE}
#' is not to display them.
#'
#' An object of class \code{"bm"} is a list containing at least the following
#' components:
#'
#' \item{coefficients}{a named vector of coefficients}
#' \item{residuals}{the working residuals, that is the residuals in the final
#' iteration of model.}
#' \item{fitted.values}{the fitted mean values, obtained by transforming the
#' linear predictors by the inverse of the log link function.}
#' \item{linear.predictors}{the linear fit on log scale.}
#' \item{deviance}{twice the absolute value of maximized log-likelihood.}
#' \item{aic}{A version of Akaike's An Information Criterion, minus twice the
#' maximized log-likelihood plus twice the number of parameters, computed by
#' the \code{aic} component of the family. For binomial model the dispersion
#' is fixed at one and the number of parameters is the number of coefficients.}
#' \item{null.deviance}{The deviance for the null model, comparable with
#' \code{deviance}. The null model will only include an intercept if there is
#' one in the model.}
#' \item{df.residual}{the residual degrees of freedom.}
#' \item{df.null}{the residual degrees of freedom for the null model.}
#' \item{response}{the response vector used in the mode.l}
#' \item{vcov}{the unscaled (\code{dispersion = 1}) estimated covariance matrix
#' of the estimated coefficients.}
#' \item{vce}{the type of standard error estimated based on the information matrix
#' (observed or espected) applied.}
#' \item{call}{the matched call.}
#' \item{na.action}{(where relevant) information returned by \code{stats::model.frame}
#' on the special handling of \code{NA}.}
#' \item{contrasts}{(where relevant) the contrasts used.}
#' \item{formula}{the formula supplied.}
#' \item{factor}{the order of factors used in response variable.}
#' \item{bvector}{the matrix of boundary vectors.}
#' \item{bv}{logical. Was the model judged to have boundary vector.}
#' \item{link}{Link function applied in the model}
#' \item{bound}{An indicator to describe the location of MLE on the boundary of
#' parameter space in the identity-link binomial model only.}
#'
#' @references
#'     Petersen, M. R. & Deddens, J. A. (2010). Maximum likelihood estimation
#'     of the log-binomial model. \emph{Communications in Statistics - Theory
#'     and Methods}, 39: 5, 874 - 883.
#'
#' @seealso
#' \code{glm}, \code{lm}, \code{lbm}.
#'
#' @examples
#' ## Two examples are from Petersen, M. R. & Deddens, J. A. (2010).
#'
#' ## Example 1.
#' x<-c(1:10)
#' y<-c(0,0,0,0,1,0,1,1,1,1)
#' data<-data.frame(x,y)
#' a<-bm(formula=y~x,data=data,link=log,vce=eim)
#'
#' ## Example 2.
#' x1<-c(1:11)
#' x2<-x1^2
#' y<-c(10,6,4,3,3,2,3,3,4,6,10)
#' dat<-cbind(x1,x2,y)
#' dat1<-apply(dat, 1, function(t) {
#'   temp<-data.frame(x1=rep(t[1],10),x2=rep(t[2],10),y=0)
#'   temp$y[1:t[3]]<-1
#'   return(temp)
#' })
#' data<-do.call(rbind, dat1)
#' a<-bm(formula=y~x1+x2,data=data)
#' summary(a)
#' a<-bm(formula=y~x1+x2,data=data,link=identity)
#' summary(a)
#'
#' @export bm

bm <- function(formula, data,contrasts = NULL,subset,na.action,link="identity",
                lfv.u=0.95,lfv.l=0.001,vce = "oim",rescode=NULL,control=bm.control(),...)
{
  scall <- match.call()
  cf<-match.call(expand.dots = FALSE)
  c<-match(c("formula", "data", "subset","na.action"), names(cf), 0L)
  # print(c)
  cf <- cf[c(1L, c)]
  cf$drop.unused.levels <- TRUE
  link<-substitute(link)
  if(!is.character(link)) link <- deparse(link)
  okLinks <- c("identity", "log")
  if (!(link %in% okLinks)) {
    stop(gettextf("link \"%s\" not available for binomial family; available links are %s",
        link, paste(sQuote(okLinks), collapse = ", ")), domain = NA)
  }
  vce<-substitute(vce)
  if(!is.character(vce)) vce <- deparse(vce)
  cf[[1L]] <- quote(stats::model.frame)
  # print(cf)
  data <- eval(cf, parent.frame())
  #data <- stats::model.frame(formula, data = data)
  na.action <- attr(data, "na.action")
  terms<-attr(data,"terms")
  y<-model.extract(data,"response")
  # names_observed<-attr(data,"names")[1]
  factor_index<-NULL
  if(is.factor(y)) {
    if(is.null(rescode))lvl<-levels(y) else lvl<-as.factor(rescode)
    if(length(lvl)!=2L) stop("Response variable has to be binary.")
    y_temp<-y
    y<-as.numeric(y)
    y[which(y_temp==lvl[1])]<-0
    y[which(y_temp==lvl[2])]<-1
    factor_index<-data.frame(Levels=lvl, Num=c(0,1))
  }
  x<-model.matrix(terms,data,contrasts)
  contrasts<-attr(x,"contrasts")
  family<-binomial(link=link)
  start<-bm.InitValue(x=x,y=y,family=family,control=control,noconstant=F)
  #start_temp<<-start
  # print(start)
  # if(any(is.na(ststart)==TRUE))
  #   stop("Please check the model. Some covariates may be 100 percent explained by others.
  #       \rIf interaction term is included, please reform the data and generate them and
  #       \rrerun the function.")
  res<-bm.NR(x=x,y=y,start=start,family=family,
              noconstant=F,control=control,...)
  # print(format(res$fitted.values,nsmall=15))
  # print(format(res$deviance,nsmall=15))
  x<-cbind(x,indicator.u=0,indicator.l=0,fitvalues=0)
  temp_null.deviance <- res$null.deviance
  deviance_temp<-res$deviance
  # x[,"indicator.u"][which(res$fitted.values > lfv & y!=0)]<-1
  # x[,"indicator.l"][which((1-res$fitted.values) > lfv & y!=1)]<-1
  x[,"indicator.u"][which(res$fitted.values > lfv.u & y==1)]<-1
  x[,"indicator.l"][which((res$fitted.values) < lfv.l & y==0)]<-1
  x[,"fitvalues"]<-res$fitted.values
  indicator<-data.frame(u=x[,"indicator.u"],l=x[,"indicator.l"])
  # print(indicator)
  start<-res$coefficients
  # cat("\nbm",start,"\n")
  b_vectors <- as.data.frame(x)[which(as.data.frame(x)$indicator.u==1 |
                                      as.data.frame(x)$indicator.l==1),]
  x<-subset.matrix(x,select=-c(indicator.u,indicator.l,fitvalues))
  # b_vectors<-subset.data.frame(b_vectors,select=-c(indicator.u,indicator.l))
  # print(b_vectors)
  # print(format(deviance_temp,nsmall=15))
  df.null <- nrow(x)-1
  df.residual <- nrow(x) - ncol(x)
  aic <- res$deviance + 2*ncol(x)
  res <- c(res, list(bv=0,aic = aic,vce=vce,bound=NULL,
          df.null = df.null,df.residual = df.residual))

  if(nrow(b_vectors)!=0) {

    # res.temp <-bm.fit(x=x,y=y, bvectors.old = b_vectors,control=control,family=family,
    #   num_bvectors = nrow(b_vectors), indicator.old = indicator,vce=vce,
    #   start=start, deviance=deviance_temp,null.dev=temp_null.deviance,...)

    res.temp <- tryCatch({
      suppressWarnings(bm.fit(x=x,y=y, bvectors.old = b_vectors,control=control,family=family,
              num_bvectors = nrow(b_vectors), indicator.old = indicator,vce=vce,
              start=start, deviance=deviance_temp,null.dev=temp_null.deviance,...))
    },
      error=function(err) {
        return(ind<-0)
    })
    if(is.list(res.temp))
      res <- c(res.temp, list(vce=vce, bv=1))
  }

  res <- c(res, list(factor=factor_index,formula = formula,call=scall,
          contrasts=contrasts,na.action=na.action,response=y, link=link))
  if(res$bv==0 & vce=="oim")
    warning("\rThe vce option works when the MLE is on the boundary of parameter space only.
            \rThe results from glm are reported for the non-boundary case.",call.=F)
  class(res)<-c("bm")
  return(res)
}

bm.fit<-function(x,y, bvectors.old = NULL,null.dev=NULL,vce,
  num_bvectors = NULL, indicator.old = NULL,control,
  start=NULL, deviance=NULL,family,...)
{
  num_var<-ncol(x)-1

  # print(bvectors.old)
  for(m in 1:2)
  {
    # if(m==2) break
    improved<-0
    if(m==1) {
      if(sum(bvectors.old$indicator.u)!=0L) cat("\nSearching MLE on the upper bound.\n")
      else {
        cat("\nNo candidate boundary vector closes to upper bound. \n")
        next
      }
    } else {
      if(family$link=="log") break
      if(sum(bvectors.old$indicator.l)!=0L) cat("\n\nSearching MLE on the lower bound.\n")
      else {
        cat("\nNo candidate boundary vector closes to lower bound. \n")
        next
      }
    }
    # m=1 for the data that MLE could be on the upper bound.
    # m=2 for the data that MLE could be on the lower bound (identity link only).
    if(m==1){
      if(sum(bvectors.old$indicator.l)!=0L) {
        bvectors<-bvectors.old[-which(bvectors.old$indicator.l==1),]
      } else bvectors<-bvectors.old
      bvectors<-subset.data.frame(bvectors,select=-c(indicator.u,indicator.l))
      # print(bvectors)
      indicator<-indicator.old$u
    } else {
      if(sum(bvectors.old$indicator.u)!=0) {
        bvectors<-bvectors.old[-which(bvectors.old$indicator.u==1),]
      } else bvectors<-bvectors.old
      bvectors<-subset.data.frame(bvectors,select=-c(indicator.u,indicator.l))
      bvectors$fitvalues<-1-bvectors$fitvalues
      indicator<-indicator.old$l
    }
    # print(bvectors)

    if(nrow(bvectors)>=1) {
      bvectors<-bvectors[order(bvectors$fitvalues, decreasing = T),]
      names<-row.names(bvectors)
      bvectors<-subset(bvectors, select=-c(`(Intercept)`,fitvalues))
      bvectors_1e<-floor(bvectors*1e+5)

      bvectors<-cbind(bvectors,ind_dup=cumsum(!duplicated(bvectors_1e)))
      bvectors_temp <- bvectors[!duplicated(bvectors_1e[,1:(ncol(bvectors_1e)-1)]), ]
      # print(bvectors_temp)
      max_num<-min(num_var,nrow(bvectors_temp))
      #if(m==2 && nrow(bvectors_temp)>max_num) bvectors_temp<-bvectors_temp[max_num,]
      for(j in 1:max_num)
      {
        cat("\nSearching for the admissible combinations of",j,"boundary vectors.\n")
        if(improved==0) ind_matrix<-t(utils::combn(nrow(bvectors_temp),j))
        else ind_matrix<-t(utils::combn(nrow(bvectors_temp)-bv_num,j-bv_num))
        for(i in 1:nrow(ind_matrix))
        {
          cat(".",sep="")
          if(i %% 10==0) cat("|",sep="")
          if(i %% 5==0 & i %% 10!=0) cat(" ",sep="")
          if(i %% 50==0) cat("\n",sep="")
          x_temp<-cbind(x, indicator=indicator)
          # print(x_temp)
          if(improved==0) temp_bvectors<-bvectors_temp[ind_matrix[i,],]
            else temp_bvectors<-rbind(bvectors_temp[1:bv_num,],
              bvectors_temp[(bv_num+1):nrow(bvectors_temp),][ind_matrix[i,],])
          x_temp[,"indicator"][!(row.names(x_temp) %in%
              names[which(sapply(bvectors$ind_dup,
                function(t){t %in% temp_bvectors$ind_dup}))])]<-0
          indicator_temp<-x_temp[,"indicator"]
          x_temp<-subset(x_temp, select=-indicator)
          temp_bvectors<-subset(temp_bvectors,select=-ind_dup)

          # print(temp_bvectors)

          conv<-tryCatch({
            bm.reform(x=x_temp,y=y, t = temp_bvectors,
              num_bvectors = nrow(temp_bvectors),bound=m,
              indicator=indicator_temp, start=start,
              vce=vce,control=control,family=family)
          },
            error=function(err){
              return(ind<-0)
            })
          # conv<-bm.reform(x=x_temp,y=y, t = temp_bvectors,
          #       num_bvectors = nrow(temp_bvectors),bound=m,
          #       indicator=indicator_temp, start=start,
          #       vce=vce,control=control,family=family)

          # if(is.list(conv)) print(1) else print(0)
          if(is.list(conv)) {
            # print(format(conv$deviance,nsmall=15))
            if(conv$deviance<=deviance) break
          }
          if(i==control$maxit) break
        }
        if(is.list(conv)) {
          if(conv$deviance<=deviance) {
            deviance<-conv$deviance
            res.temp<-conv
            # Sort the fitted probability by 1-p instead of p.
            if(m==2) conv$fitted.values<-1-conv$fitted.values
            b_vector<-temp_bvectors
            if(nrow(temp_bvectors)!=nrow(bvectors)) {
              temp_fit<-conv$fitted.values[row.names(bvectors_temp)[-which(row.names(bvectors_temp)%in%
                  row.names(temp_bvectors))]]
              # Reorder the boundary vectors left in the bvectors_temp by their fitted probabilities.
              temp_fit_name<-names(temp_fit[order(temp_fit,decreasing = T)])
              temp_fit_name<-c(row.names(temp_bvectors),temp_fit_name)
              # print(temp_fit_name)
              bvectors_temp<-bvectors_temp[temp_fit_name,]
            }
            bv_num<-j
            improved<-1
            bound<-m

            # print(format(deviance,nsmall=15))
            # print(format(res.temp$beta,nsmall=15))
            # print(b_vector)
            # start[colnames(res.temp$beta)]<-res.temp$beta
            cat("\nMinimum deviance improved. Deviance= ",deviance, "\n",sep = "")
          }
        }
      }
    }
    # if(!is.list(res.temp))
    #   if(m==1) stop("No admissible pairs of boundary vectors could be found on the upper bound.")
    # else stop("No admissible pairs of boundary vectors could be found on the lower bound.")
  }
  if(is.list(res.temp)) {
    res.temp<-bm.makeup(beta=res.temp$beta,
      vcov=res.temp$vcov,
      t.repa=res.temp$t.repa)
    res.temp$beta<-res.temp$beta[colnames(x)]
    if(family$link=="identity" && bound==1) res.temp$beta[1]<-1+res.temp$beta[1]
    res.temp$vcov<-res.temp$vcov[colnames(x),colnames(x)]
    eta <- x %*% res.temp$beta
    eta <- as.vector(round(eta,digits=15), mode = "numeric")
    if(family$link=="identity") mu<-eta else mu <- exp(eta)
    ## Working residuals
    residuals <-  (y - mu)/exp(eta)
    deviance<-bm.NR.dev(y=y,u=mu)
    df.null <- nrow(x) - 1
    df.residual <- nrow(x) - ncol(x)
    aic <- deviance + 2*ncol(x)
    res <- list(coefficients = res.temp$beta,bvector = b_vector,
      linear.predictors = eta, fitted.values = mu,
      null.deviance = null.dev, deviance = deviance,
      df.null = df.null, df.residual = df.residual,bound=bound,
      residuals=residuals,aic = aic,vcov = res.temp$vcov)
    return(res)
  } else {
    stop("No admissible combination of boundary vectors in the model.")
  }
}


bm.makeup<-function(beta,vcov,t.repa)
{
  t<-cbind(0,t.repa)
  name.vcov<-c("(Intercept)",names(t.repa))
  beta.complete<-rep(0,length(name.vcov))
  names(t)<-names(beta.complete)<-name.vcov
  vcov.complete<-matrix(0,nrow=length(name.vcov),
    ncol=length(name.vcov),
    dimnames=list(name.vcov,name.vcov))
  vcov.complete[row.names(vcov),colnames(vcov)]<-as.matrix(vcov)
  beta.complete[names(beta)]<-beta
  # Makeup estimates and variance-covariance matrix.
  for (i in nrow(t.repa):1)
  {
    cov.vector<-apply(vcov.complete,1,function(x) {
      cov<--sum(x*t[i,])
      return(cov)
    })
    var.temp<-sum(diag(t[i,])%*%vcov.complete%*%t(t[i,]))
    vcov.complete[i,]<-vcov.complete[,i]<-cov.vector
    vcov.complete[i,i]<-var.temp
    beta.complete[i]<--sum(t[i,]*beta.complete)
  }
  return(list(beta=beta.complete,
    vcov=vcov.complete))
}

#' @export
bm.control <- function(epsilon = 1e-8, maxit = 50)
{
  if(!is.numeric(epsilon) || epsilon <= 0)
    stop("value of 'epsilon' must be > 0")
  if(!is.numeric(maxit) || maxit <= 0)
    stop("maximum number of iterations must be > 0")
  list(epsilon = epsilon, maxit = maxit,trace=F)
}


#' @method summary bm
#' @export
summary.bm<-function(object, CF.lvl=0.95, RR=FALSE,...)
{
  if(0.5>=CF.lvl || CF.lvl>=1)
    stop("CF.lvl should be a value between 0.5 and 1.")
  alpha<-1-CF.lvl
  var.cf <- diag(object$vcov)
  coef<-object$coefficients
  ## calculate coef table
  s.err <- sqrt(var.cf)
  zvalue <- coef/s.err
  if(object$bv==1) {
    if(tolower(object$vce)=="oim"){
      dn<-c("Estimate", "Std.Err(OIM)")
    } else {
      dn<-c("Estimate", "Std.Err(EIM)")
    }
  } else {
    dn<-c("Estimate", "Std.Err(EIM)")
  }
  pvalue <- 2*pnorm(-abs(zvalue))
  coef.table <- cbind(coef, s.err, zvalue, pvalue)
  confint<-apply(coef.table[,1:2], 1, function(t){
    lower<-t[1]+qnorm(alpha/2)*t[2]
    upper<-t[1]+qnorm(1-alpha/2)*t[2]
    dat<-cbind(lower,upper)
  })
  coef.table<-cbind(coef.table,t(confint))
  ## Borrow idea from the Stata to present confidence interval.
  confint_name<-c(paste("[",(1-alpha)*100,"% ","Conf.",sep=""),"Interval]")
  dimnames(coef.table)<-list(names(object$coefficients),
    c(dn, "z value","Pr(>|z|)",confint_name))
  object$coefficients<-coef.table
  if(RR==TRUE) {
    if(object$link=="log"){
      RR.table<-exp(cbind(coef,t(confint)))
      dimnames(RR.table)<-list(row.names(object$coefficients),
        c("Risk Ratio", confint_name))
      object$RR<-RR.table
    } else stop("Relative risk (RR) only works for \"log\" link function.")
  }
  ans<-object
  class(ans) <- "summary.bm"
  return(ans)
}

printCoefmat.bm<-function (x, digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),signif.legend = signif.stars,
  dig.tst = max(1L, min(5L, digits-1L)), cs.ind = 1:k, tst.ind = k + 1,
  zap.ind = integer(),  P.values = NULL, has.Pvalue = nc >= 4L &&
    length(cn <- colnames(x)) && substr(cn[nc], 1L, 3L) %in% c("Pr(", "p-v"),
  eps.Pvalue = .Machine$double.eps, na.print = "NA", quote = FALSE, right = TRUE, ...)
{
  if (is.null(d <- dim(x)) || length(d) != 2L)
    stop("'x' must be coefficient matrix/data frame")
  nc <- d[2L]-2
  if (is.null(P.values)) {
    scp <- getOption("show.coef.Pvalues")
    if (!is.logical(scp) || is.na(scp)) {
      warning("option \"show.coef.Pvalues\" is invalid: assuming TRUE")
      scp <- TRUE
    }
    P.values <- has.Pvalue && scp
  } else if (P.values && !has.Pvalue)
    stop("'P.values' is TRUE, but 'has.Pvalue' is not")
  if (has.Pvalue && !P.values) {
    d <- dim(xm <- data.matrix(x[, -nc, drop = FALSE]))
    nc <- nc - 1
    has.Pvalue <- FALSE
  } else xm <- data.matrix(x)
  k <- nc - has.Pvalue - (if (missing(tst.ind))
    1
    else length(tst.ind))
  if (!missing(cs.ind) && length(cs.ind) > k)
    stop("wrong k / cs.ind")
  Cf <- array("", dim = d, dimnames = dimnames(xm))
  ok <- !(ina <- is.na(xm))
  for (i in zap.ind) xm[, i] <- zapsmall(xm[, i], digits)
  if (length(cs.ind)) {
    acs <- abs(coef.se <- xm[, cs.ind, drop = FALSE])
    if (any(ia <- is.finite(acs))) {
      digmin <- 1 + if (length(acs <- acs[ia & acs != 0]))
        floor(log10(range(acs[acs != 0], finite = TRUE)))
      else 0
      Cf[, cs.ind] <- format(round(coef.se, max(1L, digits -
          digmin)), digits = digits)
    }
  }
  if (length(tst.ind))
    Cf[, tst.ind] <- format(round(xm[, tst.ind], digits = dig.tst),
      digits = digits)
  if (any(r.ind <- !((1L:d[2L]) %in% c(cs.ind, tst.ind, if (has.Pvalue) nc))))
    for (i in which(r.ind)) Cf[, i] <- format(xm[, i], digits = digits)
  ok[, tst.ind] <- FALSE
  okP <- if (has.Pvalue)
    ok[, -nc]
  else ok
  x1 <- Cf[okP]
  dec <- getOption("OutDec")
  if (dec != ".")
    x1 <- chartr(dec, ".", x1)
  x0 <- (xm[okP] == 0) != (as.numeric(x1) == 0)
  if (length(not.both.0 <- which(x0 & !is.na(x0)))) {
    Cf[okP][not.both.0] <- format(xm[okP][not.both.0], digits = max(1L,
      digits - 1L))
  }
  if (any(ina))
    Cf[ina] <- na.print
  ## print(Cf)
  if (P.values) {
    if (!is.logical(signif.stars) || is.na(signif.stars)) {
      warning("option \"show.signif.stars\" is invalid: assuming TRUE")
      signif.stars <- TRUE
    }
    if (any(okP <- ok[, nc])) {
      pv <- as.vector(xm[, nc])
      Cf[okP, nc] <- format.pval(pv[okP], digits = dig.tst,
        eps = eps.Pvalue)
      ## signif.stars <- signif.stars && any(pv[okP] < 0.1)
      ## print(signif.stars)
      if (signif.stars) {
        Signif <- symnum(pv, corr = FALSE, na = FALSE,
          cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
          symbols = c("***", "**", "*", ".", " "))
        Cf <- cbind(Cf, format(Signif))
        ## print(Cf)
      }
    }
    else signif.stars <- FALSE
  }
  else signif.stars <- FALSE

  print.default(Cf, quote = quote, right = right, na.print = na.print,
    ...)
  if (signif.stars && signif.legend) {
    if ((w <- getOption("width")) < nchar(sleg <- attr(Signif,
      "legend")))
      sleg <- strwrap(sleg, width = w - 2, prefix = "  ")
    cat("---\nSignif. codes:  ", sleg, sep = "", fill = w +
        4 + max(nchar(sleg, "bytes") - nchar(sleg)))
  }
  invisible(x)
}


#' @method print summary.bm
#' @export
print.summary.bm <- function (x, digits = max(6L, getOption("digits")-3L),
  signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\nCall:\n")
  dput(x$call)
  cat("\nFormula:\n", paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n", sep = "")
  cat("\nLink-function: ",x$link,"\n",sep="")
  if(!is.null(x$bound)) {
    if(x$bound==1) cat("\nMLE is on the upper bound of parameter space.\n")
      else cat("\nMLE is on the lower bound of parameter space.\n")
  }
  cat("\nBoundary Vectors:\n")
  print(x$bvector, digits=digits)
  if(!is.null(x$factor)) {
    cat("\nIndex of response variable:\n")
    factor_names<-as.character(x$factor[,1])
    factor_names<-format(factor_names,
      width=max(nchar(factor_names)),
      justify="right")
    cat(factor_names[1],"$ ",x$factor[,2][1],"\n",
      factor_names[2],"$ ",x$factor[,2][2],"\n",sep="")
  }
  cat("\nCoefficients:\n")
  coefs <- x$coefficients
  printCoefmat.bm(coefs, digits=min(4L,digits-2L),
    signif.stars = signif.stars, na.print = "NA")
  if(!is.null(x$RR)) {
    cat("\nRisk Ratio Estimates:\n")
    print(x$RR, digits=digits)
  }
  cat("\n", apply(cbind(paste(format(c("Null","Residual"), justify="right"),
    "deviance:"), format(unlist(x[c("null.deviance","deviance")]),
      digits = max(5L, digits + 3L)), " on",
    format(unlist(x[c("df.null","df.residual")])),
    " degrees of freedom\n"), 1L, paste, collapse = " "), sep = "")
  if(nzchar(mess <- naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")
  cat("AIC: ", format(x$aic, digits = max(5L, digits + 3L)), "\n", sep = "")
  invisible(x)
}

#' @method print bm
#' @export
print.bm <- function(x, digits = max(3L, getOption("digits")-1L), ...)
{
  cat("\nCall:\n")
  dput(x$call)
  cat("\nFormula:\n", paste(deparse(x$formula),
    sep = "\n", collapse = "\n"), "\n", sep = "")
  cat("\nLink-function: ",x$link,"\n",sep="")
  if(!is.null(x$bound)) {
    if(x$bound==1) cat("\nMLE is on the upper bound of parameter space.\n")
    else cat("\nMLE is on the lower bound of parameter space.\n")
  }
  cat("\nBoundary Vectors:\n")
  print(x$bvector, digits=digits-2L)
  if(!is.null(x$factor)) {
    cat("\nIndex of response variable:\n")
    factor_names<-as.character(x$factor[,1])
    factor_names<-format(factor_names,
      width=max(nchar(factor_names)),
      justify="right")
    cat(factor_names[1],"$ ",x$factor[,2][1],"\n",
      factor_names[2],"$ ",x$factor[,2][2],"\n",sep="")
  }
  if(length(coef(x))) {
    cat("\nCoefficients")
    if(is.character(co <- x$contrasts))
      cat("  [contrasts: ",
        apply(cbind(names(co),co), 1L, paste, collapse = "="), "]")
    cat(":\n")
    print.default(format(x$coefficients, digits = digits),
      print.gap=2, quote = FALSE)
  } else cat("No coefficients\n\n")
  cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
    x$df.residual, "Residual\n")
  if(nzchar(mess <- naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")
  cat("     Null Deviance:",	format(signif(x$null.deviance,
    digits=max(5L, digits + 3L))),
    "\n Residual Deviance:", format(signif(x$deviance,
      digits=max(5L, digits + 3L))),
    "\tAIC:", format(signif(x$aic, digits=max(4L, digits + 3L))),"\n")
  invisible(x)
}

#' @method vcov bm
#' @export
vcov.bm <- function(object, ...) object$vcov

#' @method vcov summary.bm
#' @export
vcov.summary.bm <- function(object, ...)  object$vcov


#' @method logLik bm
#' @export
logLik.bm <- function(object, ...)  -object$deviance/2

#' @method logLik summary.bm
#' @export
logLik.summary.bm <- function(object, ...)  -object$deviance/2

