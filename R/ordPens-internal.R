coding <- function(x, constant=TRUE, splitcod=TRUE)
  {
  # Converts ordinal/categorical data into split-coding, or dummy-coding.

  # Input:
  # x: data matrix with ordinal/categorical data
  # constant: should a constant added in the regression

  # Output:
  # Split/dummy-coded data-matrix
  # Reference-Category = 1
  n <- nrow(x)
  kx <- apply(x,2,max) - 1
  xds <- matrix(0,n,sum(kx))

  # Loop to run through the columns of x
  for (j in 1:ncol(x))
		{
      j1col <- ifelse(j>1,sum(kx[1:(j-1)]),0)
  		# Loop to run through the rows of x
  		for (i in 1:n)
  			{
          if (x[i,j] > 1)
            {
              if (splitcod)
                xds[i,j1col + 1:(x[i,j]-1)] <- 1
              else
                xds[i,j1col + (x[i,j]-1)] <- 1
            }
  		  }
    }
  ## Output
  if (constant)
    return(cbind(1,xds))
  else
    return(xds)
 	}



genRidge <- function(x, y, offset, omega, lambda, model, delta=1e-6, maxit=25)
  {
    coefs <- matrix(0,ncol(x),length(lambda))
    fits <- matrix(NA,nrow(x),length(lambda))
    l <- 1
    for (lam in lambda)
      {
        if (model == "linear")
          {
            yw <- y - offset
            chollam <- chol(t(x)%*%x + lam*omega)
            coefs[,l] <- backsolve(chollam,
            backsolve(chollam, t(x)%*%yw, transpose=TRUE))
            fits[,l] <- x%*%coefs[,l] + offset
          }
        else
          {
            ## penalized logistic/poisson regression
            conv <- FALSE
    
            # start values
            if (l > 1)
              {
                b.start <- coefs[,l-1]
              }
            else
              {
                if (model == "logit")
                  {
                    gy <- rep(log(mean(y)/(1-mean(y))),length(y)) - offset
                  }
                else
                  {
                    gy <- rep(log(mean(y)),length(y)) - offset
                  }
                chollam <- chol(t(x)%*%x + lam*omega)
                b.start <- backsolve(chollam,
                backsolve(chollam, t(x)%*%gy, transpose=TRUE))
              }
    
            # fisher scoring
            b.old <- b.start
            i <- 1
            while(!conv)
              {
                eta <- x%*%b.old + offset
                if (model == "logit")
                  {
                    mu <- plogis(eta)
                    sigma <- as.numeric(mu*(1-mu))
                  }
                else
                  {
                    mu <- exp(eta)
                    sigma <- as.numeric(mu)
                  }
                score <- t(x)%*%(y-mu)
                fisher <- t(x)%*%(x*sigma)
                choli <- chol(fisher + lam*omega)
                b.new <- b.old + backsolve(choli,
                backsolve(choli, (score - lam*omega%*%b.old), transpose=TRUE))
  
                if(sqrt(sum((b.new-b.old)^2)/sum(b.old^2)) < delta | i>=maxit)
                  {
                    # check the stop criterion
                    conv <- TRUE
                  }
                b.old <- b.new
                i <- i+1
              }  # end while
            coefs[,l] <- b.old   
            fits[,l] <- mu       
          }
        l <- l+1
      }
    rownames(fits) <- NULL
    colnames(fits) <- lambda
    rownames(coefs) <- NULL
    colnames(coefs) <- lambda
    
    return(list(fitted = fits, coefficients = coefs))   
  }



cd <- function(x)
  {
    n <- length(x)

    X <- matrix(1,n,1)

    Z <- matrix(rep(x,max(x)-1),n,max(x)-1)
    Z <- Z - matrix(rep(1:(max(x)-1),n),dim(Z)[1],dim(Z)[2],byrow=T)
    Z[Z<0] <- 0
    Z[Z>1] <- 1

    res <- list(X = X, Z = Z)
    return(res)
  }



ordAOV1 <- function(x, y, type = "RLRT", nsim = 10000, null.sample = NULL ,...){

  x <- as.numeric(x)

  ## check x
  if(min(x)!=1 | length(unique(x)) != max(x))
    stop("x has to contain levels 1,...,max")
    
  if(length(x) != length(y))
    stop("x and y have to be of the same length")

  k <- length(unique(x))
  cdx <- cd(x)
  X <- cdx$X
  Z <- cdx$Z

  # RLRT
  if (type == "RLRT")
  {
  # model under the alternative
  m <- gam(y ~ Z, paraPen=list(Z=list(diag(k-1))), method="REML")

  # null model
  m0 <- gam(y ~ 1, method="REML")

  # log-likelihood
  rlogLik.m <- -summary(m)$sp.criterion
  rlogLik.m0 <- -summary(m0)$sp.criterion
  rlrt.obs <- max(0, 2*(rlogLik.m - rlogLik.m0))

  # null distribution
  if (rlrt.obs != 0) {
      if (length(null.sample)==0)
        RLRTsample <- RLRTSim(X, Z, qr(cdx$X), chol(diag(k-1)), nsim = nsim, ...)
      else
        RLRTsample <- null.sample
        
      p <- mean(rlrt.obs < RLRTsample)
    }
  else
    {
      if (length(null.sample)==0)
        RLRTsample <- NULL
      else
        RLRTsample <- null.sample

      p <- 1
    }

  # return
  RVAL <- list(statistic = c(RLRT = rlrt.obs), p.value = p,
        method = paste("simulated finite sample distribution of RLRT.\n (p-value based on",
        length(RLRTsample), "simulated values)"), sample = RLRTsample)
  }

  # LRT
  else
  {
  # model under the alternative
  m <- gam(y ~ Z, paraPen=list(Z=list(diag(k-1))), method="ML")

  # null model
  m0 <- gam(y ~ 1, method="ML")

  # log-likelihood
  logLik.m <- -summary(m)$sp.criterion
  logLik.m0 <- -summary(m0)$sp.criterion
  lrt.obs <- max(0, 2*(logLik.m - logLik.m0))

  # null distribution
  if (lrt.obs != 0) {
      if (length(null.sample)==0)
        LRTsample <- LRTSim(X, Z, q=0, chol(diag(k-1)), nsim = nsim, ...)
      else
        LRTsample <- null.sample
        
      p <- mean(lrt.obs < LRTsample)
    }
  else
    {
      if (length(null.sample)==0)
        LRTsample <- NULL
      else
        LRTsample <- null.sample

      p <- 1
    }

  # return
  RVAL <- list(statistic = c(LRT = lrt.obs), p.value = p,
        method = paste("simulated finite sample distribution of LRT.\n (p-value based on",
        length(LRTsample), "simulated values)"), sample = LRTsample)
  }

  class(RVAL) <- "htest"
  return(RVAL)
}



ordAOV2 <- function(x, y, type = "RLRT", nsim = 10000, null.sample = NULL, ...){

  n <- length(y)
  p <- ncol(x)
  k <- apply(x, 2, max)

  ## check x
  if(min(x[,1])!=1 | length(unique(x[,1])) != max(x[,1]))
    stop("x has to contain levels 1,...,max")

  if(nrow(x) != length(y))
    stop("nrow(x) and length(y) do not match")

  cdx <- cd(x[,1])
  X <- cdx$X
  ZZ <- vector("list", p)
  ZZ[[1]] <- cdx$Z
  Z <- ZZ[[1]]
  DD <- vector("list", p)
  DD[[1]] <- diag(c(rep(1,k[1]-1),rep(0,sum(k[2:p]-1))))
  for (j in 2:p)
    {

      ## check x
      if(min(x[,j])!=1 | length(unique(x[,j])) != max(x[,j]))
        stop("x has to contain levels 1,...,max")

      ZZ[[j]] <- cd(x[,j])$Z
      Z <- cbind(Z,ZZ[[j]])
      if (j < p)
        DD[[j]] <- diag(c(rep(0,sum(k[1:(j-1)]-1)),rep(1,k[j]-1),rep(0,sum(k[(j+1):p]-1))))
      else
        DD[[j]] <- diag(c(rep(0,sum(k[1:(j-1)]-1)),rep(1,k[j]-1)))
    }
  RRVAL <- vector("list", p)

  # RLRT
  if (type == "RLRT")
  {
  # model under the alternative
  mA <- gam(y ~ Z, paraPen=list(Z=DD), method="REML")

  for (j in 1:p)
  {
  # null model
  Z0 <- matrix(unlist(ZZ[-j]),n,sum(k[-j]-1))
  D0 <- DD[-j]
  if(!is.list(D0))
    D0 <- list(D0)

  if (j == 1)
    out <- 1:(k[1]-1)
  else
    out <- sum(k[1:(j-1)]-1)+(1:(k[j]-1))

  for(j0 in 1:length(D0))
    {
      D0[[j0]] <- D0[[j0]][-out,-out]
    }
  m0 <- gam(y ~ Z0, paraPen=list(Z0=D0), method="REML")

  # log-likelihood
  rlogLik.mA <- -summary(mA)$sp.criterion
  rlogLik.m0 <- -summary(m0)$sp.criterion
  rlrt.obs <- max(0, 2*(rlogLik.mA - rlogLik.m0))

  # model that only contains the variance
  # ...

  # null distribution
  Z1 <- ZZ[[j]]
  if (rlrt.obs != 0) {
      if (length(null.sample) == 0)
        RLRTsample <- RLRTSim(X, Z1, qr(X), chol(diag(k[j]-1)), nsim = nsim, ...)
      else
        {
          if (length(null.sample) == p)
            RLRTsample <- null.sample[[j]]
          else
            stop("wrong number of null.sample elements")
        }

      p <- mean(rlrt.obs < RLRTsample)
    }
  else
    {
      if (length(null.sample) == 0)
        RLRTsample <- NULL
      else
        {
          if (length(null.sample) == p)
            RLRTsample <- null.sample[[j]]
          else
            stop("wrong number of null.sample elements")
        }

      p <- 1
    }

  # return
  RVAL <- list(statistic = c(RLRT = rlrt.obs), p.value = p,
        method = paste("simulated finite sample distribution of RLRT.\n (p-value based on",
        length(RLRTsample), "simulated values)"), sample = RLRTsample)

  class(RVAL) <- "htest"
  RRVAL[[j]] <- RVAL
  }
  }

  # LRT
  else
  {
  # model under the alternative
  mA <- gam(y ~ Z, paraPen=list(Z=DD), method="ML")

  for (j in 1:p)
  {
  # null model
  Z0 <- matrix(unlist(ZZ[-j]),n,sum(k[-j]-1))
  D0 <- DD[-j]
  if(!is.list(D0))
    D0 <- list(D0)

  if (j == 1)
    out <- 1:(k[1]-1)
  else
    out <- sum(k[1:(j-1)]-1)+(1:(k[j]-1))

  for(j0 in 1:length(D0))
    {
      D0[[j0]] <- D0[[j0]][-out,-out]
    }
  m0 <- gam(y ~ Z0, paraPen=list(Z0=D0), method="ML")

  # log-likelihood
  logLik.mA <- -summary(mA)$sp.criterion
  logLik.m0 <- -summary(m0)$sp.criterion
  lrt.obs <- max(0, 2*(logLik.mA - logLik.m0))

  # model that only contains the variance
  # ...

  # null distribution
  Z1 <- ZZ[[j]]
  if (lrt.obs != 0) {
      if (length(null.sample) == 0)
      LRTsample <- LRTSim(X, Z1, q=0, chol(diag(k[j]-1)), nsim = nsim, ...)
      else
        {
          if (length(null.sample) == p)
            LRTsample <- null.sample[[j]]
          else
            stop("wrong number of null.sample elements")
        }

      p <- mean(lrt.obs < LRTsample)
    }
  else
    {
      if (length(null.sample) == 0)
        LRTsample <- NULL
      else
        {
          if (length(null.sample) == p)
            LRTsample <- null.sample[[j]]
          else
            stop("wrong number of null.sample elements")
        }

      p <- 1
    }

  # return
  RVAL <- list(statistic = c(LRT = lrt.obs), p.value = p,
        method = paste("simulated finite sample distribution of LRT.\n (p-value based on",
        length(LRTsample), "simulated values)"), sample = LRTsample)

  class(RVAL) <- "htest"
  RRVAL[[j]] <- RVAL
  }
  }

  names(RRVAL) <- colnames(x)
  return(RRVAL)
}



## Defines a new type of "spline basis" for ordered factors
## with difference penalty:

#' smooth constructor:
smooth.construct.ordinal.smooth.spec <- 
  function(object, data, knots){
    x <- data[[object$term]]
    ## stop if co is not an ordered factor:
    stopifnot(is.ordered(x))
    
    nlvls <- nlevels(x)
    
    #default to 1st order differences (penalizations of deviations from constant):
    if(is.na(object$p.order)) object$p.order <- 1
    
    # construct difference penalty
    Diffmat <- diag(nlvls)
    if(object$p.order > 0){
      for(d in 1:object$p.order) Diffmat <- diff(Diffmat)
    } 
    object$S[[1]] <- crossprod(Diffmat)
    
    # construct design matrix
    object$X <- model.matrix(~ x - 1)
    object$rank <- c(nlvls - object$p.order)
    object$null.space.dim <- object$p.order
    object$knots <- levels(x) 
    object$df <- nlvls
    
    class(object) <- "ordinal.smooth"
    object
  }
#' for predictions:
Predict.matrix.ordinal.smooth <- 
  function(object,data){
    x <- ordered(data[[object$term]], levels = object$knots)
    X <- model.matrix(~ x - 1)
  }

#' for plots:
plot.ordinal.smooth <- 
  function(x,P=NULL,data=NULL,label="",se1.mult=1,se2.mult=2,
           partial.resids=FALSE,rug=TRUE,se=TRUE,scale=-1,n=100,n2=40,n3=3,
           pers=FALSE,theta=30,phi=30,jit=FALSE,xlab=NULL,ylab=NULL,main=NULL,
           ylim=NULL,xlim=NULL,too.far=0.1,shade=FALSE,shade.col="gray80",
           shift=0,trans=I,by.resids=FALSE,scheme=0,...) {
    if (is.null(P)) { ## get plotting info
      xx <- unique(data[x$term])  
      dat <- data.frame(xx)
      names(dat) <- x$term
      X <- PredictMat(x, dat)
      
      if (is.null(xlab)) xlabel <- x$term else xlabel <- xlab
      if (is.null(ylab)) ylabel <- label else ylabel <- ylab
      return(list(X=X,scale=FALSE,se=TRUE, xlab=xlabel,ylab=ylabel,
                  raw = data[x$term], 
                  main="",x=xx,n=n, se.mult=se1.mult))
    } else { ## produce the plot
      n <- length(P$fit)
      if (se) { ## produce CI's
        if (scheme == 1) shade <- TRUE
        ul <- P$fit + P$se ## upper CL
        ll <- P$fit - P$se ## lower CL
        if (scale==0 && is.null(ylim)) { ## get scale 
          ylimit <- c(min(ll), max(ul))
          if (partial.resids) { 
            max.r <- max(P$p.resid,na.rm=TRUE)
            if (max.r> ylimit[2]) ylimit[2] <- max.r
            min.r <-  min(P$p.resid,na.rm=TRUE)
            if (min.r < ylimit[1]) ylimit[1] <- min.r
          }
        }
        if (!is.null(ylim)) ylimit <- ylim
        if (length(ylimit)==0) ylimit <- range(ul,ll)
        ## plot the smooth...
        plot(P$x, trans(P$fit+shift),  
             xlab=P$xlab, ylim=trans(ylimit+shift),
             xlim=P$xlim, ylab=P$ylab, main=P$main, ...)
        if (shade) { 
          rect(xleft=as.numeric(P$x[[1]])-.4, 
               xright=as.numeric(P$x[[1]])+.4, 
               ybottom=trans(ll+shift), 
               ytop=trans(ul+shift), col = shade.col,border = NA)
          segments(x0=as.numeric(P$x[[1]])-.4, 
                   x1=as.numeric(P$x[[1]])+.4, 
                   y0=trans(P$fit+shift), 
                   y1=trans(P$fit+shift), lwd=3)
        } else { ## ordinary plot 
          if (is.null(list(...)[["lty"]])) { 
            segments(x0=as.numeric(P$x[[1]]), 
                     x1=as.numeric(P$x[[1]]), 
                     y0=trans(ul+shift), 
                     y1=trans(ll+shift), lty=2,...)
          } else { 
            segments(x0=as.numeric(P$x[[1]]), 
                     x1=as.numeric(P$x[[1]]), 
                     y0=trans(ul+shift), 
                     y1=trans(ll+shift),...)
            segments(x0=as.numeric(P$x[[1]]), 
                     x1=as.numeric(P$x[[1]]), 
                     y0=trans(ul+shift), 
                     y1=trans(ll+shift),...)
          }
        }
      }  else {
        if (!is.null(ylim)) ylimit <- ylim
        if (is.null(ylimit)) ylimit <- range(P$fit) 
        ## plot the smooth... 
        plot(P$x, trans(P$fit+shift),  
             xlab=P$xlab, ylim=trans(ylimit+shift),
             xlim=P$xlim, ylab=P$ylab, main=P$main, ...)
      }    ## ... smooth plotted
      if (partial.resids&&(by.resids||x$by=="NA")) { ## add any partial residuals
        if (length(P$raw)==length(P$p.resid)) {
          if (is.null(list(...)[["pch"]]))
            points(P$raw,trans(P$p.resid+shift),pch=".",...) else
              points(P$raw,trans(P$p.resid+shift),...) 
        } else {
          warning("Partial residuals not working.")
        }
      } ## partial residuals finished 
      
    }
  } 


crO <- function(k, d=2){
  
  Dd <- cbind(diag(-1,k-1),0) + cbind(0,diag(1,k-1))
  if (d > 1)
  {
    Dj <- Dd
    for (j in 2:d) {
      Dj <- Dj[-1,-1]
      Dd <- Dj%*%Dd
    }
  }
  
  Om <- t(Dd)%*%Dd
  return(Om)
}


penALS <- function(H, p, lambda, qstart, crit, maxit, Ks, constr){
  
  Q <- as.matrix(H) 
  n <- nrow(Q)     
  m <- ncol(Q)      
  
  qs <- list()
  Z <- list()
  ZZO <- list()
  iZZO <- list()
  
  tracing <- c()
  
  for (j in 1:m) 
  { 
    qj <- c(Q[,j],1:Ks[j])                          
    Z[[j]] <- model.matrix(~ factor(qj) - 1)[1:n,] 
    Om <- (Ks[j]-1)*crO(Ks[j])                      
    ZZO[[j]] <- (t(Z[[j]])%*%Z[[j]] + lambda*Om)   
    iZZO[[j]] <- chol2inv(chol(ZZO[[j]]))          
    qs[[j]] <- ((1:Ks[j]) - mean(Q[,j]))/sd(Q[,j]) 
  }
  
  if(length(qstart) > 0){
    Qstart <- mapply("%*%", Z, qstart, SIMPLIFY = TRUE)
    Q <- scale(Qstart)
    qs <- qstart
  }else{
    Q <- scale(Q) 
  }
  
  iter <- 0
  conv <- FALSE
  QQ <- Q
  while(!conv & iter < maxit)
  {
    
    pca <- prcomp(Q, scale = FALSE)  
    X <- pca$x[,1:p, drop = FALSE]        
    A <- pca$rotation[,1:p, drop = FALSE]  
    
    for (j in 1:m){
      
      if (constr[j]){ 
        bvec <- c(n-1,numeric(Ks[j]-1)) 
        dvec <- t(Z[[j]])%*%(X%*%A[j,]) 
        Amat2 <- qs[[j]] %*% t(Z[[j]]) %*% Z[[j]]
        Amat3 <- cbind(0,diag(Ks[j]-1)) - cbind(diag(Ks[j]-1),0)
        Amat <- rbind( Amat2, Amat3)  
        Amat <- t(Amat) 
      }else{
        bvec <- c(n-1)   
        dvec <- t(Z[[j]])%*%(X%*%A[j,]) 
        Amat2 <- qs[[j]] %*% t(Z[[j]]) %*% Z[[j]]
        Amat <- rbind( Amat2)  
        Amat <- t(Amat)
      }

      qj <- solve.QP(Dmat=ZZO[[j]], dvec=dvec, Amat=Amat, bvec=bvec, meq=1)$solution
      
      qj <- as.numeric(qj)
      Zqj <- Z[[j]]%*%qj
      Q[,j] <- Zqj/sd(Zqj)
      qs[[j]] <- qj/sd(Zqj)   
      
    }
    
    tracing <- c(tracing, sum(pca$sdev[1:p]^2) / sum(pca$sdev^2))
    
    # convergence?
    if (sum((QQ - Q)^2)/(n*m) < crit)
      conv <- TRUE
    
    QQ <- Q
    iter <- iter + 1
  }
  
  # final pca
  pca <- prcomp(Q, scale = FALSE)  
  X <- pca$x[,1:p, drop = FALSE]  
  A <- pca$rotation[,1:p, drop = FALSE]  
  
  tracing[length(tracing)] <- sum(pca$sdev[1:p]^2) / sum(pca$sdev^2)
  
  return(list("Q" = Q, "qs" = qs, "iter" = iter, "trace" = tracing))
}



