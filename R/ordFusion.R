ordFusion <- function(x, y, u = NULL, z = NULL, offset = rep(0,length(y)),
                      lambda, model = c("linear", "logit", "poisson"),
                      scalex = TRUE, nonpenx = NULL, frac.arclength = NULL, ...)
  {
  model <- match.arg(model)
  model <- switch(model, linear="linear", logit="logit", poisson="poisson")
  
  ## Check the x matrix
  if(!(is.matrix(x) | is.numeric(x) | is.data.frame(x)))
    stop("x has to be a matrix, numeric vector or data.frame")
  
  if(any(is.na(x)))
    stop("Missing values in x are not allowed")
  
  
  ## Check the response
  if(!is.numeric(y))
    stop("y has to be numeric")
  
  if(model == "logit" & any(!is.element(y[!is.na(y)], c(0,1))))
    stop("y has to be 0/1 coded")
  
  tol <- .Machine$double.eps^0.5 
  if(model == "poisson" & 
     (any(abs(y[!is.na(y)] - round(y[!is.na(y)])) > tol) | any(y[!is.na(y)] < 0)))
    stop("y has to contain nonnegative integers")
  
  
  ## Check the other arguments
  if(length(offset) != length(y))
    stop("length(offset) not equal length(y)")  
  
  if(!is.null(nonpenx))
  {if(max(nonpenx) > ncol(as.matrix(x)))
    stop("max(nonpenx) > ncol(x)")}
  
  ## ordinal predictors  
    x <- as.matrix(x)
    px <- ncol(x)
    kx <- apply(x,2,max)
    xnames <- colnames(x)
    grp <- rep(1:px,kx-1)
    x <- coding(x, constant=FALSE)
    x <- x[!is.na(y),]
    y <- y[!is.na(y)]
    
    
    ## nominal predictors
    if (length(u) > 0)
    {
      if(!(is.matrix(u) | is.numeric(u) | is.data.frame(u)))
        stop("u has to be a matrix, numeric vector or data.frame")
      if(any(is.na(u)))
        stop("Missing values in u are not allowed")
      
      u <- as.matrix(u)
      if(nrow(u) != length(y))
        stop("u and y do not have correct dimensions")
      if(any(!apply(u,2,is.numeric)))
        stop("Entries of u have to be of type 'numeric'")
      if(any(abs(u - round(u)) > tol) | any(u < 1))
        stop("u has to contain positive integers")
      pu <- ncol(u)
      ku <- apply(u,2,max)
      unames <- colnames(u)
      grp <- c(grp,rep(max(grp)+(1:pu),ku-1))
      u <- coding(u, constant=FALSE, splitcod=FALSE)
      u <- u[!is.na(y),]
      nonpenu <- 1:pu
    }
    else
    {
      pu <- NULL
      nonpenu <- NULL
      ku <- NULL
    }
    
    ## metric predictors
    if (length(z) > 0)
    {
      if(!(is.matrix(z) | is.numeric(z) | is.data.frame(z)))
        stop("z has to be a matrix, numeric vector or data.frame")
      if(any(is.na(z)))
        stop("Missing values in z are not allowed")
      
      z <- as.matrix(z)
      if(nrow(z) != length(y))
        stop("z and y do not have correct dimensions")
      if(any(!apply(z,2,is.numeric)))
        stop("Entries of z have to be of type 'numeric'")
      pz <- ncol(z)
      znames <- colnames(z)
      grp <- c(grp,max(grp)+(1:pz))
      z <- z[!is.na(y),]
      nonpenz <- 1:pz
    }
    else
    {
      pz <- NULL
      nonpenz <- NULL
    }
    
    xuz <- cbind(x,u,z)
    offset <- offset[!is.na(y)]
    y <- y[!is.na(y)]
    
    ## fitting
    nonpen <- c(nonpenx, px+nonpenu, px+ifelse(length(pu)>0,pu,0)+nonpenz)
    if (length(nonpen) > 0)
    {
      nonpen2 <- 1:ncol(xuz)
      nonpen2[!is.element(grp,nonpen)] <- -1
      nonpen2 <- nonpen2[nonpen2 > 0]
    }
    else
    {
      nonpen2 <- NULL
    }
    
    if (model == "linear")
      {
        if (is.null(frac.arclength))
          frac.arclength <- 1
        
        flmodel <- glmpath(xuz, y, nopenalty.subset = nonpen2,
                           family = gaussian, standardize = scalex,
                           frac.arclength = frac.arclength, ...)
      }
    else if (model == "logit")
      {
        if (is.null(frac.arclength))
          frac.arclength <- 0.1
        
        flmodel <- glmpath(xuz, y, nopenalty.subset = nonpen2,
                           family = binomial, standardize = scalex,
                           frac.arclength = frac.arclength, ...)
      }
    else
      {
        if (is.null(frac.arclength))
          frac.arclength <- 0.1
      
        flmodel <- glmpath(xuz, y, nopenalty.subset = nonpen2,
                           family = poisson, standardize = scalex,
                           frac.arclength = frac.arclength, ...)
      }
    
    lmbd <- lmbs <- sort(lambda, decreasing = TRUE)
    lmbd <- c(max(flmodel$lambda),lmbd[lmbd < max(flmodel$lambda)])
    
    flx <- predict(flmodel, s = lmbd, type = "coefficients",
    mode = "lambda")
    constant <- flx[,1]
    xuzc <- t(flx[,-1])
    xgrp <- grp[1:ncol(x)]
    xc <- xuzc[1:ncol(x),]
    
    for (j in 1:max(xgrp))
      {
        xc[xgrp==j,] <- apply(xc[xgrp==j, ,drop=FALSE],2,cumsum)
      }
    if (length(xnames)==0)
      xnames <- paste("x",1:px,sep="")

    xnames <- rep(xnames,kx)
    xnames <- paste(xnames,":",sequence(kx),sep="")
    
    xcoefs <- matrix(0,length(xnames),ncol(xc))
    xcoefs[sequence(kx)>1,] <- xc

    if (length(u) > 0)
    {
      if (ncol(u) > 1)
        uc <- cbind(xuzc[ncol(x)+(1:ncol(u)),])
      else
        uc <- rbind(xuzc[ncol(x)+1,])
      
      if (length(unames)==0)
        unames <- paste("u",1:pu,sep="")
      
      unames <- rep(unames,ku)
      unames <- paste(unames,":",sequence(ku),sep="")
      
      ucoefs <- matrix(0,length(unames),ncol(uc))
      ucoefs[sequence(ku)>1,] <- uc
    }
    else
    {
      ucoefs <- NULL
      unames <- NULL
    }
    
    if (length(z) > 0)
    {
      if (ncol(z) > 1)
        zcoefs <- cbind(xuzc[length(grp)+2-(ncol(z):1),])
      else
        zcoefs <- rbind(xuzc[length(grp)+1,])
      
      if (length(znames)==0)
        znames <- paste("z",1:pz,sep="")
    }
    else
    {
      zcoefs <- NULL
      znames <- NULL
    }
    
    coefs <- rbind(constant,xcoefs,ucoefs,zcoefs)
    rownames(coefs) <- c("intercept",xnames,unames,znames)
    colnames(coefs) <- lmbd

    fits <- predict(flmodel, newx = xuz, s = lmbd, type = "response",
    mode = "lambda")
    colnames(fits) <- lmbd
    rownames(fits) <- NULL
    
    
    ladd <- sum(lmbs > max(flmodel$lambda))
    if (ladd > 0)
    {
      cadd <- matrix(coefs[,1],nrow(coefs),ladd, byrow = FALSE)
      fadd <- matrix(fits[,1],nrow(fits),ladd, byrow = FALSE)
      coefs <- cbind(cadd,coefs[,-1])
      fits <- cbind(fadd,fits[,-1])
    }
    else
    {
      coefs <- coefs[,-1]
      fits <- fits[,-1]
    }
    
    colnames(coefs) <- colnames(fits) <- lmbs
    rownames(fits) <- NULL

    out <- list(fitted = fits,
                coefficients = coefs,
                model = model,
                lambda = lmbs,
                xlevels = kx,
                ulevels = ku,
                zcovars = length(znames))
    structure(out, class = "ordPen")
  }