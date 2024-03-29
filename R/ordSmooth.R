ordSmooth <- function(x, y, u = NULL, z = NULL, offset = rep(0,length(y)), 
  lambda, model = c("linear", "logit", "poisson"), restriction = c("refcat", "effect"), penscale = identity,  
  scalex = TRUE, nonpenx = NULL, eps = 1e-3, delta = 1e-6, maxit = 25, ...)
  {
    model <- match.arg(model)
    model <- switch(model, linear="linear", logit="logit", poisson="poisson")
    restriction <- match.arg(restriction) 
    
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
      
    if(is.unsorted(rev(lambda)))
      warning("lambda values should be sorted in decreasing order")

    lmbd <- sort(lambda, decreasing = TRUE)
    
    ## ordinal predictors
    x <- as.matrix(x)
    if(nrow(x) != length(y))
      stop("x and y do not have correct dimensions")
    if(any(!apply(x,2,is.numeric)))
      stop("Entries of x have to be of type 'numeric'")
    if(any(abs(x - round(x)) > tol) | any(x < 1))
      stop("x has to contain positive integers")
    px <- ncol(x)
    kx <- apply(x,2,max)
    xnames <- colnames(x)
    grp <- rep(1:px,kx-1)
    x <- coding(x, constant=TRUE)
    x <- x[!is.na(y),]
    p0 <- ncol(x)-1
    
        if (scalex)
          {
            stdx <- apply(cbind(x),2,sd)[-1]
            stdx[stdx==0] <- 1
            stdx <- stdx + eps
          }
        else
          {
            stdx <- rep(1,ncol(cbind(x))-1)
          }
        x <- scale(x,center=FALSE,scale=c(1,stdx))

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
    omega <- rep(1,length(grp))
    omega <- omega*penscale(rep(table(grp),table(grp)))
    if (length(nonpen) > 0)
      {
        omega[is.element(grp,nonpen)] <- 0
      }
    omega <- c(0,omega)
    omega <- diag(omega)
    
    ridgemodel <- genRidge(x=xuz,y=y,offset=offset,omega=omega,lambda=lmbd,
    model=model,delta=delta,maxit=maxit)
      
    ## fitted coefficients
    constant <- as.numeric(ridgemodel$coef[1,])
    if (ncol(x) > 2)
      xc <- cbind(ridgemodel$coef[2:ncol(x),]/stdx)
    else
      xc <- rbind(ridgemodel$coef[2,]/stdx)
          
    xgrp <- grp[1:(ncol(x)-1)]

    for (j in 1:max(xgrp))
      {
        if (sum(xgrp==j) > 1)
          xc[xgrp==j,] <- apply(cbind(xc[xgrp==j,]),2,cumsum)
        else
          xc[xgrp==j,] <- apply(rbind(xc[xgrp==j,]),2,cumsum)
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
          uc <- cbind(ridgemodel$coef[ncol(x)+(1:ncol(u)),])
        else
          uc <- rbind(ridgemodel$coef[ncol(x)+1,])
        
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
          zcoefs <- cbind(ridgemodel$coef[length(grp)+2-(ncol(z):1),])
        else
          zcoefs <- rbind(ridgemodel$coef[length(grp)+1,])
          
        if (length(znames)==0)
          znames <- paste("z",1:pz,sep="")
      }
    else
      {
        zcoefs <- NULL
        znames <- NULL
      }

    
    ## Effect coding / sum to zero restriction
    if(restriction == "effect")
    {
      modelbjk <-  xc[1:p0, ,drop=F]  
      transcoef <- matrix(NA, p0, ncol(xc)) 
      transb1 <- matrix(NA, px, ncol(xc)) 
      for(j in 1:px){
        means <- drop(apply(modelbjk[grp==j,, drop=F], 2, function(x) mean(c(0,x))) )
        for(i in 1:ncol(xc)){
          transcoef[grp==j,i] <- beffcjk <- modelbjk[grp==j,i] - means[i]
          transb1[j,i] <- - sum(beffcjk)
        }   
      }
      xcoefs[sequence(kx)>1,] <- transcoef 
      xcoefs[sequence(kx)==1,] <- transb1 
    } 
    
    
    coefs <- rbind(constant,xcoefs,ucoefs,zcoefs)
    rownames(coefs) <- c("intercept",xnames,unames,znames)
    colnames(coefs) <- lmbd
    fits <- ridgemodel$fitted
    colnames(fits) <- lmbd
    rownames(fits) <- NULL
    
    ## output
    out <- list(fitted = fits,
                coefficients = coefs,
                model = model,
                restriction = restriction, 
                lambda = lmbd,
                xlevels = kx,
                ulevels = ku,
                zcovars = length(znames))
    structure(out, class = "ordPen")    
  }
