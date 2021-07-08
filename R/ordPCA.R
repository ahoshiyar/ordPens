ordPCA <- function(H, p, lambda = c(1), maxit = 100, crit = 1e-7,  
                    Qstart = NULL, Ks = apply(H,2,max), constr = rep(FALSE, ncol(H)),
                    CV = FALSE, k = 10, CVsupress = TRUE){ 
  
   
     
  
  if(length(lambda) > 5 & CVsupress == FALSE){
    ask <- askYesNo("Lists of undesireably high dimensions will be produced. Do you want to proceed?", 
                    prompts = getOption("askYesNo", gettext(c("Yes", "No", "Cancel")))
    ) 
    if(!is.na(ask) & ask == TRUE){
      CVsupress <- FALSE
      warning("Lists of undesireably high dimensions will be produced")
    }else if(is.na(ask)){
      stop("calculation canceled", call. = FALSE)
    }else{
      CVsupress <- TRUE
      message("VAF computed only")
    } 
  } 
  
  if(!(CV)){
    
    if(length(lambda) == 1){ 
      
      res <- penALS(H = H, p = p, lambda = lambda, Qstart = Qstart, Ks = Ks,
                    constr = constr, maxit = maxit, crit = crit)
      
      # final pca
      pca <- prcomp(res$Q, scale = FALSE)
      X <- pca$x[,1:p, drop = FALSE]  
      A <- pca$rotation[,1:p, drop = FALSE]  
      
      out <- list("qs" = res$qs, "Q" = res$Q, "X" = X, "A" = A, "iter" = res$iter, "pca" = pca, 
                  "VAFtrain" = NULL, "VAFtest" = NULL) 
      
    }else if(length(lambda) > 1){
      
      pcai <- list()
      Ai <- list()
      Xi <- list()
      Qi <- list() 
      iteri <- c()
      qsji <- list()
      for(j in 1:ncol(as.matrix(H))){
        qsji[[j]] <- matrix(NA, ncol = length(lambda), nrow = Ks[j])
      }
      
      for(i in 1:length(lambda)){
        
        res <- penALS(H = H, p = p, lambda = lambda[i], Qstart = Qstart, Ks = Ks,
                      constr = constr, maxit = maxit, crit = crit)
        
        # final pca
        pcai[[i]] <- pcafit <- prcomp(res$Q, scale = FALSE)
        Xi[[i]] <- pcafit$x[,1:p, drop = FALSE]  
        Ai[[i]] <- pcafit$rotation[,1:p, drop = FALSE] 
        Qi[[i]] <- Qstart <- res$Q
        iteri[i] <- res$iter
        for(j in 1:ncol(as.matrix(H))){
          qsji[[j]][,i] <- res$qs[[j]]
        }
      }  
      out <- list("qs" = qsji, "Q" = Qi, "X" = Xi, "A" = Ai, "iter" = iteri, "pca" = pcai,
                  "VAFtrain" = NULL, "VAFtest" = NULL)
    } 
    
  }else if(CV){  
    
    nh <- nrow(H) 
    mh <- ncol(H) 
    os <- sample(nh,nh)    
    
    if(length(lambda) == 1){ 
      
      TVAFtrain <- matrix(NA,k,length(lambda))
      TVAFtest <- matrix(NA,k,length(lambda))
      
      nk <- floor(nh/k)
      for (wk in 1:k)
      {
        if (wk < k)
          inds <- os[(wk-1)*nk+(1:nk)] # groups of individuals without kth fold
        else
          inds <- os[-(1:((wk-1)*nk))] # group of kth fold
        
        # prepare test obs
        Hk <- H[inds,] # data submatrix (test data)
        Zk <- list()
        for (j in 1:mh)
        {
          hkj <- c(Hk[,j],1:Ks[j])                                     # ensure all categories to be represented
          Zk[[j]] <- model.matrix( ~ factor(hkj) - 1)[1:length(inds),] # quantification submatrix for test data 
        }
        Qk <- matrix(NA,length(inds),mh) # data submatrix for test data 
        
        # Qstart <- NULL # for the first fold
        
        # ordPCA: 
        # training obs 
        res <- penALS(H = H[-inds,], p = p, lambda = lambda, Qstart = Qstart, Ks = Ks,
                      constr = constr, maxit = maxit, crit = crit)
        
        # final pca
        pca <- prcomp(res$Q, scale = FALSE)
        X <- pca$x[,1:p, drop = FALSE]  
        A <- pca$rotation[,1:p, drop = FALSE]  
        
        Qstart <- pca$Q # update Qstart for remaining folds 
        
        # VAF on training set
        TVAFtrain[wk,] <- sum((pca$sdev[1:p])^2)/sum(pca$sdev^2) 
        
        # scaled test obs
        for (j in 1:ncol(as.matrix(H[-inds,])))
          Qk[,j] <- Zk[[j]]%*%res$qs[[j]]  # Zk from kth fold, quantifications from training data
        # result to validate model chosen on training data on testing obs
        
        # pca on test obs
        pcak <- princomp(Qk, cor = TRUE) 
        # VAF on test obs
        TVAFtest[wk,] <- sum((pcak$sdev[1:p])^2)/sum(pcak$sdev^2)
        
      }
      
      VAFtrain <- apply(TVAFtrain,2,mean)
      VAFtest <- apply(TVAFtest,2,mean)
      
      out <- list("qs" = res$qs, "Q" = res$Q, "X" = X, "A" = A, "iter" = res$iter, 
                  "pca" = pca, "VAFtrain" = VAFtrain, "VAFtest" = VAFtest) 
      
    }else if(length(lambda) > 1){
      
      if(CVsupress){
        
        TVAFtrain <- matrix(NA,k,length(lambda))
        TVAFtest <- matrix(NA,k,length(lambda))
        nk <- floor(nh/k)
        for (wk in 1:k)
        {
          if (wk < k)
            inds <- os[(wk-1)*nk+(1:nk)]  
          else
            inds <- os[-(1:((wk-1)*nk))]  
          
          # test obs
          Hk <- H[inds,]  
          Zk <- list()
          for (j in 1:mh)
          {
            hkj <- c(Hk[,j],1:Ks[j])  
            Zk[[j]] <- model.matrix( ~ factor(hkj) - 1)[1:length(inds),]
          }
          Qk <- matrix(NA,length(inds),mh)  
          
          # Qstart <- NULL
          
          # ordPCA: 
          for (i in 1:length(lambda))
          {
            # ordinal pca
            # training obs
            res <- penALS(H = H[-inds,], p = p, lambda = lambda[i], Qstart = Qstart, Ks = Ks,
                          constr = constr, maxit = maxit, crit = crit)
            
            # final pca
            pcafit <- prcomp(res$Q, scale = FALSE)
            Qstart <- res$Q
            
            # VAF on training obs
            TVAFtrain[wk,i] <- sum((pcafit$sdev[1:p])^2)/sum(pcafit$sdev^2)  
            
            # scaled test obs
            for (j in 1:ncol(as.matrix(H[-inds,])))
              Qk[,j] <- Zk[[j]]%*%res$qs[[j]]  
            
            # pca on test obs
            pcak <- princomp(Qk, cor = TRUE) 
            # VAF on test obs
            TVAFtest[wk,i] <- sum((pcak$sdev[1:p])^2)/sum(pcak$sdev^2)
            
          } 
        }
        
        VAFtrain <- apply(TVAFtrain,2,mean)
        VAFtest <- apply(TVAFtest,2,mean)

        out <- list("qs" = NULL, "Q" = NULL, "X" = NULL, "A" = NULL, "iter" = NULL, 
                    "pca" = NULL, "VAFtrain" = VAFtrain, "VAFtest" = VAFtest) 
        
      }else if(!CVsupress){
        
        TVAFtrain <- matrix(NA,k,length(lambda))
        TVAFtest <- matrix(NA,k,length(lambda))
        nk <- floor(nh/k)
        for (wk in 1:k)
        {
          if (wk < k)
            inds <- os[(wk-1)*nk+(1:nk)]  
          else
            inds <- os[-(1:((wk-1)*nk))]  
          
          # test obs
          Hk <- H[inds,]  
          Zk <- list()
          for (j in 1:mh)
          {
            hkj <- c(Hk[,j],1:Ks[j])  
            Zk[[j]] <- model.matrix( ~ factor(hkj) - 1)[1:length(inds),]
          }
          Qk <- matrix(NA,length(inds),mh)  
          
          # Qstart <- NULL
          
          pcai <- list()
          Ai <- list()
          Xi <- list()
          Qi <- list() 
          iteri <- c()
          qsji <- list()
          for(j in 1:ncol(as.matrix(H))){
            qsji[[j]] <- matrix(NA, ncol = length(lambda), nrow = Ks[j])
          }
          
          # ordPCA: 
          for (i in 1:length(lambda))
          {
            # ordinal pca
            # training obs
            res <- penALS(H = H[-inds,], p=p, lambda = lambda[i], Qstart = Qstart, Ks = Ks,
                          constr = constr, maxit = maxit, crit = crit)
            
            # final pca
            pcai[[i]] <- pcafit <- prcomp(res$Q, scale = FALSE)
            Xi[[i]] <- pcafit$x[,1:p, drop = FALSE]  
            Ai[[i]] <- pcafit$rotation[,1:p, drop = FALSE] 
            Qi[[i]] <- Qstart <- res$Q
            iteri[i] <- res$iter
            for(j in 1:ncol(as.matrix(H))){ 
              qsji[[j]][,i] <- res$qs[[j]]
            }
            
            # VAF on training obs
            TVAFtrain[wk,i] <- sum((pcafit$sdev[1:p])^2)/sum(pcafit$sdev^2)  
            
            # scaled test obs
            for (j in 1:ncol(as.matrix(H)))
              Qk[,j] <- Zk[[j]]%*%res$qs[[j]]  
            
            # pca on test obs
            pcak <- princomp(Qk, cor = TRUE) 
            # VAF on test obs
            TVAFtest[wk,i] <- sum((pcak$sdev[1:p])^2)/sum(pcak$sdev^2)
            
          } 
        }
        
        VAFtrain <- apply(TVAFtrain,2,mean)
        VAFtest <- apply(TVAFtest,2,mean)
        
        out <- list("qs" = qsji, "Q" = Qi, "X" = Xi, "A" = Ai, "iter" = iteri, "pca" = pcai,
                    "VAFtrain" = VAFtrain, "VAFtest" = VAFtest)
      }
    }  
  }   
  return(out)
}
