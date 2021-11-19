

library(testthat)
library(ordPens)

context("check ordinal penalized PCA (ordPCA) function")


data(ICFCoreSetCWP) 
H <- ICFCoreSetCWP[, 1:67] + matrix(c(rep(1, 50), rep(5, 16), 1),
                                    nrow(ICFCoreSetCWP), 67,
                                    byrow = TRUE)
y <- ICFCoreSetCWP$phcs 


## checks on warning and error messages
test_that("lambda in ordPCA gets sorted properly",
          {
           
            expect_warning( 
              ordPCA(H, p = 2, lambda = sort(c(5, 0.5, 0.001)), 
                     Ks = c(rep(5,50), rep(9,16), 5), 
                     constr = c(rep(TRUE,50), rep(FALSE,16), TRUE) 
                      ),  
              "lambda values should be sorted in decreasing order" 
            )       
                      
})

test_that("CV and fit works well for large lambda vec",
          { 
            expect_warning( 
              myPCA <- ordPCA(H, p = 2, lambda = seq(6,1,by=-1), 
                              Ks = c(rep(5,50), rep(9,16), 5), 
                              constr = c(rep(TRUE,50), rep(FALSE,16), TRUE),
                              CV = TRUE, CVfit = TRUE),
              "Lists of undesireably high dimensions will be produced"
            )     
})
yes 
test_that("CV and fit error checking",
          { 
            expect_error( 
              myPCA <- ordPCA(H, p = 2, lambda = seq(6,1,by=-1), 
                              Ks = c(rep(5,50), rep(9,16), 5), 
                              constr = c(rep(TRUE,50), rep(FALSE,16), TRUE),
                              CV = TRUE, CVfit = TRUE),
              "calculation canceled"
            )     
})
cancel 
test_that("CV only for large lambda vec works properly",
          { 
            expect_message( 
              myPCA <- ordPCA(H, p = 2, lambda = seq(6,1,by=-1), 
                              Ks = c(rep(5,50), rep(9,16), 5), 
                              constr = c(rep(TRUE,50), rep(FALSE,16), TRUE),
                              CV = TRUE, CVfit = TRUE),
              "VAF computed only"
            )     
})
no
 
 
