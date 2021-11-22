
---
title: '``ordPens``: An R package for Selection, Smoothing and Principal Components
  Analysis for Ordinal Variables'
authors:
- name: Aisouda Hoshiyar
  orcid: 0000-0002-5702-130X
  affiliation: 1
affiliations:
- name: School of Economics and Social Sciences, Helmut Schmidt University, Hamburg, Germany
  index: 1
date: "\\today"
output: pdf_document
bibliography: references.bib
csl: apa-6th-edition.csl
tags:
- R
- smoothing
- ordinal ANOVA
- Nonlinear Principal Components Analysis
- fusion
- lasso
header-includes:
  \usepackage{bm}
---

# Summary


Ordinal data are a common case in applied statistics. In order to incorporate the ordinal scale level, among other things, regularization techniques are often suggested in the literature [@Tutz:2014; @Tutz:2016]. In particular, penalization approaches for smoothing and selection when dealing with Likert-type data -- which are by no means restricted to Likert scale -- are commonly proposed. 
``ordPens`` is a package in the R programming language [@R] and provides several penalty approaches for ordinal predictors in regression models and ordinal variables for principal component analysis (PCA).    
The main idea behind smoothing in the regression context is to maximize the penalized log-likelihood by introducing a penalty term and a tuning parameter controlling the amount of penalty. Different types of penalization can be considered, depending on whether to achieve smoothing, selection or clustering of variables. Smoothing only can be done by penalizing the sum of squared differences of adjacent coefficients. A modified group lasso based on a difference penalty can be used for selection. Clustering/fusion of categories can be achieved by the fused lasso penalizing absolute differences by using the $L_1$-norm. 

  
# Statement of Need 
 
As suggested by @Tutz:2014 and @Tutz:2016, selection, and/or smoothing/fusing of ordinally scaled independent variables shall be done using a modified group lasso or generalized ridge penalty when dealing with ordinally scaled predictors in regression analysis.
The penalized log-likelihood to be maximized takes the form $\bm{l}_p(\bm{\beta}) = \bm{l}(\bm{\beta}) - \lambda \bm{J}(\bm{\beta})$, with $\bm{\beta}$ corresponding to the vector of regression parameters, $\lambda$ representing the smoothing parameter and $\bm{J}(\cdot)$ being the penalty function. 
\newline
`ordPens` [@ordPens:2021] offers various tools for data analysis of ordinally scaled data. The package attacks the afore mentioned tasks and offers penalized regression for smoothing, selection and fusion. 
Specifically, the function `ordSmooth()` for smoothing only incorporates the generalized ridge penalty
$$
  \bm{J}(\bm{\beta}) = \sum_{s=1}^p \bm{\beta}_s^T \bm{D}_{d,s}^T \bm{D}_{d,s} \bm{\beta}_s,
$$
with $\bm{D}_{d,s}$ being the matrix generating differences of order $d$ and $\bm{\beta}_s^T = (\beta_{s1}, ..., \beta_{sk_{s}})$ being the parameter vector linked to the $s$th (dummy-coded) predictor with categories $1,...,k_s$. The `ordSelect()` function performs smoothing and selection by adopting a modified group lasso penalty based on differences of the form
$$
  \bm{J}(\bm{\beta}) = \sum_{s=1}^p \sqrt{k_s} \sqrt{ \bm{\beta}_s^T \bm{D}_{d,s}^T \bm{D}_{d,s} \bm{\beta}_s }.
$$
Clustering of categories is done by the function `ordFusion()`, which uses a fused lasso penalty based on differences of first order:
$$
  \bm{J}(\bm{\beta}) = \sum_{s=1}^p \sum_{j=2}^{k_s}  |  \beta_{sj} -  \beta_{s,j-1} |.
$$ 
For more information on the original group lasso (for nominal predictors and grouped variables in general), see @Meier:2008 and @Yuan:2006. For details on the fused lasso, see @Tibshirani:2005.
In the case of smoothing only, the package includes auxiliary functions such that ``mgcv::gam()`` [@Wood:2008; @Wood:2017] can be used for fitting generalized linear and additive models with first- and second-order ordinal smoothing penalty as well as built-in smoothing parameter selection. Also, ``mgcv`` tools for further statistical inference can be used, see @Gertheiss:2021 for details.
Furthermore, testing for differences in the means, known as analysis of variance (ANOVA), is provided for ordered factors by the function `ordAOV()` penalizing (squared) differences of adjacent means. Testing for differentially expressed genes, when analyzing microarrays of gene expression data, is incorporated by the function `ordGene()`. 
Technical details can be viewed from @Gertheiss:2014 and @Sweeney:2016, respectively.


If, in contrast, dimension reduction is desired in an unsupervised way, principal components analysis can be applied to ordinal data as well. However, those data are usually either treated as numeric implying linear relationships between the variables at hand, or non-linear PCA is applied where the obtained 
coefficients are sometimes hard to interpret. Note that in IBM SPSS Statistics (Version 25.0), for instance, there is an option available for smoothing quantifications by the use of spline functions, which, however, limits the type of functions that can be fitted when using a small number of knots and a suitable choice may be challenging for the (inexperienced) user. On the other hand, as splines are defined on interval scale whereas ordinal variables can only take some discrete values, the usage of spline functions may be seen as unnecessarily complex for scaling ordinal data. To incorporate the ordinal scale level, the concept of penalization can also be adapted here, as suggested in @Hoshiyar:2021. Penalized non-linear principal components analysis for ordinal variables is incorporated in the function `ordPCA()` using a second-order difference penalty. In addition, the function provides performance evaluation and selection of an optimal penalty parameter using k-fold cross-validation. Also, the option of both non-monotone effects and incorporating constraints enforcing monotonicity is provided. Penalized non-linear PCA therefore serves as an intermediate between the standard methods typically used so far (see above). The new approach offers both better interpretability as well as better performance on validation data.


A topic of future research would be the analysis of dependencies within a (high dimensional) set of ordinal variables by graphical models. A further typical approach when dealing with ordinal data is motivated by assuming a latent continuous variable linked to the ordinal variable via thresholds. The proportional odds model, which is also motivated as a latent variable approach, in combination with the ordinal penalty could be also of interest for future research. Another interesting field is found in @Huang:2021, who analyze (mixed) ordinal dependencies using a latent Gaussian copula model based on rank correlations. 
Assuming a latent continuous variable, however, may not always be desirable by the data analyst. The methods implemented in `ordPens` (up to version 1.0.0) therefore do not underly the latent variable assumption. 

# Availability 
  
The R package ``ordPens`` is publicly available on [CRAN](https://cran.r-project.org/web/packages/ordPens/index.html) and [Github](https://github.com/cran/ordPens), where issues can be opened. ``ordPens`` is licensed under the GPL-2 General Public License.
Documentation and examples are contained in the package manual, which can be found on [CRAN](https://cran.r-project.org/web/packages/ordPens/ordPens.pdf). 

To install ``ordPens``, simply run:
```r
install.packages("ordPens")
```
For penalized regression and ordinal ANOVA see also ``vignette("ordPens", package = "ordPens")``. Penalized non-linear PCA is also documented in detail and can be accessed via ``vignette("ordPCA", package = "ordPens")``. 
  
# ``ordPens`` in action

This example illustrates penalized non-linear PCA on the so-called brief ICF core set on Chronic Widespread Pain (CWP) consisting of 26 ordinally scaled variables. Details on the data can be found in @GerHogObeTut:2011; or by typing `?ICFCoreSetCWP`.
Analysis of the so-called comprehensive core set for CWP, consisting of 67 ICF variables, is found in @Hoshiyar:2021. Note that the `ordPCA` vignette also analyzes the comprehensive core set. Figure 1, generated by the following code, illustrates the estimated coefficients of selected variables for different values of the penalty parameter $\lambda$ along with cross-validation results.  
```r

library(ordPens)

# load ICF data & code adequately
data(ICFCoreSetCWP) 
H <- ICFCoreSetCWP[,1:67] + matrix(c(rep(1,50), rep(5,16), 1), 
                                   nrow(ICFCoreSetCWP), 67, byrow = TRUE)

# select brief core set variables 
brief <- c("b130","b134","b140","b147","b152","b1602","b280","b455",
           "b730","b760","d175","d230","d240","d430","d450","d640",
           "d760","d770","d850","d920","e1101","e310","e355","e410",
           "e420","e570")
           
H <- H[, brief]
xnames <- names(H)

# ordinal penalized PCA 
icf_pca1 <- ordPCA(H, p = 2, lambda = c(5, 0.5, 0.001), qstart = NULL, 
                   crit = 1e-7, maxit = 100, Ks = c(rep(5,20),rep(9,6)), 
                   constr = c(rep(TRUE,20), rep(FALSE,6)), CV = FALSE, 
                   k = 5) 

# 5-fold cross-validation
lambda <- 10^seq(4, -4, by = -0.1)
set.seed(1234) 
cvResult <- ordPCA(H, p = 2, lambda = lambda, Ks = c(rep(5,20),rep(9,6)), 
                   constr = c(rep(TRUE,20), rep(FALSE,6)), 
                   CV = TRUE, k = 5, CVfit = FALSE) 

# plotting results for selected variables
par(mfrow = c(2,3))
for(i in which(xnames=="b280"|xnames=="d450"|xnames=="e1101"|
               xnames=="e410")){
  plot(icf_pca1$qs[[i]][,3], type = "b", xlab = "category", col = 1,
       ylab = "quantification", main = xnames[i], bty = "n", xaxt = "n",
       ylim = range(icf_pca1$qs[[i]])) 
  lines(icf_pca1$qs[[i]][,2], type="b", col=2, lty=2, pch=2, lwd=2)
  lines(icf_pca1$qs[[i]][,1], type="b", col=3, lty=3, pch=3, lwd=2)
  axis(1, at = 1:length(icf_pca1$qs[[i]][,1])) 
}
  
plot(log10(lambda), apply(cvResult$VAFtrain, 2, mean), type = "l",  
     xlab = expression(log[10](lambda)), ylab = "VAF", cex.axis = 1.2,
     main = "training data", cex.lab = 1.2) 

plot(log10(lambda), apply(cvResult$VAFtest, 2, mean), type = "l",
     xlab = expression(log[10](lambda)), ylab = "VAF", cex.axis = 1.2,
     main = "validation data", cex.lab = 1.2)
abline(v = log10(lambda)[which.max(apply(cvResult$VAFtest,2,mean))], 
       lty=2)  
```
 
![Category quantifications/scores for $\lambda \to 0$ (solid black), $\lambda = 0.5$ (dashed red), $\lambda = 5$ (dotted green) (a)–(d); VAF by the first 5 principal components: (e) training data, (f) validation data with optimal $\lambda$ (dashed line).](ordPens_pca.pdf){ width=100% }

 

# Acknowledgements

This work was supported in part by Deutsche Forschungsgemeinschaft (DFG) under Grant GE2353/2-1.

I thank Jan Gertheiss and Fabian Scheipl for their contributions to the software package. 
Jan Gertheiss created initial package versions 0.1-1 up to 0.3-1 and helped discussing the manuscript. Fabian Scheipl implemented the ordinal smoothing penalty for use within ``mgcv``.

# References
