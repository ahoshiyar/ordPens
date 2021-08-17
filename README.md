
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `ordPens`: Selection and/or Smoothing and Principal Components Analysis for Ordinal Variables

<!-- badges: start -->
<!-- badges: end -->

We provide selection, and/or smoothing/fusion of ordinally scaled
independent variables using a group lasso or generalized ridge penalty.
In addition, nonlinear principal components analysis for ordinal
variables is offered, using a second-order difference penalty.

Also, ANOVA with ordered factors is provided by the function `ordAOV`;
testing for differentially expressed genes can be done using `ordGene`.
For details cf. Gertheiss (2014) and Sweeney et al. (2015),
respectively.

For smoothing, selection and fusion, details may be found in Tutz and
Gertheiss (2014, 2016). All functions are documented in detail in
`vignette("ordPens", package = "ordPens")`. For smoothing only, the
package also builds a bridge to `mgcv::gam()`, see Gertheiss et
al. (2021) for further information.

For the function implementing nonlinear principal components analysis,
`ordPCA`, details can be found in Hoshiyar et al. (2021) and
`vignette("ordPCA", package = "ordPens")`.

Version 1.0.0 is a major release with new functions:

-   `ordPCA` applies nonlinear principal components analysis for ordinal
    variables. Also, performance evaluation and selection of an optimal
    penalty parameter provided.  
-   `ordFusion` fits dummy coefficients of ordinally scaled independent
    variables with a fused lasso penalty for fusion and selection.
-   A new type of spline basis for ordered factors
    `s(..., bs = "ordinal")`is provided, such that smooth terms in the
    `mgcv::gam()` formula can be used as an alternative and extension to
    `ordSmooth()`. Additionally, generic functions for prediction and
    plotting are provided.

## References

-   Gertheiss, J. (2014). ANOVA for factors with ordered levels.
    *Journal of Agricultural, Biological and Environmental Statistics
    19*, 258-277.

-   Gertheiss, J., F. Scheipl, T. Lauer, and H. Ehrhardt (2021).
    Statistical inference for ordinal predictors in generalized linear
    and additive models with application to bronchopulmonary dysplasia.
    Preprint, available from <https://arxiv.org/abs/2102.01946>.

-   Hoshiyar, A., H.A.L. Kiers, and J. Gertheiss (2021). Penalized
    non-linear principal components analysis for ordinal variables with
    an application to international classification of functioning core
    sets, Preprint.

-   Sweeney, E., C. Crainiceanu, and J. Gertheiss (2015). Testing
    differentially expressed genes in dose-response studies and with
    ordinal phenotypes. *Statistical Applications in Genetics and
    Molecular Biology 15*, 213-235.

-   Tutz, G. and J. Gertheiss (2014). Rating scales as predictors – the
    old question of scale level and some answers. *Psychometrica 79*,
    357-376.

-   Tutz, G. and J. Gertheiss (2016). Regularized regression for
    categorical data. *Statistical Modelling 16*, 161-200.
