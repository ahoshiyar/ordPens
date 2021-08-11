\name{ordGene}
\alias{ordGene}
\title{Testing for differentially expressed genes}
\description{This function can be used to test for genes that are differentially
expressed between levels of an ordinal factor, such as dose levels or ordinal
phenotypes.}

\usage{
ordGene(xpr, lvs, type = c("RLRT", "LRT"), nsim = 1e6,
null.sample=NULL, ...)
}
\arguments{
  \item{xpr}{a matrix or data frame of gene expression data with Probe IDs as
  row names.}
  \item{lvs}{a numeric vector containing the factor levels (e.g., dose levels)
  corresponding to the columns of \code{xpr}.}
  \item{type}{the type of test to carry out: likelihood ratio ("LRT") or
  restricted likelihood ratio ("RLRT").}
  \item{nsim}{number of values to simulate from the null distribution.}
  \item{null.sample}{a vector containing values already simulated from the null 
  distribution (overrides \code{nsim})}
  \item{...}{additional arguments to \code{\link{LRTSim}} and
  \code{\link{RLRTSim}}, respectively.}
}
\details{
For each gene in the dataset, \code{\link{ordAOV}} is applied to test for
differences between levels given in \code{lvs}. See \code{\link{ordAOV}} for
further information on the testing procedure. Simulation studies by Gertheiss (2014)
suggest that a restricted likelihood test (RLRT) should rather be used than
a likelihood ratio test (LRT).

In addition to (R)LRT, results of usual one-way ANOVA (not taking the factor's
ordinal scale level into account) and a t-test assuming a linear trend across
factor levels are reported. Note that the t-test does not assume linearity in the
doses (such as 0, 0.5, 2.0, 5.0, ...), if given, but in the levels, i.e., 1, 2, 3, etc.
}

\value{
A matrix containing the raw p-values for each gene (rows) when using (R)LRT,
ANOVA or a t-test (columns).
}

\references{
Crainiceanu, C. and D. Ruppert (2004). \emph{Likelihood ratio tests in linear
mixed models with one variance component}, Journal of the Royal Statistical
Society B, 66, 165-185.

Gertheiss, J. (2014). \emph{ANOVA for factors with ordered levels}, Journal of
Agricultural, Biological and Environmental Statistics, 19, 258-277.

Gertheiss, J. and F. Oehrlein (2011). \emph{Testing relevance and linearity of
ordinal predictors}, Electronic Journal of Statistics, 5, 1935-1959.

Sweeney, E., C. Crainiceanu, and J. Gertheiss (2015). \emph{Testing
differentially expressed genes in dose-response studies and with ordinal
phenotypes}, Statistical Applications in Genetics and Molecular Biology, 15, 213-235.
}

\author{Jan Gertheiss}


\seealso{\code{\link{ordAOV}}}

\examples{
\dontrun{
# use dopamine data from package IsoGene
 
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")
BiocManager::install(c("Biobase","IsoGene")) 

require(IsoGene)
require(Biobase)  
data(dopamine)
xpr <- data.frame(exprs(dopamine))
dose <- unlist(pData(dopamine))
plot(dose,xpr[83,], col=as.factor(dose), lwd=2, ylab="expression") 

# calculate p-values
pvals <- ordGene(xpr = xpr, lvs = dose, nsim=1e6)
# compare distribution of (small) p-values
plot(ecdf(pvals[,1]), xlim=c(0,0.05), ylim=c(0,0.12),
main="dopamine", xlab="p-value", ylab="F(p-value)")
plot(ecdf(pvals[,2]), xlim=c(0,0.05), add=TRUE, col=2)
plot(ecdf(pvals[,3]), xlim=c(0,0.05), add=TRUE, col=3)
legend('topleft', colnames(pvals), col=1:3, lwd=2, lty=1) 
 
}
}


\keyword{anova}
\keyword{gene expression}
