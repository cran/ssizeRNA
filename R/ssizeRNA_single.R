#' Sample Size Calculations for Two-Sample RNA-seq Experiments with Single Set of Parameters
#' 
#' This function calculates appropriate sample sizes for two-sample RNA-seq experiments for a desired power 
#' in which mean and dispersion parameters are identical for all genes. 
#' Sample size calculations are performed at controlled false discovery rates, user-specified proportions 
#' of non-differentially expressed genes, mean counts in control group, dispersion, and log fold change. 
#' A plot of power versus sample size is generated.
#' 
#' @import edgeR limma
#' @param nG total number of genes.
#' @param pi0 a vector (or scalar) of proportions of non-differentially expressed genes.
#' @param mu a scalar of mean counts in control group from which to simulate.
#' @param disp a scalar of dispersion parameter from which to simulate.
#' @param logfc log fold change between treatment group and control group.
#' @param up proportion of up-regulated genes among all differentially expressed genes, the default value is \code{0.5}.
#' @param replace sample with or without replacement from given parameters. See Details for more information.
#' @param m pseudo sample size for generated data.
#' @param fdr the false discovery rate to be controlled.
#' @param power the desired power to be achieved.
#' @param maxN the maximum sample size used for power calculations.
#' @param side options are "two-sided", "upper", or "lower".
#' @param cex.title controls size of chart titles.
#' @param cex.legend controls size of chart legend.
#' 
#' @details If a vector is input for \code{pi0}, sample size calculations are performed for each proportion.
#' 
#' If the total number of genes is larger than length of \code{mu} or \code{disp}, 
#' \code{replace} always equals \code{TRUE}.
#' 
#' @return \item{ssize}{sample sizes (for each treatment) at which desired power is first reached.}
#' @return \item{power}{power calculations with corresponding sample sizes.}
#' @return \item{crit.vals}{critical value calculations with corresponding sample sizes.}
#' 
#' @author Ran Bi \email{biran@@iastate.edu}, Peng Liu \email{pliu@@iastate.edu}
#' 
#' @references Liu, P. and Hwang, J. T. G. (2007) Quick calculation for sample size while controlling 
#' false discovery rate with application to microarray analysis. \emph{Bioinformatics} 23(6): 739-746. 
#' 
#' Orr, M. and Liu, P. (2009) Sample size estimation while controlling false discovery rate for microarray 
#' experiments using ssize.fdr package. \emph{The R Journal}, 1, 1, May 2009, 47-53. 
#' 
#' Law, C. W., Chen, Y., Shi, W., Smyth, G. K. (2014). Voom: precision weights unlock linear model analysis tools for 
#' RNA-seq read counts. \emph{Genome Biology} 15, R29.
#' 
#' @seealso \code{\link{ssizeRNA_vary}}
#' 
#' @examples
#' p0 <- 0.8                                ## proportion of non-differentially expressed genes
#' mu <- 10                                 ## mean counts in control group for all genes
#' disp <- 0.1                              ## dispersion for all genes
#' logfc <- log(2)                          ## log fold change for up-regulated genes
#' 
#' size <- ssizeRNA_single(pi0 = p0, mu = mu, disp = disp, logfc = logfc, m = 30, maxN = 20)
#' size$ssize                               ## first sample size to reach desired power
#' size$power                               ## calculated power for each sample size
#' size$crit.vals                           ## calculated critical value for each sample size
#' 
#' @export
#' 
ssizeRNA_single <- function(nG = 10000, pi0 = 0.95, mu, disp, logfc, up = 0.5, replace = TRUE,
                            m = 200, fdr = 0.05, power = 0.8, maxN = 35, 
                            side = "two-sided", cex.title = 1.15, cex.legend = 1){ 

  arg = list(
    nG = nG,  # total number of genes
    pi0 = pi0,  # proportion of non-differentially expressed genes
    group = rep(c(1, 2), each = m)  # treatment groups
  )
  
  sim <- sim.counts(arg, mu, disp, logfc, up, replace)

  d_cpm <- DGEList(sim$counts)
  d_cpm <- calcNormFactors(d_cpm)
  design <- model.matrix(~factor(arg$group))
  y <- voom(d_cpm, design, plot = F)  # convert read counts to log2-cpm with associated weights
  fit <- lmFit(y, design)
  fit <- eBayes(fit)
  fit$logcpm <- y$E  # normalized log-cpm value
  fit$weights <- y$weights  # precision weights for each observation 
  fit$www <- sqrt(2/m * apply(y$weights[, 1:m], 1, sum) * apply(y$weights[, (m+1):(2*m)], 1, sum)
                  / apply(y$weights, 1, sum))
  fit$Delta <- fit$coef[, 2] * fit$www  # effect size defining as weighted mean difference of log-cpm values, Delta_g
  
  dm <- mean((fit$Delta * sim$de)[sim$de!=0])
  ds <- sd((fit$Delta * sim$de)[sim$de!=0])
  sig <- density(fit$sigma)$x[which.max(density(fit$sigma)$y)]
  
  ret <- ssize.twoSampVaryDelta(deltaMean = dm, deltaSE = ds, sigma = sig, fdr = fdr, power = power, pi0 = pi0, 
                                maxN = maxN, side = side, cex.title = cex.title, cex.legend = cex.legend)
  return(ret)
}