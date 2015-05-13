#' Average Power and True FDR Based on limma/voom RNAseq Analysis Pipeline
#' 
#' For the limma/voom RNAseq analysis pipeline, when we control false discovery
#' rate by using the Benjamini and Hochberg step-up procedure (1995) 
#' and/or Storey and Tibshirani's q-value procedure (Storey et al, 2004),
#' \code{check.power} calculates average power and true FDR for given sample 
#' size, user-specified proportions of non-differentially expressed genes, 
#' number of iterations, FDR level to control, mean counts in control group, 
#' dispersion, and log fold change.
#' 
#' @import qvalue
#'
#' @param nGenes total number of genes to test for differential expression.
#' @param mu mean counts across all tested genes in the control group.
#' @param dispersion mean dispersion across all genes in the control group.
#' @param group treatment group assignments for simulations.
#' @param fdr the target false discovery rate to be controlled.
#' @param pi0 proportion of non-differentially expressed genes. 
#' @param logfc minimum DE log fold change between treatment and control groups.
#' @param sims number of simulations to run when computing power and FDR.
#' 
#' @return \item{pow_bh_ave}{average power when controlling FDR 
#'                           by Benjamini and Hochberg (1995) method.}
#' @return \item{fdr_bh_ave}{true false discovery rate when controlling FDR 
#'                           by Benjamini and Hochberg (1995) method.}
#' @return \item{pow_bh_ave}{average power when controlling FDR 
#'                           by q-value procedure (Storey et al., 2004).}
#' @return \item{fdr_bh_ave}{true false discovery rate when controlling FDR
#'                           by q-value procedure (Storey et al., 2004).}
#' 
#' @author Ran Bi \email{biran@@iastate.edu}, 
#'         Peng Liu \email{pliu@@iastate.edu}, 
#'         Tim Triche \email{tim.triche@@gmail.com}
#' 
#' @references Benjamini, Y. and Hochberg, Y. (1995) 
#'             Controlling the false discovery rate: a practical and 
#'             powerful approach to multiple testing. 
#'             \emph{J. R. Stat. Soc. B}, 57, 289-300.
#' 
#'             Storey, J. D., Taylor, J. E. and Siegmund, D. (2004)
#'             Strong control, conservative point estimation and 
#'             simultaneous rates: a unified approach. 
#'             \emph{J. R. Stat. Soc. B}, 66, 187- 205.
#' 
#' @examples
#' library(limma)
#' library(qvalue)
#' nGenes <- 10000                          ## number of genes to test for DE
#' mu <- 10                                 ## mean read counts in control group
#' dispersion <- 0.1                        ## mean dispersion across all genes
#' group <- rep(c(1, 2), each = 13)         ## group assignments for simulations
#' 
#' check.power(nGenes, mu, dispersion, group)
#'
#' @export
#' 
check.power <- function(nGenes, mu, dispersion, group, 
                        fdr=0.05, pi0=0.8, logfc=1, sims=20) {

  ## empirical "power" & "fdr" function
  powerfdr.fun <- function(fdr, p){ # {{{
    V <- sum( (p < fdr) & (sim$de == FALSE))
    R <- sum( p < fdr )
    S <- R - V
    power <- S / (nGenes * (1 - pi0))
    fdr_true <- V / R
    return(c(power, fdr_true))
  } # }}}
  
  res <- list()
  pow_bh <- fdr_bh <- pow_qvalue <- fdr_qvalue <- rep(0, sims)
  for (j in 1:sims){
    message("Performing simulation ", j, "/", sims, "...")
    sim <- sim.counts(arg=list(nGenes=nGenes, pi0=pi0, group=group), 
                      mu=mu, disp=dispersion, logfc=logfc)
    cts <- sim$counts
    lib.size <- colSums(cts)
    d <- DGEList(cts, lib.size, group=group)
    d <- calcNormFactors(d)
    design <- model.matrix(~ factor(group))
    y <- voom(d, design, plot=FALSE)       # convert read counts to log2-cpm 
                                           # with associated weights
    fit <- lmFit(y, design)
    fit <- eBayes(fit)
    pvalue <- fit$p.value[, 2]             # pvalue
    
    p_bh <- p.adjust(pvalue, method="BH")  # Benjamini & Hochberg
    pow_bh[j] <- powerfdr.fun(fdr, p_bh)[1]
    fdr_bh[j] <- powerfdr.fun(fdr, p_bh)[2]
  
    p_qvalue <- qvalue(pvalue)$qvalues         # Storey and Tibshirani
    pow_qvalue[j] <- powerfdr.fun(fdr, p_qvalue)[1]
    fdr_qvalue[j] <- powerfdr.fun(fdr, p_qvalue)[2]
  }
  message("Simulations completed.")
    
  ## average power & true fdr over sims simulations
  res$pow_bh_ave <- mean(pow_bh)
  res$fdr_bh_ave <- mean(fdr_bh)
  res$pow_qvalue_ave <- mean(pow_qvalue)
  res$fdr_qvalue_ave <- mean(fdr_qvalue)
  return(res)
}
