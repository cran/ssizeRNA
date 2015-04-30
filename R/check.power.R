#' Average Power and True FDR Based on Voom and Limma Pipeline
#' 
#' For the voom and limma pipeline method, when we control false discovery rate by using both Benjamini 
#' and Hochberg (1995) method and q-value procedure (Storey et al., 2004), \code{check.power} calculates
#' average power and true FDR for given sample size, user-specified proportions of non-differentially 
#' expressed genes, number of iterations, FDR level to control,
#' mean counts in control group, dispersion, and log fold change.
#' 
#' @import qvalue
#' @param arg a list of global parameters to pass into the function, such as total number of genes, 
#' proportion of non-differentially expressed genes and treatment groups. See Details for more information.
#' @param k number of iterations.
#' @param fdr the false discovery rate to be controlled.
#' @param mu a scalar of mean counts in control group from which to simulate.
#' @param disp a scalar of dispersion parameter from which to simulate.
#' @param logfc log fold change between treatment group and control group.
#' 
#' @details \code{arg = list(nG, pi0, group)} where \code{nG} is the total number of genes, \code{pi0} 
#' is the proportion of non-differentially expressed genes, and \code{group} is the treatment groups.
#' 
#' @return \item{pow_bh_ave}{average power when controlling FDR by Benjamini 
#' and Hochberg (1995) method.}
#' @return \item{fdr_bh_ave}{true false discovery rate when controlling FDR by Benjamini 
#' and Hochberg (1995) method.}
#' @return \item{pow_bh_ave}{average power when controlling FDR by q-value procedure 
#' (Storey et al., 2004).}
#' @return \item{fdr_bh_ave}{true false discovery rate when controlling FDR by q-value procedure 
#' (Storey et al., 2004).}
#' 
#' @author Ran Bi \email{biran@@iastate.edu}, Peng Liu \email{pliu@@iastate.edu}
#' 
#' @references Benjamini, Y. and Hochberg, Y. (1995) Controlling the false discovery rate: a practical and 
#' powerful approach to multiple testing. \emph{J. R. Stat. Soc. B}, 57, 289-300.
#' 
#' Storey, J. D., Taylor, J. E. and Siegmund, D. (2004) Strong control, conservative point estimation and 
#' simultaneous rates: a unified approach. \emph{J. R. Stat. Soc. B}, 66, 187- 205.
#' 
#' @examples
#' library(edgeR)
#' library(qvalue)
#' arg = list(
#'   nG = 10000,  
#'   pi0 = 0.8,   
#'   group = rep(c(1, 2), each = 13) 
#' )
#' 
#' k <- 2                                   ## number of simulations (defined by user)
#' fdr <- 0.05                              ## the false discovery rate to be controlled
#' mu <- 10                                 ## mean counts in control group for all genes
#' disp <- 0.1                              ## dispersion for all genes
#' logfc <- log(2)                          ## log fold change for up-regulated genes
#' 
#' check.power(arg, k, fdr, mu, disp, logfc)
#'
#'  @export
#' 
check.power <- function(arg, k, fdr, mu, disp, logfc){
  res <- list()
  ## "power" & "fdr" function
  powerfdr.fun <- function(fdr, p){
    V <- sum((p < fdr) & (sim$de == FALSE))
    R <- sum(p < fdr)
    S <- R - V
    power <- S / (arg$nG * (1 - arg$pi0))
    fdr_true <- V / R
    return(c(power, fdr_true))
  }
  
  pow_bh <- fdr_bh <- pow_st <- fdr_st <- rep(0, k)
  for(j in 1:k){
    sim <- sim.counts(arg, mu, disp, logfc)
    cts <- sim$counts
    lib.size <- colSums(cts)
    d <- DGEList(cts, lib.size, group = arg$group)
    d <- calcNormFactors(d)
    design <- model.matrix(~ factor(arg$group))
    y <- voom(d, design, plot=FALSE)       # convert read counts to log2-cpm with associated weights
    fit <- lmFit(y, design)
    fit <- eBayes(fit)
    pvalue <- fit$p.value[, 2]             # pvalue
    
    p_bh <- p.adjust(pvalue, method="BH")  # Benjamini & Hochberg
    pow_bh[j] <- powerfdr.fun(fdr, p_bh)[1]
    fdr_bh[j] <- powerfdr.fun(fdr, p_bh)[2]
  
    p_st <- qvalue(pvalue)$qvalues         # Storey and Tibshirani
    pow_st[j] <- powerfdr.fun(fdr, p_st)[1]
    fdr_st[j] <- powerfdr.fun(fdr, p_st)[2]
  }
    
  ## average power & true fdr over k simulations
  res$pow_bh_ave <- mean(pow_bh)
  res$fdr_bh_ave <- mean(fdr_bh)
  res$pow_st_ave <- mean(pow_st)
  res$fdr_st_ave <- mean(fdr_st)
  return(res)
}