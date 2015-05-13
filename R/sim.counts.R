#' RNA-seq Count Data Simulation from Negative-Binomial Distribution
#' 
#' This function simulates count data from Negative-Binomial distribution
#' for RNA-seq experiments with given mean, dispersion and log fold change. 
#' A count data matrix is generated.
#' 
#' @import MASS
#' @param arg a list of global parameters to pass into the function, such as total number of genes, 
#' proportion of non-differentially expressed genes and treatment groups. See Details for more information.
#' @param mu a vector (or scalar) of mean counts in control group from which to simulate.
#' @param disp a vector (or scalar) of dispersion parameter from which to simulate.
#' @param logfc a vector (or scalar) of log fold change between treatment group and control group.
#' @param up proportion of up-regulated genes among all differentially expressed genes, the default value is \code{0.5}.
#' @param replace sample with or without replacement from given parameters. See Details for more information.
#' 
#' @details \code{arg = list(nGenes, pi0, group)} where \code{nGenes} is the total number of genes, \code{pi0} 
#' is the proportion of non-differentially expressed genes, and \code{group} is the treatment groups.
#' 
#' If the total number of genes is larger than length of \code{mu} or \code{disp}, 
#' \code{replace} always equals \code{TRUE}.
#' 
#' @return \item{counts}{RNA-seq count data matrix.}
#' @return \item{group}{treatment group vector.}
#' @return \item{lambda0}{mean counts in control group for each gene.}
#' @return \item{phi0}{dispersion parameter for each gene.}
#' @return \item{de}{differentially expressed genes indicator: \code{0} for non-differentially expressed genes, 
#' \code{1} for up-regulated genes, \code{-1} for down-regulated genes.}
#' @return \item{delta}{log fold change for each gene between treatment group and control group.}
#' 
#' @author Ran Bi \email{biran@@iastate.edu}, 
#'         Peng Liu \email{pliu@@iastate.edu}, 
#'         Tim Triche \email{tim.triche@@gmail.com}
#' 
#' @examples
#' arg <- list(nGenes = 10000,   ## total number of genes
#'             pi0 = 0.8,        ## proportion non-differentially-expressed 
#'             group = rep(c(1, 2), each = 3))  ## treatment groups
#' mu <- 10                                     ## mean counts in control group 
#' disp <- 0.1                                  ## dispersion across all genes
#' logfc <- log(2)                              ## log FC difference for DEGs 
#' 
#' sim <- sim.counts(arg, mu, disp, logfc, up = 0.5, replace = TRUE)
#' sim$counts                              ## count data matrix
#' 
#' @export
#' 
sim.counts <- function(arg, mu, disp, logfc, up = 0.5, replace = TRUE){
 
  ## expected false positives
  FP <- arg$nGenes * arg$pi0
  TP <- arg$nGenes - FP 
  
  ## types of true positives
  TP_up <- TP * up
  TP_down <- TP - TP_up 

  de <- c(rep(0, FP), rep(1, TP_up), rep(-1, TP_down))
  de <- de[sample.int(length(de))] ## resample
  
  # h = vector indicating which pseudo-genes to re-simulate
  h <- rep(T, arg$nGenes) 
  counts <- matrix(0, nrow = arg$nGenes, ncol = length(arg$group))
  
  ## log fold change, approximately half positive, half negative
  delta <- logfc * de 
  
  selected_genes <- true_means <- true_disps <- rep(0, arg$nGenes)
  left_genes <- 1:length(mu)
  lambda <- phi <- matrix(0, nrow = arg$nGenes, ncol = length(arg$group))
  
  while(any(h)){
    temp <- sample.int(length(left_genes), sum(h), replace)
    temp <- temp[order(temp)]
    selected_genes[h] <- left_genes[temp]
    if (replace == FALSE){
      left_genes <- left_genes[-temp]
    }
    
    true_means[h] <- mu[selected_genes[h]]
    true_disps[h] <- disp[selected_genes[h]]
    
    lambda[h,] <- matrix(true_means[h], ncol = 1) %*% 
                  matrix(rep(1, length(arg$group)), nrow = 1) * 
                  cbind(matrix(rep(1, sum(h) * length(arg$group) / 2), 
                               ncol = length(arg$group) / 2), 
                        matrix(rep(exp(delta[h]), length(arg$group) / 2), 
                               ncol = length(arg$group) / 2))

    ## mean of counts
    phi[h,] <- matrix(rep(true_disps[h], length(arg$group)), 
                      ncol = length(arg$group))
    ## dispersion of counts
    
    counts[h,] <- rnegbin(sum(h) * length(arg$group), lambda[h,], 1 / phi[h,])
    h <- (rowSums(cpm(counts) > 2) < 3)
    # print(sum(h))
  }
  
  if(any(rowSums(cpm(counts) > 2) < 3 ))
    print("Error: Failed to simulate data: some genes are not expressed.")
  
  list(counts = counts, 
       group = arg$group, 
       lambda0 = lambda[, 1],   # mean counts in control group
       phi0 = phi[, 1],   # dispersion
       de = de,   # DE indicator
       delta = delta  # log fold change
  )
}
