#############
# GWAS_scan #
#############

#' GWAS scan
#'
#' Perform a SNP genome wide association study using optional principal
#' components terms or a kinship matrix to correct for the genetic background.
#'
#' The function is a wrapper for the \code{GWAS} function from the \code{sommer}
#' package (Covarrubias-Pazaran, 2016). The model is fitted using the EMMA
#' algorithm proposed by Kang et al. (2008).
#'
#' By default, the kinship matrix is computed using the method of Astle and
#' Balding (2009). It is possible to compute a linkage disequilibrium adjusted
#' kinship (LDAK) kinship matrix using the method of Speed et al. (2012) by
#' introducing the weights computed with the function
#' \code{LDAK_weights()} (Not available for the moment).
#' The model can be fitted using the kinship containing all markers or removing
#' the markers of the scanned chromosome (\code{K_i = TRUE}).
#'
#' @param gp \code{gpData} object with elements geno coded 0 1 2, map with
#' marker position in bp or cM and phenotype.
#'
#' @param trait \code{Numerical} or \code{character} indicator to specify which
#' trait of the gp object should be used. Default = 1.
#'
#' @param model \code{Character} argument specifying the type of model that is
#' used to perform the GWAS analysis: 1) 'null' for a simple model with only the
#' QTL position; 2) 'cofactors' for a model including QTL and cofactors;
#' 3) 'PC' model with genetic background correction using the \code{nPC}
#' first principal components; 4) 'kinship', model with kinship
#' genetic background correction. Default = 'kinship'.
#'
#' @param thre_cof \code{Numeric} -log10(pval) above which a cofactor is
#' selected. Default = 6.
#'
#' @param win_cof \code{Numeric} value indicating the minimum distance between
#' two cofactors positions.
#'
#' @param nPC Optional number of principal components for genetic background
#' correction if the 'PC' model is selected. Default = 2.
#'
#' @param K_i \code{Logical} specifying if the kinship correction should be done
#' by removing the markers of the scanned chromosome. Default = TRUE.
#'
#' @param weights (Not available for the moment) object of class \code{LD_wgh}
#' obtained by the function LDAK_weights() representing a data.frame with two
#' columns: the marker identifier and the LD adjusted weights. These weight will
#' be used to compute a LDAK as defined by Speed et al. (2012) for the genetic
#' background adjustement. Default = NULL.
#'
#' @param power \code{Numerical} value specifying the value of the
#' parameter for marker scores standardization. The column of the marker matrix
#' (X.j) are multiplied by var(X.j)^(power/2). It correspond to alpha in the
#' formula. Default = -1.
#'
#' @param mk.sel \code{Character vector} specifying a list of marker to use
#' for the kinship matrix computation. By default, the function use all markers
#' of the \code{gp}.
#'
#' @param verbose \code{Logical} indicating if function outputs should be printed.
#' Default = FALSE.
#'
#' @return Return:
#'
#' \item{G_res}{Object of class \code{G_res} representing a data.frame with four
#' columns: marker identifier, chromosome, position in cM or bp and
#' -log10(p-value).}
#'
#' @author Vincent Garin
#'
#' @references
#'
#' Astle, W., & Balding, D. J. (2009). Population structure and cryptic
#' relatedness in genetic association studies. Statistical Science, 451-471.
#'
#' Covarrubias-Pazaran G. 2016. Genome assisted prediction of quantitative traits
#' using the R package sommer. PLoSONE 11(6):1-15.
#'
#' Kang, H. M., Zaitlen, N. A., Wade, C. M., Kirby, A., Heckerman, D., Daly,
#' M. J., & Eskin, E. (2008). Efficient control of population structure in model
#' organism association mapping. Genetics, 178(3), 1709-1723.
#'
#' Speed, D., Hemani, G., Johnson, M. R., & Balding, D. J. (2012).
#' Improved heritability estimation from genome-wide SNPs. The American Journal
#' of Human Genetics, 91(6), 1011-1021.
#'
#'@examples
#'
#' Come later
#'
#' @export
#'

# arguments

# source('G:/GWAS_ToolBox/GWASToolBox/R/check_object_format.R')
# source('G:/GWAS_ToolBox/GWASToolBox/R/test_class.R')
#
# gp <- gp
# trait <- 1
# model = 'kinship'
# thre_cof = 6
# win_cof = 70000000/3
# nPC = 2
# K_i = FALSE
# weights = NULL
# power = -1
# mk.sel = NULL
# verbose = FALSE


GWAS_scan <- function(gp, trait = 1, model = "kinship", thre_cof = 6, win_cof,
                      nPC = 2, K_i = TRUE, weights = NULL, power = -1,
                      mk.sel = NULL, verbose = FALSE){

  # check gp, trait and weights

  check_gpData(gp)

  # check_weights(weights = weights, gp = gp)

  check_trait(trait = trait, gp = gp)

  # reconstruct the map

  map <- data.frame(rownames(gp$map), gp$map[, 1:2], stringsAsFactors = FALSE)
  colnames(map) <- c("mk.id", "chr", "pos")

  # check that the selection of marker is present in the map

  if( sum(!(mk.sel %in% map[, 1])) !=0 ){

    prob.mk <- paste(mk.sel[!(mk.sel %in% map[, 1])], collapse = ", ")

    stop(paste("the following markers:", prob.mk, "are present in mk.sel but",
               "not in the map."))

  }

  ############ end checks

  # select the trait

  if(is.numeric(trait)){

    pheno <- gp$pheno[, trait, 1]

  } else {

    trait.names <- attr(gp$pheno, "dimnames")[[2]]
    pheno <- gp$pheno[, 1, which(trait %in% trait.names)]

  }


  if(model == "null"){

    Z1 <- diag(length(pheno))
    ETA <- list(list(Z=Z1))

    ans <- sommer::GWAS(Y = pheno, Z = ETA, M = gp$geno, silent = !verbose,
                        gwas.plots = FALSE)

    G_res <- data.frame(map, ans$M.scores$score[1, ], stringsAsFactors = FALSE)
    colnames(G_res) <- c("mk.id", "Chrom", "Position", "p.val")
    class(G_res) <- c("data.frame", "G_res")

  } else if (model == 'cofactors'){

    # SIM scan

    scan1 <- rep(0, dim(map)[1])

    Xc <- matrix(rep(1, length(pheno)), ncol = 1)

    for(i in 1:dim(map)[1]){

      QTL <- matrix(gp$geno[, i], ncol = 1)

      # scan1[i] <- tryCatch(-log10(anova(lm(pheno ~ QTL))[1, 5]),
      #                      error = function(e) 1)

      scan1[i] <- tryCatch(FastLmSIM(y = pheno, Xc = Xc, Xq = QTL),
                           error = function(e) 0)

    }

    # cofactor selection

    Qprof <- data.frame(map[, 1:2], 1:dim(map)[1], map[, 3], scan1)
    colnames(Qprof) <- c("mk.id", "chr", "pos.id", "pos.cM", "log10pval")
    class(Qprof) <- c("data.frame", "QTLprof")

    cofactors <- QTL_select(Qprof = Qprof, threshold = thre_cof,
                            window = win_cof)

    # CIM scan

    cof_mat <- gp$geno[, map[, 1] %in% cofactors[, 1]]

    cof.part <- apply(X = cofactors[, c(2, 4)], MARGIN = 1, FUN = test_cof,
                      map = Qprof[, 1:4], window = win_cof)

    scan2 <- rep(0, dim(map)[1])

    for(i in 1:dim(map)[1]){

      QTL <- matrix(gp$geno[, i], ncol = 1)

      Xc <- cbind(rep(1, length(pheno)), cof_mat[, cof.part[i, ]])

      # scan2.1[i] <- tryCatch(-log10(anova(lm(pheno ~ -1 + Xc +  QTL))[2, 5]),
      #                        error = function(e) 1)

      scan2[i] <- tryCatch(FastLmCIM(y = pheno, Xc = Xc, Xq = QTL),
                           error = function(e) 0)

    }

    G_res <- data.frame(map, scan2, stringsAsFactors = FALSE)
    colnames(G_res) <- c("mk.id", "Chrom", "Position", "p.val")
    class(G_res) <- c("data.frame", "G_res")


  } else if (model == 'PC'){

    # computation of the PC

    PC <- prcomp(gp$geno, center = TRUE, scale. = TRUE)
    PC <- PC$x[, 1:nPC]

    Z1 <- diag(length(pheno))
    ETA <- list(list(Z=Z1))

    ans <- sommer::GWAS(Y = pheno, X = PC, Z = ETA, M = gp$geno,
                        silent = !verbose, gwas.plots = FALSE)

    G_res <- data.frame(map, ans$M.scores$score[1, ], stringsAsFactors = FALSE)
    colnames(G_res) <- c("mk.id", "Chrom", "Position", "p.val")
    class(G_res) <- c("data.frame", "G_res")


  } else if (model == 'kinship'){


    if(!K_i){ # Use the whole genome for K.

      K <- kin_mat(gp = gp, weights = weights, power = power, mk.sel = mk.sel)

      # Could have a check that the matrix is positive semi-definite

      # eg_val <- eigen(K)
      # eg_dec <- solve(K)

      Z1 <- diag(length(pheno))
      ETA <- list( list(Z=Z1, K=K))
      ans <- sommer::GWAS(Y = pheno, Z = ETA, M = gp$geno, method = "EMMA",
                          silent = !verbose, gwas.plots = FALSE)

      G_res <- data.frame(map, ans$M.scores$score[1, ], stringsAsFactors = FALSE)
      colnames(G_res) <- c("mk.id", "Chrom", "Position", "p.val")
      class(G_res) <- c("data.frame", "G_res")

    } else { # Remove the kth chromosome for the computation of K.

      # K all markers - ith chromosome

      p.val <- c()
      n.chr <- length(unique(map[, 2]))
      chr.id <- unique(map[, 2])

      if(!is.null(mk.sel)){

        mk.sel_temp <- map[map[, 1] %in% mk.sel, ]

      } else { mk.sel_temp <- map }

        for(i in 1:n.chr){

          mk_chr_i <- map[map[, 2] == chr.id[i], 1]

          # adapt the list of selected markers

          mk.sel_i <- mk.sel_temp[mk.sel_temp[, 2] != chr.id[i], 1]

          # obtain a reduced genotype matrix

          geno_i <- gp$geno[, colnames(gp$geno) %in% mk_chr_i]

          K <- kin_mat(gp = gp, weights = weights, power = power,
                           mk.sel = mk.sel_i)

          Z1 <- diag(length(pheno))
          ETA <- list( list(Z=Z1, K=K))
          ans <- sommer::GWAS(Y = pheno, Z = ETA, M = geno_i, method = "EMMA",
                              silent = !verbose, gwas.plots = FALSE)
          p.val <- c(p.val, ans$M.scores$score[1, ])

        }

      G_res <- data.frame(map, p.val)
      colnames(G_res) <- c("mk.id", "Chrom", "Position", "p.val")
      class(G_res) <- c("data.frame", "G_res")

    }

  }


  return(G_res)

}
