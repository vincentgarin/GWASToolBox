################
# GWAS_GenABEL #
################

#' GWAS using GenAbel functions
#'
#' Perform a SNP genome wide association study using optional principal
#' components terms or a kinship matrix to correct for the genetic background.
#' The function is a wrapper for the \code{GenAbel} package functions
#' (Covarrubias-Pazaran, 2016).
#'
#' The model can be fitted using the kinship containing all markers or removing
#' the markers of the scanned chromosome (\code{K_i = TRUE}).
#'
#' @param gwaa \code{gwaa.data} object.
#'
#' @param trait \code{Numerical} or \code{character} indicator to specify which
#' trait of the gp object should be used. Default = 1.
#'
#' @param model \code{Character} argument specifying the type of model that is
#' used to perform the GWAS analysis: 1) 'null' for a simple model with only the
#' QTL position; 2) 'PC' model with genetic background correction using the
#' \code{nPC} first principal components; 3) 'kinship', model with kinship
#' genetic background correction. Default = 'kinship'.
#'
#' @param nPC Optional matrix of principal components for genetic background
#' correction. Default = NULL.
#'
#' @param K_i \code{Logical} specifying if the kinship correction should be done
#' by removing the markers of the scanned chromosome. Default = TRUE.
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
#' Aulchenko, Y. S., Ripke, S., Isaacs, A., & Van Duijn, C. M. (2007). GenABEL:
#' an R library for genome-wide association analysis. Bioinformatics, 23(10),
#' 1294-1296.
#'
#' Aulchenko, Y. S., De Koning, D. J., & Haley, C. (2007). Genomewide rapid
#' association using mixed model and regression: a fast and simple method for
#' genomewide pedigree-based quantitative trait loci association analysis.
#' Genetics, 177(1), 577-585.
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
# gwaa <- df
# trait <- 1
# model <- 'kinship'
# nPC = 2
# K_i = FALSE
# verbose = FALSE


GWAS_GenABEL <- function(gwaa, trait = 1, model = 'kinship', nPC = NULL,
                         K_i = TRUE, verbose = FALSE){

  # check gwaa

  if(!is_gwaa(gwaa)){

    stop("The object given in gwaa is not a gwaa.data object.")

  }

  # check_weights(weights = weights, gp = gp)

  check_trait_gwaa(trait = trait, gwaa = gwaa)

  # reconstruct the map

  chr <- gwaa@gtdata@chromosome
  pos <- gwaa@gtdata@map
  map <- data.frame(attr(pos, "names"), as.character(chr), pos,
                    stringsAsFactors = FALSE)
  map[, 2] <- as.numeric(map[, 2])

  colnames(map) <- c("mk.id", "chr", "pos")

  # check that the selection of marker is present in the map

  # if( sum(!(mk.sel %in% map[, 1])) !=0 ){
  #
  #   prob.mk <- paste(mk.sel[!(mk.sel %in% map[, 1])], collapse = ", ")
  #
  #   stop(paste("the following markers:", prob.mk, "are present in mk.sel but",
  #              "not in the map."))
  #
  # }

  ############ end checks

  # get the phenotypic values

  pheno <- gwaa@phdata
  pheno <- pheno[3:dim(pheno)[2]]
  ph.val <- pheno[, trait]

  if(model == 'null'){

    ans <- GenABEL::qtscore(formula = ph.val, data =  gwaa,
                            trait.type = "gaussian")

    G_res <- data.frame(map, -log10(ans@results$P1df), stringsAsFactors = FALSE)
    colnames(G_res) <- c("mk.id", "Chrom", "Position", "p.val")
    class(G_res) <- c("data.frame", "G_res")

  } else if (model == 'PC'){

    # computation of the PC

    # PC <- prcomp(gp$geno, center = TRUE, scale. = TRUE)
    # PC <- PC$x[, 1:nPC]

    gkin <- ibs(data = gwaa, w = "freq")
    diag(gkin) <- hom(gwaa)$Var

    ans <- GenABEL::egscore(formula = ph.val, data = gwaa, kinship.matrix = gkin,
                            naxes = nPC)

    G_res <- data.frame(map, -log10(ans@results$P1df), stringsAsFactors = FALSE)
    colnames(G_res) <- c("mk.id", "Chrom", "Position", "p.val")
    class(G_res) <- c("data.frame", "G_res")


  } else if (model == 'kinship'){


    if(!K_i){

      gkin <- ibs(data = gwaa, w = "freq")

      gkin <- gkin + (diag(length(ph.val)) * 0.01)
      eg_val <- eigen(gkin)
      eg_dec <- solve(gkin)

      # perform polygenic analysis

      h2ht <- polygenic(ph.val, kin = gkin, data = gwaa, trait.type = "gaussian",
                        quiet = !verbose)

      # maybe need to add a small digit to the kinship...

      ans <- grammar(polyObject = h2ht, data = gwaa, method = "gamma")

      G_res <- data.frame(map, -log10(ans@results$P1df), stringsAsFactors = FALSE)
      colnames(G_res) <- c("mk.id", "Chrom", "Position", "p.val")
      class(G_res) <- c("data.frame", "G_res")

    } else {

        p.val <- c()
        n.chr <- length(unique(map[, 2]))
        chr.id <- unique(map[, 2])

      for(i in 1:n.chr){

          mk_chr_i <- map[map[, 2] == chr.id[i], 1]

          # adapt the list of selected markers

          mk.sel_i <- map[map[, 2] != chr.id[i], 1]

          # obtain a reduced genotype matrix

          geno_i <- gp$geno[, colnames(gp$geno) %in% mk_chr_i]

          gkin <- ibs(data = gwaa, snpsubset = mk.sel_i, w = "freq")

          # perform polygenic analysis

          h2ht <- polygenic(ph.val, kin = gkin, data = gwaa[, mk_chr_i],
                            trait.type = "gaussian", quiet = !verbose)

          # maybe need to add a small digit to the kinship...

          ans <- grammar(polyObject = h2ht, data = gwaa[, mk_chr_i],
                         method = "gamma")

          p.val <- c(p.val, -log10(ans@results$P1df))

          }

        G_res <- data.frame(map, p.val, stringsAsFactors = FALSE)
        colnames(G_res) <- c("mk.id", "Chrom", "Position", "p.val")
        class(G_res) <- c("data.frame", "G_res")

    }

  }

  return(G_res)

}
