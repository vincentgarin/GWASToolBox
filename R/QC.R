###################
# Quality control #
###################

#' Quality control
#'
#' Function performing a quality control on...
#'
#' @param geno Genotype marker score matrix in HapMap format.
#'
#' @param pheno Data.frame with the first column being character genotype
#' identifier and the next column being numeric phenotypic values.
#'
#' @param Mkmiss Maximum percentage of phenotypic values in the markers.
#' Default = 0.1.
#'
#' #' @param Genmiss Maximum percentage of phenotypic values in the genotypes.
#' Default = 0.1.
#'
#' @param MAFlim Critical minor allele frequency under which markers are removed.
#' Default = 0.05.
#'
#' @param pheno_imp Method for imputation of the missing phenotypic values.
#' One of : 'no', 'random', 'mean'. Default = 'mean'.
#'
#' @param geno_imp Method for imputation of the missing marker scores values.
#' For the moment only 'random'.
#'
#' @return Return:
#'
#' \code{List} with two objects
#'
#' \item{gp}{a gpData object containing all data processed data (genotype, map,
#' phenotype) that can directly be used in the function GWAS_scan()}
#'
#' \item{genoHapMap}{The sorted genotype marker matrix in HapMap format
#' that can be used for GWAS analysis in TASSEL using the function GWAS_TASSEL().}
#'
#' @export
#'

# geno <- geno
# pheno <- pheno
# MkMiss <- 0.1
# GenMiss <- 0.1
# MAFlim <- 0.05
# pheno_imp <- 'mean' # one of: 'no', 'random', 'mean'.
# geno_imp <- 'random'

QC <- function(geno, pheno, MkMiss = 0.1, GenMiss = 0.1, MAFlim = 0.05,
               pheno_imp = 'mean', geno_imp = 'random'){

  # 1. Genotype data
  ##################

  ### 1.1 check column names of the genotype identifier

  ref_nm_geno <- c("rs#", "alleles", "chrom", "pos", "strand", "assembly#",
                   "center", "protLSID", "assayLSID", "panel", "QCcode")

  if(!identical(colnames(geno)[1:11], ref_nm_geno)){

    stop(paste('The column names 1 to 11 of the genotype files are not',
               'correct. Please use:', paste(ref_nm_geno, collapse = ', ')))

  }

  if(!is.numeric(geno$chrom)){

    stop('The chromosome indicator (chrom) is not numeric.')

  }

  if(!is.numeric(geno$pos)){

    stop('The marker position (pos) indicator is not numeric.')

  }

  ### 1.2 check genotype identifier matching

  # Keep the original list of genotypes for later

  geno_id_gn <- colnames(geno)[12:dim(geno)[2]]

  geno_id_ph <- pheno[, 1]

  if(!identical(geno_id_gn, geno_id_ph)){

    inter.geno.pheno <- intersect(geno_id_gn, geno_id_ph)

    if(length(inter.geno.pheno) == 0) {

      stop('The genotype identifiers in the geno and pheno files do not match.')

    }

    geno.only <- setdiff(geno_id_gn, geno_id_ph)
    pheno.only <- setdiff(geno_id_ph, geno_id_gn)

    if(length(geno.only) > 0){

      cat(paste('The following genotypes:', paste(geno.only, collapse = ', '),
                'are only present in the genotype file.'), '\n')
    }

    if(length(pheno.only) > 0){

      cat(paste('The following genotypes:', paste(pheno.only, collapse = ', '),
                'are only present in the phenotype file.'), '\n')
    }

    # Keep only the genotypes and phenotypes presents

    ind_geno <- sort(which(geno_id_gn %in% inter.geno.pheno) + 11)

    geno <- geno[, c(1:11, ind_geno)]
    pheno <- pheno[pheno[, 1] %in% inter.geno.pheno, ]
    rownames(pheno) <- pheno[, 1]
    pheno <- pheno[colnames(geno)[12:dim(geno)[2]], ]

  }

  ### 1.3 check that markers are scored correctly

  mk_mat <- geno[, 12:dim(geno)[2]]
  mk_mat <- as.matrix(mk_mat)

  if(!is.character(mk_mat)){

    stop('The marker scores are not characters.')

  }

  mk_sc <- unique(c(mk_mat))

  ref_sc <- c('AA', 'CC', 'GG', 'TT', 'AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA',
              'GC', 'GT', 'TA', 'TC', 'TG', 'NN')

  mk_sc_prob <- mk_sc[!(mk_sc %in% ref_sc)]

  if(!(length(mk_sc_prob) == 0)){

    tot_mk <- dim(mk_mat)[1] * dim(mk_mat)[2]

    mk_tab <- table(c(mk_mat))

    freq_mk_prob <- (mk_tab[mk_sc_prob]/tot_mk) * 100

    mk_mat[mk_mat %in% mk_sc_prob] <- 'NN'

    mk_sc_message <- paste('The following unauthorised mk scores:',
                           paste(mk_sc_prob, collapse = ', '), 'representing respectively',
                           paste(round(freq_mk_prob, 4), collapse = ', '), '% of the',
                           'marker scores, have been replaced by missing values.')


    cat(mk_sc_message, "\n")

  }

  ### 1.4 Check that the markers have maximum 2 alleles

  mk_mat_temp <- t(mk_mat)
  mk_mat_temp[mk_mat_temp == 'NN'] <- NA
  colnames(mk_mat_temp) <- geno[, 1]

  # Value to be recorded

  prob_mk <- c() # Start a list for problematic markers
  prob_gen <- c() # Start a list for problematic genotypes
  n_geno_0 <- dim(mk_mat_temp)[1]
  n_mk_0 <- dim(mk_mat_temp)[2]

  geno_err <- QC_GenotypingError(mk.mat = mk_mat_temp)

  if(length(geno_err) > 0){

    prob_mk <- c(prob_mk, geno_err)

    mk_mat_temp <- mk_mat_temp[, !(colnames(mk_mat_temp) %in% geno_err)]
    geno <- geno[!(geno[, 1]%in% geno_err), ]

  }

  ### 1.5 Check marker missingness

  mk_miss <- QC_missing(mk.mat = mk_mat_temp, threshold = MkMiss)


  if(dim(mk_miss)[1] > 0){

    prob_mk <- c(prob_mk, colnames(mk_mat_temp)[mk_miss[, 2]])

    mk_mat_temp <- mk_mat_temp[, -mk_miss[, 2]]
    geno <- geno[-mk_miss[, 2], ]

  }

  ### 1.6 Check genotype missingness

  Gen_miss <- QC_missing(mk.mat = mk_mat_temp, threshold = GenMiss, MARGIN = 1)

  if(dim(Gen_miss)[1] > 0){

    prob_gen <- c(prob_gen, rownames(mk_mat_temp)[Gen_miss[, 2]])

    mk_mat_temp <- mk_mat_temp[-Gen_miss[, 2], ]
    geno <- geno[, -c(Gen_miss[, 2] + 11)]
    pheno <- pheno[-Gen_miss[, 2], ]

  }

  ### 1.7 Check marker MAF


  mk.MAF <- QC_MAF(mk.mat = mk_mat_temp)
  prob.mk.id <- which(mk.MAF < MAFlim)

  if(length(prob.mk.id) > 0){

    prob_mk <- c(prob_mk, colnames(mk_mat_temp)[mk_miss[, 2]])

    mk_mat_temp <- mk_mat_temp[, -prob.mk.id]
    geno <- geno[-prob.mk.id, ]

  }


  # 2. Phenotype data
  ###################

  ### 2.1 check if the phenotype values are numeric

  check_num <- unlist(lapply(X = pheno[, 2:dim(pheno)[2]],
                             FUN = function(x) is.numeric(x)))

  if(any(!check_num)){

    prob_ph <- colnames(pheno)[2:dim(pheno)[2]][!check_num]

    err_msg <- paste('The phenotypic values are not completely numeric or missing.',
                     'The following traits:', paste(prob_ph, collapse = ', '),
                     'contain some non-numeric information.')

    stop(err_msg)

  }

  ### 2.2 imputation of the missing phenotype values

  if(pheno_imp == 'no'){

    not_comp <- !complete.cases(pheno)

    if(any(not_comp)){

      prob_gen <- c(prob_gen, pheno[not_comp, 1])

      comp_gen <- pheno[!not_comp, 1]

      mk_mat_temp <- mk_mat_temp[comp_gen, ]
      geno <- geno[, -c(which(not_comp) + 11)]
      pheno <- pheno[!not_comp, ]

    }

  } else if (pheno_imp == 'mean'){

    for(i in 2:ncol(pheno)){
      pheno[is.na(pheno[,i]), i] <- mean(pheno[,i], na.rm = TRUE)
    }

  } else if(pheno_imp == 'random'){

    for(i in 2:ncol(pheno)){
      n_miss <- length(pheno[is.na(pheno[,i]), i])
      min_val <- min(pheno[,i], na.rm = TRUE)
      max_val <- max(pheno[,i], na.rm = TRUE)
      pheno[is.na(pheno[,i]), i] <- runif(n = n_miss, min = min_val,
                                          max = max_val)
    }

  }


  # 4. Form the gpdata object
  ###########################

  map <- data.frame(geno[, 1], geno$chrom, geno$pos)
  colnames(map) <- c("mk.id", "chr", "bp")

  map_gp <- data.frame(map[, 2], as.numeric(map[, 3]))
  colnames(map_gp) <- c("chr", "pos")
  rownames(map_gp) <- map[, 1]

  pheno_gp <- pheno[, 2:dim(pheno)[2]]
  rownames(pheno_gp) <- pheno[, 1]

  gp <- create.gpData(pheno = pheno_gp, geno = mk_mat_temp, map = map_gp,
                      map.unit = "bp")

  # 5. gpdata object imputation
  #############################

  # percentage of missing values

  N_miss <- sum(is.na(c(mk_mat_temp)))

  Ntot <- dim(mk_mat_temp)[1] * dim(mk_mat_temp)[2]

  perc_miss <- (N_miss/Ntot) * 100

  gp <- codeGeno(gpData = gp, impute = TRUE, impute.type = geno_imp,
                 nmiss = MkMiss, maf = MAFlim, verbose = TRUE)


  mess_imp <- paste(round(perc_miss, 4),'%', 'of missing values were imputed',
                    'using', geno_imp, 'method.')

  cat(mess_imp, "\n")


  # 5. final information
  ######################

  mess_geno <- paste(dim(gp$geno)[1], 'genotypes out of', n_geno_0,
                     'remain after the quality control.')

  mess_mk <- paste(dim(gp$geno)[2], 'markers out of', n_mk_0,
                   'remain after the quality control.')

  cat(mess_geno, "\n")
  cat(mess_mk)

  # return(list(gp = gp, genoHapMap = geno))

  return(gp)

}
