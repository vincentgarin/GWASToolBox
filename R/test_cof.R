############
# test_cof #
############

test_cof <- function(x, map, window) {

  t1 <- map$chr == as.numeric(x[1])
  t2 <- abs(map$pos.cM - as.numeric(x[2])) < window
  !(t1 & t2)

}
