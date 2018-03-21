#############
# FastLmSIM #
#############

# function computing the significance of a QTL term in a SIM scan. The function
# use .lm.fit function which help to speed up the computation.

FastLmSIM <- function(y, Xc, Xq){

  index <- complete.cases(cbind(y, Xc, Xq))

  y <- y[index]
  Xc <- Xc[index, , drop = FALSE]
  Xq <- Xq[index, , drop = FALSE]

  X <- cbind(Xc, Xq)
  y.bar <- mean(y)
  n.tot <- length(y)
  n.al <- dim(Xq)[2]

  m.fit <- .lm.fit(X, y)

  n.res <- n.tot - m.fit$rank
  coefs <- m.fit$coefficients[match(1:dim(X)[2], m.fit$pivot)]
  n.al <- sum(rev(coefs)[1:n.al] != 0)

  # compute ANOVA
  ###############

  SSE <- sum(m.fit$residuals^2)

  # SSR.all

  y.hat <- X %*% coefs
  SSR.tot <- sum((y.hat - y.bar)^2)

  # SSR.c

  y.hat <- rep(.lm.fit(Xc, y)$coefficients, colSums(Xc))
  SSR.c <- sum((y.hat - y.bar)^2)


  # compute the F.stat

  pval <- pf(q =((SSR.tot - SSR.c)/n.al)/(SSE/n.res), df1 = n.al,
             df2 = n.res, lower.tail = FALSE)

  -log10(pval)

}
