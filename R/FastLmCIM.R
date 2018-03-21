#############
# FastLmCIM #
#############

# function computing the significance of a QTL term in a CIM scan. The function
# use .lm.fit function which help to speed up the computation.

FastLmCIM <- function(y, Xc, Xq){

  index <- complete.cases(cbind(y, Xc, Xq))

  y <- y[index]
  Xc <- Xc[index, , drop = FALSE]
  Xq <- Xq[index, , drop = FALSE]

  X <- cbind(Xc, Xq)
  y.bar <- mean(y)
  n.tot <- length(y)
  n.al <- dim(Xq)[2]

  m.fit <- .lm.fit(X, y)

  # re-order the coefficients

  coefs <- m.fit$coefficients[match(1:dim(X)[2], m.fit$pivot)]

  n.res <- n.tot - m.fit$rank
  n.al <- sum(rev(coefs)[1:n.al] != 0)

  # compute ANOVA
  ###############

  SSE <- sum(m.fit$residuals^2)

  # SSR.all

  y.hat <- X %*% coefs
  SSR.tot <- sum((y.hat - y.bar)^2)

  # SSR.c

  m.fit <- .lm.fit(Xc, y)
  coefs <- m.fit$coefficients[match(1:dim(Xc)[2], m.fit$pivot)]

  y.hat <- Xc %*% coefs
  SSR.c <- sum((y.hat - y.bar)^2)


  # compute the F.stat

  pval <- pf(q =((SSR.tot - SSR.c)/n.al)/(SSE/n.res), df1 = n.al,
             df2 = n.res, lower.tail = FALSE)

  -log10(pval)

}
