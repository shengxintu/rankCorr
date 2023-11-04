####Internal functions
library(rms)
scores_presid <- function (y, X, lnk){
  k <- table(X)
  N <- length(k)
  ny <- length(table(y))
  na <- ny - 1
  nb <-  N - 1
  npar <- na + nb
  mod <- orm(y ~ X, family = lnk, x = TRUE, y = TRUE, maxit=35)
  alpha <- mod$coeff[1 : na]
  beta <- c(0, mod$coeff[-(1 : na)])
  O <- tapply(y, X, function(s) outer(1 : na, s, ">="))
  Y <- tapply(y, X, I)
  psi <- matrix(0, ncol = npar, nrow = N)
  info.matrix <- matrix(0, ncol = npar, nrow = npar)
  Mu <- list()
  dpresid.dtheta <- list()
  for (i in 1 : N){
    ki <- k[i]
    g <- alpha + beta[i]
    muij <- 1 - mod$trans$cumprob(g)
    tmp <-  mod$trans$deriv(x = g, f = muij)
    tmp <- ifelse(tmp < 1e-16, 1e-16, tmp)
    dmij.dg <- - tmp
    dmij.dalpha <- diag(dmij.dg)
    dmij.dbeta <- matrix(0, ncol = nb, nrow = na)
    if(i > 1) dmij.dbeta[, i - 1] <- dmij.dg
    dmij.dtheta <- cbind(dmij.dalpha, dmij.dbeta)
    diffmuij <- diff(c(0, muij))
    diffmuij <- ifelse(diffmuij < 1e-16, 1e-16, diffmuij)
    s.inv <- matrix(1 / (1 - sum(diffmuij)), ncol = na, nrow = na) + diag(1 / diffmuij)
    s.inv[!is.finite(s.inv)] <- 1e16
    l.inv <- diag(1, ncol = na, nrow = na)
    l.inv[cbind(2:na, 1:(na-1))] <- -1
    vij.inv <- t(l.inv) %*% s.inv %*% l.inv
    tmp <- t(dmij.dtheta) %*% vij.inv
    psi[i,] <- rowSums(apply(O[[i]] - muij, 2, function(l) tmp %*% l))
    info.matrix <- info.matrix + tmp %*% dmij.dtheta * ki
    Mu[[i]] <- matrix(rep(muij, ki), nrow = ki, byrow = T)
    dlowi.dtheta <- dhii.dtheta <- matrix(0, nrow = npar, ncol = ki)
    yi <- Y[[i]]
    for(j in seq_along(yi)){
      yij <- yi[j]
      if (yij == 1) {
        dlowi.dtheta[,j] <- 0
      } else {
        dlowi.dtheta[,j] <- c(dmij.dalpha[yij - 1,], dmij.dbeta[yij - 1, ])
      }
      if (yij == ny) {
        dhii.dtheta[,j] <- 0
      } else {
        dhii.dtheta[,j] <- - c(dmij.dalpha[yij,], dmij.dbeta[yij,])
      }
    }
    dpresid.dtheta[[i]] <- dlowi.dtheta - dhii.dtheta
  }
  Mu <- do.call(rbind, Mu)
  Ny <- length(y)
  low <- cbind(0, Mu)[cbind(1 : Ny, y)]
  hi <- cbind(1 - Mu, 0)[cbind(1 : Ny, y)]
  presid <- low - hi
  result <- list(psi = psi,
                 dpsi.dtheta = - info.matrix,
                 presid = presid,
                 dpresid.dtheta = dpresid.dtheta,
                 beta = beta)
  if(length(unique(k)) != 1) result$est.mean <- predict(mod)
  return(result)
}

getCI <- function(ts, v, fisher, ci=0.95){
  if(!fisher){
    lower <- ts - abs(qnorm(0.5*(1-ci)))*sqrt(v)
    upper <- ts + abs(qnorm(0.5*(1-ci)))*sqrt(v)
  } else {
    ts_f <- log((1+ts)/(1-ts))
    v_f <- v*(2/(1-ts^2))^2
    lower_f <- ts_f - abs(qnorm(0.5*(1-ci)))*sqrt(v_f)
    upper_f <- ts_f + abs(qnorm(0.5*(1-ci)))*sqrt(v_f)
    lower <- (exp(lower_f)-1)/(1+exp(lower_f))
    upper <- (exp(upper_f)-1)/(1+exp(upper_f))
  }
  return(c(lower, upper))
}

########Within-cluster correlation
cor_rw <- function(xresid, yresid, x.psi, y.psi,
                   x.dpsi.dtheta, y.dpsi.dtheta,
                   x.dpresid.dtheta, y.dpresid.dtheta,
                   cls, wi, wij, conf.int, fisher){
  cls.xresid <- tapply(xresid, cls, mean); cls.yresid <- tapply(yresid, cls, mean)
  cls.xresid2 <- tapply(xresid ^ 2, cls, mean); cls.yresid2 <- tapply(yresid ^ 2, cls, mean)
  cls.xyresid <- tapply(xresid * yresid, cls, mean)
  mean.xresid <- sum(cls.xresid * wi); mean.yresid <- sum(cls.yresid * wi)
  mean.xyresid <- sum(cls.xyresid * wi)
  mean.xresid2 <- sum(cls.xresid2 * wi); mean.yresid2 <- sum(cls.yresid2 * wi)
  rw <- sum(wij * (xresid - mean.xresid) * (yresid - mean.yresid)) / sqrt(sum((xresid - mean.xresid)^2 * wij)* sum((yresid - mean.yresid)^2 * wij))
  rw.output <- c(rw, rep(NA, 3))
  names(rw.output) <- c("Estimate", "SE", "Lower", "Upper")
  if(conf.int){
    bigpsi <- cbind(x.psi,
                    y.psi,
                    mean.xresid - cls.xresid,
                    mean.yresid - cls.yresid,
                    mean.xyresid - cls.xyresid,
                    mean.xresid2 - cls.xresid2,
                    mean.yresid2 - cls.yresid2,
                    0)
    npar.x <- dim(x.psi)[2]
    npar.y <- dim(y.psi)[2]
    Ntheta <- npar.x + npar.y + 6
    N <- dim(x.psi)[1]
    A <- matrix(0, Ntheta, Ntheta)
    A[1 : npar.x, 1 : npar.x] <- x.dpsi.dtheta
    A[npar.x + (1 : npar.y), npar.x + (1 : npar.y)] <- y.dpsi.dtheta
    A[Ntheta - 6 + (1 : 6), Ntheta - 6 + (1 : 6)] <-  diag(N, 6)
    l.xresid <- tapply(xresid, cls, I)
    l.yresid <- tapply(yresid, cls, I)
    k <- table(cls)
    xcls.dpresid.dtheta <- matrix(0, npar.x, N)
    ycls.dpresid.dtheta <- matrix(0, npar.y, N)
    dxpresid.xresid <- dxpresid.yresid <- matrix(0, npar.x, N)
    dypresid.xresid <- dypresid.yresid <- matrix(0, npar.y, N)
    for(i in 1:N){
      xpi <- x.dpresid.dtheta[[i]]
      ypi <- y.dpresid.dtheta[[i]]
      ki <- k[i]
      xi <- matrix(l.xresid[[i]], ncol = 1)
      yi <- matrix(l.yresid[[i]], ncol = 1)
      xcls.dpresid.dtheta[, i] <- rowMeans(xpi)
      ycls.dpresid.dtheta[, i] <- rowMeans(ypi)
      dxpresid.yresid[, i] <- xpi %*% yi / ki
      dypresid.xresid[, i] <- ypi %*% xi / ki
      dxpresid.xresid[, i] <- xpi %*% (2 * xi) / ki
      dypresid.yresid[, i] <- ypi %*% (2 * yi) / ki
    }
    bigpartial <- rbind(c(rowSums(xcls.dpresid.dtheta), rep(0, npar.y)),
                        c(rep(0, npar.x), rowSums(ycls.dpresid.dtheta)),
                        c(rowSums(dxpresid.yresid), rowSums(dypresid.xresid)),
                        c(rowSums(dxpresid.xresid), rep(0, npar.y)),
                        c(rep(0, npar.x), rowSums(dypresid.yresid)))
    A[Ntheta - 6 + (1:5), 1 : (npar.x + npar.y)] <- - bigpartial
    g.rw <- mean.xyresid - mean.xresid * mean.yresid
    var.xresid <- mean.xresid2 - mean.xresid ^ 2
    var.yresid <- mean.yresid2 - mean.yresid ^ 2
    varprod <- var.xresid * var.yresid
    revsvp <- 1 / sqrt(varprod)
    drw.dvarprod <- - g.rw / 2 * revsvp ^ 3
    smallpartial <- N * c(- mean.yresid * revsvp + drw.dvarprod * (-2 * mean.xresid * var.yresid),
                          - mean.xresid * revsvp + drw.dvarprod * (-2 * mean.yresid * var.xresid),
                          revsvp, drw.dvarprod * var.yresid,
                          drw.dvarprod * var.xresid)
    A[Ntheta, Ntheta - 6 + (1 : 5)] <- - smallpartial
    SS <- solve(A, t(bigpsi), tol=1e-23)
    var.theta <- tcrossprod(SS, SS)
    var.rw <- var.theta[Ntheta, Ntheta]
    rw.output[c("Lower", "Upper")] <- getCI(ts = rw, v = var.rw, fisher = fisher, ci = conf.int)
    rw.output["SE"] <- sqrt(var.rw)
  }
  return(rw.output)
}

emp_CDF <- function(x, pij, tol = 1e-7) {
  vapply(x, function(i) {
    c1 <- x - i
    isEq <- abs(c1) < tol
    isLT <- (!isEq & c1 < 0)
    sum(isLT * pij + (isEq * pij) / 2)
  }, numeric(1))
}

########Cluster-median-based between-cluster correlation
cor_rb <- function(xbeta, ybeta, x.psi, y.psi,
                        x.dpsi.dtheta, y.dpsi.dtheta, wi, conf.int, fisher){
  xef <- emp_CDF(xbeta, wi)
  yef <- emp_CDF(ybeta, wi)
  mean.xef <- sum(wi * xef); mean.yef <- sum(wi * yef)
  rb <- sum(wi * (xef - mean.xef) * (yef - mean.yef)) / sqrt(sum((xef - mean.xef)^2 * wi)* sum((yef - mean.yef)^2 * wi))
  rb.output <- c(rb, rep(NA, 3))
  names(rb.output) <- c("Estimate", "SE", "Lower", "Upper")
  if(conf.int){
    mean.xb <- sum(wi * xbeta); mean.yb <- sum(wi * ybeta)
    mean.xb2 <- sum(wi * xbeta ^ 2); mean.yb2 <- sum(wi * ybeta ^ 2)
    sd.xb <- sqrt(mean.xb2 - mean.xb ^ 2); sd.yb <- sqrt(mean.yb2 - mean.yb ^ 2)
    xnf <- pnorm(xbeta, mean = mean.xb, sd = sd.xb)
    ynf <- pnorm(ybeta, mean = mean.yb, sd = sd.yb)
    mean.xynf <- sum(xnf * ynf * wi); mean.xnf <-  sum(xnf * wi); mean.ynf <- sum(ynf * wi)
    mean.xnf2 <- sum(xnf ^ 2 * wi); mean.ynf2 <- sum(ynf ^ 2 * wi)
    bigpsi <- cbind(x.psi, y.psi,
                    mean.xb - xbeta,
                    mean.yb - ybeta,
                    mean.xb2 - xbeta ^ 2,
                    mean.yb2 - ybeta ^ 2,
                    mean.xnf - xnf,
                    mean.ynf - ynf,
                    mean.xynf - xnf * ynf,
                    mean.xnf2 - xnf ^ 2,
                    mean.ynf2 - ynf ^ 2,
                    0)
    npar.x <- dim(x.psi)[2]
    npar.y <- dim(y.psi)[2]
    Ntheta <- npar.x + npar.y + 10
    N <- length(xbeta)
    nb <- N - 1
    na.x <- npar.x - nb
    na.y <- npar.y - nb
    A <- matrix(0, Ntheta, Ntheta)
    A[1 : npar.x, 1 : npar.x] <- x.dpsi.dtheta
    A[npar.x + (1 : npar.y), npar.x + (1 : npar.y)] <- y.dpsi.dtheta
    A[Ntheta - 10 + (1 : 10), Ntheta - 10 + (1 : 10)] <-  diag(N, 10)
    A[Ntheta - 10 + (1 : 4), na.x + (1 : nb)] <- - rbind(rep(1, nb), rep(0, nb), 2 * xbeta[-1], rep(0, nb))
    A[Ntheta - 10 + (1 : 4), npar.x + na.y + (1 : nb)] <- - rbind(rep(0, nb), rep(1, nb), rep(0, nb), 2 * ybeta[-1])
    dxnf.dxb <- dnorm(xbeta, mean = mean.xb, sd = sd.xb)
    dynf.dyb <- dnorm(ybeta, mean = mean.yb, sd = sd.yb)
    A[Ntheta - 10 + (5 : 9), na.x + (1 : nb)] <- - rbind(dxnf.dxb[-1], rep(0, nb), dxnf.dxb[-1] * ynf[-1], 2 * dxnf.dxb[-1] * xnf[-1], rep(0, nb))
    A[Ntheta - 10 + (5 : 9), npar.x + na.y + (1 : nb)] <- - rbind(rep(0, nb), dynf.dyb[-1], dynf.dyb[-1] * xnf[-1], rep(0, nb), 2 * dynf.dyb[-1] * ynf[-1])
    dxnf.dxb.s <- dnorm((xbeta - mean.xb) / sd.xb, 0, 1)
    dynf.dyb.s <- dnorm((ybeta - mean.yb) / sd.yb, 0, 1)
    dxnf.dtheta1 <- dxnf.dxb.s * (- 1 / sd.xb + (xbeta - mean.xb) * mean.xb / (sd.xb ^ 3))
    dxnf.dtheta2 <- dxnf.dxb.s * (mean.xb - xbeta) / (2 * sd.xb ^ 3)
    dynf.dtheta1 <- dynf.dyb.s * (- 1 / sd.yb + (ybeta - mean.yb) * mean.yb / (sd.yb ^ 3))
    dynf.dtheta2 <- dynf.dyb.s * (mean.yb - ybeta) / (2 * sd.yb ^ 3)
    A[Ntheta - 10 + (5 : 9), npar.x + npar.y + (1 : 4)] <- - rbind(c(sum(dxnf.dtheta1), 0, sum(dxnf.dtheta2), 0),
                                                                   c(0, sum(dynf.dtheta1), 0, sum(dynf.dtheta2)),
                                                                   c(sum(dxnf.dtheta1 * ynf), sum(dynf.dtheta1 * xnf), sum(dxnf.dtheta2 * ynf), sum(dynf.dtheta2 * xnf)),
                                                                   c(sum(dxnf.dtheta1 * xnf * 2), 0, sum(dxnf.dtheta2 * xnf * 2), 0),
                                                                   c(0, sum(dynf.dtheta1 * ynf * 2), 0, sum(dynf.dtheta2 * ynf * 2)))
    var.xf <- mean.xnf2 - mean.xnf ^ 2
    var.yf <- mean.ynf2 - mean.ynf ^ 2
    g.rw <- mean.xynf - mean.xnf * mean.ynf
    varprod <- var.xf * var.yf
    revsvp <- 1 / sqrt(varprod)
    drw.dvarprod <- - g.rw / 2 * revsvp ^ 3
    smallpartial <- N * c(- mean.ynf * revsvp + drw.dvarprod * (-2 * mean.xnf * var.yf),
                          - mean.xnf * revsvp + drw.dvarprod * (-2 * mean.ynf * var.xf),
                          revsvp, drw.dvarprod * var.yf,
                          drw.dvarprod * var.xf)
    A[Ntheta, Ntheta - 10 + (5 : 9)] <- - smallpartial
    SS <- solve(A, t(bigpsi), tol=1e-23)
    var.theta <- tcrossprod(SS, SS)
    var.rb <- var.theta[Ntheta, Ntheta]
    rb.output[c("Lower", "Upper")] <- getCI(ts = rb, v = var.rb, fisher = fisher, ci = conf.int)
    rb.output["SE"] <- sqrt(var.rb)
  }
  return(rb.output)
}

emp_CDF_xy <- function(x, y, wij, tol = 1e-7){
  n <- length(x)
  xy <- cbind(x, y)
  est <- apply(xy, 1, function(i){
    xc <- x - i[1]
    yc <- y - i[2]
    x.isEq <- abs(xc) < tol
    x.isLT <- (!x.isEq & xc < 0)
    y.isEq <- abs(yc) < tol
    y.isLT <- (!y.isEq & yc < 0)
    x.l <- x.isLT * wij + (x.isEq * wij) / 2
    y.l <- y.isLT * wij + (y.isEq * wij) / 2
    c(sum(x.l), sum(y.l), x.l / sqrt(wij), y.l / sqrt(wij))
  })
  x.est <- est[3 : (n + 2), ]
  y.est <- est[(n + 3) : (2 * n + 2), ]
  xy.est <- apply(y.est, 2, function(i) colSums(i * x.est))
  return(list(xcdf = est[1,],
              ycdf = est[2,],
              xycdf = xy.est,
              xl = x.est / sqrt(wij),
              yl = y.est / sqrt(wij)))
}

########Total correlation
cor_rt <- function(x, y, cls, wij, conf.int, fisher){
  n <- length(x)
  n.cluster <- length(unique(cls))
  CDF.est <- emp_CDF_xy(x, y, wij)
  x.cdf <- CDF.est$xcdf
  y.cdf <- CDF.est$ycdf
  x.mcdf <- sum(x.cdf * wij)
  y.mcdf <- sum(y.cdf * wij)
  cov.est <- sum(wij * (x.cdf - x.mcdf) * (y.cdf - y.mcdf))
  vx.est <- sum(wij * (x.cdf - x.mcdf) ^ 2)
  vy.est <- sum(wij * (y.cdf - y.mcdf) ^ 2)
  rt <- cov.est / sqrt(vx.est * vy.est)
  rt.output <- c(rt, rep(NA, 3))
  names(rt.output) <- c("Estimate", "SE", "Lower", "Upper")
  if(conf.int){
    xs <- rep(x, n)
    ys <- rep(y, each = n)
    xy.cdf <- CDF.est$xycdf
    h <- (2 * x.cdf * y.cdf + 2 * rowMeans(xy.cdf) + 2 * colMeans(xy.cdf) + 2 - 2 * x.cdf - 2 * y.cdf) / 6 - 1 / 4 - cov.est
    hbar <- tapply(h, cls, mean)
    hx <- (x.cdf ^ 2 + 2 * rowMeans(sweep(CDF.est$xl, 2, x.cdf, "*"))) / 3 - 1 / 4 - vx.est
    hy <- (y.cdf ^ 2 + 2 * rowMeans(sweep(CDF.est$yl, 2, y.cdf, "*"))) / 3 - 1 / 4 - vy.est
    hxbar <- tapply(hx, cls, mean)
    hybar <- tapply(hy, cls, mean)
    hbar.mat <- cbind(hbar, hxbar, hybar)
    var.h.mat <- 9 * t(hbar.mat)  %*% hbar.mat / (n.cluster ^ 2)
    g.delta <- rt * c(1 / cov.est, - 1 / (2 * vx.est), - 1 / (2 * vy.est))
    var.rt <- g.delta %*% var.h.mat %*% g.delta
    rt.output[c("Lower", "Upper")] <- getCI(ts = rt, v = var.rt, fisher = fisher, ci = conf.int)
    rt.output["SE"] <- sqrt(var.rt)
  }
  return(rt.output)
}

#########Weights of the between- and within-cluster correlations in the approximated linear relationship
linear_weights <- function(x, cluster, weights, negICC = FALSE){
  est <- rankICC::rankICC(x, cluster, weights = weights, conf.int = FALSE)
  rankicc <- est['rankICC']
  if(rankicc < 0){
    ki <- tabulate(cluster)
    n.cluster <- length(unique(cluster))
    n.obs <- length(x)
    wij <- ifelse(weights[1] == "obs", 1 / n.obs, rep(1 / n.cluster / ki, ki))
    wi <- ifelse(weights[1] == "obs", ki / n.obs, 1 / n.cluster)
    ef <- emp_CDF(x, wij)
    avg <- sum(ef * wij)
    tv <- sum((ef - avg)^2 * wij)
    avgi1 <- tapply(ef, cluster, mean)
    est2 <- sum(wi * ((avgi1 - avg)^2)) / tv
    rankicc <- c(rankicc, est2)
  }
  return(rankicc)
}

########Approx-based between-cluster correlation
cor_rb_approx <- function(x, y, cluster, rw, rt, weights,
                          conf.int = 0.95, fisher = FALSE, rb){
  irx <- linear_weights(x, cluster, weights)
  iry <- linear_weights(y, cluster, weights)
  if(length(irx) > 1 | length(iry) > 1){
    rb.approx <- (rt - sqrt((1-irx[length(irx)])*(1-iry[length(iry)])) * rw)/sqrt(irx[length(irx)]*iry[length(iry)])
  }
  else rb.approx <- (rt - sqrt((1-irx[1])*(1-iry[1])) * rw)/sqrt(irx[1]*iry[1])
  rb.approx.est <- c(rb.approx, rep(NA, 3))
  names(rb.approx.est) <- c("Estimate", "SE", "Lower", "Upper")
  rankicc <- c("x" = unname(irx[1]), "y" = unname(iry[1]))
  if(conf.int){
    rb.approx.est[c("Lower", "Upper")] <- getCI(ts = rb.approx, v = rb['SE']^2, fisher = fisher, ci = conf.int)
    rb.approx.est["SE"] <- rb['SE']
  }
  rb.approx.output <- list(rb.approx.est = rb.approx.est,
                    rankicc = rankicc)
  return(rb.approx.output)
}






