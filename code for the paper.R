########Code for Between- and Within-Cluster Spearman rank correlations
##Application 
#CD4 and CD8 data
demo.d <- read.csv("cd4cd8rat_pts_04jun2015.csv")
labs.d <- read.csv("cd4cd8ratio_labs_04jun2015.csv")
tmp <- labs.d %>% 
  filter(testName == "CD4 COUNT") %>% 
  rename(cd4Count = RESULT_NUMERIC) %>% 
  select(CFAR_PID, AGE_AT_RESULT_DATE, cd4Count)
dat <- labs.d %>% 
  filter(testName == "CD8 COUNT") %>% 
  rename(cd8Count = RESULT_NUMERIC) %>% 
  select(CFAR_PID, AGE_AT_RESULT_DATE, cd8Count) %>% 
  inner_join(tmp, by = c("CFAR_PID", "AGE_AT_RESULT_DATE")) %>% 
  mutate(cd4tcd8 = cd4Count / cd8Count, cluster = CFAR_PID) %>% 
  left_join(demo.d, by ="CFAR_PID") %>% 
  filter(sex == 'F') %>% 
  filter(CENSOR_TYPE == "END OF STUDY")
dat$cluster <- factor(dat$cluster, levels = unique(dat$cluster))
colnames(dat)[3:4] <- c("x", "y")
#Between-, within-cluster, and total Spearman rank correlations
#original scale
rankCorr::rankCorrCluster(dat$x, dat$y, dat$cluster, weights = "clusters")
#sqrt root transformation
dat$sqrtx <- sqrt(dat$x)
dat$sqrty <- sqrt(dat$y)
rankCorr::rankCorrCluster(dat$sqrtx, dat$sqrty, dat$cluster, weights = "clusters")
#log transformation
dat$logx <- log(dat$x)
dat$logy <- log(dat$y)
rankCorr::rankCorrCluster(dat$logx, dat$logy, dat$cluster, weights = "clusters")

#Between-, within-cluster, and total Pearson correlations
pi <- 1 / length(unique(dat$cluster))
k <- table(dat$cluster)
pij <- rep(pi / k, k)
#original scale
ymodr <- lmer(y ~ 1 + (1 | cluster), data = dat, REML = T)
ym <- predict(ymodr)
xmodr <- lmer(x ~ 1 + (1 | cluster), data = dat, REML = T)
xm <- predict(xmodr)
yr <- dat$y - ym
xr <- dat$x - xm
mean.xm <- sum(pij * xm); mean.ym <- sum(pij * ym)
b <- sum(pij * (xm - mean.xm) * (ym - mean.ym)) / sqrt(sum((xm - mean.xm)^2 * pij)* sum((ym - mean.ym)^2 * pij))
mean.xr <- sum(pij * xr); mean.yr <- sum(pij * yr)
w <- sum(pij * (xr - mean.xr) * (yr - mean.yr)) / sqrt(sum((xr - mean.xr)^2 * pij)* sum((yr - mean.yr)^2 * pij))
mean.x <- sum(pij * dat$x); mean.y <- sum(pij * dat$y)
t <- sum(pij * (dat$x - mean.x) * (dat$y - mean.y)) / sqrt(sum((dat$x - mean.x)^2 * pij)* sum((dat$y - mean.y)^2 * pij))
est1 <- c(b, w, t)
#sqrt root transformation
ymodr <- lmer(sqrty ~ 1 + (1 | cluster), data = dat, REML = T)
ym <- predict(ymodr)
xmodr <- lmer(sqrtx ~ 1 + (1 | cluster), data = dat, REML = T)
xm <- predict(xmodr)
yr <- dat$sqrty - ym
xr <- dat$sqrtx - xm
mean.xm <- sum(pij * xm); mean.ym <- sum(pij * ym)
b <- sum(pij * (xm - mean.xm) * (ym - mean.ym)) / sqrt(sum((xm - mean.xm)^2 * pij)* sum((ym - mean.ym)^2 * pij))
mean.xr <- sum(pij * xr); mean.yr <- sum(pij * yr)
w <- sum(pij * (xr - mean.xr) * (yr - mean.yr)) / sqrt(sum((xr - mean.xr)^2 * pij)* sum((yr - mean.yr)^2 * pij))
mean.x <- sum(pij * dat$sqrtx); mean.y <- sum(pij * dat$sqrty)
t <- sum(pij * (dat$sqrtx - mean.x) * (dat$sqrty - mean.y)) / sqrt(sum((dat$sqrtx - mean.x)^2 * pij)* sum((dat$sqrty - mean.y)^2 * pij))
est1s <- c(b, w, t)
#log transformation
ymodr <- lmer(logy ~ 1 + (1 | cluster), data = dat, REML = T)
ym <- predict(ymodr)
xmodr <- lmer(logx ~ 1 + (1 | cluster), data = dat, REML = T)
xm <- predict(xmodr)
yr <- dat$logy - ym
xr <- dat$logx - xm
mean.xm <- sum(pij * xm); mean.ym <- sum(pij * ym)
b <- sum(pij * (xm - mean.xm) * (ym - mean.ym)) / sqrt(sum((xm - mean.xm)^2 * pij)* sum((ym - mean.ym)^2 * pij))
mean.xr <- sum(pij * xr); mean.yr <- sum(pij * yr)
w <- sum(pij * (xr - mean.xr) * (yr - mean.yr)) / sqrt(sum((xr - mean.xr)^2 * pij)* sum((yr - mean.yr)^2 * pij))
mean.x <- sum(pij * dat$logx); mean.y <- sum(pij * dat$logy)
t <- sum(pij * (dat$logx - mean.x) * (dat$logy - mean.y)) / sqrt(sum((dat$logx - mean.x)^2 * pij)* sum((dat$logy - mean.y)^2 * pij))
est1l <- c(b, w, t)


#Patient Health Questionnaire-9 Score
load("hops-data.RData")
m <- hops_data %>% filter(sex.f == "Male") 
f <- hops_data %>% filter(sex.f == "Female")
lb <- c("depress", "hiv", "stigma",  "po_prop_days_nonadherent_12m", "po_prop_days_adherent_12m", "age_enrollment")
lnk <- c("logistic", "logistic", "logistic", "probit", "probit", "probit")
est <- list()
for(i in seq_along(lb)){
  df <- m[,c("couple_id", "cluster", lb[i])] %>% 
    inner_join(f[,c("couple_id", lb[i])], by = "couple_id") %>% 
    arrange(cluster)
  colnames(df) <- c("couple_id", "cluster", "x", "y")
  est[[lb[i]]] <- rankCorr::rankCorrCluster(df$x, df$y, df$cluster, link.x = lnk[i], link.y = lnk[i])
}


#Simulations
#functions to generate clustered data 
generate.data.1 <- function(seed, mu.x = 1, mu.y = -1, rho.B = 0.7, rho.W = 0.7, n.cluster = 10, 
                            cluster.limit = c(20, 30)){
  set.seed(seed)
  #generate cluster means
  sigma.xy <- matrix(c(1, rho.B, rho.B , 1), byrow = T, ncol = 2)
  xy.cluster <- mvtnorm::rmvnorm(n.cluster, c(mu.x, mu.y), sigma.xy)
  #cluster size 
  if(cluster.limit[1] == cluster.limit[2]){
    size.cluster <- rep(cluster.limit[1], n.cluster)
  }else size.cluster <- replicate(n.cluster, sample(cluster.limit[1]:cluster.limit[2], 1))
  #generate x and y for each cluster
  sgW <- matrix(c(1, rho.W, rho.W, 1), byrow = T, ncol = 2)
  dat <- list()
  cls <- pos.x <- pos.y <- list()
  for(i in 1:n.cluster){
    dat[[i]] <- mvtnorm::rmvnorm(size.cluster[i], xy.cluster[i,], sgW)
    cls[[i]] <- rep(i, each = size.cluster[i])
    pos.x[[i]] <- rep(xy.cluster[i,1], each = size.cluster[i]) 
    pos.y[[i]] <- rep(xy.cluster[i,2], each = size.cluster[i])
  }
  dat <- as.data.frame(do.call(rbind, dat))
  colnames(dat) <- c("x","y")
  dat$y <- exp(dat$y)
  dat$cluster <- as.factor(unlist(cls))
  dat$pos.x <- unlist(pos.x)
  dat$pos.y <- unlist(pos.y)
  return(dat)
} 

#exp cluster mean
generate.data.2 <- function(seed, mu.x = 1, mu.y = -1, rho.B = 0.7, rho.W = 0.7, n.cluster = 10, 
                            cluster.limit = c(20, 30)){
  set.seed(seed)
  #generate cluster means
  sigma.xy <- matrix(c(1, rho.B, rho.B , 1), byrow = T, ncol = 2)
  xy.cluster <- mvtnorm::rmvnorm(n.cluster, c(mu.x, mu.y), sigma.xy)
  #exponentiate 
  xy.cluster <- exp(xy.cluster)
  #cluster size 
  if(cluster.limit[1] == cluster.limit[2]){
    size.cluster <- rep(cluster.limit[1], n.cluster)
  }else size.cluster <- replicate(n.cluster, sample(cluster.limit[1]:cluster.limit[2], 1))
  #generate x and y for each cluster
  sgW <- matrix(c(1, rho.W, rho.W, 1), byrow = T, ncol = 2)
  dat <- list()
  cls <- pos.x <- pos.y <- list()
  for(i in 1:n.cluster){
    dat[[i]] <- mvtnorm::rmvnorm(size.cluster[i], xy.cluster[i,], sgW)
    cls[[i]] <- rep(i, each = size.cluster[i])
    pos.x[[i]] <- rep(xy.cluster[i,1], each = size.cluster[i]) 
    pos.y[[i]] <- rep(xy.cluster[i,2], each = size.cluster[i])
  }
  dat <- as.data.frame(do.call(rbind, dat))
  colnames(dat) <- c("x","y")
  dat$y <- exp(dat$y)
  dat$cluster <- as.factor(unlist(cls))
  dat$pos.x <- unlist(pos.x)
  dat$pos.y <- unlist(pos.y)
  return(dat)
} 

#Scenario I 
simulation1 <- function(iter, n.cls, n.cll, rho1, rho2, idx){
  sim_num <- 5
  ans <- list()
  trial <- (iter - 1) * sim_num
  cll <- paste(n.cll, collapse = "t")
  filename <- paste('output/S1-n',n.cls,'cls', cll, 'val', idx, '-',iter,'.RData', sep="")
  trial <- (iter - 1) * sim_num
  for(i in 1:sim_num){
      dat <- generate.data.1(seed = trial+i, mu.x = 1, mu.y = -1, rho.B = rho1, rho.W = rho2, n.cluster = n.cls, cluster.limit = n.cll)
      dat$y <- log(dat$y)
      est <- rankCorr::rankCorrCluster(df$x, df$y, df$cluster, link.x = lk.x, link.y = lk.y)
      if(n.cll[1] != n.cll[2]) est <- list(est, rankCorr::rankCorrCluster(df$x, df$y, df$cluster, link.x = lk.x, link.y = lk.y, weights = "clusters"))
      ans[[i]] <- est
      save(ans, file = filename)
  }
}

lk.x <- lk.y <- "probit"
rho <- cbind(c(0.8,0.8,0,0,0.8), c(0.7,0,0.8,0,-0.7))
#number of clusters = 100, cluster size = 10 
for(args in 1:200){
  n.cls <- 100; n.cll <- c(10, 10)
  for(i in 1:nrow(rho)){
    simulation1(args, n.cls, n.cll, rho[i,1], rho[i,2], i)
  }
}
#number of clusters = 100, cluster size = 20
for(args in 1:200){
  n.cls <- 100; n.cll <- c(20, 20)
  for(i in 1:nrow(rho)){
    simulation1(args, n.cls, n.cll, rho[i,1], rho[i,2], i)
  }
}
#number of clusters = 100, cluster size = 30
for(args in 1:200){
  n.cls <- 100; n.cll <- c(30, 30)
  for(i in 1:nrow(rho)){
    simulation1(args, n.cls, n.cll, rho[i,1], rho[i,2], i)
  }
}
#number of clusters = 100, cluster size = Unif(1,50) 
for(args in 1:200){
  n.cls <- 100; n.cll <- c(1, 50)
  for(i in 1:nrow(rho)){
    simulation1(args, n.cls, n.cll, rho[i,1], rho[i,2], i)
  }
}

#link function misspecification
lk.x <- lk.y <- "loglog"
for(args in 1:200){
  n.cls <- 100; n.cll <- c(20, 20)
  for(i in 1:nrow(rho)){
    simulation1(args, n.cls, n.cll, rho[i,1], rho[i,2], i)
  }
}
lk.x <- lk.y <- "cloglog"
for(args in 1:200){
  n.cls <- 100; n.cll <- c(20, 20)
  for(i in 1:nrow(rho)){
    simulation1(args, n.cls, n.cll, rho[i,1], rho[i,2], i)
  }
}

#Scenario II 
simulation2 <- function(iter, n.cls, n.cll, rho1, rho2, idx){
  sim_num <- 5
  ans <- list()
  trial <- (iter - 1) * sim_num
  cll <- paste(n.cll, collapse = "t")
  filename <- paste('output/S2-n',n.cls,'cls', cll, 'val', idx, '-',iter,'.RData', sep="")
  trial <- (iter - 1) * sim_num
  for(i in 1:sim_num){
    dat <- generate.data.1(seed = trial+i, mu.x = 1, mu.y = -1, rho.B = rho1, rho.W = rho2, n.cluster = n.cls, cluster.limit = n.cll)
    est <- rankCorr::rankCorrCluster(df$x, df$y, df$cluster, link.x = lk.x, link.y = lk.y)
    if(n.cll[1] != n.cll[2]) est <- list(est, rankCorr::rankCorrCluster(df$x, df$y, df$cluster, link.x = lk.x, link.y = lk.y, weights = "clusters"))
    ans[[i]] <- est
    save(ans, file = filename)
  }
}
lk.x <- lk.y <- "probit"
rho <- cbind(c(0.8,0.8,0,0,0.8), c(0.7,0,0.8,0,-0.7))
#number of clusters = 100, cluster size = 10 
for(args in 1:200){
  n.cls <- 100; n.cll <- c(10, 10)
  for(i in 1:nrow(rho)){
    simulation2(args, n.cls, n.cll, rho[i,1], rho[i,2], i)
  }
}
#number of clusters = 100, cluster size = 20
for(args in 1:200){
  n.cls <- 100; n.cll <- c(20, 20)
  for(i in 1:nrow(rho)){
    simulation2(args, n.cls, n.cll, rho[i,1], rho[i,2], i)
  }
}
#number of clusters = 100, cluster size = 30
for(args in 1:200){
  n.cls <- 100; n.cll <- c(30, 30)
  for(i in 1:nrow(rho)){
    simulation2(args, n.cls, n.cll, rho[i,1], rho[i,2], i)
  }
}
#number of clusters = 100, cluster size = Unif(1,50) 
for(args in 1:200){
  n.cls <- 100; n.cll <- c(1, 50)
  for(i in 1:nrow(rho)){
    simulation2(args, n.cls, n.cll, rho[i,1], rho[i,2], i)
  }
}

#Scenario III 
simulation3 <- function(iter, n.cls, n.cll, rho1, rho2, idx){
  sim_num <- 5
  ans <- list()
  trial <- (iter - 1) * sim_num
  cll <- paste(n.cll, collapse = "t")
  filename <- paste('output/S3-n',n.cls,'cls', cll, 'val', idx, '-',iter,'.RData', sep="")
  trial <- (iter - 1) * sim_num
  for(i in 1:sim_num){
    dat <- generate.data.2(seed = trial+i, mu.x = 1, mu.y = -1, rho.B = rho1, rho.W = rho2, n.cluster = n.cls, cluster.limit = n.cll)
    est <- rankCorr::rankCorrCluster(df$x, df$y, df$cluster, link.x = lk.x, link.y = lk.y)
    if(n.cll[1] != n.cll[2]) est <- list(est, rankCorr::rankCorrCluster(df$x, df$y, df$cluster, link.x = lk.x, link.y = lk.y, weights = "clusters"))
    ans[[i]] <- est
    save(ans, file = filename)
  }
}
lk.x <- lk.y <- "probit"
rho <- cbind(c(0.8,0.8,0,0,0.8), c(0.7,0,0.8,0,-0.7))
for(args in 1:200){
  n.cls <- 100; n.cll <- c(20, 20)
  for(i in 1:nrow(rho)){
    simulation3(args, n.cls, n.cll, rho[i,1], rho[i,2], i)
  }
}

#Ordinal data 
simulation_ordinal <- function(iter, n.cls, n.cll, rho1, rho2, idx, l){
  sim_num <- 5
  ans <- list()
  cll <- paste(n.cll, collapse = "t")
  filename <- paste('output/l',l,"n",n.cls,'cls', cll, 'val', idx, '-',iter,'.RData', sep="")
  trial <- (iter - 1) * sim_num
  for(i in 1:sim_num){
    dat <- generate.data.1(seed = trial+i, mu.x = 1, mu.y = -1, rho.B = rho1, rho.W = rho2, n.cluster = n.cls, cluster.limit = n.cll)
    dat$y <- log(dat$y)
    bi <- c(-Inf, qnorm(seq(1/l, 1-1/l, 1/l), mean = 1, sd = sqrt(2)), Inf)
    dat$x <- as.numeric(cut(dat$x, breaks = bi, labels = 1:l))  
    bi <- c(-Inf, qnorm(seq(1/l, 1-1/l, 1/l), mean = -1, sd = sqrt(2)), Inf)
    dat$y <- as.numeric(cut(dat$y, breaks = bi, labels = 1:l))  
    est <- rankCorr::rankCorrCluster(df$x, df$y, df$cluster, link.x = lk.x, link.y = lk.y)
    ans[[i]] <- est
    save(ans, file = filename)
  }
}
lk.x <- lk.y <- "logistic"
rho <- cbind(c(0.8,0.8,0,0,0.8), c(0.7,0,0.8,0,-0.7))
for(l in c(3, 5, 10)){
  for(args in 1:200){
    n.cls <- 100; n.cll <- c(20, 20)
    for(i in 1:nrow(rho)){
      simulation_ordinal(args, n.cls, n.cll, rho[i,1], rho[i,2], i, l)
    }
  }
}

#Negative rank ICC
generate.data.negICC <- function(seed, mu.x = 1, mu.y = -1, rho.B = 0.8, rho.W = 0.7, n.cluster = 10,
                                 icc.x = -0.5, icc.y = -0.5){
  set.seed(seed)
  vu <- 1; vr <- (1-icc.x)/(1+icc.x) * vu
  su <- matrix(c(1, rho.B, rho.B, 1) * vu, ncol = 2, byrow = T)
  sr <- matrix(c(1, rho.W, rho.W, 1) * vr, ncol = 2, byrow = T)
  mu.xy <- c(mu.x, mu.y)
  r.xy <- rep(0, 2)
  #generate x and y for each cluster
  dat <- list()
  cls <- list()
  for(i in 1:n.cluster){
    ui <- matrix(rep(MASS::mvrnorm(1, mu.xy, su), each=2), ncol=2)
    ri <- matrix(rep(MASS::mvrnorm(1, r.xy, sr), each = 2) * c(1, -1, 1, -1), ncol=2)
    dat[[i]] <-  ui + ri
    cls[[i]] <- rep(i, each = 2)
  }
  dat <- as.data.frame(do.call(rbind, dat))
  colnames(dat) <- c("x","y")
  dat$y <- exp(dat$y)
  dat$cluster <- as.factor(unlist(cls))
  return(dat)
}

simulation_negICC <- function(iter, n.cls, n.cll, rho1, rho2, idx){
  sim_num <- 5
  ans <- list()
  trial <- (iter - 1) * sim_num
  cll <- paste(n.cll, collapse = "t")
  filename <- paste('output/negICC-n',n.cls,'cls', cll, 'val', idx, '-',iter,'.RData', sep="")
  trial <- (iter - 1) * sim_num
  for(i in 1:sim_num){
    dat <- generate.data.negICC(seed = trial+i, mu.x = 1, mu.y = -1, rho.B = rho1, rho.W = rho2, n.cluster = n.cls)
    est <- rankCorr::rankCorrCluster(df$x, df$y, df$cluster, link.x = lk.x, link.y = lk.y)
    if(n.cll[1] != n.cll[2]) est <- list(est, rankCorr::rankCorrCluster(df$x, df$y, df$cluster, link.x = lk.x, link.y = lk.y, weights = "clusters"))
    ans[[i]] <- est
    save(ans, file = filename)
  }
}
lk.x <- lk.y <- "probit"
rho <- cbind(c(0.8, 0.8, 0, 0), c(-0.7, 0.7, 0.8, 0))
for(args in 1:200){
  n.cls <- 100; n.cll <- c(2, 2)
  for(i in 1:nrow(rho)){
    simulation_negICC(args, n.cls, n.cll, rho[i,1], rho[i,2], i)
  }
}
