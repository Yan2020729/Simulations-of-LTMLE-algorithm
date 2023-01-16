rm(list=ls())

#------- small functions ----------------
.bound <- function(x, bounds){
  x[x<min(bounds)] <- min(bounds); x[x>max(bounds)] <- max(bounds)
  return(x)
}
rowSums<-function (x, na.rm = FALSE, dims = 1) {
  
  if (is.data.frame(x)) x <- as.matrix(x)
  if (length(dn <- dim(x)) < 2) return(x)
  if (dims < 1 || dims > length(dn) - 1) stop("invalid dims")
  p <- prod(dn[-(1:dims)]); dn <- dn[1:dims]
  
  z <- if (is.complex(x)) 
    .Internal(rowSums(Re(x), prod(dn), p, na.rm)) + (0+1i) * .Internal(rowSums(Im(x), prod(dn), p, na.rm))
  else .Internal(rowSums(x, prod(dn), p, na.rm))
  if (length(dn) > 1) {dim(z) <- dn; dimnames(z) <- dimnames(x)[1:dims]
  } else names(z) <- dimnames(x)[[1]]
  z
}  

#---------------------- estimate pA and pC-----------------------
pA = function(r, long_dat) {
  p_A = NULL
  for (t in 1:2) {
    subA = long_dat[which(long_dat$t==t & long_dat$A_minus ==0),] # since A1 is before C1, so remove C==0 
    pA = glm.fit(y = subA$A, x = cbind(1, subA$w1, subA$w2, subA$L_minus), family=binomial())
    fitA <- plogis( as.matrix(cbind(1, long_dat$w1, long_dat$w2, long_dat$L_minus)) %*% coef(pA) )*a_sim2[r,t]+
      (1 - plogis( as.matrix(cbind(1, long_dat$w1, long_dat$w2, long_dat$L_minus)) %*% coef(pA)) )*(1-a_sim2[r,t])
    p_A[[t]] = fitA
  }
  return(list(pA_t1 = p_A[[1]], pA_t2 = p_A[[2]]))
}

#---------------------------------- simulation ----------------------------------
# O=(w1, w2, L0, A1, L1, C1, A2, L2, C2, Y)
source(file = "sim_data.R") 
iter = 500
IF_sim = sum_fit = NULL
ests_msm = sds = sds_c = beta_coverage = beta_c_coverage = matrix(NA, nrow = iter, ncol = 4)
seeds = sample(1:5000, iter)

a_sim2 = matrix(c(1,1,0,1,0,0), nrow = 3, byrow = TRUE)

for (j in 1:iter) {
  
  N= 5000; study_n = 100
  datl1 = datal(seed=seeds[j], N=N, study_n=study_n, study = "itt")
  id = datl1$ind_id
  long_dat = longData(data = datl1, study = "itt")
  
  pA_allt = NULL; 
  for (regim in 1:3 ) { pA_allt[[regim]] = pA(r=regim, long_dat = long_dat)}
  
  p_C = NULL
  for (t in 1:2) {
    pC <- fitC <- NULL
    
    subC = long_dat[which(long_dat$t==t ),]
    pC <- glm.fit(y = subC$C, x = cbind(1, subC$w1, subC$w2, subC$L, subC$A), family=binomial())
    fitC <- 1 - plogis(as.matrix(cbind(1, long_dat$w1, long_dat$w2, long_dat$L, long_dat$A)) %*% coef(pC))
    p_C[[t]] = fitC
  }
  
  
  #----------------------------------- cumulative weights ------------------------------------
  weights_sim_allr <- NULL
  
  for (t in 1:3) {
    weights_sim_allr[[t]] <- NULL; ws = NULL
    if (t ==1) {
      for (r in 1:3) {
        w = 1/(sapply(pA_allt[[r]]$pA_t1, as.numeric))
        ws = c(ws, w)
      }
    } else if (t ==3) {
      for (r in 1:3) {
        if (r==1) {w = 1/( exp(rowSums(log(pA_allt[[r]]$pA_t1)+log(p_C[[1]])+log(p_C[[2]] ))) )
        } else w = 1/( exp(rowSums(log(pA_allt[[r]]$pA_t1)+log(pA_allt[[r]]$pA_t2)+log(p_C[[1]])+log(p_C[[2]] ))) )
        ws = c(ws, w)
      }
    } else { # for t =2
      for (r in 1:3) {
        if (r==1) {w = 1/( exp(rowSums(log(pA_allt[[r]]$pA_t1)+log(p_C[[1]]))) )
        } else w = 1/( exp(rowSums(log(pA_allt[[r]]$pA_t1)+log(pA_allt[[r]]$pA_t2)+log(p_C[[1]]))) )
        ws = c(ws, w)
      }
    }
    weights_sim_allr[[t]] = ws
  }
  
  #######################################################################################################
  
  Ddd <- NULL;
  
  #------------------------------------------ t = 3 --------------------------------------------
  t = 3
  Q_l <- NULL; long_n = dim(long_dat)[1]
  
  for (r in 1:3) {
    
    Q_l[[r]] = rep(0, long_n)
    col_Q = c(5:8,12) # "w1" ,"w2" , "L" , "L_minus" , "cumA"
    
    subQ <- which(long_dat$t==t-1 & long_dat$C==0)
    model = glm.fit(y = long_dat[subQ,]$Y, x = cbind(1, long_dat[subQ, col_Q], long_dat$w2[subQ]* long_dat$cumA[subQ]) ) 
    Q_l[[r]] = as.matrix(cbind(1, long_dat[, 5:8], cumA = rep(sum(a_sim2[r,1:(t-1)]), long_n), long_dat$w2*rep(sum(a_sim2[r,1:(t-1)]), long_n))) %*% coef(model)
    
  }
  
  Ql_3r = unlist(Q_l) 
  
  a = min(c(long_dat$Y, Ql_3r), na.rm = TRUE); b = max(c(long_dat$Y, Ql_3r), na.rm = TRUE)
  
  Ql_3r_scaled = (Ql_3r-a)/(b-a)
  bound=0.001
  Ql_3r_scaled = .bound(Ql_3r_scaled, c(bound, 1-bound))
  Y_scaled = (long_dat$Y-a)/(b-a)
  Y_scaled = .bound(Y_scaled, c(bound, 1-bound))
  # ------------------------ update Q3; working model (w2 + cumA + w2:cumA)
  doffset = qlogis(Ql_3r_scaled)
  coe_cum3 = rep(2:0, each = long_n)
  
  stacklong_dat <- data.frame(id = rep(long_dat$id, 3), study_id= rep(long_dat$study_id, 3), t = rep(long_dat$t, 3), C = rep(long_dat$C, 3), 
                              r = rep(1:3, each = long_n), Y_scaled = rep(Y_scaled, 3), w2 = rep(long_dat$w2, 3), 
                              cumA = rep(long_dat$cumA, 3), cumA_minus = rep(long_dat$cumA_minus, 3), cumA_pl = rep(long_dat$cumA_pl, 3)) 
  stacklong_n = dim(stacklong_dat)[1] # 
  
  subjs = NULL
  for ( r in 1:3) {
    rsub = (1 + (r-1) * long_n) : (long_n +  (r-1) * long_n) 
    subs = length(rsub)*(r-1)+which(stacklong_dat[rsub,]$t==t-1 & stacklong_dat[rsub,]$C==0 & stacklong_dat[rsub,]$cumA == t-r)
    subjs = c(subjs, subs)
  }
    y = stacklong_dat[subjs,]$Y_scaled
    
    wc3 = matrix(cbind(1, stacklong_dat$w2, coe_cum3, stacklong_dat$w2*coe_cum3), nrow = stacklong_n) 
    modelu <- glm(y ~ -1 + wc3[subjs,], offset = doffset[subjs], weights = weights_sim_allr[[t]][subjs], family = quasibinomial()) 
    
    estQ_sim = rep(NA, stacklong_n)
    estQ_sim[which(stacklong_dat$t==t-1)] <- plogis( doffset[which(stacklong_dat$t==t-1)] +
                                                       wc3[which(stacklong_dat$t==t-1),] %*% coef(modelu) )
    Dds= rep(0, stacklong_n)
    Dds[subjs] = (stacklong_dat[subjs,]$Y_scaled - estQ_sim[subjs]) * weights_sim_allr[[t]][subjs] 
    mean(Dds, na.rm = TRUE)
  
  #------------------------------------------ t = 2 --------------------------------------------
  t = 2
  Q_l <- NULL; 
  
  for (r in 1:3) {
    rsub = (1 + (r-1) * long_n) : (long_n +  (r-1) * long_n)
    
    Q_l[[r]] = rep(0, long_n)
    col_Q = c(5:8, 14)
    
    subQ <- which(long_dat$t==t-1 & long_dat$C==0)
    model = glm.fit(y = estQ_sim[rsub][which(long_dat$t==t)], x = cbind(1, long_dat[subQ, col_Q]), family=quasibinomial() )
    Q_l[[r]] = plogis(as.matrix(cbind(1, long_dat[, 5:8], rep(sum(a_sim2[r,1:(t)]), long_n))) %*% coef(model))
    
  }
  
  Ql_2r = unlist(Q_l) 
  subjs = subys = NULL
  for ( r in 1:3) {
    rsub = (1 + (r-1) * long_n) : (long_n +  (r-1) * long_n) 
    subs = length(rsub)*(r-1)+which(stacklong_dat[rsub,]$t==t-1 & stacklong_dat[rsub,]$C==0 & stacklong_dat[rsub,]$cumA_pl == sum(a_sim2[r,1:(t)]))
    subjs = c(subjs, subs)
    
    suby = length(rsub)*(r-1)+which(stacklong_dat[rsub,]$t==t & stacklong_dat[rsub,]$cumA == sum(a_sim2[r,1:(t)]))
    subys = c(subys, suby)
  }
  
  y = estQ_sim[subys]
  doffset = qlogis(Ql_2r)
  
  coe_cum2 = rep(c(2:0), each = long_n)
  wc2 = matrix(cbind(1, stacklong_dat$w2, coe_cum2, stacklong_dat$w2*coe_cum2), nrow = stacklong_n) 
  modelu <- glm(y ~ -1 + wc2[subjs,], offset = doffset[subjs], weights =  weights_sim_allr[[t]][subjs], family = quasibinomial()) 
  estQ_sim[which(stacklong_dat$t==t-1)] <- plogis( doffset[which(stacklong_dat$t==t-1)] +
                                                     wc2[which(stacklong_dat$t==t-1),]%*% coef(modelu) )
 
  mean(estQ_sim[which(stacklong_dat$t==t-1)])*(b-a)+a
  Dds[subjs] = (y - estQ_sim[subjs]) * weights_sim_allr[[t]][subjs] 
 
  #------------------------------------------ t = 1 --------------------------------------------
  t = 1
  Q_l <- NULL; 
  
  for (r in 1:3) {
    rsub = (1 + (r-1) * long_n) : (long_n +  (r-1) * long_n)
    
    Q_l[[r]] = rep(0, long_n)
    col_Q = c(5:7, 14) #  "w1", "w2" , "L" , "cumA+", remove  L-
    
    subQ <- which(long_dat$t==t-1 & long_dat$C==0)
    model = glm.fit(y = estQ_sim[rsub][which(long_dat$t==t)], x = cbind(1, long_dat[subQ, col_Q]), family=quasibinomial() ) 
    
    Q_l[[r]] = plogis(as.matrix(cbind(1, long_dat[, 5:7], rep(a_sim2[r,(t)],long_n))) %*% coef(model))
    
  }
  
  Ql_1r = unlist(Q_l) 
  
  subjs = subys = NULL
  for ( r in 1:3) {
    rsub = (1 + (r-1) * long_n) : (long_n +  (r-1) * long_n) 
    subs = length(rsub)*(r-1)+which(stacklong_dat[rsub,]$t==t-1 & stacklong_dat[rsub,]$C==0 & stacklong_dat[rsub,]$cumA_pl == sum(a_sim2[r,1:(t)]))
    subjs = c(subjs, subs)
    
    suby = length(rsub)*(r-1)+which(stacklong_dat[rsub,]$t==t & stacklong_dat[rsub,]$cumA == sum(a_sim2[r,1:(t)]))
    subys = c(subys, suby)
  }
  
  y = estQ_sim[subys]
  doffset = qlogis(Ql_1r)
  
  coe_cum1 = rep(c(2:0), each = long_n)
  wc1 = matrix(cbind(1, stacklong_dat$w2, coe_cum1, stacklong_dat$w2*coe_cum1), nrow = stacklong_n) 
  modelu <- glm(y ~ -1 + wc1[subjs,], offset = doffset[subjs], weights = weights_sim_allr[[t]][subjs], family = quasibinomial()) 
  estQ_sim[which(stacklong_dat$t==t-1)] <- plogis( doffset[which(stacklong_dat$t==t-1)] +
                                                     wc1[which(stacklong_dat$t==t-1),]%*% coef(modelu) )
  
  mean(estQ_sim[which(stacklong_dat$t==t-1)])*(b-a)+a
  Dds[subjs] = (y - estQ_sim[subjs]) * weights_sim_allr[[t]][subjs] 
  
  #------------------------------ pooled ltmle for MSM ------------------------------------
  totalq = estQ_sim[which(stacklong_dat$t==0)]
  length(totalq); mean(totalq)*(b-a)+a
  
  msmA = coe_cum3[which(stacklong_dat$t==0)] 
  covmat_sim = cbind(intercep = 1, w2=stacklong_dat[which(stacklong_dat$t==0),]$w2, cumA=msmA,
                     w2_cumA = stacklong_dat[which(stacklong_dat$t==0),]$w2 * msmA)
   totalq_orig = totalq*(b-a)+a
  fit_sim = glm.fit(y = totalq_orig, x = covmat_sim)
  
  ests_msm[j,] = as.vector(coef(fit_sim))
  sum_fit[[j]] = summary.glm(fit_sim)
  
  #---------------------- sum EIF on unvaried DD
  id = datl1$ind_id
  ddss = NULL
  for ( r in 1:3) {
    rsub = (1 + (r-1) * long_n) : (long_n +  (r-1) * long_n) 
    ds = sapply(1: length(id), function(x) sum(Dds[rsub][which(stacklong_dat[rsub,]$id == id[x])]))
    ddss = c(ddss, ds)
  }
 
  #------------------- IF(a[r,])  on  DD
  
  fitt_sim = glm.fit(y = totalq, x = covmat_sim)
  diff1_sim = totalq -predict.glm(fitt_sim)
  diff2_sim = ddss 
  
  sumDif = covmat_sim*(diff1_sim+diff2_sim)
  sumD = sumDif[1:N,]+sumDif[(N+1):(2*N),]+sumDif[(2*N+1):(3*N),]
  
   # ---- normlization matrix (sum_r vv'/n)^(-1)
  covmat_p1 = t(covmat_sim[1:N,]) %*% covmat_sim[1:N,]
  covmat_p2 = t(covmat_sim[(N+1):(2*N),]) %*% covmat_sim[(N+1):(2*N),]
  covmat_p3 = t(covmat_sim[(2*N+1):(3*N),]) %*% covmat_sim[(2*N+1):(3*N),]
  sumCOv = covmat_p1+covmat_p2+covmat_p3
  norm_sim <- solve(sumCOv/length(id))
  
  IFbeta_sim <- (sumD) %*% norm_sim
  IF_sim[[j]] = IFbeta_sim
  
  # ------ pooled sd 
  sd_coef_sim <- sqrt(colSums(IFbeta_sim^2)/(length(id))^2)*(b-a)
  sds[j,] = as.vector(sd_coef_sim)
  
  ciL_coef <- ests_msm[j,]-1.96*sd_coef_sim
  ciH_coef <- ests_msm[j,]+1.96*sd_coef_sim
  beta_coverage[j,] <- as.vector(sapply(1:4, function(x) as.numeric(tru_beta[x]>=ciL_coef[x] & tru_beta[x]<=ciH_coef[x])))
  
  
  # ------ sd with clusters 
  clustDD_sim <- data.frame(study_id=datl1$study_id, IFbeta_sim)
  split_DD_sim <- split(clustDD_sim, clustDD_sim$study_id)
  clusterN_sim <- unique(datl1$study_id)# 1,...,50
  
  sd_c_coef_sim <- sqrt(rowSums( sapply(1:length(clusterN_sim), function(x) colSums(split_DD_sim[[x]][,-1])^2) )/(length(id))^2)*(b-a)
  sds_c[j,] = as.vector(sd_c_coef_sim)
  
  ciL_c_coef <- ests_msm[j,]-1.96*sd_c_coef_sim
  ciH_c_coef <- ests_msm[j,]+1.96*sd_c_coef_sim
  beta_c_coverage[j,] <- as.vector(sapply(1:4, function(x) as.numeric(tru_beta[x]>=ciL_c_coef[x] & tru_beta[x]<=ciH_c_coef[x])))
 
}

est_coefs = round(colMeans(ests_msm),4)
sd_mc = round(sapply(1:4, function(x) sd(ests_msm[,x])),4)
sd_unD = round(colMeans(sds),4)
sd_c_unD = round(colMeans(sds_c),4)
beta_cov = colSums(beta_coverage)/iter
beta_c_cov = colSums(beta_c_coverage)/iter

# source(file = "sim_true_coef.R")
tru_beta = colMeans(true_coeff_itt)
resu = round(rbind(tru_beta, est_coefs, sd_mc, sd_unD, sd_c_unD, beta_cov, beta_c_cov),2) 
#        intercep  w2  cumA w2_cumA
# tru_beta   0.06 1.95 0.59 0.20
# est_coefs  0.07 1.95 0.59 0.20
# sd_mc      0.17 0.26 0.08 0.16
# sd_unD     0.14 0.30 0.08 0.18
# sd_c_unD   0.17 0.31 0.08 0.17
# beta_cov   0.85 0.97 0.97 0.97
# beta_c_cov 0.94 0.97 0.96 0.97
