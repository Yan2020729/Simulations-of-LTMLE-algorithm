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



#---------------------------------- simulation ----------------------------------
# j=1000
# O=(w1, w2, L0, A1, L1, C1, A2, L2, C2, Y)

# datal <- function(seed, N, study_n)  {
#   
#   set.seed(seed)
#   # Scenario 1: D is independent
#   clusterN = N/study_n
#   study_id = c(rep(1:clusterN, each=study_n))
#   ind_id = c(1:N)
#   
#   # generate study-level covariate V
#   s_bar = rnorm(clusterN, mean=0.5, sd=1)
#   s = rep(s_bar, each=study_n)
#   
#   # generate 9000 individual covariates W
#   w1<-rnorm(N, mean=0.5*s, sd=0.7)
#   w2<-rbinom(N, 1, plogis(0.5+0.8*s))
#   # w2<-rbinom(N, 1, 0.5)
#   L0 = rnorm(N, mean=0.2+0.2*w1+0.3*w2, sd=0.7)
#   A0 = C0 = rep(0, N)
#   
#   # t1 = rep(1, N)
#   p_A1 = plogis(0.7*w1+0.8*w2+0.9*L0)
#   A1 = rbinom(N, 1, p_A1)
#   
#   L1 = rnorm(N, mean=0.2*w1+0.3*w2+0.5*L0 + 0.5*A1, sd=1)
#   p_C1 = plogis(-2.5+0.2*w1+0.4*w2+0.3*L1+0.2*A1)
#   C1 = rbinom(N, 1,  p_C1)
#   #A1
#   p_A2 = plogis(0.5*w1+0.3*w2+0.7*L1)
#   A2 = rbinom(N, 1, p_A2) 
#   # A2 = ifelse( A1==1, 1, A20) # itt: if A1=1, A2=1
#   
#   L20 = rnorm(N, mean=0.2*w1+0.3*w2+0.7*L1+0.4*A2, sd=0.8)
#   
#   p_C2 = plogis(-3+0.2*w1+0.4*w2+0.3*L20+0.4*A2)# ; mean(p_C2, na.rm = TRUE)
#   C2 = ifelse( C1==0,  rbinom(N, 1, p_C2), NA)
#   total_t = ifelse(C1 ==0, 2, 1)
#   
#   cumA = ifelse( C1==0, A1+A2, A1)
#   
#   err <- rnorm(N,0,1)
#   Y = ifelse(C2==0, 0.2*s+0.6*w1+1*w2+0.4*L1+0.5*L20+0.3*cumA+0.2*cumA*w2+ err, NA) 
#   
#   A2 = ifelse( C1==0, A2, NA)
#   L2 = ifelse( C1==0, L20, NA)
#   
#   datal = data.frame(ind_id, study_id, total_t, s, w1, w2, L0, A0, C0, A1, L1, C1, A2, cumA, L2, C2, Y)
#   
#   return(datal)
# }


# longData = function(data) {
#   
#   t = id = study_id = Y = s = w1 = w2 = L = L_minus = A = A_minus = A_pl = cumA = cumA_minus = cumA_pl = C = Y = NULL
#   
#   for (i in 1:dim(datl1)[1]) { 
#     ind_t=0:datl1$total_t[i]; t = c(t, ind_t)
#     id = c(id, rep(datl1$ind_id[i], length(ind_t)))
#     study_id = c(study_id, rep(datl1$study_id[i], length(ind_t)))
#     s = c(s, rep(datl1$s[i], length(ind_t)))
#     w1 = c(w1, rep(datl1$w1[i], length(ind_t)))
#     w2 = c(w2, rep(datl1$w2[i], length(ind_t)))
#     
#     if (length(ind_t)==2) {
#       ll = c(datl1$L0[i], datl1$L1[i]) 
#       llmin = c(datl1$L0[i], datl1$L0[i])
#       aa = c(datl1$A0[i], datl1$A1[i])
#       aamin = c(datl1$A0[i], datl1$A0[i])
#       aapl = c(datl1$A1[i], datl1$A1[i])
#       cc = c(datl1$C0[i], datl1$C1[i])
#       cuma = c(0, ifelse(datl1$A1[i]==1, 1, 0))
#       cumamin = c(0, 0)
#       cumapl = c(ifelse(datl1$A1[i]==1, 1, 0), ifelse(datl1$A1[i]==1, 1, 0))
#       y = c(NA, NA)
#       
#     } else {
#       ll = c(datl1$L0[i], datl1$L1[i], datl1$L2[i])
#       llmin = c(datl1$L0[i], datl1$L0[i], datl1$L1[i])
#       aa = c(datl1$A0[i], datl1$A1[i], datl1$A2[i])
#       aamin = c(datl1$A0[i], datl1$A0[i], datl1$A1[i])
#       aapl = c(datl1$A1[i], datl1$A2[i], datl1$A2[i])
#       cc = c(datl1$C0[i], datl1$C1[i], datl1$C2[i])
#       
#       if (sum(aa) ==1) {
#         if (datl1$A1[i]==1) {cuma = c(0, 1, 1); cumapl = c(1,1,1)
#         } else {cuma = c(0, 0, 1); cumapl = c(0,1,1)}
#       } else if (sum(aa) ==0) {cuma =cumapl =  c(0, 0, 0)
#       } else {cuma = c(0, 1, 2); cumapl = c(1,2,2) }
#       
#       cumamin = c(0, 0, ifelse(datl1$A1[i]==1, 1, 0))
#       y = c(NA, NA, datl1$Y[i])
#     }
#     L = c(L, ll)
#     L_minus = c(L_minus, llmin)
#     A = c(A, aa)
#     A_minus = c(A_minus, aamin)
#     A_pl = c(A_pl, aapl)
#     C = c(C, cc)
#     cumA = c(cumA, cuma)
#     cumA_minus = c(cumA_minus, cumamin)
#     cumA_pl = c(cumA_pl, cumapl)
#     Y= c(Y, y)
#   }
#   
#   long_dat = data.frame(id, study_id, t, s, w1, w2, L, L_minus, A, A_minus, A_pl, cumA, cumA_minus, cumA_pl, C, Y)
#   
#   
#   return(long_dat)
#   
# }

#---------------------- estimate pA and pC-----------------------
pA = function(r, long_dat) {
  p_A = NULL
  for (t in 1:2) {
    subA = long_dat[which(long_dat$t==t),] # remove A-==0 # diff with itt
    pA = glm.fit(y = subA$A, x = cbind(1, subA$w1, subA$w2, subA$L_minus), family=binomial())
    fitA <- plogis( as.matrix(cbind(1, long_dat$w1, long_dat$w2, long_dat$L_minus)) %*% coef(pA) )*a_asTreat[r,t]+
      (1 - plogis( as.matrix(cbind(1, long_dat$w1, long_dat$w2, long_dat$L_minus)) %*% coef(pA)) )*(1-a_asTreat[r,t])
    p_A[[t]] = fitA
  }
  return(list(pA_t1 = p_A[[1]], pA_t2 = p_A[[2]]))
}



# #-------  long format data ------------------
source(file = "./test/sim_data.R") 
# 
# datl1 = datal(seed=1234, N=5000, study_n=100, study = "at")
# head(datl1)
# View(datl1)
# 
# long_dat = longData(data = datl1, study = "at")
# head(long_dat)
# View(long_dat)

#------------------------ regimes--------------------------
a_asTreat = matrix(c(1,1,0,1,0,0,1,0), nrow = 4, byrow = TRUE)
#        [,1] [,2]
# [1,]    1    1
# [2,]    0    1
# [3,]    0    0
# [4,]    1    0


############################################################################################
iter = 500
IF_sim = sum_fit = NULL
ests_msm = sds = sds_c = beta_coverage = beta_c_coverage = matrix(NA, nrow = iter, ncol = 4)
# seeds = sample(1:5000, 1000)
load(file = "./Rda/seeds1000")
# source("./test/true_coef")

for (j in 1:iter) {
  
  N= 5000; study_n =100
  datl1 = datal(seed=seeds[200+j], N=N, study_n=study_n, study = "at")
  
  id = datl1$ind_id
  long_dat = longData(data = datl1, study = "at")
  
  
  pA_allt = NULL; 
  for (regim in 1:4 ) { pA_allt[[regim]] = pA(r=regim, long_dat = long_dat)} # diff with itt, regim =1:4
  
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
      for (r in 1:4) {# diff with itt
        w = 1/(sapply(pA_allt[[r]]$pA_t1, as.numeric))
        ws = c(ws, w)
      }
    } else if (t ==3) {
      for (r in 1:4) {# diff with itt
        w = 1/( exp(rowSums(log(pA_allt[[r]]$pA_t1)+log(pA_allt[[r]]$pA_t2)+log(p_C[[1]])+log(p_C[[2]] ))) )# diff with itt
        ws = c(ws, w)
      }
    } else { # for t =2
      for (r in 1:4) {# diff with itt
        w = 1/( exp(rowSums(log(pA_allt[[r]]$pA_t1)+log(pA_allt[[r]]$pA_t2)+log(p_C[[1]]))) )# diff with itt
        ws = c(ws, w)
      }
    }
    
    # ws = .bound(ws, quantile(ws, c(0.01, 0.99)))
    weights_sim_allr[[t]] = ws
  }
  
  
  
  #######################################################################################################
  
  Ddd <- NULL;
  
  #------------------------------------------ t = 3 --------------------------------------------
  t = 3
  Q_l <- NULL; long_n = dim(long_dat)[1]
  
  for (r in 1:4) {# diff with itt, r =1:4
    
    Q_l[[r]] = rep(0, long_n)
    col_Q = c(5:8,9,12) # "w1" ,"w2" , "L" , "L_minus" ,"A", "cumA" # diff with itt (add "A")
    
    subQ <- which(long_dat$t==t-1 & long_dat$C==0)
    model = glm.fit(y = long_dat[subQ,]$Y, x = cbind(1, long_dat[subQ, col_Q], long_dat$w2[subQ]* long_dat$cumA[subQ]) ) 
    
    Q_l[[r]] = as.matrix(cbind(1, long_dat[, 5:8], A= rep(a_asTreat[r,(t-1)], long_n), cumA = rep(sum(a_asTreat[r,1:(t-1)]), long_n), 
                               long_dat$w2*rep(sum(a_asTreat[r,1:(t-1)]), long_n))) %*% coef(model)
    
  }
  
  Ql_3r = unlist(Q_l) 
  
  # mean(Ql_3r[which(stacklong_dat$t==t-1)])
  # mean(Ql_3r[which(stacklong_dat$t==t-1)][1:(12792/3)])
  # mean(Ql_3r[which(stacklong_dat$t==t-1)][(12792/3+1):(2*12792/3)])
  # mean(Ql_3r[which(stacklong_dat$t==t-1)][(2*12792/3+1):(12792)])
  
  
  a = min(c(long_dat$Y, Ql_3r), na.rm = TRUE); b = max(c(long_dat$Y, Ql_3r), na.rm = TRUE)
  
  Ql_3r_scaled = (Ql_3r-a)/(b-a)
  bound=0.001
  Ql_3r_scaled = .bound(Ql_3r_scaled, c(bound, 1-bound))
  Y_scaled = (long_dat$Y-a)/(b-a)
  Y_scaled = .bound(Y_scaled, c(bound, 1-bound))
  
  
  # ------------------------ update Q3; working model (w2 + cumA + w2:cumA)
  doffset = qlogis(Ql_3r_scaled)
  coe_cum3 = rep(c(2,1,0,1), each = long_n) # diff with itt
  
  stacklong_dat <- data.frame(id = rep(long_dat$id, 4), study_id= rep(long_dat$study_id, 4), t = rep(long_dat$t, 4), C = rep(long_dat$C, 4), # diff with itt
                              r = rep(1:4, each = long_n), Y_scaled = rep(Y_scaled, 4), w2 = rep(long_dat$w2, 4), # diff with itt
                              A= rep(long_dat$A, 4),  A_pl= rep(long_dat$A_pl, 4), # diff with itt
                              cumA = rep(long_dat$cumA, 4), cumA_minus = rep(long_dat$cumA_minus, 4), cumA_pl = rep(long_dat$cumA_pl, 4)) # diff with itt
  stacklong_n = dim(stacklong_dat)[1] # 
  
  subjs = NULL
  for ( r in 1:4) { # diff with itt
    rsub = (1 + (r-1) * long_n) : (long_n +  (r-1) * long_n) 
    subs = length(rsub)*(r-1)+which(stacklong_dat[rsub,]$t==t-1 & stacklong_dat[rsub,]$C==0 & stacklong_dat[rsub,]$A==a_asTreat[r,(t-1)] & # diff with itt
                                      stacklong_dat[rsub,]$cumA == sum(a_asTreat[r,1:(t-1)]))
    subjs = c(subjs, subs)
  }
  
  
  y = stacklong_dat[subjs,]$Y_scaled
  # y = stacklong_dat[subjs,]$Y_scaled
  
  wc3 = matrix(cbind(1, stacklong_dat$w2, coe_cum3, stacklong_dat$w2*coe_cum3), nrow = stacklong_n) 
  modelu <- glm(y ~ -1 + wc3[subjs,], offset = doffset[subjs], weights = weights_sim_allr[[t]][subjs], family = quasibinomial()) 
  
  estQ_sim = rep(NA, stacklong_n)
  estQ_sim[which(stacklong_dat$t==t-1)] <- plogis( doffset[which(stacklong_dat$t==t-1)] +
                                                     wc3[which(stacklong_dat$t==t-1),] %*% coef(modelu) )
  
  # mean(estQ_sim[which(stacklong_dat$t==t-1)])*(b-a)+a
  # mean(estQ_sim[which(stacklong_dat$t==t-1)][1:(12792/3)])*(b-a)+a
  # mean(estQ_sim[which(stacklong_dat$t==t-1)][(12792/3+1):(2*12792/3)])*(b-a)+a
  # mean(estQ_sim[which(stacklong_dat$t==t-1)][(2*12792/31):(12792)])*(b-a)+a
  
  
  # ---------------------  compute EIF
  Dds= rep(0, stacklong_n)
  Dds[subjs] = (stacklong_dat[subjs,]$Y_scaled - estQ_sim[subjs]) * weights_sim_allr[[t]][subjs] 
  mean(Dds, na.rm = TRUE)
  
  
  
  #------------------------------------------ t = 2 --------------------------------------------
  t = 2
  Q_l <- NULL; 
  
  for (r in 1:4) { # diff with itt
    rsub = (1 + (r-1) * long_n) : (long_n +  (r-1) * long_n)
    
    Q_l[[r]] = rep(0, long_n)
    # col_Q = c(5:8, 11, 12) # # "w1" ,"w2" , "L" , "L_minus" , "A+", "cumA+" 
    col_Q = c(5:8, 11,14) # diff with itt, add "A+"
    
    subQ <- which(long_dat$t==t-1 & long_dat$C==0)
    model = glm.fit(y = estQ_sim[rsub][which(long_dat$t==t)], x = cbind(1, long_dat[subQ, col_Q]), family=quasibinomial() )
    Q_l[[r]] = plogis(as.matrix(cbind(1, long_dat[, 5:8], rep(a_asTreat[r,t], long_n), rep(sum(a_asTreat[r,1:(t)]), long_n))) %*% coef(model)) # diff with itt
    
  }
  
  Ql_2r = unlist(Q_l) 
  
  # mean(Ql_2r[which(stacklong_dat$t==t-1)])*(b-a)+a
  # mean(Ql_2r[which(stacklong_dat$t==t-1)][1:(12792/3)])*(b-a)+a
  # mean(Ql_2r[which(stacklong_dat$t==t-1)][(12792/3+1):(2*12792/3)])*(b-a)+a
  # mean(Ql_2r[which(stacklong_dat$t==t-1)][(2*12792/3+1):(12792)])*(b-a)+a
  
  # ------------------------ update Q2; 
  
  subjs = subys = NULL
  for ( r in 1:4) { # diff with itt
    rsub = (1 + (r-1) * long_n) : (long_n +  (r-1) * long_n) 
    subs = length(rsub)*(r-1)+which(stacklong_dat[rsub,]$t==t-1 & stacklong_dat[rsub,]$C==0 & stacklong_dat[rsub,]$A_pl==a_asTreat[r,t] &
                                      stacklong_dat[rsub,]$cumA_pl == sum(a_asTreat[r,1:(t)])) # diff with itt
    subjs = c(subjs, subs)
    
    suby = length(rsub)*(r-1)+which(stacklong_dat[rsub,]$t==t & stacklong_dat[rsub,]$A==a_asTreat[r,t] & 
                                      stacklong_dat[rsub,]$cumA == sum(a_asTreat[r,1:(t)])) # diff with itt
    subys = c(subys, suby)
  }
  
  y = estQ_sim[subys]
  doffset = qlogis(Ql_2r)

  wc2 = matrix(cbind(1, stacklong_dat$w2, coe_cum3, stacklong_dat$w2*coe_cum3), nrow = stacklong_n) # use the same coe_cum for t =1,2,3
  modelu <- glm(y ~ -1 + wc2[subjs,], offset = doffset[subjs], weights =  weights_sim_allr[[t]][subjs], family = quasibinomial()) 
  estQ_sim[which(stacklong_dat$t==t-1)] <- plogis( doffset[which(stacklong_dat$t==t-1)] +
                                                     wc2[which(stacklong_dat$t==t-1),]%*% coef(modelu) )
  
  mean(estQ_sim[which(stacklong_dat$t==t-1)])*(b-a)+a
  
  # mean(estQ_sim[which(stacklong_dat$t==t-1)])*(b-a)+a
  # mean(estQ_sim[which(stacklong_dat$t==t-1)][1:(12792/3)])*(b-a)+a
  # mean(estQ_sim[which(stacklong_dat$t==t-1)][(12792/3+1):(2*12792/3)])*(b-a)+a
  # mean(estQ_sim[which(stacklong_dat$t==t-1)][(2*12792/31):(12792)])*(b-a)+a
  
  
  # ---------------------  compute EIF
  
  Dds[subjs] = (y - estQ_sim[subjs]) * weights_sim_allr[[t]][subjs] 
  # mean(Dds, na.rm= TRUE)
  
  #------------------------------------------ t = 1 --------------------------------------------
  t = 1
  Q_l <- NULL; 
  
  for (r in 1:4) { # diff with itt
    rsub = (1 + (r-1) * long_n) : (long_n +  (r-1) * long_n)
    
    Q_l[[r]] = rep(0, long_n)
    col_Q = c(5:7, 14) #  "w1", "w2" , "L" ,"cumA+", remove  L-, not need "A+" since "A+"="cumA+"
    
    subQ <- which(long_dat$t==t-1 & long_dat$C==0)
    model = glm.fit(y = estQ_sim[rsub][which(long_dat$t==t)], x = cbind(1, long_dat[subQ, col_Q]), family=quasibinomial() ) 
    
    Q_l[[r]] = plogis(as.matrix(cbind(1, long_dat[, 5:7], rep(a_asTreat[r,(t)],long_n))) %*% coef(model))
    
  }
  
  Ql_1r = unlist(Q_l) 
  
  # mean(Ql_1r[which(stacklong_dat$t==t-1)])*(b-a)+a
  # mean(Ql_1r[which(stacklong_dat$t==t-1)][1:(12792/3)])*(b-a)+a
  # mean(Ql_1r[which(stacklong_dat$t==t-1)][(12792/3+1):(2*12792/3)])*(b-a)+a
  # mean(Ql_1r[which(stacklong_dat$t==t-1)][(2*12792/3+1):(12792)])*(b-a)+a
  
  # ------------------------ update Q1; 
  
  subjs = subys = NULL
  for ( r in 1:4) { # diff with itt
    rsub = (1 + (r-1) * long_n) : (long_n +  (r-1) * long_n) 
    subs = length(rsub)*(r-1)+which(stacklong_dat[rsub,]$t==t-1 & stacklong_dat[rsub,]$C==0 & stacklong_dat[rsub,]$cumA_pl == sum(a_asTreat[r,1:(t)]))
    subjs = c(subjs, subs)
    
    suby = length(rsub)*(r-1)+which(stacklong_dat[rsub,]$t==t & stacklong_dat[rsub,]$cumA == sum(a_asTreat[r,1:(t)]))
    subys = c(subys, suby)
  }
  
  # y = .bound(estQ_sim2[which(stacklong_dat$t==t)], c(bound, 1-bound))
  y = estQ_sim[subys]
  doffset = qlogis(Ql_1r)

  wc1 = matrix(cbind(1, stacklong_dat$w2, coe_cum3, stacklong_dat$w2*coe_cum3), nrow = stacklong_n) 
  modelu <- glm(y ~ -1 + wc1[subjs,], offset = doffset[subjs], weights = weights_sim_allr[[t]][subjs], family = quasibinomial()) 
  estQ_sim[which(stacklong_dat$t==t-1)] <- plogis( doffset[which(stacklong_dat$t==t-1)] +
                                                     wc1[which(stacklong_dat$t==t-1),]%*% coef(modelu) )
  
  mean(estQ_sim[which(stacklong_dat$t==t-1)])*(b-a)+a
  
  # mean(estQ_sim[which(stacklong_dat$t==t-1)])*(b-a)+a
  # mean(estQ_sim[which(stacklong_dat$t==t-1)][1:(12792/3)])*(b-a)+a
  # mean(estQ_sim[which(stacklong_dat$t==t-1)][(12792/3+1):(2*12792/3)])*(b-a)+a
  # mean(estQ_sim[which(stacklong_dat$t==t-1)][(2*12792/3+1):(12792)])*(b-a)+a
  
  
  # ---------------------  compute EIF
  Dds[subjs] = (y - estQ_sim[subjs]) * weights_sim_allr[[t]][subjs] 
  # mean(Dds, na.rm= TRUE)
  
  
  
  
  #------------------------------ pooled ltmle for MSM ------------------------------------
  totalq = estQ_sim[which(stacklong_dat$t==0)]
  length(totalq); mean(totalq)*(b-a)+a

  
  #------------ construct cumlative expo for all r
  
  msmA = coe_cum3[which(stacklong_dat$t==0)] # length(msmA) = 4*N
  covmat_sim = cbind(intercep = 1, w2=stacklong_dat[which(stacklong_dat$t==0),]$w2, cumA=msmA,
                     w2_cumA = stacklong_dat[which(stacklong_dat$t==0),]$w2 * msmA)
  # dim(covmat_sim)
  
  totalq_orig = totalq*(b-a)+a
  fit_sim = glm.fit(y = totalq_orig, x = covmat_sim)
  
  ests_msm[j,] = as.vector(coef(fit_sim))
  sum_fit[[200+j]] = summary.glm(fit_sim)
  
  #---------------------- sum EIF on unvaried DD
  id = datl1$ind_id
  ddss = NULL
  for ( r in 1:4) {
    rsub = (1 + (r-1) * long_n) : (long_n +  (r-1) * long_n) 
    ds = sapply(1: length(id), function(x) sum(Dds[rsub][which(stacklong_dat[rsub,]$id == id[x])]))
    ddss = c(ddss, ds)
  }
  # summary(ddss[15001:20000]); length(ddss)
  
  
  #------------------- IF(a[r,])  on  DD
  
  fitt_sim = glm.fit(y = totalq, x = covmat_sim)
  diff1_sim = totalq -predict.glm(fitt_sim)
  diff2_sim = ddss 
  
  sumDif = covmat_sim*(diff1_sim+diff2_sim)
  sumD = sumDif[1:N,]+sumDif[(N+1):(2*N),]+sumDif[(2*N+1):(3*N),]+sumDif[(3*N+1):(4*N),]
  
  
  # ---- normlization matrix (sum_r vv'/n)^(-1)
  covmat_p1 = t(covmat_sim[1:N,]) %*% covmat_sim[1:N,]
  covmat_p2 = t(covmat_sim[(N+1):(2*N),]) %*% covmat_sim[(N+1):(2*N),]
  covmat_p3 = t(covmat_sim[(2*N+1):(3*N),]) %*% covmat_sim[(2*N+1):(3*N),]
  covmat_p4 = t(covmat_sim[(3*N+1):(4*N),]) %*% covmat_sim[(3*N+1):(4*N),]
  sumCOv = covmat_p1+covmat_p2+covmat_p3+covmat_p4
  norm_sim <- solve(sumCOv/length(id))
  
  IFbeta_sim <- (sumD) %*% norm_sim
  IF_sim[[200+j]] = IFbeta_sim
  # summary(IFbeta_sim); dim(IFbeta_sim)
  
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

# save(ests_msm, IF_sim, sds, sds_c, sum_fit, beta_coverage, beta_c_coverage, file = "./Rda/simu_at_500") # w2 ~ Bern(p=0.5+0.8s)
# load(file = "./Rda/simu_at_500")



est_coefs = round(colMeans(ests_msm),4)
# intercep     w2      cumA   w2_cumA 
#  0.3372   1.0048   0.5688   0.2224 

sd_mc = round(sapply(1:4, function(x) sd(ests_msm[,x])),4)
# intercep   w2      cumA   w2_cumA 
# 0.1026   0.1132   0.0526   0.0582

sd_unD = round(colMeans(sds),4)
# intercep       w2     cumA  w2_cumA 
#  0.0807    0.0964    0.0004   0.0006 
sd_c_unD = round(colMeans(sds_c),4)
# intercep       w2     cumA  w2_cumA 
#  0.1000    0.0944    0.0011   0.0018 
# sqrt(colSums(sds^2)/300)
# sqrt(colSums(sds_c^2)/300)
beta_cov = colSums(beta_coverage)/iter
beta_c_cov = colSums(beta_c_coverage)/iter
# ci_L<- ests_msm-1.96*sds
# ci_H <- ests_msm+1.96*sds
# beta_cov = sapply(1:4, function(x) sum(tru_beta[x]>=ci_L[,x] & tru_beta[x]<=ci_H[,x])/iter)
# 
# ci_c_L<- ests_msm-1.96*sds_c
# ci_c_H <- ests_msm+1.96*sds_c
# beta_c_cov = sapply(1:4, function(x) sum(tru_beta[x]>=ci_c_L[,x] & tru_beta[x]<=ci_c_H[,x])/iter)

# source(file = "./test/sim_true_coef.R")
load(file = "./Rda/true_coeff_itt&at")
tru_beta = colMeans(true_coeff_at)
resu_at = round(rbind(tru_beta, est_coefs, sd_mc, sd_unD, sd_c_unD, beta_cov, beta_c_cov), 3) # w2 ~ Bern(p=0.5+0.8s)


#         intercep   w2   cumA w2_cumA
# tru_beta   0.089 1.950 0.588 0.200
# est_coefs  0.094 1.953 0.588 0.196
# sd_mc      0.162 0.242 0.095 0.165
# sd_unD     0.128 0.272 0.096 0.184
# sd_c_unD   0.164 0.284 0.094 0.182
# beta_cov   0.850 0.968 0.966 0.974
# beta_c_cov 0.944 0.976 0.958 0.972



library(xtable)
xtable(resu_at, digits = 2)

