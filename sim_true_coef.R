

true_r <- function(seed, N, study_n, r, study)  {
  
  set.seed(seed)
  # Scenario 1: D is independent
  clusterN = N/study_n
  study_id = c(rep(1:clusterN, each=study_n))
  ind_id = c(1:N)
  
  # generate study-level covariate 
  s_bar = rnorm(clusterN, mean=0.5, sd=1)
  s = rep(s_bar, each=study_n)
  
  w1<-rnorm(N, mean=0.5*s, sd=0.7)
  w2<-rbinom(N, 1, plogis(0.5+0.8*s))
  L0 = rnorm(N, mean=0.2+0.2*w1+0.3*w2, sd=0.7)
  
  # t1 = rep(1, N)
  if (study == "at") {
    if (r == 1 ) A1 = A2 = rep(1, N) else if (r == 2) {A1 = rep(0, N); A2 = rep(1,N)
    } else if (r == 3) {A1 = A2 = rep(0, N)} else if (r==4) {A1 = rep(1, N); A2 = rep(0,N)} else stop
  }
    
  if (study == "itt") {
    if (r == 1 ) A1 = A2 = rep(1, N) else if (r == 2) {A1 = rep(0, N); A2 = rep(1,N)
    } else if (r == 3) {A1 = A2 = rep(0, N)} else stop
  }
  
  L1 = rnorm(N, mean=0.2*w1+0.3*w2+0.5*L0 + 0.5*A1, sd=1)
  L2 = rnorm(N, mean=0.2*w1+0.3*w2+0.7*L1 + 0.4*A2, sd=0.8)
  
  cumA = A1+A2
  
  err <- rnorm(N,0,1)
  Y = 0.2*s+0.6*w1+1*w2+0.4*L1+0.5*L2+0.3*cumA+0.2*cumA*w2+ err
  
  return(data.frame(Y, w1, w2, L0, L1, L2, cumA))
}


tru_coeff <- function(tru_iter, study) {
  
  truY= NULL
  tru_inds = matrix(NA, nrow = tru_iter, ncol = 4)
  
  for (i in 1:tru_iter) {
    
    N=5000; study_n=100
    tr1 = true_r(seed=seeds[i], N=N, study_n=study_n, r=1, study = study)
    tr2 = true_r(seed=seeds[i], N=N, study_n=study_n, r=2, study = study)
    tr3 = true_r(seed=seeds[i], N=N, study_n=study_n, r=3, study = study)
    truY = c(tr1$Y, tr2$Y, tr3$Y)
    msmA = rep(c(2,1,0), each = N) 
    covmat = cbind(intercep = 1, w2=rep(tr1$w2, 3), cumA=msmA, w2_cumA = rep(tr1$w2, 3) * msmA)
    
    if (study == "at") {
      tr4 = true_r(seed=seeds[i], N=N, study_n=study_n, r=4, study = study)
      truY = c(tr1$Y, tr2$Y, tr3$Y, tr4$Y)
      msmA = rep(c(2,1,0,1), each = N) 
      covmat = cbind(intercep = 1, w2=rep(tr1$w2, 4), cumA=msmA, w2_cumA = rep(tr1$w2, 4) * msmA)
    }
    
    fit_tru = glm.fit(y = truY, x = covmat)
    tru_inds[i,] = as.vector(coef(fit_tru))
  }
  return(tru_inds)
  
}


n_iter=100000
seeds = sample(1:n_iter, n_iter, replace = FALSE)

true_coeff_itt = tru_coeff(tru_iter=n_iter, study="itt")
colMeans(true_coeff_itt)
# 0.05956237 1.94911048 0.58750000 0.20000000

true_coeff_at = tru_coeff(tru_iter=n_iter, study="at")
colMeans(true_coeff_at)
# 0.08872904 1.94911048 0.58750000 0.20000000

# save(true_coeff_itt, true_coeff_at, file = "./Rda/true_coeff_itt&at") 




