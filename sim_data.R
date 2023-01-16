

datal <- function(seed, N, study_n, study)  {
  
  set.seed(seed)
  # Scenario 1: D is independent
  clusterN = N/study_n
  study_id = c(rep(1:clusterN, each=study_n))
  ind_id = c(1:N)
  
  # generate study-level covariate V
  s_bar = rnorm(clusterN, mean=0.5, sd=1)
  s = rep(s_bar, each=study_n)
  
  # generate 9000 individual covariates W
  w1<-rnorm(N, mean=0.5*s, sd=0.7)
  w2<-rbinom(N, 1, plogis(0.5+0.8*s))
  # w2<-rbinom(N, 1, 0.5)
  L0 = rnorm(N, mean=0.2+0.2*w1+0.3*w2, sd=0.7)
  A0 = C0 = rep(0, N)
  
  # t1 = rep(1, N)
  p_A1 = plogis(0.7*w1+0.8*w2+0.9*L0)
  A1 = rbinom(N, 1, p_A1)
  
  L1 = rnorm(N, mean=0.2*w1+0.3*w2+0.5*L0 + 0.5*A1, sd=1)
  p_C1 = plogis(-2.5+0.2*w1+0.4*w2+0.3*L1+0.2*A1)
  C1 = rbinom(N, 1,  p_C1)
  #A1
  p_A2 = plogis(0.5*w1+0.3*w2+0.7*L1)
  
  if (study == "itt") A2 = ifelse( A1==1, 1, rbinom(N, 1, p_A2)) # itt: if A1=1, A2=1
  if (study == "at") A2 = rbinom(N, 1, p_A2) 
  
  L20 = rnorm(N, mean=0.2*w1+0.3*w2+0.7*L1+0.4*A2, sd=0.8)
  
  p_C2 = plogis(-3+0.2*w1+0.4*w2+0.3*L20+0.4*A2)# ; mean(p_C2, na.rm = TRUE)
  C2 = ifelse( C1==0,  rbinom(N, 1, p_C2), NA)
  total_t = ifelse(C1 ==0, 2, 1)
  
  cumA = ifelse( C1==0, A1+A2, A1)
  
  err <- rnorm(N,0,1)
  Y = ifelse(C2==0, 0.2*s+0.6*w1+1*w2+0.4*L1+0.5*L20+0.3*cumA+0.2*cumA*w2+ err, NA) 
  
  A2 = ifelse( C1==0, A2, NA)
  L2 = ifelse( C1==0, L20, NA)
  
  datal = data.frame(ind_id, study_id, total_t, s, w1, w2, L0, A0, C0, A1, L1, C1, A2, cumA, L2, C2, Y)
  
  return(datal)
}



# function to change the wide data to long data
longData = function(data, study) {
  
  t = id = study_id = Y = s = w1 = w2 = L = L_minus = A = A_minus = A_pl = cumA = cumA_minus = cumA_pl = C = Y = NULL
  
  for (i in 1:dim(datl1)[1]) { 
    ind_t=0:datl1$total_t[i]; t = c(t, ind_t)
    id = c(id, rep(datl1$ind_id[i], length(ind_t)))
    study_id = c(study_id, rep(datl1$study_id[i], length(ind_t)))
    s = c(s, rep(datl1$s[i], length(ind_t)))
    w1 = c(w1, rep(datl1$w1[i], length(ind_t)))
    w2 = c(w2, rep(datl1$w2[i], length(ind_t)))
    
    if (length(ind_t)==2) {
      ll = c(datl1$L0[i], datl1$L1[i]) 
      llmin = c(datl1$L0[i], datl1$L0[i])
      aa = c(datl1$A0[i], datl1$A1[i])
      aamin = c(datl1$A0[i], datl1$A0[i])
      aapl = c(datl1$A1[i], datl1$A1[i])
      cc = c(datl1$C0[i], datl1$C1[i])
      cuma = c(0, ifelse(datl1$A1[i]==1, 1, 0))
      cumamin = c(0, 0)
      cumapl = c(ifelse(datl1$A1[i]==1, 1, 0), ifelse(datl1$A1[i]==1, 1, 0))
      y = c(NA, NA)
      
    } else {
      ll = c(datl1$L0[i], datl1$L1[i], datl1$L2[i])
      llmin = c(datl1$L0[i], datl1$L0[i], datl1$L1[i])
      aa = c(datl1$A0[i], datl1$A1[i], datl1$A2[i])
      aamin = c(datl1$A0[i], datl1$A0[i], datl1$A1[i])
      aapl = c(datl1$A1[i], datl1$A2[i], datl1$A2[i])
      cc = c(datl1$C0[i], datl1$C1[i], datl1$C2[i])
      
      if (study == "itt") {
        
        if (sum(aa) ==1) {cuma = c(0, 0, 1); cumapl = c(0,1,1)
      } else if (sum(aa) ==0) {cuma =cumapl =  c(0, 0, 0)
      } else {cuma = c(0, 1, 2); cumapl = c(1,2,2) }
        
      }
      
      if (study == "at") {
        
        if (sum(aa) ==1) {
          if (datl1$A1[i]==1) {cuma = c(0, 1, 1); cumapl = c(1,1,1)
          } else {cuma = c(0, 0, 1); cumapl = c(0,1,1)}
        } else if (sum(aa) ==0) {cuma =cumapl =  c(0, 0, 0)
        } else {cuma = c(0, 1, 2); cumapl = c(1,2,2) }
        
      }
      
      cumamin = c(0, 0, ifelse(datl1$A1[i]==1, 1, 0))
      y = c(NA, NA, datl1$Y[i])
    }
    L = c(L, ll)
    L_minus = c(L_minus, llmin)
    A = c(A, aa)
    A_minus = c(A_minus, aamin)
    A_pl = c(A_pl, aapl)
    C = c(C, cc)
    cumA = c(cumA, cuma)
    cumA_minus = c(cumA_minus, cumamin)
    cumA_pl = c(cumA_pl, cumapl)
    Y= c(Y, y)
  }
  
  long_dat = data.frame(id, study_id, t, s, w1, w2, L, L_minus, A, A_minus, A_pl, cumA, cumA_minus, cumA_pl, C, Y)
  
  
  return(long_dat)
  
}






