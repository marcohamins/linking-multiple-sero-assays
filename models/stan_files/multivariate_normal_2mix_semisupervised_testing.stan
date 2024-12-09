data {
    int nP; // number of antigens/assays
    int<lower=0> N; // number of training individuals
    int<lower=0> Ntest; // number of testing individuals
    int<lower=0> Nc1; // number of individuals
    int<lower=0> Nc2; // number of individuals
    array[N] vector[nP] y1; // 1st time point
    array[N] vector[nP] y2; // 2nd time point
    array[Ntest] vector[nP] y1_test; // 1st time point
    array[Ntest] vector[nP] y2_test; // 2nd time point
    int K; // number of mixtures
    int C; // number of classes from supervised approach
}

parameters {
    simplex[K] theta; //mixing proportions
    array[C] simplex[K] r; // prob of consistency between class and mixture
    positive_ordered[nP] sigma0; // sd negatives
    positive_ordered[nP] sigma1; // sd positives
    positive_ordered[nP] beta0; // mean y1
    vector[nP] beta1; // slope positive
    array[K] cholesky_factor_corr[nP] L; // cholesky factor
}

model {
 vector[K] log_theta = log(theta);  // cache log calculation
 sigma0 ~ lognormal(0,.5);
 // sigma0 ~ uniform(0.01, 100);
 sigma1 ~ lognormal(0,.5);
 // sigma1 ~ uniform(0.01, 100);
 beta0 ~ uniform(0, 10);
 beta1 ~ uniform(-1, 1);
 r[1,1] ~ uniform(.8,1);
 r[2,2] ~ uniform(.5,1);
 theta[1] ~ uniform(.3,1);
 
 for(i in 1:K){
  L[i] ~ lkj_corr_cholesky(1); 
 }

 for (n in 1:Nc1) {
  vector[K] lps = log_theta;
  real hold = 0;
  lps[1] += multi_normal_cholesky_lpdf(y2[n]|y1[n],diag_pre_multiply(sigma0,L[1])); 
  lps[2] += multi_normal_cholesky_lpdf(y2[n]|y1[n].*beta1 + beta0,diag_pre_multiply(sigma1,L[2]));
  target += log_sum_exp(lps);
  
  for(i in 1:K){
    hold += r[1,i]*softmax(lps)[i];
  }

  target += log(hold);
  }
  
 for (n in 1:Nc2) {
  vector[K] lps = log_theta;
  real hold = 0;
  
  lps[1] += multi_normal_cholesky_lpdf(y2[Nc1 + n]|y1[Nc1 + n],diag_pre_multiply(sigma0,L[1])); 
  lps[2] += multi_normal_cholesky_lpdf(y2[Nc1 + n]|y1[Nc1 + n].*beta1 + beta0,diag_pre_multiply(sigma1,L[2]));
  target += log_sum_exp(lps);
  
  for(i in 1:K){
    hold += r[2,i]*softmax(lps)[i];
  }

  target += log(hold);
  }
}

generated quantities {
   vector[K] log_theta = log(theta);  // cache log calculation
   matrix[K, N] probabilities;
   matrix[K, Ntest] probabilitiesTest;
   
    for(n in 1:N){
      for(k in 1:K){
        probabilities[k,n] = 0;
      }
    }
    
    for(n in 1:N){
      vector[K] lps = log_theta;
      lps[1] += multi_normal_cholesky_lpdf(y2[n]|y1[n],diag_pre_multiply(sigma0,L[1])); 
      lps[2] += multi_normal_cholesky_lpdf(y2[n]|y1[n].*beta1 + beta0,diag_pre_multiply(sigma1,L[2]));
      probabilities[:, n] = softmax(lps);
    }
    
    for(n in 1:Ntest){
      for(k in 1:K){
        probabilitiesTest[k,n] = 0;
      }
    }
 
 for (n in 1:Ntest) {
  vector[K] lps = log_theta;
  lps[1] += multi_normal_cholesky_lpdf(y2_test[n]|y1_test[n],diag_pre_multiply(sigma0,L[1])); 
  lps[2] += multi_normal_cholesky_lpdf(y2_test[n]|y1_test[n].*beta1 + beta0,diag_pre_multiply(sigma1,L[2]));

  probabilitiesTest[:, n] = softmax(lps);
  }
  
}

