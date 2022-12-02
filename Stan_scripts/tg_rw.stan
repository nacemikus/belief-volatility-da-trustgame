data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsubj[N];
  int<lower=0, upper=10> transfer[N,T];  // 1: 1 - 4
  int<lower=-1, upper=1> backtransfer[N,T];  // 1,2 
  int<lower=1, upper=2> trustee[N,T];  // 1,2 

  
  int<lower=0, upper=1> sulpiride[N];
  int<lower=0, upper=1> ankk[N];
}
transformed data {
}
parameters {
  // Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters
  vector[4] mu_p;
  vector<lower=0>[4] sigma;
  matrix[4, N] z;
  cholesky_factor_corr[4] L_Omega;
  ordered[10] c;

  
  
  // fixed parameters
  real beta_sul;
  real beta_ankk;
  real beta_sul_ankk;

  real beta_sul_trustee;
  real beta_ankk_trustee;
  real beta_sul_ankk_trustee;

  real beta_noise;
  real beta_mu0;
 

  real beta_sul_ankk_noise;
  real beta_sul_ankk_mu0;
 
  
  real beta_ankk_noise;
  real beta_ankk_mu0;
  
}
transformed parameters {
  // Transform subject-level raw parameters
  vector[N] ag;
  vector[N] al;
  vector[N] noise;
  vector[N] mu0;

  matrix[4, N] r1;
  
  r1 = (diag_pre_multiply(sigma,L_Omega) * z);
 // ordered[10] c[N];

  for (i in 1:N) {
      ag[i] =  Phi_approx(mu_p[1] + r1[1,i] +   
                    (beta_sul + beta_sul_ankk*ankk[i])*sulpiride[i]+
                     beta_ankk*ankk[i] +
                  0.5*(mu_p[2] + r1[2,i] +
                   (beta_sul_trustee +beta_sul_ankk_trustee*ankk[i])*sulpiride[i]+
                   beta_ankk_trustee*ankk[i]) );
                   
                  
                   
      al[i] =   Phi_approx(mu_p[1] + r1[1,i] + 
                    (beta_sul + beta_sul_ankk*ankk[i])*sulpiride[i]+
                     beta_ankk*ankk[i] -
                    0.5*(mu_p[2] + r1[2,i] +
                    (beta_sul_trustee +beta_sul_ankk_trustee*ankk[i])*sulpiride[i]+
                    beta_ankk_trustee*ankk[i]) );
                   
                   
       
      noise[i]  =  10 + mu_p[3] + r1[3,i] + 
                  (beta_noise + beta_sul_ankk_noise*ankk[i])*sulpiride[i]+
                   beta_ankk_noise*ankk[i];

        
      mu0[i]     =Phi_approx( mu_p[4] + r1[4,i]+ 
                  (beta_mu0 + beta_sul_ankk_mu0*ankk[i])*sulpiride[i]+
                   beta_ankk_mu0*ankk[i] );
                   
    

  }
}
model {
  // Hyperparameters
  mu_p[1]  ~ normal(0, 1);
  mu_p[2]  ~ normal(0, 1);
  mu_p[3]  ~ normal(0, 1);
  mu_p[4]  ~ normal(0, 1);

  
  sigma[1]  ~ normal(0, 0.2);
  sigma[2]  ~ normal(0, 0.2);
  sigma[3]  ~ normal(0, 0.2);
  sigma[4]  ~ normal(0, 0.2); // cauchy(0,3); // 
 
  c ~ normal(0,5);
   // individual parameters
  to_vector(z) ~  normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(2);
  // om_mean_pr     ~ normal(0, 1);
  // om_diff_pr   ~ normal(0, 1);
  // noise_pr     ~ normal(0, 1);
  // // betal_pr  ~ normal(0, 1);
  // // beta2_pr  ~ normal(0, 1);
  // mu0_pr     ~ normal(0, 1);
  
  // fixed parameters
  beta_sul ~ normal(0,1);
  beta_ankk ~ normal(0,1);
  beta_sul_ankk ~ normal(0,1);

  beta_sul_trustee ~ normal(0,1);
  beta_ankk_trustee ~ normal(0,1);
  beta_sul_ankk_trustee ~ normal(0,1);

  beta_noise ~ normal(0,1);
  beta_mu0~ normal(0,1);


  beta_sul_ankk_noise ~ normal(0,1);
  beta_sul_ankk_mu0 ~ normal(0,1);

  beta_ankk_noise ~ normal(0,1);
  beta_ankk_mu0 ~ normal(0,1);
 
 
  
  for (i in 1:N) {
    // Define values
    real mu_good;   
    real mu_bad;
    real a;
    real mu;
    real da;
    int transfer_tr;
    // c_pr[i] ~ normal(0, 1);
    // print("hm, c_pr[i] is ", c_pr[i]);
    // Initialize values
    
     // print("hm, c[i] is ", c[i]);
    
    
    mu_good = mu0[i];
    mu_bad =mu0[i];
    
    
    for (t in 1:Tsubj[i])  {
      
      
      transfer_tr = transfer[i,t]+1;
      
      if(trustee[i,t]== 1) { 
     
          // mu  =(delta[i]*mu_good^gamma[i])/(delta[i]*mu_good^gamma[i] + (1-mu_good)^gamma[i]);
           mu  = mu_good;//((mu_good^gam[i]))/((mu_good^gam[i]) + (1-mu_good)^gam[i]);
          transfer_tr ~ ordered_logistic( mu*noise[i], c);

          da = backtransfer[i,t] - mu_good;
          if(da>0) {
            a = ag[i];
          } else {
            a = al[i];
          }
          mu_good = mu_good + a*da;
      } else {
        
            // mu  =(delta[i]*mu_bad^gamma[i])/(delta[i]*mu_bad^gamma[i] + (1-mu_bad)^gamma[i]);
          mu  =mu_bad; // ((mu_bad^gam[i]))/((mu_bad^gam[i]) + (1-mu_bad)^gam[i]);
          transfer_tr ~ ordered_logistic( mu*noise[i], c);
          
          da = backtransfer[i,t] - mu_bad;
          if(da>0) {
            a = ag[i];
          } else {
            a = al[i];
          }
          mu_bad = mu_bad + a*da;
      }
      
   
    } // end of t loop
    
   
    
  } // end of i loop
}

generated quantities {
  // For group level parameters
  // real<lower=0,upper=1> mu_ag;
  // real<lower=0>         mu_beta;
  // real<lower=0,upper=1> mu_al;
  // // real<lower=0>         mu_beta2;
  // 
  // real<lower=0,upper=1> mu_mu0;
  // 
  
  // For log likelihood calculation
  real log_lik[N*T];
  int count_trials;
  // For posterior predictive check
  real y_pred[N,T];
 
  // For model diagnosis
  real mu_good_vect[N,T];
  real mu_bad_vect[N,T];
  // real y_pred_step2[N,T];
  
  
  // Set all posterior predictions to 0 (avoids NULL values)
  for (i in 1:N) {
    for (t in 1:T) {
      y_pred[i,t] = -1;
      mu_good_vect[i,t] =99;
      mu_bad_vect[i,t] = 99;
      // y_pred_step2[i,t] = -1;
    }
  }
 
  // Generate group level parameter values
  // mu_ag     = Phi_approx( mu_p[1] );
  // 
  // mu_al     = Phi_approx( mu_p[2] );
  // // mu_beta2  = exp( mu_p[4] );
  // mu_mu0     = Phi_approx( mu_p[3] );
  // mu_beta  = exp( mu_p[4] );

  { // local section, this saves time and space
  count_trials =0;
  for (i in 1:N) {
    // Define values
       // Define values
    real mu_good;   
    real mu_bad;
    real a;
    real mu;
    real da;
    int transfer_tr;
    
    // Initialize values
    mu_good = mu0[i];
    mu_bad = mu0[i];
    
    
    // log_lik[i] = 0;

    for (t in 1:Tsubj[i])  {
       transfer_tr = transfer[i,t]+1;
       count_trials =count_trials +1;
      // print("lr for subject ", i, " is ", al[i])
     if(trustee[i,t]== 1) { 
     
          // mu  =(delta[i]*mu_good^gamma[i])/(delta[i]*mu_good^gamma[i] + (1-mu_good)^gamma[i]);
          mu  =mu_good ;//((mu_good^gam[i]))/((mu_good^gam[i]) + (1-mu_good)^gam[i]);
          mu_good_vect[i,t] = mu;
          // mu_bad_vect[i,t] = 99;
          log_lik[count_trials] =  ordered_logistic_lpmf(transfer_tr|noise[i]*mu, c);
          
          // generate posterior prediction for current trial
          y_pred[i,t] = ordered_logistic_rng(noise[i]*mu,c) -1;
          
        
          da = backtransfer[i,t] - mu_good;
          if(da>0) {
            a = ag[i];
          } else {
            a = al[i];
          }
          mu_good = mu_good + a*da;
      } else {
        
            // mu  =(delta[i]*mu_bad^gamma[i])/(delta[i]*mu_bad^gamma[i] + (1-mu_bad)^gamma[i]);
           mu  =mu_bad ;// ((mu_bad^gam[i]))/((mu_bad^gam[i]) + (1-mu_bad)^gam[i]);
       
          
          mu_bad_vect[i,t] = mu;
          // mu_good_vect[i,t] = -1;
            // mu  =((mu_good^gam[i]))/((mu_good^gam[i]) + (1-mu_good)^gam[i]);
       
          // mu_bad_vect[i,t] = 99;
          log_lik[count_trials] =  ordered_logistic_lpmf(transfer_tr|noise[i]*mu, c);
          
          
          // generate posterior prediction for current trial
          y_pred[i,t] = ordered_logistic_rng(noise[i]*mu,c) -1;
          
          
          da = backtransfer[i,t] - mu_bad;
          if(da>0) {
            a = ag[i];
          } else {
            a = al[i];
          }
          mu_bad = mu_bad + a*da;
          // print("back tr =", backtransfer[i,t], ", mu_bad=", mu_bad, ", da = ", da, ", investment = ",  y_pred[i,t]);
      }
       
      
      } // end of t loop
        
     
    } // end of i loop
   } // end local
}
