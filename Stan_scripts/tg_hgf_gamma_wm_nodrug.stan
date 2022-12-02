data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsubj[N];
  int<lower=0, upper=10> transfer[N,T];  // 1: 1 - 4
  int<lower=-1, upper=1> backtransfer[N,T];  // 1,2 
  int<lower=1, upper=2> trustee[N,T];  // 1,2 
  int<lower=0, upper=1> Tsubj_remove[N];
  
  // int<lower=0, upper=1> sulpiride[N];
  // int<lower=0, upper=1> ankk[N];
  real swm_error[N];
  
}
transformed data {
}
parameters {
  // Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters
  vector[5] mu_p;
  vector<lower=0>[5] sigma;
  matrix[5, N] z;
  cholesky_factor_corr[5] L_Omega;
  ordered[10] c;
 
  // fixed parameters
  // real beta_sul;
  // real beta_ankk;
  // real beta_sul_ankk;
  // // 
  // real beta_sul_trustee;
  // real beta_ankk_trustee;
  // real beta_sul_ankk_trustee;
  // 
  // real beta_noise;
  // real beta_mu0;
  // // 
  // // 
  // // 
  // real beta_sul_ankk_noise;
  // real beta_sul_ankk_mu0;
  // // 
  // // 
  // real beta_ankk_noise;
  // real beta_ankk_mu0;
  // 
  // real beta_gam;
  // real beta_sul_ankk_gam;
  // real beta_ankk_gam;
  // 
  real beta_om_swm;
  real beta_gam_swm;
  real beta_mu0_swm;
  real beta_noise_swm;
  
}
transformed parameters {
  // Transform subject-level raw parameters
  vector[N] om_good;
  vector[N] om_bad;
  vector[N] noise;
  vector[N] mu0;
  vector[N] gam;
  
  matrix[5, N] r1;
  
  r1 = (diag_pre_multiply(sigma,L_Omega) * z);
  // ordered[10] c[N];
  
  for (i in 1:N) {
    om_good[i] =  mu_p[1] + r1[1,i]  -2 + 0.5*(mu_p[2] + r1[2,i]) +  beta_om_swm*swm_error[i]; 
    
    
    
    om_bad[i] =   mu_p[1] + r1[1,i] -2 - 0.5*(mu_p[2] + r1[2,i]) + beta_om_swm*swm_error[i]; 
    
    
    
    noise[i]  =  10 + mu_p[3] + r1[3,i]  + beta_noise_swm*swm_error[i]; 
    
    
    mu0[i]     = mu_p[4] + r1[4,i]+  beta_mu0_swm*swm_error[i]; 
    // 
    gam[i]     =  exp(mu_p[5] + r1[5,i] + beta_gam_swm*swm_error[i]);              
  }
}
model {
  // Hyperparameters
  mu_p[1]  ~ normal(0, 1);
  mu_p[2]  ~ normal(0, 1);
  mu_p[3]  ~ normal(0, 1);
  mu_p[4]  ~ normal(0, 1);
  mu_p[5]  ~ normal(0, 0.5);
  
  sigma[1]  ~ normal(0, 0.2);
  sigma[2]  ~ normal(0, 0.2);
  sigma[3]  ~ normal(0, 0.2);
  sigma[4]  ~ normal(0, 0.2); // cauchy(0,3); // 
  sigma[5]  ~ normal(0, 0.2);
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
  // beta_sul ~ normal(0,1);
  // beta_ankk ~ normal(0,1);
  // beta_sul_ankk ~ normal(0,1);
  // 
  // beta_sul_trustee ~ normal(0,1);
  // beta_ankk_trustee ~ normal(0,1);
  // beta_sul_ankk_trustee ~ normal(0,1);
  // 
  // beta_noise ~ normal(0,1);
  // beta_mu0~ normal(0,1);
  // 
  // 
  // beta_sul_ankk_noise ~ normal(0,1);
  // beta_sul_ankk_mu0 ~ normal(0,1);
  // 
  // beta_ankk_noise ~ normal(0,1);
  // beta_ankk_mu0 ~ normal(0,1);
  // 
  // beta_gam ~ normal(0,1);
  // beta_sul_ankk_gam ~ normal(0,1);
  // beta_ankk_gam ~ normal(0,1);
  
 beta_om_swm~ normal(0,1);
  beta_gam_swm~ normal(0,1);
  beta_mu0_swm~ normal(0,1);
  beta_noise_swm~ normal(0,1);
  // beta_mu0 ~ normal(0,5);
  
  // beta_ami ~ normal(0,2);
  // 
  // beta_dat_nal~ normal(0,2);
  // beta_dat_ami~ normal(0,2);
  // 
  // beta_dat~ normal(0,2);
  
  
  
  
  for (i in 1:N) {
    // Define values
    real mu_good1;
    real mu_good2;
    real mu_good2_prev; 
    real pi_good1;
    real pi_good2;
    real pi_good2_prev;
    
    real mu_bad1;
    real mu_bad2;
    real mu_bad2_prev; 
    real pi_bad1;
    real pi_bad2;
    real pi_bad2_prev;
    
    real da;
    
    real transfer_temp;
    
    real ka_good_par;
    real ka_bad_par;
    
    real om_good_par;
    real om_bad_par;
    
    real mu;
    
    
    int transfer_tr;
    
    // c_pr[i] ~ normal(0, 1);
    // print("hm, c_pr[i] is ", c_pr[i]);
    // Initialize values
    
    // print("hm, c[i] is ", c[i]);
    mu_good2_prev = mu0[i];
    pi_good2_prev = 1/0.1;// 1.0/sa0_2[i];
    
    // mu_bad1 = mu0_1[i];
    // pi_bad1 = 1.0/sa0_1[i];
    mu_bad2_prev = mu0[i];
    pi_bad2_prev = 1/0.1;// 1.0/sa0_2[i];
    
    ka_good_par = 1;//ka_good[i];
    om_good_par = om_good[i];
    
    ka_bad_par = 1;//ka_bad[i];
    om_bad_par = om_bad[i];
    
    // mu_good = mu0[i];
    // mu_bad =mu0[i] ;
    
    
    for (t in 1:Tsubj[i])  {
      
      
      transfer_tr = transfer[i,t]+1;
      
      if(trustee[i,t]== 1) { 
        
        // predictions 
        mu_good2 = mu_good2_prev;
        mu_good1 = 1.0/(1.0+exp(-ka_good_par*mu_good2));
        pi_good1 = 1/(mu_good1*(1-mu_good1));
        
        // model
        
        mu  =((mu_good1^gam[i]))/((mu_good1^gam[i]) + (1-mu_good1)^gam[i]);
        
        if (Tsubj_remove[i] == 0) {
          transfer_tr ~ ordered_logistic( mu*noise[i], c);
        }
        da = backtransfer[i,t] - mu_good1;
        
        // 2nd level
        pi_good2 = 1/(1/pi_good2_prev + exp(om_good_par));
        
        // save for next trial
        pi_good2_prev = pi_good2 + (ka_good_par^2)/pi_good1;
        mu_good2_prev = mu_good2+ ka_good_par/pi_good2_prev*da;
        
        // if (pi_good2 <= 0) 
        //     pi_good2 = 0.001;
        //     // reject("pi_good2 must be positive; negative posterior precision is not allowed: ", pi_good2);
        
      } else {
        
        // predictions 
        mu_bad2 = mu_bad2_prev;
        mu_bad1 = 1.0/(1.0+exp(-ka_bad_par*mu_bad2));
        pi_bad1 = 1/(mu_bad1*(1-mu_bad1));
        
        // model
        // mu = mu_bad1;
        mu  =((mu_bad1^gam[i]))/((mu_bad1^gam[i]) + (1-mu_bad1)^gam[i]);
        
        if (Tsubj_remove[i] == 0) {
          transfer_tr ~ ordered_logistic( mu*noise[i], c);
        }
        da = backtransfer[i,t] - mu_bad1;
        
        // 2nd level
        pi_bad2 = 1/(1/pi_bad2_prev + exp(om_bad_par));
        
        // save for next trial
        pi_bad2_prev = pi_bad2 + (ka_bad_par^2)/pi_bad1;
        mu_bad2_prev = mu_bad2+ ka_bad_par/pi_bad2_prev*da;
        
        // if (pi_bad2 <= 0) 
        //    pi_bad2 = 0.001;
        //    // reject("pi_bad2 must be positive; negative posterior precision is not allowed: ", pi_bad2);
        
      }
      
      
      
    } // end of t loop
    
    
    
  } // end of i loop
}

generated quantities {
  // For group level parameters
  // real mu_om_good;
  // real mu_om_bad;
  // real mu_beta;
  // // real<lower=0,upper=1> mu_al;
  // // real<lower=0>         mu_beta2;
  // 
  // real mu_mu0;
  
  
  // For log likelihood calculation
  real log_lik[N*T];
  int count_trial;
  // For posterior predictive check
  real y_pred[N,T];
  
  // For model diagnosis
  real mu_good_vect[N,T];
  real mu_bad_vect[N,T];
  real mu_good2_vect[N,T];
  real mu_bad2_vect[N,T];
  real pi_good_vect[N,T];
  real pi_bad_vect[N,T];
  real pi_good2_vect[N,T];
  real pi_bad2_vect[N,T];
  // real y_pred_step2[N,T];
  
  
  // Set all posterior predictions to 0 (avoids NULL values)
  for (i in 1:N) {
    for (t in 1:T) {
      y_pred[i,t] = -1;
      mu_good_vect[i,t] = -1;
      mu_bad_vect[i,t] = -1;
      mu_good2_vect[i,t] = -1;
      mu_bad2_vect[i,t] = -1;
      pi_good_vect[i,t] = -1;
      pi_bad_vect[i,t] = -1;
      pi_good2_vect[i,t] = -1;
      pi_bad2_vect[i,t] = -1;
      // y_pred_step2[i,t] = -1;
    }
  }
  for (i in 1:T*N) {
    log_lik[i]=0;
  }
  // Generate group level parameter values
  // mu_om_good     = mu_p[1]-4;
  // mu_om_bad     = mu_p[2]-4;
  // mu_beta  =  mu_p[3] ;
  // // mu_beta2  = exp( mu_p[4] );
  // mu_mu0     =  mu_p[4];
  
  { // local section, this saves time and space
  count_trial= 0;
  for (i in 1:N) {
    // Define values
    // Define values
    real mu_good1;
    real mu_good2;
    real mu_good2_prev; 
    real pi_good1;
    real pi_good2;
    real pi_good2_prev;
    
    real mu_bad1;
    real mu_bad2;
    real mu_bad2_prev; 
    real pi_bad1;
    real pi_bad2;
    real pi_bad2_prev;
    
    real da;
    
    real transfer_temp;
    
    real ka_good_par;
    real ka_bad_par;
    
    real om_good_par;
    real om_bad_par;
    
    real mu;
    int transfer_tr;
    
    
    // Initialize values
    mu_good2_prev = mu0[i];
    pi_good2_prev = 1/0.1;// 1.0/sa0_2[i];
    
    // mu_bad1 = mu0_1[i];
    // pi_bad1 = 1.0/sa0_1[i];
    mu_bad2_prev = mu0[i];
    pi_bad2_prev = 1/0.1;// 1.0/sa0_2[i];
    
    ka_good_par = 1;//ka_good[i];
    om_good_par = om_good[i];
    
    ka_bad_par = 1;//ka_bad[i];
    om_bad_par = om_bad[i];
    
    
    
    // log_lik[i] = 0;
    
    for (t in 1:Tsubj[i])  {
      transfer_tr = transfer[i,t]+1;
      count_trial = count_trial +1; 
      // print("lr for subject ", i, " is ", al[i])
      if(trustee[i,t]== 1) { 
        
        // predictions 
        mu_good2 = mu_good2_prev;
        
        mu_good1 = 1.0/(1.0+exp(-ka_good_par*mu_good2));
        pi_good1 = 1/(mu_good1*(1-mu_good1));
        
        // model
        mu  =((mu_good1^gam[i]))/((mu_good1^gam[i]) + (1-mu_good1)^gam[i]);
        
      
          // log_lik[i] = log_lik[i] + ordered_logistic_lpmf(transfer_tr|noise[i]*mu, c);
          log_lik[count_trial]  = ordered_logistic_lpmf(transfer_tr|noise[i]*mu, c);
          // generate posterior prediction for current trial
          y_pred[i,t] = ordered_logistic_rng(noise[i]*mu,c) -1;
        
        
        
        
        da = backtransfer[i,t] - mu_good1;
        
        // 2nd level
        pi_good2 = 1/(1/pi_good2_prev + exp(om_good_par));
        
        
        // save vars
        mu_good_vect[i,t]= mu_good1;
        mu_good2_vect[i,t]= mu_good2;
        pi_good_vect[i,t]= pi_good1;
        pi_good2_vect[i,t]= pi_good2;
        // save for next trial
        pi_good2_prev = pi_good2 + (ka_good_par^2)/pi_good1;
        mu_good2_prev = mu_good2+ ka_good_par/pi_good2_prev*da;
        
        // print("mu_good2 is ", mu_good2)
        // if (pi_good2 <= 0) 
        //     pi_good2 = 0.001;
        
      } else {
        
        // predictions 
        mu_bad2 = mu_bad2_prev;
        
        mu_bad1 = 1.0/(1.0+exp(-ka_bad_par*mu_bad2));
        pi_bad1 = 1/(mu_bad1*(1-mu_bad1));
        
        // model
        mu  =((mu_bad1^gam[i]))/((mu_bad1^gam[i]) + (1-mu_bad1)^gam[i]);
        
        
          // log_lik[i] = log_lik[i] + ordered_logistic_lpmf(transfer_tr|noise[i]*mu, c);
        log_lik[count_trial]  = ordered_logistic_lpmf(transfer_tr|noise[i]*mu, c);   
          
          
          // generate posterior prediction for current trial
          y_pred[i,t] = ordered_logistic_rng(noise[i]*mu,c) -1;
       
        da = backtransfer[i,t] - mu_bad1;
        
        // 2nd level
        pi_bad2 = 1/(1/pi_bad2_prev + exp(om_bad_par));
        // save vars
        mu_bad_vect[i,t]= mu_bad1;
        mu_bad2_vect[i,t]= mu_bad2;
        pi_bad_vect[i,t]= pi_bad1;
        pi_bad2_vect[i,t]= pi_bad2;
        // save for next trial
        pi_bad2_prev = pi_bad2 + (ka_bad_par^2)/pi_bad1;
        mu_bad2_prev = mu_bad2+ ka_bad_par/pi_bad2_prev*da;
        
        // if (pi_bad2 <= 0) 
        //     pi_bad2 = 0.001;
        
      }
      
      
      
    } // end of t loop
    
    
  } // end of i loop
  } // end local
}

