data {
    int K_sys;  // number of strategies
    int nsub; // number of subjects 
    int nprob; // number of problems
    array[nsub, nprob] int<lower=0, upper=1> obs;  // observed choices
    array[nprob, K_sys] int<lower=0, upper=1> preds;  // strategy predictions
    }

transformed data {
  int K = K_sys + 1; 
}

parameters {
    
    // group-level parameters
	  simplex[K] weights; // strategy weights
	  
	  vector[K_sys] mu_probit_error; // mean of population-level normal distribution for parameters on probit scale
    vector<lower=0>[K_sys] sigma_probit_error;
	  
	  // subject-level parameter
    matrix[nsub, K_sys] z_error;
}

transformed parameters {
  
    matrix[nsub, K_sys] probit_errors;
    matrix<lower=0, upper=0.5>[nsub, K_sys] errors;
    
    for (sub in 1:nsub){
      for (k in 1:K_sys){
        
        probit_errors[sub, k] = mu_probit_error[k] + sigma_probit_error[k] * z_error[sub,k];
    	  errors[sub, k] = 0.5*(Phi(probit_errors[sub,k])); // trembling-hand error
    	  }
    	}
}
    	

model {
  
  vector[K] log_posterior_mixture;
  
  // priors

	weights ~ dirichlet(rep_vector(1.0, K)) ;
	
  mu_probit_error ~ normal(0,1) ; 
  sigma_probit_error ~ exponential(1) ;
  
  
  for (sub in 1:nsub){
    for (k in 1:K_sys){
      z_error[sub, k] ~ normal(0, 1);
    }
  }
  
  // likelihood   

  
	for (sub in 1:nsub) {
	  for (k in 1:K_sys) {
	    
	    log_posterior_mixture[k] = log(weights[k]); // takes the probability for the strategy ...
	  
		  for (prob in 1:nprob) {

        if (obs[sub,prob] == preds[prob,k])
          log_posterior_mixture[k] += log1m(errors[sub, k]); // and multiplies with the likelihood of the response given the strategy
        else
          log_posterior_mixture[k] += log(errors[sub, k]);
			  }
			}
			
		// guessing class
    log_posterior_mixture[K] = log(weights[K]) - nprob * log(2.0);	
			  
		target += log_sum_exp(log_posterior_mixture);
	}
}

generated quantities { 

  vector<lower=0, upper=.5>[K_sys] mu_error;
  
  for(k in 1:K_sys){
    mu_error[k] = 0.5 * Phi(mu_probit_error[k] / sqrt(1 + square(sigma_probit_error[k])));
  }
}
