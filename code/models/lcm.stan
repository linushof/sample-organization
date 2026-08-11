data {
    int K_sys;  // number of systematic strategies
    int nsub; // number of subjects 
    int nprob; // number of problems per subject
    array[nsub, nprob] int<lower=0, upper=1> obs;  // observed choices
    array[nprob, K_sys] int<lower=0, upper=1> preds;  // predictions of systematic strategies
    }

transformed data {
  int K = K_sys + 1; // adds guessing as K-th strategy
}

parameters {
    
	  simplex[K] weights; // group-level class-membership probabilities
	  
	  // reparameterized probit error
	  real mu_probit_error; // group-level mean probit error 
    real<lower=0> sigma_probit_error; // group-level standard deviation of probit errors
	  vector[nsub] z_error; // z-standardized (indiviual) probit error
}

transformed parameters {
  
    vector[nsub] probit_errors;
    vector<lower=0, upper=0.5>[nsub] errors;
    
    probit_errors = mu_probit_error + sigma_probit_error * z_error; // reverse z-standardization to obtain probit errors
    
    for (sub in 1:nsub){
    	errors[sub] = 0.5*(Phi(probit_errors[sub])); // inverse probit to obtain individual errors
    }
}

model {
  
  // priors

	weights ~ dirichlet(rep_vector(1.0, K)) ; //  rep_vector(1.0, K) is the concentration vector alpha 
	
  mu_probit_error ~ normal(0,1) ;
  sigma_probit_error ~ exponential(1) ;
  z_error ~ normal(0,1);
  
  // likelihood   
  vector[K] log_posterior_mixture;
  
	for (sub in 1:nsub) {
	  for (k in 1:K_sys) {
	    
	    // takes the log-probability for the strategy ...
	    
	    log_posterior_mixture[k] = log(weights[k]); 
	  
		  for (prob in 1:nprob) {

        // and adds the log-likelihood of the observed response given the strategy prediction
        
        if (obs[sub,prob] == preds[prob,k])
          log_posterior_mixture[k] += log1m(errors[sub]); 
        else
          log_posterior_mixture[k] += log(errors[sub]);
			  }
			}
			
		// add guessing
		
    log_posterior_mixture[K] = log(weights[K]) - nprob * log(2.0);	
			  
		target += log_sum_exp(log_posterior_mixture);
	}
}

generated quantities { 

  // get group-level mean error 
  real<lower=0, upper=.5> mu_error; 
  mu_error = 0.5 * Phi(mu_probit_error / sqrt(1 + square(sigma_probit_error)));
  
  // get posterior predicted class membership for each individual
  array[nsub] vector[K] class_prob;
  
  for (sub in 1:nsub) {
    
    vector[K] log_prob;
    
    for (k in 1:K_sys) {
      
      log_prob[k] = log(weights[k]);
      
        for (prob in 1:nprob) {

          if (obs[sub,prob] == preds[prob,k])
          
            log_prob[k] += log1m(errors[sub]);
            
          else
          
            log_prob[k] += log(errors[sub]);

    }
  }

  // guessing strategy
  log_prob[K] = log(weights[K]) - nprob * log(2);

  class_prob[sub] = softmax(log_prob);
}

}
  
