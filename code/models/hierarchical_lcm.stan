data {
    int K;  // number of mixture components (strategies)
    int nsub; // number of subjects 
    int nprob; // number of problems
    array[nsub, nprob] int<lower=0, upper=1> obs;  // observed choices
    array[nprob, K] int<lower=0, upper=1> preds;  // strategy predictions
    }

parameters {
    // hyperparameters
	  vector<lower=0>[K] concentration; // for dirichlet distribution
	  real mu_probit_error; // mean of population-level normal distribution for parameters on probit scale
    real<lower=0> sigma_probit_error;
	  
	  // individual parameters 
    array[nsub] simplex[K] weights; // strategy weights
    vector[nsub] z_error;
}

transformed parameters {
  
    vector[nsub] probit_errors;
    vector<lower=0, upper=0.5>[nsub] errors;

    probit_errors = mu_probit_error + sigma_probit_error * z_error;
    
    
    for (sub in 1:nsub){
    	errors[sub] = 0.5*(Phi(probit_errors[sub])); // trembling-hand error
    }
}

model {
  
  vector[K] log_posterior_mixture;
  
  // hyperpriors 
  
	concentration ~ gamma(2,1); 
  mu_probit_error ~ normal(0,1) ; 
  sigma_probit_error ~ normal(0,1) ;
  
  // individual priors
  z_error ~ normal(0,1);
  
  // likelihood   
  
	for (sub in 1:nsub) {
	  
		weights[sub] ~ dirichlet(concentration);
		// probit_errors[sub] ~ normal(mu_probit_error, sigma_probit_error);
		
		for (prob in 1:nprob) {
		  
			for (k in 1:K) {
			  
			  real choice_prob;
		

        if (obs[sub,prob] == preds[prob,k])
          choice_prob = 1 - errors[sub];
        else
          choice_prob = errors[sub];
			  
			  
    			log_posterior_mixture[k] = log(weights[sub, k]) + log(choice_prob);
			}
    	target += log_sum_exp(log_posterior_mixture);
		}
	}
}
