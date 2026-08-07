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
    simplex[K] mu_weights;
    real<lower=0> kappa;
  
  	// subject-level parameters 
    array[nsub] simplex[K] weights; // strategy weights

	  
	  real mu_probit_error; // mean of population-level normal distribution for parameters on probit scale
    real<lower=0> sigma_probit_error;
	  
	  // subject-level parameter
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
  
  // priors

  mu_weights ~ dirichlet(rep_vector(1.0, K));

  kappa ~ gamma(2,.5); // use stronger regularization

  for (s in 1:nsub)
    weights[s] ~ dirichlet(kappa * mu_weights);
	
  mu_probit_error ~ normal(0,1) ; 
  sigma_probit_error ~ normal(0,1);
  z_error ~ normal(0,1);
  
  // likelihood   
  for (sub in 1:nsub) {
  for (prob in 1:nprob) {

    vector[K] lp;

    for (k in 1:K_sys) {

      lp[k] = log(weights[sub][k]);

      if (obs[sub,prob] == preds[prob,k])
        lp[k] += log1m(errors[sub]);
      else
        lp[k] += log(errors[sub]);
    }

    // guessing
    lp[K] = log(weights[sub][K]) - log(2.0);

    target += log_sum_exp(lp);
  }
}
  
  
  
}

generated quantities { 

  real<lower=0, upper=.5> mu_error; 
  mu_error = 0.5 * Phi(mu_probit_error / sqrt(1 + square(sigma_probit_error)));
  
  }


