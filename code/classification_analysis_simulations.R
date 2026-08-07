# prep --------------------------------------------------------------------
rm(list = ls())

# pkgs
pacman::p_load(rstan, gtools, bayesplot, ggdist, tidybayes, tidyverse)


# stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# data
choices <- read_csv('data/choices.csv')
problems <- read_csv('data/problems.csv')

# design

nsub <- 35
nprob <- 20
K_sys <- 2
K <- K_sys + 1

preds <- problems |>
  filter(complexity=='complex') |> 
  select(short=short_option, long=long_option) |>
  as.matrix()


# models

models <- c(paste0(rep('lcm_',4), 1:4), paste0(rep('mixture_',2), 1:2))


# Simulations -------------------------------------------------------------

mu_probit_error <- -1
sigma_probit_error <- .5


probit_errors <- rnorm(
  nsub,
  mu_probit_error,
  sigma_probit_error
)

errors <- 0.5 * pnorm(probit_errors)

# SIM 1: Latent class model -------------------------------------------------------------------


weights <- c(.45, .35, .20)


obs <- matrix(0, nsub, nprob)
true_strategy <- integer(nsub)

set.seed(123)
for(sub in 1:nsub){
  
  ## choose one strategy for this subject
  
  k <- sample(
    1:K,
    size = 1,
    prob = weights
  )
  
  true_strategy[sub] <- k
  
  for(prob in 1:nprob){
    
    if(k == K){
      
      ## guessing strategy
      obs[sub,prob] <- rbinom(1,1,.5)
      
    } else {
      
      ## systematic strategy with trembling-hand error
      
      if(runif(1) < errors[sub]){
        
        obs[sub,prob] <- 1 - preds[prob,k]
        
      } else {
        
        obs[sub,prob] <- preds[prob,k]
        
      }
      
    }
  }
}

stan_data_lcm2 <- list(
  K_sys = K_sys,
  nsub  = nsub,
  nprob = nprob,
  obs   = obs,
  preds = preds
)

# Fit
fit_lcm2 <- stan('code/models/lcm_2.stan',
                 data = stan_data_lcm2, 
                 chains = 6,
                 cores=6,
                 iter = 2000,
                 warmup = 1000,
                 control=list(adapt_delta=.995))

posts <- rstan::extract(fit_lcm2)

summary(fit_lcm2)$summary |> View()



str(posts)

for(i in 1:35){
  
 print(round(colMeans(posts$class_prob[,i,]),4)) 
  
}

fit_lcm2_with_mixture <- stan('code/models/sw_version.stan',
                              data = stan_data_lcm2, 
                              chains = 6,
                              cores=6,
                              iter = 2000,
                              warmup = 1000,
                              control=list(adapt_delta=.995))

# check recovery

lcm2_params <- summary(fit_lcm2)$summary |> 
  as.data.frame() |> 
  rownames_to_column("parameter")

model_errors <- lcm2_params |> 
  filter(str_detect(parameter, "^errors\\[")) |> 
  separate_wider_regex(parameter,
                       patterns= c("^errors\\[", participant = "\\d+", "\\]")) |>
  select(participant, error = mean, lower_CI = `2.5%`, upper_CI = `97.5%` )

model_errors |> 
  mutate(true=errors) |> 
  ggplot(aes(true, error)) + 
  geom_point() + 
  geom_linerange(aes(ymin = lower_CI, ymax = upper_CI)) +
  geom_abline(intercept=0, slope=1) +
  scale_x_continuous(limits=c(0,1)) + 
  scale_y_continuous(limits=c(0,1)) +
  theme_minimal()


model_weigths <- lcm2_params |> 
  filter(str_detect(parameter, "^weights\\[")) |> 
  separate_wider_regex(parameter,
                       patterns= c("^weights\\[", strategy = "\\d+", "\\]")) |>
  select(strategy, error = mean, lower_CI = `2.5%`, upper_CI = `97.5%` )



# SIM 2: mixture model -------------------------------------------------------



mu_weights <- c(4, 3, 2)


set.seed(1233)


weights <- rdirichlet(nsub, mu_weights)
colMeans(weights)


obs <- matrix(0, nsub, nprob)

# (optional) store the strategy used on each trial
true_strategy <- matrix(NA, nsub, nprob)

for (s in 1:nsub){
  
  for (p in 1:nprob){
    
    ## Choose strategy for this problem
    k <- sample(
      1:K,
      size = 1,
      prob = weights[s,]
    )
    
    true_strategy[s,p] <- k
    
    ## Guessing strategy
    if(k == K){
      
      obs[s,p] <- rbinom(1,1,.5)
      
    } else {
      
      ## Systematic strategy with trembling-hand error
      if(runif(1) < errors[s]){
        
        obs[s,p] <- 1 - preds[p,k]
        
      } else {
        
        obs[s,p] <- preds[p,k]
        
      }
    }
  }
}


stan_data_mixture <- list(K_sys = K_sys,
                          nsub = nsub,
                          nprob = nprob,
                          obs = obs,
                          preds = preds)


# Fit
fit_mixture <- stan('code/models/sw_version.stan',
                    data = stan_data_mixture, 
                    chains = 6,
                    cores=6,
                    iter = 2000,
                    warmup = 1000,
                    control=list(adapt_delta=.995))

fit_mixture_with_lcm <- stan('code/models/lcm_2.stan',
                             data = stan_data_mixture, 
                             chains = 6,
                             cores=6,
                    iter = 2000,
                    warmup = 1000,
                    control=list(adapt_delta=.995))


posterior <- as.array(fit_mixture)

np <- nuts_params(fit_mixture)

mcmc_pairs(
  posterior,
  pars = c(
    "mu_weights[1]",
    "mu_weights[2]",
    "mu_weights[3]",
    "mu_probit_error"
  ),
  np = np
)
weights
mcmc_trace(fit, regex_pars = "mu_weights")
mcmc_rank_overlay(fit, regex_pars = "mu_weights")



load("fits/lcm_2_complex_frequent.RData")
posterior <- as.array(lcm)

np <- nuts_params(lcm)

mcmc_pairs(
  posterior,
  pars = c(
    "weights[1]",
    "weights[2]",
    "weights[3]",
    "mu_probit_error"
  ),
  np = np
)
