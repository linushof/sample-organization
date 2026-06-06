# preparation -------------------------------------------------------------

# packages
pacman::p_load(tidyverse, 
               gtools, 
               rstan)

# data 
choices <- read_csv('data/clean/choices_wide_org.csv') # option 1 == 1; option 2 = 0
problems <- read_csv('data/stimuli/observed_samples/test_summaries.csv')

# stan

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## initial values for MCMC
# set_inits <- function(chain_id = 1) {
#   list(alpha_pre=.5, gamma_pre=.5, delta_pre=.5, phi_pre=.01 )
# }
n_chains <- 6 # number of chains
#inits <- lapply(1:n_chains, function(id) set_inits(chain_id = id))


groups <- choices |> 
  distinct(complexity, switching) |> 
  arrange(complexity)

for (g in seq_len(nrow(groups))){ 
  
  cpx <- groups[[g, 'complexity']]
  switch <- groups[[g, 'switching']]
  
  obs <- choices |> 
    filter(complexity==cpx & switching==switch) |> 
    select(problem1:problem20) |> 
    as.matrix()
  
  preds <- problems |> 
    separate(problem, c('complexity', 'problem_match'), sep='_') |> 
    mutate(short = as.numeric(short_option=='option1') , 
           long = as.numeric(long_option=='option1')) |>
    #select(complexity, problem_match, short, long) |> 
    filter(complexity==cpx) |> 
    select(short, long) |>
    as.matrix()
  
  lamix_dat <- list(
    obs = obs, 
    preds = preds , 
    nsub = dim(obs)[1]  , # number of subjects in group 
    nprob = dim(obs)[2]  , # number of problems (20)
    K = dim(pred)[2] # number of possible strategies
    )
  
  lamix <- stan(file = "code/models/classification.stan" , 
                #init = inits,
                data = lamix_dat , 
                chains = n_chains , 
                cores = n_chains ,
                warmup = 2500 ,  
                iter = 7500 ,  
                refresh = 750 , 
                seed = 421886
                )
  
  save(lamix, file = "fits/lamix.RData")
  
  
}
