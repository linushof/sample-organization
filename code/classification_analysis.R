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
  
  lcm_dat <- list(
    obs = obs, 
    preds = preds , 
    nsub = dim(obs)[1]  , # number of subjects in group 
    nprob = dim(obs)[2]  , # number of problems (20)
    K = dim(pred)[2] # number of possible strategies
    )
  
  lcm <- stan(file = "code/models/classification.stan" , 
                #init = inits,
                data = lcm_dat , 
                chains = n_chains , 
                cores = n_chains ,
                warmup = 2500 ,  
                iter = 7500 ,  
                refresh = 750 , 
                seed = 421886
                )
  
  stanfile <- paste0('fits/lcm_',cpx, '_', switch, '.RData')
  save(lcm, file=stanfile)
  
}


# Extract summary matrix
model_params <- summary(lcm)$summary


weights <- model_params |> 
  as.data.frame() |> 
  rownames_to_column("parameter") |> 
  filter(str_detect(parameter, "^weights\\[")) |> # extract weights only
  separate_wider_regex(parameter,
                       patterns= c("weights\\[", participant = "\\d+", ",", strategy = "\\d+","\\]")) |>
  select(participant, strategy, weight = mean)

weights |> 
  mutate(weight1 = if_else(strategy==1,weight,NA) , 
         participant = fct_reorder(participant, weight1, .na_rm=TRUE)) |> 
  ggplot(aes(x=participant, y=weight, fill=strategy)) + 
  geom_col(position = 'fill') + 
  theme_minimal()


# z_error is the standard deviation, not the absolute error
errors <- model_params |> 
  as.data.frame() |> 
  rownames_to_column("parameter") |> 
  filter(str_detect(parameter, "^z_error\\[")) |> # extract weights only
  separate_wider_regex(parameter,
                       patterns= c("z_error\\[", participant = "\\d+", "\\]")) |>
  select(participant, error = mean, error_lower = `2.5%`, error_upper = `97.5%` )

errors |> 
  mutate(participant = fct_reorder(participant, error, .na_rm=TRUE)) |> 
  ggplot(aes(x=participant)) + 
  geom_point(aes(y=error)) +
  geom_errorbar(aes(ymin=error_lower, ymax=error_upper)) +
  theme_minimal()
  
