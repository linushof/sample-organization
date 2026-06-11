# preamble ----------------------------------------------------------------

# packages
pacman::p_load(gtools, rstan)

redo_analysis <- FALSE
current_model <- 'hierarchical_lcm' 


current_model_fits <- paste0('fits/summaries/',current_model,'.RData')
groups <- expand_grid(current_model=current_model , 
                      complexity=c('simple', 'complex') ,
                      switching=c('none', 'infrequent','frequent')) |> 
  mutate(file = paste0('fits/',current_model,'_',complexity, '_', switching, '.RData'))


if(redo_analysis){file.remove(groups$file,  current_model_fits)}

missing_fits <- anyNA(if_else(file.exists(groups$file),TRUE,NA))

# preparation -------------------------------------------------------------

  
if(missing_fits){
  
  # Stan settings
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  lcm_stan <- paste0('code/models/',current_model,'.stan') 
  n_chains <- 6 # number of parallel chains
  n_iter <- 7500
  n_disc <- 2500
  
  # Data sets
  choices <- read_csv('data/clean/choices_wide_org.csv') # option 1 == 1; option 2 = 0
  problems <- read_csv('data/stimuli/observed_samples/test_summaries.csv')
  #samples <- read_csv('data/stimuli/observed_samples/test_all.csv')
  
  for (g in seq_len(nrow(groups))){ 
    
    cpx <- groups[[g, 'complexity']]
    switch <- groups[[g, 'switching']]
    stanfile <- groups[[g, 'file']]
    
    # refit model
 
    obs <- choices |> 
      filter(complexity==cpx & switching==switch) |> 
      select(problem1:problem20) |> 
      as.matrix()
    
    preds <- problems |> 
      separate(problem, c('complexity', 'problem_match'), sep='_') |> 
      mutate(short = as.numeric(short_option=='option1') , 
             long = as.numeric(long_option=='option1')) |>
      filter(complexity==cpx) |> 
      select(short, long) |>
      as.matrix()
    
    lcm_dat <- list(obs = obs ,
                    preds = preds ,
                    nsub = dim(obs)[1] , # number of subjects in group
                    nprob = dim(obs)[2] , # number of problems (20)
                    K = dim(preds)[2]) # number of possible strategies
    
    lcm <- stan(file = lcm_stan, data = lcm_dat ,
                chains = n_chains, iter = n_iter, warmup = n_disc ,
                cores = n_chains, refresh = (n_iter/10) ,
                seed = 421886)
    save(lcm, file=stanfile)
  }
}


# prepare data sets

if(!file.exists(current_model_fits)){
  
  weights_list <- errors_list <- vector('list', length = nrow(groups))
  for (g in seq_len(nrow(groups))){ 
    
    cpx <- groups[[g, 'complexity']]
    switch <- groups[[g, 'switching']]
    stanfile <- groups[[g, 'file']]
    load(stanfile)
    
    # extract relevant parameters (strategy weights and parameters)
    model_params <- summary(lcm)$summary
    
    if(current_model=='hierarchical_lcm'){ 
      weights <- model_params |> 
        as.data.frame() |> 
        rownames_to_column("parameter") |> 
        filter(str_detect(parameter, "^weights\\[")) |> # extract weights only
        separate_wider_regex(parameter,
                             patterns= c("weights\\[", participant = "\\d+", ",", strategy = "\\d+","\\]")) |>
        select(participant, strategy, weight = mean) |> 
        mutate(complexity=cpx, 
               switching=switch)
      
      # CAVE: z_error is the standard deviation, not the absolute error
      errors <- model_params |> 
        as.data.frame() |> 
        rownames_to_column("parameter") |> 
        filter(str_detect(parameter, "^z_error\\[")) |> # extract weights only
        separate_wider_regex(parameter,
                             patterns= c("z_error\\[", participant = "\\d+", "\\]")) |>
        select(participant, error = mean, error_lower = `2.5%`, error_upper = `97.5%` ) |> 
        mutate(complexity=cpx, 
               switching=switch)
      
      } else{# pooling
        
        # continue here
        
      }
    
    weights_list[[g]] <- weights
    errors_list[[g]] <- errors
  }
  
  
  # make summary data frame
  
  all_weights <- weights_list |> 
    bind_rows() |> 
    mutate(strategy=if_else(strategy=='1', 'Roundwise', 'Summary'), 
           weight_round = if_else(strategy=='Roundwise', weight, 1-weight)) 
  
  all_errors <- errors_list |> 
    bind_rows() |> 
    mutate(error_transformed = pnorm(error)) |> 
    left_join(all_weights, by=join_by(participant, complexity, switching)) |> 
    distinct(participant, complexity, switching, .keep_all=T) |> 
    select(!weight)
  
  
  save(all_weights, all_errors, file=current_model_fits)
  
}
