# preparation -------------------------------------------------------------

# packages
pacman::p_load(tidyverse, 
               gtools, 
               rstan, patchwork)

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
    K = dim(preds)[2] # number of possible strategies
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


weights_list <- vector('list', length = nrow(groups))
errors_list <- vector('list', length = nrow(groups))

for (g in seq_len(nrow(groups))){
  
  cpx <- groups[[g, 'complexity']]
  switch <- groups[[g, 'switching']]
  
  # load stan fits
  load(paste0('fits/lcm_',cpx, '_', switch, '.RData'))
  
  # Extract summary matrix
  model_params <- summary(lcm)$summary
  
  weights <- model_params |> 
    as.data.frame() |> 
    rownames_to_column("parameter") |> 
    filter(str_detect(parameter, "^weights\\[")) |> # extract weights only
    separate_wider_regex(parameter,
                         patterns= c("weights\\[", participant = "\\d+", ",", strategy = "\\d+","\\]")) |>
    select(participant, strategy, weight = mean) |> 
    mutate(complexity=cpx, 
           switching=switch)

  
  # z_error is the standard deviation, not the absolute error
  errors <- model_params |> 
    as.data.frame() |> 
    rownames_to_column("parameter") |> 
    filter(str_detect(parameter, "^z_error\\[")) |> # extract weights only
    separate_wider_regex(parameter,
                         patterns= c("z_error\\[", participant = "\\d+", "\\]")) |>
    select(participant, error = mean, error_lower = `2.5%`, error_upper = `97.5%` ) |> 
    mutate(complexity=cpx, 
           switching=switch)
  
  weights_list[[g]] <- weights
  errors_list[[g]] <- errors
  }


all_weights <- weights_list |> 
  bind_rows() |> 
  mutate(strategy=if_else(strategy=='1', 'Roundwise', 'Summary')) 

plot_weights_simple <- all_weights |> 
  filter(complexity=='simple') |> 
  mutate(weight_round = if_else(strategy=='Roundwise', weight, 1-weight)) |> 
  ggplot(aes(x=reorder_within(participant, by=weight_round, within=switching), y=weight, fill=strategy)) + 
  geom_col(position = 'fill', col='black') + 
  scale_x_reordered() +
  facet_wrap(~switching, scales='free', drop=TRUE, nrow=3) +
  labs(title='Low Complexity', 
       x='Participants', 
       y='Strategy Weights', 
       fill='Decision Strategy') +
  theme_bw() +
  scale_y_continuous(breaks = seq(0,1,length.out=3)) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank() ) 

plot_weights_complex <- all_weights |> 
  filter(complexity=='complex') |> 
  mutate(weight_round = if_else(strategy=='Roundwise', weight, 1-weight)) |> 
  ggplot(aes(x=reorder_within(participant, by=weight_round, within=switching), y=weight, fill=strategy)) + 
  geom_col(position = 'fill', col='black') + 
  scale_x_reordered() +
  facet_wrap(~switching, scales='free', drop=TRUE, nrow=3) +
  labs(title='High Complexity',
       x='Participants', 
       y='', 
       fill='Decision Strategy') +
  theme_bw() +
  scale_y_continuous(breaks = seq(0,1,length.out=3)) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank() ) 

lcm_plot <- plot_weights_simple + plot_weights_complex +
  plot_layout(guides = 'collect')


ggsave("figures/strategy_classification.eps",
       width = 17.62, height=9/0.7, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/strategy_classification.png",
       width = 17.62, height=9/0.7, units="cm",dpi=900)

