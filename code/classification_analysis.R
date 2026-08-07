rm(list = ls())

# preamble ----------------------------------------------------------------

# packages
pacman::p_load(tidyverse, gtools, rstan)

# data
choices <- read_csv('data/choices.csv')
problems <- read_csv('data/problems.csv')

redo_analysis <- TRUE

# models <-  expand_grid(weights = c('pooled', 'hierarchical') , 
#                        errors = c('single', 'multiple') , 
#                        strategies = c('without guessing', 'with guessing'))


models <- c(paste0(rep('lcm_',4), 1:4), paste0(rep('mixture_',2), 1:2))
current_model <- models[6] 

current_model_summaries <- paste0('fits/summaries/',current_model,'.RData')
groups <- expand_grid(current_model=current_model , 
                      complexity=c('simple', 'complex') ,
                      switching=c('none', 'infrequent','frequent')) |> 
  mutate(file = paste0('fits/',current_model,'_',complexity, '_', switching, '.RData'))


if(redo_analysis){file.remove(groups$file,  current_model_summaries)}

missing_fits <- anyNA(if_else(file.exists(groups$file),TRUE,NA))

# fitting -------------------------------------------------------------

  
if(missing_fits){
  
  # Stan settings
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  lcm_stan <- paste0('code/models/',current_model,'.stan') 
  n_chains <- 6 # number of parallel chains
  n_iter <- 7500
  n_disc <- 2500
  
  for (g in seq_len(nrow(groups))){ 
    
    cpx <- groups[[g, 'complexity']]
    switch <- groups[[g, 'switching']]
    stanfile <- groups[[g, 'file']]
    
    # refit model
 
    obs <- choices |> 
      select(-problem_unique) |> 
      pivot_wider(names_from = problem, names_prefix = 'problem', values_from = choice) |> 
      filter(complexity==cpx & switching==switch) |> 
      select(problem1:problem20) |> 
      as.matrix()
    
    preds <- problems |>
      filter(complexity==cpx) |> 
      select(short=short_option, long=long_option) |>
      as.matrix()
    
    lcm_dat <- list(obs = obs ,
                    preds = preds ,
                    nsub = dim(obs)[1] , # number of subjects in group
                    nprob = dim(obs)[2] , # number of problems (20)
                    K_sys = dim(preds)[2]) # number of possible strategies
    
    lcm <- stan(file = lcm_stan, data = lcm_dat ,
                chains = n_chains, iter = n_iter, warmup = n_disc ,
                cores = n_chains, refresh = (n_iter/10) ,
                seed = 421886 # using varying seeds
                )
    save(lcm, file=stanfile)
  }
}



# Post-processing and analysis ---------------------------------------------

width <- 0.25
alpha <- .3
col_strat <- 'bamako' # bamako lapaz  tokyo  imola

## LCM1 ----------------------------------------------------------
' 
- 2 Pooled weights for roundwise / summary
- Subject-specific errors from partial pooling (hierarchy)
'

current_model <- models[1] 
current_model_summaries <- paste0('fits/summaries/',current_model,'.RData')

### post-processing ----------------------------------------------------------

if(!file.exists(current_model_summaries)){
  
  weights_list <- errors_list <- vector('list', length = nrow(groups))
  for (g in seq_len(nrow(groups))){ 
    
    cpx <- groups[[g, 'complexity']]
    switch <- groups[[g, 'switching']]
    stanfile <- groups[[g, 'file']]
    load(stanfile)
    
    # extract relevant parameters (strategy weights and parameters)
    model_params <- summary(lcm)$summary
    
    weights <- lcm |> 
      spread_draws(`weights[1]`, `weights[2]`) |>
      mutate(diff_round_sum = `weights[1]` - `weights[2]`) |> 
      gather_draws(`weights[1]`, `weights[2]`, diff_round_sum) |> 
      mean_qi(.width = c(.50,.80,.95))  |> 
      mutate(complexity=cpx, 
             switching=switch ,
             .variable=case_match(.variable, 
                                  'diff_round_sum' ~ 'Roundwise - Summary' ,
                                  'weights[1]' ~ 'Roundwise',
                                  'weights[2]' ~ 'Summary')) |> 
      select(complexity, switching, strategy=.variable, weight=.value, .width, lower_CI=.lower, upper_CI=.upper)
    
    group_error <- model_params |>
      as.data.frame() |> 
      rownames_to_column("parameter") |> 
      filter(str_detect(parameter, "^mu_error")) |>
      mutate(participant = 'Group Mean' ,
             complexity=cpx , 
             switching=switch) |> 
      select(complexity, switching, participant, error = mean, lower_CI = `2.5%`, upper_CI = `97.5%` )
    
    errors <- model_params |> 
      as.data.frame() |> 
      rownames_to_column("parameter") |> 
      filter(str_detect(parameter, "^errors\\[")) |> 
      separate_wider_regex(parameter,
                           patterns= c("^errors\\[", participant = "\\d+", "\\]")) |>
      mutate(complexity=cpx, 
             switching=switch) |> 
      select(complexity, switching, participant, error = mean, lower_CI = `2.5%`, upper_CI = `97.5%` )
    
    
    # make summary data frame
    weights_list[[g]] <- weights
    errors_list[[g]] <- bind_rows(group_error, errors) 
    
  }
  
  all_weights <- weights_list |> bind_rows() 
  all_errors <- errors_list |> bind_rows() 
  
  save(all_weights, all_errors, file=current_model_summaries)
  
}

### analysis ----------------------------------------------------------------

load(current_model_summaries, verbose = T)


#### strategy weights --------------------------------------------------------

lcm_1_weights <- all_weights |>
  filter(strategy %in% c("Roundwise", "Summary")) |>
  pivot_wider(names_from = .width, values_from = c(lower_CI, upper_CI)) |> 
  mutate(x = as.numeric(case_match(switching, 
                                   'infrequent' ~ 1,
                                   'frequent' ~ 2,
                                   'none' ~ 3)),
         offset = c(Roundwise = -width,
                    Summary    = 0,
                    Guessing   = width)[strategy]) |> 
  ggplot() +
  geom_rect(aes(xmin = x + offset - width/2 ,
                xmax = x + offset + width/2 ,
                ymin = lower_CI_0.95 ,
                ymax = upper_CI_0.95 ,
                fill = strategy) ,
            alpha = alpha) +
  geom_rect(aes(xmin = x + offset - width/2 ,
                xmax = x + offset + width/2 ,
                ymin = lower_CI_0.8 , 
                ymax = upper_CI_0.8 ,
                fill = strategy) ,
            alpha = alpha) +
  geom_rect(aes(xmin = x + offset - width/2 ,
                xmax = x + offset + width/2 ,
                ymin = lower_CI_0.5 ,
                ymax = upper_CI_0.5 ,
                fill = strategy) ,
            alpha = alpha) +
  geom_segment(aes(x = x + offset - width/2 ,
                   xend = x + offset + width/2 ,
                   y = weight ,
                   yend = weight ,
                   col = strategy) ,
               linewidth = 2) +
  scale_x_continuous(breaks = 1:3, labels=c('Infrequent', 'Frequent', 'Parallel')) +
  facet_wrap(~factor(complexity, levels = c('simple', 'complex'), 
                     labels=c('Simple', 'Complex')), drop=TRUE, nrow=1) +
  labs(x='Sample Organization', 
       y='Strategy Weight', 
       fill='Decision Strategy', 
       color = 'Decision Strategy') +
  theme_bw() +
  scale_y_continuous(limits=c(-.05, 1.05), breaks = seq(0,1,length.out=3)) +
  scale_fill_scico_d(palette=col_strat, begin=.1, end=.9) + 
  scale_color_scico_d(palette=col_strat, begin=.1, end=.9)

# posterior mean and 95% CI for differences in the strategy weights
# all_weights |> 
#   filter(strategy=='Difference') |> 
#   ggplot(aes(x=switching, y=weight, fill=factor(weight<=0, labels=c('Roundwise', 'Summary')))) + 
#   geom_col(position = 'dodge', col='black') + 
#   geom_linerange(aes(ymin=lower_CI, ymax=upper_CI, group=switching), position=position_dodge2(width = 0.9)) +
#   geom_hline(yintercept = 0, linetype='dashed') + 
#   facet_wrap(~factor(complexity, levels = c('simple', 'complex')), drop=TRUE, nrow=1) +
#   labs(x='Strategy', 
#        y='Difference in Weights', 
#        fill='Decision Strategy') +
#   theme_bw() +
#   scale_y_continuous(breaks = seq(0,1,length.out=3)) 


ggsave("figures/lcm_1_strategy_weights.eps", plot=lcm_1_weights, 
       width = 17.62, height=9/1.2, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/lcm_1_strategy_weights.png", plot=lcm_1_weights, 
       width = 17.62, height=9/1.2, units="cm",dpi=900)


#### errors ------------------------------------------------------------------


lcm_1_errors <- all_errors |>
  filter(participant=='Group Mean') |> 
  ggplot(aes(x=switching, y=error, group=switching, color=complexity)) + 
  geom_point(size=2, position=position_dodge2(width=.5)) +
  geom_linerange(aes(ymin=lower_CI, ymax=upper_CI), position=position_dodge2(width=.5)) + 
  geom_hline(yintercept = .5, linetype='dashed') + 
  scale_color_scico_d(palette='vanimo', begin=.1, end=.9) +
  scale_y_continuous(limits=c(0,.55), breaks=seq(0,.5,length.out=3)) +
  theme_bw() 



ggsave("figures/lcm_1_errors.eps", plot=lcm_1_errors, 
       width = 17.62, height=9/1.2, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/lcm_1_errors.png", plot=lcm_1_errors, 
       width = 17.62, height=9/1.2, units="cm",dpi=900)


## LCM 2 -------------------------------------------------------------------
' 
- 3 Pooled weights for roundwise / summary / guessing
- Subject-specific errors
'

current_model <- models[2] 
current_model_summaries <- paste0('fits/summaries/',current_model,'.RData')


### post-processing ----------------------------------------------------------

if(!file.exists(current_model_summaries)){
  
  weights_list <- errors_list <- vector('list', length = nrow(groups))
  for (g in seq_len(nrow(groups))){ 
    
    cpx <- groups[[g, 'complexity']]
    switch <- groups[[g, 'switching']]
    stanfile <- groups[[g, 'file']]
    load(stanfile)
    
    # extract relevant parameters (strategy weights and parameters)
    model_params <- summary(lcm)$summary
    
    weights <- lcm |> 
      spread_draws(`weights[1]`, `weights[2]`, `weights[3]`) |> 
      mutate(diff_round_sum = `weights[1]` - `weights[2]` , 
             diff_round_guess = `weights[1]` - `weights[3]` ,
             diff_sum_guess = `weights[2]` - `weights[3]` ,
             ) |> 
      gather_draws(`weights[1]`, `weights[2]`, `weights[3]`, 
                   diff_round_sum, diff_round_guess, diff_sum_guess) |> 
      mean_qi(.width = c(.50,.80,.95))  |> 
      mutate(complexity=cpx, 
             switching=switch ,
             .variable=case_match(.variable, 
                                  'diff_round_sum' ~ 'Roundwise - Summary' ,
                                  'diff_round_guess' ~ 'Roundwise - Guessing' ,
                                  'diff_sum_guess' ~ 'Summary - Guessing' ,
                                  'weights[1]' ~ 'Roundwise',
                                  'weights[2]' ~ 'Summary',
                                  'weights[3]' ~ 'Guessing')) |> 
      select(complexity, switching, strategy=.variable, weight=.value, .width, lower_CI=.lower, upper_CI=.upper)
    
      group_error <- model_params |>
        as.data.frame() |> 
        rownames_to_column("parameter") |> 
        filter(str_detect(parameter, "^mu_error")) |>
        mutate(participant = 'Group Mean' ,
               complexity=cpx , 
               switching=switch) |> 
        select(complexity, switching, participant, error = mean, lower_CI = `2.5%`, upper_CI = `97.5%` )
      
      errors <- model_params |> 
        as.data.frame() |> 
        rownames_to_column("parameter") |> 
        filter(str_detect(parameter, "^errors\\[")) |> 
        separate_wider_regex(parameter,
                             patterns= c("^errors\\[", participant = "\\d+", "\\]")) |>
        mutate(complexity=cpx, 
               switching=switch) |> 
        select(complexity, switching, participant, error = mean, lower_CI = `2.5%`, upper_CI = `97.5%` )
      
      
      # make summary data frame
      weights_list[[g]] <- weights
      errors_list[[g]] <- bind_rows(group_error, errors) 
         
  }
  
  all_weights <- weights_list |> bind_rows() 
  all_errors <- errors_list |> bind_rows() 
  # |> mutate(error_transformed = pnorm(error)) only required if error is on probit scale
  
  save(all_weights, all_errors, file=current_model_summaries)
  
}

### analysis ----------------------------------------------------------------

load(current_model_summaries, verbose = T)


#### strategy Weights --------------------------------------------------------

lcm_2_weights <- all_weights |>
  filter(strategy %in% c("Roundwise", "Summary", "Guessing")) |>
  pivot_wider(names_from = .width, values_from = c(lower_CI, upper_CI)) |> 
  mutate(x = as.numeric(case_match(switching, 
                                   'infrequent' ~ 1 ,
                                   'frequent' ~ 2 ,
                                   'none' ~ 3)) ,
         offset = c(Roundwise = -width ,
                    Summary    = 0 ,
                    Guessing   = width)[strategy]) |> 
  ggplot() +
  geom_rect(aes(xmin = x + offset - width/2 , 
                xmax = x + offset + width/2 , 
                ymin = lower_CI_0.95 , 
                ymax = upper_CI_0.95 ,
                fill = strategy) ,
            alpha = alpha) +
  geom_rect(aes(xmin = x + offset - width/2 ,
                xmax = x + offset + width/2 ,
                ymin = lower_CI_0.8 ,
                ymax = upper_CI_0.8 ,
                fill = strategy) , 
            alpha = .3) +
  geom_rect(aes(xmin = x + offset - width/2 ,
                xmax = x + offset + width/2 ,
                ymin = lower_CI_0.5 ,
                ymax = upper_CI_0.5 , 
                fill = strategy) ,
            alpha = .3) +
  geom_segment(aes(x = x + offset - width/2 ,
                   xend = x + offset + width/2 ,
                   y = weight ,
                   yend = weight ,
                   col = strategy) ,
               linewidth = 2) +
  scale_x_continuous(breaks = 1:3, labels=c('Infrequent', 'Frequent', 'Parallel')) +
  facet_wrap(~factor(complexity, levels = c('simple', 'complex'), 
                     labels=c('Simple', 'Complex')), drop=TRUE, nrow=1) +
  labs(x='Sample Organization', 
       y='Strategy Weight', 
       fill='Decision Strategy', 
       color = 'Decision Strategy') +
  theme_bw() +
  scale_y_continuous(limits=c(-.05, 1.05), breaks = seq(0,1,length.out=3)) +
  scale_fill_scico_d(palette=col_strat, begin=.1, end=.9) + 
  scale_color_scico_d(palette=col_strat, begin=.1, end=.9)


ggsave("figures/lcm_2_strategy_weights.eps", plot=lcm_2_weights, 
       width = 17.62, height=9/1.2, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/lcm_2_strategy_weights.png", plot=lcm_2_weights, 
       width = 17.62, height=9/1.2, units="cm",dpi=900)


#### errors ------------------------------------------------------------------

# estimate group level error and add uncertainty in the estimates
lcm_2_errors <- all_errors |> 
  filter(participant=='Group Mean') |> 
  ggplot(aes(x=switching, y=error, group=switching, colour = complexity)) + 
  geom_point(size=2, position = position_dodge2(width=0.5)) +
  geom_linerange(aes(ymin=lower_CI, ymax=upper_CI), position = position_dodge2(width=0.5)) + 
  geom_hline(yintercept = .5, linetype='dashed') + 
  scale_color_scico_d(palette='vanimo', begin=.1, end=.9) +
  scale_y_continuous(limits=c(0,.55), breaks=seq(0,.5,length.out=3)) +
  theme_bw() 

ggsave("figures/lcm_2_error.eps", plot=lcm_2_errors,
       width = 17.62, height=9/1.2, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/lcm_2_error.png", plot=lcm_2_errors,
       width = 17.62, height=9/1.2, units="cm",dpi=900)


## LCM 3 -------------------------------------------------------------------
' 
- 3 Pooled weights for roundwise / summary / guessing
- Subject-specific errors
'

current_model <- models[3] 
current_model_summaries <- paste0('fits/summaries/',current_model,'.RData')



### post-processing ----------------------------------------------------------

if(!file.exists(current_model_summaries)){
  
  weights_list <- group_errors_list <-  errors_list <- vector('list', length = nrow(groups))
  for (g in seq_len(nrow(groups))){ 
    
    cpx <- groups[[g, 'complexity']]
    switch <- groups[[g, 'switching']]
    stanfile <- groups[[g, 'file']]
    load(stanfile)
    
    # extract relevant parameters (strategy weights and parameters)
    model_params <- summary(lcm)$summary
    
    weights <- lcm |> 
      spread_draws(`weights[1]`, `weights[2]`, `weights[3]`) |> 
      mutate(diff_round_sum = `weights[1]` - `weights[2]` , 
             diff_round_guess = `weights[1]` - `weights[3]` ,
             diff_sum_guess = `weights[2]` - `weights[3]` ) |> 
      gather_draws(`weights[1]`, `weights[2]`, `weights[3]`, 
                   diff_round_sum, diff_round_guess, diff_sum_guess) |> 
      mean_qi(.width = c(.50,.80,.95))  |> 
      mutate(complexity=cpx, 
             switching=switch ,
             .variable=case_match(.variable, 
                                  'diff_round_sum' ~ 'Roundwise - Summary' ,
                                  'diff_round_guess' ~ 'Roundwise - Guessing' ,
                                  'diff_sum_guess' ~ 'Summary - Guessing' ,
                                  'weights[1]' ~ 'Roundwise',
                                  'weights[2]' ~ 'Summary',
                                  'weights[3]' ~ 'Guessing')) |> 
      select(complexity, switching, strategy=.variable, weight=.value, .width, lower_CI=.lower, upper_CI=.upper)
    
    
    group_errors <- lcm |>
      spread_draws(`mu_error[1]`, `mu_error[2]`) |> 
      mutate(diff_round_sum = `mu_error[1]` - `mu_error[2]`) |> 
      gather_draws(`mu_error[1]`, `mu_error[2]`, diff_round_sum) |>
      mean_qi(.width = c(.50,.80,.95))  |> 
      mutate(complexity=cpx, 
             switching=switch ,
             .variable=case_match(.variable, 
                                  'diff_round_sum' ~ 'Roundwise - Summary',
                                  'mu_error[1]' ~ 'Roundwise',
                                  'mu_error[2]' ~ 'Summary')) |> 
      select(complexity, switching, strategy=.variable, weight=.value, .width, lower_CI=.lower, upper_CI=.upper)
      
    errors <- model_params |> 
      as.data.frame() |> 
      rownames_to_column("parameter") |> 
      filter(str_detect(parameter, "^errors\\[")) |> 
      separate_wider_regex(parameter,
                           patterns= c("^errors\\[", participant = "\\d+", '\\,', strategy="\\d+", "\\]")) |>
      mutate(complexity=cpx, 
             switching=switch) |> 
      select(complexity, switching, participant, strategy, error = mean, lower_CI = `2.5%`, upper_CI = `97.5%` )
    
    
    # make summary data frame
    weights_list[[g]] <- weights
    group_errors_list[[g]] <- group_errors
    errors_list[[g]] <- errors
    
  }
  
  all_weights <- weights_list |> bind_rows() 
  all_group_errors <- group_errors_list |> bind_rows()
  all_errors <- errors_list |> bind_rows() 
  
  save(all_weights, all_group_errors, all_errors, file=current_model_summaries)
  
}

### analysis ----------------------------------------------------------------

load(current_model_summaries, verbose = T)


#### strategy Weights --------------------------------------------------------

lcm_3_weights <- all_weights |>
  filter(strategy %in% c("Roundwise", "Summary", "Guessing")) |>
  pivot_wider(names_from = .width, values_from = c(lower_CI, upper_CI)) |> 
  mutate(x = as.numeric(case_match(switching, 
                                   'infrequent' ~ 1 ,
                                   'frequent' ~ 2 ,
                                   'none' ~ 3)) ,
         offset = c(Roundwise = -width ,
                    Summary    = 0 ,
                    Guessing   = width)[strategy]) |> 
  ggplot() +
  geom_rect(aes(xmin = x + offset - width/2 , 
                xmax = x + offset + width/2 , 
                ymin = lower_CI_0.95 , 
                ymax = upper_CI_0.95 ,
                fill = strategy) ,
            alpha = alpha) +
  geom_rect(aes(xmin = x + offset - width/2 ,
                xmax = x + offset + width/2 ,
                ymin = lower_CI_0.8 ,
                ymax = upper_CI_0.8 ,
                fill = strategy) , 
            alpha = .3) +
  geom_rect(aes(xmin = x + offset - width/2 ,
                xmax = x + offset + width/2 ,
                ymin = lower_CI_0.5 ,
                ymax = upper_CI_0.5 , 
                fill = strategy) ,
            alpha = .3) +
  geom_segment(aes(x = x + offset - width/2 ,
                   xend = x + offset + width/2 ,
                   y = weight ,
                   yend = weight ,
                   col = strategy) ,
               linewidth = 2) +
  scale_x_continuous(breaks = 1:3, labels=c('Infrequent', 'Frequent', 'Parallel')) +
  facet_wrap(~factor(complexity, levels = c('simple', 'complex'), 
                     labels=c('Simple', 'Complex')), drop=TRUE, nrow=1) +
  labs(x='Sample Organization', 
       y='Strategy Weight', 
       fill='Decision Strategy', 
       color = 'Decision Strategy') +
  theme_bw() +
  scale_y_continuous(limits=c(-.05, 1.05), breaks = seq(0,1,length.out=3)) +
  scale_fill_scico_d(palette=col_strat, begin=.1, end=.9) + 
  scale_color_scico_d(palette=col_strat, begin=.1, end=.9)


ggsave("figures/lcm_3_strategy_weights.eps", plot=lcm_3_weights, 
       width = 17.62, height=9/1.2, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/lcm_3_strategy_weights.png", plot=lcm_3_weights, 
       width = 17.62, height=9/1.2, units="cm",dpi=900)


#### errors ------------------------------------------------------------------

# estimate group level error and add uncertainity in the estimates
all_group_errors |> 
  filter(strategy %in% c('Roundwise', 'Summary') & .width==.95) |>
  ggplot(aes(x=switching, y=weight, colour=strategy)) + 
  geom_point(size=3, position=position_dodge2(width = .5)) +
  geom_linerange(aes(ymin=lower_CI, ymax=upper_CI), position=position_dodge2(width = .5), linewidth=1) + 
  geom_hline(yintercept = .5, linetype='dashed') + 
  facet_wrap(~factor(complexity, levels = c('simple', 'complex')),  drop=TRUE, nrow=1) +
  scale_color_scico_d(palette='bamako') +
  scale_y_continuous(limits=c(0,.55), breaks=seq(0,.5,length.out=3)) +
  theme_bw() 

lcm_3_errors <- all_group_errors |> 
  filter(strategy %in% c('Roundwise', 'Summary') & .width==.95) |>
  ggplot(aes(x=switching, y=weight, linetype=strategy, colour=complexity)) + 
  geom_point(size=3, position=position_dodge2(width = .5)) +
  geom_linerange(aes(ymin=lower_CI, ymax=upper_CI), position=position_dodge2(width = .5), linewidth=1) + 
  geom_hline(yintercept = .5, linetype='dashed') + 
  #facet_wrap(~strategy, nrow=1) +
  scale_color_scico_d(palette='vanimo', begin=.1, end=.9) +
  scale_y_continuous(limits=c(0,.55), breaks=seq(0,.5,length.out=3)) +
  theme_bw() 


ggsave("figures/lcm_3_errors.eps",
       width = 17.62, height=9/1.2, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/lcm_3_errors.png",
       width = 17.62, height=9/1.2, units="cm",dpi=900)


## Mixture model -----------------------------------------------------------
' 
- 3 Pooled weights for roundwise / summary / guessing
- Subject-specific errors
- Each problem is predicted by a mixture of the strategies (subjects can rely on more than one strategy)
'

current_model <- models[4] 
current_model_summaries <- paste0('fits/summaries/',current_model,'.RData')


### post-processing ----------------------------------------------------------

if(!file.exists(current_model_summaries)){
  
  weights_list <- errors_list <- vector('list', length = nrow(groups))
  for (g in seq_len(nrow(groups))){ 
    
    cpx <- groups[[g, 'complexity']]
    switch <- groups[[g, 'switching']]
    stanfile <- groups[[g, 'file']]
    load(stanfile)
    
    # extract relevant parameters (strategy weights and parameters)
    model_params <- summary(lcm)$summary
    
    weights <- lcm |> 
      spread_draws(`weights[1]`, `weights[2]`, `weights[3]`) |> 
      mutate(diff_round_sum = `weights[1]` - `weights[2]` , 
             diff_round_guess = `weights[1]` - `weights[3]` ,
             diff_sum_guess = `weights[2]` - `weights[3]` ,
      ) |> 
      gather_draws(`weights[1]`, `weights[2]`, `weights[3]`, 
                   diff_round_sum, diff_round_guess, diff_sum_guess) |> 
      mean_qi(.width = c(.50,.80,.95))  |> 
      mutate(complexity=cpx, 
             switching=switch ,
             .variable=case_match(.variable, 
                                  'diff_round_sum' ~ 'Roundwise - Summary' ,
                                  'diff_round_guess' ~ 'Roundwise - Guessing' ,
                                  'diff_sum_guess' ~ 'Summary - Guessing' ,
                                  'weights[1]' ~ 'Roundwise',
                                  'weights[2]' ~ 'Summary',
                                  'weights[3]' ~ 'Guessing')) |> 
      select(complexity, switching, strategy=.variable, weight=.value, .width, lower_CI=.lower, upper_CI=.upper)
    
    group_error <- model_params |>
      as.data.frame() |> 
      rownames_to_column("parameter") |> 
      filter(str_detect(parameter, "^mu_error")) |>
      mutate(participant = 'Group Mean' ,
             complexity=cpx , 
             switching=switch) |> 
      select(complexity, switching, participant, error = mean, lower_CI = `2.5%`, upper_CI = `97.5%` )
    
    errors <- model_params |> 
      as.data.frame() |> 
      rownames_to_column("parameter") |> 
      filter(str_detect(parameter, "^errors\\[")) |> 
      separate_wider_regex(parameter,
                           patterns= c("^errors\\[", participant = "\\d+", "\\]")) |>
      mutate(complexity=cpx, 
             switching=switch) |> 
      select(complexity, switching, participant, error = mean, lower_CI = `2.5%`, upper_CI = `97.5%` )
    
    
    # make summary data frame
    weights_list[[g]] <- weights
    errors_list[[g]] <- bind_rows(group_error, errors) 
    
  }
  
  all_weights <- weights_list |> bind_rows() 
  all_errors <- errors_list |> bind_rows() 
  # |> mutate(error_transformed = pnorm(error)) only required if error is on probit scale
  
  save(all_weights, all_errors, file=current_model_summaries)
  
}

### analysis ----------------------------------------------------------------

load(current_model_summaries, verbose = T)


#### strategy Weights --------------------------------------------------------

mixture_1_weights <- all_weights |>
  filter(strategy %in% c("Roundwise", "Summary", "Guessing")) |>
  pivot_wider(names_from = .width, values_from = c(lower_CI, upper_CI)) |> 
  mutate(x = as.numeric(case_match(switching, 
                                   'infrequent' ~ 1 ,
                                   'frequent' ~ 2 ,
                                   'none' ~ 3)) ,
         offset = c(Roundwise = -width ,
                    Summary    = 0 ,
                    Guessing   = width)[strategy]) |> 
  ggplot() +
  geom_rect(aes(xmin = x + offset - width/2 , 
                xmax = x + offset + width/2 , 
                ymin = lower_CI_0.95 , 
                ymax = upper_CI_0.95 ,
                fill = strategy) ,
            alpha = alpha) +
  geom_rect(aes(xmin = x + offset - width/2 ,
                xmax = x + offset + width/2 ,
                ymin = lower_CI_0.8 ,
                ymax = upper_CI_0.8 ,
                fill = strategy) , 
            alpha = .3) +
  geom_rect(aes(xmin = x + offset - width/2 ,
                xmax = x + offset + width/2 ,
                ymin = lower_CI_0.5 ,
                ymax = upper_CI_0.5 , 
                fill = strategy) ,
            alpha = .3) +
  geom_segment(aes(x = x + offset - width/2 ,
                   xend = x + offset + width/2 ,
                   y = weight ,
                   yend = weight ,
                   col = strategy) ,
               linewidth = 2) +
  scale_x_continuous(breaks = 1:3, labels=c('Infrequent', 'Frequent', 'Parallel')) +
  facet_wrap(~factor(complexity, levels = c('simple', 'complex'), 
                     labels=c('Simple', 'Complex')), drop=TRUE, nrow=1) +
  labs(x='Sample Organization', 
       y='Strategy Weight', 
       fill='Decision Strategy', 
       color = 'Decision Strategy') +
  theme_bw() +
  scale_y_continuous(limits=c(-.05, 1.05), breaks = seq(0,1,length.out=3)) +
  scale_fill_scico_d(palette=col_strat, begin=.1, end=.9) + 
  scale_color_scico_d(palette=col_strat, begin=.1, end=.9)

mixture_1_weights


ggsave("figures/mixture_1_strategy_weights.eps", plot=mixture_1_weights, 
       width = 17.62, height=9/1.2, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/mixture_1_strategy_weights.png", plot=mixture_1_weights, 
       width = 17.62, height=9/1.2, units="cm",dpi=900)


#### errors ------------------------------------------------------------------

# estimate group level error and add uncertainty in the estimates
mixture_1_errors <- all_errors |> 
  filter(participant=='Group Mean') |> 
  ggplot(aes(x=switching, y=error, group=switching, colour = complexity)) + 
  geom_point(size=2, position = position_dodge2(width=0.5)) +
  geom_linerange(aes(ymin=lower_CI, ymax=upper_CI), position = position_dodge2(width=0.5)) + 
  geom_hline(yintercept = .5, linetype='dashed') + 
  scale_color_scico_d(palette='vanimo', begin=.1, end=.9) +
  scale_y_continuous(limits=c(0,.55), breaks=seq(0,.5,length.out=3)) +
  theme_bw() 

ggsave("figures/mixture_1_errors.eps", plot=mixture_1_errors,
       width = 17.62, height=9/1.2, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/mixture_1_errors.png", plot=mixture_1_errors,
       width = 17.62, height=9/1.2, units="cm",dpi=900)



















# other stuff -------------------------------------------------------------


## LCM 4: Hierarchical Weights ---------------------------------------------


if(current_model=='hierarchical_lcm'){
  
  # strategy 1 is roundwise
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
  
}
