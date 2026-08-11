# preamble -------------------------------------------------------------
rm(list = ls())

# packages 
pacman::p_load(tidyverse, gtools , 
               rstan, brms, tidybayes, ggdist , # bayes 
               scico, patchwork, tidytext, ggbeeswarm # plotting
               )

# data
choices <- read_csv('data/choices.csv')
problems <- read_csv('data/problems.csv')
dat <- choices |> left_join(problems, by=join_by(complexity, problem, problem_unique))

## plotting ----------------------------------------------------------------

width <- 0.25
alpha <- .3
col_strat <- 'bamako' # bamako lapaz  imola
scico(n=3, palette=col_strat, begin=.1, end=.9)

# ggplot
three_col<- c('#DF9D56', '#6EB5E1', '#572948')
three_lvl <- c('infrequent', 'frequent', 'none')
three_lab <- c('infrequent', 'frequent', 'parallel')



# data set description -------------------------------------------------------------

dat |> 
  distinct(participant, .keep_all = T) |> 
  summarise(n = n(), 
            m = mean(age, na.rm = T), 
            sd = sd(age, na.rm = T), 
            min = min(age, na.rm =T),
            max = max(age, na.rm =T))

dat |> 
  distinct(participant, .keep_all = T) |> 
  group_by(gender) |> 
  summarise(n = n())

# behavioral analysis ----------------------------------------------------

## main --------------------------------------------------------------------

lm_dat <- list(PART = dat$participant ,
               PROB = dat$problem_unique , 
               S = factor(dat$switching) , 
               CPX = factor(dat$complexity) , 
               C = if_else(dat$choice==dat$short_option,1,0) # short-term choice
               )
  
lm_f <- bf(C ~ 1 + CPX*S + (1|PROB) + (1|PART))

lm <- brm(formula=lm_f , 
          data = lm_dat ,
          family = bernoulli(link = "logit") ,
          iter = 2000 ,
          warmup = 1000 ,
          chains = 6  ,
          cores = 6 ,
          save_pars = save_pars(all = TRUE),
          file = 'fits/linear_model')

# posterior
lm_effect_samples <- lm |> 
  spread_draws(b_Intercept , b_Sinfrequent , b_Snone, # complex problems
               b_CPXsimple , `b_CPXsimple:Sinfrequent`, `b_CPXsimple:Snone`) |>  # simple problems 
  mutate(m_C_F = plogis(b_Intercept) ,
         m_C_I = plogis(b_Intercept + b_Sinfrequent) , 
         m_C_P = plogis(b_Intercept + b_Snone) ,
         m_S_F = plogis(b_Intercept + b_CPXsimple) ,
         m_S_I = plogis(b_Intercept + b_CPXsimple + b_Sinfrequent + `b_CPXsimple:Sinfrequent`) , 
         m_S_P = plogis(b_Intercept + b_CPXsimple + b_Snone + `b_CPXsimple:Snone`), 
         d_C_FI = m_C_F - m_C_I ,
         d_C_FP = m_C_F - m_C_P , 
         d_C_IP = m_C_I - m_C_P ,
         d_S_FI = m_S_F - m_S_I ,
         d_S_FP = m_S_F - m_S_P , 
         d_S_IP = m_S_I - m_S_P , 
         d_F = m_C_F - m_S_F ,
         d_I = m_C_I - m_S_I ,
         d_P = m_C_P - m_S_P 
         ) 

lm_effect_summaries <- lm_effect_samples |> 
  gather_draws(m_S_F, m_S_I, m_S_P , 
               m_C_F, m_C_I, m_C_P , 
               d_S_FI, d_S_FP, d_S_IP, 
               d_C_FI, d_C_FP, d_C_IP ,
               d_F, d_I, d_P) |> 
  mean_qi()


# plot group means (posterior summaries)


# replace with random effects of subjects


choice_prop_ind <- dat |> 
  mutate(short_choice = if_else(choice==short_option,1,0) ,
         switching = if_else(switching=='none', "parallel" , switching)) |>
  group_by(participant, complexity, switching) |> 
  summarize(choice_prop = mean(short_choice))


figure_group_means <- lm_effect_summaries |> 
  filter(!str_detect(.variable, "d")) |> 
  mutate(complexity = case_when(str_detect(.variable, "_C_") ~ "complex" ,
                                str_detect(.variable, "_S_") ~ "simple" ,
                                TRUE ~ NA_character_) ,
         switching = case_when(str_detect(.variable, "F") ~ "frequent"  ,
                               str_detect(.variable, "I") ~ "infrequent" ,
                               str_detect(.variable, "P") ~ "parallel" ,
                               TRUE ~ NA_character_)) |> 
  ggplot(aes(y=.value , 
             x=factor(switching, levels=c("infrequent", "frequent", "parallel")), 
             color=factor(switching, levels=c("infrequent", "frequent", "parallel")))) +
  facet_wrap(~factor(complexity, levels=c("simple", "complex"), labels=c('Simple', 'Complex')), nrow=1) +
  geom_beeswarm(data = choice_prop_ind, aes(y=choice_prop), 
                alpha = .3, 
                color = 'lightgray', 
                cex=1.5, 
                size=2, 
                side = -1L) +
  geom_linerange(aes(ymin = .lower, ymax = .upper), position = position_dodge(width=.5), linewidth=1) +
  geom_point(position = position_dodge(width=.5), size = 3) + 
  geom_hline(yintercept = .5, linetype='dashed', color='gray') +
  scale_color_manual(values=three_col) +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Switching", 
       y = "Short-term Winner Choice Proportion",
       color = "Switching") +
  theme_bw()
figure_group_means

ggsave("manuscript/figures/group_means.eps", plot=figure_group_means,
       width = 17.62, height=9/0.9, units="cm",dpi=600, device = cairo_ps)
ggsave("manuscript/figures/group_means.jpg", plot=figure_group_means,
       width = 17.62, height=9/0.9, units="cm",dpi=900)


## risk preference (simple problems) ---------------------------------------------------------

simple_dat <- dat |> filter(complexity=='simple')

lm_risk_dat <- list(PART = simple_dat$participant ,
                    PROB = simple_dat$problem_unique , 
                    S = factor(simple_dat$switching) , 
                    R = factor(simple_dat$short_type) , 
                    C = if_else(simple_dat$choice==simple_dat$short_option,1,0))

lm_risk_f <- bf(C ~ 1 + R*S + (1|PROB) + (1|PART))

lm_risk <- brm(formula=lm_risk_f , 
               data = lm_risk_dat ,
               family = bernoulli(link = "logit") ,
               iter = 2000 ,
               warmup = 1000 ,
               chains = 6  ,
               cores = 6 ,
               save_pars = save_pars(all = TRUE) ,
               file = 'fits/linear_model_risk')

# posteriors

lm_risk_effect_samples <- lm_risk |> 
  spread_draws(b_Intercept , b_Sinfrequent , b_Snone, # short-term is risky
               b_Rsafe , `b_Rsafe:Sinfrequent`, `b_Rsafe:Snone`) |>  # short-term is safe 
  mutate(m_R_F = plogis(b_Intercept) ,
         m_R_I = plogis(b_Intercept + b_Sinfrequent) , 
         m_R_P = plogis(b_Intercept + b_Snone) ,
         m_S_F = plogis(b_Intercept + b_Rsafe) ,
         m_S_I = plogis(b_Intercept + b_Rsafe + b_Sinfrequent + `b_Rsafe:Sinfrequent`) , 
         m_S_P = plogis(b_Intercept + b_Rsafe + b_Snone + `b_Rsafe:Snone`), 
         d_R_FI = m_R_F - m_R_I ,
         d_R_FP = m_R_F - m_R_P , 
         d_R_IP = m_R_I - m_R_P ,
         d_S_FI = m_S_F - m_S_I ,
         d_S_FP = m_S_F - m_S_P , 
         d_S_IP = m_S_I - m_S_P , 
         d_F = m_R_F - m_S_F ,
         d_I = m_R_I - m_S_I ,
         d_P = m_R_P - m_S_P) 

lm_risk_effect_summaries <- lm_risk_effect_samples |> 
  gather_draws(m_S_F, m_S_I, m_S_P , 
               m_R_F, m_R_I, m_R_P , 
               d_S_FI, d_S_FP, d_S_IP, 
               d_R_FI, d_R_FP, d_R_IP ,
               d_F, d_I, d_P) |> 
  mean_qi()

# plot means for safe vs. risky short-term winning problems


choice_prop_ind_risk <- simple_dat |> 
  mutate(short_choice = if_else(choice==short_option,1,0) ,
         switching = if_else(switching=='none', "parallel" , switching)) |>
  group_by(participant, short_type, switching) |> 
  summarize(choice_prop = mean(short_choice))


figure_risk <- lm_risk_effect_summaries |> 
  filter(!str_detect(.variable, "d")) |> 
  mutate(short_type = case_when(str_detect(.variable, "_R_") ~ "risky" ,
                                str_detect(.variable, "_S_") ~ "safe" ,
                                TRUE ~ NA_character_) ,
         switching = case_when(str_detect(.variable, "F") ~ "frequent"  ,
                               str_detect(.variable, "I") ~ "infrequent" ,
                               str_detect(.variable, "P") ~ "parallel" ,
                               TRUE ~ NA_character_)) |> 
  ggplot(aes(x=factor(switching, levels=c("infrequent", "frequent", "parallel")), y=.value, 
             color=factor(switching, levels=c("infrequent", "frequent", "parallel")))) +
  facet_wrap(~factor(short_type, levels=c('risky', 'safe'), labels=c('Safe Long-Term Winner','Risky Long-Term Winner')), nrow=1) +
  geom_beeswarm(data = choice_prop_ind_risk, aes(y=choice_prop), 
                alpha = .3, 
                color = 'lightgray', 
                cex=1.5, 
                size=2, 
                side = -1L) +
  geom_linerange(aes(ymin = .lower, ymax = .upper), position = position_dodge(width=.5), linewidth=1) +
  geom_point(position = position_dodge(width=.5), size = 3) + 
  geom_hline(yintercept = .5, linetype='dashed', color='gray') +
  scale_color_manual(values=three_col) +
  scale_y_continuous(limits = c(0,1), breaks=seq(0,1,length.out=3)) +
  labs(x = "Switching", 
       y = "Short-term Winner Choice Proportion",
       color = "Switching",
       linetype = "Short-Term Option") +
  theme_bw()

ggsave("manuscript/figures/risk_preference.eps", plot=figure_risk,
       width = 17.62, height=9/0.9, units="cm",dpi=600, device = cairo_ps)
ggsave("manuscript/figures/risk_preference.jpg", plot=figure_risk,
       width = 17.62, height=9/0.9, units="cm",dpi=900)



# computational analysis --------------------------------------------------

current_model <- 'lcm'  
groups <- expand_grid(current_model=current_model , 
                      complexity=c('simple', 'complex') ,
                      switching=c('none', 'infrequent','frequent')) |> 
  mutate(file = paste0('fits/',current_model,'_',complexity, '_', switching, '.RData'))
current_model_summaries <- paste0('fits/summaries/',current_model,'.RData')

redo_analysis <- FALSE
if(redo_analysis){file.remove(groups$file,  current_model_summaries)}
missing_fits <- anyNA(if_else(file.exists(groups$file),TRUE,NA))


## fitting -----------------------------------------------------------------

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
                control=list(adapt_delta=.995) , 
                cores = n_chains, refresh = (n_iter/10) ,
                seed = 421886)
    save(lcm, file=stanfile)
  }
}

## post-processing ----------------------------------------------------------

if(!file.exists(current_model_summaries)){
  
  weights_list <- posteriors_list <- errors_list <- vector('list', length = nrow(groups))
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
             diff_sum_guess = `weights[2]` - `weights[3]`) |> 
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
    
    
    posteriors <- model_params |> 
      as.data.frame() |> 
      rownames_to_column("parameter") |> 
      filter(str_detect(parameter, "^class_prob\\[")) |> 
      separate_wider_regex(parameter,
                           patterns= c("^class_prob\\[", participant = "\\d+", '\\,', strategy='\\d+', "\\]")) |>
      mutate(complexity=cpx, 
             switching=switch) |> 
      select(complexity, switching, participant, strategy, pprob = mean, lower_CI = `2.5%`, upper_CI = `97.5%` )
    
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
    posteriors_list[[g]] <- posteriors
    errors_list[[g]] <- bind_rows(group_error, errors) 
    
  }
  
  all_weights <- weights_list |> bind_rows() 
  all_posteriors <- posteriors_list |> bind_rows()
  all_errors <- errors_list |> bind_rows() 
  # |> mutate(error_transformed = pnorm(error)) only required if error is on probit scale
  
  save(all_weights, all_posteriors, all_errors, file=current_model_summaries)
  
} else{load(current_model_summaries, verbose = T)}


## analysis ----------------------------------------------------------------

### strategy weights --------------------------------------------------------
lcm_weights <- all_weights |>
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
                fill = factor(strategy, 
                              levels = c('Roundwise', 'Summary', 'Guessing'))) , 
            alpha = .3) +
  geom_rect(aes(xmin = x + offset - width/2 ,
                xmax = x + offset + width/2 ,
                ymin = lower_CI_0.5 ,
                ymax = upper_CI_0.5 , 
                fill = factor(strategy, 
                              levels = c('Roundwise', 'Summary', 'Guessing'))) ,
            alpha = .3) +
  geom_segment(aes(x = x + offset - width/2 ,
                   xend = x + offset + width/2 ,
                   y = weight ,
                   yend = weight ,
                   col = factor(strategy, 
                                levels = c('Roundwise', 'Summary', 'Guessing'))) ,
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

# jpg for lcm_weights as eps does not support transparency

# ggsave("manuscript/figures/lcm_weights.eps", plot=lcm_weights,
#        width = 17.62, height=9/1.2, units="cm",dpi=600, device = cairo_ps)

ggsave("manuscript/figures/lcm_weights.jpg", plot=lcm_weights, 
       width = 17.62, height=9/1.2, units="cm",dpi=900)

### individual posterior probabilities ---------------------------------------------------

nsub <- all_posteriors |> 
  distinct(complexity, switching, participant) |> 
  nrow()

part_uni <- all_posteriors |> 
  distinct(complexity, switching, participant) |> 
  mutate(part_unique=row_number())

all_posts <- all_posteriors |>
  left_join(part_uni, by=join_by(complexity, switching, participant)) |> 
  select(complexity, switching, participant, part_unique, everything())

part_order <- all_posts |> 
  filter(strategy == 2) |> 
  arrange(desc(pprob)) |> 
  pull(part_unique)

lcm_pprobs <- all_posts |> 
  mutate(complexity = str_to_title(complexity), 
         switching= str_to_title(switching)) |> 
  group_by(complexity, switching) |> 
  mutate(part_unique = factor(part_unique, levels = part_order)) |> 
  ggplot(aes(x = part_unique, y = pprob, 
             fill = factor(strategy, 
                           levels = c('1','2','3'), 
                           labels = c('Roundwise', 'Summary', 'Guessing')))) +
  ggh4x::facet_grid2(factor(switching, levels=c('Infrequent', 'Frequent', 'None'))~factor(complexity, levels=c('Simple', 'Complex')),  
                     scales='free_x', independent = 'x') +
  geom_col(position = "fill") +
  labs(x = "Participant", 
       y = "Probability of Strategy Use",
       fill='Strategy') +
  theme_bw() +
  scale_fill_manual(values = c("#0C4240", "#61790A", "#E4CA74")) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

ggsave("manuscript/figures/lcm_pprobs.eps", plot=lcm_pprobs, 
       width = 17.62, height=17.62, units="cm",dpi=600, device = cairo_ps)
ggsave("manuscript/figures/lcm_ppreds.jpg", plot=lcm_pprobs, 
       width = 17.62, height=17.62, units="cm",dpi=900)

### errors ------------------------------------------------------------------

# estimate group level error and add uncertainty in the estimates
lcm_mu_error <- all_errors |> 
  filter(participant=='Group Mean') |> 
  ggplot(aes(x=switching, y=error, group=switching, colour = complexity)) + 
  geom_point(size=2, position = position_dodge2(width=0.5)) +
  geom_linerange(aes(ymin=lower_CI, ymax=upper_CI), position = position_dodge2(width=0.5)) + 
  geom_hline(yintercept = .5, linetype='dashed') + 
  scale_color_scico_d(palette='vanimo', begin=.1, end=.9) +
  scale_y_continuous(limits=c(0,.55), breaks=seq(0,.5,length.out=3)) +
  theme_bw() 

ggsave("manuscript/figures/lcm_mu_error.eps", plot=lcm_mu_error,
       width = 17.62, height=9/1.2, units="cm",dpi=600, device = cairo_ps)
ggsave("manuscript/figures/lcm_mu_error.jpg", plot=lcm_mu_error,
       width = 17.62, height=9/1.2, units="cm",dpi=900)
