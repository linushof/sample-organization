# preamble -------------------------------------------------------------
rm(list = ls())

# packages 
pacman::p_load(tidyverse, gtools , 
               rstan, brms, tidybayes, ggdist , # bayes 
               scico, patchwork, tidytext, ggbeeswarm, # plotting
               knitr, kableExtra)

# data
choices <- read_csv('data/clean/choices.csv')
problems <- read_csv('data/clean/sample_sequences_summary.csv')
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


## tables ------------------------------------------------------------------

col_names <- c("Coef." , 
               "Mean", "2.5\\%", "97.5\\%", "Median", "SD" , 
               "$\\hat{R}$" , 
               "$\\text{ESS}_{\\text{bulk}}$" ,
               "$\\text{ESS}_{\\text{tail}}$")


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

## fitting --------------------------------------------------------------------

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


## tables -------------------------------------------------------------------

lower <- .025
upper <- .975

lm_draws <- lm |> 
  spread_draws(b_Intercept , b_Sinfrequent , b_Snone, # complex problems
               b_CPXsimple , `b_CPXsimple:Sinfrequent`, `b_CPXsimple:Snone`, # simple problems
               sd_PROB__Intercept, sd_PART__Intercept) |> # random effects
  # estimate conditional (subject-specific) means of the different conditions
  mutate(m_C_F = plogis(b_Intercept) , 
         m_C_I = plogis(b_Intercept + b_Sinfrequent) , 
         m_C_P = plogis(b_Intercept + b_Snone) ,
         m_S_F = plogis(b_Intercept + b_CPXsimple) ,
         m_S_I = plogis(b_Intercept + b_CPXsimple + b_Sinfrequent + `b_CPXsimple:Sinfrequent`) , 
         m_S_P = plogis(b_Intercept + b_CPXsimple + b_Snone + `b_CPXsimple:Snone`)) 

# group means and model coefficients

effects <- lm_draws |> 
  select(.chain, .iteration, .draw, 
         m_S_I, m_S_F, m_S_P,
         m_C_I, m_C_F, m_C_P,
         everything()) |> 
  summarise_draws('mean', 
                  ~quantile(.x, probs = lower) ,
                  ~quantile(.x, probs = upper) ,
                  'median', 'sd', 'rhat', 'ess_bulk', 'ess_tail') |> 
  rename(Coefficient = variable , 
         Mean = mean ,
         Median = median , 
         SD = sd ,
         R = rhat ,
         ESS_bulk = ess_bulk , 
         ESS_tail = ess_tail)
lm_effect_summaries <- effects

colnames(effects) <- col_names
  
effects$Coef. <- c(rep(c("Infrequent", "Frequent", "Parallel") , 2), 
                   "$\\beta_0$", "$\\beta_1$", "$\\beta_2$", "$\\beta_3$" , 
                   "$\\beta_4$", "$\\beta_5$", 
                   "$\\sigma_u$", "$\\sigma_v$")

effects |>
  kbl(format = "latex",
      booktabs = TRUE,
      caption = "Posterior Summaries of the Logistic Regression Model",
      label = 'lm',
      align = c("l", rep("r", 8)),
      escape = FALSE, 
      digits=3) |>
  pack_rows("Group Means", 1, 6, bold=T, italic=F) |>
  pack_rows("Simple Problems", 1, 3, bold=F, italic=T) |>
  pack_rows("Complex Problems", 4, 7, bold=F, italic=T) |>
  pack_rows("Model Coefficients", 7, 14, bold=T, italic=F) |>
  pack_rows("Fixed effects", 7, 12, bold=F, italic=T) |>
  pack_rows("Random effects (Hyperparameters)", 13, 14, bold=F, italic=T) |>
  footnote(general = "Group means refer to draws from the expected value of the 
           model's posterior predictive distribution when setting the
           problem- and participant-specific random effects to 0.",
           general_title = "Note. ",
           escape = FALSE,
           threeparttable = TRUE) |>  
  save_kable(paste0("manuscript/tables/lm.tex")) # Table C1 (dissertation)


# sampling contrasts

sampling_contrasts <- lm_draws |> 
  mutate(d_S_FI = m_S_F - m_S_I ,
         d_S_PI = m_S_P - m_S_I , 
         d_S_PF = m_S_P - m_S_F , 
         d_C_FI = m_C_F - m_C_I ,
         d_C_PI = m_C_P - m_C_I ,
         d_C_PF = m_C_P - m_C_F) |> 
  select(.chain, .iteration, .draw, 
         d_S_FI:d_C_PF) |> 
  summarise_draws('mean', 
                ~quantile(.x, probs = lower) ,
                ~quantile(.x, probs = upper) ,
                'median', 'sd') |> 
  rename(Coefficient = variable , 
         Mean = mean ,
         Median = median , 
         SD = sd)

colnames(sampling_contrasts) <- c("Groups", "Mean", "2.5\\%", "97.5\\%", "Median", "SD")  

sampling_contrasts$Groups <- rep( c('Frequent - Infrequent', 
                            'Parallel - Infrequent' , 
                            'Parallel - Frequent'), 2)

sampling_contrasts |>
  kbl(format = "latex",
      booktabs = TRUE,
      caption = "Pairwise Differences in the Expected Short-Term Winner Choice Probability Between Sampling Organizations",
      label = 'lm_sampling_contrasts',
      align = c("l", rep("r", 8)),
      escape = FALSE, 
      digits=3) |>
  pack_rows("Simple Problems", 1, 3, bold=F, italic=T) |>
  pack_rows("Complex Problems", 4, 6, bold=F, italic=T) |>
  footnote(general = "Values refer to differences between draws from the expected value of 
           the posterior predicive distributions of the two indicated sampling conditions
           (differences in the reported in Table \\\\ref{tab:lm}).",
           general_title = "Note. ",
           escape = FALSE,
           threeparttable = TRUE) |>  
  save_kable(paste0("manuscript/tables/lm_sampling_contrasts.tex")) # table 4.1 dissertation


# complexity contrasts

complexity_contrasts <- lm_draws |> 
  mutate(d_S_FI = m_S_F - m_S_I ,
         d_S_PI = m_S_P - m_S_I , 
         d_S_PF = m_S_P - m_S_F , 
         d_C_FI = m_C_F - m_C_I ,
         d_C_PI = m_C_P - m_C_I ,
         d_C_PF = m_C_P - m_C_F , 
         d_I = m_C_I - m_S_I ,
         d_F = m_C_F - m_S_F ,
         d_P = m_C_P - m_S_P , 
         d_FI = d_C_FI - d_S_FI , 
         d_PI = d_C_PI - d_S_PI ,
         d_PF = d_C_PF - d_S_PF) |> 
  select(.chain, .iteration, .draw, 
         d_I:d_PF) |> 
  summarise_draws('mean', 
                  ~quantile(.x, probs = lower) ,
                  ~quantile(.x, probs = upper) ,
                  'median', 'sd') |> 
  rename(Coefficient = variable , 
         Mean = mean ,
         Median = median , 
         SD = sd)

colnames(complexity_contrasts) <- c("Group(s)", "Mean", "2.5\\%", "97.5\\%", "Median", "SD")  

complexity_contrasts$`Group(s)` <- c('Infrequent', 'Frequent', 'Parallel', 
                          'Frequent - Infrequent',
                          'Parallel - Infrequent' , 
                          'Parallel - Frequent')

complexity_contrasts |>
  kbl(format = "latex",
      booktabs = TRUE,
      caption = "Pairwise Differences in the Expected Short-Term Winner Choice Probability Between Complexity Conditions",
      label = 'lm_complexity_contrasts',
      align = c("l", rep("r", 8)),
      escape = FALSE, 
      digits=3) |>
  pack_rows("Main effect", 1, 3, bold=F, italic=T) |>
  pack_rows("Interaction effect", 4, 6, bold=F, italic=T) |>
  footnote(general = "\\\\textit{Main effects} refer to differences in the draws from the expected value of 
           each sampling organization's posterior predicive distribution between simple vs. complex problems.
           \\\\textit{Interaction effects} refer to differences in the sampling contrasts 
           (see Table \\\\ref{tab:lm_sampling_contrasts}) between simple vs. complex problems: 
           Positive values indicate a larger between sampling organizations in the complex problems.
           Positive values were expected for all sampling contrasts.",
           general_title = "Note. ",
           escape = FALSE,
           threeparttable = TRUE) |>  
  save_kable(paste0("manuscript/tables/lm_complexity_contrasts.tex"))  # table 4.2 dissertation



## figure -----------------------------------------------------------------

# plot group means (posterior summaries)

choice_prop_ind <- dat |> 
  mutate(short_choice = if_else(choice==short_option,1,0) ,
         complexity = if_else(complexity=='simple', 'Simple', 'Complex'), 
         switching = case_match(switching, 
                                'none' ~ "Parallel" ,
                                'frequent' ~ 'Frequent', 
                                'infrequent' ~ 'Infrequent')) |>
  group_by(participant, complexity, switching) |> 
  summarize(choice_prop = mean(short_choice))


figure_group_means <- lm_effect_summaries |> 
  filter(str_detect(Coefficient, "m_")) |> 
  mutate(complexity = case_when(str_detect(Coefficient, "_C_") ~ "Complex" ,
                                str_detect(Coefficient, "_S_") ~ "Simple" ,
                                TRUE ~ NA_character_) ,
         switching = case_when(str_detect(Coefficient, "F") ~ "Frequent"  ,
                               str_detect(Coefficient, "I") ~ "Infrequent" ,
                               str_detect(Coefficient, "P") ~ "Parallel" ,
                               TRUE ~ NA_character_)) |> 
  ggplot(aes(y=Mean , 
             x=factor(switching, levels=c("Infrequent", "Frequent", "Parallel")), 
             color=factor(switching, levels=c("Infrequent", "Frequent", "Parallel")))) +
  facet_wrap(~factor(complexity, levels=c("Simple", "Complex")), nrow=1) +
  geom_beeswarm(data = choice_prop_ind, aes(y=choice_prop), 
                alpha = .3, 
                color = 'lightgray', 
                cex=1.5, 
                size=2, 
                side = -1L) +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), position = position_dodge(width=.5), linewidth=1) +
  geom_point(position = position_dodge(width=.5), size = 3) + 
  geom_hline(yintercept = .5, linetype='dashed', color='gray') +
  scale_color_manual(values=three_col) +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Sampling", 
       y = "Short-term Winner Choice Proportion",
       color = "Sampling") +
  theme_bw()

ggsave("manuscript/figures/group_means.eps", plot=figure_group_means,
       width = 17.62, height=9/1.1, units="cm",dpi=600, device = cairo_ps)
ggsave("manuscript/figures/group_means.jpg", plot=figure_group_means,
       width = 17.62, height=9/1.1, units="cm",dpi=900)
ggsave("manuscript/figures/group_means.pdf", plot=figure_group_means,
       width = 17.62, height=9/1.1, units="cm",dpi=900)


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



## tables ------------------------------------------------------------------

digits <- 2

error_tbl <- all_errors |> 
  filter(participant=='Group Mean') |> 
  mutate(compact = paste0(round(error,digits), ' [', round(lower_CI, digits), ',', round(upper_CI,digits),']')) |> 
  select(-c(participant, error, lower_CI, upper_CI)) |> 
  rename(Error=compact)

weight_tbl <- all_weights |> 
  filter(.width == 0.95 & strategy %in% c('Roundwise', 'Summary', 'Guessing')) |> 
  mutate(compact = paste0(round(weight,digits), ' [', round(lower_CI, digits), ',', round(upper_CI,digits),']')) |> 
  select(-c(.width, weight, lower_CI, upper_CI)) |> 
  pivot_wider(names_from = strategy, values_from = compact)

lcm_tbl <- weight_tbl |> 
  left_join(error_tbl, by=join_by(complexity, switching)) |> 
  rename(Complexity=complexity, 
         Sampling=switching) |> 
  mutate(Complexity = case_match(Complexity,
                                 'simple' ~ 'Simple', 
                                 'complex' ~ 'Complex'), 
         Sampling = case_match(Sampling,
                               'none' ~ 'Parallel', 
                               'infrequent' ~ 'Infrequent', 
                               'frequent' ~ 'Frequent')) |> 
  arrange(factor(Complexity, levels = c('Simple', 'Complex')), 
          factor(Sampling, levels = c('Infrequent', 'Frequent', 'Parallel'))) |> 
  select(-Complexity)

lcm_tbl |>
  kbl(format = "latex",
      booktabs = TRUE,
      caption = "Posterior Summaries of the Latent Classification Models",
      label = 'lcm',
      align = c(rep("l",2), rep("r", 4)),
      escape = FALSE, 
      digits=digits) |>
  pack_rows("Simple", 1, 3, bold=F, italic=T) |>
  pack_rows("Complex", 4, 6, bold=F, italic=T) |>
  footnote(general = "Values in squared brackets are the lower and upper boundary of the 95\\\\% credibile interval.",
           general_title = "Note. ",
           escape = FALSE,
           threeparttable = TRUE) |>  
  save_kable(paste0("manuscript/tables/lcm.tex")) # table 4.3 dissertation


## figure(s)  ----------------------------------------------------------------

### collection --------------------------------------------------------------

#### individual posterior probabilities ---------------------------------------------------

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
ggsave("manuscript/figures/lcm_ppreds.pdf", plot=lcm_pprobs, 
       width = 17.62, height=17.62, units="cm",dpi=900)

#### errors ------------------------------------------------------------------

# estimate group level error and add uncertainty in the estimates
lcm_mu_error <- all_errors |> 
  mutate(complexity=if_else(complexity=='simple','Simple', 'Complex'), 
         switching=case_match(switching, 
                              'none'~'Parallel',
                              'frequent'~'Frequent',
                              'infrequent'~'Infrequent')) |> 
  filter(participant=='Group Mean') |> 
  ggplot(aes(x=switching, y=error, color = complexity, group=switching)) + 
  geom_point(size=2, position = position_dodge2(width=0.5)) +
  geom_linerange(aes(ymin=lower_CI, ymax=upper_CI), position = position_dodge2(width=0.5)) + 
  geom_hline(yintercept = .5, linetype='dashed') + 
  scale_color_scico_d(palette='vanimo', begin=.1, end=.9) +
  scale_y_continuous(limits=c(0,.55), breaks=seq(0,.5,length.out=3)) +
  labs(x='Sampling', 
       y='Error', 
       color='Complexity') +
  theme_bw() 

ggsave("manuscript/figures/lcm_mu_error.eps", plot=lcm_mu_error,
       width = 17.62, height=9/1.2, units="cm",dpi=600, device = cairo_ps)
ggsave("manuscript/figures/lcm_mu_error.jpg", plot=lcm_mu_error,
       width = 17.62, height=9/1.2, units="cm",dpi=900)


#### guessing and error -------------------------------------------------------------

guess <- all_posts |> 
  filter(strategy==3) |> 
  select(!c(part_unique, lower_CI, upper_CI))

corr_dat <- all_errors |> 
  filter(!participant=='Group Mean') |>
  select(!c(lower_CI, upper_CI)) |> 
  left_join(guess, by=join_by(complexity, switching, participant))

lcm_guess_err_corr <- corr_dat |> 
  ggplot(aes(x=pprob, y=error)) + 
  #facet_grid(complexity~switching) +
  geom_point() +
  geom_smooth(method='gam', formula = y ~ s(x, bs = "cs"), color='#E4CA74', fill='#E4CA74') +
  labs(x='p(Guessing)',
       y='Error') +
  theme_bw()


#### merge ---------------------------------------------------------------------

lcm_collection <- lcm_pprobs + lcm_mu_error + lcm_guess_err_corr + guide_area() + 
  plot_layout(design=c(area(1,1,3), 
                       area(1,2),
                       area(2,2),
                       area(3,2)),
              widths=c(2,1),
              guides = 'collect') +
  plot_annotation(tag_levels = 'A') &
  theme(legend.box = 'horizontal')


ggsave("manuscript/figures/lcm_collection.pdf", plot=lcm_collection,
       width = 190, height= 190*.8, units="mm",dpi=900)

### strategy weights  --------------------------------------------------------
# not included in dissertation

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
  labs(x='Sampling', 
       y='Strategy Weight', 
       fill='Decision Strategy', 
       color = 'Decision Strategy') +
  theme_bw() +
  scale_y_continuous(limits=c(-.05, 1.05), breaks = seq(0,1,length.out=3)) +
  scale_fill_scico_d(palette=col_strat, begin=.1, end=.9) + 
  scale_color_scico_d(palette=col_strat, begin=.1, end=.9)

ggsave("manuscript/figures/lcm_weights.jpg", plot=lcm_weights, 
       width = 17.62, height=9/1.1, units="cm",dpi=900)
ggsave("manuscript/figures/lcm_weights.pdf", plot=lcm_weights, 
       width = 17.62, height=9/1.1, units="cm",dpi=900)



# Appendix: Risk Analysis -------------------------------------------------


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


### figure ------------------------------------------------------------------

# plot means for safe vs. risky short-term winning problems

choice_prop_ind_risk <- simple_dat |> 
  mutate(short_choice = if_else(choice==short_option,1,0) ,
         complexity = if_else(complexity=='simple', 'Simple', 'Complex'), 
         switching = case_match(switching, 
                                'none' ~ "Parallel" ,
                                'frequent' ~ 'Frequent', 
                                'infrequent' ~ 'Infrequent')) |>
  group_by(participant, short_type, switching) |> 
  summarize(choice_prop = mean(short_choice))


figure_risk <- lm_risk_effect_summaries |> 
  filter(!str_detect(.variable, "d")) |> 
  mutate(short_type = case_when(str_detect(.variable, "_R_") ~ "risky" ,
                                str_detect(.variable, "_S_") ~ "safe" ,
                                TRUE ~ NA_character_) ,
         switching = case_when(str_detect(.variable, "F") ~ "Frequent"  ,
                               str_detect(.variable, "I") ~ "Infrequent" ,
                               str_detect(.variable, "P") ~ "Parallel" ,
                               TRUE ~ NA_character_)) |> 
  ggplot(aes(x=factor(switching, levels=c("Infrequent", "Frequent", "Parallel")), y=.value, 
             color=factor(switching, levels=c("Infrequent", "Frequent", "Parallel")))) +
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
  labs(x = "Sampling", 
       y = "Short-term Winner Choice Proportion",
       color = "Sampling",
       linetype = "Short-Term Option") +
  theme_bw()

ggsave("manuscript/figures/risk_preference.eps", plot=figure_risk,
       width = 17.62, height=9/1.1, units="cm",dpi=600, device = cairo_ps)
ggsave("manuscript/figures/risk_preference.jpg", plot=figure_risk,
       width = 17.62, height=9/1.1, units="cm",dpi=900)
ggsave("manuscript/figures/risk_preference.pdf", plot=figure_risk,
       width = 17.62, height=9/1.1, units="cm",dpi=600, device = cairo_pdf)


