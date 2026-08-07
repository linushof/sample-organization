# preparation -------------------------------------------------------------

# packages 
pacman::p_load(tidyverse , 
               brms ,
               tidybayes,
               ggdist , 
               scico,
               patchwork,
               tidytext,
               beeswarm
               )

# data
choices <- read_csv('data/choices.csv')
problems <- read_csv('data/problems.csv')
dat <- choices |> left_join(problems, by=join_by(complexity, problem, problem_unique))

# ggplot
three_col <- c('#DF9D56', '#6EB5E1', '#572948')
three_lvl <- c('infrequent', 'frequent', 'none')
three_lab <- c('infrequent', 'frequent', 'parallel')

# descriptive -------------------------------------------------------------

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

# linear models ----------------------------------------------------

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
          file = 'fits/linear_model.rds')

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

ggsave("figures/group_means.eps", plot=figure_group_means,
       width = 17.62, height=9/0.9, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/group_means.png", plot=figure_group_means,
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
               file = 'fits/linear_model_risk.rds')

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
figure_risk

ggsave("figures/risk_preference.eps", plot=figure_risk,
       width = 17.62, height=9/0.9, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/risk_preference.png", plot=figure_risk,
       width = 17.62, height=9/0.9, units="cm",dpi=900)



# computational models --------------------------------------------------

models <- c('complete_pooling', 'pooled_lcm', 'pooled_lcm2', 'hierarchical_lcm')
current_model <- models[3]  
current_model_summaries <- paste0('fits/summaries/',current_model,'.RData') 

if(file.exists(paste0('fits/summaries/',current_model,'.RData'))){
  load(current_model_summaries, verbose=T)} else{
    print(paste0('File ', current_model_summaries, ' could not be found. First run script code/classification_analysis.R' ))}


## complete pooling --------------------------------------------------------

all_weights |> 
  ggplot(aes(x=switching, y=weight, fill=strategy)) + 
  geom_col(position = 'dodge', col='black') + 
  facet_wrap(~factor(complexity, levels = c('simple', 'complex')), scales='free', drop=TRUE, nrow=1) +
  labs(x='Strategy', 
       y='Weight', 
       fill='Decision Strategy') +
  theme_bw() +
  scale_fill_manual(values=three_col) +
  scale_y_continuous(breaks = seq(0,1,length.out=3)) 

all_errors |> 
  ggplot(aes(x=factor(complexity, levels = c('simple', 'complex')), y=error, col=switching)) + 
  geom_point(size=2) + 
  scale_color_manual(values=three_col) +
  labs(x = 'Complexity') +
  theme_bw()

## pooling -----------------------------------------------------------------

all_weights |> 
  ggplot(aes(x=switching, y=weight, fill=strategy)) + 
  geom_col(position = 'dodge', col='black') + 
  #geom_hline(yintercept = 1/3, linetype='dashed') + 
  facet_wrap(~factor(complexity, levels = c('simple', 'complex')), drop=TRUE, nrow=1) +
  labs(x='Strategy', 
       y='Weight', 
       fill='Decision Strategy') +
  theme_bw() +
  scale_y_continuous(breaks = seq(0,1,length.out=3)) 

ggsave("figures/strategy_classification_weights.eps",
       width = 17.62, height=9/1.2, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/strategy_classification_weights.png",
       width = 17.62, height=9/1.2, units="cm",dpi=900)

all_errors |> 
  ggplot(aes(x=switching, y=error, fill=switching)) + 
  geom_boxplot(alpha = .8) + 
  geom_hline(yintercept = .5, linetype='dashed') + 
  facet_wrap(~factor(complexity, levels = c('simple', 'complex')),  drop=TRUE, nrow=1) +
  scale_fill_manual(values=three_col) +
  scale_y_continuous(limits=c(0,.55), breaks=seq(0,.5,length.out=3)) +
  theme_bw() 
  
ggsave("figures/strategy_classification_error.eps",
       width = 17.62, height=9/1.2, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/strategy_classification_error.png",
       width = 17.62, height=9/1.2, units="cm",dpi=900)


## hierarchical ------------------------------------------------------------

all_errors_s <- all_errors |> filter(complexity=='simple') 
plot_weights_simple <- all_weights |> 
  filter(complexity=='simple') |> 
  ggplot(aes(x=reorder_within(participant, by=weight_round, within=switching), y=weight, fill=strategy)) + 
  geom_col(position = 'fill', col='black') + 
  geom_point(data=all_errors_s, aes(y=error_transformed)) + 
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

all_errors_c <- all_errors |> filter(complexity=='complex') 
plot_weights_complex <- all_weights |> 
  filter(complexity=='complex') |> 
  mutate(weight_round = if_else(strategy=='Roundwise', weight, 1-weight)) |> 
  ggplot(aes(x=reorder_within(participant, by=weight_round, within=switching), y=weight, fill=strategy)) + 
  geom_col(position = 'fill', col='black') + 
  geom_point(data=all_errors_c, aes(y=error_transformed)) + 
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
lcm_plot

ggsave("figures/strategy_classification.eps",
       width = 17.62, height=9/0.7, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/strategy_classification.png",
       width = 17.62, height=9/0.7, units="cm",dpi=900)



# additional stuff --------------------------------------------------


## individual differences --------------------------------------------------






|> 
  ggplot(aes(x=switching, y=choice_prop)) + 
  facet_wrap(~factor(complexity, levels=c('simple', 'complex')), nrow=1) + 
  geom_beeswarm(cex = 3, size = 1.5, color = 'gray') + 
  geom_
  #geom_hline(yintercept = .5) + 
  theme_bw()

  geom_jitter(size = 2) + 
  geom_hline(yintercept = .5 ) + 
  theme_bw()



dat |> 
  group_by(participant, complexity, switching) |> 
  summarize(choice_prop = mean(choice)) |> 
  ggplot(aes(x=choice_prop, fill=switching)) + 
  facet_wrap(~factor(complexity, levels=c('simple', 'complex')), nrow=1) + 
  geom_density(alpha=.5, position='identity') + 
  geom_vline(xintercept = .5 ) + 
  theme_bw()



dat |> 
  group_by(participant, complexity, long_type, switching) |> 
  summarize(choice_prop = mean(choice)) |> 
  ggplot(aes(x=switching, y=choice_prop, col=switching)) + 
  facet_grid(long_type~factor(complexity, levels=c('simple', 'complex'))) + 
  geom_jitter(size=2) + 
  geom_hline(yintercept = .5 ) + 
  theme_bw()


dat |> 
  group_by(participant, complexity, long_type,switching) |> 
  summarize(choice_prop = mean(choice)) |> 
  ggplot(aes(x=choice_prop, color=switching)) + 
  facet_grid(long_type~factor(complexity, levels=c('simple', 'complex'))) + 
  geom_density(alpha = .5) + 
  geom_vline(xintercept = .5 ) + 
  theme_bw()





# replace with random effects of subjects
group_means <- dat |> 
  group_by(complexity, switching) |> 
  summarize(choice_prop = mean(choice)) |> 
  mutate(switching = if_else(switching=='none', "parallel" , switching))




dat |> 
  group_by(participant, complexity, switching) |> 
  summarize(choice_prop = mean(choice)) |> 
  mutate(switching=if_else(switching=='none', 'parallel', switching)) |> 
  ggplot(aes(x=switching, y=choice_prop, group=switching), color='lightgray', alpha = .4) + 
  facet_wrap(~factor(complexity, levels=c('simple', 'complex')), nrow=1) + 
  geom_beeswarm(cex = 2) + 
  geom_point(data=group_means, aes(color=switching), size = 2) +
  geom_hline(yintercept = .5 ) + 
  theme_bw()





ggplot(aes(x=factor(complexity, levels=c("simple", "complex")), y=.value, color=factor(switching, levels=c("infrequent", "frequent", "parallel")))) +
  geom_linerange(aes(ymin = .lower, ymax = .upper), position = position_dodge(width=.5), linewidth=1) +
  geom_point(position = position_dodge(width=.5), size = 3) + 
  # geom_beeswarm(data = ind_choice_prop, aes(y=choice_prop), alpha = .5) +
  scale_color_manual(values=three_col) +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Complexity", 
       y = "Short-term Winner Choice Proportion",
       color = "Switching") +
  theme_bw()
figure_group_means





