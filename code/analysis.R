# preparation -------------------------------------------------------------

# packages 
pacman::p_load(tidyverse , 
               brms ,
               tidybayes,
               ggdist , 
               scico,
               patchwork,
               tidytext
               )

# data
dat <- read_csv('data/clean/choices_long.csv')
part_to_exclude <- c("0319_3_17_1", "0320_5_17_1", "0321_8_13_3")
dat <- dat |> filter(!(participant %in% part_to_exclude))

# ggplot
three_col <- c('#DF9D56', '#6EB5E1', '#572948')
three_lvl <- c('infrequent', 'frequent', 'none')
three_lab <- c('infrequent', 'frequent', 'parallel')


# descriptive -------------------------------------------------------------

dat |> 
  distinct(participant, .keep_all = T) |> 
  mutate(age = if_else(age > 1900, 2025-age, age)) |> 
  summarise(n = n(), 
            m = mean(age, na.rm = T), 
            sd = sd(age, na.rm = T), 
            min = min(age, na.rm =T),
            max = max(age, na.rm =T)
            )


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
               C = dat$short_term_choice)
  
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

figure_group_means <- lm_effect_summaries |> 
  filter(!str_detect(.variable, "d")) |> 
  mutate(complexity = case_when(str_detect(.variable, "_C_") ~ "complex" ,
                                str_detect(.variable, "_S_") ~ "simple" ,
                                TRUE ~ NA_character_) ,
         switching = case_when(str_detect(.variable, "F") ~ "frequent"  ,
                               str_detect(.variable, "I") ~ "infrequent" ,
                               str_detect(.variable, "P") ~ "parallel" ,
                               TRUE ~ NA_character_)) |> 
  ggplot(aes(x=factor(complexity, levels=c("simple", "complex")), y=.value, color=factor(switching, levels=c("infrequent", "frequent", "parallel")))) +
  geom_linerange(aes(ymin = .lower, ymax = .upper), position = position_dodge(width=.5), linewidth=1) +
  geom_point(position = position_dodge(width=.5), size = 3) + 
  scale_color_manual(values=three_col) +
  #scale_y_continuous(limits = c(0,1)) +
  labs(x = "Complexity", 
       y = "Short-term Winner Choice Proportion",
       color = "Switching") +
  theme_bw()
figure_group_means

ggsave("figures/group_means.eps",
       width = 17.62, height=9/0.7, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/group_means.png",
       width = 17.62, height=9/0.7, units="cm",dpi=900)


## risk preference (simple problems) ---------------------------------------------------------

simple_dat <- dat |> filter(complexity=='simple')

lm_risk_dat <- list(PART = simple_dat$participant ,
                PROB = simple_dat$problem_unique , 
                S = factor(simple_dat$switching) , 
                R = factor(simple_dat$short_type) , 
                C = simple_dat$short_term_choice)

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
         d_P = m_R_P - m_S_P 
  ) 

lm_risk_effect_summaries <- lm_risk_effect_samples |> 
  gather_draws(m_S_F, m_S_I, m_S_P , 
               m_R_F, m_R_I, m_R_P , 
               d_S_FI, d_S_FP, d_S_IP, 
               d_R_FI, d_R_FP, d_R_IP ,
               d_F, d_I, d_P) |> 
  mean_qi()

# plot means for safe vs. risky short-term winning problems

figure_risk <- lm_risk_effect_summaries |> 
  filter(!str_detect(.variable, "d")) |> 
  mutate(risk = case_when(str_detect(.variable, "_R_") ~ "risky" ,
                                str_detect(.variable, "_S_") ~ "safe" ,
                                TRUE ~ NA_character_) ,
         switching = case_when(str_detect(.variable, "F") ~ "frequent"  ,
                               str_detect(.variable, "I") ~ "infrequent" ,
                               str_detect(.variable, "P") ~ "parallel" ,
                               TRUE ~ NA_character_)) |> 
  ggplot(aes(x=factor(switching, levels=c("infrequent", "frequent", "parallel")), y=.value, 
             color=factor(switching, levels=c("infrequent", "frequent", "parallel")),
             linetype = risk)) +
  geom_linerange(aes(ymin = .lower, ymax = .upper), position = position_dodge(width=.5), linewidth=1) +
  geom_point(position = position_dodge(width=.5), size = 3) + 
  scale_color_manual(values=three_col) +
  #scale_y_continuous(limits = c(0,1)) +
  labs(x = "Switching", 
       y = "Short-term Winner Choice Proportion",
       color = "Switching",
       linetype = "Short-Term Option") +
  theme_bw()
figure_risk

ggsave("figures/risk_preference.eps",
       width = 17.62, height=9/0.7, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/risk_preferences.png",
       width = 17.62, height=9/0.7, units="cm",dpi=900)



# computational models --------------------------------------------------

current_model <- 'pooled_lcm' 
current_model_summaries <- paste0('fits/summaries/',current_model,'.RData') 

if(file.exists(paste0('fits/summaries/',current_model,'.RData'))){
  load(current_model_summaries, verbose=T)} else{
    print(paste0('File ', current_model_summaries, ' could not be found. First run script code/classification_analysis.R' ))}


## pooling -----------------------------------------------------------------

all_weights |> 
  ggplot(aes(x=switching, y=weight, fill=strategy)) + 
  geom_col(position = 'dodge', col='black') + 
  facet_wrap(~factor(complexity, levels = c('simple', 'complex')), scales='free', drop=TRUE, nrow=1) +
  labs(x='Strategy', 
       y='Weight', 
       fill='Decision Strategy') +
  theme_bw() +
  scale_y_continuous(breaks = seq(0,1,length.out=3)) 

all_errors


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

dat |> 
  filter(group==6) |> 
  ggplot(aes(x=problem_match,y=short_term_choice)) +
  geom_point() + 
  geom_line() +
  facet_wrap(~participant) 







