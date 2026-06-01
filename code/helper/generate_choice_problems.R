
# test problems -----------------------------------------------------------

## generate problems -------------------------------------------------------

source("code/helper/fun_generate_problems.R")
run <- generate_problems_beta1.1(n=20, seed=2539)


## merge data sets ------------------------------------------------------

# to obtain simple problems with the same columns/format as for complex problems

## left-skewed problems

simple_short_full <- run$simple_short %>% 
  rename(short_ev = short_ev_new , 
         ratio = ratio_new , 
         winprob = winprob_new , 
         short_o1 = short_o1_new , 
         short_p1 = short_p1_new , 
         short_p3 = short_p3_new ,
         long_o1 = long_o1_new , 
         long_p1 = long_p1_new) %>% 
  mutate(short_o2 = 0 , 
         short_p2 = 0 , 
         long_o2 = 0 ,
         long_p2 = 0 , 
         long_o3 = 0 , 
         long_p3 = 0 ) %>% 
  select(id, type, short_ev, long_ev, ratio, winprob, 
         short_o1, short_o2, short_o3, 
         short_p1, short_p2, short_p3, 
         long_o1, long_o2, long_o3, 
         long_p1, long_p2, long_p3)

## right skewed options
  
simple_long_full <- run$simple_long %>% 
  rename(long_ev = long_ev_new , 
         ratio = ratio_new , 
         winprob = winprob_new , 
         short_o1 = short_o1_new , 
         short_p1 = short_p1_new , 
         long_p1 = long_p1_new , 
         long_o3 = long_o3_new , 
         long_p3 = long_p3_new) %>% 
  mutate(short_o2 = 0 , 
         short_p2 = 0 ,
         short_o3 = 0 , 
         short_p3 = 0 , 
         long_o2 = 0 ,
         long_p2 = 0) %>% 
  select(id, type, short_ev, long_ev, ratio, winprob, 
         short_o1, short_o2, short_o3, 
         short_p1, short_p2, short_p3, 
         long_o1, long_o2, long_o3, 
         long_p1, long_p2, long_p3)                     

## merge all data sets

all_problems <- bind_rows(run$complex, simple_short_full, simple_long_full)

## store data sets ---------------------------------------------------------

# .csv
write_csv(all_problems, "data/stimuli/all_problems.csv" ) # all
write_csv(run$complex, "data/stimuli/complex_problems.csv" ) # only complex
write_csv(run$simple_short, "data/stimuli/simple_short_problems.csv" ) # only left-skewed
write_csv(run$simple_long, "data/stimuli/simple_long_problems.csv" ) # only right-skewed


# training ----------------------------------------------------------------

training <- generate_problems_beta1.1(n=3, seed=78915)


## merge data sets ------------------------------------------------------


simple_short_full <- training$simple_short %>% 
  rename(short_ev = short_ev_new , 
         ratio = ratio_new , 
         winprob = winprob_new , 
         short_o1 = short_o1_new , 
         short_p1 = short_p1_new , 
         short_p3 = short_p3_new ,
         long_o1 = long_o1_new , 
         long_p1 = long_p1_new) %>% 
  mutate(short_o2 = 0 , 
         short_p2 = 0 , 
         long_o2 = 0 ,
         long_p2 = 0 , 
         long_o3 = 0 , 
         long_p3 = 0 ) %>% 
  select(id, type, short_ev, long_ev, ratio, winprob, 
         short_o1, short_o2, short_o3, 
         short_p1, short_p2, short_p3, 
         long_o1, long_o2, long_o3, 
         long_p1, long_p2, long_p3)

## right skewed options

simple_long_full <- training$simple_long %>% 
  rename(long_ev = long_ev_new , 
         ratio = ratio_new , 
         winprob = winprob_new , 
         short_o1 = short_o1_new , 
         short_p1 = short_p1_new , 
         long_p1 = long_p1_new , 
         long_o3 = long_o3_new , 
         long_p3 = long_p3_new) %>% 
  mutate(short_o2 = 0 , 
         short_p2 = 0 ,
         short_o3 = 0 , 
         short_p3 = 0 , 
         long_o2 = 0 ,
         long_p2 = 0) %>% 
  select(id, type, short_ev, long_ev, ratio, winprob, 
         short_o1, short_o2, short_o3, 
         short_p1, short_p2, short_p3, 
         long_o1, long_o2, long_o3, 
         long_p1, long_p2, long_p3)                     

## merge all data sets

all_problems <- bind_rows(training$complex, simple_short_full, simple_long_full)

## store data sets ---------------------------------------------------------

write_csv(all_problems, "data/stimuli/training_all_problems.csv" ) # all
write_csv(training$complex, "data/stimuli/training_complex_problems.csv" ) # only complex
write_csv(training$simple_short, "data/stimuli/training_simple_short_problems.csv" ) # only left-skewed
write_csv(training$simple_long, "data/stimuli/training_simple_long_problems.csv" ) # only right-skewed

