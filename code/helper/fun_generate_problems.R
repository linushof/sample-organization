# load packages

library(tidyverse) 
library(gtools)


# problem generator -------------------------------------------------------


generate_problems_beta1.1 <- function(n=20 , # number of problems per condition
                              range=c(1, 99) , # possible outcome interval
                              pdirichlet = c(80, 120, 300) , # parameter of Dirichlet distribution (defaults create skewed distributions) 
                              seed=123 # RNG seed for reproducibility
                              ){
  
  set.seed(seed)
  
  problems <- tibble() # contains accepted proposals
  features <- c("short_o1", "short_o2", "short_o3", 
                "short_p1", "short_p2", "short_p3" ,
                "long_o1", "long_o2", "long_o3" , 
                "long_p1", "long_p2", "long_p3")
  
  while(nrow(problems) < n){ # to generate proposals until n are accepted 
    
    # for each option, draw length(pdirichlet) probabilities from a Dirichlet distribution
    short_probs <-  rdirichlet(n=1, alpha=pdirichlet) %>% sort() %>% round(3) # short-term winner: left-skewed
    long_probs <-  rdirichlet(n=1, alpha=pdirichlet) %>% sort(decreasing=T) %>% round(3) #  long-term winner: right skewed
    
    
    # draw outcomes 
    outs <- sample(range[1]:range[2], size=5, replace=F) %>% round()
    long_outs <- c(max(outs) , sample(outs[outs!=min(outs) & outs!=max(outs)], 2, replace = F)) %>% sort() # ensures that largest (smallest) outcome belongs to long (short)  
    long_ev = long_outs[1]*long_probs[1] + long_outs[2]*long_probs[2] + long_outs[3]*long_probs[3]
    short_outs <- c(outs[!outs %in% long_outs], sample(1.1*long_ev:range[2], size=1, replace = F)) %>% round() %>%  sort() ## correct logic?
    
    # to prevent dominant options -> still necessary ?? (see line above)
    if(max(short_outs) < min(long_outs)){
      newShort <- min(long_outs)
      newLong <- max(short_outs)
      short_outs[3] <- newShort
      long_outs[1] <- newLong
    }
    
    # compute joint probabilities and probability winning pairwise comparison (winprob)
    short <- tibble(short_outs, short_probs)
    long  <- tibble(long_outs, long_probs)
    joint <- expand_grid(short, long) %>% 
      mutate(winner = if_else(short_outs > long_outs, 1, 0), 
             joint = short_probs*long_probs) %>% 
      group_by(winner) %>% 
      summarize(prob = sum(joint))
    winprob <- round(joint[[2,2]],3)
    
    # proposal summary
    proposal <- tibble(feature = features , 
                       value = c(short_outs, short_probs, long_outs, long_probs)) %>% 
      pivot_wider(names_from = feature, values_from = value) %>% 
      mutate(winprob = winprob , 
             short_ev = round(short_o1*short_p1 + short_o2*short_p2 + short_o3*short_p3, 3) , # compute expected value
             long_ev = round(long_o1*long_p1 + long_o2*long_p2 + long_o3*long_p3, 3) ,
             ratio = round(long_ev/short_ev, 2)) 
    
    # check validity and add
    problem <- proposal %>% filter(ratio >= 1.15 ,
                                   winprob >= .65 , 
                                   short_p1 >= .1 , 
                                   long_p3 >= .1 , 
                                   short_o3 > long_ev , 
                                   short_ev > long_o1
                                   )
    
    problems <- bind_rows(problems, problem)
    
  } # close while loop
  
  outputs <- vector("list", length = 3)
  names(outputs) <- c("complex", "simple_short", "simple_long")
  
  
  # further processing 
  
  complex <- problems %>% 
    mutate(id = row_number() , 
           type = "complex") %>% 
    select(id, type, short_ev, long_ev, ratio, winprob, everything())
  outputs[[1]] <- complex
  
  # reduce complex case to safe-risky versions
  
  ## with left skew (risky is short-term winner, safe is long-term winner)
  
  simple_short <- complex %>% mutate(type = "simple-short" , 
                                     long_o1_new = long_ev ,
                                     long_p1_new = 1 ,
                                     short_p3_new = winprob ,  
                                     short_p1_new = 1-short_p3_new , 
                                     short_o1_new = round(  ( short_ev - (short_p3_new*short_o3) )/ short_p1_new ) ,
                                     short_ev_new = round ( short_p1_new*short_o1_new + short_p3_new*short_o3 , 3) ,
                                     ratio_new = round(long_ev/short_ev_new, 2) , 
                                     winprob_new = short_p3_new) %>% 
    select(id, type , short_ev_new, long_ev, ratio_new, winprob_new, short_o1_new, short_o3, short_p1_new, short_p3_new, long_o1_new, long_p1_new)
  outputs[[2]] <- simple_short
  
  # with right skew (safe is short-term winner, risky is long-term winner)
  
  simple_long <- complex %>% mutate(type = "simple-long" , 
                                    short_o1_new = short_ev , 
                                    short_p1_new = 1 , 
                                    long_p1_new = winprob , 
                                    long_p3_new = 1-long_p1_new , 
                                    long_o3_new = round ( ( long_ev - (long_p1_new*long_o1) )/long_p3_new ) ,
                                    long_ev_new = round ( long_p1_new*long_o1 + long_p3_new *long_o3_new , 3) , 
                                    ratio_new = round(long_ev_new/short_ev, 2) ,
                                    winprob_new = long_p1_new) %>%
    select(id, type , short_ev, long_ev_new, ratio_new, winprob_new, short_o1_new, short_p1_new, long_o1, long_o3_new, long_p1_new, long_p3_new)
  outputs[[3]] <- simple_long
  
  
  return(outputs)
  
}