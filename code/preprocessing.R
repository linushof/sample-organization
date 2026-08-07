# preparation -------------------------------------------------------------

# pkgs
pacman::p_load(tidyverse)

# data
probs <- read_csv('data/stimuli/observed_samples/test_summaries.csv')
parts <- list.files(path = 'data/raw/', pattern = "*.csv", full.names=T, include.dirs=T, recursive=T) # PsychpPy outputs
# length(parts) # 199

# exclude participants that did not follow the instructions
exclude <- c("0319_3_17_1", "0320_5_17_1", "0321_8_13_3")
valid_parts <- parts[!grepl(paste(exclude, collapse = "|"), parts)]
#length(valid_parts) # 196 

# choices -----------------------------------------------------------------

all_dat <- vector('list', length(valid_parts))
for (part in valid_parts) {
  
  data <- read_csv(part, col_types = cols())
  
  # Extract the necessary information and convert to character
  participant <- unique(data$participant)
  gender <- data$gender.response[!is.na(data$gender.response) & data$gender.response != ""][1]
  age <- data$age.text[!is.na(data$age.text) & data$age.text != ""][1]
  group <- data$group.text[!is.na(data$group.text) & data$group.text != ""][1]
  resp <- data$choice[!is.na(data$choice) & data$choice != ""] # 1 = option 1; 2 = option 2
  
  dat <- expand_grid(participant=participant, 
                     gender=gender, 
                     age=if_else(age > 1900, 2025-age, age),
                     group=group,
                     choice=if_else(resp==1, 1,0)) |> 
    mutate(complexity = case_when(group %in% c(1,2,5) ~ 'simple', 
                                  group %in% c(3,4,6) ~ 'complex') ,
           present = case_when(group %in% 1:4 ~ 'sequential' , 
                               group %in% 5:6 ~ 'parallel') , 
           switching = case_when(group %in% c(1,3) ~ 'infrequent', 
                                 group %in% c(2,4) ~ 'frequent',
                                 group %in% c(5,6) ~ 'none') , 
           problem=1:20,
           problem_unique = as.factor(if_else(complexity=='simple', problem, problem+20)) ,
           problem = as.factor(problem)) |> 
    select(participant:group, complexity, present, switching, problem, problem_unique, choice)
  
  all_dat[[part]] <- dat
  
}    

choices <- bind_rows(all_dat)
write_csv(choices, 'data/choices.csv')

# problems ---------------------------------------------

problems <- probs |> 
  separate(problem, c('complexity', 'problem'), sep='_') |>  
  mutate(problem = as.double(str_remove(problem, 'problem')) , 
         problem_unique = as.factor(if_else(complexity=='simple', problem, problem+20)) ,
         problem = as.factor(problem),
         short_option = if_else(short_option=='option1',1,0), 
         long_option = if_else(long_option=='option1',1,0)) |> 
  select(complexity, problem, problem_unique, everything())
problems

write_csv(problems, 'data/problems.csv')


# samples -----------------------------------------------------------------

# Get file list
file_list <- list.files(
  path = 'data/stimuli/samples/test', 
  pattern = '^(simple|complex)_problem\\d+_(frequent|infrequent|parallel)', 
  full.names = TRUE
)

# Read files and extract the number
samples <- tibble(file = file_list) |> 
  mutate(
    problem_match = as.numeric(str_extract(basename(file), "(?<=(simple|complex)_problem)\\d+(?=_(frequent|infrequent|parallel))")),
    complexity = str_extract(basename(file), '(complex|simple)') ,
    cond = str_extract(basename(file), '(frequent|infrequent|parallel)') , 
    present = if_else(cond=='parallel', 'parallel', 'sequential'),
    switching = if_else(present=='parallel', 'none', str_extract(basename(file), '(frequent|infrequent)') ) ,
    problem_unique = if_else(complexity=='complex', problem_match, problem_match+20) , 
    data = map(file, read_csv) ) |> 
  unnest(cols = data) |>  
  group_by(problem_match, complexity, switching) |>  
  mutate(sample = row_number()) |>  
  ungroup() |>  
  mutate(group = as.factor(case_when(complexity == 'simple' & switching == 'infrequent' ~ 1 , 
                           complexity == 'simple' & switching == 'frequent' ~ 2 ,
                           complexity == 'complex' & switching == 'infrequent' ~ 3 ,
                           complexity == 'complex' & switching == 'frequent' ~ 4 ,
                           complexity == 'simple' & switching == 'none' ~ 5 ,
                           complexity == 'complex' & switching == 'none' ~ 6))) |>  
  select(problem_match, problem_unique, present, switching, complexity, group, sample, option1, option2)

write_csv(samples, 'data/stimuli/observed_samples/test_all.csv')
